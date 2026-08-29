//! rusty-beagle: a fast Rust port of Beagle 5.5 genotype imputation.
//!
//! Pipeline driver (port of `main.Main` for the phased-target path).

mod bgzf;
mod genmap;
mod hmm;
mod impdata;
mod impibs;
mod impout;
mod impstates;
mod javautil;
mod marker;
mod par;
mod refpanel;
mod vcfio;
mod windows;

use par::Par;
use std::io::Write;
use std::time::Instant;

/// Printed in the `##source` header line.
pub const JAVA_EQUIV_PROGRAM: &str = "rusty-beagle 0.1.0 (port of beagle.27Feb25.75f.jar)";

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        println!("{}", JAVA_EQUIV_PROGRAM);
        println!("{}", par::usage());
        return;
    }
    let par = Par::new(&args);
    if par.window < 1.1 * par.overlap {
        eprintln!(
            "ERROR: The \"window\" parameter must be at least 1.1 times the \"overlap\" parameter"
        );
        std::process::exit(1);
    }
    if par.reff.is_none() {
        eprintln!(
            "ERROR: rusty-beagle currently requires a reference panel (ref=).\n\
             Phasing without a reference panel is not yet ported; use Java Beagle for that."
        );
        std::process::exit(1);
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(par.nthreads)
        .build_global()
        .expect("rayon pool");

    let start = Instant::now();
    let mut log = Log::new(&par.out);
    log.println(JAVA_EQUIV_PROGRAM);
    log.println(&format!("nthreads       : {}", par.nthreads));
    log.println(&format!("gt             : {}", par.gt));
    log.println(&format!("ref            : {}", par.reff.as_deref().unwrap()));
    log.println(&format!("out            : {}", par.out));
    if let Some(m) = &par.map {
        log.println(&format!("map            : {}", m));
    }
    log.println(&format!("seed           : {}", par.seed));

    run(&par, &mut log);

    log.println(&format!(
        "run time       : {:.2} seconds",
        start.elapsed().as_secs_f64()
    ));
    log.close();
}

fn run(par: &Par, log: &mut Log) {
    let exclude_samples = vcfio::read_exclude_file(&par.excludesamples);
    let exclude_markers = vcfio::read_exclude_file(&par.excludemarkers);
    let chrom_filter: Option<&str> = par.chrom.as_ref().map(|c| c.chrom.as_str());
    let genmap = genmap::GeneticMap::new(&par.map, chrom_filter);

    // target reader
    let mut targ_raw = vcfio::open_text(&par.gt);
    let (targ_header, targ_first) =
        vcfio::read_header(&mut targ_raw, &par.gt, &exclude_samples);
    let targ_lines = vcfio::LineSource::new(targ_raw, targ_first);
    let targ_it = windows::TargReader::new(
        targ_header,
        targ_lines,
        exclude_markers.clone(),
        par.chrom.clone(),
    );

    // reference reader
    let ref_path = par.reff.as_deref().unwrap();
    if ref_path.ends_with(".bref") || ref_path.ends_with(".bref3") {
        eprintln!(
            "ERROR: bref/bref3 reference files are not supported by rusty-beagle yet;\n\
             convert with unbref3 to a .vcf.gz file"
        );
        std::process::exit(1);
    }
    let mut ref_raw = vcfio::open_text(ref_path);
    let (ref_header, ref_first) = vcfio::read_header(&mut ref_raw, ref_path, &exclude_samples);
    let ref_lines = vcfio::LineSource::new(ref_raw, ref_first);
    let ref_reader = refpanel::RefReader::new(ref_header, ref_lines, exclude_markers);
    let ref_it = windows::FilteredRefReader::new(ref_reader, par.chrom.clone());

    let targ_samples = targ_it.header.samples.clone();
    let ref_samples = ref_it.samples();
    log.println(&format!(
        "\nReference samples: {}\nStudy     samples: {}",
        ref_samples.len(),
        targ_samples.len()
    ));

    let genmap = std::sync::Arc::new(genmap);
    let sliding = windows::SlidingWindows::new(par, genmap.clone(), targ_it, ref_it);
    let mut bg = windows::BackgroundWindows::spawn(sliding);
    let mut writer =
        impout::WindowWriter::new(&par.out, targ_samples.clone(), par.ap, par.gp);

    let timing = std::env::var("RUSTY_BEAGLE_TIMING").is_ok();
    let mut window_count = 0usize;
    loop {
        let t_read = Instant::now();
        let window = match bg.next_window() {
            Some(w) => w,
            None => break,
        };
        if timing {
            eprintln!("[timing] window read: {:.3}s", t_read.elapsed().as_secs_f64());
        }
        window_count += 1;
        let indices = &window.indices;
        let first = &window.targ_recs[0].marker;
        let last = &window.targ_recs[window.targ_recs.len() - 1].marker;
        log.println(&format!(
            "\nWindow {} ({}:{}-{})\nReference markers: {}\nStudy     markers: {}",
            window.window_index,
            first.chrom(),
            first.pos,
            last.pos,
            indices.n_markers(),
            indices.n_targ_markers()
        ));
        if !window.targ_is_phased() {
            eprintln!(
                "\nERROR: the target VCF contains unphased or missing genotypes.\n\
                 rusty-beagle currently ports Beagle's imputation stage only, which requires\n\
                 a phased, non-missing target (Java Beagle skips phasing for such input, so\n\
                 outputs are directly comparable).  Phase the target first, e.g. with:\n\
                 java -jar beagle.jar gt={} out=phased impute=false",
                par.gt
            );
            std::process::exit(1);
        }
        let impute = indices.n_markers() != indices.n_targ_markers();
        if !impute {
            let m_start = indices.targ_prev_splice;
            let m_end = indices.targ_next_splice;
            writer.print_phased(&window, m_start, m_end);
        } else {
            let t0 = Instant::now();
            let imp_data = impdata::ImpData::new(par, &window, genmap.as_ref(), &targ_samples);
            let t1 = Instant::now();
            let ibs_haps = impibs::ImpIbs::new(&imp_data);
            let t2 = Instant::now();
            let state_probs = hmm::state_probs(&imp_data, &ibs_haps);
            let t3 = Instant::now();
            let m_start = indices.prev_splice;
            let m_end = indices.next_splice;
            writer.print_imputed(&imp_data, &window, m_start, m_end, &state_probs);
            if timing {
                eprintln!(
                    "[timing] impdata: {:.3}s  ibs: {:.3}s  hmm: {:.3}s  output: {:.3}s",
                    (t1 - t0).as_secs_f64(),
                    (t2 - t1).as_secs_f64(),
                    (t3 - t2).as_secs_f64(),
                    t3.elapsed().as_secs_f64()
                );
            }
            log.println(&format!(
                "Imputation time  : {:.2} seconds",
                t0.elapsed().as_secs_f64()
            ));
        }
    }
    writer.close();
    let stats = bg.finish();
    log.println(&format!(
        "\nCumulative Statistics:\nReference markers: {}\nStudy     markers: {}\nWindows:           {}",
        stats.cum_ref_markers, stats.cum_targ_markers, window_count
    ));
}

/// Writes run info to stdout and to `<out>.log` (like `main.RunStats`).
struct Log {
    file: Option<std::fs::File>,
}

impl Log {
    fn new(out_prefix: &str) -> Log {
        let path = format!("{}.log", out_prefix);
        let file = std::fs::File::create(&path).ok();
        Log { file }
    }

    fn println(&mut self, s: &str) {
        println!("{}", s);
        if let Some(f) = &mut self.file {
            let _ = writeln!(f, "{}", s);
        }
    }

    fn close(self) {}
}
