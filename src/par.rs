//! Port of `main.Par`: Beagle-style `key=value` command-line parameters.

use std::collections::HashMap;
use std::process::exit;

#[derive(Clone, Debug)]
#[allow(dead_code)] // phasing parameters accepted for CLI compatibility
pub struct Par {
    // data parameters
    pub gt: String,
    pub reff: Option<String>,
    pub out: String,
    pub map: Option<String>,
    pub chrom: Option<ChromInterval>,
    pub excludesamples: Option<String>,
    pub excludemarkers: Option<String>,

    // phasing parameters (accepted; used once phasing is ported)
    pub burnin: i32,
    pub iterations: i32,
    pub phase_states: i32,

    // imputation parameters
    pub impute: bool,
    pub imp_states: usize,
    pub imp_segment: f32,
    pub imp_step: f32,
    pub imp_nsteps: usize,
    pub cluster: f32,
    pub ap: bool,
    pub gp: bool,

    // general parameters
    pub em: bool,
    pub initial_lr: f32,
    pub step_scale: f32,
    pub rare: f32,
    pub ne: f32,
    pub err: f32,
    pub window: f32,
    pub window_markers: usize,
    pub overlap: f32,
    pub buffer: f32,
    pub seed: i64,
    pub nthreads: usize,
    /// Number of imputation replicates averaged into the output (rusty-beagle
    /// extension; 1 = plain Beagle behavior, bit-identical to Java).
    pub ensemble: usize,

    /// Species preset applied (rusty-beagle extension; not a Java Beagle
    /// parameter). Records the parameters it filled in, for the run log.
    pub preset: Option<PresetInfo>,
}

/// Defaults applied by `preset=cattle`. The single change with a large,
/// reproducible effect is `ne`: Beagle's default of 100,000 suits human
/// panels, but cattle breeds have an effective population size near 100,
/// and the imputation HMM's state-switch rate scales with
/// `0.04 * ne / nRefHaps`, so the default massively over-switches.
///
/// ne=1000 sits on the accuracy plateau at every tested panel size. On the
/// public synbreedData panel of 500 real dairy bulls (LD-chip targets
/// imputed to the full 7,250-marker panel, sire-family-aware
/// cross-validation), mean masked-marker dosage r2 rises from 0.36 at the
/// default to 0.77; on coalescent panels simulated under the MacLeod et al.
/// (2013) Holstein demography the gain ranges from +0.05 (20,000-haplotype
/// reference) to +0.50 (800-haplotype reference). Secondary parameters
/// (err, imp-states, cluster, window) moved accuracy by at most 0.005 and
/// not consistently across datasets, so the preset leaves them alone. See
/// tests/bovine/ for the harness that reproduces all of this.
const CATTLE_PRESET: &[(&str, &str)] = &[("ne", "1000")];

/// A preset supplies defaults for parameters the user did not set
/// explicitly; explicit `key=value` arguments always win.
#[derive(Clone, Debug)]
pub struct PresetInfo {
    pub name: String,
    /// (parameter, value) pairs the preset actually applied
    pub applied: Vec<(String, String)>,
}

#[derive(Clone, Debug)]
pub struct ChromInterval {
    pub chrom: String,
    pub start: i32, // inclusive; -1 if unbounded
    pub end: i32,   // inclusive; i32::MAX if unbounded
}

impl ChromInterval {
    pub fn parse(s: &str) -> Option<ChromInterval> {
        if s.is_empty() {
            return None;
        }
        let (chrom, rest) = match s.rfind(':') {
            Some(i) => (&s[..i], Some(&s[i + 1..])),
            None => (s, None),
        };
        if chrom.is_empty() {
            return None;
        }
        match rest {
            None => Some(ChromInterval {
                chrom: chrom.to_string(),
                start: -1,
                end: i32::MAX,
            }),
            Some(r) => {
                let dash = r.find('-')?;
                let start_str = &r[..dash];
                let end_str = &r[dash + 1..];
                let start = if start_str.is_empty() {
                    -1
                } else {
                    start_str.parse::<i32>().ok()?
                };
                let end = if end_str.is_empty() {
                    i32::MAX
                } else {
                    end_str.parse::<i32>().ok()?
                };
                Some(ChromInterval {
                    chrom: chrom.to_string(),
                    start,
                    end,
                })
            }
        }
    }

    pub fn contains(&self, chrom: &str, pos: i32) -> bool {
        self.chrom == chrom && pos >= self.start && (self.end == i32::MAX || pos <= self.end)
    }
}

pub fn usage() -> String {
    format!(
        "Usage: rusty-beagle [arguments]

data parameters ...
  gt=<VCF file with GT FORMAT field>                 (required)
  preset=<species preset: cattle>                    (optional; rusty-beagle
        extension: fills in defaults tuned for the species -- currently
        ne=1000 for cattle -- for any parameter not given explicitly)
  ensemble=<imputation replicates to average>        (default=1; rusty-beagle
        extension: runs the stochastic phasing/imputation K times with
        distinct seeds and averages the allele probabilities, which
        measurably raises dosage-r2 accuracy at ~K-fold run time; 1 keeps
        output bit-identical to Java Beagle)
  ref=<VCF file with phased genotypes>               (optional)
  out=<output file prefix>                           (required)
  map=<PLINK map file with cM units>                 (optional)
  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)
  excludesamples=<file with 1 sample ID per line>    (optional)
  excludemarkers=<file with 1 marker ID per line>    (optional)

imputation parameters ...
  impute=<impute ungenotyped markers (true/false)>   (default=true)
  imp-states=<model states for imputation>           (default=1600)
  cluster=<max cM in a marker cluster>               (default=0.005)
  ap=<print posterior allele probabilities>          (default=false)
  gp=<print posterior genotype probabilities>        (default=false)

general parameters ...
  ne=<effective population size>                     (default=100000)
  err=<allele mismatch probability>                  (default: data dependent)
  window=<window length in cM>                       (default=40.0)
  window-markers=<maximum markers per window>        (default=4000000)
  overlap=<window overlap in cM>                     (default=2.0)
  seed=<random seed>                                 (default=-99999)
  nthreads=<number of threads>                       (default: machine dependent)
"
    )
}

fn get<'a>(map: &mut HashMap<String, String>, key: &str) -> Option<String> {
    map.remove(key)
}

fn parse_or_exit<T: std::str::FromStr>(key: &str, val: &str) -> T {
    val.parse::<T>().unwrap_or_else(|_| {
        eprintln!("ERROR: invalid value for parameter \"{}\": {}", key, val);
        exit(1)
    })
}

fn bool_arg(map: &mut HashMap<String, String>, key: &str, default: bool) -> bool {
    match get(map, key) {
        None => default,
        Some(v) => match v.to_lowercase().as_str() {
            "true" => true,
            "false" => false,
            _ => {
                eprintln!("ERROR: invalid value for parameter \"{}\": {}", key, v);
                exit(1)
            }
        },
    }
}

impl Par {
    pub fn new(args: &[String]) -> Par {
        let mut map: HashMap<String, String> = HashMap::new();
        for arg in args {
            match arg.split_once('=') {
                Some((k, v)) => {
                    if map.insert(k.to_string(), v.to_string()).is_some() {
                        eprintln!("ERROR: duplicate parameter: {}", arg);
                        exit(1);
                    }
                }
                None => {
                    eprintln!("ERROR: invalid parameter (missing '='): {}", arg);
                    exit(1);
                }
            }
        }
        // Species presets fill in defaults for parameters not given
        // explicitly (rusty-beagle extension; see usage()).
        let preset = match get(&mut map, "preset") {
            None => None,
            Some(name) => {
                let defaults: &[(&str, &str)] = match name.as_str() {
                    "cattle" | "bovine" => CATTLE_PRESET,
                    _ => {
                        eprintln!("ERROR: unknown preset \"{}\" (available: cattle)", name);
                        exit(1)
                    }
                };
                let mut applied = Vec::new();
                for (k, v) in defaults {
                    if !map.contains_key(*k) {
                        map.insert(k.to_string(), v.to_string());
                        applied.push((k.to_string(), v.to_string()));
                    }
                }
                Some(PresetInfo { name, applied })
            }
        };

        let gt = get(&mut map, "gt").unwrap_or_else(|| {
            eprintln!("{}", usage());
            eprintln!("ERROR: missing required parameter \"gt\"");
            exit(1)
        });
        let out = get(&mut map, "out").unwrap_or_else(|| {
            eprintln!("{}", usage());
            eprintln!("ERROR: missing required parameter \"out\"");
            exit(1)
        });
        let reff = get(&mut map, "ref");
        let mapfile = get(&mut map, "map");
        let chrom = match get(&mut map, "chrom") {
            None => None,
            Some(s) => match ChromInterval::parse(&s) {
                Some(ci) => Some(ci),
                None => {
                    eprintln!("ERROR: invalid chrom parameter: {}", s);
                    exit(1)
                }
            },
        };
        let excludesamples = get(&mut map, "excludesamples");
        let excludemarkers = get(&mut map, "excludemarkers");

        let burnin = get(&mut map, "burnin").map_or(3, |v| parse_or_exit("burnin", &v));
        let iterations =
            get(&mut map, "iterations").map_or(12, |v| parse_or_exit("iterations", &v));
        let phase_states =
            get(&mut map, "phase-states").map_or(280, |v| parse_or_exit("phase-states", &v));

        let impute = bool_arg(&mut map, "impute", true);
        let imp_states =
            get(&mut map, "imp-states").map_or(1600, |v| parse_or_exit("imp-states", &v));
        let imp_segment =
            get(&mut map, "imp-segment").map_or(6.0, |v| parse_or_exit("imp-segment", &v));
        let imp_step = get(&mut map, "imp-step").map_or(0.1, |v| parse_or_exit("imp-step", &v));
        let imp_nsteps =
            get(&mut map, "imp-nsteps").map_or(7, |v| parse_or_exit("imp-nsteps", &v));
        let cluster = get(&mut map, "cluster").map_or(0.005, |v| parse_or_exit("cluster", &v));
        let ap = bool_arg(&mut map, "ap", false);
        let gp = bool_arg(&mut map, "gp", false);

        let em = bool_arg(&mut map, "em", true);
        let initial_lr =
            get(&mut map, "initial-lr").map_or(100_000.0, |v| parse_or_exit("initial-lr", &v));
        let step_scale =
            get(&mut map, "step-scale").map_or(3.0, |v| parse_or_exit("step-scale", &v));
        let rare = get(&mut map, "rare").map_or(0.002, |v| parse_or_exit("rare", &v));
        let ne = get(&mut map, "ne").map_or(100_000.0, |v| parse_or_exit("ne", &v));
        // D_ERR = -Float.MIN_VALUE signals "data dependent"
        let err = get(&mut map, "err").map_or(-f32::MIN_POSITIVE, |v| parse_or_exit("err", &v));
        let window = get(&mut map, "window").map_or(40.0, |v| parse_or_exit("window", &v));
        let window_markers = get(&mut map, "window-markers")
            .map_or(4_000_000, |v| parse_or_exit("window-markers", &v));
        let overlap = get(&mut map, "overlap").map_or(2.0, |v| parse_or_exit("overlap", &v));
        let buffer = get(&mut map, "buffer").map_or(1.0, |v| parse_or_exit("buffer", &v));
        let seed = get(&mut map, "seed").map_or(-99999, |v| parse_or_exit("seed", &v));
        let nthreads = match get(&mut map, "nthreads") {
            None => std::thread::available_parallelism().map_or(1, |n| n.get()),
            Some(v) => parse_or_exit("nthreads", &v),
        };
        let ensemble: usize =
            get(&mut map, "ensemble").map_or(1, |v| parse_or_exit("ensemble", &v));
        if ensemble < 1 {
            eprintln!("ERROR: ensemble must be >= 1");
            exit(1);
        }
        // Accept and ignore documented phasing-only params so Beagle command
        // lines run unchanged.
        for k in ["truth", "ped"] {
            map.remove(k);
        }
        if !map.is_empty() {
            let mut keys: Vec<&String> = map.keys().collect();
            keys.sort();
            eprintln!("{}", usage());
            eprintln!("ERROR: unrecognized parameter(s): {:?}", keys);
            exit(1);
        }
        Par {
            gt,
            reff,
            out,
            map: mapfile,
            chrom,
            excludesamples,
            excludemarkers,
            burnin,
            iterations,
            phase_states,
            impute,
            imp_states,
            imp_segment,
            imp_step,
            imp_nsteps,
            cluster,
            ap,
            gp,
            em,
            initial_lr,
            step_scale,
            rare,
            ne,
            err,
            window,
            window_markers,
            overlap,
            buffer,
            seed,
            nthreads,
            ensemble,
            preset,
        }
    }

    pub fn rare(&self) -> f32 {
        self.rare
    }

    pub fn step_scale(&self) -> f32 {
        self.step_scale
    }

    pub fn initial_lr(&self) -> f32 {
        self.initial_lr
    }

    /// `Par.err(nHaps)`: explicit err parameter, or the Li&Stephens allele
    /// mismatch approximation.
    pub fn err_for(&self, n_haps: usize) -> f32 {
        if self.err >= 0.0 {
            self.err
        } else {
            li_stephens_p_mismatch(n_haps)
        }
    }
}

/// `Par.liStephensPMismatch`
pub fn li_stephens_p_mismatch(n_haps: usize) -> f32 {
    let theta = 1.0 / ((n_haps as f64).ln() + 0.5);
    (theta / (2.0 * (theta + n_haps as f64))) as f32
}
