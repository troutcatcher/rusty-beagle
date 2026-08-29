//! VCF text input: gzip/BGZF decoding, header parsing (`vcf.VcfHeader`),
//! sample lists (`vcf.Samples`), and target GT records
//! (`vcf.VcfRecGTParser` with `VcfIt.TO_LOWMEM_GT_REC` semantics).

use crate::marker::Marker;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::sync::Arc;

/// Opens a possibly-gzipped (including BGZF) text file for line reading.
pub fn open_text(path: &str) -> Box<dyn BufRead + Send> {
    let f = File::open(path).unwrap_or_else(|e| {
        eprintln!("ERROR: cannot open file {}: {}", path, e);
        std::process::exit(1)
    });
    let mut reader = BufReader::with_capacity(1 << 20, f);
    let (is_gz, is_bgzf) = {
        let buf = reader.fill_buf().unwrap_or(&[]);
        let gz = buf.len() >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;
        (gz, gz && crate::bgzf::sniff_bgzf(buf))
    };
    if is_bgzf {
        // BGZF: decompress blocks in parallel (like Java's BGZipIt)
        let bgzf = crate::bgzf::ParallelBgzfReader::new(reader);
        Box::new(BufReader::with_capacity(1 << 20, bgzf))
    } else if is_gz {
        let gz = flate2::bufread::MultiGzDecoder::new(reader);
        Box::new(BufReader::with_capacity(1 << 20, gz))
    } else {
        Box::new(reader)
    }
}

/// Port of `vcf.Samples`.
#[derive(Clone, Debug)]
pub struct Samples {
    pub ids: Arc<Vec<String>>,
    pub is_diploid: Arc<Vec<bool>>,
}

impl Samples {
    pub fn len(&self) -> usize {
        self.ids.len()
    }
}

/// Port of `vcf.VcfHeader`.
pub struct VcfHeader {
    pub src: String,
    pub samples: Samples,
    /// maps filtered sample index -> unfiltered sample index
    pub unfiltered_index: Vec<usize>,
    pub n_unfiltered_samples: usize,
}

const SAMPLE_OFFSET: usize = 9;

/// Reads meta-information lines and the header line; returns the header and
/// the first data line (if any).
pub fn read_header(
    reader: &mut dyn BufRead,
    src: &str,
    exclude_samples: &HashSet<String>,
) -> (VcfHeader, Option<String>) {
    let mut header_line: Option<String> = None;
    loop {
        let mut line = String::new();
        let n = reader.read_line(&mut line).unwrap_or_else(|e| {
            eprintln!("ERROR: I/O error reading {}: {}", src, e);
            std::process::exit(1)
        });
        if n == 0 {
            break;
        }
        while line.ends_with('\n') || line.ends_with('\r') {
            line.pop();
        }
        if line.starts_with("##") {
            continue;
        } else if line.starts_with('#') {
            header_line = Some(line);
            break;
        } else if line.is_empty() {
            continue;
        } else {
            eprintln!(
                "ERROR: missing #CHROM header line in {} (found data line first)",
                src
            );
            std::process::exit(1);
        }
    }
    let header_line = header_line.unwrap_or_else(|| {
        eprintln!("ERROR: missing #CHROM header line in {}", src);
        std::process::exit(1)
    });
    let fields: Vec<&str> = header_line.split('\t').collect();
    if fields.len() < SAMPLE_OFFSET {
        eprintln!(
            "ERROR: VCF header line has fewer than 9 fields in {}",
            src
        );
        std::process::exit(1);
    }
    let all_ids: Vec<String> = fields[SAMPLE_OFFSET..]
        .iter()
        .map(|s| s.to_string())
        .collect();

    // first data line determines haploid/diploid status per sample
    let mut first_data_line: Option<String> = None;
    loop {
        let mut line = String::new();
        let n = reader.read_line(&mut line).unwrap_or(0);
        if n == 0 {
            break;
        }
        while line.ends_with('\n') || line.ends_with('\r') {
            line.pop();
        }
        if !line.is_empty() {
            first_data_line = Some(line);
            break;
        }
    }
    let is_diploid_all: Vec<bool> = match &first_data_line {
        Some(line) => is_diploid_per_sample(line, all_ids.len(), src),
        None => vec![true; all_ids.len()],
    };

    let mut unfiltered_index = Vec::new();
    let mut ids = Vec::new();
    let mut is_diploid = Vec::new();
    for (j, id) in all_ids.iter().enumerate() {
        if !exclude_samples.contains(id) {
            unfiltered_index.push(j);
            ids.push(id.clone());
            is_diploid.push(is_diploid_all[j]);
        }
    }
    if ids.is_empty() {
        eprintln!("ERROR: all samples in {} are excluded", src);
        std::process::exit(1);
    }
    (
        VcfHeader {
            src: src.to_string(),
            samples: Samples {
                ids: Arc::new(ids),
                is_diploid: Arc::new(is_diploid),
            },
            unfiltered_index,
            n_unfiltered_samples: all_ids.len(),
        },
        first_data_line,
    )
}

/// Port of `VcfHeader.isDiploid(String)`: a sample field is diploid when it
/// contains an allele separator ('|' or '/') anywhere in the field.
fn is_diploid_per_sample(rec: &str, n_samples: usize, src: &str) -> Vec<bool> {
    let start = match nth_tab_pos(rec, SAMPLE_OFFSET) {
        Some(p) => p + 1,
        None => {
            eprintln!("ERROR: VCF record format error in {}: {}", src, rec);
            std::process::exit(1)
        }
    };
    let mut list = Vec::with_capacity(n_samples);
    let mut no_allele_sep = true;
    for &b in rec.as_bytes()[start..].iter() {
        if b == b'\t' {
            list.push(!no_allele_sep);
            no_allele_sep = true;
        } else if b == b'/' || b == b'|' {
            no_allele_sep = false;
        }
    }
    list.push(!no_allele_sep);
    if list.len() != n_samples {
        eprintln!(
            "ERROR: VCF header line has {} sample fields, but a data line has {} in {}",
            n_samples,
            list.len(),
            src
        );
        std::process::exit(1);
    }
    list
}

/// Byte offset of the n-th tab (1-based count) in the record.
pub fn nth_tab_pos(rec: &str, n: usize) -> Option<usize> {
    let mut count = 0;
    for (i, &b) in rec.as_bytes().iter().enumerate() {
        if b == b'\t' {
            count += 1;
            if count == n {
                return Some(i);
            }
        }
    }
    None
}

/// A parsed target GT record (dense allele storage; -1 = missing).
pub struct GtRec {
    pub marker: Marker,
    /// 2 slots per sample; haploid samples have the allele duplicated
    pub alleles: Vec<i16>,
    /// true iff every genotype is phased and non-missing
    pub phased: bool,
}

/// Parses the GT field of every sample in the record
/// (`VcfRecGTParser.storeAlleles` semantics: if one allele of a diploid
/// genotype is missing, both are set to missing; a record is unphased if any
/// genotype is unphased or missing).
pub fn parse_gt_rec(header: &VcfHeader, rec: &str) -> Result<GtRec, String> {
    let (marker, _) = Marker::parse(rec)?;
    let n_alleles = marker.n_alleles as i32;
    let bytes = rec.as_bytes();
    let ninth_tab = nth_tab_pos(rec, 9).ok_or_else(|| {
        format!(
            "VCF record format error (fewer than 9 fields): {}",
            crate::marker::truncate(rec, 100)
        )
    })?;
    let n_samples = header.samples.len();
    let mut alleles: Vec<i16> = vec![0; n_samples << 1];
    let mut all_phased = true;

    let mut pos = ninth_tab; // index of tab preceding next sample field
    let mut unfilt: isize = -1;
    for s in 0..n_samples {
        let next_unfiltered = header.unfiltered_index[s] as isize;
        while unfilt + 1 < next_unfiltered {
            unfilt += 1;
            pos = match find_byte(bytes, pos + 1, b'\t') {
                Some(p) => p,
                None => return Err(field_count_err(header, rec)),
            };
        }
        unfilt += 1;
        if pos >= bytes.len() {
            return Err(field_count_err(header, rec));
        }
        let al_start = pos + 1;
        let al_end1 = allele_end1(bytes, al_start);
        if al_start == al_end1 {
            return Err(format!(
                "missing data for sample {} at marker [{}:{}]",
                header.samples.ids[s],
                marker.chrom(),
                marker.pos
            ));
        }
        let al_end2 = allele_end2(bytes, al_end1);
        let is_diploid = al_end1 != al_end2;
        if is_diploid != header.samples.is_diploid[s] {
            return Err(format!(
                "Sample {} has an inconsistent number of alleles at {}:{}",
                header.samples.ids[s],
                marker.chrom(),
                marker.pos
            ));
        }
        let mut a1 = parse_allele(rec, al_start, al_end1, n_alleles)?;
        let mut a2 = if al_end1 == al_end2 {
            a1
        } else {
            parse_allele(rec, al_end1 + 1, al_end2, n_alleles)?
        };
        if (a1 == -1) != (a2 == -1) {
            a1 = -1;
            a2 = -1;
        }
        let phased_gt = !is_diploid || bytes[al_end1] == b'|';
        if !phased_gt || a1 == -1 {
            all_phased = false;
        }
        let h1 = s << 1;
        alleles[h1] = a1 as i16;
        alleles[h1 | 1] = a2 as i16;
        pos = match find_byte(bytes, al_end2, b'\t') {
            Some(p) => p,
            None => bytes.len(),
        };
    }
    Ok(GtRec {
        marker,
        alleles,
        phased: all_phased,
    })
}

fn field_count_err(header: &VcfHeader, rec: &str) -> String {
    format!(
        "VCF data line has too few fields (source: {}): {}",
        header.src,
        crate::marker::truncate(rec, 100)
    )
}

#[inline]
fn find_byte(bytes: &[u8], start: usize, byte: u8) -> Option<usize> {
    bytes[start.min(bytes.len())..]
        .iter()
        .position(|&b| b == byte)
        .map(|p| p + start)
}

/// end (exclusive) of the first allele: stops at '/', '|', '\t', or ':'
#[inline]
fn allele_end1(bytes: &[u8], start: usize) -> usize {
    let mut i = start;
    while i < bytes.len() {
        match bytes[i] {
            b'/' | b'|' | b'\t' | b':' => return i,
            _ => i += 1,
        }
    }
    i
}

/// end (exclusive) of the second allele: stops at ':' or '\t'
#[inline]
fn allele_end2(bytes: &[u8], start: usize) -> usize {
    let mut i = start;
    while i < bytes.len() {
        match bytes[i] {
            b':' | b'\t' => return i,
            _ => i += 1,
        }
    }
    i
}

fn parse_allele(rec: &str, start: usize, end: usize, n_alleles: i32) -> Result<i32, String> {
    if start == end {
        return Err(format!(
            "Missing sample allele: {}",
            crate::marker::truncate(rec, 100)
        ));
    }
    let bytes = rec.as_bytes();
    let al: i32 = if start + 1 == end {
        let c = bytes[start];
        if c == b'.' {
            return Ok(-1);
        }
        (c as i32) - ('0' as i32)
    } else {
        rec[start..end]
            .parse()
            .map_err(|_| format!("Invalid allele [{}]", &rec[start..end]))?
    };
    if al < 0 || al >= n_alleles {
        return Err(format!(
            "Invalid allele [{}] in record \"{}...\"",
            &rec[start..end],
            crate::marker::truncate(rec, 60)
        ));
    }
    Ok(al)
}

/// Reads a marker/sample exclusion file (one identifier per line).
pub fn read_exclude_file(path: &Option<String>) -> HashSet<String> {
    let mut set = HashSet::new();
    if let Some(path) = path {
        let mut reader = open_text(path);
        let mut line = String::new();
        loop {
            line.clear();
            let n = reader.read_line(&mut line).unwrap_or(0);
            if n == 0 {
                break;
            }
            let t = line.trim();
            if !t.is_empty() {
                set.insert(t.to_string());
            }
        }
    }
    set
}

/// Iterator over the data lines of a VCF file.
/// Reads VCF data lines in batches into one reusable byte buffer, handing
/// out `&str` slices into it.  Reference-panel lines are large (one field
/// per sample, so ~25 KB at 5,000 samples), which made a `String` per line
/// costly: each allocation regrew from empty through a dozen reallocations.
/// A shared buffer amortizes that to one append per line.
pub struct LineSource {
    reader: Box<dyn BufRead + Send>,
    /// raw bytes of the current batch, newlines included
    buf: Vec<u8>,
    /// byte ranges of the batch's lines within `buf`, EOL stripped
    ranges: Vec<(usize, usize)>,
    /// index of the next line to hand out from `ranges` (single-line API)
    cursor: usize,
    pending: Option<String>,
    eof: bool,
}

impl LineSource {
    pub fn new(reader: Box<dyn BufRead + Send>, first_data_line: Option<String>) -> Self {
        LineSource {
            reader,
            buf: Vec::new(),
            ranges: Vec::new(),
            cursor: 0,
            pending: first_data_line,
            eof: false,
        }
    }

    /// Reads up to `max_lines` non-blank lines into the internal buffer,
    /// replacing any previous batch.  Returns the number of lines read;
    /// they are then accessible via `line(i)`.
    pub fn fill_batch(&mut self, max_lines: usize) -> usize {
        self.buf.clear();
        self.ranges.clear();
        self.cursor = 0;
        if let Some(first) = self.pending.take() {
            self.buf.extend_from_slice(first.as_bytes());
            self.ranges.push((0, self.buf.len()));
        }
        while self.ranges.len() < max_lines && !self.eof {
            let start = self.buf.len();
            let n = self.reader.read_until(b'\n', &mut self.buf).unwrap_or_else(|e| {
                eprintln!("ERROR: I/O error reading VCF file: {}", e);
                std::process::exit(1)
            });
            if n == 0 {
                self.eof = true;
                break;
            }
            let mut end = self.buf.len();
            while end > start && (self.buf[end - 1] == b'\n' || self.buf[end - 1] == b'\r') {
                end -= 1;
            }
            if end > start {
                self.ranges.push((start, end));
            } else {
                self.buf.truncate(start); // blank line
            }
        }
        self.ranges.len()
    }

    /// The `i`-th line of the current batch.
    #[inline]
    pub fn line(&self, i: usize) -> &str {
        let (s, e) = self.ranges[i];
        std::str::from_utf8(&self.buf[s..e]).unwrap_or_else(|_| {
            eprintln!("ERROR: VCF file contains invalid UTF-8");
            std::process::exit(1)
        })
    }

    /// The current batch as `&str` slices, for parallel parsing.
    pub fn batch(&self) -> Vec<&str> {
        (0..self.ranges.len()).map(|i| self.line(i)).collect()
    }

    /// Single-line access (used for target records, which are parsed one at
    /// a time); refills the batch buffer as needed.
    pub fn next_line(&mut self) -> Option<&str> {
        if self.cursor == self.ranges.len() {
            if self.fill_batch(LINE_BATCH) == 0 {
                return None;
            }
        }
        let i = self.cursor;
        self.cursor += 1;
        Some(self.line(i))
    }
}

/// Lines per read batch; keeps the shared buffer to a few MB even for
/// reference panels with thousands of samples.
pub const LINE_BATCH: usize = 512;
