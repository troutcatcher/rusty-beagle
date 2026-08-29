//! Reference-panel input: port of `vcf.RefIt` + `bref.SeqCoder3`.
//!
//! Records are parsed into an allele-coded form (major allele + carrier
//! lists).  High-frequency records are then "sequence coded" in blocks
//! exactly as Java's `SeqCoder3` does; the resulting per-block shared
//! hap→seq map and per-record seq→allele maps are semantically relevant:
//! block boundaries force marker-cluster splits in `ImpData`, and the
//! composed representation makes per-cluster haplotype coding fast.

use crate::marker::Marker;
use crate::vcfio::{nth_tab_pos, VcfHeader};
use std::sync::Arc;

/// Allele storage for one reference record.
pub enum RefAlleles {
    /// Major allele + sorted haplotype-index lists for each non-major allele.
    AlleleCoded {
        major: u8,
        /// `carriers[a]` is empty for the major allele
        carriers: Vec<Vec<u32>>,
    },
    /// Sequence-coded record: shared per-block hap→seq map plus a
    /// per-record seq→allele map.
    SeqCoded {
        block: u32,
        hap2seq: Arc<Vec<u16>>,
        seq2allele: Vec<u8>,
    },
}

pub struct RefRec {
    pub marker: Marker,
    pub alleles: RefAlleles,
    pub n_haps: usize,
}

impl RefRec {
    #[inline]
    pub fn is_allele_coded(&self) -> bool {
        matches!(self.alleles, RefAlleles::AlleleCoded { .. })
    }

    /// Block id for sequence-coded records (`map(0)` identity in Java).
    #[inline]
    pub fn block_id(&self) -> Option<u32> {
        match &self.alleles {
            RefAlleles::SeqCoded { block, .. } => Some(*block),
            _ => None,
        }
    }

    /// `RefGTRec.get(hap)`
    pub fn allele(&self, hap: usize) -> u8 {
        match &self.alleles {
            RefAlleles::AlleleCoded { major, carriers } => {
                for (a, list) in carriers.iter().enumerate() {
                    if !list.is_empty() && list.binary_search(&(hap as u32)).is_ok() {
                        return a as u8;
                    }
                }
                *major
            }
            RefAlleles::SeqCoded {
                hap2seq,
                seq2allele,
                ..
            } => seq2allele[hap2seq[hap] as usize],
        }
    }

    /// `RefGTRec.alleleCount(a)` for a non-major allele.
    pub fn allele_count(&self, allele: usize) -> usize {
        match &self.alleles {
            RefAlleles::AlleleCoded { carriers, .. } => carriers[allele].len(),
            RefAlleles::SeqCoded {
                hap2seq,
                seq2allele,
                ..
            } => hap2seq
                .iter()
                .filter(|&&s| seq2allele[s as usize] as usize == allele)
                .count(),
        }
    }

    /// Major allele: allele with the maximal count (ties → smallest index).
    pub fn major_allele(&self) -> u8 {
        match &self.alleles {
            RefAlleles::AlleleCoded { major, .. } => *major,
            RefAlleles::SeqCoded { .. } => {
                unreachable!("major_allele is only queried on allele-coded records")
            }
        }
    }

    /// Fills `dst` (length n_haps) with this record's alleles.
    pub fn fill_alleles(&self, dst: &mut [u8]) {
        match &self.alleles {
            RefAlleles::AlleleCoded { major, carriers } => {
                dst.fill(*major);
                for (a, list) in carriers.iter().enumerate() {
                    for &h in list {
                        dst[h as usize] = a as u8;
                    }
                }
            }
            RefAlleles::SeqCoded {
                hap2seq,
                seq2allele,
                ..
            } => {
                for (d, &s) in dst.iter_mut().zip(hap2seq.iter()) {
                    *d = seq2allele[s as usize];
                }
            }
        }
    }
}

/// Parses a reference VCF record: all genotypes must be phased and
/// non-missing (`VcfRecGTParser.phasedAlleles`).
pub fn parse_ref_rec(header: &VcfHeader, rec: &str) -> Result<(Marker, Vec<Vec<u32>>), String> {
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
    let mut hap_lists: Vec<Vec<u32>> = vec![Vec::new(); n_alleles as usize];
    let mut pos = ninth_tab;
    let mut unfilt: isize = -1;
    for s in 0..n_samples {
        let next_unfiltered = header.unfiltered_index[s] as isize;
        while unfilt + 1 < next_unfiltered {
            unfilt += 1;
            pos = match memchr_from(bytes, pos + 1, b'\t') {
                Some(p) => p,
                None => return Err(format!("VCF data line has too few fields: {}", header.src)),
            };
        }
        unfilt += 1;
        if pos >= bytes.len() {
            return Err(format!("VCF data line has too few fields: {}", header.src));
        }
        let al_start = pos + 1;
        let mut al_end1 = al_start;
        while al_end1 < bytes.len() {
            match bytes[al_end1] {
                b'/' | b'|' | b'\t' | b':' => break,
                _ => al_end1 += 1,
            }
        }
        if al_start == al_end1 {
            return Err(format!(
                "missing data for reference sample {} at marker [{}:{}]",
                header.samples.ids[s],
                marker.chrom(),
                marker.pos
            ));
        }
        let mut al_end2 = al_end1;
        while al_end2 < bytes.len() {
            match bytes[al_end2] {
                b':' | b'\t' => break,
                _ => al_end2 += 1,
            }
        }
        let is_diploid = al_end1 != al_end2;
        if is_diploid != header.samples.is_diploid[s] {
            return Err(format!(
                "Reference sample {} has an inconsistent number of alleles at {}:{}",
                header.samples.ids[s],
                marker.chrom(),
                marker.pos
            ));
        }
        let a1 = parse_allele_strict(rec, al_start, al_end1, n_alleles, &marker)?;
        let a2 = if al_end1 == al_end2 {
            a1
        } else {
            parse_allele_strict(rec, al_end1 + 1, al_end2, n_alleles, &marker)?
        };
        if is_diploid && bytes[al_end1] != b'|' {
            return Err(format!(
                "unphased genotype for reference sample {} at marker [{}:{}]",
                header.samples.ids[s],
                marker.chrom(),
                marker.pos
            ));
        }
        let h1 = (s << 1) as u32;
        hap_lists[a1 as usize].push(h1);
        hap_lists[a2 as usize].push(h1 | 1);
        pos = match memchr_from(bytes, al_end2, b'\t') {
            Some(p) => p,
            None => bytes.len(),
        };
    }
    Ok((marker, hap_lists))
}

#[inline]
fn memchr_from(bytes: &[u8], start: usize, byte: u8) -> Option<usize> {
    bytes[start.min(bytes.len())..]
        .iter()
        .position(|&b| b == byte)
        .map(|p| p + start)
}

fn parse_allele_strict(
    rec: &str,
    start: usize,
    end: usize,
    n_alleles: i32,
    marker: &Marker,
) -> Result<i32, String> {
    let bytes = rec.as_bytes();
    if start == end {
        return Err("Missing sample allele".to_string());
    }
    let al: i32 = if start + 1 == end {
        let c = bytes[start];
        if c == b'.' {
            return Err(format!(
                "missing allele in reference VCF at marker [{}:{}]",
                marker.chrom(),
                marker.pos
            ));
        }
        (c as i32) - ('0' as i32)
    } else {
        rec[start..end]
            .parse()
            .map_err(|_| format!("Invalid allele [{}]", &rec[start..end]))?
    };
    if al < 0 || al >= n_alleles {
        return Err(format!(
            "Invalid allele [{}] at marker [{}:{}]",
            &rec[start..end],
            marker.chrom(),
            marker.pos
        ));
    }
    Ok(al)
}

/// hap-index lists → allele-coded record parts (major allele = maximal count,
/// ties broken by smallest allele index; the major carrier list is dropped).
pub fn to_allele_coded(mut hap_lists: Vec<Vec<u32>>) -> (u8, Vec<Vec<u32>>) {
    let mut major = 0usize;
    for a in 1..hap_lists.len() {
        if hap_lists[a].len() > hap_lists[major].len() {
            major = a;
        }
    }
    hap_lists[major] = Vec::new();
    (major as u8, hap_lists)
}

/// Port of `bref.SeqCoder3` restricted to what imputation needs:
/// the partition-refinement state, block flush boundaries, and the final
/// per-block hap→seq / per-record seq→allele maps.
pub struct SeqCoder {
    n_haps: usize,
    max_n_seq: usize,
    /// pending allele-coded records to be sequence-coded in this block
    recs: Vec<(Marker, u8, Vec<Vec<u32>>)>,
    hap2seq: Vec<u32>,
    seq2cnt: Vec<i32>,
    seq2allele_seq_map: Vec<Vec<i32>>,
}

pub const MAX_NALLELES: usize = 255;
pub const COMPRESS_FREQ_THRESHOLD: f64 = 0.995;

/// `SeqCoder3.defaultMaxNSeq`
pub fn default_max_n_seq(n_samples: usize) -> usize {
    assert!(n_samples >= 1);
    if n_samples == 1 {
        3
    } else {
        let exponent = 2.0 * (n_samples as f64).log10() + 1.0;
        let max_n_seq = (2.0f64).powf(exponent).floor();
        if max_n_seq > 65535.0 {
            65535
        } else {
            max_n_seq as usize
        }
    }
}

impl SeqCoder {
    pub fn new(n_samples: usize) -> SeqCoder {
        let n_haps = n_samples << 1;
        let max_n_seq = default_max_n_seq(n_samples);
        let mut coder = SeqCoder {
            n_haps,
            max_n_seq,
            recs: Vec::with_capacity(128),
            hap2seq: vec![0; n_haps],
            seq2cnt: Vec::new(),
            seq2allele_seq_map: Vec::new(),
        };
        coder.initialize();
        coder
    }

    pub fn max_n_seq(&self) -> usize {
        self.max_n_seq
    }

    fn initialize(&mut self) {
        self.recs.clear();
        self.seq2cnt.clear();
        self.seq2allele_seq_map.clear();
        self.hap2seq.fill(0);
        self.seq2cnt.push(self.n_haps as i32);
        self.seq2allele_seq_map.push(Vec::new());
    }

    /// `SeqCoder3.add`; returns false when the record does not fit in the
    /// current block (caller must flush and retry).
    pub fn add(&mut self, marker: Marker, major: u8, carriers: Vec<Vec<u32>>) -> bool {
        let success = self.set_allele_map(major, &carriers);
        if success {
            let n_alleles = carriers.len();
            for a in 0..n_alleles {
                if a != major as usize {
                    for &h in &carriers[a] {
                        let h = h as usize;
                        let old_seq = self.hap2seq[h] as usize;
                        let list = &self.seq2allele_seq_map[old_seq];
                        let mut index = 0;
                        while index < list.len() && list[index] != a as i32 {
                            index += 2;
                        }
                        let new_seq = list[index + 1] as usize;
                        if new_seq != old_seq {
                            while new_seq >= self.seq2cnt.len() {
                                self.seq2cnt.push(0);
                            }
                            self.hap2seq[h] = new_seq as u32;
                            self.seq2cnt[old_seq] -= 1;
                            self.seq2cnt[new_seq] += 1;
                        }
                    }
                }
            }
            self.recs.push((marker, major, carriers));
        }
        debug_assert_eq!(self.seq2cnt.len(), self.seq2allele_seq_map.len());
        success
    }

    fn set_allele_map(&mut self, major: u8, carriers: &[Vec<u32>]) -> bool {
        let n_start_seq = self.seq2cnt.len();
        let mut seq2non_major_cnt = vec![0i32; n_start_seq];
        for list in self.seq2allele_seq_map.iter_mut() {
            list.clear();
        }
        let n_alleles = carriers.len();
        for a in 0..n_alleles {
            if a != major as usize {
                for &h in &carriers[a] {
                    let seq = self.hap2seq[h as usize] as usize;
                    seq2non_major_cnt[seq] += 1;
                    let list_len = self.seq2allele_seq_map[seq].len();
                    if list_len == 0 {
                        let list = &mut self.seq2allele_seq_map[seq];
                        list.push(a as i32);
                        list.push(seq as i32);
                    } else {
                        let mut index = 0;
                        while index < list_len && self.seq2allele_seq_map[seq][index] != a as i32 {
                            index += 2;
                        }
                        if index == list_len {
                            let next = self.seq2allele_seq_map.len() as i32;
                            let list = &mut self.seq2allele_seq_map[seq];
                            list.push(a as i32);
                            list.push(next);
                            self.seq2allele_seq_map.push(Vec::with_capacity(4));
                        }
                    }
                }
            }
        }
        self.add_major_allele(&seq2non_major_cnt, major);
        if self.seq2allele_seq_map.len() > self.max_n_seq {
            self.seq2allele_seq_map.truncate(n_start_seq);
            false
        } else {
            true
        }
    }

    fn add_major_allele(&mut self, seq2non_major_cnt: &[i32], major: u8) {
        for seq in 0..seq2non_major_cnt.len() {
            if seq2non_major_cnt[seq] < self.seq2cnt[seq] {
                let list_len = self.seq2allele_seq_map[seq].len();
                if list_len == 0 {
                    let list = &mut self.seq2allele_seq_map[seq];
                    list.push(major as i32);
                    list.push(seq as i32);
                } else {
                    // assign the major allele the existing sequence index
                    let next = self.seq2allele_seq_map.len() as i32;
                    let list = &mut self.seq2allele_seq_map[seq];
                    let first_allele = list[0];
                    debug_assert_eq!(list[1], seq as i32);
                    list.push(first_allele);
                    list.push(next);
                    list[0] = major as i32;
                    self.seq2allele_seq_map.push(Vec::with_capacity(4));
                }
            }
        }
    }

    /// `SeqCoder3.getCompressedList`: finalizes the current block and
    /// returns its sequence-coded records; resets the coder.
    pub fn flush(&mut self, block_id: u32, n_haps: usize) -> Vec<RefRec> {
        if self.recs.is_empty() {
            self.initialize();
            return Vec::new();
        }
        let n_seq = self.seq2allele_seq_map.len();
        // seq -> first haplotype carrying it
        let mut seq2hap: Vec<i32> = vec![-1; n_seq];
        for (h, &seq) in self.hap2seq.iter().enumerate() {
            if seq2hap[seq as usize] == -1 {
                seq2hap[seq as usize] = h as i32;
            }
        }
        let hap2seq: Arc<Vec<u16>> =
            Arc::new(self.hap2seq.iter().map(|&s| s as u16).collect());
        let mut out = Vec::with_capacity(self.recs.len());
        for (marker, major, carriers) in self.recs.drain(..) {
            let mut seq2allele = vec![major; n_seq];
            for (a, list) in carriers.iter().enumerate() {
                if a != major as usize {
                    for &h in list {
                        // note: only sequences whose first hap carries allele a
                        // matter here; we translate rec.get(seq2Hap[s]) directly
                        let _ = h;
                    }
                }
            }
            // seq2allele[s] = rec.get(seq2hap[s])
            for s in 0..n_seq {
                let h = seq2hap[s];
                if h >= 0 {
                    let mut allele = major;
                    for (a, list) in carriers.iter().enumerate() {
                        if a != major as usize && list.binary_search(&(h as u32)).is_ok() {
                            allele = a as u8;
                            break;
                        }
                    }
                    seq2allele[s] = allele;
                }
            }
            out.push(RefRec {
                marker,
                alleles: RefAlleles::SeqCoded {
                    block: block_id,
                    hap2seq: hap2seq.clone(),
                    seq2allele,
                },
                n_haps,
            });
        }
        self.initialize();
        out
    }
}

/// Port of `vcf.RefIt`: streams reference records, sequence-coding
/// high-frequency records in blocks.
pub struct RefReader {
    pub header: VcfHeader,
    lines: crate::vcfio::LineSource,
    coder: SeqCoder,
    max_seq_coded_alleles: usize,
    max_seq_coding_major_cnt: i64,
    /// finalized records ready to be consumed
    out: std::collections::VecDeque<RefRec>,
    /// pending records of the current block (None = slot for seq-coded record)
    pending: Vec<Option<RefRec>>,
    /// parallel with block flush: parts for seq-coded slots
    n_haps: usize,
    last_chrom: i32,
    next_block_id: u32,
    exclude_markers: std::collections::HashSet<String>,
    done: bool,
}

impl RefReader {
    pub fn new(
        header: VcfHeader,
        lines: crate::vcfio::LineSource,
        exclude_markers: std::collections::HashSet<String>,
    ) -> RefReader {
        let n_samples = header.samples.len();
        let n_haps = n_samples << 1;
        let coder = SeqCoder::new(n_samples);
        let max_seq_coded_alleles = coder.max_n_seq().min(MAX_NALLELES);
        let max_seq_coding_major_cnt =
            ((n_haps as f64) * COMPRESS_FREQ_THRESHOLD - 1.0).floor() as i64;
        RefReader {
            header,
            lines,
            coder,
            max_seq_coded_alleles,
            max_seq_coding_major_cnt,
            out: std::collections::VecDeque::new(),
            pending: Vec::new(),
            n_haps,
            last_chrom: -1,
            next_block_id: 0,
            exclude_markers,
            done: false,
        }
    }

    pub fn n_haps(&self) -> usize {
        self.n_haps
    }

    fn apply_seq_coding(&self, marker: &Marker, major: u8, carriers: &[Vec<u32>]) -> bool {
        if marker.n_alleles as usize > self.max_seq_coded_alleles {
            return false;
        }
        let mut maj_cnt = self.n_haps as i64;
        for (a, list) in carriers.iter().enumerate() {
            if a != major as usize {
                maj_cnt -= list.len() as i64;
            }
        }
        maj_cnt <= self.max_seq_coding_major_cnt
    }

    fn flush(&mut self) {
        let block_id = self.next_block_id;
        self.next_block_id += 1;
        let mut coded = self.coder.flush(block_id, self.n_haps).into_iter();
        for slot in self.pending.drain(..) {
            match slot {
                Some(rec) => self.out.push_back(rec),
                None => self.out.push_back(coded.next().expect("seq-coded slot")),
            }
        }
        debug_assert!(coded.next().is_none());
    }

    /// Parses lines until at least one record is finalized or input ends.
    fn fill(&mut self) {
        while self.out.is_empty() && !self.done {
            match self.lines.next_line() {
                None => {
                    self.done = true;
                    self.flush();
                }
                Some(line) => {
                    let (marker, hap_lists) = parse_ref_rec(&self.header, &line)
                        .unwrap_or_else(|e| {
                            eprintln!("ERROR: {}", e);
                            std::process::exit(1)
                        });
                    if let Some(id) = &marker.id {
                        if self.exclude_markers.contains(id.as_ref()) {
                            continue;
                        }
                    }
                    let (major, carriers) = to_allele_coded(hap_lists);
                    let chrom = marker.chrom_idx as i32;
                    if self.last_chrom == -1 {
                        self.last_chrom = chrom;
                    }
                    if chrom != self.last_chrom {
                        self.flush();
                        self.last_chrom = chrom;
                    }
                    if !self.apply_seq_coding(&marker, major, &carriers) {
                        self.pending.push(Some(RefRec {
                            marker,
                            alleles: RefAlleles::AlleleCoded { major, carriers },
                            n_haps: self.n_haps,
                        }));
                    } else {
                        let ok = self.coder.add(marker.clone(), major, carriers.clone());
                        if !ok {
                            self.flush();
                            let ok2 = self.coder.add(marker, major, carriers);
                            assert!(ok2, "SeqCoder add failed after flush");
                        }
                        self.pending.push(None);
                    }
                }
            }
        }
    }

    pub fn next(&mut self) -> Option<Arc<RefRec>> {
        if self.out.is_empty() {
            self.fill();
        }
        self.out.pop_front().map(Arc::new)
    }
}
