//! Port of `imp.RefHapHash`, `imp.ImputedVcfWriter`, `imp.ImputedRecBuilder`,
//! `vcf.VcfWriter` meta lines, and `main.WindowWriter`.

use crate::bgzf;
use crate::hmm::StateProbs;
use crate::impdata::ImpData;
use crate::javautil::{java_rint, JavaRandom};
use crate::refpanel::{RefAlleles, RefRec};
use crate::vcfio::Samples;
use crate::windows::Window;
use rayon::prelude::*;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

/// Formatting tables shared by all output records.
pub struct Fmt {
    pub ds_vals: Vec<String>,
    pub r2_vals: Vec<String>,
}

impl Fmt {
    pub fn new() -> Fmt {
        Fmt {
            ds_vals: crate::javautil::ds_vals(),
            r2_vals: crate::javautil::r2_vals(),
        }
    }
}

/// Per-thread buffers reused across the marker clusters of a window.
#[derive(Default)]
struct OutScratch {
    hap_to_index: Vec<u32>,
    seen: Vec<u64>,
}

/// Port of `imp.RefHapHash`.
///
/// The per-haplotype ALT-allele lists are held in one flat CSR buffer rather
/// than a `Vec<Vec<_>>`: this is rebuilt for every marker cluster, and a
/// separate growing `Vec` per haplotype cost thousands of allocations per
/// cluster.
struct RefHapHash<'a> {
    ref_recs: &'a [Arc<RefRec>],
    i2hap: Vec<i32>,
    i2hash: Vec<i32>,
    /// `alt_data[alt_start[i]..alt_start[i + 1]]` = (marker offset, ALT
    /// allele) for haplotype `i`
    alt_start: Vec<u32>,
    alt_data: Vec<(u32, u32)>,
    /// reference haplotype -> its index in `i2hap`; every state of every
    /// target haplotype is looked up here, which made a binary search over
    /// `i2hap` (a dozen scattered probes) the cost of the accumulation loop.
    /// Borrowed per thread and never cleared between clusters: only the
    /// entries this cluster just wrote are ever read back.
    hap_to_index: &'a mut Vec<u32>,
    start: usize,
    end: usize,
}

impl<'a> RefHapHash<'a> {
    fn new(
        state_probs: &[StateProbs],
        targ_cluster: usize,
        ref_recs: &'a [Arc<RefRec>],
        n_ref_haps: usize,
        hap_to_index: &'a mut Vec<u32>,
        seen: &mut Vec<u64>,
        start: usize,
        end: usize,
    ) -> RefHapHash<'a> {
        assert!(start < end);
        // Every state of every target haplotype names a reference haplotype,
        // but only a small fraction are distinct (millions of states, tens of
        // thousands of haplotypes at 200k reference samples). Mark first
        // sightings in a bitset so only the distinct haplotypes are collected
        // and sorted, rather than sorting the full multiset.
        let n_words = (n_ref_haps + 63) >> 6;
        if seen.len() < n_words {
            seen.resize(n_words, 0);
        }
        let mut list: Vec<i32> = Vec::new();
        for sp in state_probs {
            for st in sp.states(targ_cluster) {
                let h = st.hap as usize;
                let bit = 1u64 << (h & 63);
                let word = &mut seen[h >> 6];
                if *word & bit == 0 {
                    *word |= bit;
                    list.push(st.hap);
                }
            }
        }
        list.sort_unstable();
        // clearing whole words is safe: every bit set above belongs to a
        // haplotype in `list`, so this leaves the bitset all-zero again
        for &h in &list {
            seen[h as usize >> 6] = 0;
        }
        let i2hap = list;
        let n = i2hap.len();
        let _ = targ_cluster;
        if hap_to_index.len() < n_ref_haps {
            hap_to_index.resize(n_ref_haps, 0);
        }
        for (i, &h) in i2hap.iter().enumerate() {
            hap_to_index[h as usize] = i as u32;
        }
        let mut hash = RefHapHash {
            ref_recs,
            i2hash: vec![0; n],
            alt_start: Vec::new(),
            alt_data: Vec::new(),
            hap_to_index,
            i2hap,
            start,
            end,
        };
        hash.set_hash_and_alt_alleles();
        hash
    }

    fn set_hash_and_alt_alleles(&mut self) {
        let mut rand = JavaRandom::new(self.start as i64);
        let n = self.i2hap.len();
        // (haplotype index, marker offset, ALT allele), later bucketed by
        // haplotype into the CSR arrays
        let mut staged: Vec<(u32, u32, u32)> = Vec::with_capacity(4 * n);
        let mut counts = vec![0u32; n + 1];
        let mut allele_hash: Vec<i32> = Vec::new();
        // block sequence of each `i2hap` entry, valid while `cached_block` holds
        let mut seq_of: Vec<u16> = Vec::with_capacity(n);
        let mut cached_block: Option<u32> = None;
        for m in self.start..self.end {
            let rec = &self.ref_recs[m];
            let marker_offset = (m - self.start) as u32;
            match &rec.alleles {
                RefAlleles::AlleleCoded { major, carriers } if *major == 0 => {
                    // low-ALT-frequency update: walk the (short) carrier lists
                    let n_alleles = carriers.len();
                    for al in 1..n_alleles {
                        let hash = rand.next_int();
                        let n_copies = carriers[al].len();
                        if n < n_copies {
                            for (i, &h) in self.i2hap.iter().enumerate() {
                                if carriers[al].binary_search(&(h as u32)).is_ok() {
                                    self.i2hash[i] = self.i2hash[i].wrapping_add(hash);
                                    staged.push((i as u32, marker_offset, al as u32));
                                    counts[i + 1] += 1;
                                }
                            }
                        } else {
                            for &hap in &carriers[al] {
                                if let Ok(i) = self.i2hap.binary_search(&(hap as i32)) {
                                    self.i2hash[i] = self.i2hash[i].wrapping_add(hash);
                                    staged.push((i as u32, marker_offset, al as u32));
                                    counts[i + 1] += 1;
                                }
                            }
                        }
                    }
                }
                RefAlleles::SeqCoded {
                    block,
                    hap2seq,
                    seq2allele,
                } => {
                    // hap2seq is shared by every marker of a sequence-coded
                    // block, so resolve each haplotype's sequence once per
                    // block instead of once per marker: the per-marker lookup
                    // then hits the small seq2allele table rather than the
                    // panel-sized hap2seq.
                    if cached_block != Some(*block) {
                        seq_of.clear();
                        seq_of.extend(self.i2hap.iter().map(|&h| hap2seq[h as usize]));
                        cached_block = Some(*block);
                    }
                    let n_alleles = rec.marker.n_alleles as usize;
                    allele_hash.clear();
                    allele_hash.push(0);
                    for _ in 1..n_alleles {
                        allele_hash.push(rand.next_int());
                    }
                    for (i, &s) in seq_of.iter().enumerate() {
                        let allele = seq2allele[s as usize] as usize;
                        if allele != 0 {
                            self.i2hash[i] = self.i2hash[i].wrapping_add(allele_hash[allele]);
                            staged.push((i as u32, marker_offset, allele as u32));
                            counts[i + 1] += 1;
                        }
                    }
                }
                _ => {
                    let n_alleles = rec.marker.n_alleles as usize;
                    allele_hash.clear();
                    allele_hash.push(0);
                    for _ in 1..n_alleles {
                        allele_hash.push(rand.next_int());
                    }
                    for (i, &h) in self.i2hap.iter().enumerate() {
                        let allele = rec.allele(h as usize) as usize;
                        if allele != 0 {
                            self.i2hash[i] = self.i2hash[i].wrapping_add(allele_hash[allele]);
                            staged.push((i as u32, marker_offset, allele as u32));
                            counts[i + 1] += 1;
                        }
                    }
                }
            }
        }
        // prefix-sum the per-haplotype counts, then scatter into CSR order
        for i in 0..n {
            counts[i + 1] += counts[i];
        }
        let total = counts[n] as usize;
        let mut alt_data = vec![(0u32, 0u32); total];
        let mut cursor = counts.clone();
        for &(i, off, al) in &staged {
            let slot = &mut cursor[i as usize];
            alt_data[*slot as usize] = (off, al);
            *slot += 1;
        }
        self.alt_start = counts;
        self.alt_data = alt_data;
    }

    #[inline]
    fn hap2index(&self, hap: i32) -> usize {
        self.hap_to_index[hap as usize] as usize
    }

    #[inline]
    fn hash(&self, index: usize) -> i32 {
        self.i2hash[index]
    }

    /// `RefHapHash.setAlleles`
    fn set_alleles(&self, index: usize, alleles: &mut [i32]) {
        alleles[..self.end - self.start].fill(0);
        let lo = self.alt_start[index] as usize;
        let hi = self.alt_start[index + 1] as usize;
        for &(off, al) in &self.alt_data[lo..hi] {
            alleles[off as usize] = al as i32;
        }
    }
}

/// Port of `imp.ImputedRecBuilder`.
struct ImputedRecBuilder<'a> {
    marker: &'a crate::marker::Marker,
    n_alleles: usize,
    n_input_targ_haps: usize,
    ap: bool,
    gp: bool,
    sum_al_probs: Vec<f32>,
    sum_al_probs2: Vec<f32>,
    sample_data: String,
    hap_cnt: usize,
    fmt: &'a Fmt,
    hom_ref_field: &'a str, // "" when not applicable
}

impl<'a> ImputedRecBuilder<'a> {
    fn new(
        marker: &'a crate::marker::Marker,
        n_input_targ_haps: usize,
        ap: bool,
        gp: bool,
        fmt: &'a Fmt,
        hom_ref_fields: &'a [String],
    ) -> ImputedRecBuilder<'a> {
        let n_alleles = marker.n_alleles as usize;
        let hom_ref_field = if n_alleles < hom_ref_fields.len() {
            &hom_ref_fields[n_alleles]
        } else {
            ""
        };
        ImputedRecBuilder {
            marker,
            n_alleles,
            n_input_targ_haps,
            ap,
            gp,
            sum_al_probs: vec![0.0; n_alleles],
            sum_al_probs2: vec![0.0; n_alleles],
            sample_data: String::with_capacity(200 + n_input_targ_haps * 5),
            hap_cnt: 0,
            fmt,
            hom_ref_field,
        }
    }

    /// diploid sample
    fn add_sample_data2(&mut self, a1: &mut [f32], a2: &mut [f32]) {
        self.hap_cnt += 2;
        if a1[0] == 1.0 && a2[0] == 1.0 && !self.hom_ref_field.is_empty() {
            self.sample_data.push_str(self.hom_ref_field);
        } else {
            scale(&mut a1[..self.n_alleles]);
            scale(&mut a2[..self.n_alleles]);
            self.sample_data.push('\t');
            let m1 = max_index(&a1[..self.n_alleles]);
            let m2 = max_index(&a2[..self.n_alleles]);
            push_usize(&mut self.sample_data, m1);
            self.sample_data.push('|');
            push_usize(&mut self.sample_data, m2);
            for a in 1..self.n_alleles {
                let dose = a1[a] + a2[a];
                let dose2 = a1[a] * a1[a] + a2[a] * a2[a];
                self.sum_al_probs[a] += dose;
                self.sum_al_probs2[a] += dose2;
                self.sample_data.push(if a == 1 { ':' } else { ',' });
                self.sample_data
                    .push_str(&self.fmt.ds_vals[ds_index(dose)]);
            }
            if self.ap {
                for a in 1..self.n_alleles {
                    self.sample_data.push(if a == 1 { ':' } else { ',' });
                    self.sample_data
                        .push_str(&self.fmt.ds_vals[ds_index(a1[a])]);
                }
                for a in 1..self.n_alleles {
                    self.sample_data.push(if a == 1 { ':' } else { ',' });
                    self.sample_data
                        .push_str(&self.fmt.ds_vals[ds_index(a2[a])]);
                }
            }
            if self.gp {
                for i2 in 0..self.n_alleles {
                    for i1 in 0..=i2 {
                        let mut prob = a1[i1] * a2[i2];
                        if i1 != i2 {
                            prob += a1[i2] * a2[i1];
                        }
                        self.sample_data.push(if i2 == 0 { ':' } else { ',' });
                        self.sample_data
                            .push_str(&self.fmt.ds_vals[ds_index(prob)]);
                    }
                }
            }
        }
    }

    /// haploid sample
    fn add_sample_data1(&mut self, a1: &mut [f32]) {
        self.hap_cnt += 1;
        scale(&mut a1[..self.n_alleles]);
        self.sample_data.push('\t');
        let m1 = max_index(&a1[..self.n_alleles]);
        push_usize(&mut self.sample_data, m1);
        for a in 1..self.n_alleles {
            let dose = a1[a];
            let dose2 = a1[a] * a1[a];
            self.sum_al_probs[a] += dose;
            self.sum_al_probs2[a] += dose2;
            self.sample_data.push(if a == 1 { ':' } else { ',' });
            self.sample_data
                .push_str(&self.fmt.ds_vals[ds_index(dose)]);
        }
        if self.ap {
            for a in 1..self.n_alleles {
                self.sample_data.push(if a == 1 { ':' } else { ',' });
                self.sample_data
                    .push_str(&self.fmt.ds_vals[ds_index(a1[a])]);
            }
        }
    }

    fn print_rec(&self, out: &mut String, is_imputed: bool) {
        assert_eq!(
            self.hap_cnt, self.n_input_targ_haps,
            "inconsistent data in ImputedRecBuilder"
        );
        // marker fields
        out.push_str(&self.marker.chrom());
        out.push('\t');
        let _ = write!(out, "{}", self.marker.pos);
        out.push('\t');
        out.push_str(self.marker.id_str());
        out.push('\t');
        out.push_str(&self.marker.alleles);
        out.push('\t');
        out.push('.'); // QUAL
        out.push('\t');
        out.push_str("PASS"); // FILTER
        out.push('\t');
        self.print_info_field(out, is_imputed); // INFO
        out.push('\t');
        out.push_str("GT:DS"); // FORMAT
        if self.ap {
            out.push_str(":AP1:AP2");
        }
        if self.gp {
            out.push_str(":GP");
        }
        out.push_str(&self.sample_data);
        out.push('\n');
    }

    fn print_info_field(&self, out: &mut String, is_imputed: bool) {
        if self.n_alleles == 1 {
            if is_imputed {
                out.push_str("IMP");
            }
        } else {
            for a in 1..self.n_alleles {
                out.push_str(if a == 1 { "DR2=" } else { "," });
                let idx = java_rint((100.0f32 * self.r2(a)) as f64) as usize;
                out.push_str(&self.fmt.r2_vals[idx]);
            }
            for a in 1..self.n_alleles {
                out.push_str(if a == 1 { ";AF=" } else { "," });
                let af = self.sum_al_probs[a] / self.n_input_targ_haps as f32;
                let _ = write!(out, "{:.4}", af as f64);
            }
            if let Some(end_subfield) = self.marker.end_subfield() {
                out.push(';');
                out.push_str(&end_subfield);
            }
            if is_imputed {
                out.push_str(";IMP");
            }
        }
    }

    fn r2(&self, allele: usize) -> f32 {
        let sum = self.sum_al_probs[allele];
        if sum == 0.0 {
            0.0
        } else {
            let sum2 = self.sum_al_probs2[allele];
            let mean_term = sum * sum / self.n_input_targ_haps as f32;
            let num = sum2 - mean_term;
            let den = sum - mean_term;
            if num <= 0.0 {
                0.0
            } else {
                num / den
            }
        }
    }
}

#[inline]
fn ds_index(dose: f32) -> usize {
    java_rint((100.0f32 * dose) as f64) as usize
}

#[inline]
fn scale(fa: &mut [f32]) {
    let mut sum = 0.0f32;
    for &f in fa.iter() {
        sum += f;
    }
    for f in fa.iter_mut() {
        *f /= sum;
    }
}

#[inline]
fn max_index(fa: &[f32]) -> usize {
    let mut max_index = 0;
    for j in 1..fa.len() {
        if fa[j] > fa[max_index] {
            max_index = j;
        }
    }
    max_index
}

fn push_usize(s: &mut String, v: usize) {
    let _ = write!(s, "{}", v);
}

/// `ImputedRecBuilder.defaultHomRefFields` /  `homRefFields(ap, gp)`:
/// index k = the shared sample field for a hom-ref sample at a marker with
/// k alleles (only defined for k < 5).
fn hom_ref_fields(ap: bool, gp: bool, fmt: &Fmt) -> Vec<String> {
    let mut sa = vec![String::new(); 5];
    sa[1] = "\t0|0".to_string();
    sa[2] = "\t0|0:0".to_string();
    for j in 3..5 {
        sa[j] = format!("{},0", sa[j - 1]);
    }
    if ap || gp {
        for n_al in 1..5usize {
            let mut sb = sa[n_al].clone();
            if ap {
                for a in 1..n_al {
                    sb.push(if a == 1 { ':' } else { ',' });
                    sb.push_str(&fmt.ds_vals[0]);
                }
                for a in 1..n_al {
                    sb.push(if a == 1 { ':' } else { ',' });
                    sb.push_str(&fmt.ds_vals[0]);
                }
            }
            if gp {
                sb.push(':');
                sb.push_str(&fmt.ds_vals[100]);
                for i2 in 1..n_al {
                    for _i1 in 0..=i2 {
                        sb.push(',');
                        sb.push_str(&fmt.ds_vals[0]);
                    }
                }
            }
            sa[n_al] = sb;
        }
    }
    sa
}

/// Port of `imp.ImputedVcfWriter`: writes the VCF records for one target
/// cluster into a `String`.
#[allow(clippy::too_many_arguments)]
fn append_records(
    imp_data: &ImpData,
    window: &Window,
    samples: &Samples,
    state_probs: &[StateProbs],
    win_ref_start: usize,
    win_ref_end: usize,
    targ_cluster: usize,
    fmt: &Fmt,
    hom_ref: &[String],
    scratch_out: &mut OutScratch,
) -> String {
    // bounds from the ImputedVcfWriter constructor
    let ref_start = if targ_cluster == 0 {
        win_ref_start
    } else {
        win_ref_start.max(imp_data.ref_cluster_start[targ_cluster])
    };
    let (clust_end, ref_end) = if targ_cluster < imp_data.n_clusters - 1 {
        let tmp_clust_end = win_ref_start.max(imp_data.ref_cluster_end[targ_cluster]);
        (
            tmp_clust_end.min(win_ref_end),
            imp_data.ref_cluster_start[targ_cluster + 1].min(win_ref_end),
        )
    } else {
        (win_ref_end, win_ref_end)
    };
    let mut out = String::new();
    if ref_start >= ref_end {
        return out;
    }
    let ref_recs = &window.ref_recs;
    let tq0 = std::time::Instant::now();
    let rhh = RefHapHash::new(
        state_probs,
        targ_cluster,
        ref_recs,
        imp_data.n_ref_haps,
        &mut scratch_out.hap_to_index,
        &mut scratch_out.seen,
        ref_start,
        ref_end,
    );
    out_add(2, tq0.elapsed().as_nanos() as u64);
    let n_markers = ref_end - ref_start;
    let mut rec_builders: Vec<ImputedRecBuilder> = (ref_start..ref_end)
        .map(|m| {
            ImputedRecBuilder::new(
                &ref_recs[m].marker,
                imp_data.n_input_targ_haps,
                imp_data.ap,
                imp_data.gp,
                fmt,
                hom_ref,
            )
        })
        .collect();
    // allele probabilities for the two haplotypes of the current sample,
    // flat with per-marker offsets: a `Vec` per marker meant an allocation
    // per marker per cluster, and a pointer chase in the accumulation loop
    let mut al_off: Vec<u32> = Vec::with_capacity(n_markers + 1);
    let mut total = 0u32;
    for m in ref_start..ref_end {
        al_off.push(total);
        total += ref_recs[m].marker.n_alleles as u32;
    }
    al_off.push(total);
    let mut a1_probs: Vec<f32> = vec![0.0; total as usize];
    let mut a2_probs: Vec<f32> = vec![0.0; total as usize];
    let is_imputed: Vec<bool> = (ref_start..ref_end)
        .map(|m| window.indices.marker_to_targ_marker[m] == -1)
        .collect();

    let mut scratch = SetAlProbsScratch::new(n_markers);
    let n = state_probs.len();
    let mut h = 0;
    while h < n {
        let is_diploid = samples.is_diploid[h >> 1];
        let tp0 = std::time::Instant::now();
        set_al_probs(
            imp_data,
            &state_probs[h],
            &rhh,
            targ_cluster,
            ref_start,
            clust_end,
            ref_end,
            &mut a1_probs,
            &al_off,
            &mut scratch,
        );
        set_al_probs(
            imp_data,
            &state_probs[h + 1],
            &rhh,
            targ_cluster,
            ref_start,
            clust_end,
            ref_end,
            &mut a2_probs,
            &al_off,
            &mut scratch,
        );
        out_add(3, tp0.elapsed().as_nanos() as u64);
        let tb0 = std::time::Instant::now();
        if is_diploid {
            for m in 0..n_markers {
                if !is_imputed[m] {
                    set_to_obs_alleles(
                        imp_data, window, ref_start, m, h, &mut a1_probs, &mut a2_probs, &al_off,
                    );
                }
                let (lo, hi) = (al_off[m] as usize, al_off[m + 1] as usize);
                rec_builders[m].add_sample_data2(&mut a1_probs[lo..hi], &mut a2_probs[lo..hi]);
                a1_probs[lo..hi].fill(0.0);
                a2_probs[lo..hi].fill(0.0);
            }
        } else {
            for m in 0..n_markers {
                if !is_imputed[m] {
                    set_to_obs_alleles(
                        imp_data, window, ref_start, m, h, &mut a1_probs, &mut a2_probs, &al_off,
                    );
                }
                let (lo, hi) = (al_off[m] as usize, al_off[m + 1] as usize);
                rec_builders[m].add_sample_data1(&mut a1_probs[lo..hi]);
                a1_probs[lo..hi].fill(0.0);
                a2_probs[lo..hi].fill(0.0);
            }
        }
        out_add(4, tb0.elapsed().as_nanos() as u64);
        h += 2;
    }
    for (m, rb) in rec_builders.iter().enumerate() {
        rb.print_rec(&mut out, is_imputed[m]);
    }
    out
}

struct SetAlProbsScratch {
    indices: Vec<usize>,
    hashes: Vec<i32>,
    seq_probs: Vec<f32>,
    seq_probs_p1: Vec<f32>,
    alleles: Vec<i32>,
}

impl SetAlProbsScratch {
    fn new(n_markers: usize) -> SetAlProbsScratch {
        SetAlProbsScratch {
            indices: Vec::with_capacity(4),
            hashes: Vec::with_capacity(4),
            seq_probs: Vec::with_capacity(4),
            seq_probs_p1: Vec::with_capacity(4),
            alleles: vec![0; n_markers],
        }
    }
}

/// `ImputedVcfWriter.setAlProbs` (both overloads).
#[allow(clippy::too_many_arguments)]
fn set_al_probs(
    imp_data: &ImpData,
    state_probs: &StateProbs,
    rhh: &RefHapHash,
    targ_cluster: usize,
    ref_start: usize,
    clust_end: usize,
    ref_end: usize,
    al_probs: &mut [f32],
    al_off: &[u32],
    scratch: &mut SetAlProbsScratch,
) {
    scratch.indices.clear();
    scratch.hashes.clear();
    scratch.seq_probs.clear();
    scratch.seq_probs_p1.clear();
    for st in state_probs.states(targ_cluster) {
        let (val, val_p1) = (st.prob, st.prob_p1);
        let index = rhh.hap2index(st.hap);
        let hash = rhh.hash(index);
        let mut i = 0;
        while i < scratch.hashes.len() && scratch.hashes[i] != hash {
            i += 1;
        }
        if i == scratch.hashes.len() {
            scratch.indices.push(index);
            scratch.hashes.push(hash);
            scratch.seq_probs.push(val);
            scratch.seq_probs_p1.push(val_p1);
        } else {
            scratch.seq_probs[i] += val;
            scratch.seq_probs_p1[i] += val_p1;
        }
    }
    let n_seq = scratch.seq_probs.len();
    if n_seq == 1 {
        let index = scratch.indices[0];
        rhh.set_alleles(index, &mut scratch.alleles);
        for m in ref_start..ref_end {
            let mm = m - ref_start;
            al_probs[al_off[mm] as usize + scratch.alleles[mm] as usize] = 1.0;
        }
    } else {
        for j in 0..n_seq {
            let index = scratch.indices[j];
            rhh.set_alleles(index, &mut scratch.alleles);
            let prob = scratch.seq_probs[j];
            let prob_p1 = scratch.seq_probs_p1[j];
            for m in ref_start..clust_end {
                let mm = m - ref_start;
                al_probs[al_off[mm] as usize + scratch.alleles[mm] as usize] += prob;
            }
            for m in clust_end..ref_end {
                // Java: float += double  =>  add in f64, then narrow to f32
                let wt = imp_data.weight[m] as f64;
                let mm = m - ref_start;
                let slot = &mut al_probs[al_off[mm] as usize + scratch.alleles[mm] as usize];
                *slot = (*slot as f64 + (wt * prob as f64 + (1.0 - wt) * prob_p1 as f64)) as f32;
            }
        }
    }
}

/// `ImputedVcfWriter.setToObsAlleles`
fn set_to_obs_alleles(
    imp_data: &ImpData,
    window: &Window,
    ref_start: usize,
    m: usize,
    targ_hap: usize,
    a1_probs: &mut [f32],
    a2_probs: &mut [f32],
    al_off: &[u32],
) {
    let (lo, hi) = (al_off[m] as usize, al_off[m + 1] as usize);
    a1_probs[lo..hi].fill(0.0);
    a2_probs[lo..hi].fill(0.0);
    let t = window.indices.marker_to_targ_marker[ref_start + m] as usize;
    let a1 = imp_data.targ_alleles[t][targ_hap] as usize;
    let a2 = imp_data.targ_alleles[t][targ_hap + 1] as usize;
    a1_probs[lo + a1] = 1.0;
    a2_probs[lo + a2] = 1.0;
}

/// Port of `main.WindowWriter`.
#[allow(dead_code)]
pub struct WindowWriter {
    file: File,
    samples: Samples,
    ap: bool,
    gp: bool,
    fmt: Fmt,
    hom_ref: Vec<String>,
}

impl WindowWriter {
    pub fn new(out_prefix: &str, samples: Samples, ap: bool, gp: bool) -> WindowWriter {
        let fmt = Fmt::new();
        let hom_ref = hom_ref_fields(ap, gp, &fmt);
        let path = format!("{}.vcf.gz", out_prefix);
        let mut file = File::create(&path).unwrap_or_else(|e| {
            eprintln!("ERROR: cannot create output file {}: {}", path, e);
            std::process::exit(1)
        });
        let header = meta_lines(&samples.ids, crate::JAVA_EQUIV_PROGRAM, true, ap, gp, false);
        let compressed = bgzf::compress(header.as_bytes());
        file.write_all(&compressed).expect("write VCF header");
        WindowWriter {
            file,
            samples,
            ap,
            gp,
            fmt,
            hom_ref,
        }
    }

    /// `WindowWriter.printImputed`
    pub fn print_imputed(
        &mut self,
        imp_data: &ImpData,
        window: &Window,
        start: usize,
        end: usize,
        state_probs: &[StateProbs],
    ) {
        let chunks: Vec<Vec<u8>> = (0..imp_data.n_clusters)
            .into_par_iter()
            .map_init(OutScratch::default, |scratch_out, c| {
                let t0 = std::time::Instant::now();
                let text = append_records(
                    imp_data,
                    window,
                    &self.samples,
                    state_probs,
                    start,
                    end,
                    c,
                    &self.fmt,
                    &self.hom_ref,
                    scratch_out,
                );
                let t1 = std::time::Instant::now();
                let out = if text.is_empty() {
                    Vec::new()
                } else {
                    bgzf::compress(text.as_bytes())
                };
                out_add(0, (t1 - t0).as_nanos() as u64);
                out_add(1, t1.elapsed().as_nanos() as u64);
                out
            })
            .collect();
        if std::env::var("RUSTY_BEAGLE_TIMING2").is_ok() {
            eprintln!(
                "[timing2/cpu-total] format: {:.3}s  deflate: {:.3}s  (hash: {:.3}s  alprobs: {:.3}s  build: {:.3}s)",
                OUT_NANOS[0].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
                OUT_NANOS[1].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
                OUT_NANOS[2].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
                OUT_NANOS[3].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
                OUT_NANOS[4].load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9,
            );
        }
        for chunk in chunks {
            if !chunk.is_empty() {
                self.file.write_all(&chunk).expect("write VCF records");
            }
        }
    }

    /// `WindowWriter.printPhased` (no-imputation windows: GT-only records).
    /// `phased[m][hap]` supplies the phased alleles for target marker m.
    pub fn print_phased(
        &mut self,
        window: &Window,
        phased: &[Vec<i16>],
        start: usize,
        end: usize,
    ) {
        let step = 100usize;
        let mut chunk_ranges = Vec::new();
        let mut m = start;
        while m < end {
            let e = (m + step).min(end);
            chunk_ranges.push((m, e));
            m = e;
        }
        let chunks: Vec<Vec<u8>> = chunk_ranges
            .into_par_iter()
            .map(|(s, e)| {
                let mut text = String::new();
                for t in s..e {
                    let rec = &window.targ_recs[t];
                    phased_rec(&rec.marker, &phased[t], &self.samples, &mut text);
                }
                bgzf::compress(text.as_bytes())
            })
            .collect();
        for chunk in chunks {
            self.file.write_all(&chunk).expect("write VCF records");
        }
    }

    pub fn close(mut self) {
        let eof = bgzf::eof_block();
        self.file.write_all(&eof).expect("write BGZF EOF block");
        self.file.flush().expect("flush output file");
    }
}

/// `VcfWriter.appendRecords` for one phased target record (GT only).
fn phased_rec(
    marker: &crate::marker::Marker,
    alleles: &[i16],
    samples: &Samples,
    out: &mut String,
) {
    out.push_str(&marker.chrom());
    out.push('\t');
    let _ = write!(out, "{}", marker.pos);
    out.push('\t');
    out.push_str(marker.id_str());
    out.push('\t');
    out.push_str(&marker.alleles);
    out.push('\t');
    out.push('.'); // QUAL
    out.push('\t');
    out.push('.'); // FILTER
    out.push('\t');
    match marker.end_subfield() {
        Some(e) => out.push_str(&e),
        None => out.push('.'),
    }
    out.push('\t');
    out.push_str("GT");
    for (s, &diploid) in samples.is_diploid.iter().enumerate() {
        let h1 = s << 1;
        out.push('\t');
        let a1 = alleles[h1];
        let _ = write!(out, "{}", a1);
        if diploid {
            out.push('|');
            let _ = write!(out, "{}", alleles[h1 | 1]);
        }
    }
    out.push('\n');
}

/// `VcfWriter.writeMetaLines`
pub fn meta_lines(
    sample_ids: &[String],
    source: &str,
    ds: bool,
    ap: bool,
    gp: bool,
    gl: bool,
    ) -> String {
    let mut out = String::with_capacity(1024 + 8 * sample_ids.len());
    out.push_str("##fileformat=VCFv4.2\n");
    out.push_str("##filedate=");
    out.push_str(&file_date());
    out.push('\n');
    out.push_str("##source=\"");
    out.push_str(source);
    out.push_str("\"\n");
    if ds {
        out.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated ALT Allele Frequencies\">\n");
        out.push_str("##INFO=<ID=DR2,Number=A,Type=Float,Description=\"Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">\n");
        out.push_str("##INFO=<ID=IMP,Number=0,Type=Flag,Description=\"Imputed marker\">\n");
    }
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    if ds {
        out.push_str("##FORMAT=<ID=DS,Number=A,Type=Float,Description=\"estimated ALT dose [P(RA) + 2*P(AA)]\">\n");
    }
    if ap {
        out.push_str("##FORMAT=<ID=AP1,Number=A,Type=Float,Description=\"estimated ALT dose on first haplotype\">\n");
        out.push_str("##FORMAT=<ID=AP2,Number=A,Type=Float,Description=\"estimated ALT dose on second haplotype\">\n");
    }
    if gp {
        out.push_str("##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Estimated Genotype Probability\">\n");
    }
    if gl {
        out.push_str("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled Genotype Likelihood\">\n");
    }
    out.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for id in sample_ids {
        out.push('\t');
        out.push_str(id);
    }
    out.push('\n');
    out
}

/// `yyyyMMdd` for the current local date (UTC is used, matching common
/// server configurations; the field is excluded from output comparisons).
fn file_date() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let days = (secs / 86_400) as i64;
    let (y, m, d) = civil_from_days(days);
    format!("{:04}{:02}{:02}", y, m, d)
}

/// Howard Hinnant's civil-from-days algorithm.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}

static OUT_NANOS: [std::sync::atomic::AtomicU64; 5] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];

#[inline]
fn out_add(idx: usize, nanos: u64) {
    OUT_NANOS[idx].fetch_add(nanos, std::sync::atomic::Ordering::Relaxed);
}
