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

/// Port of `imp.RefHapHash`.
struct RefHapHash<'a> {
    ref_recs: &'a [Arc<RefRec>],
    i2hap: Vec<i32>,
    i2hash: Vec<i32>,
    alt_alleles: Vec<Vec<i32>>, // (marker offset, ALT allele) pairs
    start: usize,
    end: usize,
}

impl<'a> RefHapHash<'a> {
    fn new(
        state_probs: &[StateProbs],
        targ_cluster: usize,
        ref_recs: &'a [Arc<RefRec>],
        start: usize,
        end: usize,
    ) -> RefHapHash<'a> {
        assert!(start < end);
        let mut list: Vec<i32> = Vec::with_capacity(10 * state_probs.len());
        for sp in state_probs {
            for k in 0..sp.n_states(targ_cluster) {
                list.push(sp.ref_hap(targ_cluster, k));
            }
        }
        list.sort_unstable();
        list.dedup();
        let i2hap = list;
        let n = i2hap.len();
        let _ = targ_cluster;
        let mut hash = RefHapHash {
            ref_recs,
            i2hash: vec![0; n],
            alt_alleles: vec![Vec::new(); n],
            i2hap,
            start,
            end,
        };
        hash.set_hash_and_alt_alleles();
        hash
    }

    fn set_hash_and_alt_alleles(&mut self) {
        let mut rand = JavaRandom::new(self.start as i64);
        for m in self.start..self.end {
            let rec = &self.ref_recs[m];
            let marker_offset = (m - self.start) as i32;
            match &rec.alleles {
                RefAlleles::AlleleCoded { major, carriers } if *major == 0 => {
                    self.low_alt_freq_update(carriers, marker_offset, &mut rand);
                }
                _ => {
                    self.standard_update(rec, marker_offset, &mut rand);
                }
            }
        }
    }

    fn low_alt_freq_update(
        &mut self,
        carriers: &[Vec<u32>],
        marker_offset: i32,
        rand: &mut JavaRandom,
    ) {
        let n_alleles = carriers.len();
        for al in 1..n_alleles {
            let hash = rand.next_int();
            let n_copies = carriers[al].len();
            if self.i2hap.len() < n_copies {
                for i in 0..self.i2hap.len() {
                    let hap = self.i2hap[i] as u32;
                    if carriers[al].binary_search(&hap).is_ok() {
                        self.i2hash[i] = self.i2hash[i].wrapping_add(hash);
                        self.alt_alleles[i].push(marker_offset);
                        self.alt_alleles[i].push(al as i32);
                    }
                }
            } else {
                for &hap in &carriers[al] {
                    if let Ok(i) = self.i2hap.binary_search(&(hap as i32)) {
                        self.i2hash[i] = self.i2hash[i].wrapping_add(hash);
                        self.alt_alleles[i].push(marker_offset);
                        self.alt_alleles[i].push(al as i32);
                    }
                }
            }
        }
    }

    fn standard_update(&mut self, rec: &RefRec, marker_offset: i32, rand: &mut JavaRandom) {
        let n_alleles = rec.marker.n_alleles as usize;
        let mut allele_hash = vec![0i32; n_alleles];
        for a in allele_hash.iter_mut().skip(1) {
            *a = rand.next_int();
        }
        for i in 0..self.i2hap.len() {
            let allele = rec.allele(self.i2hap[i] as usize) as usize;
            if allele != 0 {
                self.i2hash[i] = self.i2hash[i].wrapping_add(allele_hash[allele]);
                self.alt_alleles[i].push(marker_offset);
                self.alt_alleles[i].push(allele as i32);
            }
        }
    }

    #[inline]
    fn hap2index(&self, hap: i32) -> usize {
        self.i2hap.binary_search(&hap).expect("state hap in hash")
    }

    #[inline]
    fn hash(&self, index: usize) -> i32 {
        self.i2hash[index]
    }

    /// `RefHapHash.setAlleles`
    fn set_alleles(&self, index: usize, alleles: &mut [i32]) {
        alleles[..self.end - self.start].fill(0);
        let il = &self.alt_alleles[index];
        let mut j = 0;
        while j < il.len() {
            alleles[il[j] as usize] = il[j + 1];
            j += 2;
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
    let rhh = RefHapHash::new(state_probs, targ_cluster, ref_recs, ref_start, ref_end);
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
    let mut a1_probs: Vec<Vec<f32>> = (ref_start..ref_end)
        .map(|m| vec![0.0f32; ref_recs[m].marker.n_alleles as usize])
        .collect();
    let mut a2_probs: Vec<Vec<f32>> = a1_probs.clone();
    let is_imputed: Vec<bool> = (ref_start..ref_end)
        .map(|m| window.indices.marker_to_targ_marker[m] == -1)
        .collect();

    let mut scratch = SetAlProbsScratch::new(n_markers);
    let n = state_probs.len();
    let mut h = 0;
    while h < n {
        let is_diploid = samples.is_diploid[h >> 1];
        set_al_probs(
            imp_data,
            &state_probs[h],
            &rhh,
            targ_cluster,
            ref_start,
            clust_end,
            ref_end,
            &mut a1_probs,
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
            &mut scratch,
        );
        if is_diploid {
            for m in 0..n_markers {
                if !is_imputed[m] {
                    set_to_obs_alleles(imp_data, window, ref_start, m, h, &mut a1_probs, &mut a2_probs);
                }
                let (a1, a2) = (&mut a1_probs[m], &mut a2_probs[m]);
                rec_builders[m].add_sample_data2(a1, a2);
                a1.fill(0.0);
                a2.fill(0.0);
            }
        } else {
            for m in 0..n_markers {
                if !is_imputed[m] {
                    set_to_obs_alleles(imp_data, window, ref_start, m, h, &mut a1_probs, &mut a2_probs);
                }
                rec_builders[m].add_sample_data1(&mut a1_probs[m]);
                a1_probs[m].fill(0.0);
                a2_probs[m].fill(0.0);
            }
        }
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
    al_probs: &mut [Vec<f32>],
    scratch: &mut SetAlProbsScratch,
) {
    scratch.indices.clear();
    scratch.hashes.clear();
    scratch.seq_probs.clear();
    scratch.seq_probs_p1.clear();
    for j in 0..state_probs.n_states(targ_cluster) {
        let hap = state_probs.ref_hap(targ_cluster, j);
        let val = state_probs.probs(targ_cluster, j);
        let val_p1 = state_probs.probs_p1(targ_cluster, j);
        let index = rhh.hap2index(hap);
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
            al_probs[mm][scratch.alleles[mm] as usize] = 1.0;
        }
    } else {
        for j in 0..n_seq {
            let index = scratch.indices[j];
            rhh.set_alleles(index, &mut scratch.alleles);
            let prob = scratch.seq_probs[j];
            let prob_p1 = scratch.seq_probs_p1[j];
            for m in ref_start..clust_end {
                let mm = m - ref_start;
                al_probs[mm][scratch.alleles[mm] as usize] += prob;
            }
            for m in clust_end..ref_end {
                // Java: float += double  =>  add in f64, then narrow to f32
                let wt = imp_data.weight[m] as f64;
                let mm = m - ref_start;
                let slot = &mut al_probs[mm][scratch.alleles[mm] as usize];
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
    a1_probs: &mut [Vec<f32>],
    a2_probs: &mut [Vec<f32>],
) {
    a1_probs[m].fill(0.0);
    a2_probs[m].fill(0.0);
    let t = window.indices.marker_to_targ_marker[ref_start + m] as usize;
    let a1 = imp_data.targ_alleles[t][targ_hap] as usize;
    let a2 = imp_data.targ_alleles[t][targ_hap + 1] as usize;
    a1_probs[m][a1] = 1.0;
    a2_probs[m][a2] = 1.0;
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
            .map(|c| {
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
                );
                if text.is_empty() {
                    Vec::new()
                } else {
                    bgzf::compress(text.as_bytes())
                }
            })
            .collect();
        for chunk in chunks {
            if !chunk.is_empty() {
                self.file.write_all(&chunk).expect("write VCF records");
            }
        }
    }

    /// `WindowWriter.printPhased` (no-imputation windows: GT-only records).
    pub fn print_phased(&mut self, window: &Window, start: usize, end: usize) {
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
                    phased_rec(rec, &self.samples, &mut text);
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
fn phased_rec(rec: &crate::vcfio::GtRec, samples: &Samples, out: &mut String) {
    let marker = &rec.marker;
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
        let a1 = rec.alleles[h1];
        let _ = write!(out, "{}", a1);
        if diploid {
            out.push('|');
            let _ = write!(out, "{}", rec.alleles[h1 | 1]);
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
