//! Ports of `phase.HmmUpdater`, `phase.BasicPhaseStates`, `phase.PhaseBaum2`,
//! `phase.HmmParamData`, `phase.ParamEstimates`, and the `phase.PhaseLS`
//! stage-1 driver.

use crate::bits::BitArray;
use crate::codedsteps::CodedSteps;
use crate::javautil::{IntIntMap, JavaOrd, JavaPriorityQueue, JavaRandom};
use crate::par::Par;
use crate::phaseibs::PbwtPhaseIbs;
use crate::phasedata::{
    swap_rate_increment, EstPhase, FixedPhaseData, MarkerClusterInfo, PhaseData, SamplePhase,
    CLUST_MASKED_HET, CLUST_MISSING_GT, CLUST_UNPHASED_HET,
};
use crate::xref::XRefGT;
use std::sync::Arc;

const NIL: i32 = -103;

/// Port of `beagleutil.CompHapSegment` (with `startMarker`).
#[derive(Clone)]
pub struct CompHapSegment {
    pub hap: i32,
    pub start_marker: usize,
    pub last_ibs_step: i32,
    pub comp_hap_index: usize,
}

impl JavaOrd for CompHapSegment {
    #[inline]
    fn compare_to(&self, other: &Self) -> std::cmp::Ordering {
        self.last_ibs_step.cmp(&other.last_ibs_step)
    }
}

// ---------------------------------------------------------------------------
// HmmUpdater

/// `HmmUpdater.fwdUpdate`
pub fn fwd_update(
    fwd: &mut [f32],
    fwd_sum: f32,
    p_switch: f32,
    p_mismatch: &[f32; 2],
    mismatch: &[u8],
    n_states: usize,
) -> f32 {
    let shift = p_switch / n_states as f32;
    let scale = (1.0f32 - p_switch) / fwd_sum;
    let mut sum = 0.0f32;
    for k in 0..n_states {
        fwd[k] = p_mismatch[mismatch[k] as usize] * (scale * fwd[k] + shift);
        sum += fwd[k];
    }
    sum
}

/// `HmmUpdater.bwdUpdate`
pub fn bwd_update(
    bwd: &mut [f32],
    p_switch: f32,
    p_mismatch: &[f32; 2],
    mismatch: &[u8],
    n_states: usize,
) {
    let mut sum = 0.0f32;
    for k in 0..n_states {
        bwd[k] *= p_mismatch[mismatch[k] as usize];
        sum += bwd[k];
    }
    let shift = p_switch / n_states as f32;
    let scale = (1.0f32 - p_switch) / sum;
    for b in bwd[..n_states].iter_mut() {
        *b = scale * *b + shift;
    }
}

// ---------------------------------------------------------------------------
// BasicPhaseStates

pub struct BasicPhaseStates<'a> {
    ibs_haps: &'a PbwtPhaseIbs,
    fpd: &'a FixedPhaseData,
    all_haps: &'a XRefGT,
    n_markers: usize,
    max_states: usize,
    min_steps: usize,
    hap_to_last_ibs_step: IntIntMap,
    q: JavaPriorityQueue<CompHapSegment>,
    comp_haps: Vec<BitArray>,
    it_seed: i64,
}

impl<'a> BasicPhaseStates<'a> {
    pub fn new(
        pd: &PhaseData,
        ibs_haps: &'a PbwtPhaseIbs,
        all_haps: &'a XRefGT,
        fpd: &'a FixedPhaseData,
        max_states: usize,
    ) -> BasicPhaseStates<'a> {
        let n_markers = all_haps.n_markers();
        let phase_step = fpd.ibs_step;
        let min_steps = std::cmp::max(200, (1.0f32 / phase_step).ceil() as usize);
        let n_bits = all_haps.layout.total_bits();
        BasicPhaseStates {
            ibs_haps,
            fpd,
            all_haps,
            n_markers,
            max_states,
            min_steps,
            hap_to_last_ibs_step: IntIntMap::with_capacity(max_states),
            q: JavaPriorityQueue::new(),
            comp_haps: (0..max_states).map(|_| BitArray::new(n_bits)).collect(),
            it_seed: pd.it_seed(),
        }
    }

    /// `BasicPhaseStates.ibsStates(mc, refAtMissingGT, nMismatches)`
    /// (three-row variant used by PhaseBaum2)
    pub fn ibs_states_clusters(
        &mut self,
        sp: &SamplePhase,
        mc: &MarkerClusterInfo,
        ref_at_missing_gt: &mut Vec<Vec<i32>>,
        mismatches: &mut [Vec<u8>; 3],
        row_len: usize,
    ) -> usize {
        let n_comp_haps = self.set_comp_ref_haps(sp.sample);
        self.copy_data_clusters(sp, mc, n_comp_haps, ref_at_missing_gt, mismatches, row_len);
        n_comp_haps
    }

    /// `BasicPhaseStates.ibsStates(sample, nMismatches)` (two-row variant used
    /// by HmmParamData)
    pub fn ibs_states_markers(
        &mut self,
        sample: usize,
        mismatches: &mut [Vec<u8>; 2],
        row_len: usize,
    ) -> usize {
        let n_comp_haps = self.set_comp_ref_haps(sample);
        let h1 = sample << 1;
        let h2 = h1 | 1;
        for m in 0..self.n_markers {
            let a1 = self.all_haps.allele(m, h1);
            let a2 = self.all_haps.allele(m, h2);
            let row = m * row_len;
            for j in 0..n_comp_haps {
                let ref_allele = self.all_haps.layout.allele(&self.comp_haps[j], m);
                mismatches[0][row + j] = if ref_allele == a1 { 0 } else { 1 };
                mismatches[1][row + j] = if ref_allele == a2 { 0 } else { 1 };
            }
        }
        n_comp_haps
    }

    fn set_comp_ref_haps(&mut self, sample: usize) -> usize {
        let h1 = sample << 1;
        let h2 = h1 | 1;
        self.q.clear();
        self.hap_to_last_ibs_step.clear();
        for step in 0..self.fpd.stage1_steps.size() {
            let ibs_hap1 = self.ibs_haps.ibs_hap(h1, step);
            if ibs_hap1 >= 0 {
                self.add_ibs_hap(ibs_hap1, step as i32);
            }
            let ibs_hap2 = self.ibs_haps.ibs_hap(h2, step);
            if ibs_hap2 >= 0 {
                self.add_ibs_hap(ibs_hap2, step as i32);
            }
        }
        if self.q.is_empty() {
            self.fill_q_with_random_haps(sample);
        }
        self.copy_final_ref_segs()
    }

    fn add_ibs_hap(&mut self, ibs_hap: i32, step: i32) {
        if self.hap_to_last_ibs_step.get(ibs_hap, NIL) == NIL {
            self.update_head_of_q();
            let recycle = self.q.len() == self.max_states
                || (!self.q.is_empty()
                    && step - self.q.peek().unwrap().last_ibs_step >= self.min_steps as i32);
            if recycle {
                let head = self.q.poll().unwrap();
                let index = head.comp_hap_index;
                let prev_hap = head.hap;
                let prev_start = head.start_marker;
                let next_start = self
                    .fpd
                    .stage1_steps
                    .start((((head.last_ibs_step + step) as u32) >> 1) as usize);
                self.hap_to_last_ibs_step.remove(head.hap);
                self.all_haps.copy_to(
                    prev_hap as usize,
                    prev_start,
                    next_start,
                    &mut self.comp_haps[index],
                );
                self.q.offer(CompHapSegment {
                    hap: ibs_hap,
                    start_marker: next_start,
                    last_ibs_step: step,
                    comp_hap_index: index,
                });
            } else {
                let index = self.q.len();
                self.q.offer(CompHapSegment {
                    hap: ibs_hap,
                    start_marker: 0,
                    last_ibs_step: step,
                    comp_hap_index: index,
                });
            }
        }
        self.hap_to_last_ibs_step.put(ibs_hap, step);
    }

    fn update_head_of_q(&mut self) {
        if let Some(head) = self.q.peek() {
            let mut last_ibs_step = self.hap_to_last_ibs_step.get(head.hap, NIL);
            let mut head_step = head.last_ibs_step;
            while head_step != last_ibs_step {
                let mut head = self.q.poll().unwrap();
                head.last_ibs_step = last_ibs_step;
                self.q.offer(head);
                let head_ref = self.q.peek().unwrap();
                head_step = head_ref.last_ibs_step;
                last_ibs_step = self.hap_to_last_ibs_step.get(head_ref.hap, NIL);
            }
        }
    }

    fn copy_final_ref_segs(&mut self) -> usize {
        let n_comp_haps = self.q.len();
        while let Some(head) = self.q.poll() {
            let index = head.comp_hap_index;
            self.all_haps.copy_to(
                head.hap as usize,
                head.start_marker,
                self.n_markers,
                &mut self.comp_haps[index],
            );
        }
        n_comp_haps
    }

    fn copy_data_clusters(
        &self,
        sp: &SamplePhase,
        mc: &MarkerClusterInfo,
        n_comp_haps: usize,
        ref_at_missing_gt: &mut [Vec<i32>],
        mismatches: &mut [Vec<u8>; 3],
        row_len: usize,
    ) {
        let layout = &self.all_haps.layout;
        let hap1 = &sp.hap1;
        let hap2 = &sp.hap2;
        let mut miss_index = 0usize;
        let n_clusters = mc.n_clusters();
        for c in 0..n_clusters {
            let row = c * row_len;
            for r in mismatches.iter_mut() {
                r[row..row + n_comp_haps].fill(0);
            }
            let m_start = mc.cluster_start(c);
            let m_end = mc.cluster_end(c);
            let ct = sp.clust_type_at(c);
            if ct == CLUST_MISSING_GT || ct == CLUST_MASKED_HET {
                debug_assert_eq!(m_end - m_start, 1);
                let ref_alleles = &mut ref_at_missing_gt[miss_index];
                miss_index += 1;
                for j in 0..n_comp_haps {
                    ref_alleles[j] = layout.allele(&self.comp_haps[j], m_start);
                }
            } else {
                let b_start = layout.sum_hap_bits[m_start] as usize;
                let b_end = layout.sum_hap_bits[m_end] as usize;
                if hap1.equal(hap2, b_start, b_end) {
                    for j in 0..n_comp_haps {
                        if !hap1.equal(&self.comp_haps[j], b_start, b_end) {
                            mismatches[0][row + j] = 1;
                            mismatches[1][row + j] = 1;
                            mismatches[2][row + j] = 1;
                        }
                    }
                } else {
                    // cluster contains a heterozygote genotype
                    for j in 0..n_comp_haps {
                        if !hap1.equal(&self.comp_haps[j], b_start, b_end) {
                            mismatches[1][row + j] = 1;
                        }
                        if !hap2.equal(&self.comp_haps[j], b_start, b_end) {
                            mismatches[2][row + j] = 1;
                        }
                    }
                }
            }
        }
    }

    fn fill_q_with_random_haps(&mut self, sample: usize) {
        debug_assert!(self.q.is_empty());
        let n_haps = self.all_haps.n_haps();
        let n_states = std::cmp::min(n_haps - 2, self.max_states);
        if n_states == 0 {
            eprintln!("ERROR: there is only one sample");
            std::process::exit(1);
        }
        let mut rand = JavaRandom::new(self.it_seed.wrapping_add(sample as i64));
        let ibs_step = 0i32;
        let start_marker = 0usize;
        let mut comp_hap_index = 0usize;
        for _ in 0..n_states {
            let mut h = rand.next_int_bound(n_haps as i32);
            while (h >> 1) as usize == sample {
                h = rand.next_int_bound(n_haps as i32);
            }
            if self.hap_to_last_ibs_step.get(h, NIL) == NIL {
                self.q.offer(CompHapSegment {
                    hap: h,
                    start_marker,
                    last_ibs_step: ibs_step,
                    comp_hap_index,
                });
                comp_hap_index += 1;
                self.hap_to_last_ibs_step.put(h, start_marker as i32);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// PhaseBaum2

pub struct PhaseBaum2<'a> {
    pd: &'a PhaseData,
    par: &'a Par,
    burnin: bool,
    lr_threshold: f32,
    mask_trailing_hets: bool,
    est_phase: &'a EstPhase,
    #[allow(dead_code)]
    n_markers: usize,
    ref_alleles: Vec<Vec<i32>>,
    mismatches: [Vec<u8>; 3],
    p_mismatch: f32,
    em_probs: [f32; 2],
    max_states: usize,
    states: BasicPhaseStates<'a>,
    n_states: usize,
    fwd: [Vec<f32>; 3],
    bwd: [Vec<f32>; 3],
    fwd_sums: [f32; 3],
    bwd_miss1: Vec<Vec<f32>>,
    bwd_miss2: Vec<Vec<f32>>,
    bwd_het1: Vec<Vec<f32>>,
    bwd_het2: Vec<Vec<f32>>,
    swap_haps: bool,
    n_swaps: usize,
}

impl<'a> PhaseBaum2<'a> {
    pub fn new(
        pd: &'a PhaseData,
        par: &'a Par,
        ibs: &'a PbwtPhaseIbs,
        all_haps: &'a XRefGT,
    ) -> PhaseBaum2<'a> {
        let fpd = pd.fpd();
        let burnin = pd.it < par.burnin as usize;
        let lr_threshold = pd.lr_threshold;
        let lr_mask_threshold = 50f32;
        let max_states = par.phase_states as usize;
        let n_markers = fpd.n_stage1_markers();
        let p_mismatch = pd.p_mismatch;
        PhaseBaum2 {
            pd,
            par,
            burnin,
            lr_threshold,
            mask_trailing_hets: lr_threshold < lr_mask_threshold,
            est_phase: &pd.est_phase,
            n_markers,
            ref_alleles: Vec::new(),
            mismatches: [Vec::new(), Vec::new(), Vec::new()],
            p_mismatch,
            em_probs: [1.0 - p_mismatch, p_mismatch],
            max_states,
            states: BasicPhaseStates::new(pd, ibs, all_haps, fpd, max_states),
            n_states: 0,
            fwd: [Vec::new(), Vec::new(), Vec::new()],
            bwd: [Vec::new(), Vec::new(), Vec::new()],
            fwd_sums: [0.0; 3],
            bwd_miss1: Vec::new(),
            bwd_miss2: Vec::new(),
            bwd_het1: Vec::new(),
            bwd_het2: Vec::new(),
            swap_haps: false,
            n_swaps: 0,
        }
    }

    /// `PhaseBaum2.phase(sample)`
    pub fn phase(&mut self, sample: usize) {
        let _ = self.par;
        let mut sp = self.est_phase.take(sample);
        if self.mask_trailing_hets {
            sp.mask_trailing_unphased_hets(&self.pd.fpd().stage1_positions);
        }
        let n_unph_hets = sp.n_unphased();
        let n_masked_hets = sp.n_masked();
        let n_missing_or_masked = sp.n_missing() + n_masked_hets;
        if n_missing_or_masked > 0 || n_unph_hets > 0 {
            self.n_swaps = 0;
            self.swap_haps = false;
            let mc = MarkerClusterInfo::new(&sp, &self.pd.p_recomb);
            self.ensure_capacity(n_unph_hets, n_missing_or_masked, mc.n_clusters());
            self.n_states = self.states.ibs_states_clusters(
                &sp,
                &mc,
                &mut self.ref_alleles,
                &mut self.mismatches,
                self.max_states,
            );
            self.bwd_alg(&sp, &mc);
            self.fwd_alg(&mut sp, &mc);
            swap_rate_increment(n_unph_hets, self.n_swaps);
        }
        self.est_phase.put(sample, sp);
    }

    fn ensure_capacity(&mut self, n_unph: usize, n_miss: usize, n_clusters: usize) {
        while self.ref_alleles.len() < n_miss {
            self.ref_alleles.push(vec![0i32; self.max_states]);
            self.bwd_miss1.push(vec![0f32; self.max_states]);
            self.bwd_miss2.push(vec![0f32; self.max_states]);
        }
        while self.bwd_het1.len() < n_unph {
            self.bwd_het1.push(vec![0f32; self.max_states]);
            self.bwd_het2.push(vec![0f32; self.max_states]);
        }
        let need = n_clusters * self.max_states;
        for r in self.mismatches.iter_mut() {
            if r.len() < need {
                r.resize(need, 0);
            }
        }
        for f in self.fwd.iter_mut() {
            if f.len() < self.max_states {
                f.resize(self.max_states, 0.0);
            }
        }
        for b in self.bwd.iter_mut() {
            if b.len() < self.max_states {
                b.resize(self.max_states, 0.0);
            }
        }
    }

    fn bwd_alg(&mut self, sp: &SamplePhase, mc: &MarkerClusterInfo) {
        let n = self.n_states;
        let mut miss_index = (sp.n_missing() + sp.n_masked()) as isize - 1;
        let mut unph_index = sp.n_unphased() as isize - 1;
        // initializeBwdFields
        self.bwd[0][..n].fill(1.0 / n as f32);
        let (b0, rest) = self.bwd.split_at_mut(1);
        rest[0][..n].copy_from_slice(&b0[0][..n]);
        rest[1][..n].copy_from_slice(&b0[0][..n]);
        let last_cluster = mc.n_clusters() - 1;
        if is_missing_or_masked(sp, last_cluster) {
            self.bwd_miss1[miss_index as usize][..n].copy_from_slice(&self.bwd[0][..n]);
            self.bwd_miss2[miss_index as usize][..n].copy_from_slice(&self.bwd[0][..n]);
            miss_index -= 1;
        }
        for c in (0..last_cluster).rev() {
            self.bwd_step(mc, c);
            if is_missing_or_masked(sp, c) {
                self.bwd_miss1[miss_index as usize][..n].copy_from_slice(&self.bwd[1][..n]);
                self.bwd_miss2[miss_index as usize][..n].copy_from_slice(&self.bwd[2][..n]);
                miss_index -= 1;
            }
            if sp.clust_type_at(c + 1) == CLUST_UNPHASED_HET {
                self.bwd_het1[unph_index as usize][..n].copy_from_slice(&self.bwd[1][..n]);
                self.bwd_het2[unph_index as usize][..n].copy_from_slice(&self.bwd[2][..n]);
                let (b0, rest) = self.bwd.split_at_mut(1);
                rest[0][..n].copy_from_slice(&b0[0][..n]);
                rest[1][..n].copy_from_slice(&b0[0][..n]);
                unph_index -= 1;
            }
        }
        debug_assert_eq!(miss_index, -1);
        debug_assert_eq!(unph_index, -1);
    }

    fn bwd_step(&mut self, mc: &MarkerClusterInfo, cluster: usize) {
        let c_p1 = cluster + 1;
        let p_rec = mc.p_recomb[c_p1];
        let mut clust_em = (mc.cluster_end(c_p1) - mc.cluster_start(c_p1)) as f32 * self.p_mismatch;
        if clust_em >= 0.5 {
            clust_em = 0.5;
        }
        self.em_probs[1] = clust_em;
        self.em_probs[0] = 1.0 - clust_em;
        let row = c_p1 * self.max_states;
        let n = self.n_states;
        for i in 0..3 {
            bwd_update(
                &mut self.bwd[i],
                p_rec,
                &self.em_probs,
                &self.mismatches[i][row..row + n],
                n,
            );
        }
    }

    fn fwd_alg(&mut self, sp: &mut SamplePhase, mc: &MarkerClusterInfo) {
        let n = self.n_states;
        let mut miss_index = 0usize;
        let mut unph_het_index = 0usize;
        // initializeFwdFields
        self.fwd[0][..n].fill(1.0 / n as f32);
        let (f0, rest) = self.fwd.split_at_mut(1);
        rest[0][..n].copy_from_slice(&f0[0][..n]);
        rest[1][..n].copy_from_slice(&f0[0][..n]);
        self.fwd_sums = [1.0, 1.0, 1.0];
        for c in 0..mc.n_clusters() {
            if sp.clust_type_at(c) == CLUST_UNPHASED_HET {
                self.phase_het(sp, unph_het_index, c);
                unph_het_index += 1;
                if self.swap_haps {
                    let swap_end = if unph_het_index < mc.unph_het_clusters.len() {
                        mc.unph_het_clusters[unph_het_index]
                    } else {
                        mc.n_clusters()
                    };
                    self.do_swap_haps(sp, mc, c, swap_end);
                }
                let (f0, rest) = self.fwd.split_at_mut(1);
                rest[0][..n].copy_from_slice(&f0[0][..n]);
                rest[1][..n].copy_from_slice(&f0[0][..n]);
                self.fwd_sums[1] = self.fwd_sums[0];
                self.fwd_sums[2] = self.fwd_sums[0];
            }
            self.fwd_step(mc, c);
            if is_missing_or_masked(sp, c) {
                self.impute_alleles(sp, mc, c, miss_index);
                miss_index += 1;
            }
        }
    }

    fn fwd_step(&mut self, mc: &MarkerClusterInfo, cluster: usize) {
        let p_rec = mc.p_recomb[cluster];
        let mut clust_em =
            (mc.cluster_end(cluster) - mc.cluster_start(cluster)) as f32 * self.p_mismatch;
        if clust_em >= 0.5 {
            clust_em = 0.5;
        }
        self.em_probs[1] = clust_em;
        self.em_probs[0] = 1.0 - clust_em;
        let row = cluster * self.max_states;
        let n = self.n_states;
        for i in 0..3 {
            self.fwd_sums[i] = fwd_update(
                &mut self.fwd[i],
                self.fwd_sums[i],
                p_rec,
                &self.em_probs,
                &self.mismatches[i][row..row + n],
                n,
            );
        }
    }

    fn do_swap_haps(
        &mut self,
        sp: &mut SamplePhase,
        mc: &MarkerClusterInfo,
        start_clust: usize,
        end_clust: usize,
    ) {
        for c in start_clust..end_clust {
            let row = c * self.max_states;
            let (m1, m2) = {
                let (a, b) = self.mismatches.split_at_mut(2);
                (&mut a[1], &mut b[0])
            };
            for j in 0..self.n_states {
                std::mem::swap(&mut m1[row + j], &mut m2[row + j]);
            }
        }
        let layout = &self.pd.fpd().stage1_layout;
        sp.swap_haps(layout, mc.cluster_start(start_clust), mc.cluster_end(end_clust - 1));
    }

    fn impute_alleles(
        &mut self,
        sp: &mut SamplePhase,
        mc: &MarkerClusterInfo,
        cluster: usize,
        miss_index: usize,
    ) {
        debug_assert_eq!(mc.cluster_end(cluster) - mc.cluster_start(cluster), 1);
        let n = self.n_states;
        let (probs1, probs2) = if self.swap_haps {
            (&mut self.bwd_miss2[miss_index], &mut self.bwd_miss1[miss_index])
        } else {
            (&mut self.bwd_miss1[miss_index], &mut self.bwd_miss2[miss_index])
        };
        let ref_al = &self.ref_alleles[miss_index];
        for k in 0..n {
            probs1[k] *= self.fwd[1][k];
            probs2[k] *= self.fwd[2][k];
        }
        let marker = mc.cluster_start(cluster);
        let layout = &self.pd.fpd().stage1_layout;
        let n_alleles = layout.n_alleles[marker] as usize;
        let mut al_freq1 = vec![0.0f32; n_alleles];
        let mut al_freq2 = vec![0.0f32; n_alleles];
        for k in 0..n {
            al_freq1[ref_al[k] as usize] += probs1[k];
            al_freq2[ref_al[k] as usize] += probs2[k];
        }
        let ct = sp.clust_type_at(cluster);
        if ct == CLUST_MISSING_GT {
            let mut a1 = 0usize;
            let mut a2 = 0usize;
            for j in 1..n_alleles {
                if al_freq1[j] > al_freq1[a1] {
                    a1 = j;
                }
                if al_freq2[j] > al_freq2[a2] {
                    a2 = j;
                }
            }
            sp.set_allele1(layout, marker, a1 as i32);
            sp.set_allele2(layout, marker, a2 as i32);
        } else if ct == CLUST_MASKED_HET {
            let a1 = sp.allele1(layout, marker);
            let a2 = sp.allele2(layout, marker);
            debug_assert_ne!(a1, a2);
            let p_no_switch = al_freq1[a1 as usize] * al_freq2[a2 as usize];
            let p_switch = al_freq1[a2 as usize] * al_freq2[a1 as usize];
            if p_switch > p_no_switch {
                sp.set_allele1(layout, marker, a2);
                sp.set_allele2(layout, marker, a1);
                if p_switch >= self.lr_threshold * p_no_switch {
                    sp.mark_masked_het_cluster_as_phased(cluster);
                }
            } else if p_no_switch >= self.lr_threshold * p_switch {
                sp.mark_masked_het_cluster_as_phased(cluster);
            }
        }
    }

    fn phase_het(&mut self, sp: &mut SamplePhase, unph_het_index: usize, cluster: usize) {
        let b1 = &self.bwd_het1[unph_het_index];
        let b2 = &self.bwd_het2[unph_het_index];
        let mut p11 = 0.0f32;
        let mut p12 = 0.0f32;
        let mut p21 = 0.0f32;
        let mut p22 = 0.0f32;
        for k in 0..self.n_states {
            p11 += self.fwd[1][k] * b1[k];
            p12 += self.fwd[1][k] * b2[k];
            p21 += self.fwd[2][k] * b1[k];
            p22 += self.fwd[2][k] * b2[k];
        }
        let num = p11 * p22;
        let den = p12 * p21;
        let last_swap_haps = self.swap_haps;
        self.swap_haps = num < den;
        if self.swap_haps != last_swap_haps {
            self.n_swaps += 1;
        }
        if !self.burnin
            && (num >= den * self.lr_threshold
                || (self.swap_haps && den >= num * self.lr_threshold))
        {
            sp.mark_unphased_het_cluster_as_phased(cluster);
        }
    }
}

#[inline]
fn is_missing_or_masked(sp: &SamplePhase, cluster: usize) -> bool {
    let ct = sp.clust_type_at(cluster);
    ct == CLUST_MISSING_GT || ct == CLUST_MASKED_HET
}

// ---------------------------------------------------------------------------
// ParamEstimates / HmmParamData

#[derive(Default)]
pub struct ParamEstimates {
    switch_data: std::sync::Mutex<Vec<(f64, f64)>>, // (genDistance, switchProb)
    mismatch_data: std::sync::Mutex<Vec<(i64, f64)>>, // (markerCnt, pMismatchSum)
}

impl ParamEstimates {
    pub fn add_mismatch_data(&self, marker_cnt: i64, p_mismatch_sum: f64) {
        if marker_cnt > 0 && p_mismatch_sum > 0.0 && p_mismatch_sum.is_finite() {
            self.mismatch_data
                .lock()
                .unwrap()
                .push((marker_cnt, p_mismatch_sum));
        }
    }

    pub fn add_switch_data(&self, gen_distances: f64, switch_probs: f64) {
        if gen_distances > 0.0
            && switch_probs > 0.0
            && gen_distances.is_finite()
            && switch_probs.is_finite()
        {
            self.switch_data
                .lock()
                .unwrap()
                .push((gen_distances, switch_probs));
        }
    }

    /// `ParamEstimates.pMismatch()` (sorted before summation)
    pub fn p_mismatch(&self) -> f32 {
        let mut mda = self.mismatch_data.lock().unwrap().clone();
        // MismatchData.compareTo: by pMismatchSum, then markerCnt
        mda.sort_by(|a, b| {
            a.1.partial_cmp(&b.1)
                .unwrap()
                .then(a.0.cmp(&b.0))
        });
        let mut sum_markers = 0i64;
        let mut sum_p_mismatch = 0f64;
        for (cnt, sum) in mda {
            sum_markers += cnt;
            sum_p_mismatch += sum;
        }
        if sum_markers == 0 {
            f32::NAN
        } else {
            (sum_p_mismatch / sum_markers as f64) as f32
        }
    }

    /// `ParamEstimates.recombIntensity()` (sorted before summation)
    pub fn recomb_intensity(&self) -> f32 {
        let mut rda = self.switch_data.lock().unwrap().clone();
        // RecombData.compareTo: by genDistance, then switchProb
        rda.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap()
                .then(a.1.partial_cmp(&b.1).unwrap())
        });
        let mut sum_switches = 0f64;
        let mut sum_distances = 0f64;
        for (dist, sw) in rda {
            sum_switches += sw;
            sum_distances += dist;
        }
        if sum_distances == 0.0 {
            f32::NAN
        } else {
            (sum_switches / sum_distances) as f32
        }
    }
}

/// Port of `phase.HmmParamData`.
pub struct HmmParamData<'a> {
    fpd: &'a FixedPhaseData,
    n_markers: usize,
    p_recomb: &'a [f32],
    al_match: [Vec<u8>; 2],
    states: BasicPhaseStates<'a>,
    fwd: Vec<f32>,
    bwd: Vec<f32>,
    saved_bwd: Vec<f32>, // nMarkers x maxStates
    em_probs: [f32; 2],
    max_states: usize,
    mismatch_cnt: i64,
    sum_mismatch_prob: f64,
    sum_gen_dist: f64,
    sum_switch_prob: f64,
}

impl<'a> HmmParamData<'a> {
    pub fn new(
        pd: &'a PhaseData,
        par: &'a Par,
        ibs: &'a PbwtPhaseIbs,
        all_haps: &'a XRefGT,
    ) -> HmmParamData<'a> {
        let fpd = pd.fpd();
        let max_states = par.phase_states as usize;
        let n_markers = fpd.n_stage1_markers();
        let p_mismatch = pd.p_mismatch;
        HmmParamData {
            fpd,
            n_markers,
            p_recomb: &pd.p_recomb,
            al_match: [
                vec![0u8; n_markers * max_states],
                vec![0u8; n_markers * max_states],
            ],
            states: BasicPhaseStates::new(pd, ibs, all_haps, fpd, max_states),
            fwd: vec![0f32; max_states],
            bwd: vec![0f32; max_states],
            saved_bwd: vec![0f32; n_markers * max_states],
            em_probs: [1.0 - p_mismatch, p_mismatch],
            max_states,
            mismatch_cnt: 0,
            sum_mismatch_prob: 0.0,
            sum_gen_dist: 0.0,
            sum_switch_prob: 0.0,
        }
    }

    pub fn sum_switch_probs(&self) -> f64 {
        self.sum_switch_prob
    }

    pub fn add_estimation_data(&mut self, param_est: &ParamEstimates) {
        param_est.add_mismatch_data(self.mismatch_cnt, self.sum_mismatch_prob);
        param_est.add_switch_data(self.sum_gen_dist, self.sum_switch_prob);
        self.mismatch_cnt = 0;
        self.sum_mismatch_prob = 0.0;
        self.sum_gen_dist = 0.0;
        self.sum_switch_prob = 0.0;
    }

    pub fn update(&mut self, sample: usize) {
        let n_states =
            self.states
                .ibs_states_markers(sample, &mut self.al_match, self.max_states);
        if n_states > 1 {
            self.get_param_data(0, n_states);
            self.get_param_data(1, n_states);
        }
    }

    fn get_param_data(&mut self, which: usize, n_states: usize) {
        let n_markers = self.n_markers;
        self.bwd[..n_states].fill(1.0);
        let last_row = (n_markers - 1) * self.max_states;
        self.saved_bwd[last_row..last_row + n_states].fill(1.0);
        for m in (0..n_markers - 1).rev() {
            let m_p1 = m + 1;
            let row_p1 = m_p1 * self.max_states;
            bwd_update(
                &mut self.bwd,
                self.p_recomb[m_p1],
                &self.em_probs,
                &self.al_match[which][row_p1..row_p1 + n_states],
                n_states,
            );
            let row = m * self.max_states;
            self.saved_bwd[row..row + n_states].copy_from_slice(&self.bwd[..n_states]);
        }
        let h_factor = n_states as f32 / (n_states as f32 - 1.0);
        self.fwd[..n_states].fill(1.0 / n_states as f32);
        let mut sum = 1.0f32;
        for m in 0..n_markers {
            sum = self.fwd_update_est(m, which, n_states, sum, h_factor);
        }
    }

    fn fwd_update_est(
        &mut self,
        m: usize,
        which: usize,
        n_states: usize,
        last_sum: f32,
        h_factor: f32,
    ) -> f32 {
        let p_switch = self.p_recomb[m];
        let shift = p_switch / n_states as f32;
        let scale = (1.0f32 - p_switch) / last_sum;
        let no_switch_scale = ((1.0f32 - p_switch) + shift) / last_sum;
        let mut joint_state_sum = 0.0f32;
        let mut state_sum = 0.0f32;
        let row = m * self.max_states;
        let bwd_m = &self.saved_bwd[row..row + n_states];
        let al_discord = &self.al_match[which][row..row + n_states];
        let mut fwd_sum = 0.0f32;
        let mut mismatch_sum = 0.0f32;
        for k in 0..n_states {
            let em = self.em_probs[al_discord[k] as usize];
            joint_state_sum += bwd_m[k] * em * no_switch_scale * self.fwd[k];
            self.fwd[k] = em * (scale * self.fwd[k] + shift);
            fwd_sum += self.fwd[k];
            let state_prob = self.fwd[k] * bwd_m[k];
            state_sum += state_prob;
            if al_discord[k] > 0 {
                mismatch_sum += state_prob;
            }
        }
        self.mismatch_cnt += 1;
        self.sum_mismatch_prob += (mismatch_sum / state_sum) as f64;
        let switch_prob = (h_factor * (1.0f32 - joint_state_sum / state_sum)) as f64;
        if switch_prob > 0.0 {
            self.sum_gen_dist += self.fpd.stage1_map.gen_dist[m] as f64;
            self.sum_switch_prob += switch_prob;
        }
        fwd_sum
    }
}

// ---------------------------------------------------------------------------
// PhaseLS stage-1 driver

/// Port of `PhaseLS.runStage1`.
pub fn run_stage1(pd: &mut PhaseData, par: &Par) {
    let timing2 = std::env::var("RUSTY_BEAGLE_TIMING2").is_ok();
    let t0 = std::time::Instant::now();
    let use_bwd = pd.it & 1 == 0;
    let coded_steps = CodedSteps::new(&pd.est_phase, par.nthreads);
    let t1 = std::time::Instant::now();
    let phase_ibs = PbwtPhaseIbs::new(pd, par, &coded_steps, use_bwd);
    let t2 = std::time::Instant::now();
    let all_haps = coded_steps.all_haps.clone();
    if par.em {
        let mut rand = JavaRandom::new(pd.it_seed());
        if pd.it == 0 {
            initialize_parameters(pd, par, &phase_ibs, &all_haps, &mut rand);
        } else if pd.it < par.burnin as usize {
            update_parameters(pd, par, &phase_ibs, &all_haps, &mut rand);
        }
    }
    let t3 = std::time::Instant::now();
    // phase every sample (work-stealing across threads; per-sample results
    // are order-independent because each sample's phase is stored separately)
    let n_samples = pd.fpd().n_targ_samples;
    let counter = std::sync::atomic::AtomicUsize::new(0);
    let pd_ref: &PhaseData = pd;
    let n_threads = par.nthreads;
    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            let counter = &counter;
            let phase_ibs = &phase_ibs;
            let all_haps: &XRefGT = &all_haps;
            scope.spawn(move || {
                let mut baum = PhaseBaum2::new(pd_ref, par, phase_ibs, all_haps);
                loop {
                    let s = counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if s >= n_samples {
                        break;
                    }
                    baum.phase(s);
                }
            });
        }
    });
    if timing2 {
        eprintln!(
            "[timing2] it {}: coded {:.3}s  ibs {:.3}s  em {:.3}s  baum {:.3}s",
            pd.it,
            (t1 - t0).as_secs_f64(),
            (t2 - t1).as_secs_f64(),
            (t3 - t2).as_secs_f64(),
            t3.elapsed().as_secs_f64()
        );
    }
}

fn initialize_parameters(
    pd: &mut PhaseData,
    par: &Par,
    phase_ibs: &PbwtPhaseIbs,
    all_haps: &Arc<XRefGT>,
    rand: &mut JavaRandom,
) {
    let mut prev_rec_int = pd.recomb_intensity;
    let max_initial_its = 15;
    for _ in 0..max_initial_its {
        update_parameters(pd, par, phase_ibs, all_haps, rand);
        let rec_int = pd.recomb_intensity;
        if (rec_int - prev_rec_int).abs() <= 0.1 * prev_rec_int {
            break;
        }
        prev_rec_int = rec_int;
    }
}

fn update_parameters(
    pd: &mut PhaseData,
    par: &Par,
    phase_ibs: &PbwtPhaseIbs,
    all_haps: &Arc<XRefGT>,
    rand: &mut JavaRandom,
) {
    let param_est = get_param_est(pd, par, phase_ibs, all_haps, rand);
    let prev_p_mismatch = pd.p_mismatch;
    let p_mismatch = param_est.p_mismatch();
    let recomb_intensity = param_est.recomb_intensity();
    if p_mismatch.is_finite() && p_mismatch > prev_p_mismatch {
        pd.update_p_mismatch(p_mismatch);
    }
    if recomb_intensity.is_finite() && recomb_intensity > 0.0 {
        pd.update_recomb_intensity(recomb_intensity);
    }
}

fn get_param_est(
    pd: &PhaseData,
    par: &Par,
    phase_ibs: &PbwtPhaseIbs,
    all_haps: &Arc<XRefGT>,
    rand: &mut JavaRandom,
) -> ParamEstimates {
    let param_est = ParamEstimates::default();
    let sample_indices = samples_to_analyze(pd, rand);
    let n_threads = par.nthreads.min(sample_indices.len());
    let max_sum = 20000.0 / n_threads as f64;
    let min_indices = 50usize;
    let counter = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..n_threads {
            let counter = &counter;
            let param_est = &param_est;
            let sample_indices = &sample_indices;
            let all_haps: &XRefGT = all_haps;
            scope.spawn(move || {
                let mut hpd = HmmParamData::new(pd, par, phase_ibs, all_haps);
                let mut index = counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                while (hpd.sum_switch_probs() < max_sum || index < min_indices)
                    && index < sample_indices.len()
                {
                    hpd.update(sample_indices[index] as usize);
                    index = counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    hpd.add_estimation_data(param_est);
                }
            });
        }
    });
    param_est
}

fn samples_to_analyze(pd: &PhaseData, rand: &mut JavaRandom) -> Vec<i32> {
    let max_samples_to_analyze = 500usize;
    let n_targ_samples = pd.fpd().n_targ_samples;
    let mut ia: Vec<i32> = (0..n_targ_samples as i32).collect();
    if n_targ_samples <= max_samples_to_analyze {
        ia
    } else {
        for j in 0..max_samples_to_analyze {
            let x = rand.next_int_bound((ia.len() - j) as i32) as usize;
            ia.swap(j, j + x);
        }
        ia.truncate(max_samples_to_analyze);
        ia
    }
}
