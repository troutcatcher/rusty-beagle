//! Ports of `phase.LowFreqPhaseStates`, `phase.HmmStateProbs`,
//! `phase.Stage2Haps`, `phase.Stage2Baum`, and the `phase.PhaseLS`
//! stage-2 driver.

use crate::codedsteps::CodedSteps;
use crate::javautil::{IntIntMap, JavaPriorityQueue, JavaRandom};
use crate::par::Par;
use crate::phasebaum::CompHapSegment;
use crate::phaseibs::LowFreqPhaseIbs;
use crate::phasedata::{Carriers, FixedPhaseData, PhaseData};
use crate::xref::XRefGT;
use std::sync::Mutex;

const NIL: i32 = -103;

// ---------------------------------------------------------------------------
// LowFreqPhaseStates

pub struct LowFreqPhaseStates<'a> {
    ibs_haps: &'a LowFreqPhaseIbs,
    fpd: &'a FixedPhaseData,
    all_haps: &'a XRefGT,
    #[allow(dead_code)]
    n_markers: usize,
    max_states: usize,
    min_steps: usize,
    it_seed: i64,
    hap_to_last_ibs_step: IntIntMap,
    q: JavaPriorityQueue<CompHapSegment>,
    comp_hap_hap: Vec<Vec<i32>>,
    comp_hap_end: Vec<Vec<i32>>,
    segment_index: Vec<usize>,
    comp_hap_to_hap: Vec<i32>,
    comp_hap_to_end: Vec<i32>,
}

impl<'a> LowFreqPhaseStates<'a> {
    pub fn new(
        pd: &PhaseData,
        ibs_haps: &'a LowFreqPhaseIbs,
        all_haps: &'a XRefGT,
        fpd: &'a FixedPhaseData,
        max_states: usize,
    ) -> LowFreqPhaseStates<'a> {
        let phase_step = fpd.ibs_step;
        let min_steps = std::cmp::max(200, (1.0f32 / phase_step).ceil() as usize);
        LowFreqPhaseStates {
            ibs_haps,
            fpd,
            all_haps,
            n_markers: all_haps.n_markers(),
            max_states,
            min_steps,
            it_seed: pd.it_seed(),
            hap_to_last_ibs_step: IntIntMap::with_capacity(max_states),
            q: JavaPriorityQueue::new(),
            comp_hap_hap: vec![Vec::new(); max_states],
            comp_hap_end: vec![Vec::new(); max_states],
            segment_index: vec![0; max_states],
            comp_hap_to_hap: vec![0; max_states],
            comp_hap_to_end: vec![0; max_states],
        }
    }

    /// `LowFreqPhaseStates.ibsStates(targHap, haps, nMismatches)`
    pub fn ibs_states(
        &mut self,
        targ_hap: usize,
        haps: &mut [i32],      // nMarkers x maxStates
        mismatches: &mut [u8], // nMarkers x maxStates
        row_len: usize,
    ) -> usize {
        let n_comp_haps = self.set_comp_ref_haps(targ_hap);
        self.copy_data(targ_hap, n_comp_haps, haps, mismatches, row_len);
        n_comp_haps
    }

    fn set_comp_ref_haps(&mut self, targ_hap: usize) -> usize {
        self.q.clear();
        self.hap_to_last_ibs_step.clear();
        for j in 0..self.max_states {
            self.comp_hap_hap[j].clear();
            self.comp_hap_end[j].clear();
        }
        for step in 0..self.fpd.stage1_steps.size() {
            self.add_ibs_hap(self.ibs_haps.fwd.ibs_hap(targ_hap, step), step as i32);
            self.add_ibs_hap(self.ibs_haps.bwd.ibs_hap(targ_hap, step), step as i32);
        }
        if self.q.is_empty() {
            self.fill_q_with_random_haps(targ_hap);
        }
        self.set_final_ref_segs()
    }

    fn add_ibs_hap(&mut self, ibs_hap: i32, step: i32) {
        if ibs_hap < 0 {
            return;
        }
        if self.hap_to_last_ibs_step.get(ibs_hap, NIL) == NIL {
            self.update_head_of_q();
            let recycle = self.q.len() == self.max_states
                || (!self.q.is_empty()
                    && step - self.q.peek().unwrap().last_ibs_step >= self.min_steps as i32);
            if recycle {
                let head = self.q.poll().unwrap();
                let index = head.comp_hap_index;
                let prev_hap = head.hap;
                let next_start = self
                    .fpd
                    .stage1_steps
                    .start((((head.last_ibs_step + step) as u32) >> 1) as usize);
                self.hap_to_last_ibs_step.remove(prev_hap);
                self.comp_hap_hap[index].push(ibs_hap); // hap of new segment
                self.comp_hap_end[index].push(next_start as i32); // end of old segment
                self.q.offer(CompHapSegment {
                    hap: ibs_hap,
                    start_marker: next_start,
                    last_ibs_step: step,
                    comp_hap_index: index,
                });
            } else {
                let index = self.q.len();
                self.comp_hap_hap[index].push(ibs_hap); // hap of new segment
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

    fn set_final_ref_segs(&mut self) -> usize {
        let n_comp_haps = self.q.len();
        while let Some(head) = self.q.poll() {
            let comp_hap = head.comp_hap_index;
            self.comp_hap_end[comp_hap].push(self.n_markers as i32); // add missing end of last segment
            self.segment_index[comp_hap] = 0;
            self.comp_hap_to_hap[comp_hap] = self.comp_hap_hap[comp_hap][0];
            self.comp_hap_to_end[comp_hap] = self.comp_hap_end[comp_hap][0];
        }
        n_comp_haps
    }

    fn copy_data(
        &mut self,
        targ_hap: usize,
        n_comp_haps: usize,
        haps: &mut [i32],
        mismatches: &mut [u8],
        row_len: usize,
    ) {
        for m in 0..self.n_markers {
            let obs_allele = self.all_haps.allele(m, targ_hap);
            let row = m * row_len;
            for j in 0..n_comp_haps {
                if m as i32 == self.comp_hap_to_end[j] {
                    self.segment_index[j] += 1;
                    self.comp_hap_to_hap[j] = self.comp_hap_hap[j][self.segment_index[j]];
                    self.comp_hap_to_end[j] = self.comp_hap_end[j][self.segment_index[j]];
                }
                let ref_hap = self.comp_hap_to_hap[j];
                haps[row + j] = ref_hap;
                mismatches[row + j] = if self.all_haps.allele(m, ref_hap as usize) == obs_allele
                {
                    0
                } else {
                    1
                };
            }
        }
    }

    fn fill_q_with_random_haps(&mut self, hap: usize) {
        debug_assert!(self.q.is_empty());
        let n_haps = self.all_haps.n_haps();
        let n_states = std::cmp::min(n_haps - 2, self.max_states);
        if n_states == 0 {
            eprintln!("ERROR: there is only one sample");
            std::process::exit(1);
        }
        let mut rand = JavaRandom::new(self.it_seed.wrapping_add(hap as i64));
        let sample = hap >> 1;
        let ibs_step = 0i32;
        let start_marker = 0usize;
        for j in 0..n_states {
            let mut h = rand.next_int_bound(n_haps as i32);
            while (h >> 1) as usize == sample {
                h = rand.next_int_bound(n_haps as i32);
            }
            self.comp_hap_hap[self.q.len()].push(h);
            self.q.offer(CompHapSegment {
                hap: h,
                start_marker,
                last_ibs_step: ibs_step,
                comp_hap_index: j,
            });
        }
    }
}

// ---------------------------------------------------------------------------
// HmmStateProbs

pub struct HmmStateProbs<'a> {
    states: LowFreqPhaseStates<'a>,
    p_recomb: &'a [f32],
    mismatch: Vec<u8>, // nMarkers x maxStates
    bwd: Vec<f32>,
    p_mismatch: [f32; 2],
    n_markers: usize,
    max_states: usize,
}

impl<'a> HmmStateProbs<'a> {
    pub fn new(
        pd: &'a PhaseData,
        par: &Par,
        ibs: &'a LowFreqPhaseIbs,
        all_haps: &'a XRefGT,
    ) -> HmmStateProbs<'a> {
        let fpd = pd.fpd();
        let n_markers = fpd.n_stage1_markers();
        let max_states = par.phase_states as usize / 2;
        let p_miss = pd.p_mismatch;
        HmmStateProbs {
            states: LowFreqPhaseStates::new(pd, ibs, all_haps, fpd, max_states),
            p_recomb: &pd.p_recomb,
            mismatch: vec![0u8; n_markers * max_states],
            bwd: vec![0f32; max_states],
            p_mismatch: [1.0 - p_miss, p_miss],
            n_markers,
            max_states,
        }
    }

    pub fn max_states(&self) -> usize {
        self.max_states
    }

    /// `HmmStateProbs.run(targHap, refHaps, stateProbs)`; both output arrays
    /// are (nMarkers x maxStates) flattened.
    pub fn run(
        &mut self,
        targ_hap: usize,
        ref_haps: &mut [i32],
        state_probs: &mut [f32],
    ) -> usize {
        let n_states =
            self.states
                .ibs_states(targ_hap, ref_haps, &mut self.mismatch, self.max_states);
        self.run_fwd(state_probs, n_states);
        self.run_bwd(state_probs, n_states);
        n_states
    }

    fn run_fwd(&mut self, probs: &mut [f32], n_states: usize) {
        let mut last_sum = 0.0f32;
        for j in 0..n_states {
            probs[j] = self.p_mismatch[self.mismatch[j] as usize];
            last_sum += probs[j];
        }
        for m in 1..self.n_markers {
            let m_m1_row = (m - 1) * self.max_states;
            let row = m * self.max_states;
            let p_rec = self.p_recomb[m];
            let shift = p_rec / n_states as f32;
            let scale = (1.0f32 - p_rec) / last_sum;
            last_sum = 0.0f32;
            for j in 0..n_states {
                let em = self.p_mismatch[self.mismatch[row + j] as usize];
                probs[row + j] = em * (scale * probs[m_m1_row + j] + shift);
                last_sum += probs[row + j];
            }
        }
    }

    fn run_bwd(&mut self, probs: &mut [f32], n_states: usize) {
        let incl_end = self.n_markers - 1;
        self.bwd[..n_states].fill(1.0 / n_states as f32);
        for m in (0..incl_end).rev() {
            let m_p1 = m + 1;
            let row_p1 = m_p1 * self.max_states;
            let mut sum = 0.0f32;
            for j in 0..n_states {
                self.bwd[j] *= self.p_mismatch[self.mismatch[row_p1 + j] as usize];
                sum += self.bwd[j];
            }
            let p_rec = self.p_recomb[m_p1];
            let scale = (1.0f32 - p_rec) / sum;
            let shift = p_rec / n_states as f32;
            let row = m * self.max_states;
            let mut sum2 = 0.0f32;
            for j in 0..n_states {
                self.bwd[j] = scale * self.bwd[j] + shift;
                probs[row + j] *= self.bwd[j];
                sum2 += probs[row + j];
            }
            for j in 0..n_states {
                probs[row + j] /= sum2;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Stage2Haps

/// Thread-safe store of stage-2 phased genotypes. Phased rare-allele
/// carriers are recorded per (marker, allele) and sorted on read-out,
/// making the result independent of thread scheduling (as in Java).
pub struct Stage2Haps {
    #[allow(dead_code)]
    pub n_targ_samples: usize,
    /// rare carrier hap lists indexed by sumAlleles(marker)+allele;
    /// None = high-frequency allele
    rare_carriers: Vec<Option<Mutex<Vec<i32>>>>,
}

impl Stage2Haps {
    pub fn new(fpd: &FixedPhaseData) -> Stage2Haps {
        let size = *fpd.sum_alleles.last().unwrap() as usize;
        let mut rare_carriers: Vec<Option<Mutex<Vec<i32>>>> = Vec::with_capacity(size);
        for m in 0..fpd.targ.len() {
            let n_alleles = fpd.targ_markers[m].n_alleles as usize;
            for al in 0..n_alleles {
                if matches!(fpd.carriers[m][al], Carriers::HighFreq) {
                    rare_carriers.push(None);
                } else {
                    rare_carriers.push(Some(Mutex::new(Vec::new())));
                }
            }
        }
        Stage2Haps {
            n_targ_samples: fpd.n_targ_samples,
            rare_carriers,
        }
    }

    /// `Stage2Haps.setPhasedGT(marker, sample, a1, a2)`
    pub fn set_phased_gt(
        &self,
        fpd: &FixedPhaseData,
        marker: usize,
        sample: usize,
        a1: i32,
        a2: i32,
    ) {
        let offset = fpd.sum_alleles[marker] as usize;
        if let Some(list) = &self.rare_carriers[offset + a1 as usize] {
            list.lock().unwrap().push((sample << 1) as i32);
        }
        if let Some(list) = &self.rare_carriers[offset + a2 as usize] {
            list.lock().unwrap().push(((sample << 1) | 1) as i32);
        }
    }

    /// Returns the phased alleles for output marker `m`
    /// (`Stage2Haps.gtRec`): stage-1 markers come from the stage-1 phase,
    /// others from the sorted rare-carrier lists.
    pub fn alleles_at(
        &self,
        fpd: &FixedPhaseData,
        stage1_alleles: &dyn Fn(usize) -> Vec<i16>,
        m: usize,
    ) -> Vec<i16> {
        let m1 = fpd.prev_stage1_marker[m];
        let m2 = fpd.stage1_to2[m1];
        if m == m2 {
            stage1_alleles(m1)
        } else {
            self.stage2_alleles(fpd, m)
        }
    }

    fn stage2_alleles(&self, fpd: &FixedPhaseData, m: usize) -> Vec<i16> {
        let n_haps = fpd.n_targ_haps;
        let al_start = fpd.sum_alleles[m] as usize;
        let al_end = fpd.sum_alleles[m + 1] as usize;
        // determine major allele (no carrier list) and collect sorted lists
        let mut major_allele: isize = -1;
        let mut lists: Vec<Option<Vec<i32>>> = Vec::with_capacity(al_end - al_start);
        for j in 0..(al_end - al_start) {
            match &self.rare_carriers[al_start + j] {
                Some(list) => {
                    let mut v = list.lock().unwrap().clone();
                    v.sort_unstable();
                    lists.push(Some(v));
                }
                None => {
                    major_allele = j as isize;
                    lists.push(None);
                }
            }
        }
        if major_allele == -1 {
            // can occur if all alleles are rare due to high missing rate
            let mut maj = 0usize;
            for j in 1..lists.len() {
                if lists[j].as_ref().unwrap().len() > lists[maj].as_ref().unwrap().len() {
                    maj = j;
                }
            }
            major_allele = maj as isize;
            lists[maj] = None;
        }
        let mut alleles = vec![major_allele as i16; n_haps];
        for (j, list) in lists.iter().enumerate() {
            if let Some(v) = list {
                for &h in v {
                    alleles[h as usize] = j as i16;
                }
            }
        }
        alleles
    }
}

// ---------------------------------------------------------------------------
// Stage2Baum

pub struct Stage2Baum<'a> {
    fpd: &'a FixedPhaseData,
    pd: &'a PhaseData,
    state_probs: HmmStateProbs<'a>,
    n_states: [usize; 2],
    states: [Vec<i32>; 2],
    probs: [Vec<f32>; 2],
    n_targ_haps: usize,
    n_stage1_markers: usize,
    stage2_haps: &'a Stage2Haps,
    rand: JavaRandom,
}

impl<'a> Stage2Baum<'a> {
    pub fn new(
        pd: &'a PhaseData,
        par: &Par,
        ibs: &'a LowFreqPhaseIbs,
        all_haps: &'a XRefGT,
        stage2_haps: &'a Stage2Haps,
    ) -> Stage2Baum<'a> {
        let fpd = pd.fpd();
        let n_stage1_markers = fpd.n_stage1_markers();
        let state_probs = HmmStateProbs::new(pd, par, ibs, all_haps);
        let max_states = state_probs.max_states();
        Stage2Baum {
            fpd,
            pd,
            state_probs,
            n_states: [0, 0],
            states: [
                vec![0i32; n_stage1_markers * max_states],
                vec![0i32; n_stage1_markers * max_states],
            ],
            probs: [
                vec![0f32; n_stage1_markers * max_states],
                vec![0f32; n_stage1_markers * max_states],
            ],
            n_targ_haps: fpd.n_targ_haps,
            n_stage1_markers,
            stage2_haps,
            rand: JavaRandom::new(0),
        }
    }

    /// `Stage2Baum.phase(targSample)`
    pub fn phase(&mut self, targ_sample: usize) {
        self.rand
            .set_seed(self.pd.it_seed().wrapping_add(targ_sample as i64));
        let h1 = targ_sample << 1;
        let h2 = h1 | 1;
        let (s0, s1) = {
            let max_states = self.state_probs.max_states();
            let _ = max_states;
            let n0 = self
                .state_probs
                .run(h1, &mut self.states[0], &mut self.probs[0]);
            let n1 = self
                .state_probs
                .run(h2, &mut self.states[1], &mut self.probs[1]);
            (n0, n1)
        };
        self.n_states = [s0, s1];

        let mut start = 0usize;
        for j in 0..self.n_stage1_markers {
            let end = self.fpd.stage1_to2[j];
            self.impute_interval(targ_sample, start, end);
            start = end + 1;
        }
        self.impute_interval(targ_sample, start, self.fpd.targ.len());
    }

    fn impute_interval(&mut self, sample: usize, start: usize, end: usize) {
        let hap1 = sample << 1;
        let hap2 = hap1 | 1;
        for m in start..end {
            let mut a1 = self.fpd.targ[m][hap1] as i32;
            let mut a2 = self.fpd.targ[m][hap2] as i32;
            if a1 >= 0 && a2 >= 0 {
                if a1 != a2 {
                    let al_probs1 = self.unscaled_al_probs(m, 0, a1, a2);
                    let al_probs2 = self.unscaled_al_probs(m, 1, a1, a2);
                    let p1 = al_probs1[a1 as usize] * al_probs2[a2 as usize];
                    let p2 = al_probs1[a2 as usize] * al_probs2[a1 as usize];
                    let switch_alleles = p1 < p2 || (p1 == p2 && self.rand.next_boolean());
                    if switch_alleles {
                        std::mem::swap(&mut a1, &mut a2);
                    }
                }
            } else {
                a1 = self.impute_allele(m, 0);
                a2 = self.impute_allele(m, 1);
            }
            self.stage2_haps
                .set_phased_gt(self.fpd, m, sample, a1, a2);
        }
    }

    fn unscaled_al_probs(&self, m: usize, hap_bit: usize, a1: i32, a2: i32) -> Vec<f32> {
        let n_alleles = self.fpd.targ_markers[m].n_alleles as usize;
        let mut al_probs = vec![0.0f32; n_alleles];
        let rare1 = self.fpd.is_low_freq(m, a1 as usize);
        let rare2 = self.fpd.is_low_freq(m, a2 as usize);
        let mkr_a = self.fpd.prev_stage1_marker[m];
        let mkr_b = (mkr_a + 1).min(self.n_stage1_markers - 1);
        let max_states = self.state_probs.max_states();
        let row_a = mkr_a * max_states;
        let row_b = mkr_b * max_states;
        let states_a = &self.states[hap_bit][row_a..];
        let probs_a = &self.probs[hap_bit][row_a..];
        let probs_b = &self.probs[hap_bit][row_b..];
        let wt = self.fpd.prev_stage1_wt[m];
        for j in 0..self.n_states[hap_bit] {
            let hap = states_a[j];
            let b1 = self.allele(m, hap as usize);
            let b2 = self.allele(m, (hap ^ 1) as usize);
            if b1 >= 0 && b2 >= 0 {
                let prob = wt * probs_a[j] + (1.0f32 - wt) * probs_b[j];
                if b1 == b2 {
                    al_probs[b1 as usize] += prob;
                } else {
                    let match1 = rare1 && (a1 == b1 || a1 == b2);
                    let match2 = rare2 && (a2 == b1 || a2 == b2);
                    if match1 ^ match2 {
                        if match1 {
                            al_probs[a1 as usize] += prob;
                        } else {
                            al_probs[a2 as usize] += prob;
                        }
                    }
                }
            }
        }
        al_probs
    }

    fn impute_allele(&self, m: usize, hap_bit: usize) -> i32 {
        let n_alleles = self.fpd.targ_markers[m].n_alleles as usize;
        let mut al_probs = vec![0.0f32; n_alleles];
        let mkr_a = self.fpd.prev_stage1_marker[m];
        let mkr_b = (mkr_a + 1).min(self.n_stage1_markers - 1);
        let max_states = self.state_probs.max_states();
        let row_a = mkr_a * max_states;
        let row_b = mkr_b * max_states;
        let states_a = &self.states[hap_bit][row_a..];
        let probs_a = &self.probs[hap_bit][row_a..];
        let probs_b = &self.probs[hap_bit][row_b..];
        for j in 0..self.n_states[hap_bit] {
            let wt = self.fpd.prev_stage1_wt[m];
            let prob = wt * probs_a[j] + (1.0f32 - wt) * probs_b[j];
            let hap = states_a[j];
            let b1 = self.allele(m, hap as usize);
            let b2 = self.allele(m, (hap ^ 1) as usize);
            if b1 >= 0 && b2 >= 0 {
                if b1 == b2 || hap as usize >= self.n_targ_haps {
                    al_probs[b1 as usize] += prob;
                } else {
                    let is_rare1 = self.fpd.is_low_freq(m, b1 as usize);
                    let is_rare2 = self.fpd.is_low_freq(m, b2 as usize);
                    // Java: float += double  =>  add in f64, narrow to f32
                    if is_rare1 ^ is_rare2 {
                        let (w1, w2) = if is_rare1 { (0.55f64, 0.45f64) } else { (0.45f64, 0.55f64) };
                        let s1 = &mut al_probs[b1 as usize];
                        *s1 = (*s1 as f64 + w1 * prob as f64) as f32;
                        let s2 = &mut al_probs[b2 as usize];
                        *s2 = (*s2 as f64 + w2 * prob as f64) as f32;
                    } else {
                        let s1 = &mut al_probs[b1 as usize];
                        *s1 = (*s1 as f64 + 0.5f64 * prob as f64) as f32;
                        let s2 = &mut al_probs[b2 as usize];
                        *s2 = (*s2 as f64 + 0.5f64 * prob as f64) as f32;
                    }
                }
            }
        }
        let mut max_index = 0usize;
        for j in 1..al_probs.len() {
            if al_probs[j] > al_probs[max_index] {
                max_index = j;
            }
        }
        max_index as i32
    }

    /// target-then-reference allele lookup over ALL window markers
    fn allele(&self, marker: usize, hap: usize) -> i32 {
        if hap < self.n_targ_haps {
            self.fpd.targ[marker][hap] as i32
        } else {
            self.fpd.restrict_ref[marker].allele(hap - self.n_targ_haps) as i32
        }
    }
}

// ---------------------------------------------------------------------------
// PhaseLS.runStage2

/// Port of `PhaseLS.runStage2`.
pub fn run_stage2(pd: &PhaseData, par: &Par) -> (Stage2Haps, Vec<Vec<i16>>) {
    let fpd = pd.fpd();
    let coded_steps = CodedSteps::new(&pd.est_phase, par.nthreads);
    let phase_ibs = LowFreqPhaseIbs::new(pd, par, &coded_steps);
    let all_haps = coded_steps.all_haps.clone();
    let stage2_haps = Stage2Haps::new(fpd);
    let n_samples = fpd.n_targ_samples;
    let counter = std::sync::atomic::AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..par.nthreads {
            let counter = &counter;
            let phase_ibs = &phase_ibs;
            let all_haps: &XRefGT = &all_haps;
            let stage2_haps = &stage2_haps;
            scope.spawn(move || {
                let mut baum = Stage2Baum::new(pd, par, phase_ibs, all_haps, stage2_haps);
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
    // stage-1 phased alleles, marker-major over stage1 markers
    let stage1_alleles = stage1_marker_alleles(pd);
    (stage2_haps, stage1_alleles)
}

/// marker-major stage-1 phased alleles from the current EstPhase
pub fn stage1_marker_alleles(pd: &PhaseData) -> Vec<Vec<i16>> {
    let fpd = pd.fpd();
    let n_markers = fpd.n_stage1_markers();
    let n_haps = fpd.n_targ_haps;
    let layout = &fpd.stage1_layout;
    let mut result = vec![vec![0i16; n_haps]; n_markers];
    for s in 0..fpd.n_targ_samples {
        pd.est_phase.with(s, |sp| {
            let h1 = s << 1;
            for (m, row) in result.iter_mut().enumerate() {
                row[h1] = layout.allele(&sp.hap1, m) as i16;
                row[h1 | 1] = layout.allele(&sp.hap2, m) as i16;
            }
        });
    }
    result
}
