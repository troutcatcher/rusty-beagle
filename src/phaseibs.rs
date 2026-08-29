//! Ports of `phase.PbwtIbsData`, `phase.PbwtPhaseIbs`,
//! `phase.LowFreqPbwtPhaseIbs`, and `phase.LowFreqPhaseIbs`.

use crate::codedsteps::CodedSteps;
use crate::javautil::JavaRandom;
use crate::par::Par;
use crate::pbwt::PbwtDivUpdater;
use crate::phasedata::{Carriers, FixedPhaseData, PhaseData};
use rayon::prelude::*;

const BURNIN_CANDIDATES: usize = 100;
const MAX_PHASE_CANDIDATES: usize = 90;
const MIN_PHASE_CANDIDATES: usize = 5;
const STAGE2_CANDIDATES: usize = 10;
const MAX_BACKOFF_CM: f32 = 0.3;

/// Port of `phase.PbwtIbsData` (parameters shared by the IBS selectors).
#[allow(dead_code)]
pub struct PbwtIbsData {
    pub n_haps: usize,
    pub n_targ_haps: usize,
    pub n_candidates: usize,
    pub n_steps: usize,
    pub n_overlap_steps: usize,
    pub max_backoff_steps: usize,
    pub steps_per_batch: usize,
    pub n_batches: usize,
}

impl PbwtIbsData {
    pub fn new(pd: &PhaseData, par: &Par) -> PbwtIbsData {
        let fpd = pd.fpd();
        let n_threads = par.nthreads;
        let n_its = (par.burnin + par.iterations) as usize;
        let n_candidates = if pd.it < n_its {
            n_candidates1(pd, par)
        } else {
            STAGE2_CANDIDATES.min(fpd.n_haps)
        };
        let n_steps = fpd.stage1_steps.size();
        let n_overlap_steps =
            java_rint_f32(par.buffer / fpd.ibs_step) as usize;
        let max_backoff_steps = java_rint_f32(MAX_BACKOFF_CM / fpd.ibs_step) as usize;
        let steps_per_batch = (n_steps + n_threads - 1) / n_threads;
        let n_batches = (n_steps + steps_per_batch - 1) / steps_per_batch;
        PbwtIbsData {
            n_haps: fpd.n_haps,
            n_targ_haps: fpd.n_targ_haps,
            n_candidates,
            n_steps,
            n_overlap_steps,
            max_backoff_steps,
            steps_per_batch,
            n_batches,
        }
    }

    pub fn start_step(&self, batch: usize) -> usize {
        batch * self.steps_per_batch
    }

    pub fn end_step(&self, batch: usize) -> usize {
        ((batch + 1) * self.steps_per_batch).min(self.n_steps)
    }

    pub fn buffer_start_step(&self, start_step: usize) -> usize {
        start_step.saturating_sub(self.n_overlap_steps)
    }

    pub fn buffer_end_step(&self, end_step: usize) -> usize {
        (end_step + self.n_overlap_steps).min(self.n_steps)
    }
}

/// `Math.rint(float)` after f32->f64 promotion.
fn java_rint_f32(x: f32) -> f64 {
    (x as f64).round_ties_even()
}

fn n_candidates1(pd: &PhaseData, par: &Par) -> usize {
    let mut n_candidates = BURNIN_CANDIDATES;
    let it = pd.it;
    if it >= par.burnin as usize {
        let n_its_remaining = (par.burnin + par.iterations) as f64 - it as f64;
        let p = n_its_remaining / par.iterations as f64;
        n_candidates = (p * MAX_PHASE_CANDIDATES as f64).round() as usize;
        n_candidates = n_candidates.max(MIN_PHASE_CANDIDATES);
    }
    n_candidates.min(pd.fpd().n_haps)
}

/// Port of `phase.PbwtPhaseIbs` (stage-1 IBS selection).
pub struct PbwtPhaseIbs {
    pub ibs_haps: Vec<Vec<i32>>, // [step][targ hap]
}

impl PbwtPhaseIbs {
    pub fn new(
        pd: &PhaseData,
        par: &Par,
        coded_steps: &CodedSteps,
        use_bwd: bool,
    ) -> PbwtPhaseIbs {
        let data = PbwtIbsData::new(pd, par);
        let fpd = pd.fpd();
        let ibs_haps: Vec<Vec<i32>> = (0..data.n_batches)
            .into_par_iter()
            .flat_map_iter(|batch| {
                if use_bwd {
                    bwd_ibs_haps(pd, fpd, coded_steps, &data, batch)
                } else {
                    fwd_ibs_haps(pd, fpd, coded_steps, &data, batch)
                }
            })
            .collect();
        PbwtPhaseIbs { ibs_haps }
    }

    #[inline]
    pub fn ibs_hap(&self, hap: usize, step: usize) -> i32 {
        self.ibs_haps[step][hap]
    }
}

fn bwd_ibs_haps(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    coded_steps: &CodedSteps,
    data: &PbwtIbsData,
    batch: usize,
) -> Vec<Vec<i32>> {
    let start_step = data.start_step(batch);
    let end_step = data.end_step(batch);
    let buffer_end_step = data.buffer_end_step(end_step);
    let n_haps = data.n_haps;
    let mut pbwt = PbwtDivUpdater::new(n_haps);
    let mut a: Vec<i32> = (0..n_haps as i32).collect();
    let mut d: Vec<i32> = vec![(buffer_end_step as i32) - 1; n_haps + 1];
    let mut out: Vec<Vec<i32>> = vec![Vec::new(); end_step - start_step];
    for j in (end_step..buffer_end_step).rev() {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.bwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
    }
    for j in (start_step..end_step).rev() {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.bwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
        out[j - start_step] = get_bwd_ibs_haps(pd, fpd, j, &a, &mut d, data);
    }
    out
}

fn fwd_ibs_haps(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    coded_steps: &CodedSteps,
    data: &PbwtIbsData,
    batch: usize,
) -> Vec<Vec<i32>> {
    let start_step = data.start_step(batch);
    let end_step = data.end_step(batch);
    let buffer_start_step = data.buffer_start_step(start_step);
    let n_haps = data.n_haps;
    let mut pbwt = PbwtDivUpdater::new(n_haps);
    let mut a: Vec<i32> = (0..n_haps as i32).collect();
    let mut d: Vec<i32> = vec![buffer_start_step as i32; n_haps + 1];
    let mut out: Vec<Vec<i32>> = vec![Vec::new(); end_step - start_step];
    for j in buffer_start_step..start_step {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.fwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
    }
    for j in start_step..end_step {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.fwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
        out[j - start_step] = get_fwd_ibs_haps(pd, fpd, j, &a, &mut d, data);
    }
    out
}

fn get_bwd_ibs_haps(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    step: usize,
    a: &[i32],
    d: &mut [i32],
    data: &PbwtIbsData,
) -> Vec<i32> {
    let mut rand = JavaRandom::new(pd.it_seed().wrapping_add(step as i64));
    let m_start = fpd.stage1_steps.start(step);
    let m_incl_end = fpd.stage1_steps.end(step) - 1;
    let mut selected = vec![0i32; data.n_targ_haps];
    let ibs2 = &fpd.stage1_ibs2;
    let step_i = step as i32;
    d[0] = step_i - 2;
    d[a.len()] = step_i - 2;
    for i in 0..a.len() {
        if (a[i] as usize) < data.n_targ_haps {
            let hap = a[i];
            let s1 = (hap >> 1) as usize;
            let mut u = i; // inclusive start
            let mut v = i + 1; // exclusive end
            let mut u_next_match_end = d[u];
            let mut v_next_match_end = d[v];
            while v - u < data.n_candidates
                && (step_i <= u_next_match_end || step_i <= v_next_match_end)
            {
                if u_next_match_end <= v_next_match_end {
                    v += 1;
                    v_next_match_end = d[v].min(v_next_match_end);
                } else {
                    u -= 1;
                    u_next_match_end = d[u].min(u_next_match_end);
                }
            }
            let n = v - u;
            selected[hap as usize] = -1;
            if n > 1 {
                let mut index = u + rand.next_int_bound(n as i32) as usize;
                for _ in 0..n {
                    if index == v {
                        index = u;
                    }
                    if index != i
                        && !ibs2.are_ibs2_range(
                            s1,
                            (a[index] >> 1) as usize,
                            m_start,
                            m_incl_end,
                        )
                    {
                        selected[hap as usize] = a[index];
                        break;
                    }
                    index += 1;
                }
            }
        }
    }
    selected
}

fn get_fwd_ibs_haps(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    step: usize,
    a: &[i32],
    d: &mut [i32],
    data: &PbwtIbsData,
) -> Vec<i32> {
    let mut rand = JavaRandom::new(pd.it_seed().wrapping_add(step as i64));
    let m_start = fpd.stage1_steps.start(step);
    let m_incl_end = fpd.stage1_steps.end(step) - 1;
    let mut selected = vec![0i32; data.n_targ_haps];
    let ibs2 = &fpd.stage1_ibs2;
    let step_i = step as i32;
    d[0] = step_i + 2;
    d[a.len()] = step_i + 2;
    for i in 0..a.len() {
        if (a[i] as usize) < data.n_targ_haps {
            let hap = a[i];
            let s1 = (hap >> 1) as usize;
            let mut u = i;
            let mut v = i + 1;
            let mut u_next_match_start = d[u];
            let mut v_next_match_start = d[v];
            while v - u < data.n_candidates
                && (u_next_match_start <= step_i || v_next_match_start <= step_i)
            {
                if v_next_match_start <= u_next_match_start {
                    v += 1;
                    v_next_match_start = d[v].max(v_next_match_start);
                } else {
                    u -= 1;
                    u_next_match_start = d[u].max(u_next_match_start);
                }
            }
            let n = v - u;
            selected[hap as usize] = -1;
            if n > 1 {
                let mut index = u + rand.next_int_bound(n as i32) as usize;
                for _ in 0..n {
                    if index == v {
                        index = u;
                    }
                    if index != i
                        && !ibs2.are_ibs2_range(
                            s1,
                            (a[index] >> 1) as usize,
                            m_start,
                            m_incl_end,
                        )
                    {
                        selected[hap as usize] = a[index];
                        break;
                    }
                    index += 1;
                }
            }
        }
    }
    selected
}

// ---------------------------------------------------------------------------
// LowFreqPbwtPhaseIbs / LowFreqPhaseIbs (stage-2 IBS selection)

pub struct LowFreqPbwtPhaseIbs {
    pub ibs_haps: Vec<Vec<i32>>, // [step][targ hap]
}

impl LowFreqPbwtPhaseIbs {
    pub fn new(
        pd: &PhaseData,
        par: &Par,
        coded_steps: &CodedSteps,
        use_bwd: bool,
    ) -> LowFreqPbwtPhaseIbs {
        let data = PbwtIbsData::new(pd, par);
        let fpd = pd.fpd();
        let ibs_haps: Vec<Vec<i32>> = (0..data.n_batches)
            .into_par_iter()
            .flat_map_iter(|batch| {
                if use_bwd {
                    lf_bwd_ibs_haps(pd, fpd, coded_steps, &data, batch)
                } else {
                    lf_fwd_ibs_haps(pd, fpd, coded_steps, &data, batch)
                }
            })
            .collect();
        LowFreqPbwtPhaseIbs { ibs_haps }
    }

    #[inline]
    pub fn ibs_hap(&self, hap: usize, step: usize) -> i32 {
        self.ibs_haps[step][hap]
    }
}

fn lf_bwd_ibs_haps(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    coded_steps: &CodedSteps,
    data: &PbwtIbsData,
    batch: usize,
) -> Vec<Vec<i32>> {
    let start_step = data.start_step(batch);
    let end_step = data.end_step(batch);
    let buffer_end_step = data.buffer_end_step(end_step);
    let n_haps = data.n_haps;
    let mut pbwt = PbwtDivUpdater::new(n_haps);
    let mut a: Vec<i32> = (0..n_haps as i32).collect();
    let mut d: Vec<i32> = vec![(buffer_end_step as i32) - 1; n_haps + 1];
    let mut a_inv = vec![0i32; n_haps];
    let mut i_to_prev_i = vec![0i32; n_haps];
    let mut i_to_next_i = vec![0i32; n_haps];
    let mut out: Vec<Vec<i32>> = vec![Vec::new(); end_step - start_step];
    for j in (end_step..buffer_end_step).rev() {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.bwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
    }
    for j in (start_step..end_step).rev() {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.bwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
        set_inv(&a, &mut a_inv);
        set_i_to_prev_next_i(fpd, j, &a_inv, &mut i_to_prev_i, &mut i_to_next_i);
        out[j - start_step] =
            lf_get_bwd(pd, fpd, j, &a, &mut d, &i_to_prev_i, &i_to_next_i, data);
    }
    out
}

fn lf_fwd_ibs_haps(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    coded_steps: &CodedSteps,
    data: &PbwtIbsData,
    batch: usize,
) -> Vec<Vec<i32>> {
    let start_step = data.start_step(batch);
    let end_step = data.end_step(batch);
    let buffer_start_step = data.buffer_start_step(start_step);
    let n_haps = data.n_haps;
    let mut pbwt = PbwtDivUpdater::new(n_haps);
    let mut a: Vec<i32> = (0..n_haps as i32).collect();
    let mut d: Vec<i32> = vec![buffer_start_step as i32; n_haps + 1];
    let mut a_inv = vec![0i32; n_haps];
    let mut i_to_prev_i = vec![0i32; n_haps];
    let mut i_to_next_i = vec![0i32; n_haps];
    let mut out: Vec<Vec<i32>> = vec![Vec::new(); end_step - start_step];
    for j in buffer_start_step..start_step {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.fwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
    }
    for j in start_step..end_step {
        let (h2s, vs) = coded_steps.get(j);
        pbwt.fwd_update(|h| h2s[h], *vs as usize, j as i32, &mut a, &mut d);
        set_inv(&a, &mut a_inv);
        set_i_to_prev_next_i(fpd, j, &a_inv, &mut i_to_prev_i, &mut i_to_next_i);
        out[j - start_step] =
            lf_get_fwd(pd, fpd, j, &a, &mut d, &i_to_prev_i, &i_to_next_i, data);
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn lf_get_bwd(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    step: usize,
    a: &[i32],
    d: &mut [i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> Vec<i32> {
    let mut rand = JavaRandom::new(pd.it_seed().wrapping_add(step as i64));
    let m_start = fpd.stage1_steps.start(step);
    let m_incl_end = fpd.stage1_steps.end(step) - 1;
    let mut selected = vec![0i32; data.n_targ_haps];
    let step_i = step as i32;
    d[a.len()] = step_i - 1; // set sentinel
    for i in 0..a.len() {
        if (a[i] as usize) < data.n_targ_haps {
            let best_i = lf_best_bwd(fpd, step_i, m_start, m_incl_end, i, a, d, i_to_prev_i, i_to_next_i, data);
            if best_i >= 0 {
                selected[a[i] as usize] = a[best_i as usize];
            } else {
                let mut u = i;
                let mut v = i + 1;
                let mut u_next_match_end = d[u];
                let mut v_next_match_end = d[v];
                while v - u < data.n_candidates
                    && (step_i <= u_next_match_end || step_i <= v_next_match_end)
                {
                    if u_next_match_end <= v_next_match_end {
                        v += 1;
                        v_next_match_end = d[v].min(v_next_match_end);
                    } else {
                        u -= 1;
                        u_next_match_end = d[u].min(u_next_match_end);
                    }
                }
                selected[a[i] as usize] =
                    lf_get_match(fpd, m_start, m_incl_end, i, u, v, a, &mut rand);
            }
        }
    }
    selected
}

#[allow(clippy::too_many_arguments)]
fn lf_get_fwd(
    pd: &PhaseData,
    fpd: &FixedPhaseData,
    step: usize,
    a: &[i32],
    d: &mut [i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> Vec<i32> {
    let mut rand = JavaRandom::new(pd.it_seed().wrapping_add(step as i64));
    let m_start = fpd.stage1_steps.start(step);
    let m_incl_end = fpd.stage1_steps.end(step) - 1;
    let mut selected = vec![0i32; data.n_targ_haps];
    let step_i = step as i32;
    d[a.len()] = step_i + 1; // set sentinel
    for i in 0..a.len() {
        if (a[i] as usize) < data.n_targ_haps {
            let best_i = lf_best_fwd(fpd, step_i, m_start, m_incl_end, i, a, d, i_to_prev_i, i_to_next_i, data);
            if best_i >= 0 {
                selected[a[i] as usize] = a[best_i as usize];
            } else {
                let mut u = i;
                let mut v = i + 1;
                let mut u_next_match_start = d[u];
                let mut v_next_match_start = d[v];
                while v - u < data.n_candidates
                    && (u_next_match_start <= step_i || v_next_match_start <= step_i)
                {
                    if v_next_match_start <= u_next_match_start {
                        v += 1;
                        v_next_match_start = d[v].max(v_next_match_start);
                    } else {
                        u -= 1;
                        u_next_match_start = d[u].max(u_next_match_start);
                    }
                }
                selected[a[i] as usize] =
                    lf_get_match(fpd, m_start, m_incl_end, i, u, v, a, &mut rand);
            }
        }
    }
    selected
}

#[allow(clippy::too_many_arguments)]
fn lf_best_fwd(
    fpd: &FixedPhaseData,
    step: i32,
    m_start: usize,
    m_incl_end: usize,
    i: usize,
    a: &[i32],
    d: &[i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> i32 {
    let ibs2 = &fpd.stage1_ibs2;
    let mut best_prev_match: i32 = -1;
    let mut best_next_match: i32 = -1;
    let mut prev_match_start: i32 = 0;
    let mut next_match_start: i32 = 0;
    let min_match_start = if i + 1 < a.len() {
        d[i].min(d[i + 1])
    } else {
        d[i]
    };
    let d_max = (min_match_start + data.max_backoff_steps as i32).min(step);
    let mut prev_i = i_to_prev_i[i];
    while prev_i > i32::MIN
        && ibs2.are_ibs2_range(
            (a[i] >> 1) as usize,
            (a[prev_i as usize] >> 1) as usize,
            m_start,
            m_incl_end,
        )
    {
        prev_i = i_to_prev_i[prev_i as usize];
    }
    if prev_i > i32::MIN {
        let mut u = i;
        while (u as i32 - 1) != prev_i && d[u] <= d_max {
            prev_match_start = prev_match_start.max(d[u]);
            u -= 1;
        }
        if (u as i32 - 1) == prev_i && d[u] <= d_max {
            prev_match_start = prev_match_start.max(d[u]);
            best_prev_match = prev_i;
        }
    }
    let mut next_i = i_to_next_i[i];
    while next_i < i32::MAX
        && ibs2.are_ibs2_range(
            (a[i] >> 1) as usize,
            (a[next_i as usize] >> 1) as usize,
            m_start,
            m_incl_end,
        )
    {
        next_i = i_to_next_i[next_i as usize];
    }
    if next_i < i32::MAX {
        let mut v = i;
        while (v as i32 + 1) != next_i && d[v + 1] <= d_max {
            v += 1;
            next_match_start = next_match_start.max(d[v]);
        }
        if (v as i32 + 1) == next_i && d[v + 1] <= d_max {
            v += 1;
            next_match_start = next_match_start.max(d[v]);
            best_next_match = next_i;
        }
    }
    if prev_match_start < next_match_start && best_prev_match != -1 {
        best_prev_match
    } else {
        best_next_match
    }
}

#[allow(clippy::too_many_arguments)]
fn lf_best_bwd(
    fpd: &FixedPhaseData,
    step: i32,
    m_start: usize,
    m_incl_end: usize,
    i: usize,
    a: &[i32],
    d: &[i32],
    i_to_prev_i: &[i32],
    i_to_next_i: &[i32],
    data: &PbwtIbsData,
) -> i32 {
    let ibs2 = &fpd.stage1_ibs2;
    let n_steps_m1 = (fpd.stage1_steps.size() - 1) as i32;
    let mut best_prev_match: i32 = -1;
    let mut best_next_match: i32 = -1;
    let mut prev_match_incl_end: i32 = n_steps_m1;
    let mut next_match_incl_end: i32 = n_steps_m1;
    let max_match_start = if i + 1 < a.len() {
        d[i].max(d[i + 1])
    } else {
        d[i]
    };
    let d_min = (max_match_start - data.max_backoff_steps as i32).max(step);
    let mut prev_i = i_to_prev_i[i];
    while prev_i > i32::MIN
        && ibs2.are_ibs2_range(
            (a[i] >> 1) as usize,
            (a[prev_i as usize] >> 1) as usize,
            m_start,
            m_incl_end,
        )
    {
        prev_i = i_to_prev_i[prev_i as usize];
    }
    if prev_i > i32::MIN {
        let mut u = i;
        while (u as i32 - 1) != prev_i && d[u] >= d_min {
            prev_match_incl_end = prev_match_incl_end.min(d[u]);
            u -= 1;
        }
        if (u as i32 - 1) == prev_i && d[u] >= d_min {
            prev_match_incl_end = prev_match_incl_end.min(d[u]);
            best_prev_match = prev_i;
        }
    }
    let mut next_i = i_to_next_i[i];
    while next_i < i32::MAX
        && ibs2.are_ibs2_range(
            (a[i] >> 1) as usize,
            (a[next_i as usize] >> 1) as usize,
            m_start,
            m_incl_end,
        )
    {
        next_i = i_to_next_i[next_i as usize];
    }
    if next_i < i32::MAX {
        let mut v = i;
        while (v as i32 + 1) != next_i && d[v + 1] >= d_min {
            v += 1;
            next_match_incl_end = next_match_incl_end.min(d[v]);
        }
        if (v as i32 + 1) == next_i && d[v + 1] >= d_min {
            v += 1;
            next_match_incl_end = next_match_incl_end.min(d[v]);
            best_next_match = next_i;
        }
    }
    if prev_match_incl_end > next_match_incl_end && best_prev_match != -1 {
        best_prev_match
    } else {
        best_next_match
    }
}

#[allow(clippy::too_many_arguments)]
fn lf_get_match(
    fpd: &FixedPhaseData,
    m_start: usize,
    m_incl_end: usize,
    i: usize,
    i_start: usize,
    i_end: usize,
    a: &[i32],
    rand: &mut JavaRandom,
) -> i32 {
    let i_length = i_end - i_start;
    if i_length == 1 {
        return -1;
    }
    let ibs2 = &fpd.stage1_ibs2;
    let mut m: i32 = -1;
    let mut index = i_start + rand.next_int_bound(i_length as i32) as usize;
    let mut j = 0;
    while j < i_length && m == -1 {
        if !ibs2.are_ibs2_range(
            (a[i] >> 1) as usize,
            (a[index] >> 1) as usize,
            m_start,
            m_incl_end,
        ) {
            m = a[index];
        }
        index += 1;
        if index == i_end {
            index = i_start;
        }
        j += 1;
    }
    m
}

fn set_inv(a: &[i32], a_inv: &mut [i32]) {
    for (j, &h) in a.iter().enumerate() {
        a_inv[h as usize] = j as i32;
    }
}

fn set_i_to_prev_next_i(
    fpd: &FixedPhaseData,
    step: usize,
    inv_a: &[i32],
    i_to_prev_i: &mut [i32],
    i_to_next_i: &mut [i32],
) {
    i_to_prev_i.fill(i32::MIN);
    i_to_next_i.fill(i32::MAX);
    let steps = &fpd.stage1_steps;
    let hi_freq_indices = &fpd.stage1_to2;
    let start = if step == 0 {
        0
    } else {
        hi_freq_indices[steps.start(step)]
    };
    let end = if step + 1 < steps.size() {
        hi_freq_indices[steps.start(step + 1)]
    } else {
        fpd.targ.len()
    };
    for m in start..end {
        let n_alleles = fpd.targ_markers[m].n_alleles as usize;
        for al in 0..n_alleles {
            if let Carriers::List(carriers) = &fpd.carriers[m][al] {
                if carriers.len() > 1 {
                    // hap list: both haps of every carrier sample
                    let mut idx: Vec<i32> = Vec::with_capacity(carriers.len() * 2);
                    for &sample in carriers.iter() {
                        let h1 = (sample as usize) << 1;
                        idx.push(inv_a[h1]);
                        idx.push(inv_a[h1 | 1]);
                    }
                    idx.sort_unstable();
                    for k in 1..idx.len() {
                        let i0 = idx[k - 1];
                        let i1 = idx[k];
                        if i0 > i_to_prev_i[i1 as usize] {
                            i_to_prev_i[i1 as usize] = i0;
                        }
                        if i1 < i_to_next_i[i0 as usize] {
                            i_to_next_i[i0 as usize] = i1;
                        }
                    }
                }
            }
        }
    }
}

/// Port of `phase.LowFreqPhaseIbs`.
pub struct LowFreqPhaseIbs {
    pub fwd: LowFreqPbwtPhaseIbs,
    pub bwd: LowFreqPbwtPhaseIbs,
}

impl LowFreqPhaseIbs {
    pub fn new(pd: &PhaseData, par: &Par, coded_steps: &CodedSteps) -> LowFreqPhaseIbs {
        LowFreqPhaseIbs {
            fwd: LowFreqPbwtPhaseIbs::new(pd, par, coded_steps, false),
            bwd: LowFreqPbwtPhaseIbs::new(pd, par, coded_steps, true),
        }
    }
}
