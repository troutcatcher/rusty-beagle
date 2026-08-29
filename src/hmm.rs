//! Port of `imp.ImpLSBaum`, `imp.StateProbsFactory`, `imp.StateProbs`,
//! and the `imp.ImpLS` driver.

use crate::impdata::ImpData;
use crate::impibs::ImpIbs;
use crate::impstates::ImpStates;
use rayon::prelude::*;

/// Sparse per-cluster HMM state probabilities for one target haplotype
/// (CSR layout; port of `imp.StateProbs`).
pub struct StateProbs {
    offsets: Vec<u32>, // n_clusters + 1
    haps: Vec<i32>,
    probs: Vec<f32>,
    probs_p1: Vec<f32>,
}

impl StateProbs {
    #[inline]
    pub fn n_states(&self, marker: usize) -> usize {
        (self.offsets[marker + 1] - self.offsets[marker]) as usize
    }

    #[inline]
    pub fn ref_hap(&self, marker: usize, index: usize) -> i32 {
        self.haps[self.offsets[marker] as usize + index]
    }

    #[inline]
    pub fn probs(&self, marker: usize, index: usize) -> f32 {
        self.probs[self.offsets[marker] as usize + index]
    }

    #[inline]
    pub fn probs_p1(&self, marker: usize, index: usize) -> f32 {
        self.probs_p1[self.offsets[marker] as usize + index]
    }
}

/// Per-thread reusable buffers (port of `imp.ImpLSBaum`).
struct Baum {
    states: ImpStates,
    n_markers: usize,
    words_per_row: usize,
    alleles_match: Vec<u64>, // n_markers x words_per_row bitset
    fwd_val: Vec<f32>,       // n_markers x max_states
    bwd_val: Vec<f32>,       // max_states
    max_states: usize,
}

impl Baum {
    fn new(imp_data: &ImpData) -> Baum {
        let n_markers = imp_data.n_clusters;
        let max_states = imp_data.imp_states;
        let words_per_row = (max_states + 63) >> 6;
        Baum {
            states: ImpStates::new(imp_data),
            n_markers,
            words_per_row,
            alleles_match: vec![0; n_markers * words_per_row],
            fwd_val: vec![0.0; n_markers * max_states],
            bwd_val: vec![0.0; max_states],
            max_states,
        }
    }

    /// `ImpLSBaum.impute(targHap)`
    fn impute(&mut self, imp_data: &ImpData, ibs_haps: &ImpIbs, targ_hap: usize) -> StateProbs {
        let last_marker = imp_data.n_clusters - 1;
        let t0 = std::time::Instant::now();
        let n_states = self.states.ibs_states(
            imp_data,
            ibs_haps,
            targ_hap,
            &mut self.alleles_match,
            self.words_per_row,
        );
        let t1 = std::time::Instant::now();
        self.set_fwd_values(imp_data, n_states);
        let t2 = std::time::Instant::now();
        self.bwd_val[..n_states].fill(1.0 / n_states as f32);
        let mut last_sum = 1.0f32;
        for m in (0..=last_marker).rev() {
            last_sum = self.set_bwd_value(imp_data, m, n_states, last_sum);
        }
        let t3 = std::time::Instant::now();
        let r = self.state_probs(n_states);
        phase_add(0, (t1 - t0).as_nanos() as u64);
        phase_add(1, (t2 - t1).as_nanos() as u64);
        phase_add(2, (t3 - t2).as_nanos() as u64);
        phase_add(3, t3.elapsed().as_nanos() as u64);
        r
    }

    fn set_fwd_values(&mut self, imp_data: &ImpData, n_states: usize) {
        let mut last_sum = 1.0f32;
        for m in 0..self.n_markers {
            let p_recomb = imp_data.p_recomb[m];
            let p_err = imp_data.err_prob[m];
            let p_no_err = 1.0f32 - p_err;
            let shift = p_recomb / n_states as f32;
            let scale = (1.0f32 - p_recomb) / last_sum;
            let mut sum = 0.0f32;
            let row = m * self.max_states;
            let bits = &self.alleles_match[m * self.words_per_row..];
            if m == 0 {
                for j in 0..n_states {
                    let is_match = (bits[j >> 6] >> (j & 63)) & 1 != 0;
                    let em = if is_match { p_no_err } else { p_err };
                    self.fwd_val[row + j] = em;
                    sum += em;
                }
            } else {
                let (prev_rows, cur_rows) = self.fwd_val.split_at_mut(row);
                let prev = &prev_rows[row - self.max_states..];
                let cur = &mut cur_rows[..n_states];
                for j in 0..n_states {
                    let is_match = (bits[j >> 6] >> (j & 63)) & 1 != 0;
                    let em = if is_match { p_no_err } else { p_err };
                    let v = em * (scale * prev[j] + shift);
                    cur[j] = v;
                    sum += v;
                }
            }
            last_sum = sum;
        }
    }

    fn set_bwd_value(
        &mut self,
        imp_data: &ImpData,
        m: usize,
        n_states: usize,
        last_sum: f32,
    ) -> f32 {
        let m_p1 = m + 1;
        let p_recomb = if m_p1 < self.n_markers {
            imp_data.p_recomb[m_p1]
        } else {
            0.0f32
        };
        let p_err = imp_data.err_prob[m];
        let p_no_err = 1.0f32 - p_err;
        let scale = (1.0f32 - p_recomb) / last_sum;
        let shift = p_recomb / n_states as f32;
        let mut bwd_val_sum = 0.0f32;
        let mut state_sum = 0.0f32;
        let row = m * self.max_states;
        let bits = &self.alleles_match[m * self.words_per_row..];
        let fwd = &mut self.fwd_val[row..row + n_states];
        let bwd = &mut self.bwd_val[..n_states];
        for j in 0..n_states {
            let b = scale * bwd[j] + shift; // finish calculating bwd value
            let f = fwd[j] * b; // store state probabilities in fwd_val[m]
            fwd[j] = f;
            state_sum += f;

            let is_match = (bits[j >> 6] >> (j & 63)) & 1 != 0;
            let em = if is_match { p_no_err } else { p_err };
            let b2 = b * em;
            bwd[j] = b2;
            bwd_val_sum += b2;
        }
        for f in fwd.iter_mut() {
            *f /= state_sum; // normalize state probabilities
        }
        bwd_val_sum
    }

    /// `StateProbsFactory.stateProbs`: sparsify, keeping states whose
    /// probability at marker m or m+1 exceeds `min(0.005, 0.9999/nStates)`.
    /// State haplotype indices are re-derived by replaying the composite
    /// haplotype segments.
    fn state_probs(&mut self, n_states: usize) -> StateProbs {
        let n_markers = self.n_markers;
        let n_markers_m1 = n_markers - 1;
        let threshold = (0.005f32).min(0.9999f32 / n_states as f32);
        let mut offsets: Vec<u32> = Vec::with_capacity(n_markers + 1);
        let mut haps: Vec<i32> = Vec::new();
        let mut probs: Vec<f32> = Vec::new();
        let mut probs_p1: Vec<f32> = Vec::new();
        offsets.push(0);
        let fwd_val = &self.fwd_val;
        let max_states = self.max_states;
        self.states.replay(n_states, |m, state_haps| {
            let row = m * max_states;
            let row_p1 = if m < n_markers_m1 {
                row + max_states
            } else {
                row
            };
            for j in 0..n_states {
                if fwd_val[row + j] > threshold || fwd_val[row_p1 + j] > threshold {
                    haps.push(state_haps[j]);
                    probs.push(fwd_val[row + j]);
                    probs_p1.push(fwd_val[row_p1 + j]);
                }
            }
            offsets.push(haps.len() as u32);
        });
        StateProbs {
            offsets,
            haps,
            probs,
            probs_p1,
        }
    }
}

/// Port of `imp.ImpLS.stateProbs`: HMM state probabilities for every target
/// haplotype, computed in parallel.
pub fn state_probs(imp_data: &ImpData, ibs_haps: &ImpIbs) -> Vec<StateProbs> {
    let n_targ_haps = imp_data.n_targ_haps;
    let timing = std::env::var("RUSTY_BEAGLE_TIMING2").is_ok();
    let result: Vec<StateProbs> = (0..n_targ_haps)
        .into_par_iter()
        .map_init(
            || Baum::new(imp_data),
            |baum, hap| baum.impute(imp_data, ibs_haps, hap),
        )
        .collect();
    if timing {
        let t: Vec<u64> = PHASE_NANOS
            .iter()
            .map(|a| a.load(std::sync::atomic::Ordering::Relaxed))
            .collect();
        eprintln!(
            "[timing2/cpu-total] states: {:.3}s (qbuild: {:.3}s) fwd: {:.3}s bwd: {:.3}s sparsify: {:.3}s",
            t[0] as f64 / 1e9,
            t[4] as f64 / 1e9,
            t[1] as f64 / 1e9,
            t[2] as f64 / 1e9,
            t[3] as f64 / 1e9
        );
    }
    result
}

static PHASE_NANOS: [std::sync::atomic::AtomicU64; 5] = [
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
    std::sync::atomic::AtomicU64::new(0),
];

#[inline]
pub fn phase_add(idx: usize, nanos: u64) {
    PHASE_NANOS[idx].fetch_add(nanos, std::sync::atomic::Ordering::Relaxed);
}
