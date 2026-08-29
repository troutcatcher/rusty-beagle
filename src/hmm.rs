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
    /// number of retained states in the previous haplotype's `StateProbs`;
    /// sparsity is similar between haplotypes, so this preallocates the CSR
    /// arrays instead of growing them by repeated doubling
    size_hint: usize,
    /// one cluster's worth of retained states, so the filter writes into a
    /// small L1-resident buffer and each cluster appends in one memcpy
    scr_haps: Vec<i32>,
    scr_probs: Vec<f32>,
    scr_probs_p1: Vec<f32>,
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
            size_hint: 0,
            scr_haps: vec![0; max_states],
            scr_probs: vec![0.0; max_states],
            scr_probs_p1: vec![0.0; max_states],
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
                let cur = &mut self.fwd_val[row..row + n_states];
                let mut word = 0u64;
                for (j, c) in cur.iter_mut().enumerate() {
                    if j & 63 == 0 {
                        word = bits[j >> 6];
                    }
                    let em = if word & 1 != 0 { p_no_err } else { p_err };
                    word >>= 1;
                    *c = em;
                    sum += em;
                }
            } else {
                let (prev_rows, cur_rows) = self.fwd_val.split_at_mut(row);
                let prev = &prev_rows[row - self.max_states..];
                let cur = &mut cur_rows[..n_states];
                let mut word = 0u64;
                for (j, (c, &p)) in cur.iter_mut().zip(prev.iter()).enumerate() {
                    if j & 63 == 0 {
                        word = bits[j >> 6];
                    }
                    let em = if word & 1 != 0 { p_no_err } else { p_err };
                    word >>= 1;
                    let v = em * (scale * p + shift);
                    *c = v;
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
        let mut word = 0u64;
        for (j, (f_slot, b_slot)) in fwd.iter_mut().zip(bwd.iter_mut()).enumerate() {
            let b = scale * *b_slot + shift; // finish calculating bwd value
            let f = *f_slot * b; // store state probabilities in fwd_val[m]
            *f_slot = f;
            state_sum += f;

            if j & 63 == 0 {
                word = bits[j >> 6];
            }
            let em = if word & 1 != 0 { p_no_err } else { p_err };
            word >>= 1;
            let b2 = b * em;
            *b_slot = b2;
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
        let cap = self.size_hint;
        let mut haps: Vec<i32> = Vec::with_capacity(cap);
        let mut probs: Vec<f32> = Vec::with_capacity(cap);
        let mut probs_p1: Vec<f32> = Vec::with_capacity(cap);
        offsets.push(0);
        let fwd_val = &self.fwd_val;
        let max_states = self.max_states;
        let mut scr_haps = std::mem::take(&mut self.scr_haps);
        let mut scr_probs = std::mem::take(&mut self.scr_probs);
        let mut scr_probs_p1 = std::mem::take(&mut self.scr_probs_p1);
        self.states.replay(n_states, |m, state_haps| {
            let row = m * max_states;
            let row_p1 = if m < n_markers_m1 {
                row + max_states
            } else {
                row
            };
            let cur = &fwd_val[row..row + n_states];
            let nxt = &fwd_val[row_p1..row_p1 + n_states];
            let mut k = 0;
            for j in 0..n_states {
                let a = cur[j];
                let b = nxt[j];
                if a > threshold || b > threshold {
                    scr_haps[k] = state_haps[j];
                    scr_probs[k] = a;
                    scr_probs_p1[k] = b;
                    k += 1;
                }
            }
            haps.extend_from_slice(&scr_haps[..k]);
            probs.extend_from_slice(&scr_probs[..k]);
            probs_p1.extend_from_slice(&scr_probs_p1[..k]);
            offsets.push(haps.len() as u32);
        });
        self.scr_haps = scr_haps;
        self.scr_probs = scr_probs;
        self.scr_probs_p1 = scr_probs_p1;
        self.size_hint = haps.len();
        // `size_hint` only approximates this haplotype's sparsity, so the rows
        // can retain spare capacity. Every target haplotype's rows are held at
        // once, so on large cohorts that waste dominates peak memory; reclaim
        // it when it is worth a copy.
        if haps.capacity() > haps.len() + (haps.len() >> 3) {
            haps.shrink_to_fit();
            probs.shrink_to_fit();
            probs_p1.shrink_to_fit();
        }
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
