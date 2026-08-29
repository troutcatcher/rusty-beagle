//! Port of `imp.ImpLSBaum`, `imp.StateProbsFactory`, `imp.StateProbs`,
//! and the `imp.ImpLS` driver.

use crate::impdata::ImpData;
use crate::impibs::ImpIbs;
use crate::impstates::ImpStates;
use rayon::prelude::*;

/// One retained HMM state: its reference haplotype and its probability at
/// the cluster and at the next cluster (the pair interpolation needs).
#[derive(Clone, Copy)]
pub struct StateProb {
    pub hap: i32,
    pub prob: f32,
    pub prob_p1: f32,
}

/// Sparse per-cluster HMM state probabilities for one target haplotype
/// (CSR layout; port of `imp.StateProbs`).
///
/// The three values are interleaved rather than held in parallel arrays:
/// they are written together when sparsifying and read together when
/// building output records, so one stream costs a third of the cache misses
/// and a third of the allocations -- and these rows are retained for every
/// target haplotype at once.
pub struct StateProbs {
    offsets: Vec<u32>, // n_clusters + 1
    data: Vec<StateProb>,
}

impl StateProbs {
    /// The retained states of `marker`.
    #[inline]
    pub fn states(&self, marker: usize) -> &[StateProb] {
        let lo = self.offsets[marker] as usize;
        let hi = self.offsets[marker + 1] as usize;
        &self.data[lo..hi]
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
    /// per-cluster bitset of states clearing the sparsification threshold,
    /// filled by the backward pass (which already has the values in hand)
    keep_mask: Vec<u64>,
    max_states: usize,
    /// one cluster's worth of retained states, so the filter writes into a
    /// small L1-resident buffer and each cluster appends in one memcpy
    scratch: Vec<StateProb>,
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
            keep_mask: vec![0; n_markers * words_per_row],
            max_states,
            scratch: vec![
                StateProb {
                    hap: 0,
                    prob: 0.0,
                    prob_p1: 0.0
                };
                max_states
            ],
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
        let threshold = (0.005f32).min(0.9999f32 / n_states as f32);
        for m in (0..=last_marker).rev() {
            last_sum = self.set_bwd_value(imp_data, m, n_states, last_sum, threshold);
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
        threshold: f32,
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
        // Normalize the state probabilities, and while each one is in hand
        // record whether it clears the sparsification threshold here or at
        // the next cluster. The backward pass runs from the last cluster
        // down, so cluster m+1 is already final; recording the mask here
        // saves sparsification a second pass over every state.
        let max_states = self.max_states;
        let words_per_row = self.words_per_row;
        let (head, tail) = self.fwd_val.split_at_mut(row + max_states);
        let fwd = &mut head[row..row + n_states];
        let next_row: &[f32] = if m_p1 < self.n_markers {
            &tail[..n_states]
        } else {
            &[]
        };
        let mask_row = &mut self.keep_mask[m * words_per_row..(m + 1) * words_per_row];
        for (wi, word) in mask_row.iter_mut().enumerate() {
            let base = wi * 64;
            let hi = (base + 64).min(n_states);
            let mut bits = 0u64;
            for j in base..hi {
                let v = fwd[j] / state_sum; // normalize state probabilities
                fwd[j] = v;
                // the last cluster compares against itself, as before
                let nv = if next_row.is_empty() { v } else { next_row[j] };
                bits |= (((v > threshold) | (nv > threshold)) as u64) << (j - base);
            }
            *word = bits;
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
        let mut offsets: Vec<u32> = Vec::with_capacity(n_markers + 1);
        // The backward pass already recorded which states are retained, so
        // the exact final size is a popcount away. Sizing the buffer up front
        // avoids both the repeated growth and the shrink-back that a guessed
        // capacity needs -- and since every target haplotype's rows stay live,
        // each byte over-allocated is a page this run has to fault in.
        let total: usize = self.keep_mask[..n_markers * self.words_per_row]
            .iter()
            .map(|w| w.count_ones() as usize)
            .sum();
        let mut data: Vec<StateProb> = Vec::with_capacity(total);
        offsets.push(0);
        let fwd_val = &self.fwd_val;
        let keep_mask = &self.keep_mask;
        let words_per_row = self.words_per_row;
        let max_states = self.max_states;
        let mut scratch = std::mem::take(&mut self.scratch);
        self.states.replay(n_states, |m, state_haps| {
            let row = m * max_states;
            let row_p1 = if m < n_markers_m1 {
                row + max_states
            } else {
                row
            };
            let cur = &fwd_val[row..row + n_states];
            let nxt = &fwd_val[row_p1..row_p1 + n_states];
            // the backward pass already marked which states are retained, so
            // this walks only the set bits
            let mut k = 0;
            for (wi, &word) in keep_mask[m * words_per_row..(m + 1) * words_per_row]
                .iter()
                .enumerate()
            {
                let mut bits = word;
                while bits != 0 {
                    let j = wi * 64 + bits.trailing_zeros() as usize;
                    scratch[k] = StateProb {
                        hap: state_haps[j],
                        prob: cur[j],
                        prob_p1: nxt[j],
                    };
                    k += 1;
                    bits &= bits - 1;
                }
            }
            data.extend_from_slice(&scratch[..k]);
            offsets.push(data.len() as u32);
        });
        self.scratch = scratch;
        debug_assert_eq!(data.len(), total);
        StateProbs { offsets, data }
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
