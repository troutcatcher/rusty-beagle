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
    max_states: usize,
    hap_indices: Vec<i32>,   // n_markers x max_states
    alleles_match: Vec<bool>, // n_markers x max_states
    fwd_val: Vec<f32>,       // n_markers x max_states
    bwd_val: Vec<f32>,       // max_states
}

impl Baum {
    fn new(imp_data: &ImpData) -> Baum {
        let n_markers = imp_data.n_clusters;
        let max_states = imp_data.imp_states;
        Baum {
            states: ImpStates::new(imp_data),
            n_markers,
            max_states,
            hap_indices: vec![0; n_markers * max_states],
            alleles_match: vec![false; n_markers * max_states],
            fwd_val: vec![0.0; n_markers * max_states],
            bwd_val: vec![0.0; max_states],
        }
    }

    /// `ImpLSBaum.impute(targHap)`
    fn impute(&mut self, imp_data: &ImpData, ibs_haps: &ImpIbs, targ_hap: usize) -> StateProbs {
        let last_marker = imp_data.n_clusters - 1;
        let n_states = self.states.ibs_states(
            imp_data,
            ibs_haps,
            targ_hap,
            &mut self.hap_indices,
            &mut self.alleles_match,
        );
        self.set_fwd_values(imp_data, n_states);
        self.bwd_val[..n_states].fill(1.0 / n_states as f32);
        let mut last_sum = 1.0f32;
        for m in (0..=last_marker).rev() {
            last_sum = self.set_bwd_value(imp_data, m, n_states, last_sum);
        }
        self.state_probs(n_states)
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
            let prev_row = row.wrapping_sub(self.max_states);
            for j in 0..n_states {
                let em = if self.alleles_match[row + j] {
                    p_no_err
                } else {
                    p_err
                };
                let v = if m == 0 {
                    em
                } else {
                    em * (scale * self.fwd_val[prev_row + j] + shift)
                };
                self.fwd_val[row + j] = v;
                sum += v;
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
        for j in 0..n_states {
            self.bwd_val[j] = scale * self.bwd_val[j] + shift; // finish calculating bwd value
            self.fwd_val[row + j] *= self.bwd_val[j]; // store state probabilities in fwd_val[m]
            state_sum += self.fwd_val[row + j];

            let em = if self.alleles_match[row + j] {
                p_no_err
            } else {
                p_err
            };
            self.bwd_val[j] *= em;
            bwd_val_sum += self.bwd_val[j];
        }
        for j in 0..n_states {
            self.fwd_val[row + j] /= state_sum; // normalize state probabilities
        }
        bwd_val_sum
    }

    /// `StateProbsFactory.stateProbs`: sparsify, keeping states whose
    /// probability at marker m or m+1 exceeds `min(0.005, 0.9999/nStates)`.
    fn state_probs(&self, n_states: usize) -> StateProbs {
        let n_markers = self.n_markers;
        let n_markers_m1 = n_markers - 1;
        let threshold = (0.005f32).min(0.9999f32 / n_states as f32);
        let mut offsets: Vec<u32> = Vec::with_capacity(n_markers + 1);
        let mut haps: Vec<i32> = Vec::new();
        let mut probs: Vec<f32> = Vec::new();
        let mut probs_p1: Vec<f32> = Vec::new();
        offsets.push(0);
        for m in 0..n_markers {
            let row = m * self.max_states;
            let row_p1 = if m < n_markers_m1 {
                row + self.max_states
            } else {
                row
            };
            for j in 0..n_states {
                if self.fwd_val[row + j] > threshold || self.fwd_val[row_p1 + j] > threshold {
                    haps.push(self.hap_indices[row + j]);
                    probs.push(self.fwd_val[row + j]);
                    probs_p1.push(self.fwd_val[row_p1 + j]);
                }
            }
            offsets.push(haps.len() as u32);
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
    (0..n_targ_haps)
        .into_par_iter()
        .map_init(
            || Baum::new(imp_data),
            |baum, hap| baum.impute(imp_data, ibs_haps, hap),
        )
        .collect()
}
