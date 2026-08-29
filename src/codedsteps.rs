//! Port of `phase.CodedSteps`: indexes unique allele sequences within each
//! step interval, over the combined (target-first, then reference) haplotypes.

use crate::phasedata::{EstPhase, Steps};
use crate::xref::XRefGT;
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;

pub struct CodedSteps {
    pub all_haps: Arc<XRefGT>,
    pub coded_steps: Vec<(Vec<u32>, u32)>, // (hap->seq, valueSize)
}

#[allow(dead_code)]
impl CodedSteps {
    /// `new CodedSteps(estPhase)`
    pub fn new(est_phase: &EstPhase, n_threads: usize) -> CodedSteps {
        let fpd = &est_phase.fpd;
        let targ_haps = XRefGT::from_est_phase(est_phase);
        let all_haps = match &fpd.stage1_xref {
            Some(ref_haps) => Arc::new(XRefGT::combine(targ_haps, ref_haps)),
            None => Arc::new(targ_haps),
        };
        let steps = &fpd.stage1_steps;
        let coded_steps = coded_steps(&all_haps, steps, n_threads);
        CodedSteps {
            all_haps,
            coded_steps,
        }
    }

    #[inline]
    pub fn get(&self, step: usize) -> &(Vec<u32>, u32) {
        &self.coded_steps[step]
    }

    pub fn n_steps(&self) -> usize {
        self.coded_steps.len()
    }
}

fn coded_steps(gt: &XRefGT, steps: &Steps, n_threads: usize) -> Vec<(Vec<u32>, u32)> {
    let max_steps_per_batch = 512usize;
    let mut n_steps_per_batch = (steps.size() + n_threads - 1) / n_threads;
    while n_steps_per_batch > max_steps_per_batch {
        n_steps_per_batch = (n_steps_per_batch + 1) >> 1;
    }
    let steps_per_batch = n_steps_per_batch;
    let n_batches = (steps.size() + steps_per_batch - 1) / steps_per_batch;
    (0..n_batches)
        .into_par_iter()
        .flat_map_iter(|batch| coded_steps_batch(gt, steps, batch, steps_per_batch))
        .collect()
}

fn coded_steps_batch(
    gt: &XRefGT,
    steps: &Steps,
    batch: usize,
    batch_size: usize,
) -> Vec<(Vec<u32>, u32)> {
    let start_step = batch * batch_size;
    let end_step = (start_step + batch_size).min(steps.size());
    let n_steps = end_step - start_step;
    let n_haps = gt.n_haps();
    let mut hap_to_seq: Vec<Vec<u32>> = vec![vec![0u32; n_haps]; n_steps];
    let mut seq_cnt = vec![0u32; n_steps];
    let mut seq_map: Vec<HashMap<i32, u32>> = (0..n_steps).map(|_| HashMap::new()).collect();
    for h in 0..n_haps {
        let mut m_start = steps.start(start_step);
        for j in 0..n_steps {
            let m_end = steps.end(start_step + j);
            let key = gt.hash(h, m_start, m_end);
            let seq_index = match seq_map[j].get(&key) {
                Some(&v) => v,
                None => {
                    let v = seq_cnt[j];
                    seq_cnt[j] += 1;
                    seq_map[j].insert(key, v);
                    v
                }
            };
            hap_to_seq[j][h] = seq_index;
            m_start = m_end;
        }
    }
    hap_to_seq
        .into_iter()
        .zip(seq_cnt)
        .map(|(v, c)| (v, c))
        .collect()
}
