//! Port of `imp.CodedSteps` and `imp.ImpIbs`.

use crate::impdata::ImpData;
use crate::javautil::{java_shuffle, JavaRandom};
use rayon::prelude::*;
use std::sync::Arc;

/// Compact per-step hap→sequence map (Java packs bits; we pick the smallest
/// unsigned width that fits `value_size`).
pub enum PackedSeq {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
}

impl PackedSeq {
    fn new(vals: &[u32], value_size: u32) -> PackedSeq {
        if value_size <= u8::MAX as u32 + 1 {
            PackedSeq::U8(vals.iter().map(|&v| v as u8).collect())
        } else if value_size <= u16::MAX as u32 + 1 {
            PackedSeq::U16(vals.iter().map(|&v| v as u16).collect())
        } else {
            PackedSeq::U32(vals.to_vec())
        }
    }

    #[inline]
    pub fn get(&self, i: usize) -> u32 {
        match self {
            PackedSeq::U8(v) => v[i] as u32,
            PackedSeq::U16(v) => v[i] as u32,
            PackedSeq::U32(v) => v[i],
        }
    }
}

/// Port of `imp.CodedSteps`.
pub struct CodedSteps {
    step_starts: Vec<usize>,
    coded_steps: Vec<(PackedSeq, u32)>, // (hap→seq, valueSize)
}

impl CodedSteps {
    pub fn new(imp_data: &ImpData) -> CodedSteps {
        let step_starts = step_starts(imp_data);
        let coded_steps: Vec<(PackedSeq, u32)> = (0..step_starts.len())
            .into_par_iter()
            .map(|j| code_step(imp_data, &step_starts, j))
            .collect();
        CodedSteps {
            step_starts,
            coded_steps,
        }
    }

    pub fn n_steps(&self) -> usize {
        self.step_starts.len()
    }

    pub fn step_start(&self, step: usize) -> usize {
        self.step_starts[step]
    }

    pub fn get(&self, step: usize) -> &(PackedSeq, u32) {
        &self.coded_steps[step]
    }
}

fn step_starts(imp_data: &ImpData) -> Vec<usize> {
    let pos = &imp_data.pos;
    let step = imp_data.imp_step as f64;
    let mut indices = Vec::with_capacity(pos.len() / 10 + 1);
    indices.push(0usize);
    let mut next_pos = pos[0] + step / 2.0; // first step is half-length
    let mut index = next_index(pos, 0, next_pos);
    while index < pos.len() {
        indices.push(index);
        next_pos = pos[index] + step;
        index = next_index(pos, index, next_pos);
    }
    indices
}

/// `Arrays.binarySearch(pos, start, pos.length, targetPos)` insertion point.
fn next_index(pos: &[f64], start: usize, target_pos: f64) -> usize {
    start + pos[start..].partition_point(|&p| p < target_pos)
}

fn code_step(imp_data: &ImpData, starts: &[usize], start_index: usize) -> (PackedSeq, u32) {
    let n_ref_haps = imp_data.n_ref_haps;
    let n_haps = imp_data.n_haps;
    let mut hap_to_seq = vec![1u32; n_haps];
    let start = starts[start_index];
    let end = if start_index + 1 < starts.len() {
        starts[start_index + 1]
    } else {
        imp_data.n_clusters
    };

    let mut seq_cnt: u32 = 2; // seq 0 is reserved for sequences not found in target
    for m in start..end {
        let h2s = &imp_data.hap_to_seq[m];
        let n_alleles = h2s.value_size();
        let mut seq_map = vec![0u32; (seq_cnt * n_alleles) as usize];
        seq_cnt = 1;
        for h in n_ref_haps..n_haps {
            let index = (n_alleles * hap_to_seq[h] + h2s.get(h, n_ref_haps)) as usize;
            if seq_map[index] == 0 {
                seq_map[index] = seq_cnt;
                seq_cnt += 1;
            }
            hap_to_seq[h] = seq_map[index];
        }
        for h in 0..n_ref_haps {
            if hap_to_seq[h] != 0 {
                let index = (hap_to_seq[h] * n_alleles + h2s.get(h, n_ref_haps)) as usize;
                hap_to_seq[h] = seq_map[index];
            }
        }
    }
    (PackedSeq::new(&hap_to_seq, seq_cnt), seq_cnt)
}

/// Port of `imp.ImpIbs`.
#[allow(dead_code)]
pub struct ImpIbs {
    pub coded_steps: CodedSteps,
    n_ref_haps: usize,
    seed: i64,
    n_steps: usize,       // imp_nsteps
    n_haps_per_step: usize,
    /// `ibs_haps[step][targ_hap]` = shared list of IBS reference haplotypes
    ibs_haps: Vec<Vec<Arc<Vec<i32>>>>,
}

impl ImpIbs {
    pub fn new(imp_data: &ImpData) -> ImpIbs {
        let coded_steps = CodedSteps::new(imp_data);
        let n_ref_haps = imp_data.n_ref_haps;
        let seed = imp_data.seed;
        let n_steps = imp_data.imp_nsteps;
        let n_steps_per_segment =
            (imp_data.imp_segment / imp_data.imp_step).round() as usize;
        let n_haps_per_step = imp_data.imp_states / n_steps_per_segment;

        let helper = IbsHelper {
            imp_data,
            coded_steps: &coded_steps,
            n_ref_haps,
            seed,
            n_steps,
            n_haps_per_step,
        };
        let ibs_haps: Vec<Vec<Arc<Vec<i32>>>> = (0..coded_steps.n_steps())
            .into_par_iter()
            .map(|j| helper.get_ibs_haps(j))
            .collect();
        ImpIbs {
            coded_steps,
            n_ref_haps,
            seed,
            n_steps,
            n_haps_per_step,
            ibs_haps,
        }
    }

    /// `ImpIbs.ibsHaps(hap, step)`
    #[inline]
    pub fn ibs_haps(&self, hap: usize, step: usize) -> &[i32] {
        &self.ibs_haps[step][hap]
    }

    #[allow(dead_code)]
    pub fn n_haps_per_step(&self) -> usize {
        self.n_haps_per_step
    }
}

struct IbsHelper<'a> {
    imp_data: &'a ImpData,
    coded_steps: &'a CodedSteps,
    n_ref_haps: usize,
    seed: i64,
    n_steps: usize,
    n_haps_per_step: usize,
}

impl<'a> IbsHelper<'a> {
    fn get_ibs_haps(&self, index: usize) -> Vec<Arc<Vec<i32>>> {
        let n_targ_haps = self.imp_data.n_targ_haps;
        let empty: Arc<Vec<i32>> = Arc::new(Vec::new());
        let mut results: Vec<Arc<Vec<i32>>> = vec![empty; n_targ_haps];
        let n_steps_to_merge = self.n_steps.min(self.coded_steps.n_steps() - index);
        let mut children = self.init_partition(self.coded_steps.get(index));
        let mut next_parents: Vec<Vec<i32>> = Vec::with_capacity(children.len());
        self.init_update_results(&mut children, &mut next_parents, &mut results);
        for i in 1..n_steps_to_merge {
            let parents = std::mem::take(&mut next_parents);
            let coded_step = self.coded_steps.get(index + i);
            for parent in parents {
                let mut children = self.partition(&parent, coded_step);
                self.update_results(&parent, &mut children, &mut next_parents, &mut results);
            }
        }
        self.final_update_results(next_parents, &mut results);
        results
    }

    fn init_partition(&self, coded_step: &(PackedSeq, u32)) -> Vec<Vec<i32>> {
        let (hap2seq, value_size) = coded_step;
        let n_haps = self.imp_data.n_haps;
        let mut seq_to_child: Vec<i32> = vec![-1; *value_size as usize];
        let mut children: Vec<Vec<i32>> = Vec::new();
        for h in self.n_ref_haps..n_haps {
            let seq = hap2seq.get(h) as usize;
            if seq_to_child[seq] < 0 {
                seq_to_child[seq] = children.len() as i32;
                children.push(Vec::new());
            }
        }
        for h in 0..n_haps {
            let seq = hap2seq.get(h) as usize;
            let c = seq_to_child[seq];
            if c >= 0 {
                children[c as usize].push(h as i32);
            }
        }
        children
    }

    fn partition(&self, parent: &[i32], coded_step: &(PackedSeq, u32)) -> Vec<Vec<i32>> {
        let (hap2seq, value_size) = coded_step;
        let mut seq_to_child: Vec<i32> = vec![-1; *value_size as usize];
        let mut children: Vec<Vec<i32>> = Vec::new();
        let targ_start = ins_pt(parent, self.n_ref_haps as i32);
        for &hap in &parent[targ_start..] {
            let seq = hap2seq.get(hap as usize) as usize;
            if seq_to_child[seq] < 0 {
                seq_to_child[seq] = children.len() as i32;
                children.push(Vec::new());
            }
        }
        for &hap in parent {
            let seq = hap2seq.get(hap as usize) as usize;
            let c = seq_to_child[seq];
            if c >= 0 {
                children[c as usize].push(hap);
            }
        }
        children
    }

    fn init_update_results(
        &self,
        children: &mut Vec<Vec<i32>>,
        next_parents: &mut Vec<Vec<i32>>,
        results: &mut [Arc<Vec<i32>>],
    ) {
        for child in children.drain(..) {
            let n_ref = ins_pt(&child, self.n_ref_haps as i32);
            if n_ref <= self.n_haps_per_step {
                let ibs_list: Vec<i32> = child[..n_ref].to_vec();
                self.set_result(&child, n_ref, Arc::new(ibs_list), results);
            } else {
                next_parents.push(child);
            }
        }
    }

    fn update_results(
        &self,
        parent: &[i32],
        children: &mut Vec<Vec<i32>>,
        next_ibs: &mut Vec<Vec<i32>>,
        results: &mut [Arc<Vec<i32>>],
    ) {
        for child in children.drain(..) {
            let n_child_ref = ins_pt(&child, self.n_ref_haps as i32);
            if n_child_ref <= self.n_haps_per_step {
                let ibs_list = self.ibs_haps_list(parent, &child, n_child_ref);
                self.set_result(&child, n_child_ref, Arc::new(ibs_list), results);
            } else {
                next_ibs.push(child);
            }
        }
    }

    /// `ImpIbs.ibsHaps(parent, child, nChildRef)`
    fn ibs_haps_list(&self, parent: &[i32], child: &[i32], n_child_ref: usize) -> Vec<i32> {
        let mut combined: Vec<i32> = Vec::with_capacity(self.n_haps_per_step);
        combined.extend_from_slice(&child[..n_child_ref]);
        let size = self.n_haps_per_step - n_child_ref;
        let mut rand = JavaRandom::new(self.seed.wrapping_add(parent[0] as i64));
        let uniq_to_parent = self.uniq_to_parent(parent, child, n_child_ref);
        let rand_subset = random_subset(uniq_to_parent, size, &mut rand);
        combined.extend_from_slice(&rand_subset);
        combined.sort_unstable();
        combined
    }

    fn uniq_to_parent(&self, parent: &[i32], child: &[i32], n_child_ref: usize) -> Vec<i32> {
        let n_child_ref_m1 = n_child_ref as isize - 1;
        let n_parent_ref = ins_pt(parent, self.n_ref_haps as i32);
        let mut uniq = Vec::with_capacity(parent.len());
        let mut c: isize = 0;
        // Note: Java reads child.get(0) even when nChildRef==0; the child
        // list always contains at least one (target) haplotype.
        let mut c_val = child[0];
        for &p_val in &parent[..n_parent_ref] {
            while c_val < p_val && c < n_child_ref_m1 {
                c += 1;
                c_val = child[c as usize];
            }
            if p_val != c_val {
                uniq.push(p_val);
            }
        }
        uniq
    }

    fn final_update_results(
        &self,
        children: Vec<Vec<i32>>,
        results: &mut [Arc<Vec<i32>>],
    ) {
        for child in children {
            let n_ref = ins_pt(&child, self.n_ref_haps as i32);
            let mut ibs_list: Vec<i32> = child[..n_ref].to_vec();
            if self.n_haps_per_step < ibs_list.len() {
                let mut rand = JavaRandom::new(self.seed.wrapping_add(child[0] as i64));
                java_shuffle(&mut ibs_list, &mut rand);
                ibs_list.truncate(self.n_haps_per_step);
                ibs_list.sort_unstable();
            }
            self.set_result(&child, n_ref, Arc::new(ibs_list), results);
        }
    }

    fn set_result(
        &self,
        child: &[i32],
        first_targ_index: usize,
        ibs_haps: Arc<Vec<i32>>,
        results: &mut [Arc<Vec<i32>>],
    ) {
        for &h in &child[first_targ_index..] {
            results[h as usize - self.n_ref_haps] = ibs_haps.clone();
        }
    }
}

/// `IntList.binarySearch(nRefHaps)` insertion point on a sorted list.
#[inline]
fn ins_pt(list: &[i32], n_ref_haps: i32) -> usize {
    list.partition_point(|&h| h < n_ref_haps)
}

/// `ImpIbs.randomSubset`
fn random_subset(mut list: Vec<i32>, mut size: usize, rand: &mut JavaRandom) -> Vec<i32> {
    if list.len() < size {
        size = list.len();
    }
    for j in 0..size {
        let x = rand.next_int_bound((list.len() - j) as i32) as usize;
        list.swap(j, j + x);
    }
    list.truncate(size);
    list
}
