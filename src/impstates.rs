//! Port of `imp.ImpStates` and `beagleutil.CompHapSegment`.

use crate::impdata::ImpData;
use crate::impibs::ImpIbs;
use crate::javautil::{IntIntMap, JavaOrd, JavaPriorityQueue, JavaRandom};

const NIL: i32 = -103;

/// Port of `beagleutil.CompHapSegment` (ordered by `lastIbsStep` only; ties
/// are resolved by the priority-queue's heap layout, replicated in
/// `JavaPriorityQueue`).
#[derive(Clone)]
struct CompHapSegment {
    hap: i32,
    last_ibs_step: i32,
    comp_hap_index: usize,
}

impl JavaOrd for CompHapSegment {
    #[inline]
    fn compare_to(&self, other: &Self) -> std::cmp::Ordering {
        self.last_ibs_step.cmp(&other.last_ibs_step)
    }
}

/// Port of `imp.ImpStates`: builds the composite reference haplotypes for
/// one target haplotype and copies per-cluster state data.
pub struct ImpStates {
    n_clusters: usize,
    max_states: usize,
    hap_to_last_ibs_step: IntIntMap,
    q: JavaPriorityQueue<CompHapSegment>,
    comp_hap_hap: Vec<Vec<i32>>,
    comp_hap_end: Vec<Vec<i32>>,
    comp_hap_to_list_index: Vec<usize>,
    comp_hap_to_hap: Vec<i32>,
    comp_hap_to_end: Vec<i32>,
}

impl ImpStates {
    pub fn new(imp_data: &ImpData) -> ImpStates {
        let max_states = imp_data.imp_states;
        ImpStates {
            n_clusters: imp_data.n_clusters,
            max_states,
            hap_to_last_ibs_step: IntIntMap::with_capacity(max_states),
            q: JavaPriorityQueue::new(),
            comp_hap_hap: vec![Vec::new(); max_states],
            comp_hap_end: vec![Vec::new(); max_states],
            comp_hap_to_list_index: vec![0; max_states],
            comp_hap_to_hap: vec![0; max_states],
            comp_hap_to_end: vec![0; max_states],
        }
    }

    /// `ImpStates.ibsStates`: fills `hap_indices[m*max_states + j]` and
    /// `al_match[m*max_states + j]`; returns the number of states.
    pub fn ibs_states(
        &mut self,
        imp_data: &ImpData,
        ibs_haps: &ImpIbs,
        targ_hap: usize,
        hap_indices: &mut [i32],
        al_match: &mut [bool],
    ) -> usize {
        self.initialize_fields();
        let n_steps = ibs_haps.coded_steps.n_steps();
        for j in 0..n_steps {
            let ibs = ibs_haps.ibs_haps(targ_hap, j);
            for k in 0..ibs.len() {
                let hap = ibs[k];
                self.update_fields(ibs_haps, hap, j as i32);
            }
        }
        if self.q.is_empty() {
            self.fill_q_with_random_haps(imp_data, targ_hap);
        }
        self.copy_data(imp_data, targ_hap, hap_indices, al_match)
    }

    fn initialize_fields(&mut self) {
        self.hap_to_last_ibs_step.clear();
        for j in 0..self.q.len() {
            self.comp_hap_hap[j].clear();
            self.comp_hap_end[j].clear();
        }
        self.q.clear();
    }

    fn update_fields(&mut self, ibs_haps: &ImpIbs, hap: i32, step: i32) {
        if self.hap_to_last_ibs_step.get(hap, NIL) == NIL {
            // hap not currently in q
            self.update_head_of_q();
            if self.q.len() == self.max_states {
                let head = self.q.poll().unwrap();
                let start_marker = ibs_haps
                    .coded_steps
                    .step_start((((head.last_ibs_step + step) as u32) >> 1) as usize)
                    as i32;
                self.hap_to_last_ibs_step.remove(head.hap);
                self.comp_hap_hap[head.comp_hap_index].push(hap); // hap of new segment
                self.comp_hap_end[head.comp_hap_index].push(start_marker); // end of previous segment
                let new_head = CompHapSegment {
                    hap,
                    last_ibs_step: step,
                    comp_hap_index: head.comp_hap_index,
                };
                self.q.offer(new_head);
            } else {
                let comp_hap_index = self.q.len();
                self.comp_hap_hap[comp_hap_index].push(hap); // hap of new segment
                self.q.offer(CompHapSegment {
                    hap,
                    last_ibs_step: step,
                    comp_hap_index,
                });
            }
        }
        self.hap_to_last_ibs_step.put(hap, step);
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

    fn copy_data(
        &mut self,
        imp_data: &ImpData,
        targ_hap: usize,
        hap_indices: &mut [i32],
        al_match: &mut [bool],
    ) -> usize {
        let n_comp_haps = self.q.len();
        let shifted_targ_hap = imp_data.n_ref_haps + targ_hap;
        let max_states = self.max_states;
        // initializeCopy
        for j in 0..n_comp_haps {
            self.comp_hap_end[j].push(self.n_clusters as i32); // add missing end of last segment
            self.comp_hap_to_list_index[j] = 0;
            self.comp_hap_to_hap[j] = self.comp_hap_hap[j][0];
            self.comp_hap_to_end[j] = self.comp_hap_end[j][0];
        }
        let n_ref = imp_data.n_ref_haps;
        for m in 0..self.n_clusters {
            let targ_allele = imp_data.allele(m, shifted_targ_hap);
            let row = m * max_states;
            let m_i32 = m as i32;
            for j in 0..n_comp_haps {
                if m_i32 == self.comp_hap_to_end[j] {
                    self.comp_hap_to_list_index[j] += 1;
                    self.comp_hap_to_hap[j] =
                        self.comp_hap_hap[j][self.comp_hap_to_list_index[j]];
                    self.comp_hap_to_end[j] =
                        self.comp_hap_end[j][self.comp_hap_to_list_index[j]];
                }
            }
            let haps = &self.comp_hap_to_hap[..n_comp_haps];
            let hap_out = &mut hap_indices[row..row + n_comp_haps];
            let match_out = &mut al_match[row..row + n_comp_haps];
            match &imp_data.hap_to_seq[m] {
                crate::impdata::ClusterCoding::Composed {
                    hap2seq,
                    seq1_to_seq2,
                    targ,
                    ..
                } => {
                    for j in 0..n_comp_haps {
                        let hap = haps[j];
                        let h = hap as usize;
                        let v = if h < n_ref {
                            seq1_to_seq2[hap2seq[h] as usize]
                        } else {
                            targ[h - n_ref]
                        };
                        hap_out[j] = hap;
                        match_out[j] = v == targ_allele;
                    }
                }
                crate::impdata::ClusterCoding::Direct {
                    coded_ref, targ, ..
                } => {
                    for j in 0..n_comp_haps {
                        let hap = haps[j];
                        let h = hap as usize;
                        let v = if h < n_ref {
                            coded_ref[h]
                        } else {
                            targ[h - n_ref]
                        };
                        hap_out[j] = hap;
                        match_out[j] = v == targ_allele;
                    }
                }
            }
        }
        n_comp_haps
    }

    fn fill_q_with_random_haps(&mut self, imp_data: &ImpData, targ_hap: usize) {
        debug_assert!(self.q.is_empty());
        let n_ref_haps = imp_data.n_ref_haps;
        let n_states = n_ref_haps.min(self.max_states);
        let mut rand = JavaRandom::new(targ_hap as i64);
        let ibs_step = 0;
        for j in 0..n_states {
            let mut h = rand.next_int_bound(n_ref_haps as i32);
            while h == targ_hap as i32 {
                h = rand.next_int_bound(n_ref_haps as i32);
            }
            self.comp_hap_hap[j].push(h); // hap of new segment
            self.q.offer(CompHapSegment {
                hap: h,
                last_ibs_step: ibs_step,
                comp_hap_index: j,
            });
        }
    }
}
