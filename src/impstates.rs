//! Port of `imp.ImpStates` and `beagleutil.CompHapSegment`.

use crate::impdata::ImpData;
use crate::impibs::ImpIbs;
use crate::javautil::{IntIntMap, JavaOrd, JavaPriorityQueue, JavaRandom};

const NIL: i32 = -103;

/// Packs per-state match booleans into a bitset row (whole-word writes, so
/// no prior clearing is needed; bits at or above `haps.len()` are zero).
#[inline]
fn fill_match_bits<F: Fn(usize) -> bool>(haps: &[i32], out: &mut [u64], is_match: F) {
    let mut w = 0u64;
    let mut word_idx = 0;
    for (j, &hap) in haps.iter().enumerate() {
        if is_match(hap as usize) {
            w |= 1u64 << (j & 63);
        }
        if j & 63 == 63 {
            out[word_idx] = w;
            word_idx += 1;
            w = 0;
        }
    }
    if !haps.is_empty() && haps.len() & 63 != 0 {
        out[word_idx] = w;
    }
}

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
    comp_hap_to_hap: Vec<i32>,
    /// segment transitions as `(cluster, comp_hap_index, new_hap)`, sorted by
    /// cluster; see `build_events`
    events: Vec<(u32, u32, i32)>,
    cache_seq: Vec<u16>,
    dirty: Vec<u32>,
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
            comp_hap_to_hap: vec![0; max_states],
            events: Vec::new(),
            cache_seq: vec![0; max_states],
            dirty: Vec::new(),
        }
    }

    /// `ImpStates.ibsStates`: fills the allele-match bitset
    /// (`al_match[m*words_per_row + (j>>6)]` bit `j&63`); returns the number
    /// of states.  Per-state haplotype indices are re-derived afterwards via
    /// `replay` (this avoids materializing the nClusters x maxStates
    /// haplotype matrix that Java builds).
    pub fn ibs_states(
        &mut self,
        imp_data: &ImpData,
        ibs_haps: &ImpIbs,
        targ_hap: usize,
        al_match: &mut [u64],
        words_per_row: usize,
    ) -> usize {
        let t0 = std::time::Instant::now();
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
        crate::hmm::phase_add(4, t0.elapsed().as_nanos() as u64);
        self.copy_data(imp_data, targ_hap, al_match, words_per_row)
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
        al_match: &mut [u64],
        words_per_row: usize,
    ) -> usize {
        let n_comp_haps = self.q.len();
        let shifted_targ_hap = imp_data.n_ref_haps + targ_hap;
        // initializeCopy
        for j in 0..n_comp_haps {
            self.comp_hap_end[j].push(self.n_clusters as i32); // add missing end of last segment
        }
        self.build_events(n_comp_haps);
        self.reset_comp_haps(n_comp_haps);
        let events = std::mem::take(&mut self.events);
        let mut ev = 0usize;
        let n_ref = imp_data.n_ref_haps;
        // Per-state cache of the seq-coded block sequence index
        // (hap2seq[state hap]): valid while the state's haplotype and the
        // cluster's block are unchanged.  This avoids the random-access
        // hap2seq lookup for every (cluster, state) pair.
        let mut cache_block: i64 = -1;
        if self.cache_seq.len() < n_comp_haps {
            self.cache_seq.resize(n_comp_haps, 0);
        }
        self.dirty.clear();
        for m in 0..self.n_clusters {
            let targ_allele = imp_data.allele(m, shifted_targ_hap);
            let row = m * words_per_row;
            while ev < events.len() && events[ev].0 as usize == m {
                let (_, j, new_hap) = events[ev];
                self.comp_hap_to_hap[j as usize] = new_hap;
                self.dirty.push(j);
                ev += 1;
            }
            let haps = &self.comp_hap_to_hap[..n_comp_haps];
            let match_out = &mut al_match[row..row + words_per_row];
            match &imp_data.hap_to_seq[m] {
                crate::impdata::ClusterCoding::Composed {
                    block,
                    hap2seq,
                    seq1_to_seq2,
                    targ,
                    ..
                } => {
                    if *block as i64 != cache_block {
                        for j in 0..n_comp_haps {
                            let h = haps[j] as usize;
                            if h < n_ref {
                                self.cache_seq[j] = hap2seq[h];
                            }
                        }
                        cache_block = *block as i64;
                    } else if !self.dirty.is_empty() {
                        for &j in &self.dirty {
                            let h = haps[j as usize] as usize;
                            if h < n_ref {
                                self.cache_seq[j as usize] = hap2seq[h];
                            }
                        }
                    }
                    self.dirty.clear();
                    let cache_seq = &self.cache_seq;
                    let mut w = 0u64;
                    let mut word_idx = 0;
                    for j in 0..n_comp_haps {
                        let h = haps[j] as usize;
                        let v = if h < n_ref {
                            seq1_to_seq2[cache_seq[j] as usize]
                        } else {
                            targ[h - n_ref]
                        };
                        if v == targ_allele {
                            w |= 1u64 << (j & 63);
                        }
                        if j & 63 == 63 {
                            match_out[word_idx] = w;
                            word_idx += 1;
                            w = 0;
                        }
                    }
                    if n_comp_haps & 63 != 0 {
                        match_out[word_idx] = w;
                    }
                }
                crate::impdata::ClusterCoding::Direct {
                    coded_ref, targ, ..
                } => {
                    fill_match_bits(haps, match_out, |h| {
                        let v = if h < n_ref {
                            coded_ref[h]
                        } else {
                            targ[h - n_ref]
                        };
                        v == targ_allele
                    });
                }
            }
        }
        self.events = events;
        n_comp_haps
    }

    /// Resets each composite haplotype to its first segment.
    fn reset_comp_haps(&mut self, n_comp_haps: usize) {
        for j in 0..n_comp_haps {
            self.comp_hap_to_hap[j] = self.comp_hap_hap[j][0];
        }
    }

    /// Precomputes the cluster at which each composite haplotype switches to
    /// its next segment, replacing the per-cluster scan over all states with
    /// a cursor over the (far shorter) list of actual transitions.
    ///
    /// Java advances a composite haplotype at most once per cluster, testing
    /// only its current segment end, so a segment whose end is not strictly
    /// past the previous transition can never fire and pins the composite
    /// haplotype on its current segment for the rest of the window; the
    /// `end <= last_fired` break reproduces that exactly.
    fn build_events(&mut self, n_comp_haps: usize) {
        let mut events = std::mem::take(&mut self.events);
        events.clear();
        let n_clusters = self.n_clusters as i64;
        for j in 0..n_comp_haps {
            let haps = &self.comp_hap_hap[j];
            let ends = &self.comp_hap_end[j];
            let mut idx = 0usize;
            let mut last_fired: i64 = -1;
            loop {
                let end = ends[idx] as i64;
                if end >= n_clusters || end <= last_fired {
                    break;
                }
                idx += 1;
                events.push((end as u32, j as u32, haps[idx]));
                last_fired = end;
            }
        }
        events.sort_unstable_by_key(|e| e.0);
        self.events = events;
    }

    /// Re-walks the composite-haplotype segments from the last `ibs_states`
    /// call, handing the per-cluster state haplotypes to `f(m, haps)`.
    pub fn replay<F: FnMut(usize, &[i32])>(&mut self, n_comp_haps: usize, mut f: F) {
        self.reset_comp_haps(n_comp_haps);
        let events = std::mem::take(&mut self.events);
        let mut ev = 0usize;
        for m in 0..self.n_clusters {
            while ev < events.len() && events[ev].0 as usize == m {
                let (_, j, new_hap) = events[ev];
                self.comp_hap_to_hap[j as usize] = new_hap;
                ev += 1;
            }
            f(m, &self.comp_hap_to_hap[..n_comp_haps]);
        }
        self.events = events;
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
