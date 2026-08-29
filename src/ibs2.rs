//! Port of `phase.Ibs2`, `phase.Ibs2Markers`, `phase.Ibs2Sets`,
//! and `phase.SampleSeg`.

use crate::phasedata::MarkerMap;
use rayon::prelude::*;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct SampleSeg {
    pub sample: u32,
    pub start: u32,
    pub incl_end: u32,
}

const MIN_IBS2_CM: f64 = 2.0;
const MAX_IBD_GAP_CM: f64 = 4.0;

pub struct Ibs2 {
    #[allow(dead_code)]
    n_markers: usize,
    sample_segs: Vec<Vec<SampleSeg>>, // [targ sample][segment]
}

#[allow(dead_code)]
impl Ibs2 {
    pub fn new(targ: &[Vec<i16>], map: &MarkerMap, maf: &[f32]) -> Ibs2 {
        let n_markers = targ.len();
        let n_samples = if n_markers == 0 {
            0
        } else {
            targ[0].len() >> 1
        };
        let ibs2_markers = Ibs2Markers::new(targ, map, maf);
        let ibs2_sets = Ibs2Sets::new(targ, &ibs2_markers, n_samples);
        let gen_pos = &map.gen_pos;
        let sample_segs: Vec<Vec<SampleSeg>> = (0..n_samples)
            .into_par_iter()
            .map(|s| ibs2_segments(targ, gen_pos, &ibs2_sets, s, n_markers))
            .collect();
        Ibs2 {
            n_markers,
            sample_segs,
        }
    }

    /// `Ibs2.areIbs2(targSample, otherSample, start, inclEnd)`
    pub fn are_ibs2_range(
        &self,
        targ_sample: usize,
        other_sample: usize,
        start: usize,
        incl_end: usize,
    ) -> bool {
        debug_assert!(start <= incl_end);
        let same = targ_sample == other_sample;
        if self.sample_segs[targ_sample].is_empty() || same {
            return same;
        }
        for ss in &self.sample_segs[targ_sample] {
            if ss.sample as usize == other_sample
                && start <= ss.incl_end as usize
                && ss.start as usize <= incl_end
            {
                return true;
            }
        }
        false
    }

    pub fn n_markers(&self) -> usize {
        self.n_markers
    }
}

fn ibs2_segments(
    targ: &[Vec<i16>],
    gen_pos: &[f64],
    ibs2_sets: &Ibs2Sets,
    sample: usize,
    n_markers: usize,
) -> Vec<SampleSeg> {
    let mut seg_list = ibs2_sets.seg_list(sample);
    // sampleComp: by sample, then start, then inclEnd
    seg_list.sort_by(|a, b| {
        a.sample
            .cmp(&b.sample)
            .then(a.start.cmp(&b.start))
            .then(a.incl_end.cmp(&b.incl_end))
    });
    let merged = merge_segments(&seg_list, gen_pos);
    let extended: Vec<SampleSeg> = merged
        .iter()
        .map(|ss| extend(targ, sample, *ss, n_markers))
        .collect();
    let merged = merge_segments(&extended, gen_pos);
    merged
        .into_iter()
        .filter(|ss| gen_pos[ss.incl_end as usize] - gen_pos[ss.start as usize] >= MIN_IBS2_CM)
        .collect()
}

fn merge_segments(list: &[SampleSeg], gen_pos: &[f64]) -> Vec<SampleSeg> {
    if list.len() < 2 {
        return list.to_vec();
    }
    let mut merged = Vec::new();
    let mut prev = list[0];
    for &next in &list[1..] {
        if prev.sample == next.sample
            && gen_pos[next.start as usize] - gen_pos[prev.incl_end as usize] <= MAX_IBD_GAP_CM
        {
            prev = SampleSeg {
                sample: prev.sample,
                start: prev.start,
                incl_end: next.incl_end,
            };
        } else {
            merged.push(prev);
            prev = next;
        }
    }
    merged.push(prev);
    merged
}

fn extend(targ: &[Vec<i16>], sample: usize, ss: SampleSeg, n_markers: usize) -> SampleSeg {
    let sample2 = ss.sample as usize;
    let mut incl_start = ss.start as usize;
    let mut excl_end = ss.incl_end as usize + 1;
    while incl_start > 0 && ibs2(targ, incl_start - 1, sample, sample2) {
        incl_start -= 1;
    }
    while excl_end < n_markers && ibs2(targ, excl_end, sample, sample2) {
        excl_end += 1;
    }
    SampleSeg {
        sample: ss.sample,
        start: incl_start as u32,
        incl_end: (excl_end - 1) as u32,
    }
}

fn ibs2(targ: &[Vec<i16>], m: usize, s1: usize, s2: usize) -> bool {
    let row = &targ[m];
    let h1 = s1 << 1;
    let h2 = s2 << 1;
    let a1 = row[h1];
    let a2 = row[h1 | 1];
    let b1 = row[h2];
    let b2 = row[h2 | 1];
    phase_consistent(a1, a2, b1, b2) || phase_consistent(a1, a2, b2, b1)
}

#[inline]
fn phase_consistent(a1: i16, a2: i16, b1: i16, b2: i16) -> bool {
    (a1 < 0 || b1 < 0 || a1 == b1) && (a2 < 0 || b2 < 0 || a2 == b2)
}

// ---------------------------------------------------------------------------
// Ibs2Markers

struct Ibs2Markers {
    use_marker: Vec<bool>,
    step_starts: Vec<usize>,
}

const MAX_MISS_FREQ: f64 = 0.1;
const MIN_MINOR_FREQ: f32 = 0.1;
const MIN_MARKER_CNT: usize = 50;
const MIN_INTERMARKER_CM: f64 = 0.02;

impl Ibs2Markers {
    fn new(targ: &[Vec<i16>], map: &MarkerMap, maf: &[f32]) -> Ibs2Markers {
        let n_markers = targ.len();
        let n_haps = if n_markers == 0 { 0 } else { targ[0].len() };
        let max_miss = (MAX_MISS_FREQ * n_haps as f64).ceil() as usize;
        let mut use_marker: Vec<bool> = (0..n_markers)
            .into_par_iter()
            .map(|m| {
                if maf[m] >= MIN_MINOR_FREQ {
                    let miss_cnt = targ[m].iter().filter(|&&a| a < 0).count();
                    miss_cnt <= max_miss
                } else {
                    false
                }
            })
            .collect();
        let step_starts = step_starts(&mut use_marker, map);
        Ibs2Markers {
            use_marker,
            step_starts,
        }
    }

    fn markers(&self, start: usize, end: usize) -> Vec<usize> {
        (start..end).filter(|&m| self.use_marker[m]).collect()
    }
}

fn step_starts(use_marker: &mut [bool], map: &MarkerMap) -> Vec<usize> {
    let gen_pos = &map.gen_pos;
    let n_markers = gen_pos.len();
    let mut indices = Vec::new();
    let mut last_start = 0usize;
    let mut next = next_start(gen_pos, last_start, use_marker);
    // following code combines the last two steps
    while next < n_markers {
        indices.push(last_start);
        last_start = next;
        next = next_start(gen_pos, next, use_marker);
    }
    indices
}

fn next_start(gen_pos: &[f64], start: usize, use_marker: &mut [bool]) -> usize {
    let mut cm_pos = gen_pos[start];
    let mut min_cm_pos = cm_pos + MIN_INTERMARKER_CM;
    let mut next_start = start + 1;
    let mut mkr_cnt = 0usize;
    while next_start < use_marker.len() && mkr_cnt < MIN_MARKER_CNT {
        if use_marker[next_start] {
            cm_pos = gen_pos[next_start];
            if cm_pos < min_cm_pos {
                use_marker[next_start] = false;
            } else {
                mkr_cnt += 1;
                min_cm_pos = cm_pos + MIN_INTERMARKER_CM;
            }
        }
        next_start += 1;
    }
    next_start
}

// ---------------------------------------------------------------------------
// Ibs2Sets

struct Ibs2Sets {
    n_markers_m1: usize,
    window_starts: Vec<usize>,
    ibs2_sets: Vec<Vec<Vec<u32>>>, // [window][targ_sample][ibs2 samples]
}

const MAX_MISS_STEP_FREQ: f64 = 0.1;

impl Ibs2Sets {
    fn new(targ: &[Vec<i16>], ibs2_markers: &Ibs2Markers, n_samples: usize) -> Ibs2Sets {
        let n_markers = targ.len();
        let window_starts = ibs2_markers.step_starts.clone();
        let ibs2_sets: Vec<Vec<Vec<u32>>> = (0..window_starts.len())
            .into_par_iter()
            .map(|w| {
                let start = window_starts[w];
                let end = if w + 1 < window_starts.len() {
                    window_starts[w + 1]
                } else {
                    n_markers
                };
                let step_markers = ibs2_markers.markers(start, end);
                ibs2_sets_for_window(targ, &step_markers, n_samples)
            })
            .collect();
        Ibs2Sets {
            n_markers_m1: n_markers - 1,
            window_starts,
            ibs2_sets,
        }
    }

    fn seg_list(&self, sample: usize) -> Vec<SampleSeg> {
        let mut list = Vec::new();
        for w in 0..self.ibs2_sets.len() {
            let ia = &self.ibs2_sets[w][sample];
            if !ia.is_empty() {
                let start = self.window_starts[w];
                let incl_end = if w + 1 < self.window_starts.len() {
                    self.window_starts[w + 1] - 1
                } else {
                    self.n_markers_m1
                };
                for &s2 in ia {
                    if s2 as usize != sample {
                        list.push(SampleSeg {
                            sample: s2,
                            start: start as u32,
                            incl_end: incl_end as u32,
                        });
                    }
                }
            }
        }
        list
    }
}

struct SampClust {
    samples: Vec<u32>,
    is_homozygous: bool,
}

fn ibs2_sets_for_window(
    targ: &[Vec<i16>],
    step_markers: &[usize],
    n_samples: usize,
) -> Vec<Vec<u32>> {
    // initCluster: exclude samples with too many missing genotypes
    let mut miss_cnt = vec![0usize; n_samples];
    for &m in step_markers {
        let row = &targ[m];
        for (s, mc) in miss_cnt.iter_mut().enumerate() {
            let h1 = s << 1;
            if row[h1] == -1 || row[h1 | 1] == -1 {
                *mc += 1;
            }
        }
    }
    let max_miss = (MAX_MISS_STEP_FREQ * step_markers.len() as f64).floor() as usize;
    let init: Vec<u32> = (0..n_samples as u32)
        .filter(|&s| miss_cnt[s as usize] <= max_miss)
        .collect();
    let mut partition = vec![SampClust {
        samples: init,
        is_homozygous: true,
    }];
    for &m in step_markers {
        let mut next_partition = Vec::with_capacity(partition.len());
        for parent in &partition {
            partition_cluster(targ, parent, m, &mut next_partition);
        }
        partition = next_partition;
    }
    // results
    let empty: Vec<u32> = Vec::new();
    let mut results: Vec<Vec<u32>> = vec![empty; n_samples];
    for clust in &partition {
        if !clust.is_homozygous {
            debug_assert!(clust.samples.len() > 1);
            for &s in &clust.samples {
                let entry = &mut results[s as usize];
                if entry.is_empty() {
                    *entry = clust.samples.clone();
                } else {
                    // sample can be in >1 list due to missing genotypes
                    let mut combined: Vec<u32> = entry
                        .iter()
                        .chain(clust.samples.iter())
                        .copied()
                        .collect();
                    combined.sort_unstable();
                    combined.dedup();
                    *entry = combined;
                }
            }
        }
    }
    results
}

fn partition_cluster(
    targ: &[Vec<i16>],
    parent: &SampClust,
    m: usize,
    out: &mut Vec<SampClust>,
) {
    let row = &targ[m];
    // number of possible genotypes given the alleles present is bounded by
    // observed genotype indices; use a map keyed by gt index
    let mut gt_to_list: std::collections::HashMap<i32, usize> = std::collections::HashMap::new();
    let mut lists: Vec<(i32, Vec<u32>)> = Vec::new();
    let mut missing: Vec<u32> = Vec::new();
    for &s in &parent.samples {
        let h1 = (s as usize) << 1;
        let a1 = row[h1];
        let a2 = row[h1 | 1];
        let gt_index: i32 = if a1 < 0 || a2 < 0 {
            -1
        } else if a1 <= a2 {
            (((a2 as i32) * (a2 as i32 + 1)) >> 1) + a1 as i32
        } else {
            (((a1 as i32) * (a1 as i32 + 1)) >> 1) + a2 as i32
        };
        if gt_index < 0 {
            missing.push(s);
            for (_, list) in lists.iter_mut() {
                list.push(s);
            }
        } else {
            match gt_to_list.get(&gt_index) {
                Some(&idx) => lists[idx].1.push(s),
                None => {
                    let mut list = missing.clone();
                    list.push(s);
                    gt_to_list.insert(gt_index, lists.len());
                    lists.push((gt_index, list));
                }
            }
        }
    }
    // stream order in Java: gtToList array indexed by gt index ascending
    lists.sort_by_key(|(gt, _)| *gt);
    for (gt, list) in lists {
        if list.len() > 1 {
            let a_hom = is_hom_gt(gt);
            out.push(SampClust {
                samples: list,
                is_homozygous: parent.is_homozygous && a_hom,
            });
        }
    }
}

/// true iff the genotype index corresponds to a homozygous genotype
#[inline]
fn is_hom_gt(gt: i32) -> bool {
    // gt index for (a,a) is a*(a+1)/2 + a = a*(a+3)/2
    let mut a = 0i32;
    loop {
        let hom = (a * (a + 3)) >> 1;
        if hom == gt {
            return true;
        }
        if hom > gt {
            return false;
        }
        a += 1;
    }
}
