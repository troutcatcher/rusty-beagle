//! Ports of `phase.PbwtRecPhaser`, `phase.FwdPbwtPhaser`,
//! `phase.RevPbwtPhaser`, and `phase.PbwtPhaser` (initial PBWT phasing).

use crate::javautil::JavaRandom;
use crate::par::Par;
use crate::pbwt::PbwtUpdater;
use crate::phasedata::{EstPhase, FixedPhaseData, SamplePhase};
use rayon::prelude::*;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// PbwtRecPhaser

struct PbwtRecPhaser<'a> {
    fpd: &'a FixedPhaseData,
    phased_overlap: usize,
    n_targ_haps: usize,
    n_targ_samples: usize,
    #[allow(dead_code)]
    n_haps: usize,
    a: Vec<i32>,
    inv_a: Vec<i32>,
    pbwt: PbwtUpdater,
    ref_buf: Vec<u8>,
}

impl<'a> PbwtRecPhaser<'a> {
    fn new(fpd: &'a FixedPhaseData) -> PbwtRecPhaser<'a> {
        let n_targ_haps = fpd.n_targ_haps;
        let n_haps = fpd.n_haps;
        PbwtRecPhaser {
            fpd,
            phased_overlap: fpd.stage1_overlap,
            n_targ_haps,
            n_targ_samples: n_targ_haps >> 1,
            n_haps,
            a: (0..n_haps as i32).collect(),
            inv_a: (0..n_haps as i32).collect(),
            pbwt: PbwtUpdater::new(n_haps),
            ref_buf: vec![0; fpd.n_ref_haps],
        }
    }

    /// `PbwtRecPhaser.phase(currentMkr, alleles, nextMkr, missing, unphHet)`
    fn phase(
        &mut self,
        current_mkr: isize,
        alleles: &mut [i32],
        next_mkr: usize,
        missing: &mut [bool],
        unph_het: &mut [bool],
    ) -> Vec<i32> {
        if current_mkr != -1 {
            let n_alleles = self.fpd.stage1_markers[current_mkr as usize].n_alleles as usize;
            self.pbwt.update(alleles, n_alleles, &mut self.a);
        }
        let al_cnts = self.set_alleles(next_mkr, alleles, unph_het, missing);
        if next_mkr >= self.phased_overlap {
            self.phase_inner(alleles, unph_het);
        }
        al_cnts
    }

    fn phase_inner(&mut self, alleles: &mut [i32], unph_het: &mut [bool]) {
        for j in 0..self.a.len() {
            self.inv_a[self.a[j] as usize] = j as i32;
        }
        let mut threshold = 2i32;
        let mut change_made = true;
        while threshold > 0 || change_made {
            change_made = false;
            for s in 0..self.n_targ_samples {
                if unph_het[s] {
                    change_made |= self.phase_sample(s, threshold, alleles, unph_het);
                } else {
                    let h1 = s << 1;
                    let h2 = h1 | 1;
                    if alleles[h1] == -1 {
                        alleles[h1] =
                            self.impute(alleles, unph_het, self.inv_a[h1] as usize);
                        change_made |= alleles[h1] >= 0;
                    }
                    if alleles[h2] == -1 {
                        alleles[h2] =
                            self.impute(alleles, unph_het, self.inv_a[h2] as usize);
                        change_made |= alleles[h2] >= 0;
                    }
                }
            }
            if !change_made {
                threshold -= 1;
            }
        }
    }

    fn phase_sample(
        &self,
        s: usize,
        threshold: i32,
        alleles: &mut [i32],
        unph_het: &mut [bool],
    ) -> bool {
        let h1 = s << 1;
        let h2 = h1 | 1;
        let a1 = alleles[h1];
        let a2 = alleles[h2];
        debug_assert!(a1 >= 0 && a2 >= 0 && a1 != a2);
        let cnt1 = self.phase_cnt(alleles, unph_het, self.inv_a[h1] as usize, a1, a2);
        let cnt2 = self.phase_cnt(alleles, unph_het, self.inv_a[h2] as usize, a2, a1);
        let cnt = cnt1 + cnt2;
        if cnt >= threshold {
            unph_het[s] = false;
            return true;
        }
        if cnt <= -threshold {
            alleles[h1] = a2;
            alleles[h2] = a1;
            unph_het[s] = false;
            return true;
        }
        false
    }

    fn phase_cnt(
        &self,
        alleles: &[i32],
        unphased_het: &[bool],
        ai: usize,
        a1: i32,
        a2: i32,
    ) -> i32 {
        let mut phase_cnt = 0;
        if ai > 0 {
            let h = self.a[ai - 1] as usize;
            let s = h >> 1;
            if s >= unphased_het.len() || !unphased_het[s] {
                phase_cnt += adj_cnt(alleles[h], a1, a2);
            }
        }
        if ai + 1 < alleles.len() {
            let h = self.a[ai + 1] as usize;
            let s = h >> 1;
            if s >= unphased_het.len() || !unphased_het[s] {
                phase_cnt += adj_cnt(alleles[h], a1, a2);
            }
        }
        phase_cnt
    }

    fn impute(&self, input_alleles: &[i32], unphased_het: &[bool], ai: usize) -> i32 {
        let mut prev = -1;
        let mut next = -1;
        if ai > 0 {
            let h = self.a[ai - 1] as usize;
            let s = h >> 1;
            if s >= unphased_het.len() || !unphased_het[s] {
                prev = input_alleles[h];
            }
        }
        if ai + 1 < self.a.len() {
            let h = self.a[ai + 1] as usize;
            let s = h >> 1;
            if s >= unphased_het.len() || !unphased_het[s] {
                next = input_alleles[h];
            }
        }
        if prev >= 0 && (prev == next || next < 0) {
            prev
        } else if prev < 0 && next >= 0 {
            next
        } else {
            -1
        }
    }

    /// `PbwtRecPhaser.setAlleles`: returns allele-count CDF
    fn set_alleles(
        &mut self,
        m: usize,
        input_alleles: &mut [i32],
        unph_het: &mut [bool],
        missing: &mut [bool],
    ) -> Vec<i32> {
        let n_alleles = self.fpd.stage1_markers[m].n_alleles as usize;
        let mut al_cnts = vec![0i32; n_alleles];
        let row = &self.fpd.stage1_targ[m];
        for s in 0..self.n_targ_samples {
            let h1 = s << 1;
            let h2 = h1 | 1;
            let a1 = row[h1] as i32;
            let a2 = row[h2] as i32;
            input_alleles[h1] = a1;
            input_alleles[h2] = a2;
            unph_het[s] = m >= self.phased_overlap && (a1 >= 0 && a2 >= 0 && a1 != a2);
            missing[s] = a1 < 0 || a2 < 0;
            if a1 >= 0 {
                al_cnts[a1 as usize] += 1;
            }
            if a2 >= 0 {
                al_cnts[a2 as usize] += 1;
            }
        }
        if !self.fpd.stage1_ref.is_empty() {
            let rec = &self.fpd.stage1_ref[m];
            rec.fill_alleles(&mut self.ref_buf);
            for (i, &al) in self.ref_buf.iter().enumerate() {
                input_alleles[self.n_targ_haps + i] = al as i32;
                al_cnts[al as usize] += 1;
            }
        }
        for j in 1..al_cnts.len() {
            al_cnts[j] += al_cnts[j - 1];
        }
        al_cnts
    }
}

#[inline]
fn adj_cnt(adjacent_allele: i32, a1: i32, a2: i32) -> i32 {
    if adjacent_allele == a1 {
        1
    } else if adjacent_allele == a2 {
        -1
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Per-marker packed phased alleles (Java stores BitArray[nTargHaps*bits])

struct PackedAlleles {
    #[allow(dead_code)]
    bits_per_allele: Vec<u8>,
    rows: Vec<Vec<u8>>, // [marker - start][targ hap] allele (simplified storage)
}

impl PackedAlleles {
    #[inline]
    fn allele(&self, index: usize, hap: usize) -> i32 {
        self.rows[index][hap] as i32
    }
}

// ---------------------------------------------------------------------------
// RevPbwtPhaser

struct RevPbwtPhaser {
    start: usize,
    packed: PackedAlleles,
}

impl RevPbwtPhaser {
    fn new(fpd: &FixedPhaseData, start: usize, end: usize, seed: i64) -> RevPbwtPhaser {
        let mut rand = JavaRandom::new(seed);
        let overlap = fpd.stage1_overlap;
        let n_targ_samples = fpd.n_targ_haps >> 1;
        let mut rec_phaser = PbwtRecPhaser::new(fpd);
        let mut missing_gt = vec![false; n_targ_samples];
        let mut unph_het = vec![false; n_targ_samples];
        let mut alleles = vec![0i32; fpd.n_haps];
        let mut rows: Vec<Vec<u8>> = vec![Vec::new(); end - start];
        let mut last_m: isize = -1;
        for m in (start..end).rev() {
            let allele_cdf =
                rec_phaser.phase(last_m, &mut alleles, m, &mut missing_gt, &mut unph_het);
            if m >= overlap {
                Self::finish_phasing(&mut alleles, &mut unph_het, &allele_cdf, &mut rand);
            }
            rows[m - start] = alleles[..fpd.n_targ_haps]
                .iter()
                .map(|&a| a as u8)
                .collect();
            last_m = m as isize;
        }
        let bits_per_allele: Vec<u8> = (start..end)
            .map(|m| {
                let n = fpd.stage1_markers[m].n_alleles as u32;
                (32 - (n - 1).leading_zeros().min(32)) as u8
            })
            .collect();
        RevPbwtPhaser {
            start,
            packed: PackedAlleles {
                bits_per_allele,
                rows,
            },
        }
    }

    fn finish_phasing(
        alleles: &mut [i32],
        unph_het: &mut [bool],
        allele_cdf: &[i32],
        rand: &mut JavaRandom,
    ) {
        for s in 0..unph_het.len() {
            let h1 = s << 1;
            let h2 = h1 | 1;
            if unph_het[s] {
                let a1 = alleles[h1];
                let a2 = alleles[h2];
                if rand.next_boolean() {
                    alleles[h1] = a2;
                    alleles[h2] = a1;
                }
                unph_het[s] = false;
            } else {
                if alleles[h1] == -1 {
                    alleles[h1] = impute_allele(allele_cdf, rand);
                }
                if alleles[h2] == -1 {
                    alleles[h2] = impute_allele(allele_cdf, rand);
                }
            }
        }
    }

    #[inline]
    fn allele(&self, marker: usize, hap: usize) -> i32 {
        self.packed.allele(marker - self.start, hap)
    }

    #[inline]
    #[allow(dead_code)]
    fn bits_per_allele(&self, marker: usize) -> u8 {
        self.packed.bits_per_allele[marker - self.start]
    }
}

fn impute_allele(allele_cdf: &[i32], rand: &mut JavaRandom) -> i32 {
    let bound = allele_cdf[allele_cdf.len() - 1];
    if bound == 0 {
        0
    } else {
        let r = rand.next_int_bound(bound);
        let mut allele = 0usize;
        while r >= allele_cdf[allele] {
            allele += 1;
        }
        allele as i32
    }
}

// ---------------------------------------------------------------------------
// FwdPbwtPhaser

struct FwdPbwtPhaser {
    start: usize,
    #[allow(dead_code)]
    end: usize,
    packed: PackedAlleles,
}

impl FwdPbwtPhaser {
    fn new(fpd: &FixedPhaseData, start: usize, end: usize, seed: i64) -> FwdPbwtPhaser {
        let overlap = fpd.stage1_overlap;
        let n_targ_samples = fpd.n_targ_haps >> 1;
        let mut rec_phaser = PbwtRecPhaser::new(fpd);
        let rev_pbwt = RevPbwtPhaser::new(fpd, start, end, seed);
        let mut missing_gt = vec![false; n_targ_samples];
        let mut unph_het = vec![false; n_targ_samples];
        let mut last_het = vec![-1isize; n_targ_samples];
        let mut alleles = vec![0i32; fpd.n_haps];

        let mut rows: Vec<Vec<u8>> = vec![Vec::new(); end - start];
        let mut last_m: isize = -1;
        for m in start..end {
            rec_phaser.phase(last_m, &mut alleles, m, &mut missing_gt, &mut unph_het);
            if m >= overlap {
                Self::finish_phasing(
                    &rows, &rev_pbwt, start, m, &mut alleles, &last_het, &mut unph_het,
                );
            }
            rows[m - start] = alleles[..fpd.n_targ_haps]
                .iter()
                .map(|&a| a as u8)
                .collect();
            for s in 0..n_targ_samples {
                let h1 = s << 1;
                if !missing_gt[s] && alleles[h1] != alleles[h1 | 1] {
                    last_het[s] = m as isize;
                }
            }
            last_m = m as isize;
        }
        let bits_per_allele: Vec<u8> = (start..end)
            .map(|m| {
                let n = fpd.stage1_markers[m].n_alleles as u32;
                (32 - (n - 1).leading_zeros().min(32)) as u8
            })
            .collect();
        FwdPbwtPhaser {
            start,
            end,
            packed: PackedAlleles {
                bits_per_allele,
                rows,
            },
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_phasing(
        rows: &[Vec<u8>],
        rev_pbwt: &RevPbwtPhaser,
        start: usize,
        m: usize,
        alleles: &mut [i32],
        last_het: &[isize],
        unph_het: &mut [bool],
    ) {
        for s in 0..unph_het.len() {
            let h1 = s << 1;
            let h2 = h1 | 1;
            if unph_het[s] {
                let prev_het = last_het[s];
                if prev_het >= 0 {
                    let prev_het = prev_het as usize;
                    let a1 = rev_pbwt.allele(prev_het, h1);
                    let a2 = rev_pbwt.allele(prev_het, h2);
                    let b1 = rev_pbwt.allele(m, h1);
                    let b2 = rev_pbwt.allele(m, h2);
                    let rev_same_phase = (a1 < a2) == (b1 < b2);
                    let c1 = rows[prev_het - start][h1] as i32;
                    let c2 = rows[prev_het - start][h2] as i32;
                    let fwd_same_phase = (c1 < c2) == (alleles[h1] < alleles[h2]);
                    if rev_same_phase != fwd_same_phase {
                        alleles.swap(h1, h2);
                    }
                }
                unph_het[s] = false;
            } else {
                if alleles[h1] == -1 {
                    alleles[h1] =
                        Self::impute_allele(rows, rev_pbwt, start, last_het[s], m, h1);
                }
                if alleles[h2] == -1 {
                    alleles[h2] =
                        Self::impute_allele(rows, rev_pbwt, start, last_het[s], m, h2);
                }
            }
        }
    }

    fn impute_allele(
        rows: &[Vec<u8>],
        rev_pbwt: &RevPbwtPhaser,
        start: usize,
        last_het: isize,
        m: usize,
        hap: usize,
    ) -> i32 {
        if last_het < 0 {
            return rev_pbwt.allele(m, hap);
        }
        let last_het = last_het as usize;
        let comp_hap = hap ^ 1;
        let a1 = rev_pbwt.allele(last_het, hap);
        let a2 = rev_pbwt.allele(last_het, comp_hap);
        let b1 = rows[last_het - start][hap] as i32;
        let b2 = rows[last_het - start][comp_hap] as i32;
        if (a1 < a2) == (b1 < b2) {
            rev_pbwt.allele(m, hap)
        } else {
            rev_pbwt.allele(m, comp_hap)
        }
    }

    #[inline]
    fn allele(&self, marker: usize, hap: usize) -> i32 {
        self.packed.allele(marker - self.start, hap)
    }
}

// ---------------------------------------------------------------------------
// PbwtPhaser (initPhase)

/// Port of `PbwtPhaser.initPhase`.
pub fn init_phase(fpd: &Arc<FixedPhaseData>, par: &Par, seed: i64) -> EstPhase {
    let windows = hi_freq_windows(fpd, par);
    let phasers: Vec<FwdPbwtPhaser> = windows
        .par_iter()
        .enumerate()
        .map(|(j, &(start, end))| {
            FwdPbwtPhaser::new(fpd, start, end, seed.wrapping_add(j as i64))
        })
        .collect();

    let n_samples = fpd.n_targ_samples;
    let n_threads = par.nthreads;
    let max_step_size = 128usize;
    let step_size = ((n_samples + n_threads - 1) / n_threads).min(max_step_size);
    let n_steps = (n_samples + step_size - 1) / step_size;
    let mut phase: Vec<std::sync::Mutex<Option<SamplePhase>>> =
        (0..n_samples).map(|_| std::sync::Mutex::new(None)).collect();
    let results: Vec<(usize, SamplePhase)> = (0..n_steps)
        .into_par_iter()
        .flat_map_iter(|step| set_sample_phase(fpd, &phasers, step, step_size))
        .collect();
    for (s, sp) in results {
        phase[s] = std::sync::Mutex::new(Some(sp));
    }
    EstPhase {
        fpd: fpd.clone(),
        phase,
    }
}

struct SampleIndices {
    miss_indices: Vec<usize>,
    het_indices: Vec<usize>,
}

fn set_sample_phase(
    fpd: &FixedPhaseData,
    ppa: &[FwdPbwtPhaser],
    step: usize,
    step_size: usize,
) -> Vec<(usize, SamplePhase)> {
    let n_markers = fpd.n_stage1_markers();
    let s_start = step * step_size;
    let s_end = (s_start + step_size).min(fpd.n_targ_samples);
    let indices = sample_indices(fpd, s_start, s_end);

    let n_local_haps = (s_end - s_start) << 1;
    let mut haps: Vec<Vec<i32>> = vec![vec![0i32; n_markers]; n_local_haps];
    // copyHaps for each phaser window with alignment at overlaps;
    // the row swap (Java's alignedHaps) is scoped to a single phaser window
    let mut overlap_end = 0usize;
    for (pi, phaser) in ppa.iter().enumerate() {
        let mut hap_slot: Vec<usize> = (0..n_local_haps).collect();
        if pi > 0 {
            overlap_end = ppa[pi - 1].end;
        }
        let copy_start = (phaser.start + overlap_end) >> 1;
        if phaser.start > 0 {
            for s in s_start..s_end {
                let ss = s - s_start;
                let hh1 = ss << 1;
                let hh2 = hh1 | 1;
                let align_het = alignment_het(
                    &indices[ss].het_indices,
                    phaser.start,
                    copy_start,
                    overlap_end,
                );
                if let Some(het) = align_het {
                    let h1 = s << 1;
                    let h2 = h1 | 1;
                    let a1 = haps[hap_slot[hh1]][het];
                    let a2 = haps[hap_slot[hh2]][het];
                    let b1 = phaser.allele(het, h1);
                    let b2 = phaser.allele(het, h2);
                    if a1 == b2 && a2 == b1 {
                        hap_slot.swap(hh1, hh2);
                    }
                }
            }
        }
        for m in copy_start..phaser.end {
            for s in s_start..s_end {
                let h1 = s << 1;
                let h2 = h1 | 1;
                let hh1 = (s - s_start) << 1;
                let hh2 = hh1 | 1;
                haps[hap_slot[hh1]][m] = phaser.allele(m, h1);
                haps[hap_slot[hh2]][m] = phaser.allele(m, h2);
            }
        }
    }

    let gen_pos = &fpd.stage1_map.gen_pos;
    (s_start..s_end)
        .map(|s| {
            let ss = s - s_start;
            let hh1 = ss << 1;
            let hh2 = hh1 | 1;
            let sp = SamplePhase::new(
                s,
                &fpd.stage1_layout,
                gen_pos,
                &haps[hh1],
                &haps[hh2],
                &indices[ss].het_indices,
                &indices[ss].miss_indices,
                &fpd.stage1_positions,
            );
            (s, sp)
        })
        .collect()
}

/// `PbwtPhaser.alignmentHet`: returns None if no alignment het exists.
fn alignment_het(
    het_list: &[usize],
    start: usize,
    copy_start: usize,
    overlap_end: usize,
) -> Option<usize> {
    if het_list.is_empty() {
        return None;
    }
    let mut index = het_list.partition_point(|&h| h < copy_start);
    if index == het_list.len() || (het_list[index] >= overlap_end && index > 0) {
        index -= 1;
    }
    let het = het_list[index];
    if start <= het && het < overlap_end {
        Some(het)
    } else {
        None
    }
}

fn sample_indices(fpd: &FixedPhaseData, s_start: usize, s_end: usize) -> Vec<SampleIndices> {
    let overlap = fpd.stage1_overlap;
    let n_markers = fpd.n_stage1_markers();
    let len = s_end - s_start;
    let mut miss: Vec<Vec<usize>> = vec![Vec::new(); len];
    let mut hets: Vec<Vec<usize>> = vec![Vec::new(); len];
    let mut not_first_het = vec![false; len];
    for m in 0..n_markers {
        let row = &fpd.stage1_targ[m];
        for s in s_start..s_end {
            let ss = s - s_start;
            let h1 = s << 1;
            let a1 = row[h1];
            let a2 = row[h1 | 1];
            if a1 < 0 || a2 < 0 {
                miss[ss].push(m);
            } else if a1 != a2 {
                if m >= overlap && not_first_het[ss] {
                    hets[ss].push(m);
                } else {
                    not_first_het[ss] = true;
                }
            }
        }
    }
    miss.into_iter()
        .zip(hets)
        .map(|(miss_indices, het_indices)| SampleIndices {
            miss_indices,
            het_indices,
        })
        .collect()
}

/// `PbwtPhaser.hiFreqWindows`
fn hi_freq_windows(fpd: &FixedPhaseData, par: &Par) -> Vec<(usize, usize)> {
    let gen_pos = &fpd.stage1_map.gen_pos;
    let n_markers = gen_pos.len();
    let n_threads = par.nthreads;
    let total_cm = gen_pos[n_markers - 1] - gen_pos[0];
    let overlap_cm = 0.5f64;
    let advance_cm = (total_cm / n_threads as f64).max(4.0 * overlap_cm);
    let mut windows = Vec::with_capacity(n_threads);
    let mut from = 0usize;
    let mut to = to_index(gen_pos, gen_pos[from] + advance_cm);
    while to < n_markers {
        windows.push((from, to));
        from = from_index(gen_pos, gen_pos[to] - overlap_cm);
        to = to_index(gen_pos, gen_pos[to] + advance_cm);
    }
    windows.push((from, to));
    windows
}

fn from_index(gen_pos: &[f64], pos: f64) -> usize {
    // Java binarySearch: found -> index; else insertion point
    match gen_pos.binary_search_by(|p| p.partial_cmp(&pos).unwrap()) {
        Ok(i) => i,
        Err(i) => i,
    }
}

fn to_index(gen_pos: &[f64], pos: f64) -> usize {
    match gen_pos.binary_search_by(|p| p.partial_cmp(&pos).unwrap()) {
        Ok(i) => i + 1,
        Err(i) => i,
    }
}
