//! Phasing data structures: ports of `vcf.Markers` (bit layout),
//! `vcf.MarkerMap`, `vcf.Steps`, `phase.SamplePhase`, `phase.EstPhase`,
//! `phase.FixedPhaseData`, `phase.PhaseData`, `phase.MarkerCluster`,
//! and `phase.SwapRate`.

use crate::bits::BitArray;
use crate::genmap::GeneticMap;
use crate::javautil::JavaRandom;
use crate::marker::Marker;
use crate::par::Par;
use crate::refpanel::RefRec;
use crate::windows::Window;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Bit layout over a marker list (part of Java's vcf.Markers)

pub struct HapLayout {
    pub bits_per_allele: Vec<u8>,
    pub sum_hap_bits: Vec<u32>, // length nMarkers+1
    pub n_alleles: Vec<u16>,
}

#[allow(dead_code)]
impl HapLayout {
    pub fn new(n_alleles: Vec<u16>) -> HapLayout {
        let bits_per_allele: Vec<u8> = n_alleles
            .iter()
            .map(|&n| (32 - ((n as u32) - 1).leading_zeros().min(32)) as u8)
            .collect();
        let mut sum = Vec::with_capacity(n_alleles.len() + 1);
        let mut acc = 0u32;
        sum.push(0);
        for &b in &bits_per_allele {
            acc += b as u32;
            sum.push(acc);
        }
        HapLayout {
            bits_per_allele,
            sum_hap_bits: sum,
            n_alleles,
        }
    }

    #[inline]
    pub fn n_markers(&self) -> usize {
        self.bits_per_allele.len()
    }

    #[inline]
    pub fn total_bits(&self) -> usize {
        *self.sum_hap_bits.last().unwrap() as usize
    }

    /// `Markers.allele(bits, marker)`
    #[inline]
    pub fn allele(&self, bits: &BitArray, m: usize) -> i32 {
        let start = self.sum_hap_bits[m] as usize;
        let end = self.sum_hap_bits[m + 1] as usize;
        if end == start + 1 {
            return if bits.get(start) { 1 } else { 0 };
        }
        let mut allele = 0i32;
        let mut mask = 1i32;
        for j in start..end {
            if bits.get(j) {
                allele |= mask;
            }
            mask <<= 1;
        }
        allele
    }

    /// `Markers.setAllele(marker, allele, bits)`
    pub fn set_allele(&self, m: usize, allele: i32, bits: &mut BitArray) {
        let start = self.sum_hap_bits[m] as usize;
        let end = self.sum_hap_bits[m + 1] as usize;
        let mut mask = 1i32;
        for j in start..end {
            if allele & mask == mask {
                bits.set(j);
            } else {
                bits.clear_bit(j);
            }
            mask <<= 1;
        }
    }

    /// `Markers.allelesToBits(alleles, bits)` (bits assumed clear)
    pub fn alleles_to_bits(&self, alleles: &[i32], bits: &mut BitArray) {
        for (m, &al) in alleles.iter().enumerate() {
            let start = self.sum_hap_bits[m] as usize;
            let end = self.sum_hap_bits[m + 1] as usize;
            let mut mask = 1i32;
            for j in start..end {
                if al & mask == mask {
                    bits.set(j);
                }
                mask <<= 1;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// MarkerMap / Steps (ports of vcf.MarkerMap, vcf.Steps)

pub struct MarkerMap {
    pub gen_pos: Vec<f64>,
    pub gen_dist: Vec<f32>,
}

impl MarkerMap {
    /// `MarkerMap.create(genMap, markers)` over the given markers.
    pub fn create(genmap: &GeneticMap, markers: &[Arc<Marker>]) -> MarkerMap {
        let a = &markers[0];
        let b = &markers[markers.len() - 1];
        if a.pos == b.pos {
            eprintln!(
                "ERROR: Window has only one position: CHROM={} POS={}",
                a.chrom(),
                a.pos
            );
            std::process::exit(1);
        }
        let mean_single_base = ((genmap.gen_pos_marker(b) - genmap.gen_pos_marker(a)).abs()
            / ((b.pos - a.pos).abs() as f64))
            .max(1e-8);
        // GeneticMap.genPos(genMap, minGenDist, markers)
        let mut gp = vec![0.0f64; markers.len()];
        gp[0] = genmap.gen_pos_marker(&markers[0]);
        let mut last_map_pos = gp[0];
        for j in 1..gp.len() {
            let map_pos = genmap.gen_pos_marker(&markers[j]);
            let dist = (map_pos - last_map_pos).max(mean_single_base);
            gp[j] = gp[j - 1] + dist;
            last_map_pos = map_pos;
        }
        MarkerMap::from_gen_pos(gp)
    }

    pub fn from_gen_pos(gen_pos: Vec<f64>) -> MarkerMap {
        let mut gen_dist = vec![0.0f32; gen_pos.len()];
        for j in 1..gen_pos.len() {
            gen_dist[j] = (gen_pos[j] - gen_pos[j - 1]) as f32;
            if gen_dist[j] <= 0.0 {
                eprintln!("ERROR: Nonpositive genetic distance: dist[{}]={}", j, gen_dist[j]);
                std::process::exit(1);
            }
        }
        MarkerMap { gen_pos, gen_dist }
    }

    pub fn restrict(&self, indices: &[usize]) -> MarkerMap {
        let gp: Vec<f64> = indices.iter().map(|&i| self.gen_pos[i]).collect();
        MarkerMap::from_gen_pos(gp)
    }

    /// `MarkerMap.pRecomb(recombIntensity)`
    pub fn p_recomb(&self, recomb_intensity: f32) -> Vec<f32> {
        let c = -(recomb_intensity as f64);
        self.gen_dist
            .iter()
            .map(|&d| (-f64::exp_m1(c * d as f64)) as f32)
            .collect()
    }
}

/// Port of `vcf.Steps`.
pub struct Steps {
    step_ends: Vec<usize>,
}

impl Steps {
    pub fn new(map: &MarkerMap, min_step: f32) -> Steps {
        let min_step = min_step as f64;
        let gen_pos = &map.gen_pos;
        let n_markers = gen_pos.len();
        let mut indices = Vec::with_capacity(n_markers >> 1);
        let mut end = 0usize;
        while end < n_markers {
            let min_gen_pos = gen_pos[end] + min_step;
            end += 1;
            while end < n_markers && gen_pos[end] < min_gen_pos {
                end += 1;
            }
            indices.push(end);
        }
        Steps { step_ends: indices }
    }

    #[inline]
    pub fn size(&self) -> usize {
        self.step_ends.len()
    }

    #[inline]
    pub fn start(&self, step: usize) -> usize {
        if step == 0 {
            0
        } else {
            self.step_ends[step - 1]
        }
    }

    #[inline]
    pub fn end(&self, step: usize) -> usize {
        self.step_ends[step]
    }
}

// ---------------------------------------------------------------------------
// SamplePhase (port of phase.SamplePhase)

pub const CLUST_MISSING_GT: u8 = 0;
pub const CLUST_MASKED_HET: u8 = 1;
pub const CLUST_HOMOZYGOUS_GT: u8 = 2;
pub const CLUST_PHASED_HET: u8 = 3;
pub const CLUST_UNPHASED_HET: u8 = 4;

pub struct SamplePhase {
    pub sample: usize,
    pub hap1: BitArray,
    pub hap2: BitArray,
    pub clust_size: Vec<u8>,
    pub clust_type: Vec<u8>,
    pub clust_type_cnt: [usize; 5],
}

#[allow(dead_code)]
impl SamplePhase {
    pub fn new(
        sample: usize,
        layout: &HapLayout,
        gen_pos: &[f64],
        hap1: &[i32],
        hap2: &[i32],
        unphased_hets: &[usize],
        missing_gts: &[usize],
        stage1_pos: &[i32],
    ) -> SamplePhase {
        let _ = stage1_pos;
        let n_markers = layout.n_markers();
        debug_assert_eq!(gen_pos.len(), n_markers);
        let mut b1 = BitArray::new(layout.total_bits());
        let mut b2 = BitArray::new(layout.total_bits());
        layout.alleles_to_bits(hap1, &mut b1);
        layout.alleles_to_bits(hap2, &mut b2);
        let max_cluster_cm = 0.005f32;
        let mut clust_type_list: Vec<u8> = Vec::new();
        let mut clust_size_list: Vec<u8> = Vec::new();
        let mut clust_type_cnt = [0usize; 5];
        set_clusters(
            hap1,
            hap2,
            missing_gts,
            unphased_hets,
            gen_pos,
            max_cluster_cm,
            &mut clust_type_list,
            &mut clust_type_cnt,
            &mut clust_size_list,
        );
        SamplePhase {
            sample,
            hap1: b1,
            hap2: b2,
            clust_size: clust_size_list,
            clust_type: clust_type_list,
            clust_type_cnt,
        }
    }

    pub fn n_clusters(&self) -> usize {
        self.clust_size.len()
    }

    pub fn clust_size_at(&self, c: usize) -> usize {
        self.clust_size[c] as usize
    }

    pub fn clust_type_at(&self, c: usize) -> u8 {
        self.clust_type[c]
    }

    pub fn n_unphased(&self) -> usize {
        self.clust_type_cnt[CLUST_UNPHASED_HET as usize]
    }

    pub fn n_masked(&self) -> usize {
        self.clust_type_cnt[CLUST_MASKED_HET as usize]
    }

    pub fn n_missing(&self) -> usize {
        self.clust_type_cnt[CLUST_MISSING_GT as usize]
    }

    pub fn clust_ends(&self) -> Vec<usize> {
        let mut ends = Vec::with_capacity(self.clust_size.len());
        let mut cum = 0usize;
        for &s in &self.clust_size {
            cum += s as usize;
            ends.push(cum);
        }
        ends
    }

    /// `SamplePhase.maskTrailingUnphasedHets` (needs stage1 marker base
    /// positions for the 3000bp rule).
    pub fn mask_trailing_unphased_hets(&mut self, positions: &[i32]) {
        let max_unph_het_clusters = 3usize;
        let max_masked_base_pairs = 3000i32;
        let mut unph_het_markers: Vec<usize> = Vec::new();
        let mut unph_het_clusters: Vec<usize> = Vec::new();
        let mut start_marker = 0usize;
        for c in 0..self.clust_type.len() {
            let ct = self.clust_type[c];
            if ct == CLUST_PHASED_HET {
                if 2 <= unph_het_clusters.len() && unph_het_clusters.len() <= max_unph_het_clusters
                {
                    self.mask_trailing(&unph_het_clusters, &unph_het_markers, max_masked_base_pairs, positions);
                }
                unph_het_markers.clear();
                unph_het_clusters.clear();
            } else if ct == CLUST_UNPHASED_HET {
                unph_het_markers.push(start_marker);
                unph_het_clusters.push(c);
            }
            start_marker += self.clust_size[c] as usize;
        }
        if 2 <= unph_het_clusters.len() && unph_het_clusters.len() <= max_unph_het_clusters {
            self.mask_trailing(&unph_het_clusters, &unph_het_markers, max_masked_base_pairs, positions);
        }
    }

    fn mask_trailing(
        &mut self,
        unph_het_clusters: &[usize],
        unph_het_markers: &[usize],
        max_masked_base_pairs: i32,
        positions: &[i32],
    ) {
        let last_masked_index = unph_het_clusters.len() - 2;
        if last_masked_index == 0 {
            self.mask_het_cluster(unph_het_clusters[0]);
        } else {
            let start_pos = positions[unph_het_markers[0]];
            let end_pos = positions[unph_het_markers[last_masked_index]];
            if end_pos - start_pos <= max_masked_base_pairs {
                for &c in &unph_het_clusters[..=last_masked_index] {
                    self.mask_het_cluster(c);
                }
            }
        }
    }

    pub fn mask_het_cluster(&mut self, cluster: usize) {
        debug_assert_eq!(self.clust_type[cluster], CLUST_UNPHASED_HET);
        self.clust_type[cluster] = CLUST_MASKED_HET;
        self.clust_type_cnt[CLUST_UNPHASED_HET as usize] -= 1;
        self.clust_type_cnt[CLUST_MASKED_HET as usize] += 1;
    }

    pub fn mark_unphased_het_cluster_as_phased(&mut self, cluster: usize) {
        debug_assert_eq!(self.clust_type[cluster], CLUST_UNPHASED_HET);
        self.clust_type[cluster] = CLUST_PHASED_HET;
        self.clust_type_cnt[CLUST_UNPHASED_HET as usize] -= 1;
        self.clust_type_cnt[CLUST_PHASED_HET as usize] += 1;
    }

    pub fn mark_masked_het_cluster_as_phased(&mut self, cluster: usize) {
        debug_assert_eq!(self.clust_type[cluster], CLUST_MASKED_HET);
        self.clust_type[cluster] = CLUST_PHASED_HET;
        self.clust_type_cnt[CLUST_MASKED_HET as usize] -= 1;
        self.clust_type_cnt[CLUST_PHASED_HET as usize] += 1;
    }

    pub fn allele1(&self, layout: &HapLayout, m: usize) -> i32 {
        layout.allele(&self.hap1, m)
    }

    pub fn allele2(&self, layout: &HapLayout, m: usize) -> i32 {
        layout.allele(&self.hap2, m)
    }

    pub fn set_allele1(&mut self, layout: &HapLayout, m: usize, allele: i32) {
        layout.set_allele(m, allele, &mut self.hap1);
    }

    pub fn set_allele2(&mut self, layout: &HapLayout, m: usize, allele: i32) {
        layout.set_allele(m, allele, &mut self.hap2);
    }

    /// `SamplePhase.swapHaps(start, end)`
    pub fn swap_haps(&mut self, layout: &HapLayout, start: usize, end: usize) {
        let start_bit = layout.sum_hap_bits[start] as usize;
        let end_bit = layout.sum_hap_bits[end] as usize;
        BitArray::swap_bits(&mut self.hap1, &mut self.hap2, start_bit, end_bit);
    }
}

#[allow(clippy::too_many_arguments)]
fn set_clusters(
    hap1: &[i32],
    hap2: &[i32],
    missing_gt: &[usize],
    unph_hets: &[usize],
    gen_pos: &[f64],
    max_cm: f32,
    clust_type: &mut Vec<u8>,
    clust_type_cnt: &mut [usize; 5],
    clust_size_list: &mut Vec<u8>,
) {
    let n_markers = gen_pos.len();
    let mut max_clust_end = gen_pos[0] + max_cm as f64;
    let mut prev_is_missing_or_het = false;
    let mut last_end = 0usize;
    let mut miss_index = 0usize;
    let mut unph_index = 0usize;
    let mut next_miss: isize = if miss_index < missing_gt.len() {
        let v = missing_gt[miss_index] as isize;
        miss_index += 1;
        v
    } else {
        -1
    };
    let mut next_unph: isize = if unph_index < unph_hets.len() {
        let v = unph_hets[unph_index] as isize;
        unph_index += 1;
        v
    } else {
        -1
    };
    let mut prev_type = CLUST_HOMOZYGOUS_GT;
    for m in 0..n_markers {
        let size = m - last_end;
        let t = clust_type_of(m as isize == next_miss, m as isize == next_unph, hap1[m], hap2[m]);
        if t == CLUST_MISSING_GT {
            next_miss = if miss_index < missing_gt.len() {
                let v = missing_gt[miss_index] as isize;
                miss_index += 1;
                v
            } else {
                -1
            };
        } else if t == CLUST_UNPHASED_HET {
            next_unph = if unph_index < unph_hets.len() {
                let v = unph_hets[unph_index] as isize;
                unph_index += 1;
                v
            } else {
                -1
            };
        }
        let is_missing_or_het =
            t == CLUST_MISSING_GT || t == CLUST_UNPHASED_HET || t == CLUST_PHASED_HET;
        if is_missing_or_het || prev_is_missing_or_het || gen_pos[m] > max_clust_end || size == 255
        {
            if m > 0 {
                clust_type.push(prev_type);
                clust_type_cnt[prev_type as usize] += 1;
                clust_size_list.push(size as u8);
                max_clust_end = gen_pos[m] + max_cm as f64;
                last_end = m;
            }
            prev_type = t;
        }
        prev_is_missing_or_het = is_missing_or_het;
    }
    clust_type.push(prev_type);
    clust_type_cnt[prev_type as usize] += 1;
    clust_size_list.push((n_markers - last_end) as u8);
}

fn clust_type_of(is_missing: bool, is_unphased: bool, a1: i32, a2: i32) -> u8 {
    if is_missing {
        CLUST_MISSING_GT
    } else if a1 == a2 {
        CLUST_HOMOZYGOUS_GT
    } else if is_unphased {
        CLUST_UNPHASED_HET
    } else {
        CLUST_PHASED_HET
    }
}

// ---------------------------------------------------------------------------
// SwapRate (port of phase.SwapRate)

pub static SWAP_N_SWAPS: AtomicU64 = AtomicU64::new(0);
pub static SWAP_N_UNPH_HETS: AtomicU64 = AtomicU64::new(0);

pub fn swap_rate_increment(n_unph_hets: usize, n_swaps: usize) {
    SWAP_N_SWAPS.fetch_add(n_swaps as u64, Ordering::Relaxed);
    SWAP_N_UNPH_HETS.fetch_add(n_unph_hets as u64, Ordering::Relaxed);
}

pub fn swap_rate_get_and_reset() -> f64 {
    let swaps = SWAP_N_SWAPS.swap(0, Ordering::Relaxed);
    let hets = SWAP_N_UNPH_HETS.swap(0, Ordering::Relaxed);
    swaps as f64 / hets as f64
}

// ---------------------------------------------------------------------------
// Carrier lists (port of vcf.Window.carriers)

#[derive(Clone, PartialEq)]
pub enum Carriers {
    ZeroFreq,
    HighFreq,
    List(Arc<Vec<u32>>), // sample indices: target samples, then ref samples shifted
}

#[allow(dead_code)]
impl Carriers {
    pub fn is_low_freq(&self) -> bool {
        !matches!(self, Carriers::HighFreq)
    }

    pub fn len(&self) -> usize {
        match self {
            Carriers::List(v) => v.len(),
            _ => 0,
        }
    }
}

// ---------------------------------------------------------------------------
// FixedPhaseData (port of phase.FixedPhaseData)

#[allow(dead_code)]
pub struct FixedPhaseData {
    pub window_index: usize,
    pub n_targ_samples: usize,
    pub n_targ_haps: usize,
    pub n_haps: usize,
    pub overlap: usize,

    /// all-window target alleles (with phased overlap spliced in),
    /// marker-major, -1 = missing
    pub targ: Vec<Vec<i16>>,
    pub targ_markers: Vec<Arc<Marker>>,
    /// per-marker, per-allele carrier lists
    pub carriers: Vec<Vec<Carriers>>,
    pub map: MarkerMap,
    /// prefix sums of nAlleles over all target markers
    pub sum_alleles: Vec<u32>,

    // stage 1 data
    pub stage1_to2: Vec<usize>,
    pub stage1_targ: Vec<Vec<i16>>,
    pub stage1_markers: Vec<Arc<Marker>>,
    pub stage1_ref: Vec<Arc<RefRec>>,
    pub stage1_layout: Arc<HapLayout>,
    pub stage1_map: MarkerMap,
    pub ibs_step: f32,
    pub stage1_steps: Steps,
    pub stage1_overlap: usize,
    pub stage1_positions: Vec<i32>,

    pub prev_stage1_marker: Vec<usize>,
    pub prev_stage1_wt: Vec<f32>,

    pub n_ref_haps: usize,
    /// reference records restricted to all window target markers
    pub restrict_ref: Vec<Arc<RefRec>>,
    pub stage1_ibs2: crate::ibs2::Ibs2,
    /// hap-major phased reference genotypes at stage1 markers
    pub stage1_xref: Option<Arc<crate::xref::XRefGT>>,
}

const MAX_HIFREQ_PROP: f32 = 0.75;

impl FixedPhaseData {
    pub fn new(
        par: &Par,
        genmap: &GeneticMap,
        window: &Window,
        phased_overlap: Option<&PhasedOverlap>,
    ) -> FixedPhaseData {
        let n_targ_markers = window.targ_recs.len();
        let targ_samples = window.targ_recs[0].alleles.len() >> 1;
        let n_targ_haps = targ_samples << 1;
        let restrict_ref = window.restrict_ref();
        let n_ref_haps = if restrict_ref.is_empty() {
            0
        } else {
            restrict_ref[0].n_haps
        };
        let n_ref_samples = n_ref_haps >> 1;
        let n_haps = n_targ_haps + n_ref_haps;

        // marker-major target alleles with overlap splice
        let mut targ: Vec<Vec<i16>> = window
            .targ_recs
            .iter()
            .map(|r| r.alleles.clone())
            .collect();
        let overlap = match phased_overlap {
            None => 0,
            Some(ov) => {
                for (m, row) in ov.alleles.iter().enumerate() {
                    targ[m].copy_from_slice(row);
                }
                ov.alleles.len()
            }
        };
        let targ_markers: Vec<Arc<Marker>> = window
            .targ_recs
            .iter()
            .map(|r| Arc::new(r.marker.clone()))
            .collect();
        let map = MarkerMap::create(genmap, &targ_markers);

        // carriers and hi-freq indices
        let max_carriers = std::cmp::max(
            3,
            ((targ_samples + n_ref_samples) as f64 * par.rare() as f64).floor() as usize,
        );
        let carriers = carriers(&targ, &targ_markers, &restrict_ref, targ_samples, max_carriers);
        let mut hi_freq_ind: Vec<usize> = (0..n_targ_markers)
            .filter(|&m| {
                carriers[m]
                    .iter()
                    .filter(|c| matches!(c, Carriers::HighFreq))
                    .count()
                    > 1
            })
            .collect();
        if hi_freq_ind.len() < 2
            || hi_freq_ind.len() as f32 > MAX_HIFREQ_PROP * n_targ_markers as f32
        {
            hi_freq_ind = (0..n_targ_markers).collect();
        }
        let all_markers_are_stage1 = hi_freq_ind.len() == n_targ_markers;
        let carriers = if all_markers_are_stage1 {
            // ignoreLowFreqCarriers: every allele treated as high-frequency
            carriers
                .iter()
                .map(|row| vec![Carriers::HighFreq; row.len()])
                .collect()
        } else {
            carriers
        };

        // stage-1 restrictions
        let stage1_map = if all_markers_are_stage1 {
            MarkerMap::from_gen_pos(map.gen_pos.clone())
        } else {
            map.restrict(&hi_freq_ind)
        };
        let ibs_step = par.step_scale() * median_diff(&stage1_map.gen_pos);
        let stage1_steps = Steps::new(&stage1_map, ibs_step);
        let stage1_targ: Vec<Vec<i16>> =
            hi_freq_ind.iter().map(|&m| targ[m].clone()).collect();
        let stage1_markers: Vec<Arc<Marker>> =
            hi_freq_ind.iter().map(|&m| targ_markers[m].clone()).collect();
        let stage1_ref: Vec<Arc<RefRec>> =
            hi_freq_ind.iter().map(|&m| restrict_ref[m].clone()).collect();
        let stage1_layout = Arc::new(HapLayout::new(
            stage1_markers.iter().map(|m| m.n_alleles).collect(),
        ));
        let stage1_positions: Vec<i32> = stage1_markers.iter().map(|m| m.pos).collect();
        let stage1_overlap = match phased_overlap {
            None => 0,
            Some(ov) => hi_freq_ind.partition_point(|&m| m < ov.alleles.len()),
        };
        let stage1_xref = if n_ref_haps > 0 {
            Some(Arc::new(crate::xref::XRefGT::from_ref_recs(
                &stage1_ref,
                stage1_layout.clone(),
            )))
        } else {
            None
        };

        let prev_stage1_marker = prev_stage1_marker(n_targ_markers, &hi_freq_ind);
        let prev_stage1_wt = prev_wt(&map, &hi_freq_ind);

        // stage1 minor allele frequencies
        let stage1_maf = maf(&stage1_targ, &stage1_markers, &stage1_ref, n_targ_haps, n_ref_haps, 10_000, par.seed);
        let stage1_ibs2 = crate::ibs2::Ibs2::new(&stage1_targ, &stage1_map, &stage1_maf);

        let mut sum_alleles = Vec::with_capacity(n_targ_markers + 1);
        let mut acc = 0u32;
        sum_alleles.push(0);
        for m in &targ_markers {
            acc += m.n_alleles as u32;
            sum_alleles.push(acc);
        }

        FixedPhaseData {
            window_index: window.window_index,
            n_targ_samples: targ_samples,
            n_targ_haps,
            n_haps,
            overlap,
            targ,
            targ_markers,
            carriers,
            map,
            sum_alleles,
            stage1_to2: hi_freq_ind,
            stage1_targ,
            stage1_markers,
            stage1_ref,
            stage1_layout,
            stage1_map,
            ibs_step,
            stage1_steps,
            stage1_overlap,
            stage1_positions,
            prev_stage1_marker,
            prev_stage1_wt,
            n_ref_haps,
            restrict_ref,
            stage1_ibs2,
            stage1_xref,
        }
    }

    pub fn n_stage1_markers(&self) -> usize {
        self.stage1_to2.len()
    }

    pub fn is_low_freq(&self, marker: usize, allele: usize) -> bool {
        self.carriers[marker][allele].is_low_freq()
    }
}

/// A phased-overlap block from the previous window: marker-major alleles
/// for the target haplotypes.
pub struct PhasedOverlap {
    pub alleles: Vec<Vec<i16>>,
}

/// Port of `vcf.Window.carriers(maxCarriers)`.
fn carriers(
    targ: &[Vec<i16>],
    targ_markers: &[Arc<Marker>],
    restrict_ref: &[Arc<RefRec>],
    n_targ_samples: usize,
    max_carriers: usize,
) -> Vec<Vec<Carriers>> {
    use rayon::prelude::*;
    let mut ref_alleles_buf: Vec<u8> = Vec::new();
    let _ = &mut ref_alleles_buf;
    (0..targ.len())
        .into_par_iter()
        .map_init(Vec::new, |ref_buf: &mut Vec<u8>, m| {
            let n_alleles = targ_markers[m].n_alleles as usize;
            let mut lists: Vec<Vec<u32>> = vec![Vec::new(); n_alleles];
            let row = &targ[m];
            for s in 0..n_targ_samples {
                let h1 = s << 1;
                let a1 = row[h1];
                let a2 = row[h1 | 1];
                if a1 >= 0 && lists[a1 as usize].len() <= max_carriers {
                    lists[a1 as usize].push(s as u32);
                }
                if a2 >= 0 && a2 != a1 && lists[a2 as usize].len() <= max_carriers {
                    lists[a2 as usize].push(s as u32);
                }
            }
            if !restrict_ref.is_empty() {
                let rec = &restrict_ref[m];
                ref_buf.clear();
                ref_buf.resize(rec.n_haps, 0);
                rec.fill_alleles(ref_buf);
                let n_ref_samples = rec.n_haps >> 1;
                for s in 0..n_ref_samples {
                    let h1 = s << 1;
                    let a1 = ref_buf[h1] as usize;
                    let a2 = ref_buf[h1 | 1] as usize;
                    if lists[a1].len() <= max_carriers {
                        lists[a1].push((n_targ_samples + s) as u32);
                    }
                    if a2 != a1 && lists[a2].len() <= max_carriers {
                        lists[a2].push((n_targ_samples + s) as u32);
                    }
                }
            }
            lists
                .into_iter()
                .map(|l| {
                    if l.is_empty() {
                        Carriers::ZeroFreq
                    } else if l.len() <= max_carriers {
                        Carriers::List(Arc::new(l))
                    } else {
                        Carriers::HighFreq
                    }
                })
                .collect()
        })
        .collect()
}

/// `FixedPhaseData.medianDiff`
fn median_diff(gen_pos: &[f64]) -> f32 {
    let mut diffs: Vec<f64> = (1..gen_pos.len())
        .map(|j| gen_pos[j] - gen_pos[j - 1])
        .collect();
    diffs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = diffs.len();
    0.5f32 * ((diffs[(n - 1) >> 1] + diffs[n >> 1]) as f32)
}

fn prev_stage1_marker(n_markers: usize, stage1_indices: &[usize]) -> Vec<usize> {
    let mut mkr = vec![0usize; n_markers];
    let n_hi_freq = stage1_indices.len();
    let mut start = stage1_indices[1];
    for j in 2..n_hi_freq {
        let end = stage1_indices[j];
        for v in mkr[start..end].iter_mut() {
            *v = j - 1;
        }
        start = end;
    }
    for v in mkr[start..n_markers].iter_mut() {
        *v = n_hi_freq - 1;
    }
    mkr
}

fn prev_wt(map: &MarkerMap, marker_indices: &[usize]) -> Vec<f32> {
    let gen_pos = &map.gen_pos;
    let mut prev_wt = vec![0.0f32; gen_pos.len()];
    for v in prev_wt[0..marker_indices[0]].iter_mut() {
        *v = 1.0;
    }
    let mut start = marker_indices[0];
    for &end in &marker_indices[1..] {
        let pos_a = gen_pos[start];
        let pos_b = gen_pos[end];
        let d = pos_b - pos_a;
        prev_wt[start] = 1.0;
        for m in start + 1..end {
            prev_wt[m] = ((pos_b - gen_pos[m]) / d) as f32;
        }
        start = end;
    }
    for v in prev_wt[start..].iter_mut() {
        *v = 1.0;
    }
    prev_wt
}

/// `FixedPhaseData.maf` over stage1 markers.
fn maf(
    stage1_targ: &[Vec<i16>],
    stage1_markers: &[Arc<Marker>],
    stage1_ref: &[Arc<RefRec>],
    n_targ_haps: usize,
    n_ref_haps: usize,
    max_haps: usize,
    seed: i64,
) -> Vec<f32> {
    let mut rand = JavaRandom::new(seed);
    let targ_haps = rand_haps(n_targ_haps, max_haps, &mut rand);
    let ref_haps: Vec<i32> = if targ_haps.len() < max_haps && n_ref_haps > 0 {
        rand_haps(n_ref_haps, max_haps - targ_haps.len(), &mut rand)
    } else {
        Vec::new()
    };
    use rayon::prelude::*;
    (0..stage1_targ.len())
        .into_par_iter()
        .map_init(Vec::new, |ref_buf: &mut Vec<u8>, m| {
            let n_alleles = stage1_markers[m].n_alleles as usize;
            let mut mod_cnts = vec![0u32; n_alleles + 1];
            let row = &stage1_targ[m];
            for &h in &targ_haps {
                let a = row[h as usize];
                mod_cnts[(a + 1) as usize] += 1;
            }
            if !ref_haps.is_empty() {
                let rec = &stage1_ref[m];
                ref_buf.clear();
                ref_buf.resize(rec.n_haps, 0);
                rec.fill_alleles(ref_buf);
                for &h in &ref_haps {
                    mod_cnts[ref_buf[h as usize] as usize + 1] += 1;
                }
            }
            mod_cnts[0] = 0; // zero-out missing count
            mod_cnts.sort_unstable();
            let den: u32 = mod_cnts[1..].iter().sum();
            if den == 0 {
                0.0f32
            } else {
                (mod_cnts[mod_cnts.len() - 2] as f64 / den as f64) as f32
            }
        })
        .collect()
}

fn rand_haps(n_haps: usize, max_haps: usize, rand: &mut JavaRandom) -> Vec<i32> {
    let mut ia: Vec<i32> = (0..n_haps as i32).collect();
    if n_haps > max_haps {
        // Utilities.shuffle(ia, maxHaps, rand)
        for j in 0..max_haps {
            let x = rand.next_int_bound((ia.len() - j) as i32) as usize;
            ia.swap(j, j + x);
        }
        ia.truncate(max_haps);
        ia.sort_unstable();
    }
    ia
}

// ---------------------------------------------------------------------------
// EstPhase / PhaseData (ports of phase.EstPhase, phase.PhaseData)

#[allow(dead_code)]
pub struct EstPhase {
    pub fpd: Arc<FixedPhaseData>,
    pub phase: Vec<std::sync::Mutex<Option<SamplePhase>>>,
}

impl EstPhase {
    #[allow(dead_code)]
    pub fn get_clone_haps(&self, sample: usize) -> (BitArray, BitArray) {
        let guard = self.phase[sample].lock().unwrap();
        let sp = guard.as_ref().unwrap();
        (sp.hap1.clone(), sp.hap2.clone())
    }

    /// takes the SamplePhase out for exclusive processing
    pub fn take(&self, sample: usize) -> SamplePhase {
        self.phase[sample].lock().unwrap().take().unwrap()
    }

    pub fn put(&self, sample: usize, sp: SamplePhase) {
        *self.phase[sample].lock().unwrap() = Some(sp);
    }

    pub fn with<R>(&self, sample: usize, f: impl FnOnce(&SamplePhase) -> R) -> R {
        let guard = self.phase[sample].lock().unwrap();
        f(guard.as_ref().unwrap())
    }
}

pub struct PhaseData {
    pub est_phase: Arc<EstPhase>,
    pub seed: i64,
    pub it: usize,
    pub lr_threshold: f32,
    pub recomb_intensity: f32,
    pub p_recomb: Arc<Vec<f32>>,
    pub p_mismatch: f32,
}

impl PhaseData {
    pub fn new(fpd: Arc<FixedPhaseData>, par: &Par, seed: i64) -> PhaseData {
        let est_phase = Arc::new(crate::initphase::init_phase(&fpd, par, seed));
        let recomb_intensity = 0.04f32 * par.ne / fpd.n_haps as f32;
        let p_recomb = Arc::new(fpd.stage1_map.p_recomb(recomb_intensity));
        let p_mismatch = crate::par::li_stephens_p_mismatch(fpd.n_haps);
        PhaseData {
            est_phase,
            seed,
            it: 0,
            lr_threshold: lr_threshold(par, 0),
            recomb_intensity,
            p_recomb,
            p_mismatch,
        }
    }

    pub fn fpd(&self) -> &Arc<FixedPhaseData> {
        &self.est_phase.fpd
    }

    pub fn update_recomb_intensity(&mut self, recomb_intensity: f32) {
        assert!(recomb_intensity > 0.0 && recomb_intensity.is_finite());
        self.recomb_intensity = recomb_intensity;
        self.p_recomb = Arc::new(self.est_phase.fpd.stage1_map.p_recomb(recomb_intensity));
    }

    pub fn update_p_mismatch(&mut self, p_mismatch: f32) {
        assert!((0.0..=1.0).contains(&p_mismatch) && p_mismatch.is_finite());
        self.p_mismatch = p_mismatch;
    }

    pub fn increment_it(&mut self, par: &Par) {
        self.it += 1;
        self.lr_threshold = lr_threshold(par, self.it);
    }

    pub fn advance_to_first_phasing_it(&mut self, par: &Par) {
        let n_burnin = par.burnin as usize;
        if self.it < n_burnin {
            self.it = n_burnin;
            self.lr_threshold = lr_threshold(par, self.it);
        }
    }

    /// `PhaseData.seed()`: iteration-dependent seed
    pub fn it_seed(&self) -> i64 {
        self.seed.wrapping_add(self.it as i64)
    }
}

/// `PhaseData.lrThreshold(par, it)`
pub fn lr_threshold(par: &Par, it: usize) -> f32 {
    let n_burnin_its = par.burnin as usize;
    let n_its_m1 = par.iterations as usize - 1;
    if it < n_burnin_its {
        f32::INFINITY
    } else if it == n_its_m1 + n_burnin_its {
        1.0
    } else {
        let last_val = 4.0f64;
        let exp = (n_its_m1 - (it - n_burnin_its)) as f64 / n_its_m1 as f64;
        let base = par.initial_lr() as f64 / last_val;
        (last_val * base.powf(exp)) as f32
    }
}

// ---------------------------------------------------------------------------
// MarkerCluster (port of phase.MarkerCluster); operates on a SamplePhase

pub struct MarkerClusterInfo {
    pub cluster_to_end: Vec<usize>,
    pub unph_het_clusters: Vec<usize>,
    pub p_recomb: Vec<f32>,
}

impl MarkerClusterInfo {
    pub fn new(sp: &SamplePhase, marker_p_recomb: &[f32]) -> MarkerClusterInfo {
        let cluster_to_end = sp.clust_ends();
        let unph_het_clusters: Vec<usize> = (0..sp.n_clusters())
            .filter(|&c| sp.clust_type_at(c) == CLUST_UNPHASED_HET)
            .collect();
        let n_clusters = cluster_to_end.len();
        let mut p_clust_recomb = vec![0.0f32; n_clusters];
        let mut start = cluster_to_end[0];
        for j in 1..n_clusters {
            let end = cluster_to_end[j];
            let mut p_no_recomb = 1.0f32;
            for k in start..end {
                p_no_recomb *= 1.0f32 - marker_p_recomb[k];
            }
            p_clust_recomb[j] = 1.0f32 - p_no_recomb;
            start = end;
        }
        MarkerClusterInfo {
            cluster_to_end,
            unph_het_clusters,
            p_recomb: p_clust_recomb,
        }
    }

    #[inline]
    pub fn n_clusters(&self) -> usize {
        self.cluster_to_end.len()
    }

    #[inline]
    pub fn cluster_start(&self, index: usize) -> usize {
        if index == 0 {
            0
        } else {
            self.cluster_to_end[index - 1]
        }
    }

    #[inline]
    pub fn cluster_end(&self, index: usize) -> usize {
        self.cluster_to_end[index]
    }
}
