//! Port of `imp.ImpData` and `imp.HaplotypeCoder`.

use crate::genmap::GeneticMap;
use crate::par::Par;
use crate::refpanel::{RefAlleles, RefRec};
use crate::windows::Window;
use rayon::prelude::*;
use std::sync::Arc;

const MIN_CM_DIST: f64 = 1e-7;

/// Per-cluster map from haplotype index (ref haps first, then target haps)
/// to allele-sequence index; `value_size` is the number of distinct
/// sequence indices (index 0 = ref-only sequences not present in the target).
pub enum ClusterCoding {
    /// all cluster markers share one sequence-coded block:
    /// ref hap h -> seq1_to_seq2[hap2seq[h]]
    Composed {
        block: u32,
        hap2seq: Arc<Vec<u16>>,
        seq1_to_seq2: Vec<u32>,
        targ: Vec<u32>,
        value_size: u32,
        /// `seq1_to_seq2` transposed into one bitset per allele-sequence
        /// index, keyed by block sequence
        match_bits: MatchBits,
    },
    /// explicit per-ref-hap codes
    Direct {
        coded_ref: Vec<u32>,
        targ: Vec<u32>,
        value_size: u32,
        /// `coded_ref` transposed into one bitset per allele-sequence index,
        /// keyed by reference haplotype
        match_bits: MatchBits,
    },
}

/// One bitset per allele-sequence index, marking which entries of a code
/// array carry that index.
///
/// The HMM's per-cluster inner loop only asks "does this state carry the
/// target's allele sequence?", so it can read a single bit instead of
/// gathering a `u32` code.  Clusters have very few distinct sequences (~3 in
/// practice), which makes one row far smaller than the code array — for
/// `Direct` clusters, 1.25 KB per row against a 40 KB `coded_ref` at 10,000
/// reference haplotypes — and so much likelier to stay in L1.
pub struct MatchBits {
    words: Vec<u64>,
    stride: usize,
    n_codes: usize,
}

impl MatchBits {
    fn build(codes: &[u32], n_codes: usize) -> MatchBits {
        let stride = (codes.len() + 63) >> 6;
        let mut words = vec![0u64; n_codes * stride];
        for (i, &c) in codes.iter().enumerate() {
            debug_assert!((c as usize) < n_codes);
            words[(c as usize) * stride + (i >> 6)] |= 1u64 << (i & 63);
        }
        MatchBits {
            words,
            stride,
            n_codes,
        }
    }

    /// The bitset for `code`, or `None` when no entry carries it.
    #[inline]
    pub fn row(&self, code: u32) -> Option<&[u64]> {
        let c = code as usize;
        if c < self.n_codes {
            Some(&self.words[c * self.stride..(c + 1) * self.stride])
        } else {
            None
        }
    }
}

/// Tests bit `i` of a `MatchBits` row.
#[inline]
pub fn match_bit(row: Option<&[u64]>, i: usize) -> bool {
    match row {
        Some(r) => (r[i >> 6] >> (i & 63)) & 1 != 0,
        None => false,
    }
}

impl ClusterCoding {
    #[inline]
    pub fn get(&self, hap: usize, n_ref_haps: usize) -> u32 {
        match self {
            ClusterCoding::Composed {
                hap2seq,
                seq1_to_seq2,
                targ,
                ..
            } => {
                if hap < n_ref_haps {
                    seq1_to_seq2[hap2seq[hap] as usize]
                } else {
                    targ[hap - n_ref_haps]
                }
            }
            ClusterCoding::Direct {
                coded_ref, targ, ..
            } => {
                if hap < n_ref_haps {
                    coded_ref[hap]
                } else {
                    targ[hap - n_ref_haps]
                }
            }
        }
    }

    #[inline]
    pub fn value_size(&self) -> u32 {
        match self {
            ClusterCoding::Composed { value_size, .. } => *value_size,
            ClusterCoding::Direct { value_size, .. } => *value_size,
        }
    }
}

/// Port of `imp.ImpData`.
#[allow(dead_code)] // parity fields kept for the phasing port
pub struct ImpData {
    pub imp_states: usize,
    pub imp_step: f32,
    pub imp_segment: f32,
    pub imp_nsteps: usize,
    pub seed: i64,
    pub ap: bool,
    pub gp: bool,
    pub n_threads: usize,

    pub n_clusters: usize,
    pub n_ref_haps: usize,
    pub n_targ_haps: usize,
    pub n_input_targ_haps: usize,
    pub n_haps: usize,

    pub targ_clust_start_end: Vec<usize>,
    pub ref_cluster_start: Vec<usize>,
    pub ref_cluster_end: Vec<usize>,
    pub hap_to_seq: Vec<ClusterCoding>,
    pub err_prob: Vec<f32>,
    pub pos: Vec<f64>,
    pub p_recomb: Vec<f32>,
    pub weight: Vec<f32>,

    /// target alleles, marker-major: `targ_alleles[t][hap]`
    pub targ_alleles: Vec<Vec<u8>>,
    pub targ_is_diploid: Vec<bool>,
}

impl ImpData {
    pub fn new(
        par: &Par,
        window: &Window,
        genmap: &GeneticMap,
        targ_samples: &crate::vcfio::Samples,
        targ_alleles: Vec<Vec<u8>>,
    ) -> ImpData {
        let indices = &window.indices;
        let targ_to_ref = &indices.targ_marker_to_marker;
        let n_targ_markers = targ_to_ref.len();
        let ref_recs = &window.ref_recs;
        let n_targ_haps = if n_targ_markers > 0 {
            targ_alleles[0].len()
        } else {
            0
        };
        let n_ref_haps = ref_recs[0].n_haps;
        let n_haps = n_ref_haps + n_targ_haps;
        let n_input_targ_haps = targ_samples
            .is_diploid
            .iter()
            .map(|&d| if d { 2 } else { 1 })
            .sum();

        // cumPos over target markers
        let targ_pos = cum_pos_targ(window, genmap);
        let block_end = targ_block_end(ref_recs, targ_to_ref);
        let targ_clust_start_end =
            targ_clust_start_end(&targ_pos, &block_end, par.cluster);
        let pos = mid_pos(&targ_pos, &targ_clust_start_end);
        let n_clusters = targ_clust_start_end.len() - 1;

        let restrict_ref = window.restrict_ref();
        let hap_to_seq: Vec<ClusterCoding> = (0..n_clusters)
            .into_par_iter()
            .map(|c| {
                code_cluster(
                    &restrict_ref,
                    &targ_alleles,
                    targ_clust_start_end[c],
                    targ_clust_start_end[c + 1],
                    n_ref_haps,
                )
            })
            .collect();

        let ref_cluster_start: Vec<usize> = (0..n_clusters)
            .map(|j| targ_to_ref[targ_clust_start_end[j]])
            .collect();
        let ref_cluster_end: Vec<usize> = (1..=n_clusters)
            .map(|j| targ_to_ref[targ_clust_start_end[j] - 1] + 1)
            .collect();

        let err_rate = par.err_for(n_haps);
        let err_prob = err(err_rate, &targ_clust_start_end);
        let p_recomb = p_recomb(par.ne, n_ref_haps, &pos);
        let weight = wts(ref_recs, &ref_cluster_start, &ref_cluster_end, genmap);

        ImpData {
            imp_states: par.imp_states,
            imp_step: par.imp_step,
            imp_segment: par.imp_segment,
            imp_nsteps: par.imp_nsteps,
            seed: par.seed,
            ap: par.ap,
            gp: par.gp,
            n_threads: par.nthreads,
            n_clusters,
            n_ref_haps,
            n_targ_haps,
            n_input_targ_haps,
            n_haps,
            targ_clust_start_end,
            ref_cluster_start,
            ref_cluster_end,
            hap_to_seq,
            err_prob,
            pos,
            p_recomb,
            weight,
            targ_alleles,
            targ_is_diploid: targ_samples.is_diploid.as_ref().clone(),
        }
    }

    /// `ImpData.allele(cluster, hap)`
    #[inline]
    pub fn allele(&self, cluster: usize, hap: usize) -> u32 {
        self.hap_to_seq[cluster].get(hap, self.n_ref_haps)
    }

    #[allow(dead_code)]
    pub fn targ_cluster_start(&self, cluster: usize) -> usize {
        self.targ_clust_start_end[cluster]
    }

    #[allow(dead_code)]
    pub fn targ_cluster_end(&self, cluster: usize) -> usize {
        self.targ_clust_start_end[cluster + 1]
    }
}

/// `ImpData.cumPos` over the target markers of the window.
fn cum_pos_targ(window: &Window, map: &GeneticMap) -> Vec<f64> {
    let markers: Vec<&crate::marker::Marker> =
        window.targ_recs.iter().map(|r| &r.marker).collect();
    let mut cum = vec![0.0f64; markers.len()];
    let mut last = map.gen_pos_marker(markers[0]);
    for j in 1..cum.len() {
        let gen_pos = map.gen_pos_marker(markers[j]);
        let dist = (gen_pos - last).abs().max(MIN_CM_DIST);
        cum[j] = cum[j - 1] + dist;
        last = gen_pos;
    }
    cum
}

/// `ImpData.cumPos` over reference markers.
fn cum_pos_ref(ref_recs: &[Arc<RefRec>], map: &GeneticMap) -> Vec<f64> {
    let mut cum = vec![0.0f64; ref_recs.len()];
    let mut last = map.gen_pos_marker(&ref_recs[0].marker);
    for j in 1..cum.len() {
        let gen_pos = map.gen_pos_marker(&ref_recs[j].marker);
        let dist = (gen_pos - last).abs().max(MIN_CM_DIST);
        cum[j] = cum[j - 1] + dist;
        last = gen_pos;
    }
    cum
}

/// `ImpData.targBlockEnd`: cluster-split boundaries where the sequence-coded
/// block (Java `map(0)` identity) changes.
fn targ_block_end(ref_recs: &[Arc<RefRec>], targ_to_ref: &[usize]) -> Vec<usize> {
    let mut list = Vec::new();
    let mut last_block: Option<u32> = None;
    for (j, &ref_index) in targ_to_ref.iter().enumerate() {
        if let Some(bid) = ref_recs[ref_index].block_id() {
            if last_block != Some(bid) {
                if last_block.is_some() {
                    list.push(j);
                }
                last_block = Some(bid);
            }
        }
    }
    list.push(targ_to_ref.len());
    list
}

/// `ImpData.targClustStartEnd`
fn targ_clust_start_end(raw_pos: &[f64], targ_block_end: &[usize], cluster_dist: f32) -> Vec<usize> {
    let cluster_dist = cluster_dist as f64;
    let mut clust_start_end = vec![0usize; raw_pos.len() + 1];
    let mut size = 1usize; // clust_start_end[0] = 0
    for &block_end in targ_block_end {
        let clust_start = clust_start_end[size - 1];
        let mut start_pos = raw_pos[clust_start];
        for m in clust_start + 1..block_end {
            let pos = raw_pos[m];
            if pos - start_pos > cluster_dist {
                clust_start_end[size] = m;
                size += 1;
                start_pos = pos;
            }
        }
        clust_start_end[size] = block_end;
        size += 1;
    }
    clust_start_end.truncate(size);
    clust_start_end
}

/// `ImpData.midPos`
fn mid_pos(pos: &[f64], start_end: &[usize]) -> Vec<f64> {
    (1..start_end.len())
        .map(|j| (pos[start_end[j - 1]] + pos[start_end[j] - 1]) / 2.0)
        .collect()
}

/// `ImpData.err`
fn err(err_rate: f32, start_end: &[usize]) -> Vec<f32> {
    let max_err_prob = 0.5f32;
    (0..start_end.len() - 1)
        .map(|j| {
            let e = err_rate * (start_end[j + 1] - start_end[j]) as f32;
            if e > max_err_prob {
                max_err_prob
            } else {
                e
            }
        })
        .collect()
}

/// `ImpData.pRecomb`: note `nHaps` here is the number of *reference*
/// haplotypes (Java passes `refGT.nHaps()`).
fn p_recomb(ne: f32, n_ref_haps: usize, pos: &[f64]) -> Vec<f32> {
    let mut p = vec![0.0f32; pos.len()];
    let c = -((0.04f64 * ne as f64) / n_ref_haps as f64); // 0.04 = 4/(100 cM/M)
    for j in 1..p.len() {
        p[j] = (-f64::exp_m1(c * (pos[j] - pos[j - 1]))) as f32;
    }
    p
}

/// `ImpData.wts`
fn wts(
    ref_recs: &[Arc<RefRec>],
    ref_cluster_start: &[usize],
    ref_cluster_end: &[usize],
    map: &GeneticMap,
) -> Vec<f32> {
    let cum_pos = cum_pos_ref(ref_recs, map);
    let n_targ_markers_m1 = ref_cluster_start.len() - 1;
    let mut wts = vec![f32::NAN; cum_pos.len()];
    for j in 0..n_targ_markers_m1 {
        let end = ref_cluster_end[j];
        let next_start = ref_cluster_start[j + 1];
        let next_start_pos = cum_pos[next_start];
        let total_length = next_start_pos - cum_pos[end - 1];
        for m in end..next_start {
            wts[m] = ((cum_pos[next_start] - cum_pos[m]) / total_length) as f32;
        }
    }
    wts
}

/// Port of `imp.HaplotypeCoder.run(start, end)` for one marker cluster.
fn code_cluster(
    restrict_ref: &[Arc<RefRec>],
    targ_alleles: &[Vec<u8>],
    start: usize,
    end: usize,
    n_ref_haps: usize,
) -> ClusterCoding {
    debug_assert!(start < end);
    // is_hap_coded: every record sequence-coded with the same block
    let hap_coded_block = {
        match restrict_ref[start].block_id() {
            None => None,
            Some(bid) => {
                if restrict_ref[start + 1..end]
                    .iter()
                    .all(|r| r.block_id() == Some(bid))
                {
                    Some(bid)
                } else {
                    None
                }
            }
        }
    };

    // codeTarg: index the target allele sequences and build seqMaps
    let n_targ_haps = targ_alleles.get(start).map_or(0, |v| v.len());
    let mut hap_to_seq = vec![1u32; n_targ_haps];
    let mut seq_maps: Vec<Vec<u32>> = Vec::with_capacity(end - start);
    let mut seq_cnt: u32 = 2;
    for m in start..end {
        let n_alleles = restrict_ref[m].marker.n_alleles as u32;
        let mut seq_map = vec![0u32; (seq_cnt * n_alleles) as usize];
        seq_cnt = 1;
        let alleles = &targ_alleles[m];
        for h in 0..n_targ_haps {
            let index = (n_alleles * hap_to_seq[h] + alleles[h] as u32) as usize;
            if seq_map[index] == 0 {
                seq_map[index] = seq_cnt;
                seq_cnt += 1;
            }
            hap_to_seq[h] = seq_map[index];
        }
        seq_maps.push(seq_map);
    }
    let coded_targ = hap_to_seq;
    let value_size = seq_cnt;

    match hap_coded_block {
        Some(block_id) => {
            // codeSeqCodedRef: compose through the block's seq alphabet
            let (block_hap2seq, first_seq2allele) = match &restrict_ref[start].alleles {
                RefAlleles::SeqCoded {
                    hap2seq,
                    seq2allele,
                    ..
                } => (hap2seq.clone(), seq2allele.len()),
                _ => unreachable!(),
            };
            let mut seq1_to_seq2 = vec![1u32; first_seq2allele];
            for (j, m) in (start..end).enumerate() {
                let n_alleles = restrict_ref[m].marker.n_alleles as u32;
                let seq2allele = match &restrict_ref[m].alleles {
                    RefAlleles::SeqCoded { seq2allele, .. } => seq2allele,
                    _ => unreachable!(),
                };
                let seq_map = &seq_maps[j];
                for s in 0..seq1_to_seq2.len() {
                    if seq1_to_seq2[s] > 0 {
                        let index =
                            (seq1_to_seq2[s] * n_alleles + seq2allele[s] as u32) as usize;
                        seq1_to_seq2[s] = seq_map[index];
                    }
                }
            }
            let match_bits = MatchBits::build(&seq1_to_seq2, value_size as usize);
            ClusterCoding::Composed {
                block: block_id,
                hap2seq: block_hap2seq,
                seq1_to_seq2,
                targ: coded_targ,
                value_size,
                match_bits,
            }
        }
        None => {
            // codeSeq: run every ref haplotype through the seqMaps
            let mut coded_ref = vec![1u32; n_ref_haps];
            let mut patches: Vec<(u32, u32)> = Vec::new();
            for (j, m) in (start..end).enumerate() {
                let n_alleles = restrict_ref[m].marker.n_alleles as u32;
                let seq_map = &seq_maps[j];
                match &restrict_ref[m].alleles {
                    RefAlleles::AlleleCoded { major, carriers } => {
                        // compute carrier updates from the old codes first
                        patches.clear();
                        for (a, list) in carriers.iter().enumerate() {
                            if a != *major as usize {
                                for &h in list {
                                    let c = coded_ref[h as usize];
                                    if c > 0 {
                                        let idx = (c * n_alleles + a as u32) as usize;
                                        patches.push((h, seq_map[idx]));
                                    }
                                }
                            }
                        }
                        let major = *major as u32;
                        for c in coded_ref.iter_mut() {
                            if *c > 0 {
                                *c = seq_map[(*c * n_alleles + major) as usize];
                            }
                        }
                        for &(h, v) in &patches {
                            coded_ref[h as usize] = v;
                        }
                    }
                    RefAlleles::SeqCoded {
                        hap2seq,
                        seq2allele,
                        ..
                    } => {
                        for h in 0..n_ref_haps {
                            let c = coded_ref[h];
                            if c > 0 {
                                let allele = seq2allele[hap2seq[h] as usize] as u32;
                                coded_ref[h] = seq_map[(c * n_alleles + allele) as usize];
                            }
                        }
                    }
                }
            }
            let match_bits = MatchBits::build(&coded_ref, value_size as usize);
            ClusterCoding::Direct {
                coded_ref,
                targ: coded_targ,
                value_size,
                match_bits,
            }
        }
    }
}
