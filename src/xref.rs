//! Port of `vcf.XRefGT`: phased genotypes in haplotype-major bit-packed form.
//!
//! Unlike the Java class, the combined (target + reference) view keeps the
//! reference haplotypes behind an `Arc` so each phasing iteration reuses the
//! same reference bit arrays instead of cloning them.

use crate::bits::BitArray;
use crate::phasedata::HapLayout;
use crate::refpanel::RefRec;
use rayon::prelude::*;
use std::sync::Arc;

pub struct XRefGT {
    pub layout: Arc<HapLayout>,
    own: Vec<BitArray>,
    shared: Option<Arc<XRefGT>>, // haplotypes appended after `own`
}

impl XRefGT {
    /// Builds hap-major bit arrays for the given (restricted) reference
    /// records. Equivalent to `XRefGT.fromPhasedGT`.
    pub fn from_ref_recs(recs: &[Arc<RefRec>], layout: Arc<HapLayout>) -> XRefGT {
        let n_haps = if recs.is_empty() { 0 } else { recs[0].n_haps };
        let n_bits = layout.total_bits();
        // build per-marker allele buffers once, then scatter into haps in
        // parallel batches of haplotypes
        let allele_rows: Vec<Vec<u8>> = recs
            .par_iter()
            .map(|rec| {
                let mut buf = vec![0u8; n_haps];
                rec.fill_alleles(&mut buf);
                buf
            })
            .collect();
        let haps: Vec<BitArray> = (0..n_haps)
            .into_par_iter()
            .map(|h| {
                let mut bits = BitArray::new(n_bits);
                for (m, row) in allele_rows.iter().enumerate() {
                    let allele = row[h] as i32;
                    if allele != 0 {
                        set_allele_bits(&layout, m, allele, &mut bits);
                    }
                }
                bits
            })
            .collect();
        XRefGT {
            layout,
            own: haps,
            shared: None,
        }
    }

    /// Builds an XRefGT for the target haplotypes from per-sample phase data.
    /// Equivalent to `XRefGT.from(samples, phase)`.
    pub fn from_est_phase(est: &crate::phasedata::EstPhase) -> XRefGT {
        let layout = est.fpd.stage1_layout.clone();
        let n_samples = est.phase.len();
        let haps: Vec<BitArray> = (0..n_samples << 1)
            .into_par_iter()
            .map(|h| {
                let guard = est.phase[h >> 1].lock().unwrap();
                let sp = guard.as_ref().unwrap();
                if h & 1 == 0 {
                    sp.hap1.clone()
                } else {
                    sp.hap2.clone()
                }
            })
            .collect();
        XRefGT {
            layout,
            own: haps,
            shared: None,
        }
    }

    /// `XRefGT.combine(first, second)`: haplotypes of `first` then `second`.
    /// The second (reference) haplotypes are shared, not copied.
    pub fn combine(first: XRefGT, second: &Arc<XRefGT>) -> XRefGT {
        debug_assert!(first.shared.is_none() && second.shared.is_none());
        XRefGT {
            layout: first.layout,
            own: first.own,
            shared: Some(second.clone()),
        }
    }

    #[inline]
    fn hap(&self, hap: usize) -> &BitArray {
        if hap < self.own.len() {
            &self.own[hap]
        } else {
            &self.shared.as_ref().unwrap().own[hap - self.own.len()]
        }
    }

    #[inline]
    pub fn n_haps(&self) -> usize {
        self.own.len() + self.shared.as_ref().map_or(0, |s| s.own.len())
    }

    #[inline]
    pub fn n_markers(&self) -> usize {
        self.layout.n_markers()
    }

    #[inline]
    pub fn allele(&self, m: usize, hap: usize) -> i32 {
        self.layout.allele(self.hap(hap), m)
    }

    /// `XRefGT.hash(hap, start, end)`
    #[inline]
    pub fn hash(&self, hap: usize, start: usize, end: usize) -> i32 {
        let start_bit = self.layout.sum_hap_bits[start] as usize;
        let end_bit = self.layout.sum_hap_bits[end] as usize;
        self.hap(hap).hash(start_bit, end_bit)
    }

    /// `XRefGT.copyTo(hap, start, end, bitList)`
    #[inline]
    pub fn copy_to(&self, hap: usize, start: usize, end: usize, dst: &mut BitArray) {
        let start_bit = self.layout.sum_hap_bits[start] as usize;
        let end_bit = self.layout.sum_hap_bits[end] as usize;
        dst.copy_from(self.hap(hap), start_bit, end_bit);
    }
}

#[inline]
fn set_allele_bits(layout: &HapLayout, m: usize, allele: i32, bits: &mut BitArray) {
    let start = layout.sum_hap_bits[m] as usize;
    let end = layout.sum_hap_bits[m + 1] as usize;
    let mut mask = 1i32;
    for j in start..end {
        if allele & mask == mask {
            bits.set(j);
        }
        mask <<= 1;
    }
}
