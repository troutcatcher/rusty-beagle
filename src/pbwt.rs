// SPDX-License-Identifier: GPL-3.0-or-later
//
// rusty-beagle - a Rust port of Beagle 5.5 genotype phasing and imputation.
// Copyright (C) 2026 The rusty-beagle authors
//
// This file is part of a Rust port of Beagle 5.5 (release
// beagle.27Feb25.75f), Copyright (C) 2014-2024 Brian L. Browning, and is
// distributed as a modified version of that GPL-licensed work.  The module
// documentation below names the upstream Java class(es) this file
// corresponds to; docs/PORT_NOTES.md records the full source-to-source
// mapping and the places where this port deviates from the Java.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

//! Ports of `beagleutil.PbwtUpdater` and `beagleutil.PbwtDivUpdater`.

/// `beagleutil.PbwtUpdater`: prefix-array update (no divergence).
pub struct PbwtUpdater {
    n_haps: usize,
    a: Vec<Vec<i32>>,
}

impl PbwtUpdater {
    pub fn new(n_haps: usize) -> PbwtUpdater {
        PbwtUpdater {
            n_haps,
            a: (0..4).map(|_| Vec::new()).collect(),
        }
    }

    /// `update(alleles, nAlleles, prefix)`
    pub fn update(&mut self, alleles: &[i32], n_alleles: usize, prefix: &mut [i32]) {
        debug_assert_eq!(alleles.len(), self.n_haps);
        debug_assert_eq!(prefix.len(), self.n_haps);
        if n_alleles > self.a.len() {
            self.a.resize_with(n_alleles, Vec::new);
        }
        for l in self.a[..n_alleles].iter_mut() {
            l.clear();
        }
        for &h in prefix.iter() {
            let allele = alleles[h as usize] as usize;
            self.a[allele].push(h);
        }
        let mut start = 0;
        for al in 0..n_alleles {
            let size = self.a[al].len();
            prefix[start..start + size].copy_from_slice(&self.a[al]);
            start += size;
        }
        debug_assert_eq!(start, self.n_haps);
    }
}

/// `beagleutil.PbwtDivUpdater`: prefix + divergence array updates.
pub struct PbwtDivUpdater {
    n_haps: usize,
    a: Vec<Vec<i32>>,
    d: Vec<Vec<i32>>,
    p: Vec<i32>,
}

impl PbwtDivUpdater {
    pub fn new(n_haps: usize) -> PbwtDivUpdater {
        PbwtDivUpdater {
            n_haps,
            a: (0..4).map(|_| Vec::new()).collect(),
            d: (0..4).map(|_| Vec::new()).collect(),
            p: vec![0; 4],
        }
    }

    fn ensure_capacity(&mut self, n_alleles: usize) {
        if n_alleles > self.a.len() {
            self.a.resize_with(n_alleles, Vec::new);
            self.d.resize_with(n_alleles, Vec::new);
            self.p.resize(n_alleles, 0);
        }
    }

    /// `fwdUpdate(rec, nAlleles, marker, prefix, div)`;
    /// `get(h)` returns the allele of haplotype h.
    pub fn fwd_update<F: Fn(usize) -> u32>(
        &mut self,
        get: F,
        n_alleles: usize,
        marker: i32,
        prefix: &mut [i32],
        div: &mut [i32],
    ) {
        debug_assert_eq!(prefix.len(), self.n_haps);
        self.ensure_capacity(n_alleles);
        for j in 0..n_alleles {
            self.a[j].clear();
            self.d[j].clear();
            self.p[j] = marker + 1;
        }
        for i in 0..self.n_haps {
            let h = prefix[i];
            let allele = get(h as usize) as usize;
            debug_assert!(allele < n_alleles);
            let di = div[i];
            for j in 0..n_alleles {
                if di > self.p[j] {
                    self.p[j] = di;
                }
            }
            self.a[allele].push(h);
            self.d[allele].push(self.p[allele]);
            self.p[allele] = i32::MIN;
        }
        self.update_prefix_and_div(n_alleles, prefix, div);
    }

    /// `bwdUpdate(rec, nAlleles, marker, prefix, div)`
    pub fn bwd_update<F: Fn(usize) -> u32>(
        &mut self,
        get: F,
        n_alleles: usize,
        marker: i32,
        prefix: &mut [i32],
        div: &mut [i32],
    ) {
        debug_assert_eq!(prefix.len(), self.n_haps);
        self.ensure_capacity(n_alleles);
        for j in 0..n_alleles {
            self.a[j].clear();
            self.d[j].clear();
            self.p[j] = marker - 1;
        }
        for i in 0..self.n_haps {
            let h = prefix[i];
            let allele = get(h as usize) as usize;
            debug_assert!(allele < n_alleles);
            let di = div[i];
            for j in 0..n_alleles {
                if di < self.p[j] {
                    self.p[j] = di;
                }
            }
            self.a[allele].push(h);
            self.d[allele].push(self.p[allele]);
            self.p[allele] = i32::MAX;
        }
        self.update_prefix_and_div(n_alleles, prefix, div);
    }

    fn update_prefix_and_div(&mut self, n_alleles: usize, prefix: &mut [i32], div: &mut [i32]) {
        let mut start = 0usize;
        for al in 0..n_alleles {
            let size = self.a[al].len();
            prefix[start..start + size].copy_from_slice(&self.a[al]);
            div[start..start + size].copy_from_slice(&self.d[al]);
            start += size;
        }
        debug_assert_eq!(start, self.n_haps);
    }
}
