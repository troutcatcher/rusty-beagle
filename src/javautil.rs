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

//! Ports of Java runtime utilities whose exact behavior Beagle's output
//! depends on: `java.util.Random`, `java.util.PriorityQueue` ordering,
//! `blbutil.Utilities.shuffle`, and `DecimalFormat` value tables.

/// Port of `java.util.Random`: 48-bit linear congruential generator.
#[derive(Clone)]
pub struct JavaRandom {
    seed: u64,
}

const MULTIPLIER: u64 = 0x5DEECE66D;
const ADDEND: u64 = 0xB;
const MASK: u64 = (1u64 << 48) - 1;

impl JavaRandom {
    pub fn new(seed: i64) -> Self {
        JavaRandom {
            seed: (seed as u64 ^ MULTIPLIER) & MASK,
        }
    }

    #[inline]
    fn next(&mut self, bits: u32) -> i32 {
        self.seed = self.seed.wrapping_mul(MULTIPLIER).wrapping_add(ADDEND) & MASK;
        (self.seed >> (48 - bits)) as i32
    }

    /// `Random.nextInt()`
    #[inline]
    pub fn next_int(&mut self) -> i32 {
        self.next(32)
    }

    /// `Random.setSeed(seed)`
    pub fn set_seed(&mut self, seed: i64) {
        self.seed = (seed as u64 ^ MULTIPLIER) & MASK;
    }

    /// `Random.nextLong()`
    pub fn next_long(&mut self) -> i64 {
        ((self.next(32) as i64) << 32).wrapping_add(self.next(32) as i64)
    }

    /// `Random.nextBoolean()`
    pub fn next_boolean(&mut self) -> bool {
        self.next(1) != 0
    }

    /// `Random.nextInt(bound)` for `bound > 0`.
    pub fn next_int_bound(&mut self, bound: i32) -> i32 {
        debug_assert!(bound > 0);
        let m = bound - 1;
        let mut r = self.next(31);
        if bound & m == 0 {
            // power of two
            r = (((bound as i64).wrapping_mul(r as i64)) >> 31) as i32;
        } else {
            let mut u = r;
            loop {
                r = u % bound;
                if u.wrapping_sub(r).wrapping_add(m) >= 0 {
                    break;
                }
                u = self.next(31);
            }
        }
        r
    }
}

/// `blbutil.Utilities.shuffle(int[], Random)`
pub fn java_shuffle(ia: &mut [i32], rand: &mut JavaRandom) {
    for j in 0..ia.len() {
        let x = rand.next_int_bound((ia.len() - j) as i32) as usize;
        ia.swap(j, j + x);
    }
}

/// Trait mirroring `Comparable.compareTo` used by the priority queue.
pub trait JavaOrd {
    fn compare_to(&self, other: &Self) -> std::cmp::Ordering;
}

/// Port of `java.util.PriorityQueue`: an unbalanced binary min-heap whose
/// sift-up/sift-down order determines the placement of elements that
/// compare as equal.  Beagle's composite-haplotype construction depends on
/// this tie-breaking behavior, so we replicate it exactly.
pub struct JavaPriorityQueue<T: JavaOrd + Clone> {
    queue: Vec<T>,
}

impl<T: JavaOrd + Clone> JavaPriorityQueue<T> {
    pub fn new() -> Self {
        JavaPriorityQueue { queue: Vec::new() }
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.queue.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.queue.is_empty()
    }

    pub fn clear(&mut self) {
        self.queue.clear();
    }

    pub fn peek(&self) -> Option<&T> {
        self.queue.first()
    }

    pub fn offer(&mut self, x: T) {
        let k = self.queue.len();
        self.queue.push(x.clone()); // reserve the slot
        if k > 0 {
            self.sift_up(k, x);
        }
    }

    pub fn poll(&mut self) -> Option<T> {
        if self.queue.is_empty() {
            return None;
        }
        let result = self.queue[0].clone();
        let x = self.queue.pop().unwrap();
        if !self.queue.is_empty() {
            self.sift_down(0, x);
        }
        Some(result)
    }

    fn sift_up(&mut self, mut k: usize, x: T) {
        while k > 0 {
            let parent = (k - 1) >> 1;
            if x.compare_to(&self.queue[parent]) != std::cmp::Ordering::Less {
                break;
            }
            self.queue[k] = self.queue[parent].clone();
            k = parent;
        }
        self.queue[k] = x;
    }

    fn sift_down(&mut self, mut k: usize, x: T) {
        let size = self.queue.len();
        let half = size >> 1;
        while k < half {
            let mut child = (k << 1) + 1;
            let right = child + 1;
            if right < size
                && self.queue[child].compare_to(&self.queue[right])
                    == std::cmp::Ordering::Greater
            {
                child = right;
            }
            if x.compare_to(&self.queue[child]) != std::cmp::Ordering::Greater {
                break;
            }
            self.queue[k] = self.queue[child].clone();
            k = child;
        }
        self.queue[k] = x;
    }
}

/// Small open-addressing map from i32 keys to i32 values
/// (replacement for `ints.IntIntMap`; only lookup semantics matter).
pub struct IntIntMap {
    keys: Vec<i32>,
    vals: Vec<i32>,
    mask: usize,
    size: usize,
}

const EMPTY_KEY: i32 = i32::MIN;

impl IntIntMap {
    pub fn with_capacity(capacity: usize) -> Self {
        let cap = (capacity * 4).next_power_of_two().max(16);
        IntIntMap {
            keys: vec![EMPTY_KEY; cap],
            vals: vec![0; cap],
            mask: cap - 1,
            size: 0,
        }
    }

    #[inline]
    fn index_of(&self, key: i32) -> usize {
        let mut idx = (key as u32 as usize).wrapping_mul(0x9E3779B9) & self.mask;
        loop {
            let k = self.keys[idx];
            if k == key || k == EMPTY_KEY {
                return idx;
            }
            idx = (idx + 1) & self.mask;
        }
    }

    #[inline]
    pub fn get(&self, key: i32, default: i32) -> i32 {
        let idx = self.index_of(key);
        if self.keys[idx] == key {
            self.vals[idx]
        } else {
            default
        }
    }

    pub fn put(&mut self, key: i32, value: i32) {
        debug_assert!(key != EMPTY_KEY);
        let idx = self.index_of(key);
        if self.keys[idx] != key {
            self.keys[idx] = key;
            self.size += 1;
            debug_assert!(self.size * 2 <= self.keys.len());
        }
        self.vals[idx] = value;
    }

    /// Removal with backward-shift deletion to keep probe chains valid.
    pub fn remove(&mut self, key: i32) {
        let mut idx = self.index_of(key);
        if self.keys[idx] != key {
            return;
        }
        self.keys[idx] = EMPTY_KEY;
        self.size -= 1;
        // re-insert the rest of the cluster
        let mut next = (idx + 1) & self.mask;
        while self.keys[next] != EMPTY_KEY {
            let k = self.keys[next];
            let v = self.vals[next];
            self.keys[next] = EMPTY_KEY;
            let new_idx = self.index_of(k);
            self.keys[new_idx] = k;
            self.vals[new_idx] = v;
            idx = next;
            next = (idx + 1) & self.mask;
        }
    }

    pub fn clear(&mut self) {
        if self.size > 0 {
            self.keys.fill(EMPTY_KEY);
            self.size = 0;
        }
    }
}

/// `DS_VALS` from `imp.ImputedRecBuilder`: `DecimalFormat("#.##")` of j/100
/// for j in 0..=200.
pub fn ds_vals() -> Vec<String> {
    (0..=200)
        .map(|j: u32| {
            if j % 100 == 0 {
                format!("{}", j / 100)
            } else if j % 10 == 0 {
                format!("{}.{}", j / 100, (j % 100) / 10)
            } else {
                format!("{}.{:02}", j / 100, j % 100)
            }
        })
        .collect()
}

/// `R2_VALS` from `imp.ImputedRecBuilder`: two-decimal formatting of i/100.
pub fn r2_vals() -> Vec<String> {
    let ds = ds_vals();
    (0..=100usize)
        .map(|i| {
            if ds[i].len() != 4 {
                format!("{}.{:02}", i / 100, i % 100)
            } else {
                ds[i].clone()
            }
        })
        .collect()
}

/// `Math.rint`: round half to even.
#[inline]
pub fn java_rint(x: f64) -> f64 {
    x.round_ties_even()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_random_matches_reference() {
        // Reference values generated by java.util.Random(42):
        // nextInt() x5: -1170105035, 234785527, -1360544799, 205897768, 1325939940
        let mut r = JavaRandom::new(42);
        assert_eq!(r.next_int(), -1170105035);
        assert_eq!(r.next_int(), 234785527);
        assert_eq!(r.next_int(), -1360544799);
        assert_eq!(r.next_int(), 205897768);
        assert_eq!(r.next_int(), 1325939940);
        // java.util.Random(42).nextInt(100) x5: 30, 63, 48, 84, 70
        let mut r = JavaRandom::new(42);
        assert_eq!(r.next_int_bound(100), 30);
        assert_eq!(r.next_int_bound(100), 63);
        assert_eq!(r.next_int_bound(100), 48);
        assert_eq!(r.next_int_bound(100), 84);
        assert_eq!(r.next_int_bound(100), 70);
        // java.util.Random(-99999).nextInt() x3
        let mut r = JavaRandom::new(-99999);
        assert_eq!(r.next_int(), 1971967714);
        assert_eq!(r.next_int(), -411896953);
        assert_eq!(r.next_int(), -249951563);
    }

    #[test]
    fn ds_val_table() {
        let ds = ds_vals();
        assert_eq!(ds[0], "0");
        assert_eq!(ds[5], "0.05");
        assert_eq!(ds[50], "0.5");
        assert_eq!(ds[100], "1");
        assert_eq!(ds[155], "1.55");
        assert_eq!(ds[200], "2");
        let r2 = r2_vals();
        assert_eq!(r2[0], "0.00");
        assert_eq!(r2[10], "0.10");
        assert_eq!(r2[42], "0.42");
        assert_eq!(r2[100], "1.00");
    }
}
