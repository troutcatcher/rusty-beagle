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

//! Port of `blbutil.BitArray`: a fixed-size bit list backed by u64 words.
//! The XOR-based range hash depends on the exact word layout, so all
//! operations replicate Java's bit order (bit i lives in word i>>6 at
//! position i&63).

#[derive(Clone)]
pub struct BitArray {
    words: Vec<u64>,
    size: usize,
}

const WORD_MASK: u64 = u64::MAX;

/// Java `WORD_MASK << from` / `WORD_MASK >>> -to` use the shift count mod 64.
#[inline]
fn shl(mask: u64, n: usize) -> u64 {
    mask << (n & 63)
}

#[inline]
fn shr_neg(mask: u64, to: usize) -> u64 {
    // Java: WORD_MASK >>> -to  ==  WORD_MASK >>> ((64 - to) & 63)
    mask >> ((64usize.wrapping_sub(to)) & 63)
}

impl BitArray {
    pub fn new(size: usize) -> BitArray {
        BitArray {
            words: vec![0; (size + 63) >> 6],
            size,
        }
    }

    #[inline]
    #[allow(dead_code)]
    pub fn size(&self) -> usize {
        self.size
    }

    #[inline]
    pub fn get(&self, index: usize) -> bool {
        debug_assert!(index < self.size);
        (self.words[index >> 6] >> (index & 63)) & 1 != 0
    }

    #[inline]
    pub fn set(&mut self, index: usize) {
        debug_assert!(index < self.size);
        self.words[index >> 6] |= 1u64 << (index & 63);
    }

    #[inline]
    pub fn clear_bit(&mut self, index: usize) {
        debug_assert!(index < self.size);
        self.words[index >> 6] &= !(1u64 << (index & 63));
    }

    /// `BitArray.hash(from, to)`: XOR of the masked words, folded to i32.
    pub fn hash(&self, from: usize, to: usize) -> i32 {
        debug_assert!(from <= to && to <= self.size);
        if from == to {
            return 0;
        }
        let start_word = from >> 6;
        let end_word = (to - 1) >> 6;
        let start_word_mask = shl(WORD_MASK, from);
        let end_word_mask = shr_neg(WORD_MASK, to);
        let long_hash = if start_word == end_word {
            self.words[start_word] & (start_word_mask & end_word_mask)
        } else {
            let mut h = self.words[start_word] & start_word_mask;
            for j in start_word + 1..end_word {
                h ^= self.words[j];
            }
            h ^= self.words[end_word] & end_word_mask;
            h
        };
        long_hash_code(long_hash)
    }

    /// `BitArray.equal(other, from, to)`: bit-range equality.
    pub fn equal(&self, other: &BitArray, from: usize, to: usize) -> bool {
        debug_assert!(from <= to && to <= self.size);
        if from == to {
            return true;
        }
        let start_word = from >> 6;
        let end_word = (to - 1) >> 6;
        let start_word_mask = shl(WORD_MASK, from);
        let end_word_mask = shr_neg(WORD_MASK, to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            (self.words[start_word] ^ other.words[start_word]) & mask == 0
        } else {
            if (self.words[start_word] ^ other.words[start_word]) & start_word_mask != 0 {
                return false;
            }
            for j in start_word + 1..end_word {
                if self.words[j] != other.words[j] {
                    return false;
                }
            }
            (self.words[end_word] ^ other.words[end_word]) & end_word_mask == 0
        }
    }

    /// `BitArray.copyFrom(src, from, to)`: replaces bits [from, to) in `self`
    /// with the corresponding bits of `src`.
    pub fn copy_from(&mut self, src: &BitArray, from: usize, to: usize) {
        debug_assert!(from <= to && to <= self.size && to <= src.size);
        if from == to {
            return;
        }
        let start_word = from >> 6;
        let end_word = (to - 1) >> 6;
        let start_word_mask = shl(WORD_MASK, from);
        let end_word_mask = shr_neg(WORD_MASK, to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            self.words[start_word] ^= (self.words[start_word] ^ src.words[start_word]) & mask;
        } else {
            self.words[start_word] ^=
                (self.words[start_word] ^ src.words[start_word]) & start_word_mask;
            self.words[start_word + 1..end_word]
                .copy_from_slice(&src.words[start_word + 1..end_word]);
            self.words[end_word] ^= (self.words[end_word] ^ src.words[end_word]) & end_word_mask;
        }
    }

    /// `BitArray.swapBits(a, b, from, to)`
    pub fn swap_bits(a: &mut BitArray, b: &mut BitArray, from: usize, to: usize) {
        debug_assert_eq!(a.size, b.size);
        debug_assert!(from <= to && to <= a.size);
        if from == to {
            return;
        }
        let start_word = from >> 6;
        let end_word = (to - 1) >> 6;
        let start_word_mask = shl(WORD_MASK, from);
        let end_word_mask = shr_neg(WORD_MASK, to);
        if start_word == end_word {
            let mask = start_word_mask & end_word_mask;
            let diff = (a.words[start_word] ^ b.words[start_word]) & mask;
            a.words[start_word] ^= diff;
            b.words[start_word] ^= diff;
        } else {
            let diff = (a.words[start_word] ^ b.words[start_word]) & start_word_mask;
            a.words[start_word] ^= diff;
            b.words[start_word] ^= diff;
            for j in start_word + 1..end_word {
                std::mem::swap(&mut a.words[j], &mut b.words[j]);
            }
            let diff = (a.words[end_word] ^ b.words[end_word]) & end_word_mask;
            a.words[end_word] ^= diff;
            b.words[end_word] ^= diff;
        }
    }
}

#[inline]
pub fn long_hash_code(value: u64) -> i32 {
    (value ^ (value >> 32)) as u32 as i32
}
