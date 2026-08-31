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

//! Minimal BGZF (blocked gzip) writer, equivalent in framing to
//! `blbutil.BGZIPOutputStream`.  Output is a series of independent gzip
//! members with the BC extra field; readers (including plain gunzip)
//! decompress the concatenation transparently.

use flate2::write::DeflateEncoder;
use flate2::Compression;
use std::io::Write;

/// Maximum uncompressed bytes per BGZF block (kept comfortably below the
/// 65535-byte compressed-block limit even for incompressible data).
const BLOCK_SIZE: usize = 0xff00;

/// Compresses `data` into a sequence of BGZF blocks.
pub fn compress(data: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(data.len() / 3 + 64);
    for chunk in data.chunks(BLOCK_SIZE.max(1)) {
        write_block(chunk, &mut out);
    }
    out
}

/// The 28-byte BGZF end-of-file marker (an empty block).
pub fn eof_block() -> Vec<u8> {
    let mut out = Vec::with_capacity(28);
    write_block(&[], &mut out);
    out
}

fn write_block(chunk: &[u8], out: &mut Vec<u8>) {
    let mut deflated = Vec::with_capacity(chunk.len() / 2 + 32);
    {
        let mut enc = DeflateEncoder::new(&mut deflated, Compression::default());
        enc.write_all(chunk).expect("deflate write");
        enc.finish().expect("deflate finish");
    }
    let bsize = 12 + 6 + deflated.len() + 8; // header + XLEN(BC) + data + footer
    assert!(bsize <= 65536, "BGZF block too large");
    let bsize_m1 = (bsize - 1) as u16;
    // gzip header with FEXTRA
    out.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, // magic, deflate, FLG.FEXTRA
        0, 0, 0, 0, // MTIME
        0, 0xff, // XFL, OS=unknown
        6, 0, // XLEN = 6
        b'B', b'C', 2, 0, // BC subfield, length 2
    ]);
    out.extend_from_slice(&bsize_m1.to_le_bytes());
    out.extend_from_slice(&deflated);
    let crc = crc32(chunk);
    out.extend_from_slice(&crc.to_le_bytes());
    out.extend_from_slice(&(chunk.len() as u32).to_le_bytes());
}

fn crc32(data: &[u8]) -> u32 {
    let mut hasher = flate2::Crc::new();
    hasher.update(data);
    hasher.sum()
}

// ---------------------------------------------------------------------------
// Parallel BGZF reader (equivalent to Java's BGZipIt): BGZF files record the
// compressed block size in a BC extra subfield, so blocks can be located
// without inflating and decompressed in parallel.

use std::io::Read;
use std::sync::mpsc::{Receiver, SyncSender};

/// Returns true when `header` starts a BGZF member (gzip + FEXTRA + BC).
pub fn sniff_bgzf(header: &[u8]) -> bool {
    header.len() >= 18
        && header[0] == 0x1f
        && header[1] == 0x8b
        && header[2] == 8
        && (header[3] & 4) != 0
        && find_bc(&header[12..12 + header[10] as usize + ((header[11] as usize) << 8)])
            .is_some()
}

/// Finds the BC subfield in a gzip FEXTRA payload; returns BSIZE (total
/// member length - 1).
fn find_bc(extra: &[u8]) -> Option<usize> {
    let mut i = 0;
    while i + 4 <= extra.len() {
        let slen = extra[i + 2] as usize + ((extra[i + 3] as usize) << 8);
        if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 && i + 6 <= extra.len() {
            return Some(extra[i + 4] as usize + ((extra[i + 5] as usize) << 8));
        }
        i += 4 + slen;
    }
    None
}

/// A `Read` implementation that decompresses BGZF blocks in parallel on the
/// rayon pool, keeping a bounded number of inflated blocks in flight.
pub struct ParallelBgzfReader {
    rx: Receiver<std::io::Result<Vec<u8>>>,
    cur: Vec<u8>,
    pos: usize,
    done: bool,
}

const BLOCKS_PER_BATCH: usize = 64;

impl ParallelBgzfReader {
    pub fn new<R: Read + Send + 'static>(mut src: R) -> ParallelBgzfReader {
        let (tx, rx) = std::sync::mpsc::sync_channel::<std::io::Result<Vec<u8>>>(4);
        std::thread::spawn(move || {
            let _ = Self::pump(&mut src, &tx);
        });
        ParallelBgzfReader {
            rx,
            cur: Vec::new(),
            pos: 0,
            done: false,
        }
    }

    fn pump<R: Read>(
        src: &mut R,
        tx: &SyncSender<std::io::Result<Vec<u8>>>,
    ) -> Result<(), ()> {
        use rayon::prelude::*;
        loop {
            // read up to BLOCKS_PER_BATCH raw blocks
            let mut raw_blocks: Vec<Vec<u8>> = Vec::with_capacity(BLOCKS_PER_BATCH);
            for _ in 0..BLOCKS_PER_BATCH {
                match read_raw_block(src) {
                    Ok(Some(b)) => raw_blocks.push(b),
                    Ok(None) => break,
                    Err(e) => {
                        let _ = tx.send(Err(e));
                        return Err(());
                    }
                }
            }
            if raw_blocks.is_empty() {
                return Ok(()); // EOF; dropping tx closes the channel
            }
            let inflated: Vec<std::io::Result<Vec<u8>>> = raw_blocks
                .par_iter()
                .map(|b| inflate_block(b))
                .collect();
            for block in inflated {
                let is_err = block.is_err();
                if tx.send(block).is_err() || is_err {
                    return Err(());
                }
            }
        }
    }
}

/// Reads one complete BGZF member (header through footer), or None at EOF.
fn read_raw_block<R: Read>(src: &mut R) -> std::io::Result<Option<Vec<u8>>> {
    let mut header = [0u8; 12];
    let mut got = 0;
    while got < header.len() {
        let n = src.read(&mut header[got..])?;
        if n == 0 {
            if got == 0 {
                return Ok(None);
            }
            return Err(bad_bgzf("truncated BGZF header"));
        }
        got += n;
    }
    if header[0] != 0x1f || header[1] != 0x8b || header[2] != 8 || (header[3] & 4) == 0 {
        return Err(bad_bgzf("not a BGZF block"));
    }
    let xlen = header[10] as usize + ((header[11] as usize) << 8);
    let mut extra = vec![0u8; xlen];
    src.read_exact(&mut extra).map_err(|_| bad_bgzf("truncated BGZF extra field"))?;
    let bsize = find_bc(&extra).ok_or_else(|| bad_bgzf("missing BGZF BC subfield"))? + 1;
    let rest_len = bsize
        .checked_sub(12 + xlen)
        .ok_or_else(|| bad_bgzf("invalid BGZF BSIZE"))?;
    let mut block = Vec::with_capacity(bsize);
    block.extend_from_slice(&header);
    block.extend_from_slice(&extra);
    let start = block.len();
    block.resize(bsize, 0);
    src.read_exact(&mut block[start..])
        .map_err(|_| bad_bgzf("truncated BGZF block"))?;
    let _ = rest_len;
    Ok(Some(block))
}

fn inflate_block(block: &[u8]) -> std::io::Result<Vec<u8>> {
    let xlen = block[10] as usize + ((block[11] as usize) << 8);
    let data_start = 12 + xlen;
    let data_end = block.len() - 8;
    let isize_bytes: [u8; 4] = block[block.len() - 4..].try_into().unwrap();
    let isize = u32::from_le_bytes(isize_bytes) as usize;
    let mut out = vec![0u8; isize];
    if isize > 0 {
        let mut d = flate2::Decompress::new(false);
        d.decompress(
            &block[data_start..data_end],
            &mut out,
            flate2::FlushDecompress::Finish,
        )
        .map_err(|_| bad_bgzf("corrupt BGZF block"))?;
        if d.total_out() as usize != isize {
            return Err(bad_bgzf("BGZF block length mismatch"));
        }
    }
    Ok(out)
}

fn bad_bgzf(msg: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, msg.to_string())
}

impl Read for ParallelBgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        loop {
            if self.pos < self.cur.len() {
                let n = (self.cur.len() - self.pos).min(buf.len());
                buf[..n].copy_from_slice(&self.cur[self.pos..self.pos + n]);
                self.pos += n;
                return Ok(n);
            }
            if self.done {
                return Ok(0);
            }
            match self.rx.recv() {
                Ok(Ok(block)) => {
                    self.cur = block;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => {
                    self.done = true;
                    return Ok(0);
                }
            }
        }
    }
}
