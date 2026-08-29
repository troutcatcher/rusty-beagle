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
