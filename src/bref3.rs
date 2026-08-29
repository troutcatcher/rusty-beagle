//! Reads a bref3 (binary reference format version 3) file as a reference
//! panel. Port of `bref.Bref3Header` + `bref.Bref3Reader` + `bref.Bref3It`.
//!
//! Unlike VCF input, a bref3 file already encodes its own sequence-coding
//! block structure (shared per-block hap->seq map, per-record seq->allele
//! map or, for records that didn't compress, explicit per-allele haplotype
//! lists) exactly as chosen by the `bref3` conversion tool. rusty-beagle
//! therefore does not re-run `SeqCoder3` for bref3 input: it decodes each
//! on-disk block verbatim and assigns it a fresh block id, which reproduces
//! Java's `rec.map(0) != hapToSeq` block-boundary identity check (every
//! marker decoded from one `readBlock()` call shares one block id, exactly
//! as they share one `hapToSeq` object reference in Java).
//!
//! The bref3 index/footer (written for random access by other tools) is
//! never consulted: `Bref3It` in Java reads blocks strictly sequentially
//! from the start of the file until the zero-record end-of-data sentinel,
//! and so do we.

use crate::marker::{ChromIds, Marker};
use crate::refpanel::{RefAlleles, RefRec, RefSource};
use crate::vcfio::Samples;
use std::collections::{HashSet, VecDeque};
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::sync::Arc;

const MAGIC_NUMBER_V3: i32 = 2055763188;

/// 24 permutations of "A","C","G","T" in the same (lexicographic) order
/// produced by `Bref3Reader.snvPerms()`'s recursive generator, used to
/// decode the compact SNV allele coding (`Bref3Reader.alleleString`).
const SNV_PERMS: [[&str; 4]; 24] = [
    ["A", "C", "G", "T"],
    ["A", "C", "T", "G"],
    ["A", "G", "C", "T"],
    ["A", "G", "T", "C"],
    ["A", "T", "C", "G"],
    ["A", "T", "G", "C"],
    ["C", "A", "G", "T"],
    ["C", "A", "T", "G"],
    ["C", "G", "A", "T"],
    ["C", "G", "T", "A"],
    ["C", "T", "A", "G"],
    ["C", "T", "G", "A"],
    ["G", "A", "C", "T"],
    ["G", "A", "T", "C"],
    ["G", "C", "A", "T"],
    ["G", "C", "T", "A"],
    ["G", "T", "A", "C"],
    ["G", "T", "C", "A"],
    ["T", "A", "C", "G"],
    ["T", "A", "G", "C"],
    ["T", "C", "A", "G"],
    ["T", "C", "G", "A"],
    ["T", "G", "A", "C"],
    ["T", "G", "C", "A"],
];

fn io_exit(path: &str, err: io::Error) -> ! {
    eprintln!("ERROR reading file {}: {}", path, err);
    std::process::exit(1)
}

/// Decodes Java's "modified UTF-8" (`DataInput.readUTF`'s payload format):
/// like UTF-8 except NUL is coded as the two bytes 0xC0,0x80, and
/// characters outside the BMP are coded as a pair of independently
/// 3-byte-encoded UTF-16 surrogates instead of one 4-byte sequence.
fn decode_modified_utf8(bytes: &[u8]) -> String {
    let mut out = String::with_capacity(bytes.len());
    let mut i = 0;
    let mut pending_high: Option<u32> = None;
    while i < bytes.len() {
        let b0 = bytes[i] as u32;
        let (cp, len) = if b0 & 0x80 == 0 {
            (b0, 1)
        } else if b0 & 0xE0 == 0xC0 && i + 1 < bytes.len() {
            (((b0 & 0x1F) << 6) | (bytes[i + 1] as u32 & 0x3F), 2)
        } else if b0 & 0xF0 == 0xE0 && i + 2 < bytes.len() {
            (
                ((b0 & 0x0F) << 12)
                    | ((bytes[i + 1] as u32 & 0x3F) << 6)
                    | (bytes[i + 2] as u32 & 0x3F),
                3,
            )
        } else {
            (0xFFFD, 1)
        };
        i += len;
        if (0xD800..0xDC00).contains(&cp) {
            pending_high = Some(cp);
            continue;
        } else if (0xDC00..0xE000).contains(&cp) {
            if let Some(hi) = pending_high.take() {
                let combined = 0x10000 + ((hi - 0xD800) << 10) + (cp - 0xDC00);
                if let Some(c) = char::from_u32(combined) {
                    out.push(c);
                }
                continue;
            }
        }
        pending_high = None;
        if let Some(c) = char::from_u32(cp) {
            out.push(c);
        }
    }
    out
}

/// `RefGTRec.alleleToHaps` doc: joins REF then ALT alleles (comma-separated,
/// "." when there is no ALT) the way `Marker`'s "REF\tALT" field is stored.
fn ref_alt_fields(alleles: &[String]) -> Arc<str> {
    let mut s = String::with_capacity(8 + alleles.len() * 2);
    s.push_str(&alleles[0]);
    s.push('\t');
    if alleles.len() == 1 {
        s.push('.');
    } else {
        for (i, a) in alleles[1..].iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(a);
        }
    }
    Arc::from(s.as_str())
}

/// Reads one bref3 file sequentially, applying sample/marker filters and
/// exposing decoded records as `RefRec`s (the same type the VCF reference
/// reader in `refpanel.rs` produces).
pub struct Bref3RefReader {
    r: BufReader<File>,
    path: String,
    /// original (unfiltered) haplotype index for each kept haplotype, in
    /// increasing order; `Bref3Header.filteredHapIndices()`
    included_hap_indices: Vec<u32>,
    /// original haplotype index -> kept haplotype index, or -1 if excluded;
    /// `Bref3Header.invfilteredHapIndices()`
    inv_included_hap_indices: Vec<i32>,
    samples: Samples,
    n_haps_filtered: usize,
    raw_hap2seq_buf: Vec<u8>,
    pending: VecDeque<Arc<RefRec>>,
    next_block_id: u32,
    exclude_markers: HashSet<String>,
    eof: bool,
}

impl Bref3RefReader {
    pub fn new(
        path: &str,
        exclude_samples: &HashSet<String>,
        exclude_markers: HashSet<String>,
    ) -> Bref3RefReader {
        let f = File::open(path).unwrap_or_else(|e| {
            eprintln!("ERROR: cannot open file {}: {}", path, e);
            std::process::exit(1)
        });
        let mut r = BufReader::with_capacity(1 << 20, f);
        let magic = read_i32(&mut r).unwrap_or_else(|e| io_exit(path, e));
        if magic != MAGIC_NUMBER_V3 {
            eprintln!(
                "ERROR: Unrecognized input file.  Was input file created\n\
                 with a different version of the bref program?\n\
                 File: {}",
                path
            );
            std::process::exit(1);
        }
        let _program = read_utf(&mut r).unwrap_or_else(|e| io_exit(path, e));
        let sample_ids = read_string_array(&mut r).unwrap_or_else(|e| io_exit(path, e));

        let mut included_hap_indices = Vec::new();
        let mut ids = Vec::new();
        for (j, id) in sample_ids.iter().enumerate() {
            if !exclude_samples.contains(id) {
                included_hap_indices.push((2 * j) as u32);
                included_hap_indices.push((2 * j + 1) as u32);
                ids.push(id.clone());
            }
        }
        if ids.is_empty() {
            eprintln!(
                "ERROR: All samples in the bref3 file have been excluded\n\
                 Bref3 file :  {}",
                path
            );
            std::process::exit(1);
        }
        let n_haps_unfiltered = sample_ids.len() << 1;
        let mut inv_included_hap_indices = vec![-1i32; n_haps_unfiltered];
        for (k, &h) in included_hap_indices.iter().enumerate() {
            inv_included_hap_indices[h as usize] = k as i32;
        }
        let n_haps_filtered = ids.len() << 1;
        let samples = Samples {
            is_diploid: Arc::new(vec![true; ids.len()]),
            ids: Arc::new(ids),
        };

        Bref3RefReader {
            r,
            path: path.to_string(),
            included_hap_indices,
            inv_included_hap_indices,
            samples,
            n_haps_filtered,
            raw_hap2seq_buf: vec![0u8; 2 * n_haps_unfiltered],
            pending: VecDeque::new(),
            next_block_id: 0,
            exclude_markers,
            eof: false,
        }
    }

    pub fn samples(&self) -> Samples {
        self.samples.clone()
    }

    fn accept_marker(&self, marker: &Marker) -> bool {
        match &marker.id {
            Some(id) => !self.exclude_markers.contains(id.as_ref()),
            None => true,
        }
    }

    pub fn next(&mut self) -> Option<Arc<RefRec>> {
        while self.pending.is_empty() && !self.eof {
            self.read_block();
        }
        self.pending.pop_front()
    }

    /// Reads one bref3 data block (or the end-of-data sentinel) and pushes
    /// any markers that pass the marker filter onto `pending`.
    fn read_block(&mut self) {
        let n_recs = self.read_i32();
        if n_recs == 0 {
            self.eof = true;
            return;
        }
        let n_recs = n_recs as usize;
        let chrom = self.read_utf();
        let chrom_idx = ChromIds::instance().get_index(&chrom);
        let n_seq = self.read_u16() as usize;
        self.r
            .read_exact(&mut self.raw_hap2seq_buf)
            .unwrap_or_else(|e| io_exit(&self.path, e));
        let hap2seq: Arc<Vec<u16>> = Arc::new(
            self.included_hap_indices
                .iter()
                .map(|&h| {
                    let off = (h as usize) << 1;
                    ((self.raw_hap2seq_buf[off] as u16) << 8)
                        | (self.raw_hap2seq_buf[off + 1] as u16)
                })
                .collect(),
        );
        let block_id = self.next_block_id;
        self.next_block_id += 1;
        for _ in 0..n_recs {
            let marker = self.read_marker(chrom_idx);
            let flag = self.read_i8();
            let alleles = match flag {
                0 => {
                    let mut seq2allele = vec![0u8; n_seq];
                    self.r
                        .read_exact(&mut seq2allele)
                        .unwrap_or_else(|e| io_exit(&self.path, e));
                    RefAlleles::SeqCoded {
                        block: block_id,
                        hap2seq: hap2seq.clone(),
                        seq2allele,
                    }
                }
                1 => self.read_allele_coded_rec(marker.n_alleles as usize),
                _ => {
                    eprintln!("ERROR reading file {}: unrecognized record flag", self.path);
                    std::process::exit(1)
                }
            };
            if self.accept_marker(&marker) {
                self.pending.push_back(Arc::new(RefRec {
                    marker,
                    alleles,
                    n_haps: self.n_haps_filtered,
                }));
            }
        }
    }

    fn read_allele_coded_rec(&mut self, n_alleles: usize) -> RefAlleles {
        let mut carriers: Vec<Vec<u32>> = Vec::with_capacity(n_alleles);
        let mut major = 0u8;
        for a in 0..n_alleles {
            let len = self.read_i32();
            if len == -1 {
                major = a as u8;
                carriers.push(Vec::new());
            } else {
                // The file stores each allele's carrier haplotypes in
                // increasing original-hap-index order; `inv_included_hap_indices`
                // is order-preserving on its included domain (filtered samples
                // keep their relative original order), so the remapped list
                // stays sorted without an explicit re-sort.
                let mut list = Vec::with_capacity(len as usize);
                for _ in 0..len {
                    let hap = self.read_i32() as usize;
                    let filtered = self.inv_included_hap_indices[hap];
                    if filtered >= 0 {
                        list.push(filtered as u32);
                    }
                }
                carriers.push(list);
            }
        }
        RefAlleles::AlleleCoded { major, carriers }
    }

    fn read_marker(&mut self, chrom_idx: u16) -> Marker {
        let pos = self.read_i32();
        let id = self.read_marker_id();
        let allele_code = self.read_i8() as i32;
        if allele_code == -1 {
            let str_alleles = self.read_string_array();
            let n_alleles = str_alleles.len() as u16;
            let alleles = ref_alt_fields(&str_alleles);
            let end_val = self.read_i32();
            let end = if end_val >= 0 {
                Some(Arc::from(end_val.to_string().as_str()))
            } else {
                None
            };
            Marker {
                chrom_idx,
                pos,
                id,
                alleles,
                n_alleles,
                end,
            }
        } else {
            let n_alleles = 1 + (allele_code & 0b11) as u16;
            let perm_index = (allele_code >> 2) as usize;
            let str_alleles: Vec<String> = SNV_PERMS[perm_index][..n_alleles as usize]
                .iter()
                .map(|s| s.to_string())
                .collect();
            let alleles = ref_alt_fields(&str_alleles);
            Marker {
                chrom_idx,
                pos,
                id,
                alleles,
                n_alleles,
                end: None,
            }
        }
    }

    fn read_marker_id(&mut self) -> Option<Arc<str>> {
        let n_ids = self.read_u8();
        if n_ids == 0 {
            return None;
        }
        let mut s = String::new();
        for i in 0..n_ids {
            if i > 0 {
                s.push(';');
            }
            s.push_str(&self.read_utf());
        }
        Some(Arc::from(s.as_str()))
    }

    fn read_string_array(&mut self) -> Vec<String> {
        let len = self.read_i32();
        if len <= 0 {
            return Vec::new();
        }
        (0..len).map(|_| self.read_utf()).collect()
    }

    fn read_i8(&mut self) -> i8 {
        let mut b = [0u8; 1];
        self.r
            .read_exact(&mut b)
            .unwrap_or_else(|e| io_exit(&self.path, e));
        b[0] as i8
    }

    fn read_u8(&mut self) -> u8 {
        let mut b = [0u8; 1];
        self.r
            .read_exact(&mut b)
            .unwrap_or_else(|e| io_exit(&self.path, e));
        b[0]
    }

    fn read_u16(&mut self) -> u16 {
        let mut b = [0u8; 2];
        self.r
            .read_exact(&mut b)
            .unwrap_or_else(|e| io_exit(&self.path, e));
        u16::from_be_bytes(b)
    }

    fn read_i32(&mut self) -> i32 {
        read_i32(&mut self.r).unwrap_or_else(|e| io_exit(&self.path, e))
    }

    fn read_utf(&mut self) -> String {
        read_utf(&mut self.r).unwrap_or_else(|e| io_exit(&self.path, e))
    }
}

impl RefSource for Bref3RefReader {
    fn next_rec(&mut self) -> Option<Arc<RefRec>> {
        self.next()
    }

    fn samples(&self) -> Samples {
        self.samples()
    }
}

fn read_i32(r: &mut impl Read) -> io::Result<i32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_be_bytes(b))
}

fn read_utf(r: &mut impl Read) -> io::Result<String> {
    let mut len_b = [0u8; 2];
    r.read_exact(&mut len_b)?;
    let len = u16::from_be_bytes(len_b) as usize;
    let mut buf = vec![0u8; len];
    r.read_exact(&mut buf)?;
    Ok(decode_modified_utf8(&buf))
}

fn read_string_array(r: &mut impl Read) -> io::Result<Vec<String>> {
    let len = read_i32(r)?;
    if len <= 0 {
        return Ok(Vec::new());
    }
    (0..len).map(|_| read_utf(r)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn modified_utf8_ascii() {
        assert_eq!(decode_modified_utf8(b"chr21"), "chr21");
        assert_eq!(decode_modified_utf8(b""), "");
    }

    #[test]
    fn modified_utf8_embedded_nul() {
        // Java encodes U+0000 as the two bytes 0xC0,0x80 rather than 0x00.
        assert_eq!(decode_modified_utf8(&[b'a', 0xC0, 0x80, b'b']), "a\u{0}b");
    }

    #[test]
    fn modified_utf8_two_byte() {
        // U+00E9 ('e' with acute accent), standard 2-byte UTF-8.
        assert_eq!(decode_modified_utf8(&[0xC3, 0xA9]), "\u{E9}");
    }

    #[test]
    fn modified_utf8_surrogate_pair() {
        // U+1F600 (an emoji outside the BMP) is coded as a CESU-8 pair: the
        // high and low UTF-16 surrogates, each independently 3-byte-encoded,
        // rather than one 4-byte UTF-8 sequence.
        let bytes = [0xED, 0xA0, 0xBD, 0xED, 0xB8, 0x80];
        assert_eq!(decode_modified_utf8(&bytes), "\u{1F600}");
    }

    #[test]
    fn snv_perms_match_java_generator() {
        // Cross-check the hardcoded table against Bref3Reader.snvPerms()'s
        // recursive "pick next remaining element in order" generator, ported
        // directly (rather than trusting the by-hand transcription alone).
        fn permute(
            start: &mut Vec<&'static str>,
            end: &[&'static str],
            out: &mut Vec<[&'static str; 4]>,
        ) {
            if end.is_empty() {
                let mut arr = ["", "", "", ""];
                arr[..start.len()].copy_from_slice(start);
                out.push(arr);
            } else {
                for j in 0..end.len() {
                    start.push(end[j]);
                    let mut new_end: Vec<&'static str> = end[..j].to_vec();
                    new_end.extend_from_slice(&end[j + 1..]);
                    permute(start, &new_end, out);
                    start.pop();
                }
            }
        }
        let mut generated = Vec::with_capacity(24);
        permute(&mut Vec::new(), &["A", "C", "G", "T"], &mut generated);
        assert_eq!(generated.len(), 24);
        assert_eq!(generated, SNV_PERMS);
    }
}
