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

//! Port of `vcf.Marker` / `vcf.MarkerParser` (with `storeId=true`,
//! `storeQual=storeFilter=storeInfo=false`, INFO/END always kept) and the
//! global chromosome-id registry (`beagleutil.ChromIds`).

use std::sync::{Arc, Mutex, OnceLock};

/// Global registry of chromosome identifiers (order of first appearance).
pub struct ChromIds {
    inner: Mutex<ChromIdsInner>,
}

struct ChromIdsInner {
    ids: Vec<Arc<str>>,
    map: std::collections::HashMap<Arc<str>, u16>,
}

static CHROM_IDS: OnceLock<ChromIds> = OnceLock::new();

impl ChromIds {
    pub fn instance() -> &'static ChromIds {
        CHROM_IDS.get_or_init(|| ChromIds {
            inner: Mutex::new(ChromIdsInner {
                ids: Vec::new(),
                map: std::collections::HashMap::new(),
            }),
        })
    }

    pub fn get_index(&self, chrom: &str) -> u16 {
        let mut inner = self.inner.lock().unwrap();
        if let Some(&idx) = inner.map.get(chrom) {
            return idx;
        }
        let idx = inner.ids.len() as u16;
        let arc: Arc<str> = Arc::from(chrom);
        inner.ids.push(arc.clone());
        inner.map.insert(arc, idx);
        idx
    }

    pub fn id(&self, index: u16) -> Arc<str> {
        let inner = self.inner.lock().unwrap();
        inner.ids[index as usize].clone()
    }
}

/// Port of `vcf.Marker`.  Equality and ordering use
/// (chromIndex, pos, alleles, END value), matching Java.
#[derive(Clone, Debug)]
pub struct Marker {
    pub chrom_idx: u16,
    pub pos: i32,
    /// VCF ID field, `None` when missing (".")
    pub id: Option<Arc<str>>,
    /// tab-separated REF and ALT fields, exactly as in the input record
    pub alleles: Arc<str>,
    pub n_alleles: u16,
    /// INFO/END subfield value (without the "END=" prefix), if present
    pub end: Option<Arc<str>>,
}

impl PartialEq for Marker {
    fn eq(&self, other: &Self) -> bool {
        self.chrom_idx == other.chrom_idx
            && self.pos == other.pos
            && self.alleles == other.alleles
            && self.end == other.end
    }
}
impl Eq for Marker {}

impl Marker {
    /// Parse the first 8 fields of a VCF record.
    /// Returns the marker and the byte offset of the tab that ends field 8
    /// (i.e. the tab before the FORMAT field), or the record length if there
    /// are exactly 8 fields.
    pub fn parse(rec: &str) -> Result<(Marker, usize), String> {
        let bytes = rec.as_bytes();
        let mut tabs = [0usize; 8];
        let mut n = 0;
        for (i, &b) in bytes.iter().enumerate() {
            if b == b'\t' {
                tabs[n] = i;
                n += 1;
                if n == 8 {
                    break;
                }
            }
        }
        if n < 8 {
            return Err(format!(
                "VCF record does not contain at least 8 tab characters: {}",
                truncate(rec, 200)
            ));
        }
        let chrom = &rec[..tabs[0]];
        if chrom.is_empty() || chrom.contains(char::is_whitespace) {
            return Err(format!("invalid CHROM field: {}", truncate(rec, 80)));
        }
        let chrom_idx = ChromIds::instance().get_index(chrom);
        let pos: i32 = rec[tabs[0] + 1..tabs[1]]
            .parse()
            .map_err(|_| format!("invalid POS field: {}", truncate(rec, 80)))?;
        let id_field = &rec[tabs[1] + 1..tabs[2]];
        let id = if id_field == "." {
            None
        } else {
            Some(Arc::from(id_field))
        };
        let alleles_str = &rec[tabs[2] + 1..tabs[4]]; // "REF\tALT"
        let alt = &rec[tabs[3] + 1..tabs[4]];
        if alt.is_empty() {
            return Err(format!("missing ALT field: {}", truncate(rec, 80)));
        }
        let n_alleles: usize = if alt == "." {
            1
        } else {
            2 + alt.bytes().filter(|&b| b == b',').count()
        };
        if n_alleles > 255 {
            return Err(format!("more than 255 alleles: {}", truncate(rec, 80)));
        }
        let info = &rec[tabs[6] + 1..tabs[7]];
        let end = extract_end_value(info);
        Ok((
            Marker {
                chrom_idx,
                pos,
                id,
                alleles: Arc::from(alleles_str),
                n_alleles: n_alleles as u16,
                end,
            },
            tabs[7],
        ))
    }

    pub fn chrom(&self) -> Arc<str> {
        ChromIds::instance().id(self.chrom_idx)
    }

    pub fn id_str(&self) -> &str {
        match &self.id {
            Some(s) => s,
            None => ".",
        }
    }

    /// The INFO field as printed in imputed output: only the END subfield
    /// is retained (`MarkerParser` is constructed with `storeInfo=false`).
    pub fn end_subfield(&self) -> Option<String> {
        self.end.as_ref().map(|v| format!("END={}", v))
    }
}

/// First INFO/END subfield of the INFO field (split on ';').
fn extract_end_value(info: &str) -> Option<Arc<str>> {
    if info == "." {
        return None;
    }
    for field in info.split(';') {
        if let Some(v) = field.strip_prefix("END=") {
            return Some(Arc::from(v));
        }
    }
    None
}

pub fn truncate(s: &str, max: usize) -> &str {
    if s.len() <= max {
        s
    } else {
        let mut end = max;
        while !s.is_char_boundary(end) {
            end -= 1;
        }
        &s[..end]
    }
}
