//! Port of `vcf.GeneticMap`, `vcf.PositionMap`, and `vcf.PlinkGenMap`.

use crate::marker::{ChromIds, Marker};
use crate::vcfio::open_text;
use std::io::BufRead;

pub enum GeneticMap {
    /// map=null: genetic position = 1e-6 * base position
    Position { scale: f64 },
    Plink(PlinkGenMap),
}

pub struct PlinkGenMap {
    /// indexed by chromosome index; empty vectors for unmapped chromosomes
    base_pos: Vec<Vec<i32>>,
    gen_pos: Vec<Vec<f64>>,
}

const MIN_END_CM_DIST: f64 = 5.0;

impl GeneticMap {
    pub fn new(map_file: &Option<String>, chrom_filter: Option<&str>) -> GeneticMap {
        match map_file {
            None => GeneticMap::Position { scale: 1e-6 },
            Some(path) => GeneticMap::Plink(PlinkGenMap::from_file(path, chrom_filter)),
        }
    }

    pub fn gen_pos_marker(&self, marker: &Marker) -> f64 {
        self.gen_pos(marker.chrom_idx, marker.pos)
    }

    pub fn gen_pos(&self, chrom: u16, base_position: i32) -> f64 {
        match self {
            GeneticMap::Position { scale } => scale * base_position as f64,
            GeneticMap::Plink(m) => m.gen_pos(chrom, base_position),
        }
    }

    pub fn base_pos(&self, chrom: u16, genetic_position: f64) -> i32 {
        match self {
            GeneticMap::Position { scale } => {
                let pos = (genetic_position / scale).round();
                if pos > i32::MAX as f64 {
                    eprintln!(
                        "ERROR: An estimated base position exceeds the maximum integer value\n\
                         Is the window parameter in cM units?"
                    );
                    std::process::exit(1);
                }
                pos as i32
            }
            GeneticMap::Plink(m) => m.base_pos(chrom, genetic_position),
        }
    }
}

impl PlinkGenMap {
    fn from_file(path: &str, chrom_filter: Option<&str>) -> PlinkGenMap {
        let mut base_pos: Vec<Vec<i32>> = Vec::new();
        let mut gen_pos: Vec<Vec<f64>> = Vec::new();
        let mut reader = open_text(path);
        let mut line = String::new();
        loop {
            line.clear();
            let n = reader.read_line(&mut line).unwrap_or(0);
            if n == 0 {
                break;
            }
            let fields: Vec<&str> = line.split_whitespace().collect();
            if fields.is_empty() {
                continue;
            }
            if fields.len() < 4 {
                eprintln!("ERROR: Map file format error: {}", line.trim_end());
                std::process::exit(1);
            }
            let chrom = fields[0];
            if let Some(cf) = chrom_filter {
                if chrom != cf {
                    continue;
                }
            }
            let chrom_index = ChromIds::instance().get_index(chrom) as usize;
            while chrom_index >= base_pos.len() {
                base_pos.push(Vec::new());
                gen_pos.push(Vec::new());
            }
            let bp: i32 = fields[3].parse().unwrap_or_else(|_| {
                eprintln!("ERROR: invalid base position in map file: {}", line.trim_end());
                std::process::exit(1)
            });
            let gp: f64 = fields[2].parse().unwrap_or_else(|_| {
                eprintln!("ERROR: invalid cM position in map file: {}", line.trim_end());
                std::process::exit(1)
            });
            if !gp.is_finite() {
                eprintln!("ERROR: invalid map position: {}", line.trim_end());
                std::process::exit(1);
            }
            let bv = &mut base_pos[chrom_index];
            let gv = &mut gen_pos[chrom_index];
            if let Some(&last_bp) = bv.last() {
                if bp == last_bp {
                    eprintln!(
                        "ERROR: duplicate base position in genetic map: {}",
                        line.trim_end()
                    );
                    std::process::exit(1);
                }
                if bp < last_bp || gp < *gv.last().unwrap() {
                    eprintln!(
                        "ERROR: genetic map positions are not sorted in ascending order: {}",
                        line.trim_end()
                    );
                    std::process::exit(1);
                }
            }
            bv.push(bp);
            gv.push(gp);
        }
        for (c, gv) in gen_pos.iter().enumerate() {
            if gv.len() >= 2 && gv.first() == gv.last() {
                eprintln!(
                    "ERROR: all base positions on chromosome {} have the same genetic map position",
                    ChromIds::instance().id(c as u16)
                );
                std::process::exit(1);
            }
        }
        PlinkGenMap { base_pos, gen_pos }
    }

    fn check_chrom(&self, chrom: u16) {
        let c = chrom as usize;
        if c >= self.base_pos.len() || self.base_pos[c].is_empty() {
            eprintln!(
                "ERROR: missing genetic map for chromosome {}",
                ChromIds::instance().id(chrom)
            );
            std::process::exit(1);
        }
    }

    /// Java `Arrays.binarySearch` on a sorted slice: returns index if found,
    /// otherwise `-(insertionPoint) - 1`.
    fn binary_search_i32(a: &[i32], key: i32) -> isize {
        match a.binary_search(&key) {
            Ok(i) => i as isize,
            Err(i) => -(i as isize) - 1,
        }
    }

    fn binary_search_f64(a: &[f64], key: f64) -> isize {
        // Java compares doubles via Double.compare; positions are finite here.
        let mut low: isize = 0;
        let mut high: isize = a.len() as isize - 1;
        while low <= high {
            let mid = (low + high) >> 1;
            let v = a[mid as usize];
            if v < key {
                low = mid + 1;
            } else if v > key {
                high = mid - 1;
            } else {
                return mid;
            }
        }
        -(low) - 1
    }

    pub fn gen_pos(&self, chrom: u16, base_position: i32) -> f64 {
        self.check_chrom(chrom);
        let c = chrom as usize;
        let bp = &self.base_pos[c];
        let gp = &self.gen_pos[c];
        let map_size_m1 = bp.len() - 1;
        assert!(map_size_m1 > 0, "genetic map for chromosome has < 2 positions");
        let index = Self::binary_search_i32(bp, base_position);
        if index >= 0 {
            return gp[index as usize];
        }
        let ins_pt = (-index - 1) as usize;
        let mut a_index = ins_pt as isize - 1;
        let mut b_index = ins_pt;
        if a_index == map_size_m1 as isize {
            let mut ins = Self::binary_search_f64(gp, gp[map_size_m1] - MIN_END_CM_DIST);
            if ins < 0 {
                ins = -ins - 2;
            }
            a_index = ins.max(0);
            b_index = map_size_m1;
        } else if b_index == 0 {
            let mut ins = Self::binary_search_f64(gp, gp[0] + MIN_END_CM_DIST);
            if ins < 0 {
                ins = -ins - 1;
            }
            a_index = 0;
            b_index = (ins as usize).min(map_size_m1);
        }
        let x = base_position;
        let a = bp[a_index as usize];
        let b = bp[b_index];
        let fa = gp[a_index as usize];
        let fb = gp[b_index];
        fa + (((x - a) as f64) / ((b - a) as f64)) * (fb - fa)
    }

    pub fn base_pos(&self, chrom: u16, genetic_position: f64) -> i32 {
        self.check_chrom(chrom);
        let c = chrom as usize;
        let bp = &self.base_pos[c];
        let gp = &self.gen_pos[c];
        let map_size_m1 = gp.len() - 1;
        assert!(map_size_m1 > 0, "genetic map for chromosome has < 2 positions");
        let index = Self::binary_search_f64(gp, genetic_position);
        if index >= 0 {
            return bp[index as usize];
        }
        let ins_pt = (-index - 1) as usize;
        let mut a_index = ins_pt as isize - 1;
        let mut b_index = ins_pt;
        if a_index == map_size_m1 as isize {
            let mut ins = Self::binary_search_f64(gp, gp[map_size_m1] - MIN_END_CM_DIST);
            if ins < 0 {
                ins = -ins - 2;
            }
            a_index = ins.max(0);
            b_index = map_size_m1;
        } else if b_index == 0 {
            let mut ins = Self::binary_search_f64(gp, gp[0] + MIN_END_CM_DIST);
            if ins < 0 {
                ins = -ins - 1;
            }
            a_index = 0;
            b_index = (ins as usize).min(map_size_m1);
        }
        let x = genetic_position;
        let a = gp[a_index as usize];
        let b = gp[b_index];
        let fa = bp[a_index as usize];
        let fb = bp[b_index];
        let interp = fa as f64 + ((x - a) / (b - a)) * ((fb - fa) as f64);
        if interp >= i32::MAX as f64 {
            eprintln!(
                "ERROR: An estimated base position exceeds the maximum integer value\n\
                 Are the window parameter and the genetic map in cM units?"
            );
            std::process::exit(1);
        }
        interp.round() as i32
    }
}
