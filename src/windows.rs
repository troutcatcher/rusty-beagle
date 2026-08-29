//! Port of `vcf.MarkerIndices`, `vcf.Window`, and
//! `vcf.RefTargSlidingWindow` (ref+target sliding windows).

use crate::genmap::GeneticMap;
use crate::marker::Marker;
use crate::par::Par;
use crate::refpanel::{RefReader, RefRec};
use crate::vcfio::{GtRec, LineSource, VcfHeader};
use std::collections::HashSet;
use std::sync::Arc;

/// Port of `vcf.MarkerIndices`.
pub struct MarkerIndices {
    pub prev_splice: usize,
    pub overlap_end: usize,
    pub overlap_start: usize,
    pub next_splice: usize,

    pub targ_marker_to_marker: Vec<usize>,
    pub marker_to_targ_marker: Vec<isize>,

    pub targ_prev_splice: usize,
    pub targ_overlap_end: usize,
    pub targ_overlap_start: usize,
    pub targ_next_splice: usize,
}

impl MarkerIndices {
    pub fn new(in_targ: &[bool], overlap_end: usize, overlap_start: usize) -> MarkerIndices {
        let n_markers = in_targ.len();
        let prev_splice = overlap_end >> 1;
        let next_splice = (n_markers + overlap_start) >> 1;
        let targ_marker_to_marker: Vec<usize> = in_targ
            .iter()
            .enumerate()
            .filter(|(_, &b)| b)
            .map(|(i, _)| i)
            .collect();
        let mut marker_to_targ_marker: Vec<isize> = vec![-1; n_markers];
        for (t, &m) in targ_marker_to_marker.iter().enumerate() {
            marker_to_targ_marker[m] = t as isize;
        }
        let targ_index = |marker: usize| -> usize {
            targ_marker_to_marker.partition_point(|&m| m < marker)
        };
        MarkerIndices {
            prev_splice,
            overlap_end,
            overlap_start,
            next_splice,
            targ_prev_splice: targ_index(prev_splice),
            targ_overlap_end: targ_index(overlap_end),
            targ_overlap_start: targ_index(overlap_start),
            targ_next_splice: targ_index(next_splice),
            targ_marker_to_marker,
            marker_to_targ_marker,
        }
    }

    pub fn n_markers(&self) -> usize {
        self.marker_to_targ_marker.len()
    }

    pub fn n_targ_markers(&self) -> usize {
        self.targ_marker_to_marker.len()
    }
}

/// One sliding window of reference + target records.
pub struct Window {
    pub window_index: usize,
    pub last_window: bool,
    pub indices: MarkerIndices,
    pub ref_recs: Vec<Arc<RefRec>>,
    pub targ_recs: Vec<Arc<GtRec>>,
}

impl Window {
    pub fn chrom_index(&self) -> u16 {
        self.targ_recs[0].marker.chrom_idx
    }

    /// true iff every target genotype in the window is phased & non-missing
    pub fn targ_is_phased(&self) -> bool {
        self.targ_recs.iter().all(|r| r.phased)
    }

    /// reference records restricted to target markers
    pub fn restrict_ref(&self) -> Vec<Arc<RefRec>> {
        self.indices
            .targ_marker_to_marker
            .iter()
            .map(|&m| self.ref_recs[m].clone())
            .collect()
    }
}

/// Streaming target-record reader with chrom-interval and marker filters.
pub struct TargReader {
    pub header: VcfHeader,
    lines: LineSource,
    exclude_markers: HashSet<String>,
    chrom_int: Option<crate::par::ChromInterval>,
}

impl TargReader {
    pub fn new(
        header: VcfHeader,
        lines: LineSource,
        exclude_markers: HashSet<String>,
        chrom_int: Option<crate::par::ChromInterval>,
    ) -> TargReader {
        TargReader {
            header,
            lines,
            exclude_markers,
            chrom_int,
        }
    }

    fn accept(&self, marker: &Marker) -> bool {
        if let Some(id) = &marker.id {
            if self.exclude_markers.contains(id.as_ref()) {
                return false;
            }
        }
        if let Some(ci) = &self.chrom_int {
            if !ci.contains(&marker.chrom(), marker.pos) {
                return false;
            }
        }
        true
    }

    pub fn next(&mut self) -> Option<Arc<GtRec>> {
        loop {
            let line = self.lines.next_line()?;
            let rec = crate::vcfio::parse_gt_rec(&self.header, &line).unwrap_or_else(|e| {
                eprintln!("ERROR: {}", e);
                std::process::exit(1)
            });
            if self.accept(&rec.marker) {
                return Some(Arc::new(rec));
            }
        }
    }
}

/// Reference reader wrapper applying the chrom-interval filter.
pub struct FilteredRefReader {
    inner: RefReader,
    chrom_int: Option<crate::par::ChromInterval>,
}

impl FilteredRefReader {
    pub fn new(inner: RefReader, chrom_int: Option<crate::par::ChromInterval>) -> Self {
        FilteredRefReader { inner, chrom_int }
    }

    pub fn n_haps(&self) -> usize {
        self.inner.n_haps()
    }

    pub fn samples(&self) -> crate::vcfio::Samples {
        self.inner.header.samples.clone()
    }

    pub fn next(&mut self) -> Option<Arc<RefRec>> {
        loop {
            let rec = self.inner.next()?;
            match &self.chrom_int {
                None => return Some(rec),
                Some(ci) => {
                    if ci.contains(&rec.marker.chrom(), rec.marker.pos) {
                        return Some(rec);
                    }
                }
            }
        }
    }
}

/// Port of `vcf.RefTargSlidingWindow.Reader`.
pub struct SlidingWindows {
    genmap: Arc<GeneticMap>,
    window_cm: f64,
    window_markers: usize,
    overlap_cm: f64,
    overlap_markers: usize,
    impute: bool,

    targ_it: TargReader,
    ref_it: FilteredRefReader,

    targ_overlap: Vec<Arc<GtRec>>,
    ref_overlap: Vec<Arc<RefRec>>,
    in_targ_overlap: Vec<bool>,
    targ_recs: Vec<Arc<GtRec>>,
    ref_recs: Vec<Arc<RefRec>>,
    in_targ: Vec<bool>,
    next_targ_rec: Option<Arc<GtRec>>,
    next_ref_rec: Option<Arc<RefRec>>,

    window_index: usize,
    started: bool,
    finished: bool,

    pub cum_targ_markers: usize,
    pub cum_ref_markers: usize,
}

impl SlidingWindows {
    pub fn new(
        par: &Par,
        genmap: Arc<GeneticMap>,
        targ_it: TargReader,
        ref_it: FilteredRefReader,
    ) -> SlidingWindows {
        SlidingWindows {
            genmap,
            window_cm: par.window as f64,
            window_markers: par.window_markers,
            overlap_cm: par.overlap as f64,
            overlap_markers: par.window_markers >> 2,
            impute: par.impute,
            targ_it,
            ref_it,
            targ_overlap: Vec::new(),
            ref_overlap: Vec::new(),
            in_targ_overlap: Vec::new(),
            targ_recs: Vec::new(),
            ref_recs: Vec::new(),
            in_targ: Vec::new(),
            next_targ_rec: None,
            next_ref_rec: None,
            window_index: 0,
            started: false,
            finished: false,
            cum_targ_markers: 0,
            cum_ref_markers: 0,
        }
    }

    pub fn targ_samples(&self) -> crate::vcfio::Samples {
        self.targ_it.header.samples.clone()
    }

    /// `Reader.run` loop body: produces the next window, or None.
    pub fn next_window(&mut self) -> Option<Window> {
        if self.finished {
            return None;
        }
        if !self.started {
            self.started = true;
            self.next_targ_rec = self.targ_it.next();
            self.next_ref_rec = self.ref_it.next();
            if self.next_targ_rec.is_none() || self.next_ref_rec.is_none() {
                eprintln!("ERROR: no genotype data");
                std::process::exit(1);
            }
        }
        if self.next_targ_rec.is_none() || self.next_ref_rec.is_none() {
            self.finished = true;
            return None;
        }
        let chrom_index = self.next_targ_rec.as_ref().unwrap().marker.chrom_idx;
        self.advance_ref_it_to_chrom(chrom_index);
        let next_ref_marker = match &self.next_ref_rec {
            Some(r) => r.marker.clone(),
            None => {
                self.finished = true;
                return None;
            }
        };
        let end_cm = self.next_end_cm(&next_ref_marker);
        let end_pos = self.genmap.base_pos(chrom_index, end_cm);
        self.window_index += 1;
        let window = self.read_window(chrom_index, end_pos, self.window_index);

        // save overlap for next window
        let targ_overlap_start = window.indices.targ_overlap_start;
        let ref_overlap_start = window.indices.overlap_start;
        self.targ_overlap
            .extend_from_slice(&self.targ_recs[targ_overlap_start..]);
        self.ref_overlap
            .extend_from_slice(&self.ref_recs[ref_overlap_start..]);
        self.in_targ_overlap
            .extend_from_slice(&self.in_targ[ref_overlap_start..]);

        if window.last_window {
            self.finished = true;
        }
        self.cum_targ_markers += window.indices.n_targ_markers() - window.indices.targ_overlap_end;
        self.cum_ref_markers += window.indices.n_markers() - window.indices.overlap_end;
        Some(window)
    }

    fn next_end_cm(&self, next_ref_marker: &Marker) -> f64 {
        let mut end_cm = self.genmap.gen_pos_marker(next_ref_marker);
        if self.ref_overlap.is_empty() {
            end_cm += self.window_cm;
        } else {
            end_cm += self.window_cm - self.overlap_cm;
        }
        end_cm
    }

    fn read_window(&mut self, chrom_index: u16, end_pos: i32, window_index: usize) -> Window {
        let ref_overlap_end = self.ref_overlap.len();
        self.reset_lists();
        while let Some(targ_rec) = self.next_targ_rec.clone() {
            if targ_rec.marker.chrom_idx != chrom_index
                || targ_rec.marker.pos >= end_pos
                || self.ref_recs.len() >= self.window_markers
            {
                break;
            }
            let targ_marker = &targ_rec.marker;
            let targ_pos = targ_marker.pos;
            while let Some(ref_rec) = self.next_ref_rec.clone() {
                if ref_rec.marker.chrom_idx == chrom_index
                    && (ref_rec.marker.pos < targ_pos
                        || (ref_rec.marker.pos == targ_pos && &ref_rec.marker != targ_marker))
                {
                    if self.impute {
                        self.ref_recs.push(ref_rec);
                        self.in_targ.push(false);
                    }
                    self.next_ref_rec = self.ref_it.next();
                } else {
                    break;
                }
            }
            if let Some(ref_rec) = self.next_ref_rec.clone() {
                if ref_rec.marker == *targ_marker {
                    self.targ_recs.push(targ_rec.clone());
                    self.ref_recs.push(ref_rec);
                    self.in_targ.push(true);
                    self.next_ref_rec = self.ref_it.next();
                }
            }
            self.next_targ_rec = self.targ_it.next();
        }
        if self.impute {
            while let Some(ref_rec) = self.next_ref_rec.clone() {
                if ref_rec.marker.chrom_idx == chrom_index
                    && ref_rec.marker.pos < end_pos
                    && self.ref_recs.len() < self.window_markers
                {
                    self.ref_recs.push(ref_rec);
                    self.in_targ.push(false);
                    self.next_ref_rec = self.ref_it.next();
                } else {
                    break;
                }
            }
        }
        self.build_window(ref_overlap_end, window_index, chrom_index, end_pos)
    }

    fn reset_lists(&mut self) {
        self.targ_recs.clear();
        self.ref_recs.clear();
        self.in_targ.clear();
        self.targ_recs.append(&mut self.targ_overlap);
        self.ref_recs.append(&mut self.ref_overlap);
        self.in_targ.append(&mut self.in_targ_overlap);
    }

    fn build_window(
        &mut self,
        ref_overlap_end: usize,
        window_index: usize,
        chrom_index: u16,
        end_pos: i32,
    ) -> Window {
        if self.targ_recs.is_empty() || self.ref_recs.is_empty() {
            if self.ref_recs.is_empty() {
                eprintln!(
                    "ERROR: The window ending at {}:{} contains no reference markers\n\
                     Do the reference and target VCF files contain the same\n\
                     chromosomes in the same order?",
                    crate::marker::ChromIds::instance().id(chrom_index),
                    end_pos
                );
            } else {
                eprintln!(
                    "ERROR: The reference and target VCF files contain no markers in common in \
                     the window: {}:{}-{}\n\
                     Do both VCF files share any markers in this window?\n\
                     Do both VCF files contain the same chromosomes in the same order?",
                    crate::marker::ChromIds::instance().id(chrom_index),
                    self.ref_recs[0].marker.pos,
                    end_pos
                );
            }
            std::process::exit(1);
        }
        let last_window = self.next_targ_rec.is_none() || self.next_ref_rec.is_none();
        let chrom_end = last_window
            || self.ref_recs[0].marker.chrom_idx
                != self.next_ref_rec.as_ref().unwrap().marker.chrom_idx;
        let ref_overlap_start = self.overlap_start(chrom_end, end_pos);
        let indices = MarkerIndices::new(&self.in_targ, ref_overlap_end, ref_overlap_start);
        Window {
            window_index,
            last_window,
            indices,
            ref_recs: self.ref_recs.clone(),
            targ_recs: self.targ_recs.clone(),
        }
    }

    fn overlap_start(&self, chrom_end: bool, end_pos: i32) -> usize {
        let n_markers = self.ref_recs.len();
        if chrom_end {
            return n_markers;
        }
        let n_markers_m1 = n_markers - 1;
        let chrom_index = self.ref_recs[n_markers_m1].marker.chrom_idx;
        let end_gen_pos = self.genmap.gen_pos(chrom_index, end_pos - 1);
        let start_gen_pos = end_gen_pos - self.overlap_cm;
        let key = self.genmap.base_pos(chrom_index, start_gen_pos);
        let mut low = n_markers.saturating_sub(self.overlap_markers) as isize;
        let mut high = n_markers_m1 as isize;
        while low <= high {
            let mid = (low + high) >> 1;
            let mid_pos = self.ref_recs[mid as usize].marker.pos;
            if mid_pos < key {
                low = mid + 1;
            } else if mid_pos > key {
                high = mid - 1;
            } else {
                return self.first_index_with_pos(mid as usize);
            }
        }
        debug_assert!(high < low);
        self.first_index_with_pos(high.max(0) as usize)
    }

    fn first_index_with_pos(&self, mut index: usize) -> usize {
        let pos = self.ref_recs[index].marker.pos;
        while index > 0 && self.ref_recs[index - 1].marker.pos == pos {
            index -= 1;
        }
        index
    }

    fn advance_ref_it_to_chrom(&mut self, chrom_index: u16) {
        while let Some(rec) = &self.next_ref_rec {
            if rec.marker.chrom_idx != chrom_index {
                self.next_ref_rec = self.ref_it.next();
            } else {
                break;
            }
        }
    }
}

/// Cumulative statistics reported when the background reader finishes.
pub struct WindowStats {
    pub cum_targ_markers: usize,
    pub cum_ref_markers: usize,
}

/// Runs the sliding-window reader on a background thread (like Java's
/// daemon `Reader` with an `ArrayBlockingQueue(1)`), so the next window is
/// decompressed and parsed while the current one is being imputed.
pub struct BackgroundWindows {
    rx: std::sync::mpsc::Receiver<Window>,
    stats_rx: std::sync::mpsc::Receiver<WindowStats>,
    handle: Option<std::thread::JoinHandle<()>>,
}

impl BackgroundWindows {
    pub fn spawn(mut sliding: SlidingWindows) -> BackgroundWindows {
        let (tx, rx) = std::sync::mpsc::sync_channel::<Window>(1);
        let (stats_tx, stats_rx) = std::sync::mpsc::sync_channel::<WindowStats>(1);
        let handle = std::thread::spawn(move || {
            while let Some(window) = sliding.next_window() {
                if tx.send(window).is_err() {
                    return; // consumer dropped
                }
            }
            let _ = stats_tx.send(WindowStats {
                cum_targ_markers: sliding.cum_targ_markers,
                cum_ref_markers: sliding.cum_ref_markers,
            });
        });
        BackgroundWindows {
            rx,
            stats_rx,
            handle: Some(handle),
        }
    }

    pub fn next_window(&mut self) -> Option<Window> {
        self.rx.recv().ok()
    }

    pub fn finish(mut self) -> WindowStats {
        let stats = self.stats_rx.recv().unwrap_or(WindowStats {
            cum_targ_markers: 0,
            cum_ref_markers: 0,
        });
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
        stats
    }
}
