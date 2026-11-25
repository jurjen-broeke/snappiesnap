//! Sequence view geometry: which residues go on each line, and the feature bars under them.
//! Lines are computed on demand so only visible lines cost anything, whatever the sequence
//! length.

use snappiesnap_doc::{Molecule, Range};

use crate::{Arrow, overlaps};

/// A feature piece under one line of the sequence view.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrackBar {
    pub feature: usize,
    /// Absolute 0-based positions, within the line.
    pub start: usize,
    pub end: usize,
    pub lane: usize,
    /// Arrowhead, only on the piece where the feature really ends.
    pub arrow: Arrow,
    /// The first piece of the feature on this line (where a label goes).
    pub label: bool,
}

/// Number of lines for `len` residues at `per_line` residues per line.
pub fn line_count(len: usize, per_line: usize) -> usize {
    len.div_ceil(per_line.max(1))
}

/// The residue range of line `index`.
pub fn line_range(len: usize, per_line: usize, index: usize) -> (usize, usize) {
    let per_line = per_line.max(1);
    let start = index.saturating_mul(per_line).min(len);
    (start, start.saturating_add(per_line).min(len))
}

/// Feature pieces under the line covering `start..end`, in lanes that don't overlap.
pub fn line_tracks(m: &Molecule, start: usize, end: usize) -> Vec<TrackBar> {
    let len = m.len();
    if start >= end {
        return Vec::new();
    }
    let line = Range::new(start, end);
    let mut bars: Vec<TrackBar> = Vec::new();
    for (feature, f) in m.features().iter().enumerate() {
        let n = f.segments.len();
        for (i, seg) in f.segments.iter().enumerate() {
            if seg.range.is_empty() || !overlaps(seg.range, line, len) {
                continue;
            }
            let (a, b) = seg.range.segments(len);
            for piece in std::iter::once(a).chain(b) {
                let (s, e) = (piece.start.max(start), piece.end.min(end));
                if s >= e {
                    continue;
                }
                let last = i + 1 == n;
                let arrow = match Arrow::for_strand(f.strand) {
                    Arrow::End if last && e == seg_end(seg.range, len) => Arrow::End,
                    Arrow::Start if last && s == seg.range.start => Arrow::Start,
                    _ => Arrow::None,
                };
                let label = !bars.iter().any(|b| b.feature == feature);
                bars.push(TrackBar {
                    feature,
                    start: s,
                    end: e,
                    lane: 0,
                    arrow,
                    label,
                });
            }
        }
    }
    // Lanes: first fit, in order of start.
    let mut order: Vec<usize> = (0..bars.len()).collect();
    order.sort_by_key(|&i| bars.get(i).map(|b| (b.start, std::cmp::Reverse(b.end))));
    let mut lane_ends: Vec<usize> = Vec::new();
    for i in order {
        let Some(bar) = bars.get_mut(i) else { continue };
        let lane = lane_ends
            .iter()
            .position(|&e| e <= bar.start)
            .unwrap_or(lane_ends.len());
        if lane == lane_ends.len() {
            lane_ends.push(bar.end);
        } else if let Some(e) = lane_ends.get_mut(lane) {
            *e = bar.end;
        }
        bar.lane = lane;
    }
    bars
}

/// Where a range ends as a position in `0..=len` (a wrapping range ends at `end`).
fn seg_end(r: Range, len: usize) -> usize {
    if r.is_wrapping() {
        r.end
    } else {
        r.end.min(len)
    }
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, NewFeature, Segment, Seq, Strand, Topology};

    use super::*;

    fn molecule() -> Molecule {
        let seq = Seq::new(Alphabet::Dna, "A".repeat(250)).unwrap();
        let mut m = Molecule::new("p", seq, Topology::Circular).unwrap();
        for (name, start, end, strand) in [
            ("long", 50, 180, Strand::Forward),
            ("inner", 60, 70, Strand::Reverse),
            ("wrap", 240, 20, Strand::Forward),
        ] {
            m.add_feature(NewFeature {
                color: None,
                name: name.into(),
                kind: "misc_feature".into(),
                strand,
                segments: vec![Segment {
                    range: Range::new(start, end),
                }],
                qualifiers: vec![],
            })
            .unwrap();
        }
        m
    }

    #[test]
    fn lines() {
        assert_eq!(line_count(250, 100), 3);
        assert_eq!(line_count(0, 100), 0);
        assert_eq!(line_count(5, 0), 5);
        assert_eq!(line_range(250, 100, 2), (200, 250));
        assert_eq!(line_range(250, 100, 9), (250, 250));
    }

    #[test]
    fn tracks_split_lanes_and_arrows() {
        let m = molecule();
        let first = line_tracks(&m, 0, 100);
        // "wrap" ends on this line (0..20), "long" starts here, "inner" sits inside it.
        let wrap = first.iter().find(|b| b.feature == 2).unwrap();
        assert_eq!((wrap.start, wrap.end, wrap.arrow), (0, 20, Arrow::End));
        let long = first.iter().find(|b| b.feature == 0).unwrap();
        assert_eq!((long.start, long.end, long.arrow), (50, 100, Arrow::None));
        let inner = first.iter().find(|b| b.feature == 1).unwrap();
        assert_eq!(inner.arrow, Arrow::Start);
        assert_ne!(inner.lane, long.lane);
        let second = line_tracks(&m, 100, 200);
        let long2 = second.iter().find(|b| b.feature == 0).unwrap();
        assert_eq!((long2.end, long2.arrow), (180, Arrow::End));
        let last = line_tracks(&m, 200, 250);
        let wrap3 = last.iter().find(|b| b.feature == 2).unwrap();
        assert_eq!(
            (wrap3.start, wrap3.end, wrap3.arrow),
            (240, 250, Arrow::None)
        );
        assert!(line_tracks(&m, 10, 10).is_empty());
    }
}
