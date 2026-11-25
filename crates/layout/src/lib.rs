//! SnappieSnap `layout` (layer L2).
//!
//! Geometry for map views: where features, ticks and labels go. Pure functions of the document,
//! independent of any drawing library, so they are unit-tested here and shared by the egui view
//! and the offscreen renderer.
//!
//! Circular maps use unit coordinates: the backbone is the circle of radius 1 around the origin,
//! angles are in radians, clockwise from the top (position 0 at 12 o'clock). Linear maps use
//! `x` in `0.0..=1.0` along the backbone.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

mod circular;
mod linear;
mod sequence;
mod ticks;

pub use circular::{CircularLayout, CircularOptions, FeatureArc, LabelTarget, MapLabel, circular};
pub use linear::{FeatureBar, LinearLabel, LinearLayout, linear};
pub use sequence::{TrackBar, line_count, line_range, line_tracks};
pub use ticks::{Tick, tick_step, ticks};

use snappiesnap_doc::{Molecule, Range, Strand};

/// A position to label on a map, such as a restriction site's cut.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteMark {
    /// Residues before the mark (a cut position, `0..=len`).
    pub pos: usize,
    pub text: String,
    /// Keep this label when the map runs out of room: the caller has singled it out (the
    /// enzymes named for a cloning step), so dropping it would hide the answer.
    pub keep: bool,
}

/// Which way an arrowhead points along the sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Arrow {
    None,
    /// At the end of the span (forward strand).
    End,
    /// At the start of the span (reverse strand).
    Start,
}

impl Arrow {
    pub fn for_strand(strand: Strand) -> Self {
        match strand {
            Strand::Forward => Arrow::End,
            Strand::Reverse => Arrow::Start,
            Strand::None | Strand::Both => Arrow::None,
        }
    }
}

/// One drawable piece of a feature: a segment with its lane.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Span {
    pub feature: usize,
    pub range: Range,
    /// The feature's last segment (5'→3'), the one that carries the arrowhead.
    pub last: bool,
}

/// Every feature segment, in feature order.
pub(crate) fn spans(m: &Molecule) -> Vec<Span> {
    m.features()
        .iter()
        .enumerate()
        .flat_map(|(feature, f)| {
            let n = f.segments.len();
            f.segments.iter().enumerate().map(move |(i, s)| Span {
                feature,
                range: s.range,
                last: i + 1 == n,
            })
        })
        .filter(|s| !s.range.is_empty())
        .collect()
}

/// Assign each span the lowest lane where it overlaps no earlier span (lane 0 first). Longer
/// spans are placed first so big features stay on the backbone.
pub(crate) fn lanes(spans: &[Span], len: usize) -> Vec<usize> {
    let mut order: Vec<usize> = (0..spans.len()).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(spans.get(i).map_or(0, |s| s.range.len_on(len))));
    let mut lane_of = vec![0usize; spans.len()];
    let mut placed: Vec<Vec<Range>> = Vec::new();
    for i in order {
        let Some(span) = spans.get(i) else { continue };
        let lane = placed
            .iter()
            .position(|ranges| !ranges.iter().any(|r| overlaps(*r, span.range, len)))
            .unwrap_or(placed.len());
        if lane == placed.len() {
            placed.push(Vec::new());
        }
        if let Some(ranges) = placed.get_mut(lane) {
            ranges.push(span.range);
        }
        if let Some(slot) = lane_of.get_mut(i) {
            *slot = lane;
        }
    }
    lane_of
}

/// The most feature lanes any view of this molecule can need: an upper bound for every line of
/// the sequence view, so it can use one row height.
pub fn max_lanes(m: &Molecule) -> usize {
    let spans = spans(m);
    lanes(&spans, m.len())
        .into_iter()
        .max()
        .map_or(0, |l| l.saturating_add(1))
}

/// Whether two non-empty ranges share a position (wrap-aware).
pub(crate) fn overlaps(a: Range, b: Range, len: usize) -> bool {
    a.contains(b.start, len) || b.contains(a.start, len)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlap_and_lanes() {
        let r = Range::new;
        assert!(overlaps(r(0, 10), r(5, 15), 100));
        assert!(!overlaps(r(0, 10), r(10, 15), 100));
        assert!(overlaps(r(90, 5), r(2, 3), 100));
        assert!(!overlaps(r(90, 5), r(5, 89), 100));
        let spans = [
            Span {
                feature: 0,
                range: r(0, 50),
                last: true,
            },
            Span {
                feature: 1,
                range: r(10, 20),
                last: true,
            },
            Span {
                feature: 2,
                range: r(60, 70),
                last: true,
            },
            Span {
                feature: 3,
                range: r(95, 15),
                last: true,
            },
        ];
        let l = lanes(&spans, 100);
        assert_eq!(l[0], 0);
        assert_eq!(l[2], 0);
        assert_ne!(l[1], l[0]);
        assert_ne!(l[3], l[0]);
        assert_ne!(l[3], l[1]);
    }
}
