//! Linear map geometry: `x` in `0.0..=1.0` along the backbone, lanes stacked above it.

use snappiesnap_doc::Molecule;

use crate::{Arrow, Tick, lanes, spans, ticks};

/// A feature segment as a bar. A segment that wraps past the origin of a circular molecule
/// drawn linearly becomes two bars.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeatureBar {
    pub feature: usize,
    pub x0: f32,
    pub x1: f32,
    pub lane: usize,
    pub arrow: Arrow,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LinearLabel {
    pub feature: usize,
    pub text: String,
    /// Centre of the label, in backbone units.
    pub x: f32,
    pub lane: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LinearLayout {
    pub len: usize,
    pub bars: Vec<FeatureBar>,
    pub ticks: Vec<(Tick, f32)>,
    pub labels: Vec<LinearLabel>,
    pub lanes: usize,
}

impl LinearLayout {
    pub fn x(&self, pos: f32) -> f32 {
        if self.len == 0 {
            0.0
        } else {
            pos / self.len as f32
        }
    }

    /// The position under `x` (backbone units), in `0..len`.
    pub fn position_at(&self, x: f32) -> usize {
        if self.len == 0 {
            return 0;
        }
        let pos = (x.clamp(0.0, 1.0) * self.len as f32).floor() as usize;
        pos.min(self.len.saturating_sub(1))
    }
}

pub fn linear(m: &Molecule, tick_target: usize) -> LinearLayout {
    let len = m.len();
    let spans = spans(m);
    let lane_of = lanes(&spans, len);
    let features = m.features();
    let x = |p: usize| if len == 0 { 0.0 } else { p as f32 / len as f32 };
    let mut bars = Vec::new();
    for (s, &lane) in spans.iter().zip(&lane_of) {
        let arrow = features
            .get(s.feature)
            .filter(|_| s.last)
            .map_or(Arrow::None, |f| Arrow::for_strand(f.strand));
        let (first, second) = s.range.segments(len);
        match second {
            None => bars.push(FeatureBar {
                feature: s.feature,
                x0: x(first.start),
                x1: x(first.end),
                lane,
                arrow,
            }),
            Some(second) => {
                // Arrowheads only where the feature really ends.
                let (a1, a2) = match arrow {
                    Arrow::End => (Arrow::None, Arrow::End),
                    Arrow::Start => (Arrow::Start, Arrow::None),
                    Arrow::None => (Arrow::None, Arrow::None),
                };
                bars.push(FeatureBar {
                    feature: s.feature,
                    x0: x(first.start),
                    x1: x(first.end),
                    lane,
                    arrow: a1,
                });
                bars.push(FeatureBar {
                    feature: s.feature,
                    x0: x(second.start),
                    x1: x(second.end),
                    lane,
                    arrow: a2,
                });
            }
        }
    }
    let mut labels = Vec::new();
    for (index, f) in features.iter().enumerate() {
        if let Some(bar) = bars
            .iter()
            .filter(|b| b.feature == index)
            .max_by(|a, b| (a.x1 - a.x0).total_cmp(&(b.x1 - b.x0)))
        {
            labels.push(LinearLabel {
                feature: index,
                text: f.name.clone(),
                x: (bar.x0 + bar.x1) / 2.0,
                lane: bar.lane,
            });
        }
    }
    let lanes = lane_of.iter().max().map_or(0, |m| m.saturating_add(1));
    let ticks = ticks(len, tick_target)
        .into_iter()
        .map(|t| (t, x(t.pos)))
        .collect();
    LinearLayout {
        len,
        bars,
        ticks,
        labels,
        lanes,
    }
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, NewFeature, Range, Segment, Seq, Strand, Topology};

    use super::*;

    #[test]
    fn bars_split_at_the_origin() {
        let seq = Seq::new(Alphabet::Dna, "A".repeat(100)).unwrap();
        let mut m = Molecule::new("p", seq, Topology::Circular).unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "wrap".into(),
            kind: "misc_feature".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(90, 10),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        let l = linear(&m, 10);
        assert_eq!(l.bars.len(), 2);
        assert_eq!((l.bars[0].x0, l.bars[0].x1), (0.9, 1.0));
        assert_eq!(l.bars[0].arrow, Arrow::None);
        assert_eq!((l.bars[1].x0, l.bars[1].x1), (0.0, 0.1));
        assert_eq!(l.bars[1].arrow, Arrow::End);
        assert_eq!(l.labels.len(), 1);
        assert_eq!(l.position_at(0.5), 50);
        assert_eq!(l.position_at(2.0), 99);
    }
}
