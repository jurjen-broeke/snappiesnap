//! Circular map geometry in unit coordinates (backbone radius 1, clockwise from the top).

use std::f32::consts::{FRAC_PI_2, TAU};

use snappiesnap_doc::Molecule;

use crate::{Arrow, SiteMark, Tick, lanes, spans, ticks};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CircularOptions {
    /// Radius of lane 0 (on the backbone).
    pub lane0: f32,
    /// Radial distance between lanes (lanes go inwards).
    pub lane_step: f32,
    /// Radius at which labels sit.
    pub label_radius: f32,
    /// Minimum vertical gap between labels, in unit coordinates.
    pub label_gap: f32,
    /// How far above and below the centre label text may sit, in unit coordinates (the
    /// canvas's half-height, less a margin, over the radius). Labels that don't fit are left
    /// out, site labels before feature labels.
    pub label_limit: f32,
    /// Roughly how many ticks to draw.
    pub ticks: usize,
    /// Whether features are laid out (with their labels); off leaves only the site marks.
    pub features: bool,
}

impl Default for CircularOptions {
    fn default() -> Self {
        CircularOptions {
            lane0: 0.92,
            lane_step: 0.09,
            label_radius: 1.3,
            label_gap: 0.075,
            label_limit: 1.45,
            ticks: 10,
            features: true,
        }
    }
}

/// A feature segment as an arc.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FeatureArc {
    pub feature: usize,
    /// Angle of the segment start.
    pub start: f32,
    /// Clockwise sweep in radians (0 < sweep ≤ 2π).
    pub sweep: f32,
    pub lane: usize,
    pub radius: f32,
    pub arrow: Arrow,
}

/// What a map label names.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LabelTarget {
    /// A feature, by index.
    Feature(usize),
    /// A site mark, by index into the marks passed to the layout.
    Site(usize),
}

/// A label: where the leader line starts (on the feature or the backbone), and where the
/// text goes.
#[derive(Debug, Clone, PartialEq)]
pub struct MapLabel {
    pub target: LabelTarget,
    pub text: String,
    /// Dropped last when there is no room (see [`crate::SiteMark::keep`]).
    pub keep: bool,
    /// Point on the feature (unit coordinates).
    pub anchor: (f32, f32),
    /// Where the text is anchored (unit coordinates).
    pub at: (f32, f32),
    /// True if the text extends to the right of `at`, false if to the left.
    pub right: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct CircularLayout {
    pub len: usize,
    pub arcs: Vec<FeatureArc>,
    pub ticks: Vec<(Tick, f32)>,
    pub labels: Vec<MapLabel>,
    pub lanes: usize,
    /// Labels left out because they didn't fit the height (see `label_limit`).
    pub hidden_labels: usize,
}

impl CircularLayout {
    /// The angle of a position (clockwise from the top).
    pub fn angle(&self, pos: f32) -> f32 {
        angle(pos, self.len)
    }

    /// The position under an angle, in `0..len`.
    pub fn position_at(&self, angle: f32) -> usize {
        if self.len == 0 {
            return 0;
        }
        let turn = ((angle + FRAC_PI_2).rem_euclid(TAU)) / TAU;
        let pos = (turn * self.len as f32).floor() as usize;
        pos.min(self.len.saturating_sub(1))
    }
}

/// Angle (radians, screen coordinates: x right, y down) of a position, 0 at the top.
pub fn angle(pos: f32, len: usize) -> f32 {
    if len == 0 {
        return -FRAC_PI_2;
    }
    pos / len as f32 * TAU - FRAC_PI_2
}

pub(crate) fn point(angle: f32, radius: f32) -> (f32, f32) {
    (radius * angle.cos(), radius * angle.sin())
}

/// Lay out `m` with labels for its features and for `sites` (restriction sites and other
/// position marks), all spread together so none overlap.
pub fn circular(m: &Molecule, opts: &CircularOptions, sites: &[SiteMark]) -> CircularLayout {
    let len = m.len();
    let spans = if opts.features { spans(m) } else { Vec::new() };
    let lane_of = lanes(&spans, len);
    let features = m.features();
    let arcs: Vec<FeatureArc> = spans
        .iter()
        .zip(&lane_of)
        .map(|(s, &lane)| {
            let span = s.range.len_on(len) as f32;
            FeatureArc {
                feature: s.feature,
                start: angle(s.range.start as f32, len),
                sweep: if len == 0 {
                    0.0
                } else {
                    span / len as f32 * TAU
                },
                lane,
                radius: opts.lane0 - opts.lane_step * lane as f32,
                arrow: features
                    .get(s.feature)
                    .filter(|_| s.last)
                    .map_or(Arrow::None, |f| Arrow::for_strand(f.strand)),
            }
        })
        .collect();
    let lanes = lane_of.iter().max().map_or(0, |m| m.saturating_add(1));

    // One label per feature, at the middle of its longest arc.
    let mut labels: Vec<MapLabel> = Vec::new();
    for (index, f) in features.iter().enumerate().filter(|_| opts.features) {
        let Some(arc) = arcs
            .iter()
            .filter(|a| a.feature == index)
            .max_by(|a, b| a.sweep.total_cmp(&b.sweep))
        else {
            continue;
        };
        let mid = arc.start + arc.sweep / 2.0;
        let (x, y) = point(mid, opts.label_radius);
        labels.push(MapLabel {
            target: LabelTarget::Feature(index),
            text: f.name.clone(),
            anchor: point(mid, arc.radius),
            at: (x, y),
            right: x >= 0.0,
            keep: false,
        });
    }
    for (index, site) in sites.iter().enumerate() {
        let a = angle(site.pos as f32, len);
        let (x, y) = point(a, opts.label_radius);
        labels.push(MapLabel {
            target: LabelTarget::Site(index),
            text: site.text.clone(),
            anchor: point(a, 1.0),
            at: (x, y),
            right: x >= 0.0,
            keep: site.keep,
        });
    }
    let hidden_labels = spread_labels(&mut labels, opts.label_gap, opts.label_limit);

    let ticks = ticks(len, opts.ticks)
        .into_iter()
        .map(|t| (t, angle(t.pos as f32, len)))
        .collect();
    CircularLayout {
        len,
        arcs,
        ticks,
        labels,
        lanes,
        hidden_labels,
    }
}

/// Move labels vertically (per side) so consecutive labels are at least `gap` apart and all
/// lie within `-limit..=limit`, keeping each as close to its natural height as possible. Text
/// then sits just outside the circle at the new height. When a side holds more labels than
/// fit, the ones that would end up farthest from their natural height are left out, site
/// labels before feature labels. Returns how many were left out.
fn spread_labels(labels: &mut Vec<MapLabel>, gap: f32, limit: f32) -> usize {
    let gap = if gap.is_finite() && gap > 0.0 {
        gap
    } else {
        0.075
    };
    let limit = if limit.is_finite() {
        limit.max(gap)
    } else {
        1.45
    };
    // Labels that fit in one column: n labels need (n - 1) gaps.
    let fit = ((2.0 * limit / gap).floor() as usize).saturating_add(1);
    let mut drop: Vec<bool> = vec![false; labels.len()];
    for right in [true, false] {
        let mut side: Vec<usize> = (0..labels.len())
            .filter(|&i| labels.get(i).is_some_and(|l| l.right == right))
            .collect();
        let excess = side.len().saturating_sub(fit);
        if excess == 0 {
            continue;
        }
        // How far each label would move in a column without a top limit: the crowded
        // middle of a cluster moves least, so its labels stay.
        let natural: Vec<f32> = side
            .iter()
            .map(|&i| labels.get(i).map_or(0.0, |l| l.at.1))
            .collect();
        let mut ys = natural.clone();
        let order = top_down(&natural);
        place(&mut ys, &order, gap, f32::NEG_INFINITY, limit);
        let moved = |k: usize| match (ys.get(k), natural.get(k)) {
            (Some(y), Some(n)) => (y - n).abs(),
            _ => 0.0,
        };
        // Site labels go first, then the ones that moved most (earlier label on ties).
        let mut worst: Vec<usize> = (0..side.len()).collect();
        let is_site = |k: usize| {
            side.get(k)
                .and_then(|&i| labels.get(i))
                .is_some_and(|l| matches!(l.target, LabelTarget::Site(_)))
        };
        // A label the caller singled out goes last of all: it is the one being looked for.
        let keep = |k: usize| {
            side.get(k)
                .and_then(|&i| labels.get(i))
                .is_some_and(|l| l.keep)
        };
        worst.sort_by(|&a, &b| {
            keep(a)
                .cmp(&keep(b))
                .then(is_site(b).cmp(&is_site(a)))
                .then(moved(b).total_cmp(&moved(a)))
                .then(a.cmp(&b))
        });
        for k in worst.into_iter().take(excess) {
            if let Some(&i) = side.get(k)
                && let Some(d) = drop.get_mut(i)
            {
                *d = true;
            }
        }
        side.retain(|&i| !drop.get(i).copied().unwrap_or(false));
    }
    let mut keep = drop.iter().map(|d| !d);
    labels.retain(|_| keep.next().unwrap_or(true));
    let hidden = drop.iter().filter(|d| **d).count();
    for right in [true, false] {
        let side: Vec<usize> = (0..labels.len())
            .filter(|&i| labels.get(i).is_some_and(|l| l.right == right))
            .collect();
        let mut ys: Vec<f32> = side
            .iter()
            .map(|&i| labels.get(i).map_or(0.0, |l| l.at.1))
            .collect();
        let order = top_down(&ys);
        place(&mut ys, &order, gap, -limit, limit);
        for (k, &i) in side.iter().enumerate() {
            let (Some(l), Some(&y)) = (labels.get_mut(i), ys.get(k)) else {
                continue;
            };
            l.at.1 = y;
            // Horizontal position: just outside the circle at that height (outside the tick
            // labels, which sit at about 1.1).
            let yc = y.clamp(-1.35, 1.35);
            let r: f32 = 1.24;
            let x = (r * r - yc * yc).max(0.0).sqrt().max(0.25);
            l.at.0 = if right { x } else { -x };
        }
    }
    hidden
}

/// Indices of `ys`, top (smallest) first.
fn top_down(ys: &[f32]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..ys.len()).collect();
    order.sort_by(|&a, &b| match (ys.get(a), ys.get(b)) {
        (Some(x), Some(y)) => x.total_cmp(y),
        _ => std::cmp::Ordering::Equal,
    });
    order
}

/// Space the heights `ys`, taken in `order` (top to bottom), at least `gap` apart within
/// `top..=bottom`: push down from the top, then pull back up from the bottom. If they don't
/// all fit, the top ones end up above `top`.
fn place(ys: &mut [f32], order: &[usize], gap: f32, top: f32, bottom: f32) {
    let mut prev: Option<f32> = None;
    for &k in order {
        let Some(y) = ys.get_mut(k) else { continue };
        let min = prev.map_or(top, |p| p + gap);
        if *y < min {
            *y = min;
        }
        prev = Some(*y);
    }
    let mut next: Option<f32> = None;
    for &k in order.iter().rev() {
        let Some(y) = ys.get_mut(k) else { continue };
        let max = next.map_or(bottom, |n| n - gap);
        if *y > max {
            *y = max;
        }
        next = Some(*y);
    }
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, NewFeature, Range, Segment, Seq, Strand, Topology};

    use super::*;

    fn plasmid(features: &[(&str, usize, usize, Strand)]) -> Molecule {
        let seq = Seq::new(Alphabet::Dna, "A".repeat(1000)).unwrap();
        let mut m = Molecule::new("p", seq, Topology::Circular).unwrap();
        for (name, start, end, strand) in features {
            m.add_feature(NewFeature {
                color: None,
                name: name.to_string(),
                kind: "misc_feature".into(),
                strand: *strand,
                segments: vec![Segment {
                    range: Range::new(*start, *end),
                }],
                qualifiers: vec![],
            })
            .unwrap();
        }
        m
    }

    #[test]
    fn angles_and_positions() {
        let m = plasmid(&[]);
        let l = circular(&m, &CircularOptions::default(), &[]);
        assert!((l.angle(0.0) + FRAC_PI_2).abs() < 1e-6);
        assert!((l.angle(250.0) - 0.0).abs() < 1e-5);
        assert_eq!(l.position_at(l.angle(0.0)), 0);
        assert_eq!(l.position_at(l.angle(250.5)), 250);
        assert_eq!(l.position_at(l.angle(999.5)), 999);
        assert_eq!(l.ticks.len(), 9);
    }

    #[test]
    fn arcs_lanes_arrows_and_labels() {
        let m = plasmid(&[
            ("big", 100, 600, Strand::Forward),
            ("inside", 200, 300, Strand::Reverse),
            ("wrap", 900, 50, Strand::None),
        ]);
        let l = circular(&m, &CircularOptions::default(), &[]);
        assert_eq!(l.arcs.len(), 3);
        let big = l.arcs[0];
        assert!((big.sweep - TAU / 2.0).abs() < 1e-5);
        assert_eq!(big.lane, 0);
        assert_eq!(big.arrow, Arrow::End);
        assert_eq!(l.arcs[1].lane, 1);
        assert_eq!(l.arcs[1].arrow, Arrow::Start);
        let wrap = l.arcs[2];
        assert_eq!(wrap.lane, 0);
        assert!((wrap.sweep - 150.0 / 1000.0 * TAU).abs() < 1e-5);
        assert_eq!(l.lanes, 2);
        assert_eq!(l.labels.len(), 3);
    }

    #[test]
    fn crowded_labels_are_spread() {
        let features: Vec<(String, usize, usize, Strand)> = (0..12)
            .map(|i| (format!("f{i}"), 100 + i * 5, 104 + i * 5, Strand::None))
            .collect();
        let refs: Vec<(&str, usize, usize, Strand)> = features
            .iter()
            .map(|(n, a, b, s)| (n.as_str(), *a, *b, *s))
            .collect();
        let opts = CircularOptions::default();
        let l = circular(&plasmid(&refs), &opts, &[]);
        let mut ys: Vec<f32> = l
            .labels
            .iter()
            .filter(|x| x.right)
            .map(|x| x.at.1)
            .collect();
        ys.sort_by(f32::total_cmp);
        for w in ys.windows(2) {
            assert!(w[1] - w[0] >= opts.label_gap - 1e-5, "{ys:?}");
        }
    }

    /// More labels than the height holds: the column stays within the limits, site labels are
    /// left out before feature labels, and how many is reported.
    #[test]
    fn labels_never_leave_the_canvas() {
        let m = plasmid(&[
            ("marker", 950, 990, Strand::None),
            ("gene", 600, 700, Strand::None),
        ]);
        // Many sites near the top-left (just before the origin), as in a crowded map.
        let sites: Vec<SiteMark> = (0..60)
            .map(|i| SiteMark {
                pos: 700 + i * 5,
                text: format!("E{i}I"),
                keep: false,
            })
            .collect();
        for limit in [1.45, 1.8, 3.0] {
            let opts = CircularOptions {
                label_limit: limit,
                ..CircularOptions::default()
            };
            let l = circular(&m, &opts, &sites);
            for side in [true, false] {
                let mut ys: Vec<f32> = l
                    .labels
                    .iter()
                    .filter(|x| x.right == side)
                    .map(|x| x.at.1)
                    .collect();
                ys.sort_by(f32::total_cmp);
                for y in &ys {
                    assert!((-limit - 1e-4..=limit + 1e-4).contains(y), "{limit} {ys:?}");
                }
                for w in ys.windows(2) {
                    assert!(w[1] - w[0] >= opts.label_gap - 1e-4, "{ys:?}");
                }
            }
            assert_eq!(l.labels.len() + l.hidden_labels, 62);
            for name in ["marker", "gene"] {
                assert!(
                    l.labels.iter().any(|x| x.text == name),
                    "{name} kept at {limit}"
                );
            }
        }
        let tight = CircularOptions {
            label_limit: 1.45,
            ..CircularOptions::default()
        };
        assert!(circular(&m, &tight, &sites).hidden_labels > 0);
        // A site the caller singled out survives a crowd that drops most of the others.
        let mut kept = sites.clone();
        if let Some(one) = kept.get_mut(7) {
            one.keep = true;
            one.text = "KeepMeI".to_string();
        }
        let l = circular(&m, &tight, &kept);
        assert!(l.hidden_labels > 0, "the crowd is still too big");
        assert!(
            l.labels.iter().any(|x| x.text == "KeepMeI"),
            "a kept label was dropped: {:?}",
            l.labels.iter().map(|x| &x.text).collect::<Vec<_>>()
        );
        // Hostile limits and gaps don't panic.
        for (gap, limit) in [
            (0.0, 0.0),
            (f32::NAN, f32::INFINITY),
            (-1.0, -5.0),
            (1e-9, 1e9),
        ] {
            let opts = CircularOptions {
                label_gap: gap,
                label_limit: limit,
                ..CircularOptions::default()
            };
            let _ = circular(&m, &opts, &sites);
        }
    }

    #[test]
    fn site_labels_join_the_spread() {
        let m = plasmid(&[("f", 100, 140, Strand::None)]);
        let sites: Vec<SiteMark> = (0..6)
            .map(|i| SiteMark {
                pos: 110 + i * 2,
                text: format!("E{i}I ({})", 110 + i * 2),
                keep: false,
            })
            .collect();
        let opts = CircularOptions::default();
        let l = circular(&m, &opts, &sites);
        assert_eq!(l.labels.len(), 7);
        assert_eq!(l.labels[1].target, LabelTarget::Site(0));
        let anchor = l.labels[1].anchor;
        assert!(((anchor.0 * anchor.0 + anchor.1 * anchor.1).sqrt() - 1.0).abs() < 1e-5);
        let mut ys: Vec<f32> = l.labels.iter().map(|x| x.at.1).collect();
        ys.sort_by(f32::total_cmp);
        for w in ys.windows(2) {
            assert!(w[1] - w[0] >= opts.label_gap - 1e-5, "{ys:?}");
        }
        // Positions past the end (a stale mark) don't break anything.
        let stale = [SiteMark {
            pos: 5000,
            text: "X".into(),
            keep: false,
        }];
        assert_eq!(circular(&m, &opts, &stale).labels.len(), 2);
    }

    #[test]
    fn features_can_be_left_out() {
        let m = plasmid(&[("f", 100, 140, Strand::None)]);
        let sites = [SiteMark {
            pos: 500,
            text: "E".into(),
            keep: false,
        }];
        let opts = CircularOptions {
            features: false,
            ..CircularOptions::default()
        };
        let l = circular(&m, &opts, &sites);
        assert!(l.arcs.is_empty());
        assert_eq!(l.labels.len(), 1);
        assert_eq!(l.labels[0].target, LabelTarget::Site(0));
    }

    #[test]
    fn empty_molecule_is_fine() {
        let seq = Seq::new(Alphabet::Dna, "").unwrap();
        let m = Molecule::new("e", seq, Topology::Circular).unwrap();
        let l = circular(&m, &CircularOptions::default(), &[]);
        assert!(l.arcs.is_empty() && l.ticks.is_empty());
        assert_eq!(l.position_at(1.0), 0);
    }

    #[test]
    fn only_the_last_segment_has_the_arrowhead() {
        let seq = Seq::new(Alphabet::Dna, "A".repeat(1000)).unwrap();
        let mut m = Molecule::new("p", seq, Topology::Circular).unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "join".into(),
            kind: "CDS".into(),
            strand: Strand::Reverse,
            segments: vec![
                Segment {
                    range: Range::new(500, 600),
                },
                Segment {
                    range: Range::new(100, 200),
                },
            ],
            qualifiers: vec![],
        })
        .unwrap();
        let l = circular(&m, &CircularOptions::default(), &[]);
        assert_eq!(l.arcs[0].arrow, Arrow::None);
        assert_eq!(l.arcs[1].arrow, Arrow::Start);
    }
}
