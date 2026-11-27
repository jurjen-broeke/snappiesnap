//! The navigator under the sequence view when it shows the whole molecule on one line: the
//! molecule drawn to the width of the view, its features as plain rectangles in their map
//! colours, and a box marking the part of the line on screen. Clicking or dragging in it moves
//! the view there, and hovering names the features under the pointer.
//!
//! Under the features it can show reads aligned to the molecule (the Align with pane): one
//! thin lane each, in the strand's colour, with their differences marked.

use eframe::egui::{Pos2, Rect, Sense, Stroke, StrokeKind, Ui, pos2, vec2};
use snappiesnap_engine::doc::{Molecule, Range};
use snappiesnap_layout::{self as layout, TrackBar};

use crate::theme::Tokens;

/// Feature lanes drawn; features in deeper lanes share the last one.
const LANES: usize = 4;
const LANE_H: f32 = 5.0;
const LANE_GAP: f32 = 2.0;
/// Space above and below the content, and on either side of the molecule.
const PAD_Y: f32 = 6.0;
const PAD_X: f32 = 8.0;
const BACKBONE_H: f32 = 3.0;
/// The narrowest a feature or the view box is drawn, so a 20 bp feature on a chromosome is
/// still there to find.
const MIN_FEATURE_W: f32 = 2.0;
const MIN_WINDOW_W: f32 = 6.0;
/// Read lanes drawn; further reads share the last one.
const READ_LANES: usize = 6;

/// A read aligned to the molecule, as the navigator draws it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ReadSpan {
    pub name: String,
    /// The stretches of the molecule it covers: two when it runs across the origin.
    pub pieces: Vec<(usize, usize)>,
    pub forward: bool,
    /// Molecule positions where it differs.
    pub differences: Vec<usize>,
}

/// Navigator state that lives across frames.
#[derive(Debug, Default, Clone)]
pub struct Navigator {
    /// Every feature piece of the molecule, laid out in lanes, for this document revision.
    tracks: Option<(u64, Vec<TrackBar>)>,
    /// Where in the view box a drag took hold, in residues from its middle.
    grab: Option<f32>,
    /// Centre the view on this position on the next frame: the user clicked or dragged here.
    pub target: Option<f32>,
    /// Reads to draw under the features (empty: none).
    pub reads: Vec<ReadSpan>,
    /// Feature lanes drawn, as `prepare` worked out.
    feature_lanes: usize,
}

impl Navigator {
    /// Lay out the molecule's features (once per revision) and return the navigator's height.
    pub fn prepare(&mut self, m: &Molecule, revision: u64, features: bool) -> f32 {
        if self.tracks.as_ref().is_none_or(|(r, _)| *r != revision) {
            self.tracks = Some((revision, layout::line_tracks(m, 0, m.len())));
        }
        let lanes = if features {
            self.bars()
                .iter()
                .map(|b| b.lane + 1)
                .max()
                .unwrap_or(0)
                .min(LANES)
        } else {
            0
        };
        self.feature_lanes = lanes;
        height(lanes) + self.reads.len().min(READ_LANES) as f32 * (LANE_H + LANE_GAP)
    }

    fn bars(&self) -> &[TrackBar] {
        self.tracks.as_ref().map_or(&[], |(_, b)| b.as_slice())
    }

    /// Draw the navigator across `ui`'s width. `window` is the part of the line on screen, in
    /// residues (fractional at the edges).
    #[allow(clippy::too_many_arguments)]
    pub fn show(
        &mut self,
        ui: &mut Ui,
        t: &Tokens,
        m: &Molecule,
        features: bool,
        height: f32,
        window: Option<(f32, f32)>,
        selection: Option<Range>,
    ) {
        let (rect, response) =
            ui.allocate_exact_size(vec2(ui.available_width(), height), Sense::click_and_drag());
        let len = m.len();
        if len == 0 || !ui.is_rect_visible(rect) {
            return;
        }
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, t.panel);
        painter.hline(rect.x_range(), rect.top(), Stroke::new(1.0, t.border));
        let strip = Strip::new(rect, len);

        // The molecule.
        let backbone_y = rect.top() + PAD_Y;
        painter.rect_filled(
            Rect::from_min_max(
                pos2(strip.left, backbone_y),
                pos2(strip.right, backbone_y + BACKBONE_H),
            ),
            1.0,
            t.backbone,
        );
        // The selection on it, in the sequence view's selection colour.
        if let Some(sel) = selection.filter(|r| !r.is_empty()) {
            let (a, b) = sel.segments(len);
            for piece in std::iter::once(a).chain(b) {
                let (x0, x1) = strip.span(piece.start, piece.end, MIN_FEATURE_W);
                painter.rect_filled(
                    Rect::from_min_max(
                        pos2(x0, backbone_y - 2.0),
                        pos2(x1, backbone_y + BACKBONE_H + 2.0),
                    ),
                    0.0,
                    t.selection,
                );
            }
        }

        // Features: rectangles in their colours, lane by lane.
        let lanes_top = backbone_y + BACKBONE_H + LANE_GAP + 1.0;
        let feature_rect = |b: &TrackBar| {
            let (x0, x1) = strip.span(b.start, b.end, MIN_FEATURE_W);
            let y = lanes_top + b.lane.min(LANES - 1) as f32 * (LANE_H + LANE_GAP);
            Rect::from_min_max(pos2(x0, y), pos2(x1, y + LANE_H))
        };
        if features {
            for b in self.bars() {
                painter.rect_filled(feature_rect(b), 1.0, t.color_at(m, b.feature));
            }
        }
        // Reads, a lane each under the features, their differences as ticks across them.
        let reads_top = lanes_top + self.feature_lanes as f32 * (LANE_H + LANE_GAP);
        let read_lane = |k: usize| {
            let y = reads_top + k.min(READ_LANES - 1) as f32 * (LANE_H + LANE_GAP);
            (y, y + LANE_H)
        };
        for (k, read) in self.reads.iter().enumerate() {
            let (y0, y1) = read_lane(k);
            let colour = t.feature_color(if read.forward { 0 } else { 1 });
            for &(a, b) in &read.pieces {
                let (x0, x1) = strip.span(a, b, MIN_FEATURE_W);
                painter.rect_filled(
                    Rect::from_min_max(pos2(x0, y0 + 1.0), pos2(x1, y1 - 1.0)),
                    1.0,
                    colour,
                );
            }
            for &d in &read.differences {
                let x = strip.x(d as f32 + 0.5);
                painter.vline(x, y0 - 1.0..=y1 + 1.0, Stroke::new(1.0, t.error));
            }
        }

        // The part of the line on screen.
        let shown = window.map(|(a, b)| {
            let (x0, x1) = (strip.x(a), strip.x(b));
            let mid = (x0 + x1) / 2.0;
            let half = ((x1 - x0) / 2.0).max(MIN_WINDOW_W / 2.0);
            Rect::from_min_max(
                pos2(mid - half, rect.top() + 2.0),
                pos2(mid + half, rect.bottom() - 2.0),
            )
        });
        if let Some(r) = shown {
            painter.rect_filled(r, 2.0, t.accent.gamma_multiply(0.18));
            painter.rect_stroke(r, 2.0, Stroke::new(1.5, t.accent), StrokeKind::Inside);
        }

        // Click or drag: move the view. A drag that takes hold of the box keeps its grip.
        if let Some(p) = response.interact_pointer_pos() {
            let at = strip.pos_at(p.x);
            if response.drag_started() || response.clicked() {
                self.grab = match (window, shown) {
                    (Some((a, b)), Some(r)) if response.drag_started() && r.contains(p) => {
                        Some(at - (a + b) / 2.0)
                    }
                    _ => None,
                };
            }
            if response.clicked() || response.dragged() {
                self.target = Some(at - self.grab.unwrap_or(0.0));
                ui.ctx().request_repaint();
            }
        }
        if response.drag_stopped() {
            self.grab = None;
        }

        // Hover: the position, and the features under the pointer.
        if let Some(p) = response.hover_pos() {
            let at = (strip.pos_at(p.x) as usize).min(len.saturating_sub(1)) + 1;
            let mut names: Vec<String> = Vec::new();
            if features {
                for b in self.bars() {
                    if hit(feature_rect(b), p)
                        && let Some(f) = m.features().get(b.feature)
                    {
                        let name = if f.name.is_empty() { &f.kind } else { &f.name };
                        if !names.contains(name) {
                            names.push(name.clone());
                        }
                    }
                }
            }
            for (k, read) in self.reads.iter().enumerate() {
                let (y0, y1) = read_lane(k);
                let on = read.pieces.iter().any(|&(a, b)| {
                    let (x0, x1) = strip.span(a, b, MIN_FEATURE_W);
                    hit(Rect::from_min_max(pos2(x0, y0), pos2(x1, y1)), p)
                });
                if on && !names.contains(&read.name) {
                    names.push(read.name.clone());
                }
            }
            let mut text = format!("Position {at}");
            for name in names.iter().take(8) {
                text.push('\n');
                text.push_str(name);
            }
            if names.len() > 8 {
                text.push_str(&format!("\n… and {} more", names.len() - 8));
            }
            response.on_hover_text_at_pointer(text);
        }
    }
}

/// The navigator's height with `lanes` feature lanes.
fn height(lanes: usize) -> f32 {
    PAD_Y * 2.0 + BACKBONE_H + LANE_GAP + 1.0 + lanes as f32 * (LANE_H + LANE_GAP)
}

/// Whether `p` falls on `r`, widened to a few pixels so thin features can be pointed at.
fn hit(r: Rect, p: Pos2) -> bool {
    r.expand2(vec2(((4.0 - r.width()) / 2.0).max(0.0), 1.0))
        .contains(p)
}

/// The molecule's positions spread across the navigator.
struct Strip {
    left: f32,
    right: f32,
    len: f32,
}

impl Strip {
    fn new(rect: Rect, len: usize) -> Self {
        let left = rect.left() + PAD_X;
        Strip {
            left,
            right: (rect.right() - PAD_X).max(left + 1.0),
            len: len.max(1) as f32,
        }
    }

    /// x of the boundary before residue `pos`.
    fn x(&self, pos: f32) -> f32 {
        self.left + (pos / self.len).clamp(0.0, 1.0) * (self.right - self.left)
    }

    /// The residue position under `x`, within the molecule.
    fn pos_at(&self, x: f32) -> f32 {
        ((x - self.left) / (self.right - self.left) * self.len).clamp(0.0, self.len)
    }

    /// The x range of residues `start..end`, at least `min` wide.
    fn span(&self, start: usize, end: usize, min: f32) -> (f32, f32) {
        let (x0, x1) = (self.x(start as f32), self.x(end as f32));
        if x1 - x0 >= min {
            (x0, x1)
        } else {
            let mid = (x0 + x1) / 2.0;
            (mid - min / 2.0, mid + min / 2.0)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positions_spread_across_the_strip_and_back() {
        let strip = Strip::new(
            Rect::from_min_max(pos2(0.0, 0.0), pos2(1016.0, 40.0)),
            10_000,
        );
        assert_eq!(strip.x(0.0), 8.0);
        assert_eq!(strip.x(10_000.0), 1008.0);
        assert_eq!(strip.x(5_000.0), 508.0);
        // Past either end stays on the strip.
        assert_eq!(strip.x(-5.0), 8.0);
        assert_eq!(strip.x(20_000.0), 1008.0);
        assert_eq!(strip.pos_at(508.0), 5_000.0);
        assert_eq!(strip.pos_at(-100.0), 0.0);
        assert_eq!(strip.pos_at(5_000.0), 10_000.0);
        // A single residue on a long molecule is still a few pixels wide.
        let (a, b) = strip.span(5_000, 5_001, MIN_FEATURE_W);
        assert!((b - a - MIN_FEATURE_W).abs() < 1e-3);
    }

    #[test]
    fn hostile_sizes_stay_finite() {
        // Empty molecule, and a rect too narrow for its padding.
        for (w, len) in [(0.0, 0), (4.0, 1), (1.0, usize::MAX)] {
            let strip = Strip::new(Rect::from_min_max(pos2(0.0, 0.0), pos2(w, 10.0)), len);
            assert!(strip.right > strip.left);
            for x in [f32::NEG_INFINITY, -1.0, 0.0, 3.0, 1e9] {
                assert!(strip.pos_at(x).is_finite());
            }
            assert!(strip.x(1e30).is_finite());
        }
    }

    #[test]
    fn height_grows_with_lanes() {
        assert!(height(0) > 0.0);
        assert_eq!(height(2) - height(1), LANE_H + LANE_GAP);
    }
}
