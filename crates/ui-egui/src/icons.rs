//! Icons drawn in code (AGENTS.md › Assets and brand): simple geometric glyphs of our own
//! design, so there are no image assets to license.

use eframe::egui::{Color32, Painter, Pos2, Rect, Stroke, pos2, vec2};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Icon {
    Open,
    Save,
    Undo,
    Redo,
    ReverseComplement,
    Theme,
    Sidebar,
    /// The Experiments menu.
    Experiments,
    /// The Preferences window.
    Preferences,
    /// Map layers.
    Features,
    Enzymes,
    Primers,
    Positions,
    Ruler,
    Orfs,
    /// The sequence on one line, scrolled sideways, rather than wrapped into rows.
    SingleLine,
    /// Sequences aligned under the document ("Align with").
    AlignWith,
}

/// Draw `icon` inside `rect` with `color`.
pub fn paint(painter: &Painter, rect: Rect, icon: Icon, color: Color32) {
    let s = rect.width().min(rect.height());
    let c = rect.center();
    let stroke = Stroke::new((s / 11.0).max(1.2), color);
    let p = |x: f32, y: f32| pos2(c.x + x * s / 2.0, c.y + y * s / 2.0);
    let poly = |pts: &[Pos2]| {
        for w in pts.windows(2) {
            if let [a, b] = w {
                painter.line_segment([*a, *b], stroke);
            }
        }
    };
    match icon {
        Icon::Open => {
            // A folder.
            poly(&[
                p(-0.8, 0.75),
                p(-0.8, -0.6),
                p(-0.3, -0.6),
                p(-0.15, -0.4),
                p(0.8, -0.4),
                p(0.8, 0.75),
                p(-0.8, 0.75),
            ]);
            poly(&[p(-0.8, -0.1), p(0.8, -0.1)]);
        }
        Icon::Save => {
            // A tray with a downward arrow.
            poly(&[p(-0.8, 0.2), p(-0.8, 0.8), p(0.8, 0.8), p(0.8, 0.2)]);
            poly(&[p(0.0, -0.85), p(0.0, 0.4)]);
            poly(&[p(-0.4, 0.0), p(0.0, 0.4), p(0.4, 0.0)]);
        }
        Icon::Undo | Icon::Redo => {
            let dir = if icon == Icon::Undo { -1.0 } else { 1.0 };
            // An arc over the top, ending in an arrowhead on the left (undo) or right (redo).
            let r = s * 0.42;
            let centre = c + vec2(0.0, s * 0.12);
            let steps = 12;
            let points: Vec<Pos2> = (0..=steps)
                .map(|i| {
                    let t = std::f32::consts::PI * (i as f32 / steps as f32);
                    centre + vec2(dir * r * t.cos(), -r * t.sin())
                })
                .collect();
            poly(&points);
            if let Some(&tip) = points.first() {
                // A downward arrowhead at the end of the arc.
                poly(&[
                    tip + vec2(-s * 0.2, -s * 0.2),
                    tip,
                    tip + vec2(s * 0.2, -s * 0.2),
                ]);
            }
        }
        Icon::ReverseComplement => {
            // Two opposed arrows (both strands).
            poly(&[p(-0.8, -0.35), p(0.8, -0.35)]);
            poly(&[p(0.45, -0.65), p(0.8, -0.35), p(0.45, -0.05)]);
            poly(&[p(0.8, 0.4), p(-0.8, 0.4)]);
            poly(&[p(-0.45, 0.1), p(-0.8, 0.4), p(-0.45, 0.7)]);
        }
        Icon::Sidebar => {
            // A window with a panel on the left.
            poly(&[
                p(-0.8, -0.7),
                p(0.8, -0.7),
                p(0.8, 0.7),
                p(-0.8, 0.7),
                p(-0.8, -0.7),
            ]);
            poly(&[p(-0.25, -0.7), p(-0.25, 0.7)]);
        }
        Icon::Theme => {
            // A circle, half filled.
            painter.circle_stroke(c, s * 0.4, stroke);
            let half: Vec<Pos2> = (0..=16)
                .map(|i| {
                    let t = std::f32::consts::PI * (i as f32 / 16.0) - std::f32::consts::FRAC_PI_2;
                    c + vec2(s * 0.4 * t.cos(), s * 0.4 * t.sin())
                })
                .collect();
            painter.add(eframe::egui::Shape::convex_polygon(
                half,
                color,
                Stroke::NONE,
            ));
        }
        Icon::Features => {
            // A feature arrow: a bar ending in a point.
            poly(&[
                p(-0.8, -0.35),
                p(0.35, -0.35),
                p(0.8, 0.0),
                p(0.35, 0.35),
                p(-0.8, 0.35),
                p(-0.8, -0.35),
            ]);
        }
        Icon::Enzymes => {
            // Two strands with a staggered cut through them.
            poly(&[p(-0.85, -0.45), p(0.85, -0.45)]);
            poly(&[p(-0.85, 0.45), p(0.85, 0.45)]);
            poly(&[p(-0.3, -0.8), p(-0.3, 0.0), p(0.3, 0.0), p(0.3, 0.8)]);
        }
        Icon::Primers => {
            // A template strand with a primer above it: a half arrow.
            poly(&[p(-0.85, 0.45), p(0.85, 0.45)]);
            poly(&[p(-0.6, -0.1), p(0.6, -0.1), p(0.25, -0.45)]);
        }
        Icon::Positions => {
            // A number in brackets.
            painter.text(
                c,
                eframe::egui::Align2::CENTER_CENTER,
                "(12)",
                eframe::egui::FontId::proportional(s * 0.5),
                color,
            );
        }
        Icon::Orfs => {
            // Two reading frames: one each way, each from a start tick to an arrowhead.
            poly(&[p(-0.8, -0.6), p(-0.8, -0.2)]);
            poly(&[p(-0.8, -0.4), p(0.7, -0.4)]);
            poly(&[p(0.4, -0.65), p(0.75, -0.4), p(0.4, -0.15)]);
            poly(&[p(0.8, 0.2), p(0.8, 0.6)]);
            poly(&[p(0.8, 0.4), p(-0.7, 0.4)]);
            poly(&[p(-0.4, 0.15), p(-0.75, 0.4), p(-0.4, 0.65)]);
        }
        Icon::Experiments => {
            // A conical flask: neck, sloping sides, flat base, and a liquid line.
            poly(&[
                p(-0.25, -0.85),
                p(-0.25, -0.3),
                p(-0.8, 0.8),
                p(0.8, 0.8),
                p(0.25, -0.3),
                p(0.25, -0.85),
            ]);
            poly(&[p(-0.38, -0.85), p(0.38, -0.85)]);
            poly(&[p(-0.55, 0.3), p(0.55, 0.3)]);
        }
        Icon::Preferences => {
            // A gear: a ring with eight teeth around it.
            painter.circle_stroke(c, s * 0.28, stroke);
            for i in 0..8 {
                let a = std::f32::consts::TAU * (i as f32) / 8.0;
                let (sin, cos) = a.sin_cos();
                painter.line_segment(
                    [
                        c + vec2(cos, sin) * (s * 0.28),
                        c + vec2(cos, sin) * (s * 0.46),
                    ],
                    Stroke::new(stroke.width * 1.6, color),
                );
            }
        }
        Icon::Ruler => {
            // A scale with ticks.
            poly(&[p(-0.85, 0.3), p(0.85, 0.3)]);
            for (x, h) in [(-0.8, 0.6), (-0.4, 0.3), (0.0, 0.6), (0.4, 0.3), (0.8, 0.6)] {
                poly(&[p(x, 0.3), p(x, 0.3 - h)]);
            }
        }
        Icon::SingleLine => {
            // One long line running out past both edges, with arrowheads either way.
            poly(&[p(-0.9, 0.0), p(0.9, 0.0)]);
            poly(&[p(-0.5, -0.35), p(-0.9, 0.0), p(-0.5, 0.35)]);
            poly(&[p(0.5, -0.35), p(0.9, 0.0), p(0.5, 0.35)]);
        }
        Icon::AlignWith => {
            // A long line with shorter ones staggered under it, as reads lie under a construct.
            poly(&[p(-0.9, -0.6), p(0.9, -0.6)]);
            poly(&[p(-0.7, -0.05), p(0.3, -0.05)]);
            poly(&[p(-0.2, 0.5), p(0.8, 0.5)]);
        }
    }
}
