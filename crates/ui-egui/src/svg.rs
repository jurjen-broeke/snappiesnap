//! Vector images (SVG) of what the app draws. The drawing code runs offscreen as for
//! `ui.render` ([`crate::offscreen::shapes`]), and egui's shapes are written as SVG elements:
//! lines, paths, circles, ellipses, rectangles, Bézier curves, and text as text (one `<text>`
//! per run of a row, starting where egui placed it and stretched to the width egui gave it with
//! `textLength`, so other fonts take the same room), so labels stay editable in vector editors. Meshes, bands and paint callbacks (not used by the maps) are left out.

use std::fmt::Write as _;

use eframe::egui::{Color32, Pos2, Rect, Vec2};
use eframe::epaint::{self, ClippedShape, ColorMode, PathStroke, Shape, Stroke, TextShape};

/// `#rrggbb` and an opacity for a (premultiplied) egui colour; `None` when invisible.
fn paint(c: Color32) -> Option<(String, f32)> {
    let [r, g, b, a] = c.to_srgba_unmultiplied();
    (a > 0).then(|| (format!("#{r:02x}{g:02x}{b:02x}"), f32::from(a) / 255.0))
}

fn fill_attr(c: Color32) -> String {
    match paint(c) {
        Some((hex, 1.0)) => format!(r#" fill="{hex}""#),
        Some((hex, a)) => format!(r#" fill="{hex}" fill-opacity="{a:.3}""#),
        None => r#" fill="none""#.to_string(),
    }
}

fn stroke_attr(width: f32, color: Color32) -> String {
    match paint(color) {
        Some((hex, a)) if width > 0.0 => {
            let mut s = format!(r#" stroke="{hex}" stroke-width="{width:.2}""#);
            if a < 1.0 {
                let _ = write!(s, r#" stroke-opacity="{a:.3}""#);
            }
            s
        }
        _ => String::new(),
    }
}

fn path_stroke_attr(s: &PathStroke) -> String {
    match &s.color {
        ColorMode::Solid(c) => stroke_attr(s.width, *c),
        // Gradient strokes aren't used by the maps: drawn neutral grey.
        ColorMode::UV(_) => stroke_attr(s.width, Color32::GRAY),
    }
}

fn plain_stroke_attr(s: &Stroke) -> String {
    stroke_attr(s.width, s.color)
}

fn points(ps: &[Pos2]) -> String {
    let mut out = String::new();
    for (i, p) in ps.iter().enumerate() {
        if i > 0 {
            out.push(' ');
        }
        let _ = write!(out, "{:.2},{:.2}", p.x, p.y);
    }
    out
}

/// Text escaped for XML.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            c if c.is_control() => {}
            c => out.push(c),
        }
    }
    out
}

/// The text of a text shape: per row, runs of glyphs of one section (one font and colour),
/// each a `<text>` with every character's own x position.
fn text(out: &mut String, t: &TextShape) {
    let job = &t.galley.job;
    // Glyphs come in text order (line breaks have no glyph): walk the text alongside them to
    // know each glyph's section.
    let mut chars = job
        .text
        .char_indices()
        .filter(|(_, c)| *c != '\n')
        .peekable();
    let section_at = |byte: usize| {
        job.sections
            .iter()
            .position(|s| (s.byte_range.start.0..s.byte_range.end.0).contains(&byte))
            .unwrap_or(0)
    };
    let rotate = if t.angle != 0.0 {
        format!(
            r#" transform="rotate({:.3} {:.2} {:.2})""#,
            t.angle.to_degrees(),
            t.pos.x,
            t.pos.y
        )
    } else {
        String::new()
    };
    for row in &t.galley.rows {
        let mut run: Vec<(char, f32, f32)> = Vec::new();
        let mut run_section = usize::MAX;
        let mut baseline = 0.0f32;
        let flush = |out: &mut String, run: &mut Vec<(char, f32, f32)>, section: usize, y: f32| {
            if run.is_empty() {
                return;
            }
            let Some(s) = job.sections.get(section) else {
                run.clear();
                return;
            };
            let mut color = t.override_text_color.unwrap_or(s.format.color);
            if color == Color32::PLACEHOLDER {
                color = t.fallback_color;
            }
            if t.opacity_factor < 1.0 {
                color = color.gamma_multiply(t.opacity_factor);
            }
            let family = match s.format.font_id.family {
                epaint::FontFamily::Monospace => "'Ubuntu Mono', Menlo, Consolas, monospace",
                _ => "Ubuntu, 'Helvetica Neue', Arial, sans-serif",
            };
            let (Some(first), Some(last)) = (run.first(), run.last()) else {
                return;
            };
            let x = first.1;
            let length = (last.1 + last.2 - x).max(0.0);
            let content: String = run.iter().map(|(c, _, _)| *c).collect();
            let _ = writeln!(
                out,
                r#"<text x="{x:.2}" y="{y:.2}" font-family="{family}" font-size="{:.2}" textLength="{length:.2}" lengthAdjust="spacingAndGlyphs"{}{rotate} xml:space="preserve">{}</text>"#,
                s.format.font_id.size,
                fill_attr(color),
                escape(&content),
            );
            run.clear();
        };
        for g in &row.row.glyphs {
            let byte = match chars.next() {
                Some((b, _)) => b,
                None => 0,
            };
            let section = section_at(byte);
            let x = t.pos.x + row.pos.x + g.pos.x;
            let y = t.pos.y + row.pos.y + g.pos.y;
            if section != run_section || (y - baseline).abs() > 0.01 {
                flush(out, &mut run, run_section, baseline);
                run_section = section;
                baseline = y;
            }
            run.push((g.chr, x, g.advance_width));
        }
        flush(out, &mut run, run_section, baseline);
    }
}

fn shape(out: &mut String, s: &Shape) {
    match s {
        Shape::Noop | Shape::Mesh(_) | Shape::Callback(_) | Shape::Band(_) => {}
        Shape::Vec(v) => {
            for s in v {
                shape(out, s);
            }
        }
        Shape::Circle(c) => {
            let _ = writeln!(
                out,
                r#"<circle cx="{:.2}" cy="{:.2}" r="{:.2}"{}{}/>"#,
                c.center.x,
                c.center.y,
                c.radius,
                fill_attr(c.fill),
                plain_stroke_attr(&c.stroke)
            );
        }
        Shape::Ellipse(e) => {
            let _ = writeln!(
                out,
                r#"<ellipse cx="{:.2}" cy="{:.2}" rx="{:.2}" ry="{:.2}" transform="rotate({:.3} {:.2} {:.2})"{}{}/>"#,
                e.center.x,
                e.center.y,
                e.radius.x,
                e.radius.y,
                e.angle.to_degrees(),
                e.center.x,
                e.center.y,
                fill_attr(e.fill),
                plain_stroke_attr(&e.stroke)
            );
        }
        Shape::LineSegment {
            points: [a, b],
            stroke,
        } => {
            let _ = writeln!(
                out,
                r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}"{}/>"#,
                a.x,
                a.y,
                b.x,
                b.y,
                plain_stroke_attr(stroke)
            );
        }
        Shape::Path(p) => {
            if p.points.len() < 2 {
                return;
            }
            let tag = if p.closed { "polygon" } else { "polyline" };
            let fill = if p.closed {
                fill_attr(p.fill)
            } else {
                r#" fill="none""#.to_string()
            };
            let _ = writeln!(
                out,
                r#"<{tag} points="{}"{fill}{} stroke-linejoin="round"/>"#,
                points(&p.points),
                path_stroke_attr(&p.stroke)
            );
        }
        Shape::Rect(r) => {
            let radius = f32::from(r.corner_radius.nw);
            let _ = writeln!(
                out,
                r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" rx="{radius:.2}"{}{}/>"#,
                r.rect.min.x,
                r.rect.min.y,
                r.rect.width(),
                r.rect.height(),
                fill_attr(r.fill),
                plain_stroke_attr(&r.stroke)
            );
        }
        Shape::Text(t) => text(out, t),
        Shape::QuadraticBezier(q) => {
            let [a, b, c] = q.points;
            let close = if q.closed { " Z" } else { "" };
            let fill = if q.closed {
                fill_attr(q.fill)
            } else {
                r#" fill="none""#.to_string()
            };
            let _ = writeln!(
                out,
                r#"<path d="M {:.2} {:.2} Q {:.2} {:.2} {:.2} {:.2}{close}"{fill}{}/>"#,
                a.x,
                a.y,
                b.x,
                b.y,
                c.x,
                c.y,
                path_stroke_attr(&q.stroke)
            );
        }
        Shape::CubicBezier(cb) => {
            let [a, b, c, d] = cb.points;
            let close = if cb.closed { " Z" } else { "" };
            let fill = if cb.closed {
                fill_attr(cb.fill)
            } else {
                r#" fill="none""#.to_string()
            };
            let _ = writeln!(
                out,
                r#"<path d="M {:.2} {:.2} C {:.2} {:.2} {:.2} {:.2} {:.2} {:.2}{close}"{fill}{}/>"#,
                a.x,
                a.y,
                b.x,
                b.y,
                c.x,
                c.y,
                d.x,
                d.y,
                path_stroke_attr(&cb.stroke)
            );
        }
    }
}

/// An SVG document of `shapes` on a `size` canvas filled with `background`. Shapes clipped to
/// a smaller rectangle than the canvas keep that clip.
pub fn to_svg(shapes: &[ClippedShape], size: Vec2, background: Color32, title: &str) -> String {
    let mut out = String::new();
    let _ = writeln!(out, r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    let _ = writeln!(
        out,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0}" height="{h:.0}" viewBox="0 0 {w:.2} {h:.2}">"#,
        w = size.x,
        h = size.y
    );
    let _ = writeln!(out, "<title>{}</title>", escape(title));
    let _ = writeln!(
        out,
        r#"<rect x="0" y="0" width="{:.2}" height="{:.2}"{}/>"#,
        size.x,
        size.y,
        fill_attr(background)
    );
    let canvas = Rect::from_min_size(Pos2::ZERO, size);
    let mut clips: Vec<Rect> = Vec::new();
    let mut defs = String::new();
    let mut body = String::new();
    for c in shapes {
        let mut s = String::new();
        shape(&mut s, &c.shape);
        if s.is_empty() {
            continue;
        }
        if c.clip_rect.contains_rect(canvas) || !c.clip_rect.is_finite() {
            body.push_str(&s);
            continue;
        }
        let id = match clips.iter().position(|r| *r == c.clip_rect) {
            Some(i) => i,
            None => {
                let r = c.clip_rect.intersect(canvas);
                let _ = writeln!(
                    defs,
                    r#"<clipPath id="c{}"><rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}"/></clipPath>"#,
                    clips.len(),
                    r.min.x,
                    r.min.y,
                    r.width().max(0.0),
                    r.height().max(0.0)
                );
                clips.push(c.clip_rect);
                clips.len() - 1
            }
        };
        let _ = writeln!(body, r#"<g clip-path="url(#c{id})">"#);
        body.push_str(&s);
        body.push_str("</g>\n");
    }
    if !defs.is_empty() {
        let _ = writeln!(out, "<defs>\n{defs}</defs>");
    }
    out.push_str(&body);
    out.push_str("</svg>\n");
    out
}

#[cfg(test)]
mod tests {
    use eframe::egui::{self, FontId, pos2, vec2};

    use super::*;

    #[test]
    fn shapes_become_svg_elements() {
        let shapes = crate::offscreen::shapes(
            200,
            100,
            |_| {},
            |ui| {
                let p = ui.painter();
                p.line_segment(
                    [pos2(1.0, 2.0), pos2(30.0, 40.0)],
                    Stroke::new(2.0, Color32::RED),
                );
                p.circle_filled(pos2(50.0, 50.0), 10.0, Color32::from_rgb(0, 128, 255));
                p.add(egui::Shape::convex_polygon(
                    vec![pos2(0.0, 0.0), pos2(10.0, 0.0), pos2(5.0, 8.0)],
                    Color32::GREEN,
                    Stroke::NONE,
                ));
                p.text(
                    pos2(100.0, 20.0),
                    egui::Align2::LEFT_TOP,
                    "EcoRI & <BamHI>",
                    FontId::proportional(12.0),
                    Color32::BLACK,
                );
            },
        )
        .unwrap();
        let svg = to_svg(&shapes, vec2(200.0, 100.0), Color32::WHITE, "pTest");
        assert!(svg.starts_with("<?xml"));
        assert!(svg.contains("<title>pTest</title>"));
        assert!(svg.contains(r##"<line x1="1.00" y1="2.00" x2="30.00" y2="40.00" stroke="#ff0000" stroke-width="2.00"/>"##), "{svg}");
        assert!(
            svg.contains(r##"<circle cx="50.00" cy="50.00" r="10.00" fill="#0080ff"/>"##),
            "{svg}"
        );
        assert!(
            svg.contains("<polygon points=\"0.00,0.00 10.00,0.00 5.00,8.00\" fill=\"#00ff00\""),
            "{svg}"
        );
        assert!(svg.contains(">EcoRI &amp; &lt;BamHI&gt;</text>"), "{svg}");
        assert!(svg.contains(r#"font-size="12.00""#));
        assert!(svg.ends_with("</svg>\n"));
    }

    #[test]
    fn escaping_and_invisible_colours() {
        assert_eq!(escape("a\"b\u{7}c"), "a&quot;bc");
        assert_eq!(fill_attr(Color32::TRANSPARENT), r#" fill="none""#);
        assert_eq!(stroke_attr(0.0, Color32::RED), "");
        assert!(
            fill_attr(Color32::from_rgba_unmultiplied(255, 0, 0, 128)).contains("fill-opacity")
        );
    }
}
