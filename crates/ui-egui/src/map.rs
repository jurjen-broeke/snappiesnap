//! The map view: circular or linear, drawn from `snappiesnap_layout` geometry with the theme
//! tokens. It only reads the document; clicks become [`MapAction`]s that the shell turns into
//! engine commands. Which layers are drawn (features, restriction sites, primers, positions,
//! ruler) is up to the caller.

use std::f32::consts::TAU;

use eframe::egui::{self, Align2, FontId, Pos2, Rect, Sense, Shape, Stroke, Ui, Vec2, pos2, vec2};
use snappiesnap_engine::doc::{Molecule, Range, Strand, Topology};
use snappiesnap_layout::{self as layout, Arrow, CircularOptions, LabelTarget, SiteMark};

use crate::enzymes::{MapSite, Marks, OrfMark, Pairing, PrimerMark, SiteDetail};
use crate::state::MapStyle;
use crate::theme::Tokens;

/// What the user asked for by interacting with the map.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MapAction {
    /// Select this range (`None` clears the selection).
    Select(Option<Range>),
    /// A restriction site's label was double-clicked: its index in the sites given. The
    /// cloning experiment adds that enzyme to the cut list; the document views ignore it.
    SiteChosen(usize),
    /// A primer's label was double-clicked: its index in the primer marks given. The cloning
    /// experiment puts it in the forward or reverse field, by the strand it binds.
    PrimerChosen(usize),
}

/// Which layers the map draws.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Layers {
    pub features: bool,
    pub enzymes: bool,
    pub primers: bool,
    /// Cut positions in restriction site labels.
    pub positions: bool,
    /// Ticks with coordinates.
    pub ruler: bool,
    /// Open reading frames.
    pub orfs: bool,
    /// Name labels (off for thumbnails).
    pub labels: bool,
}

impl Default for Layers {
    fn default() -> Self {
        Layers {
            features: true,
            enzymes: true,
            primers: true,
            positions: true,
            ruler: true,
            orfs: false,
            labels: true,
        }
    }
}

/// A labelled mark on the map: a restriction site or a primer binding site.
struct Mark {
    pos: usize,
    text: String,
    range: Range,
    kind: MarkKind,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum MarkKind {
    Site(usize),
    Primer(usize),
}

/// The marks for the layers that are on, in position order.
fn marks(m: &Molecule, sites: &[MapSite], primers: &[PrimerMark], layers: Layers) -> Vec<Mark> {
    let mut out = Vec::new();
    if layers.enzymes {
        for (i, s) in sites.iter().enumerate() {
            out.push(Mark {
                pos: s.mark.pos,
                text: if layers.positions {
                    s.mark.text.clone()
                } else {
                    s.detail.label.clone()
                },
                range: s.range,
                kind: MarkKind::Site(i),
            });
        }
    }
    if layers.primers {
        for (i, p) in primers.iter().enumerate() {
            out.push(Mark {
                pos: p.mid(m.len()),
                text: p.name.clone(),
                range: p.range,
                kind: MarkKind::Primer(i),
            });
        }
    }
    out.sort_by_key(|k| k.pos);
    out
}

/// Draw the map; returns an action if the user clicked something.
#[allow(clippy::too_many_arguments)]
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    m: &Molecule,
    selection: Option<Range>,
    style: MapStyle,
    all: &Marks,
    layers: Layers,
) -> Option<MapAction> {
    let (sites, primers) = (all.sites, all.primers);
    let orfs: &[OrfMark] = if layers.orfs { all.orfs } else { &[] };
    let circular = match style {
        MapStyle::Auto => m.topology() == Topology::Circular,
        MapStyle::Circular => true,
        MapStyle::Linear => false,
    };
    let marks = marks(m, sites, primers, layers);
    let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::click());
    let painter = ui.painter_at(rect);
    let hover = response.hover_pos();
    let click = response
        .clicked()
        .then(|| response.interact_pointer_pos())
        .flatten();
    let view = View {
        t,
        m,
        marks: &marks,
        sites,
        primers,
        orfs,
        layers,
        hover,
        click,
    };
    let (hovered, clicked_range, clicked_empty) = if circular {
        draw_circular(&painter, rect, &view, selection)
    } else {
        draw_linear(&painter, rect, &view, selection)
    };
    match hovered {
        // ORFs aren't labelled; their hover targets are numbered after the marks.
        Some(LabelTarget::Site(i)) if i >= marks.len() => {
            if let Some(o) = orfs.get(i - marks.len()) {
                response
                    .clone()
                    .on_hover_ui_at_pointer(|ui| orf_card(ui, t, m, &o.orf));
            }
        }
        Some(LabelTarget::Site(i)) => match marks.get(i).map(|k| k.kind) {
            Some(MarkKind::Site(j)) => {
                if let Some(site) = sites.get(j) {
                    response
                        .clone()
                        .on_hover_ui_at_pointer(|ui| site_card(ui, t, &site.detail));
                }
            }
            Some(MarkKind::Primer(j)) => {
                if let Some(p) = primers.get(j) {
                    response
                        .clone()
                        .on_hover_ui_at_pointer(|ui| primer_card(ui, t, p, m.len()));
                }
            }
            None => {}
        },
        Some(LabelTarget::Feature(index)) => {
            if let Some(f) = m.features().get(index) {
                let loc: Vec<String> = f
                    .segments
                    .iter()
                    .map(|s| {
                        format!(
                            "{}..{}",
                            s.range.start.saturating_add(1),
                            end1(s.range, m.len())
                        )
                    })
                    .collect();
                let text = format!("{}\n{} · {}", f.name, f.kind, loc.join(", "));
                response.clone().on_hover_text_at_pointer(text);
            }
        }
        None => {}
    }
    // A double-click on a site's label chooses that enzyme. The label under the pointer is
    // the one that was hovered, so no second hit test is needed.
    if response.double_clicked()
        && let Some(LabelTarget::Site(i)) = hovered
    {
        match marks.get(i).map(|k| k.kind) {
            Some(MarkKind::Site(j)) => return Some(MapAction::SiteChosen(j)),
            Some(MarkKind::Primer(j)) => return Some(MapAction::PrimerChosen(j)),
            None => {}
        }
    }
    if let Some(r) = clicked_range {
        return Some(MapAction::Select(Some(r)));
    }
    if clicked_empty {
        return Some(MapAction::Select(None));
    }
    None
}

/// Everything the drawing functions read.
struct View<'a> {
    t: &'a Tokens,
    m: &'a Molecule,
    marks: &'a [Mark],
    sites: &'a [MapSite],
    primers: &'a [PrimerMark],
    orfs: &'a [OrfMark],
    layers: Layers,
    hover: Option<Pos2>,
    click: Option<Pos2>,
}

impl View<'_> {
    /// The range a label stands for.
    fn label_range(&self, target: LabelTarget) -> Option<Range> {
        match target {
            LabelTarget::Feature(i) => feature_range(self.m, i),
            LabelTarget::Site(i) => self.marks.get(i).map(|k| k.range),
        }
    }

    /// Whether this label is a restriction site the caller marked as chosen.
    fn site_chosen(&self, target: LabelTarget) -> bool {
        let LabelTarget::Site(i) = target else {
            return false;
        };
        match self.marks.get(i).map(|k| k.kind) {
            Some(MarkKind::Site(j)) => self.sites.get(j).is_some_and(|s| s.chosen),
            _ => false,
        }
    }

    fn primer_of(&self, target: LabelTarget) -> Option<&PrimerMark> {
        match target {
            LabelTarget::Site(i) => match self.marks.get(i)?.kind {
                MarkKind::Primer(j) => self.primers.get(j),
                MarkKind::Site(_) => None,
            },
            LabelTarget::Feature(_) => None,
        }
    }

    /// Primer colour: solid for primary sites, faint for secondary ones, accent when hot.
    fn primer_color(&self, p: &PrimerMark, hot: bool) -> egui::Color32 {
        if hot {
            self.t.accent
        } else if p.primary {
            self.t.text
        } else {
            self.t.text_muted.gamma_multiply(0.6)
        }
    }
}

/// 1-based inclusive end of a range.
fn end1(r: Range, len: usize) -> usize {
    if r.end == 0 { len } else { r.end }
}

/// What the pointer is over, what was clicked, and whether a click hit nothing.
type Hits = (Option<LabelTarget>, Option<Range>, bool);

/// Room kept beside a circular map for labels (px); longer labels are shortened.
const LABEL_ROOM: f32 = 160.0;

/// One line of text, cut with an ellipsis to fit `width`.
fn elided(
    painter: &egui::Painter,
    text: &str,
    font: FontId,
    width: f32,
    color: egui::Color32,
) -> std::sync::Arc<egui::Galley> {
    let mut job = egui::text::LayoutJob::simple_singleline(text.to_string(), font, color);
    job.wrap = egui::text::TextWrapping {
        max_width: width,
        max_rows: 1,
        break_anywhere: true,
        overflow_character: Some('\u{2026}'),
    };
    painter.layout_job(job)
}

/// Gap between the backbone and the first primer lane, and between lanes (circular, px).
const PRIMER_GAP: f32 = 9.0;
const PRIMER_LANE: f32 = 6.0;

fn draw_circular(painter: &egui::Painter, rect: Rect, v: &View, selection: Option<Range>) -> Hits {
    let (t, m) = (v.t, v.m);
    let c = rect.center();
    // Leave room beside the circle for the labels (up to a typical width; longer ones are
    // shortened at the canvas's edge), as they sit about 1.24 radii from the centre at
    // mid-height.
    let widest = if v.layers.labels {
        let site_font = FontId::proportional(11.0);
        let feature_font = FontId::proportional(12.0);
        let texts = v.marks.iter().map(|k| (k.text.as_str(), &site_font));
        let names = m
            .features()
            .iter()
            .filter(|_| v.layers.features)
            .map(|f| (f.name.as_str(), &feature_font));
        texts
            .chain(names)
            .map(|(text, font)| {
                painter
                    .layout_no_wrap(text.to_string(), font.clone(), t.text)
                    .size()
                    .x
            })
            .fold(0.0_f32, f32::max)
    } else {
        0.0
    };
    let beside = (rect.width() / 2.0 - widest.min(LABEL_ROOM) - 16.0) / 1.24;
    let radius = (rect.width().min(rect.height()) * 0.30)
        .min(beside)
        .max(40.0);
    // Labels need at least a text line of vertical space, whatever the map's size, and their
    // text stays inside the canvas: half a text line plus a little room at either edge.
    let label_limit = |room: f32| ((rect.height() / 2.0 - room) / radius).max(1.0);
    let mut opts = CircularOptions {
        label_gap: (15.0 / radius).max(CircularOptions::default().label_gap),
        label_limit: label_limit(12.0),
        features: v.layers.features,
        ..CircularOptions::default()
    };
    let site_marks: Vec<SiteMark> = v
        .marks
        .iter()
        .enumerate()
        .map(|(i, k)| SiteMark {
            pos: k.pos,
            text: k.text.clone(),
            // An enzyme named for this step keeps its label when the map runs out of room.
            keep: v.site_chosen(LabelTarget::Site(i)),
        })
        .collect();
    let mut l = layout::circular(m, &opts, &site_marks);
    // When labels are left out, make room under them for the line saying so.
    if v.layers.labels && l.hidden_labels > 0 {
        opts.label_limit = label_limit(30.0);
        l = layout::circular(m, &opts, &site_marks);
    }
    let at = |u: (f32, f32)| c + vec2(u.0 * radius, u.1 * radius);
    let thick = (radius * 0.07).clamp(7.0, 18.0);
    let len = l.len.max(1) as f32;

    // Selection band under everything else.
    if let Some(sel) = selection.filter(|s| !s.is_empty()) {
        let start = l.angle(sel.start as f32);
        let sweep = sel.len_on(l.len) as f32 / len * TAU;
        painter.add(Shape::line(
            arc_points(c, radius * 0.92, start, sweep),
            Stroke::new(thick * 2.6, t.selection),
        ));
    }

    // Backbone: a double line.
    for r in [radius, radius + 4.0] {
        painter.circle_stroke(c, r, Stroke::new(1.4, t.backbone));
    }
    // Origin marker and ticks.
    painter.line_segment(
        [
            at((0.0, -1.0)) - vec2(0.0, 6.0),
            at((0.0, -1.0)) + vec2(0.0, 6.0),
        ],
        Stroke::new(2.0, t.backbone),
    );
    // Tick numbers sit outside the primer lanes.
    let primer_lanes = v
        .marks
        .iter()
        .filter_map(|k| match k.kind {
            MarkKind::Primer(j) => v.primers.get(j).map(|p| p.lane + 1),
            MarkKind::Site(_) => None,
        })
        .max()
        .unwrap_or(0);
    let tick_label = radius
        + if primer_lanes > 0 {
            4.0 + PRIMER_GAP + PRIMER_LANE * primer_lanes as f32 + 10.0
        } else {
            22.0
        };
    if v.layers.ruler {
        for (tick, angle) in &l.ticks {
            let dir = vec2(angle.cos(), angle.sin());
            painter.line_segment(
                [c + dir * (radius + 4.0), c + dir * (radius + 10.0)],
                Stroke::new(1.0, t.backbone),
            );
            painter.text(
                c + dir * tick_label,
                Align2::CENTER_CENTER,
                tick.label.to_string(),
                FontId::proportional(10.0),
                t.text_muted,
            );
        }
    }

    // Features.
    let mut hovered = None;
    let mut clicked = None;
    for arc in &l.arcs {
        let color = t.color_at(v.m, arc.feature);
        let r = radius * arc.radius;
        let head = ((thick * 1.3) / r).min(arc.sweep * 0.6);
        let (body_start, body_sweep) = match arc.arrow {
            Arrow::End => (arc.start, arc.sweep - head),
            Arrow::Start => (arc.start + head, arc.sweep - head),
            Arrow::None => (arc.start, arc.sweep),
        };
        if body_sweep * r >= 1.0 {
            painter.add(Shape::line(
                arc_points(c, r, body_start, body_sweep),
                Stroke::new(thick, color),
            ));
        }
        // Arrowheads under 2 px would be near-degenerate triangles, which anti-aliasing turns
        // into spikes; such features are too small to show direction anyway.
        let tip = match arc.arrow {
            _ if head * r < 2.0 => None,
            Arrow::End => Some((arc.start + arc.sweep, arc.start + arc.sweep - head)),
            Arrow::Start => Some((arc.start, arc.start + head)),
            Arrow::None => None,
        };
        if let Some((tip_angle, base_angle)) = tip {
            let dir = |a: f32| vec2(a.cos(), a.sin());
            painter.add(Shape::convex_polygon(
                vec![
                    c + dir(base_angle) * (r + thick * 0.85),
                    c + dir(tip_angle) * r,
                    c + dir(base_angle) * (r - thick * 0.85),
                ],
                color,
                Stroke::NONE,
            ));
        }
        let on_arc = |p: Pos2| {
            let d = p - c;
            let within_ring = (d.length() - r).abs() <= thick * 0.75;
            let a = d.y.atan2(d.x);
            within_ring && in_sweep(a, arc.start, arc.sweep)
        };
        if v.hover.is_some_and(on_arc) {
            hovered = Some(LabelTarget::Feature(arc.feature));
        }
        if v.click.is_some_and(on_arc) {
            clicked = feature_range(m, arc.feature);
        }
    }

    // ORFs: thin arrows inside the feature lanes, one lane each (hover targets after the marks).
    let orf_base =
        (radius * (opts.lane0 - opts.lane_step * l.lanes.max(1) as f32) - thick * 0.6 - 6.0)
            .max(radius * 0.35);
    for (j, o) in v.orfs.iter().enumerate() {
        let r = (orf_base - 5.0 * o.lane as f32).max(radius * 0.2);
        let start = l.angle(o.orf.range.start as f32);
        let sweep = o.orf.range.len_on(l.len) as f32 / len * TAU;
        let target = LabelTarget::Site(v.marks.len() + j);
        let on = |q: Pos2| {
            let d = q - c;
            (d.length() - r).abs() <= 3.0 && in_sweep(d.y.atan2(d.x), start, sweep)
        };
        if v.hover.is_some_and(on) {
            hovered = Some(target);
        }
        if v.click.is_some_and(on) {
            clicked = Some(o.orf.range);
        }
        let color = if hovered == Some(target) {
            t.accent
        } else {
            t.text_muted
        };
        draw_orf_arc(
            painter,
            c,
            r,
            start,
            sweep,
            o.orf.strand == Strand::Reverse,
            color,
        );
    }

    // Primer binding sites: thin arrows outside the backbone, one lane each.
    let primer_radius = |p: &PrimerMark| radius + 4.0 + PRIMER_GAP + PRIMER_LANE * p.lane as f32;
    for (i, k) in v.marks.iter().enumerate() {
        let MarkKind::Primer(j) = k.kind else {
            continue;
        };
        let Some(p) = v.primers.get(j) else { continue };
        let r = primer_radius(p);
        let start = l.angle(p.range.start as f32);
        let sweep = (p.annealed as f32 / len * TAU).max(0.004);
        let on = |q: Pos2| {
            let d = q - c;
            (d.length() - r).abs() <= 4.0 && in_sweep(d.y.atan2(d.x), start, sweep)
        };
        let target = LabelTarget::Site(i);
        if v.hover.is_some_and(on) {
            hovered = Some(target);
        }
        if v.click.is_some_and(on) {
            clicked = Some(p.range);
        }
        draw_primer_arc(
            painter,
            c,
            r,
            start,
            sweep,
            p.forward,
            v.primer_color(p, false),
        );
    }

    // Labels with leader lines: features from their arc, sites from the backbone, primers
    // from their arrow. The hovered label is outlined and its line drawn in the accent colour.
    let labels: &[_] = if v.layers.labels { &l.labels } else { &[] };
    if v.layers.labels && l.hidden_labels > 0 {
        let n = l.hidden_labels;
        painter.text(
            rect.left_bottom() + vec2(8.0, -8.0),
            Align2::LEFT_BOTTOM,
            format!(
                "{n} label{} hidden for lack of room: enlarge the window or zoom in",
                if n == 1 { "" } else { "s" }
            ),
            FontId::proportional(11.0),
            t.text_muted,
        );
    }
    for label in labels {
        let site = matches!(label.target, LabelTarget::Site(_));
        let primer = v.primer_of(label.target);
        let from = match primer {
            Some(p) => c + (at(label.anchor) - c) * ((primer_radius(p) + 2.0) / radius),
            None if site => c + (at(label.anchor) - c) * ((radius + 4.0) / radius),
            None => at(label.anchor),
        };
        let to = at(label.at);
        let (align, offset) = if label.right {
            (Align2::LEFT_CENTER, vec2(4.0, 0.0))
        } else {
            (Align2::RIGHT_CENTER, vec2(-4.0, 0.0))
        };
        // A named enzyme is drawn larger and in the accent colour. There is no bold family
        // in the fonts we ship, and asking for one by name draws every glyph as a missing
        // character box - which is exactly what it did.
        let chosen_label = v.site_chosen(label.target);
        let base = if site { 11.0 } else { 12.0 };
        let font = FontId::proportional(if chosen_label { base + 2.0 } else { base });
        // Text that would run past the canvas's edge is shortened with an ellipsis (the
        // hover card names it in full).
        let room = if label.right {
            rect.right() - (to.x + offset.x) - 4.0
        } else {
            (to.x + offset.x) - rect.left() - 4.0
        };
        let galley = elided(painter, &label.text, font, room.max(24.0), t.text);
        let r = align.anchor_size(to + offset, galley.size());
        let hit = r.expand(1.0);
        if v.hover.is_some_and(|p| hit.contains(p)) {
            hovered = Some(label.target);
        }
        if v.click.is_some_and(|p| hit.contains(p)) {
            clicked = v.label_range(label.target);
        }
        let hot = hovered == Some(label.target);
        let (line, width) = if hot && site {
            (t.accent, 1.4)
        } else if hot {
            (t.text, 0.8)
        } else {
            (t.text_muted.gamma_multiply(0.6), 0.8)
        };
        painter.line_segment([from, to], Stroke::new(width, line));
        let chosen = chosen_label;
        let color = match primer {
            Some(p) if !hot => v.primer_color(p, false),
            // An enzyme named in the cut list stands out from the ones merely offered.
            _ if site && chosen && !hot => t.accent,
            _ if site && !hot => t.text_muted,
            _ => t.text,
        };
        if hot && site {
            painter.rect_filled(r.expand2(vec2(3.0, 1.0)), 3.0, t.canvas);
            painter.rect_stroke(
                r.expand2(vec2(3.0, 1.0)),
                3.0,
                Stroke::new(1.0, t.accent),
                egui::StrokeKind::Outside,
            );
            if let Some(p) = primer {
                let start = l.angle(p.range.start as f32);
                let sweep = (p.annealed as f32 / len * TAU).max(0.004);
                draw_primer_arc(
                    painter,
                    c,
                    primer_radius(p),
                    start,
                    sweep,
                    p.forward,
                    t.accent,
                );
            }
        }
        painter.galley_with_override_text_color(r.min, galley, color);
    }

    // Name and length in the middle; a long name is shortened to fit inside the backbone.
    let name_font = FontId::proportional(17.0);
    painter.text(
        c - vec2(0.0, 9.0),
        Align2::CENTER_CENTER,
        elide(painter, m.name(), &name_font, radius * 1.5),
        name_font,
        t.text,
    );
    painter.text(
        c + vec2(0.0, 11.0),
        Align2::CENTER_CENTER,
        format!("{} bp", m.len()),
        FontId::proportional(12.0),
        t.text_muted,
    );
    let empty = v.click.is_some() && clicked.is_none();
    (hovered, clicked, empty)
}

/// `text`, shortened with an ellipsis if it is wider than `max` pixels in `font`.
fn elide(painter: &egui::Painter, text: &str, font: &FontId, max: f32) -> String {
    let width = |s: &str| {
        painter
            .layout_no_wrap(s.to_string(), font.clone(), egui::Color32::WHITE)
            .size()
            .x
    };
    if width(text) <= max {
        return text.to_string();
    }
    let chars: Vec<char> = text.chars().collect();
    let (mut lo, mut hi) = (0, chars.len());
    // Longest prefix that fits with the ellipsis (binary search).
    while lo < hi {
        let mid = (lo + hi).div_ceil(2);
        let candidate: String = chars.iter().take(mid).collect::<String>() + "\u{2026}";
        if width(&candidate) <= max {
            lo = mid;
        } else {
            hi = mid - 1;
        }
    }
    chars
        .iter()
        .take(lo)
        .collect::<String>()
        .trim_end()
        .to_string()
        + "\u{2026}"
}

/// An ORF as a thin arc with a full arrowhead at its end (start for the reverse strand).
fn draw_orf_arc(
    painter: &egui::Painter,
    c: Pos2,
    r: f32,
    start: f32,
    sweep: f32,
    reverse: bool,
    color: egui::Color32,
) {
    painter.add(Shape::line(
        arc_points(c, r, start, sweep),
        Stroke::new(1.6, color),
    ));
    let head = (6.0 / r).min(sweep * 0.5);
    if head * r < 2.0 {
        return;
    }
    let (tip, base) = if reverse {
        (start, start + head)
    } else {
        (start + sweep, start + sweep - head)
    };
    let dir = |a: f32| vec2(a.cos(), a.sin());
    painter.add(Shape::convex_polygon(
        vec![
            c + dir(tip) * r,
            c + dir(base) * (r + 3.0),
            c + dir(base) * (r - 3.0),
        ],
        color,
        Stroke::NONE,
    ));
}

/// The hover card for an ORF, laid out like the other cards: position and frame, then the
/// protein's length, mass, pI and extinction coefficient.
fn orf_card(ui: &mut Ui, t: &Tokens, m: &Molecule, o: &snappiesnap_engine::bio::translate::Orf) {
    use snappiesnap_engine::bio::translate::{
        STANDARD, code, protein_properties, strand_residues, translate,
    };
    let len = m.len();
    let end = if o.range.end == 0 { len } else { o.range.end };
    let dir = if o.strand == Strand::Reverse {
        "Reverse"
    } else {
        "Forward"
    };
    let protein = code(STANDARD)
        .and_then(|c| strand_residues(m, o.range, o.strand).map(|nt| translate(&nt, c, true)));
    let summary = match &protein {
        Some(p) => {
            let pp = protein_properties(p);
            format!(
                "{} aa  \u{b7}  {:.1} kDa  \u{b7}  pI {:.2}  \u{b7}  \u{3b5}280 {}",
                pp.length,
                pp.mass / 1000.0,
                pp.isoelectric_point,
                pp.extinction_oxidised
            )
        }
        None => format!("{} aa", o.length),
    };
    card_header(
        ui,
        t,
        "Open reading frame",
        &format!("({}..{end})", o.range.start + 1),
        &format!("{dir}, frame {}", o.frame),
        &summary,
    );
    if let Some(p) = protein {
        ui.separator();
        let shown: String = p.iter().take(60).map(|&b| b as char).collect::<String>()
            + if p.len() > 60 { "\u{2026}" } else { "" };
        ui.label(egui::RichText::new(shown).monospace().color(t.text));
    }
}

/// A primer as a thin arc with a half arrowhead (one barb, outwards) at its 3' end.
fn draw_primer_arc(
    painter: &egui::Painter,
    c: Pos2,
    r: f32,
    start: f32,
    sweep: f32,
    forward: bool,
    color: egui::Color32,
) {
    painter.add(Shape::line(
        arc_points(c, r, start, sweep),
        Stroke::new(2.2, color),
    ));
    let head = (7.0 / r).min(sweep * 0.5).max(3.0 / r);
    let (tip, base) = if forward {
        (start + sweep, start + sweep - head)
    } else {
        (start, start + head)
    };
    let dir = |a: f32| vec2(a.cos(), a.sin());
    painter.add(Shape::convex_polygon(
        vec![
            c + dir(tip) * r,
            c + dir(base) * (r + 4.5),
            c + dir(base) * r,
        ],
        color,
        Stroke::NONE,
    ));
}

fn draw_linear(painter: &egui::Painter, rect: Rect, v: &View, selection: Option<Range>) -> Hits {
    let (t, m) = (v.t, v.m);
    let l = layout::linear(m, 10);
    let margin = 48.0;
    let left = rect.left() + margin;
    let width = (rect.width() - 2.0 * margin).max(10.0);
    let thick = 14.0;
    let lane_gap = thick + 26.0;
    let baseline = rect.top() + rect.height() * 0.62;
    let x = |u: f32| left + u * width;
    let primer_lanes = v
        .marks
        .iter()
        .filter_map(|k| match k.kind {
            MarkKind::Primer(j) => v.primers.get(j).map(|p| p.lane + 1),
            MarkKind::Site(_) => None,
        })
        .max()
        .unwrap_or(0);
    // Features sit above the primer lanes.
    let feature_base = baseline - 16.0 - 6.0 * primer_lanes as f32;

    if let Some(sel) = selection.filter(|s| !s.is_empty()) {
        let (a, b) = sel.segments(l.len);
        for seg in std::iter::once(a).chain(b) {
            let r = Rect::from_min_max(
                pos2(x(l.x(seg.start as f32)), rect.top() + 8.0),
                pos2(x(l.x(seg.end as f32)), baseline + 8.0),
            );
            painter.rect_filled(r, 0.0, t.selection);
        }
    }

    painter.line_segment(
        [pos2(left, baseline), pos2(left + width, baseline)],
        Stroke::new(2.0, t.backbone),
    );
    if v.layers.ruler {
        for (tick, u) in
            std::iter::once((layout::Tick { pos: 0, label: 1 }, 0.0)).chain(l.ticks.iter().copied())
        {
            let px = x(u);
            painter.line_segment(
                [pos2(px, baseline), pos2(px, baseline + 7.0)],
                Stroke::new(1.0, t.backbone),
            );
            painter.text(
                pos2(px, baseline + 18.0),
                Align2::CENTER_CENTER,
                tick.label.to_string(),
                FontId::proportional(10.0),
                t.text_muted,
            );
        }
    }
    painter.text(
        pos2(left + width / 2.0, rect.bottom() - 20.0),
        Align2::CENTER_CENTER,
        format!("{}  ·  {} bp", m.name(), m.len()),
        FontId::proportional(14.0),
        t.text,
    );

    let mut hovered = None;
    let mut clicked = None;
    let feature_lanes = if v.layers.features { l.lanes } else { 0 };
    if v.layers.features {
        for bar in &l.bars {
            let color = t.color_at(v.m, bar.feature);
            let y = feature_base - lane_gap * bar.lane as f32;
            let (x0, x1) = (x(bar.x0), x(bar.x1).max(x(bar.x0) + 2.0));
            let head = (thick * 0.9).min((x1 - x0) * 0.6);
            let (b0, b1) = match bar.arrow {
                Arrow::End => (x0, x1 - head),
                Arrow::Start => (x0 + head, x1),
                Arrow::None => (x0, x1),
            };
            if b1 > b0 {
                painter.rect_filled(
                    Rect::from_min_max(pos2(b0, y - thick / 2.0), pos2(b1, y + thick / 2.0)),
                    2.0,
                    color,
                );
            }
            let tri = match bar.arrow {
                Arrow::End => Some([
                    pos2(x1 - head, y - thick * 0.8),
                    pos2(x1, y),
                    pos2(x1 - head, y + thick * 0.8),
                ]),
                Arrow::Start => Some([
                    pos2(x0 + head, y - thick * 0.8),
                    pos2(x0, y),
                    pos2(x0 + head, y + thick * 0.8),
                ]),
                Arrow::None => None,
            };
            if let Some(tri) = tri {
                painter.add(Shape::convex_polygon(tri.to_vec(), color, Stroke::NONE));
            }
            let hit = Rect::from_min_max(pos2(x0, y - thick), pos2(x1, y + thick));
            if v.hover.is_some_and(|p| hit.contains(p)) {
                hovered = Some(LabelTarget::Feature(bar.feature));
            }
            if v.click.is_some_and(|p| hit.contains(p)) {
                clicked = feature_range(m, bar.feature);
            }
        }
        // Labels above their bars, skipping any that would overlap an earlier one in the lane.
        let mut taken: Vec<(usize, f32, f32)> = Vec::new();
        let labels: &[_] = if v.layers.labels { &l.labels } else { &[] };
        for label in labels {
            let y = feature_base - lane_gap * label.lane as f32 - thick / 2.0 - 9.0;
            let w = label.text.chars().count() as f32 * 6.5;
            let cx = x(label.x);
            let (a, b) = (cx - w / 2.0, cx + w / 2.0);
            if taken
                .iter()
                .any(|&(lane, ta, tb)| lane == label.lane && a < tb && ta < b)
            {
                continue;
            }
            taken.push((label.lane, a, b));
            painter.text(
                pos2(cx, y),
                Align2::CENTER_CENTER,
                &label.text,
                FontId::proportional(12.0),
                t.text,
            );
        }
    }

    // ORFs: thin arrows below the ruler, one lane each (hover targets after the marks).
    let orf_top = baseline + if v.layers.ruler { 32.0 } else { 12.0 };
    for (j, o) in v.orfs.iter().enumerate() {
        let y = orf_top + 5.0 * o.lane as f32;
        let target = LabelTarget::Site(v.marks.len() + j);
        let start = o.orf.range.start;
        let end = start + o.orf.range.len_on(l.len);
        let pieces = if end <= l.len {
            vec![(x(l.x(start as f32)), x(l.x(end as f32)))]
        } else {
            vec![
                (x(l.x(start as f32)), x(1.0)),
                (x(0.0), x(l.x((end - l.len) as f32))),
            ]
        };
        for &(a, b) in &pieces {
            let hit = Rect::from_min_max(pos2(a, y - 3.0), pos2(b.max(a + 2.0), y + 3.0));
            if v.hover.is_some_and(|q| hit.contains(q)) {
                hovered = Some(target);
            }
            if v.click.is_some_and(|q| hit.contains(q)) {
                clicked = Some(o.orf.range);
            }
        }
        let color = if hovered == Some(target) {
            t.accent
        } else {
            t.text_muted
        };
        for &(a, b) in &pieces {
            painter.line_segment(
                [pos2(a, y), pos2(b.max(a + 1.0), y)],
                Stroke::new(1.6, color),
            );
        }
        let reverse = o.orf.strand == Strand::Reverse;
        let (tip, dir) = if reverse {
            (pieces.first().map_or(left, |p| p.0), 1.0)
        } else {
            (pieces.last().map_or(left, |p| p.1), -1.0)
        };
        painter.add(Shape::convex_polygon(
            vec![
                pos2(tip, y),
                pos2(tip + dir * 5.0, y - 3.0),
                pos2(tip + dir * 5.0, y + 3.0),
            ],
            color,
            Stroke::NONE,
        ));
    }

    // Primers: thin arrows just above the backbone, one lane each.
    let primer_y = |p: &PrimerMark| baseline - 6.0 - 6.0 * p.lane as f32;
    let primer_segments = |p: &PrimerMark| -> Vec<(f32, f32)> {
        let end = p.range.start.saturating_add(p.annealed);
        if end <= l.len {
            vec![(x(l.x(p.range.start as f32)), x(l.x(end as f32)))]
        } else {
            vec![
                (x(l.x(p.range.start as f32)), x(1.0)),
                (x(0.0), x(l.x((end - l.len) as f32))),
            ]
        }
    };
    let draw_primer = |p: &PrimerMark, color: egui::Color32| {
        let y = primer_y(p);
        let segs = primer_segments(p);
        for (a, b) in &segs {
            painter.line_segment(
                [pos2(*a, y), pos2(b.max(a + 1.0), y)],
                Stroke::new(2.2, color),
            );
        }
        let (tip, dir) = if p.forward {
            (segs.last().map_or(left, |s| s.1), -1.0)
        } else {
            (segs.first().map_or(left, |s| s.0), 1.0)
        };
        painter.add(Shape::convex_polygon(
            vec![
                pos2(tip, y),
                pos2(tip + dir * 6.0, y - 4.5),
                pos2(tip + dir * 6.0, y),
            ],
            color,
            Stroke::NONE,
        ));
    };
    for (i, k) in v.marks.iter().enumerate() {
        let MarkKind::Primer(j) = k.kind else {
            continue;
        };
        let Some(p) = v.primers.get(j) else { continue };
        let y = primer_y(p);
        let target = LabelTarget::Site(i);
        for (a, b) in primer_segments(p) {
            let hit = Rect::from_min_max(pos2(a, y - 4.0), pos2(b.max(a + 2.0), y + 4.0));
            if v.hover.is_some_and(|q| hit.contains(q)) {
                hovered = Some(target);
            }
            if v.click.is_some_and(|q| hit.contains(q)) {
                clicked = Some(p.range);
            }
        }
        draw_primer(p, v.primer_color(p, false));
    }

    // Marks: rows of labels above the features, each with a line down to its position. Lines
    // are drawn first and labels over them on the canvas colour, so lines never cross text.
    let row_h = 14.0;
    let first_row = feature_base - lane_gap * feature_lanes as f32 - thick / 2.0 - 16.0;
    let mut rows: Vec<f32> = Vec::new(); // right edge of the last label per row
    let mut placed = Vec::new();
    for (index, k) in v.marks.iter().enumerate() {
        let px = x(l.x(k.pos as f32));
        let galley =
            painter.layout_no_wrap(k.text.clone(), FontId::proportional(11.0), t.text_muted);
        let w = galley.size().x;
        let lx = (px - w / 2.0).clamp(rect.left() + 2.0, (rect.right() - w - 2.0).max(rect.left()));
        let row = match rows.iter().position(|&end| end + 6.0 < lx) {
            Some(r) => r,
            None => {
                rows.push(f32::NEG_INFINITY);
                rows.len() - 1
            }
        };
        let y = first_row - row_h * row as f32;
        if y < rect.top() + 6.0 {
            continue;
        }
        if let Some(end) = rows.get_mut(row) {
            *end = lx + w;
        }
        let r = Rect::from_min_size(pos2(lx, y - row_h / 2.0), galley.size());
        let target = LabelTarget::Site(index);
        if v.hover.is_some_and(|p| r.expand(1.0).contains(p)) {
            hovered = Some(target);
        }
        if v.click.is_some_and(|p| r.expand(1.0).contains(p)) {
            clicked = Some(k.range);
        }
        let foot = match k.kind {
            MarkKind::Primer(j) => v
                .primers
                .get(j)
                .map_or(baseline - 2.0, |p| primer_y(p) - 2.0),
            MarkKind::Site(_) => baseline - 2.0,
        };
        placed.push((target, px, foot, r, k.text.clone()));
    }
    let label_color = |target: LabelTarget| match v.primer_of(target) {
        Some(p) => v.primer_color(p, false),
        None => t.text_muted,
    };
    for (target, px, foot, r, _) in &placed {
        if hovered == Some(*target) {
            continue;
        }
        painter.line_segment(
            [pos2(*px, *foot), pos2(*px, r.center().y)],
            Stroke::new(0.8, t.text_muted.gamma_multiply(0.5)),
        );
    }
    let font = FontId::proportional(11.0);
    for (target, _, _, r, text) in &placed {
        if hovered == Some(*target) {
            continue;
        }
        painter.rect_filled(r.expand2(vec2(2.0, 0.0)), 0.0, t.canvas);
        painter.text(
            r.min,
            Align2::LEFT_TOP,
            text,
            font.clone(),
            label_color(*target),
        );
    }
    // The hovered mark on top: accent line and outlined label (and the primer itself).
    if let Some((target, px, foot, r, text)) = placed.iter().find(|p| hovered == Some(p.0)) {
        if let Some(p) = v.primer_of(*target) {
            draw_primer(p, t.accent);
        }
        painter.line_segment(
            [pos2(*px, *foot), pos2(*px, r.max.y)],
            Stroke::new(1.4, t.accent),
        );
        let frame = r.expand2(vec2(3.0, 1.0));
        painter.rect_filled(frame, 3.0, t.canvas);
        painter.rect_stroke(
            frame,
            3.0,
            Stroke::new(1.0, t.accent),
            egui::StrokeKind::Outside,
        );
        painter.text(r.min, Align2::LEFT_TOP, text, font, t.text);
    } else if let Some(p) = hovered.and_then(|h| v.primer_of(h)) {
        draw_primer(p, t.accent);
    }
    let empty = v.click.is_some() && clicked.is_none();
    (hovered, clicked, empty)
}

/// The top of a hover card: a title, a muted position and a badge on the right; then a
/// muted summary line. Shared by the restriction site and primer cards.
fn card_header(ui: &mut Ui, t: &Tokens, title: &str, position: &str, badge: &str, summary: &str) {
    ui.set_min_width(230.0);
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(title).strong().color(t.text));
        ui.label(egui::RichText::new(position).color(t.text_muted));
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            egui::Frame::new()
                .stroke(Stroke::new(1.0, t.accent))
                .corner_radius(4)
                .inner_margin(egui::Margin::symmetric(6, 1))
                .show(ui, |ui| {
                    ui.label(egui::RichText::new(badge).small().color(t.text));
                });
        });
    });
    ui.label(egui::RichText::new(summary).color(t.text_muted));
}

/// A canvas for rows of bases in a hover card, with the same metrics in every card.
struct StrandCanvas {
    painter: egui::Painter,
    rect: Rect,
    left: f32,
    cw: f32,
    font: FontId,
}

/// Height of a row of bases, and the gap between the two strands.
const CARD_ROW: f32 = 17.0;
const CARD_GAP: f32 = 8.0;
/// Room for the 5'/3' labels at either side.
const CARD_SIDE: f32 = 22.0;

impl StrandCanvas {
    fn new(ui: &mut Ui, t: &Tokens, columns: usize, height: f32) -> Self {
        let font = FontId::monospace(13.0);
        let cw = ui
            .painter()
            .layout_no_wrap("G".to_string(), font.clone(), t.text)
            .size()
            .x;
        let size = vec2(columns as f32 * cw + 2.0 * CARD_SIDE, height);
        let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
        StrandCanvas {
            painter: ui.painter_at(rect.expand(2.0)),
            rect,
            left: rect.left() + CARD_SIDE,
            cw,
            font,
        }
    }

    /// Left edge of a column.
    fn x(&self, col: usize) -> f32 {
        self.left + col as f32 * self.cw
    }

    fn base(&self, col: usize, y: f32, ch: char, color: egui::Color32) {
        self.painter.text(
            pos2(self.x(col) + self.cw / 2.0, y),
            Align2::CENTER_CENTER,
            ch.to_string(),
            self.font.clone(),
            color,
        );
    }

    /// A strand with its end labels at the card edges; `bright` picks the highlighted columns.
    fn strand(
        &self,
        t: &Tokens,
        y: f32,
        text: &str,
        ends: (&str, &str),
        bright: impl Fn(usize) -> bool,
    ) {
        for (i, ch) in text.chars().enumerate() {
            self.base(i, y, ch, if bright(i) { t.text } else { t.text_muted });
        }
        self.end_label(
            t,
            pos2(self.rect.left() + 4.0, y),
            Align2::LEFT_CENTER,
            ends.0,
        );
        self.end_label(
            t,
            pos2(self.rect.right() - 4.0, y),
            Align2::RIGHT_CENTER,
            ends.1,
        );
    }

    fn end_label(&self, t: &Tokens, at: Pos2, align: Align2, text: &str) {
        self.painter
            .text(at, align, text, FontId::proportional(10.0), t.text_muted);
    }
}

/// The hover card for a restriction site: enzymes and cut, overhang, and both strands of the
/// sequence with the cut drawn through them.
fn site_card(ui: &mut Ui, t: &Tokens, d: &SiteDetail) {
    let sites = if d.count == 1 { "site" } else { "sites" };
    card_header(
        ui,
        t,
        &d.names,
        &format!("({})", d.position),
        &d.overhang_text(),
        &format!(
            "{} {sites}  ·  {} bp  ·  {}",
            d.count, d.specific, d.notation
        ),
    );
    let Some(w) = &d.window else {
        return;
    };
    ui.separator();
    let c = StrandCanvas::new(
        ui,
        t,
        w.top.chars().count(),
        2.0 * CARD_ROW + CARD_GAP + 16.0,
    );
    let top_y = c.rect.top() + 8.0 + CARD_ROW / 2.0;
    let bottom_y = top_y + CARD_ROW + CARD_GAP;
    let in_site = |i: usize| (w.site.0..w.site.1).contains(&i);
    c.strand(t, top_y, &w.top, ("5'", "3'"), in_site);
    c.strand(t, bottom_y, &w.bottom, ("3'", "5'"), in_site);
    // The cut: down through the top strand, along between the strands, down through the bottom.
    let (xt, xb) = (c.x(w.cut_top), c.x(w.cut_bottom));
    let mid = (top_y + bottom_y) / 2.0;
    c.painter.add(Shape::line(
        vec![
            pos2(xt, top_y - CARD_ROW / 2.0),
            pos2(xt, mid),
            pos2(xb, mid),
            pos2(xb, bottom_y + CARD_ROW / 2.0),
        ],
        Stroke::new(1.5, t.accent),
    ));
    let tri = 4.0;
    let y = top_y - CARD_ROW / 2.0 - 1.0;
    c.painter.add(Shape::convex_polygon(
        vec![
            pos2(xt - tri, y - tri * 1.4),
            pos2(xt + tri, y - tri * 1.4),
            pos2(xt, y),
        ],
        t.accent,
        Stroke::NONE,
    ));
    let y = bottom_y + CARD_ROW / 2.0 + 1.0;
    c.painter.add(Shape::convex_polygon(
        vec![
            pos2(xb, y),
            pos2(xb + tri, y + tri * 1.4),
            pos2(xb - tri, y + tri * 1.4),
        ],
        t.accent,
        Stroke::NONE,
    ));
}

/// The primer's structure: hairpin and self-dimer ΔG, then each warning in the error colour.
pub(crate) fn structure_lines(
    ui: &mut Ui,
    t: &Tokens,
    hairpin: Option<f64>,
    dimer: Option<f64>,
    checks: &[String],
) {
    let fmt = crate::enzymes::format_dg;
    let color = if checks.is_empty() {
        t.text_muted
    } else {
        t.error
    };
    ui.label(
        egui::RichText::new(format!(
            "Hairpin \u{394}G {}  \u{b7}  self-dimer \u{394}G {} kcal/mol",
            fmt(hairpin),
            fmt(dimer)
        ))
        .color(color),
    );
    for w in checks {
        ui.label(egui::RichText::new(format!("! {w}")).color(t.error));
    }
}

/// "Tm 52.2 °C", with the whole primer's Tm too when it has a tail or mismatches.
fn tm_text(p: &PrimerMark) -> Option<String> {
    let site = p.tm_site?;
    let mut text = format!("Tm {site:.1} \u{b0}C");
    if let Some(full) = p.tm_full
        && (full - site).abs() >= 0.05
    {
        let what = if p.mismatches > 0 {
            "paired part"
        } else {
            "annealed"
        };
        text = format!("Tm {site:.1} \u{b0}C {what}  \u{b7}  {full:.1} \u{b0}C whole primer");
    }
    Some(text)
}

/// The hover card for a primer binding site, laid out like the restriction site card: the
/// primer drawn against both strands, paired bases ticked, mismatches marked, the 5' tail
/// dimmed and an arrow at the 3' end, where extension starts.
fn primer_card(ui: &mut Ui, t: &Tokens, p: &PrimerMark, len: usize) {
    let end = if p.range.end == 0 { len } else { p.range.end };
    let mut summary = format!("{} nt  ·  {} annealed", p.seq.chars().count(), p.annealed);
    let tail = p.seq.chars().count().saturating_sub(p.annealed);
    if tail > 0 {
        summary.push_str(&format!("  ·  {tail} nt 5' tail"));
    }
    if p.mismatches > 0 {
        let s = if p.mismatches == 1 { "" } else { "es" };
        summary.push_str(&format!("  ·  {} mismatch{s}", p.mismatches));
    }
    let badge = match (p.forward, p.primary) {
        (true, true) => "Forward",
        (false, true) => "Reverse",
        (true, false) => "Forward, secondary",
        (false, false) => "Reverse, secondary",
    };
    card_header(
        ui,
        t,
        &p.name,
        &format!("({}..{end})", p.range.start + 1),
        badge,
        &summary,
    );
    if let Some(text) = tm_text(p) {
        ui.label(egui::RichText::new(text).color(t.text));
    }
    structure_lines(ui, t, p.hairpin_dg, p.dimer_dg, &p.checks);
    let Some(w) = &p.window else {
        return;
    };
    ui.separator();
    let tick = 6.0;
    let height = 3.0 * CARD_ROW + CARD_GAP + tick + 16.0;
    let c = StrandCanvas::new(ui, t, w.top.chars().count(), height);
    // Forward primers sit above the top strand, reverse primers below the bottom strand.
    let first = c.rect.top() + 8.0 + CARD_ROW / 2.0;
    let (primer_y, top_y) = if p.forward {
        (first, first + CARD_ROW + tick)
    } else {
        (0.0, first)
    };
    let bottom_y = top_y + CARD_ROW + CARD_GAP;
    let primer_y = if p.forward {
        primer_y
    } else {
        bottom_y + CARD_ROW + tick
    };
    let annealed = |i: usize| {
        w.primer
            .iter()
            .any(|&(col, _, pairing)| col == i && pairing != Pairing::Tail)
    };
    c.strand(t, top_y, &w.top, ("5'", "3'"), annealed);
    c.strand(t, bottom_y, &w.bottom, ("3'", "5'"), annealed);
    // The primer, with a tick to the strand it pairs with for every paired base.
    let partner_y = if p.forward {
        top_y - CARD_ROW / 2.0
    } else {
        bottom_y + CARD_ROW / 2.0
    };
    let near_y = if p.forward {
        primer_y + CARD_ROW / 2.0
    } else {
        primer_y - CARD_ROW / 2.0
    };
    for &(col, base, pairing) in &w.primer {
        let color = match pairing {
            Pairing::Paired => t.accent,
            Pairing::Mismatch => t.error,
            Pairing::Tail => t.text_muted,
        };
        c.base(col, primer_y, base, color);
        if pairing == Pairing::Paired {
            let x = c.x(col) + c.cw / 2.0;
            c.painter.line_segment(
                [pos2(x, near_y), pos2(x, partner_y)],
                Stroke::new(1.0, t.accent.gamma_multiply(0.7)),
            );
        }
    }
    // 5' end label (or an ellipsis for a long tail) and the 3' arrow.
    let (Some(&(first_col, ..)), Some(&(last_col, ..))) = (w.primer.first(), w.primer.last())
    else {
        return;
    };
    let (five, three_x, dir) = if p.forward {
        (c.x(first_col), c.x(last_col + 1), 1.0)
    } else {
        (c.x(last_col + 1), c.x(first_col), -1.0)
    };
    let five_text = if w.hidden_5prime > 0 { "…" } else { "5'" };
    let (five_at, five_align) = if p.forward {
        (pos2(five - 3.0, primer_y), Align2::RIGHT_CENTER)
    } else {
        (pos2(five + 3.0, primer_y), Align2::LEFT_CENTER)
    };
    c.end_label(t, five_at, five_align, five_text);
    let tip = pos2(three_x + dir * 7.0, primer_y);
    c.painter.add(Shape::convex_polygon(
        vec![
            tip,
            pos2(three_x + dir * 1.0, primer_y - 4.0),
            pos2(three_x + dir * 1.0, primer_y + 4.0),
        ],
        t.accent,
        Stroke::NONE,
    ));
}

/// The range from the start of a feature's first segment to the end of its last (5'→3').
fn feature_range(m: &Molecule, index: usize) -> Option<Range> {
    let f = m.features().get(index)?;
    let first = f.segments.first()?.range;
    let last = f.segments.last()?.range;
    if f.segments.len() == 1 {
        return Some(first);
    }
    // For reverse-strand features the segments run 5'→3' on the reverse strand (descending).
    let (a, b) = if first.start <= last.start {
        (first, last)
    } else {
        (last, first)
    };
    Some(Range::new(a.start, b.end))
}

/// Points along a clockwise arc, roughly every 3 px. Very short arcs get just their two ends:
/// near-coincident points make egui's line joins spike.
fn arc_points(c: Pos2, r: f32, start: f32, sweep: f32) -> Vec<Pos2> {
    let length = sweep.abs() * r;
    let steps = (length / 3.0).ceil().clamp(1.0, 720.0) as usize;
    (0..=steps)
        .map(|i| {
            let a = start + sweep * i as f32 / steps as f32;
            c + Vec2::new(a.cos(), a.sin()) * r
        })
        .collect()
}

/// Whether angle `a` lies within the clockwise sweep from `start`.
fn in_sweep(a: f32, start: f32, sweep: f32) -> bool {
    let d = (a - start).rem_euclid(TAU);
    d <= sweep
}

/// A small map of a construct: features and the ruler, nothing to click.
pub(crate) fn thumbnail_sized(
    ui: &mut Ui,
    t: &Tokens,
    m: &snappiesnap_engine::doc::Molecule,
    size: egui::Vec2,
    labels: bool,
) {
    let marks = crate::enzymes::Marks {
        sites: &[],
        primers: &[],
        orfs: &[],
    };
    let layers = Layers {
        enzymes: false,
        primers: false,
        orfs: false,
        labels,
        ..Layers::default()
    };
    ui.allocate_ui(size, |ui| {
        ui.set_min_size(size);
        let _ = show(ui, t, m, None, crate::state::MapStyle::Auto, &marks, layers);
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sweep_membership_wraps() {
        assert!(in_sweep(0.1, 0.0, 0.5));
        assert!(!in_sweep(0.6, 0.0, 0.5));
        assert!(in_sweep(0.05, TAU - 0.1, 0.3));
        assert!(!in_sweep(3.0, TAU - 0.1, 0.3));
    }

    #[test]
    fn arc_points_cover_the_sweep() {
        let p = arc_points(pos2(0.0, 0.0), 10.0, 0.0, TAU / 4.0);
        let first = p.first().copied().unwrap_or_default();
        let last = p.last().copied().unwrap_or_default();
        assert!((first - pos2(10.0, 0.0)).length() < 1e-3);
        assert!((last - pos2(0.0, 10.0)).length() < 1e-3);
    }

    /// Hover cards for a 5' overhang, a blunt cut and an outside cutter on the reverse strand
    /// draw without trouble. Set `SNAPPIESNAP_CARD_PNG=/path.png` to look at them.
    #[test]
    fn site_cards_draw() {
        use snappiesnap_engine::doc::{Alphabet, Seq};
        use snappiesnap_engine::enzymes::EnzymeDb;
        let s = format!("TTTGAATTCTTTCCCGGGTTT{}GAGACCAAAA", "T".repeat(20));
        let m = Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Linear).unwrap();
        let sites = crate::enzymes::map_sites(
            &EnzymeDb::default(),
            &m,
            &snappiesnap_engine::enzymes::SiteQuery::default(),
        );
        let pick = |n: &str| {
            sites
                .iter()
                .find(|s| s.detail.names.split(", ").any(|x| x == n))
                .unwrap()
        };
        let cards = [pick("EcoRI"), pick("SmaI"), pick("BsaI")];
        let mut pm = Molecule::new(
            "p",
            Seq::new(
                Alphabet::Dna,
                "AAAAAAAAAAACGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAAAAAAAAAA",
            )
            .unwrap(),
            Topology::Linear,
        )
        .unwrap();
        for (n, p) in [
            ("fwd-tail", "gaattcCGTTGCATGCAAACTT"),
            ("rev", "GTAAAACGACGGCCAGTGCCA"),
        ] {
            pm.add_primer(n, Seq::new(Alphabet::Dna, p).unwrap(), "")
                .unwrap();
        }
        let primer_marks = crate::enzymes::primer_marks(&pm, &Default::default());
        assert_eq!(primer_marks.len(), 2);
        assert_eq!(cards[1].detail.overhang_text(), "Blunt");
        assert!(cards[2].detail.window.is_some());
        let t = Tokens::dark();
        let image = crate::offscreen::render(
            460,
            760,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::popup(ui.style()).show(ui, |ui| {
                    for c in &cards {
                        site_card(ui, &t, &c.detail);
                        ui.add_space(12.0);
                    }
                    for p in &primer_marks {
                        primer_card(ui, &t, p, pm.len());
                        ui.add_space(12.0);
                    }
                });
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_CARD_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }
}
