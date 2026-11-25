//! The cloning schematic: what the reaction is about to do, drawn as lines under the map.
//!
//! The vector is the canvas foreground ([`Tokens::backbone`] — black on a light theme) and
//! everything that goes into it is red ([`Tokens::insert`]), which is how these diagrams are
//! drawn at the bench. Lengths are to scale whenever they are known, so a 762 bp insert going
//! into a 3.2 kb vector looks like one; a piece whose length isn't known yet is drawn with a
//! dashed outline rather than a guessed size.

use eframe::egui::{self, Align2, Color32, Rect, RichText, Stroke, Ui, pos2, vec2};

use snappiesnap_engine::bio::cloning::Fit;

use crate::state::{CloningTechnique, Experiment, ExperimentStep, SideUse};
use crate::theme::Tokens;

/// One piece of a row: how long it is, what it is called, and whether it is an insert.
#[derive(Debug, Clone, PartialEq)]
pub struct Piece {
    /// Length in bp; `None` when the construct hasn't been read yet.
    pub length: Option<u64>,
    pub label: String,
    /// Red (a piece going in) rather than the vector's black.
    pub insert: bool,
    /// Cut out and thrown away rather than kept: the stuffer between two sites of a vector,
    /// the flanks a cut insert loses. Drawn apart from the piece that is kept.
    pub removed: bool,
}

impl Piece {
    /// A piece that is kept: the vector's backbone, or what goes into it.
    fn kept(length: Option<u64>, label: impl Into<String>, insert: bool) -> Piece {
        Piece {
            length,
            label: label.into(),
            insert,
            removed: false,
        }
    }
}

/// How a row is drawn: a straight line, or a circle, whole or cut open. A plasmid is a circle
/// at the bench, and the point of a step is often that the circle opens, so the schematic
/// draws it as one.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Ring {
    /// A straight line, left to right.
    #[default]
    None,
    /// A closed circle: a plasmid, or the product.
    Closed,
    /// A circle with a gap where it was cut. Pieces marked `removed` are drawn in the gap,
    /// between the two cut ends, which is where they came from.
    Open,
}

/// A labelled row of the schematic: the vector, what goes into it, or the product.
#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    pub title: String,
    pub pieces: Vec<Piece>,
    /// The ends of the row, named by the enzymes or sites that make them.
    pub left: Option<String>,
    pub right: Option<String>,
    /// Closed on itself (a plasmid): drawn with its ends joined and said in words. `None`
    /// when the construct hasn't been read, so nothing is claimed either way.
    pub circular: Option<bool>,
    /// Drawn as a circle rather than a line. Only set when the lengths are known, so a ring
    /// never implies a size that hasn't been read.
    pub ring: Ring,
}

/// The whole plan: the rows, a sentence saying what the reaction does, and anything the
/// user should know before running it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Plan {
    pub rows: Vec<Row>,
    pub caption: String,
    /// What might not work, or might work in a way that wasn't intended. Shown under the
    /// schematic, in the warning colour.
    pub warnings: Vec<String>,
}

/// What a side contributes, in words: the enzymes it is cut with, or how else it is used.
fn ends(
    enzymes: &[String],
    technique: CloningTechnique,
    using: SideUse,
) -> (Option<String>, Option<String>) {
    match (technique, using) {
        (CloningTechnique::GatewayLr, _) => (Some("attL1".into()), Some("attL2".into())),
        (CloningTechnique::GatewayBp, _) => (Some("attB1".into()), Some("attB2".into())),
        (_, SideUse::Pcr) => (Some("PCR end".into()), Some("PCR end".into())),
        _ => match enzymes {
            [] => (None, None),
            [one] => (Some(one.clone()), Some(one.clone())),
            [first, .., last] => (Some(first.clone()), Some(last.clone())),
        },
    }
}

/// What has been read about the constructs so far. Every field is `None` until it is known,
/// so the schematic can draw what it knows and say nothing about the rest.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Measured {
    pub vector_len: Option<u64>,
    pub insert_len: Option<u64>,
    /// The previewed product's length.
    pub product_len: Option<u64>,
    pub vector_circular: Option<bool>,
    pub insert_circular: Option<bool>,
    /// What the constructs are called, when they are open tabs rather than files.
    pub vector_name: Option<String>,
    pub insert_name: Option<String>,
    /// Whether the piece goes into the backbone, and whether its direction is settled.
    /// `None` when that hasn't been worked out (no digest, or a technique that doesn't cut).
    pub fits: Option<Fit>,
    /// What cutting each side actually gives: the fragment lengths in order along the
    /// molecule, from the digest itself rather than from a guess. Empty until the construct
    /// has been read and the enzymes named.
    pub vector_fragments: Vec<u64>,
    pub insert_fragments: Vec<u64>,
}

/// The plan for the experiment as it stands.
pub fn plan(
    exp: &Experiment,
    m: &Measured,
    vector_enzymes: &[String],
    insert_enzymes: &[String],
) -> Plan {
    let Measured {
        vector_len,
        insert_len,
        product_len,
        vector_circular,
        insert_circular,
        ref vector_fragments,
        ref insert_fragments,
        ref vector_name,
        ref insert_name,
        fits,
    } = *m;
    let t = exp.technique;
    let (vl, vr) = ends(vector_enzymes, t, exp.vector.using);
    // An insert with no enzymes of its own is cut with the vector's.
    let cuts: &[String] = if insert_enzymes.is_empty() {
        vector_enzymes
    } else {
        insert_enzymes
    };
    let (il, ir) = ends(cuts, t, exp.insert.using);
    let cuts = cuts.to_vec();

    // What to call each side: its file, or the name of the tab it is open in.
    let named = |file: &str, open: &Option<String>, fallback: &str| -> String {
        for candidate in [file, open.as_deref().unwrap_or("")] {
            if !candidate.trim().is_empty() {
                return candidate.to_string();
            }
        }
        fallback.to_string()
    };
    let vector_label = named(exp.vector.file_name(), vector_name, "vector");

    let mut rows = Vec::new();
    let opened = match t {
        CloningTechnique::Mutagenesis => "the plasmid, copied from the primers",
        CloningTechnique::GatewayLr | CloningTechnique::GatewayBp => "the backbone, att to att",
        _ if exp.vector.using == SideUse::Whole && t.is_assembly() => "the whole vector",
        _ => "the backbone, cut to cut",
    };
    rows.push(Row {
        title: "Vector".to_string(),
        pieces: vec![Piece::kept(vector_len, vector_label, false)],
        left: vl.clone(),
        right: vr.clone(),
        // The vector is cut open for the reaction, so the row shows what it was.
        circular: vector_circular,
        ring: ring_for(vector_circular, vector_len),
    });

    if t != CloningTechnique::Mutagenesis {
        let mut pieces = vec![Piece::kept(
            insert_len,
            named(exp.insert.file_name(), insert_name, "insert"),
            true,
        )];
        for (i, extra) in exp.extra.iter().filter(|s| s.is_set()).enumerate() {
            pieces.push(Piece::kept(
                None,
                {
                    let n = extra.file_name().to_string();
                    if n.is_empty() {
                        format!("part {}", i + 3)
                    } else {
                        n
                    }
                },
                true,
            ));
        }
        // One piece only can be a circle; several parts side by side are lines.
        let ring = if pieces.len() == 1 {
            ring_for(insert_circular, insert_len)
        } else {
            Ring::None
        };
        rows.push(Row {
            title: "Insert".to_string(),
            pieces,
            left: il,
            right: ir,
            circular: insert_circular,
            ring,
        });
    }

    // The product: the vector's backbone with everything else inside it, closed into a circle.
    // Where the digest has been run, the backbone is the fragment the reaction keeps rather
    // than the whole vector, so the two arcs are in the proportion the product really has.
    let backbone_len = vector_fragments.iter().max().copied().or(vector_len);
    let going_in = insert_fragments.iter().max().copied().or(insert_len);
    // A piece is named with its size only where that size came from the digest; elsewhere
    // it is the whole construct's length, which is not what ends up in the product.
    let size_of = |len: Option<u64>, measured: bool| match (len, measured) {
        (Some(n), true) => format!("{n} bp"),
        _ => String::new(),
    };
    let mut product = vec![Piece::kept(
        backbone_len,
        size_of(backbone_len, !vector_fragments.is_empty()),
        false,
    )];
    if t != CloningTechnique::Mutagenesis {
        product.push(Piece::kept(
            going_in,
            size_of(going_in, !insert_fragments.is_empty()),
            true,
        ));
    }
    // The product is drawn as a circle only when its pieces are known and they add up to the
    // product the reaction reports (allowing for the overhangs a digest counts at both ends).
    // Where they don't -- a Gateway backbone is not the whole vector, for instance -- the
    // arcs would be in the wrong proportion, so the row stays a line.
    let sum: u64 = product.iter().filter_map(|p| p.length).sum();
    let adds_up = match product_len {
        Some(n) if n > 0 => sum.abs_diff(n) * 20 <= n,
        _ => false,
    };
    let product_ring = if product.iter().all(|p| p.length.is_some()) && adds_up {
        Ring::Closed
    } else {
        Ring::None
    };
    rows.push(Row {
        title: match product_len {
            Some(n) => format!("Product · {n} bp"),
            None => "Product".to_string(),
        },
        pieces: product,
        left: None,
        right: None,
        // Every reaction here gives a plasmid, except an assembly told not to close.
        circular: Some(true),
        ring: product_ring,
    });

    // Only this step's work: what is being prepared and what that gives. The whole reaction
    // belongs to the last step, where the product is.
    let (rows, caption) = match exp.step {
        ExperimentStep::Vector => (
            vec![
                rows.first().cloned().unwrap_or_else(blank),
                prepared(
                    prepared_title(t, exp.vector.using, vector_enzymes, true),
                    vl,
                    vr,
                    false,
                    vector_fragments,
                    vector_circular,
                ),
            ],
            format!("{}: this step prepares {opened}.", t.label()),
        ),
        ExperimentStep::Insert if t == CloningTechnique::Mutagenesis => (
            vec![rows.first().cloned().unwrap_or_else(blank)],
            format!(
                "{}: this step takes the two primers that copy the plasmid; what differs in \
                 their 5' parts becomes the change.",
                t.label()
            ),
        ),
        ExperimentStep::Insert => {
            let insert_row = rows.get(1).cloned().unwrap_or_else(blank);
            let (il, ir) = (insert_row.left.clone(), insert_row.right.clone());
            (
                vec![
                    insert_row,
                    prepared(
                        prepared_title(t, exp.insert.using, &cuts, false),
                        il,
                        ir,
                        true,
                        insert_fragments,
                        insert_circular,
                    ),
                ],
                format!(
                    "{}: this step prepares the piece that goes in, and says whether it fits.",
                    t.label()
                ),
            )
        }
        _ => (rows, format!("{}: {opened}.", t.label())),
    };

    // What the schematic can't draw: whether the piece goes in at all, and whether the way
    // round is settled. Only once the step is about the piece that goes in.
    let mut warnings = Vec::new();
    // More than one cut-out piece means the reaction's rule -- keep the largest fragment --
    // is a choice, and it may not be the backbone that was meant.
    if matches!(exp.step, ExperimentStep::Vector | ExperimentStep::Product)
        && vector_fragments.len() > 2
    {
        warnings.push(format!(
            "These enzymes cut the vector into {} pieces; the reaction keeps the largest. \
             Check that is the backbone you want.",
            vector_fragments.len()
        ));
    }
    if exp.step != ExperimentStep::Vector {
        match fits {
            Some(Fit::No) => warnings.push(
                "The ends don't match: as cut, this piece won't go into the backbone. Check \
                 the enzymes on both sides."
                    .to_string(),
            ),
            Some(Fit::EitherWay) => warnings.push(
                "Both ends of the piece take the same cut, so it can go in either way round: \
                 expect both orientations and screen the colonies for the one you want."
                    .to_string(),
            ),
            _ => {}
        }
    }

    Plan {
        rows,
        caption,
        warnings,
    }
}

/// A row for a side that hasn't been chosen yet.
fn blank() -> Row {
    Row {
        title: String::new(),
        pieces: Vec::new(),
        left: None,
        right: None,
        circular: None,
        ring: Ring::None,
    }
}

/// A construct drawn as a circle only when it is known to be one and its size has been read.
fn ring_for(circular: Option<bool>, length: Option<u64>) -> Ring {
    match (circular, length) {
        (Some(true), Some(_)) => Ring::Closed,
        _ => Ring::None,
    }
}

/// What this step does to a side, in a few words.
fn prepared_title(
    t: CloningTechnique,
    using: SideUse,
    enzymes: &[String],
    is_vector: bool,
) -> String {
    match (t, using) {
        (CloningTechnique::GatewayLr | CloningTechnique::GatewayBp, _) => {
            "Between the att sites".to_string()
        }
        (CloningTechnique::Mutagenesis, _) => "Copied from the primers".to_string(),
        (
            CloningTechnique::TopoTa
            | CloningTechnique::TopoBlunt
            | CloningTechnique::TopoDirectional,
            _,
        ) if is_vector => "Opened".to_string(),
        (_, SideUse::Pcr) => "The PCR product".to_string(),
        (_, SideUse::Whole) if t.is_assembly() && is_vector => "Used whole".to_string(),
        _ if enzymes.is_empty() => "Cut (name the enzymes)".to_string(),
        _ => format!("Cut with {}", enzymes.join(" + ")),
    }
}

/// The pieces this step makes. With the digest run, these are the real fragments: the one the
/// reaction uses, and the ones cut away from it. Without it (no enzymes yet, nothing read, or
/// a technique that doesn't cut) the ends are known but the lengths are not, so the piece is
/// drawn as an outline at no particular size.
fn prepared(
    title: String,
    left: Option<String>,
    right: Option<String>,
    insert: bool,
    fragments: &[u64],
    was_circular: Option<bool>,
) -> Row {
    if fragments.is_empty() {
        return Row {
            title,
            pieces: vec![Piece::kept(None, "", insert)],
            left,
            right,
            // Cutting a plasmid opens it; a piece made this way is linear.
            circular: Some(false),
            ring: Ring::None,
        };
    }
    // The reaction keeps the largest fragment (`bio::cloning::restriction_clone` takes the
    // vector's largest; an insert's largest is the piece a user cuts out, the rest flanks).
    let keep = fragments
        .iter()
        .enumerate()
        .max_by_key(|(_, n)| **n)
        .map(|(i, _)| i);
    let pieces = fragments
        .iter()
        .enumerate()
        .map(|(i, &n)| Piece {
            length: Some(n),
            label: format!("{n} bp"),
            insert: insert && Some(i) == keep,
            removed: Some(i) != keep,
        })
        .collect();
    Row {
        title,
        pieces,
        left,
        right,
        circular: Some(false),
        // A circle that was cut is drawn as a circle with a gap: what the cut did is the
        // point of the step. A linear construct stays a line.
        ring: if was_circular == Some(true) {
            Ring::Open
        } else {
            Ring::None
        },
    }
}

/// Height of a line row, of a circle row, and the width the schematic needs.
const ROW_H: f32 = 46.0;
const RING_H: f32 = 104.0;
/// Radius of the largest circle drawn; smaller constructs are drawn smaller, to scale.
const MAX_R: f32 = 32.0;
/// The gap left where a circle was cut, in degrees, at the top of the circle.
const GAP_DEG: f32 = 52.0;
pub const MIN_WIDTH: f32 = 260.0;

/// How tall a row is: a circle needs more room than a line.
fn row_height(row: &Row) -> f32 {
    if row.ring == Ring::None {
        ROW_H
    } else {
        RING_H
    }
}

/// The space the plan needs, so a caller can leave room for it.
pub fn height(plan: &Plan) -> f32 {
    plan.rows.iter().map(row_height).sum::<f32>() + 8.0
}

/// Draw the plan. Returns the space it used.
pub fn show(ui: &mut Ui, t: &Tokens, plan: &Plan, width: f32) -> egui::Response {
    ui.label(RichText::new(&plan.caption).small().color(t.text_muted));
    let width = width.max(MIN_WIDTH);
    let (rect, response) = ui.allocate_exact_size(vec2(width, height(plan)), egui::Sense::hover());
    let painter = ui.painter_at(rect);
    let font = egui::FontId::proportional(10.0);
    // The longest row sets the scale, so rows are comparable with each other.
    let longest = plan
        .rows
        .iter()
        .map(|r| r.pieces.iter().filter_map(|p| p.length).sum::<u64>())
        .max()
        .unwrap_or(0)
        .max(1);
    // The column fits the longest row name, so a name like "Cut with EcoRI + BamHI" can't
    // run over the track; past 45% of the width it is shortened instead.
    let longest_title = plan
        .rows
        .iter()
        .map(|r| {
            painter
                .layout_no_wrap(r.title.clone(), font.clone(), t.text_muted)
                .size()
                .x
        })
        .fold(0.0_f32, f32::max);
    let label_w = (longest_title + 8.0).clamp(60.0, width * 0.45);
    // Room at the right for "circular"/"linear".
    let track = (rect.width() - label_w - 60.0).max(40.0);

    let mut top = rect.top() + 4.0;
    for row in &plan.rows {
        let h = row_height(row);
        let mid = if row.ring == Ring::None {
            // The line sits below the end labels and above the piece names, so the three
            // never land on top of each other.
            top + ROW_H * 0.55
        } else {
            top + h / 2.0
        };
        painter.text(
            pos2(rect.left(), mid),
            Align2::LEFT_CENTER,
            fit(&painter, &font, &row.title, label_w - 8.0),
            font.clone(),
            t.text_muted,
        );
        let colour = |p: &Piece| {
            if p.removed {
                t.text_muted
            } else if p.insert {
                t.insert
            } else {
                t.backbone
            }
        };
        if row.ring != Ring::None {
            draw_ring(
                &painter,
                t,
                &font,
                row,
                pos2(rect.left() + label_w + MAX_R, mid),
                longest,
            );
        } else {
            draw_line_row(
                &painter,
                t,
                &font,
                row,
                rect.left() + label_w,
                mid,
                track,
                longest,
                &colour,
            );
        }
        // Circular or linear, in words: the lines alone can't say it, and a reaction that
        // needs a linear insert or an opened vector depends on it.
        if let Some(circular) = row.circular {
            painter.text(
                pos2(rect.right(), mid),
                Align2::RIGHT_CENTER,
                if circular { "circular" } else { "linear" },
                font.clone(),
                t.text_muted,
            );
        }
        top += h;
    }
    for warning in &plan.warnings {
        ui.label(RichText::new(warning).small().color(t.warning));
    }
    response
}

/// A row drawn as a circle: the pieces that are kept as arcs round it, in proportion, and
/// anything cut away drawn in the gap between the two cut ends, where it came from.
fn draw_ring(
    painter: &egui::Painter,
    t: &Tokens,
    font: &egui::FontId,
    row: &Row,
    centre: egui::Pos2,
    longest: u64,
) {
    let kept: Vec<&Piece> = row.pieces.iter().filter(|p| !p.removed).collect();
    let total: u64 = kept.iter().filter_map(|p| p.length).sum::<u64>().max(1);
    // To scale with the other rows, but never so small that it stops reading as a circle.
    let r = (MAX_R * (total as f32 / longest as f32).sqrt()).clamp(12.0, MAX_R);
    let gap = if row.ring == Ring::Open { GAP_DEG } else { 0.0 };
    let span = 360.0 - gap;
    let point = |deg: f32| {
        let a = deg.to_radians();
        pos2(centre.x + r * a.cos(), centre.y + r * a.sin())
    };
    // The gap sits at the top, so the piece cut out can be drawn in it.
    let mut at = -90.0 + gap / 2.0;
    for piece in &kept {
        let share = piece.length.unwrap_or(0) as f32 / total as f32;
        let sweep = span * share;
        let steps = ((sweep / 4.0).ceil() as usize).max(2);
        let points: Vec<egui::Pos2> = (0..=steps)
            .map(|k| point(at + sweep * k as f32 / steps as f32))
            .collect();
        let colour = if piece.insert { t.insert } else { t.backbone };
        painter.add(egui::Shape::line(points, Stroke::new(3.0, colour)));
        at += sweep;
    }
    // What is cut away: short straight pieces in the gap, to the same scale as the arcs, so
    // a 40 bp stuffer looks like the sliver of the plasmid it is.
    let per_bp = (std::f32::consts::TAU * r * (span / 360.0)) / total as f32;
    let mut y = centre.y - r - 7.0;
    for piece in row.pieces.iter().filter(|p| p.removed) {
        // Never thinner than a few points: a 10 bp stuffer is a speck at this scale,
        // and a speck reads as a smudge rather than as a piece of DNA.
        let w = (piece.length.unwrap_or(0) as f32 * per_bp).clamp(7.0, 2.0 * r);
        painter.rect_filled(
            Rect::from_min_size(pos2(centre.x - w / 2.0, y - 1.5), vec2(w, 3.0)),
            1.0,
            t.text_muted,
        );
        if !piece.label.is_empty() {
            painter.text(
                pos2(centre.x + w / 2.0 + 5.0, y),
                Align2::LEFT_CENTER,
                &piece.label,
                font.clone(),
                t.text_muted,
            );
        }
        y -= 9.0;
    }
    // The names of the pieces kept, under the circle.
    let names: Vec<&str> = kept
        .iter()
        .map(|p| p.label.as_str())
        .filter(|l| !l.is_empty())
        .collect();
    if !names.is_empty() {
        painter.text(
            pos2(centre.x, centre.y + r + 5.0),
            Align2::CENTER_TOP,
            names.join(" + "),
            font.clone(),
            t.text_muted,
        );
    }
}

/// A row drawn as a line: the pieces end to end, with a break before anything cut away.
#[allow(clippy::too_many_arguments)]
fn draw_line_row(
    painter: &egui::Painter,
    t: &Tokens,
    font: &egui::FontId,
    row: &Row,
    start: f32,
    mid: f32,
    track: f32,
    longest: u64,
    colour: &dyn Fn(&Piece) -> Color32,
) {
    let known: u64 = row.pieces.iter().filter_map(|p| p.length).sum();
    // Pieces of unknown length share what is left, so nothing is drawn to a made-up size.
    let unknown = row.pieces.iter().filter(|p| p.length.is_none()).count();
    let used = if known == 0 {
        0.0
    } else {
        track * (known as f32 / longest as f32)
    };
    let spare = (track - used).max(0.0);
    let each_unknown = if unknown == 0 {
        0.0
    } else {
        (spare / unknown as f32).min(60.0)
    };
    let mut x = start;
    for (i, piece) in row.pieces.iter().enumerate() {
        // A break before a piece that is cut away, and before the piece after it: the gaps
        // are where the cuts are.
        if i > 0 && (piece.removed || row.pieces.get(i - 1).is_some_and(|p| p.removed)) {
            x += 5.0;
        }
        let w = match piece.length {
            Some(n) => (track * (n as f32 / longest as f32)).max(2.0),
            None => each_unknown.max(24.0),
        };
        let c = colour(piece);
        let line = Rect::from_min_size(pos2(x, mid - 2.0), vec2(w, 4.0));
        if piece.length.is_some() {
            painter.rect_filled(line, 1.0, c);
        } else {
            // Not read yet: an outline, so its width is clearly not a measurement.
            painter.rect_stroke(line, 1.0, Stroke::new(1.0, c), egui::StrokeKind::Middle);
        }
        if !piece.label.is_empty() {
            painter.text(
                pos2(x + w / 2.0, mid + 6.0),
                Align2::CENTER_TOP,
                &piece.label,
                font.clone(),
                c,
            );
        }
        x += w;
    }
    // Ends, named by what makes them.
    let tick = |at: f32, text: &Option<String>, align: Align2| {
        if let Some(name) = text {
            painter.line_segment(
                [pos2(at, mid - 7.0), pos2(at, mid + 7.0)],
                Stroke::new(1.0, t.text_muted),
            );
            painter.text(pos2(at, mid - 9.0), align, name, font.clone(), t.text_muted);
        }
    };
    tick(start, &row.left, Align2::LEFT_BOTTOM);
    // On a short piece the two names would print over each other, so the right-hand one
    // moves outside the line instead of being drawn on top of the left-hand one.
    let crowded = x - start < 90.0;
    tick(
        x,
        &row.right,
        if crowded {
            Align2::LEFT_BOTTOM
        } else {
            Align2::RIGHT_BOTTOM
        },
    );
    if row.circular == Some(true) {
        // A closed plasmid whose size isn't known: the ends meet, drawn as a hoop over the
        // whole row.
        let (a, b) = (start, x);
        // Under the line, where the end labels are not.
        let under = mid + 7.0;
        for end in [a, b] {
            painter.line_segment(
                [pos2(end, mid), pos2(end, under)],
                Stroke::new(1.0, t.backbone),
            );
        }
        painter.line_segment(
            [pos2(a, under), pos2(b, under)],
            Stroke::new(1.0, t.backbone),
        );
    }
}

/// `text`, shortened with an ellipsis until it fits `max` points.
fn fit(painter: &egui::Painter, font: &egui::FontId, text: &str, max: f32) -> String {
    let width = |s: &str| {
        painter
            .layout_no_wrap(s.to_string(), font.clone(), Color32::WHITE)
            .size()
            .x
    };
    if width(text) <= max {
        return text.to_string();
    }
    let mut cut: Vec<char> = text.chars().collect();
    while !cut.is_empty() {
        cut.pop();
        let candidate: String = cut.iter().collect::<String>() + "\u{2026}";
        if width(&candidate) <= max {
            return candidate;
        }
    }
    String::new()
}

/// A legend line, so the colours need no explaining. "Cut away" is named only when the plan
/// actually has a piece that is thrown away.
pub fn legend(ui: &mut Ui, t: &Tokens, plan: &Plan) {
    let any_removed = plan.rows.iter().any(|r| r.pieces.iter().any(|p| p.removed));
    ui.horizontal(|ui| {
        let dot = |ui: &mut Ui, colour: Color32, text: &str| {
            let (rect, _) = ui.allocate_exact_size(vec2(14.0, 4.0), egui::Sense::hover());
            ui.painter().rect_filled(rect, 1.0, colour);
            ui.label(RichText::new(text).small().color(t.text_muted));
        };
        dot(ui, t.backbone, "vector");
        dot(ui, t.insert, "what goes in");
        if any_removed {
            dot(ui, t.text_muted, "cut away");
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::{CloningSide, ExperimentKind};

    fn experiment(technique: CloningTechnique) -> Experiment {
        let mut e = Experiment::new(ExperimentKind::Cloning, "Cloning", []);
        e.technique = technique;
        e.vector = CloningSide {
            file: "pSample.gb".into(),
            enzymes: "EcoRI, BamHI".into(),
            ..Default::default()
        };
        e.insert = CloningSide {
            file: "demoInsert.fa".into(),
            ..Default::default()
        };
        e
    }

    /// The plan names the ends by the enzymes that make them, and an insert with no enzymes
    /// of its own is cut with the vector's, exactly as the reaction does it.
    #[test]
    fn the_plan_follows_the_reaction() {
        let e = experiment(CloningTechnique::Restriction);
        let p = plan(
            &e,
            &Measured {
                vector_len: Some(3200),
                insert_len: Some(762),
                product_len: Some(3950),
                vector_circular: Some(true),
                insert_circular: Some(false),
                ..Measured::default()
            },
            &["EcoRI".into(), "BamHI".into()],
            &[],
        );
        assert_eq!(p.rows.len(), 3);
        assert_eq!(p.rows[0].left.as_deref(), Some("EcoRI"));
        assert_eq!(p.rows[0].right.as_deref(), Some("BamHI"));
        // The insert takes the vector's enzymes when it names none.
        assert_eq!(p.rows[1].left.as_deref(), Some("EcoRI"));
        assert_eq!(p.rows[1].right.as_deref(), Some("BamHI"));
        assert_eq!(p.rows[1].pieces[0].label, "demoInsert.fa");
        assert!(p.rows[1].pieces[0].insert, "the insert is the red line");
        assert!(!p.rows[0].pieces[0].insert, "the vector is the black line");
        assert_eq!(p.rows[2].title, "Product · 3950 bp");
        assert_eq!(
            p.rows[2].circular,
            Some(true),
            "the product closes into a plasmid"
        );
        // A plasmid cut open and a linear insert: the schematic says which is which.
        assert_eq!(p.rows[0].circular, Some(true));
        assert_eq!(p.rows[1].circular, Some(false));
    }

    /// Each step's plan shows that step's work and no more: the whole reaction belongs to
    /// the last step, where the product is.
    #[test]
    fn the_plan_shows_only_this_steps_work() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        let measured = Measured {
            vector_len: Some(3200),
            insert_len: Some(762),
            product_len: Some(3950),
            vector_circular: Some(true),
            insert_circular: Some(false),
            ..Measured::default()
        };
        let enzymes = ["EcoRI".to_string(), "BamHI".to_string()];

        // Step 2: the vector as it stands, and what cutting it gives.
        e.step = ExperimentStep::Vector;
        let p = plan(&e, &measured, &enzymes, &[]);
        let titles: Vec<&str> = p.rows.iter().map(|r| r.title.as_str()).collect();
        assert_eq!(titles, vec!["Vector", "Cut with EcoRI + BamHI"]);
        assert_eq!(p.rows[0].circular, Some(true), "the plasmid as it stands");
        assert_eq!(p.rows[1].circular, Some(false), "cutting opens it");
        assert!(p.rows[1].pieces[0].length.is_none(), "no guessed length");
        assert!(
            !p.rows[1].pieces[0].insert,
            "the vector is still the vector"
        );
        assert!(p.caption.contains("this step"), "{}", p.caption);

        // Step 3: the insert, and the piece that goes in. Nothing about the product yet.
        e.step = ExperimentStep::Insert;
        let p = plan(&e, &measured, &enzymes, &[]);
        let titles: Vec<&str> = p.rows.iter().map(|r| r.title.as_str()).collect();
        assert_eq!(titles, vec!["Insert", "Cut with EcoRI + BamHI"]);
        assert!(p.rows.iter().all(|r| !r.title.contains("Product")));
        assert!(p.rows[1].pieces[0].insert, "what goes in is red");

        // Step 4: the whole reaction.
        e.step = ExperimentStep::Product;
        let p = plan(&e, &measured, &enzymes, &[]);
        assert_eq!(p.rows.len(), 3);
        assert_eq!(p.rows[2].title, "Product \u{b7} 3950 bp");
    }

    /// Without enzymes the step still says what it will do, and asks for them.
    #[test]
    fn a_step_without_enzymes_asks_for_them() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        e.vector.enzymes = String::new();
        e.step = ExperimentStep::Vector;
        let p = plan(&e, &Measured::default(), &[], &[]);
        assert_eq!(p.rows[1].title, "Cut (name the enzymes)");
    }

    /// Gateway names att sites, not enzymes, and mutagenesis has no insert row at all.
    #[test]
    fn each_technique_tells_its_own_story() {
        let e = experiment(CloningTechnique::GatewayLr);
        let p = plan(
            &e,
            &Measured {
                vector_len: Some(5000),
                insert_len: Some(2000),
                vector_circular: Some(true),
                insert_circular: Some(true),
                ..Measured::default()
            },
            &[],
            &[],
        );
        assert_eq!(p.rows[1].left.as_deref(), Some("attL1"));
        assert_eq!(p.rows[2].title, "Product");

        let e = experiment(CloningTechnique::Mutagenesis);
        let p = plan(
            &e,
            &Measured {
                vector_len: Some(5000),
                vector_circular: Some(true),
                ..Measured::default()
            },
            &[],
            &[],
        );
        assert_eq!(p.rows.len(), 2, "no insert row: {p:?}");
        // Row titles stay short enough for their column; what the reaction does to the
        // vector is said once, in the caption.
        assert_eq!(p.rows[0].title, "Vector");
        assert!(p.caption.contains("primers"), "{}", p.caption);
    }

    /// The pieces a cut really makes are what the step draws: the backbone the reaction keeps
    /// as a circle with a gap in it, and the stuffer between the two sites in that gap. This
    /// is what the step must show: what the action does to the DNA.
    #[test]
    fn a_cut_plasmid_is_a_broken_circle_with_its_stuffer_in_the_gap() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        e.step = ExperimentStep::Vector;
        let measured = Measured {
            vector_len: Some(3200),
            vector_circular: Some(true),
            // EcoRI and BamHI, six bases apart on a 3200 bp plasmid.
            vector_fragments: vec![10, 3198],
            ..Measured::default()
        };
        let p = plan(&e, &measured, &["EcoRI".into(), "BamHI".into()], &[]);
        assert_eq!(p.rows[0].ring, Ring::Closed, "the plasmid as it stands");
        let cut = &p.rows[1];
        assert_eq!(cut.ring, Ring::Open, "cutting it leaves a gap");
        // The backbone is kept, the small piece is cut away, and both are to their real size.
        let kept: Vec<&Piece> = cut.pieces.iter().filter(|x| !x.removed).collect();
        let gone: Vec<&Piece> = cut.pieces.iter().filter(|x| x.removed).collect();
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].length, Some(3198));
        assert_eq!(gone.len(), 1);
        assert_eq!(gone[0].length, Some(10));
        assert_eq!(
            gone[0].label, "10 bp",
            "a piece cut away says how big it is"
        );

        // One cut: nothing is taken out, so the circle is simply broken open.
        let once = Measured {
            vector_fragments: vec![3204],
            ..measured.clone()
        };
        let p = plan(&e, &once, &["EcoRI".into()], &[]);
        assert_eq!(p.rows[1].ring, Ring::Open);
        assert_eq!(p.rows[1].pieces.len(), 1);
        assert!(
            p.rows[1].pieces.iter().all(|x| !x.removed),
            "a single cut throws nothing away"
        );
    }

    /// A linear insert cut at both ends keeps the middle and loses its flanks; it is a line,
    /// not a circle, because that is what it is.
    #[test]
    fn a_cut_insert_keeps_the_middle() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        e.step = ExperimentStep::Insert;
        let p = plan(
            &e,
            &Measured {
                insert_len: Some(762),
                insert_circular: Some(false),
                insert_fragments: vec![5, 760, 5],
                ..Measured::default()
            },
            &["EcoRI".into(), "BamHI".into()],
            &[],
        );
        let cut = &p.rows[1];
        assert_eq!(cut.ring, Ring::None, "a linear piece stays a line");
        let kept: Vec<&Piece> = cut.pieces.iter().filter(|x| !x.removed).collect();
        assert_eq!(kept.len(), 1);
        assert_eq!(kept[0].length, Some(760));
        assert!(kept[0].insert, "what goes in is the red piece");
        assert!(
            cut.pieces.iter().filter(|x| x.removed).all(|x| !x.insert),
            "the flanks are not drawn as the insert"
        );
    }

    /// Without a digest nothing is drawn to a size, and nothing is drawn as a circle: the
    /// schematic says only what has been read.
    #[test]
    fn without_a_digest_no_shape_is_claimed() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        e.step = ExperimentStep::Vector;
        let p = plan(
            &e,
            &Measured {
                vector_len: Some(3200),
                vector_circular: Some(true),
                ..Measured::default()
            },
            &["EcoRI".into()],
            &[],
        );
        assert_eq!(p.rows[1].ring, Ring::None);
        assert!(p.rows[1].pieces.iter().all(|x| x.length.is_none()));
        assert!(p.rows[1].pieces.iter().all(|x| !x.removed));
    }

    /// The plan says what the drawing can't: that a piece won't go in, or that it will go in
    /// either way round. Both belong to the steps about the piece, not to the vector's own.
    #[test]
    fn the_plan_warns_about_fit_and_direction() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        let measured = |fits| Measured {
            vector_len: Some(3200),
            insert_len: Some(762),
            vector_circular: Some(true),
            insert_circular: Some(false),
            fits: Some(fits),
            ..Measured::default()
        };
        e.step = ExperimentStep::Insert;
        let p = plan(&e, &measured(Fit::EitherWay), &["EcoRI".into()], &[]);
        assert_eq!(p.warnings.len(), 1);
        assert!(p.warnings[0].contains("either way round"), "{p:?}");

        let p = plan(&e, &measured(Fit::No), &["EcoRI".into()], &[]);
        assert!(p.warnings[0].contains("ends don't match"), "{p:?}");

        // A clone whose direction is settled says nothing.
        let p = plan(&e, &measured(Fit::OneWay), &["EcoRI".into()], &[]);
        assert!(p.warnings.is_empty(), "{p:?}");
        // Nor does one where it hasn't been worked out.
        let p = plan(&e, &Measured::default(), &["EcoRI".into()], &[]);
        assert!(p.warnings.is_empty(), "{p:?}");

        // The product step carries the warning too: it is the last chance to see it.
        e.step = ExperimentStep::Product;
        let p = plan(&e, &measured(Fit::EitherWay), &["EcoRI".into()], &[]);
        assert_eq!(p.warnings.len(), 1);
        // The vector's own step is about the backbone.
        e.step = ExperimentStep::Vector;
        let p = plan(&e, &measured(Fit::EitherWay), &["EcoRI".into()], &[]);
        assert!(p.warnings.is_empty(), "{p:?}");

        // A vector cut into more than two pieces: which one is the backbone is a choice the
        // reaction makes by size, so the step says so.
        let many = Measured {
            vector_fragments: vec![3000, 150, 50],
            ..measured(Fit::OneWay)
        };
        let p = plan(&e, &many, &["EcoRI".into(), "BamHI".into()], &[]);
        assert_eq!(p.warnings.len(), 1, "{p:?}");
        assert!(p.warnings[0].contains("into 3 pieces"), "{p:?}");
        // Two pieces is the ordinary case and needs no warning.
        let two = Measured {
            vector_fragments: vec![3190, 10],
            ..measured(Fit::OneWay)
        };
        assert!(
            plan(&e, &two, &["EcoRI".into(), "BamHI".into()], &[])
                .warnings
                .is_empty()
        );
    }

    /// The product is drawn as a circle only when its arcs add up to the product the reaction
    /// reports; where the pieces are an approximation, it stays a line.
    #[test]
    fn the_product_is_a_circle_only_when_its_pieces_add_up() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        e.step = ExperimentStep::Product;
        let p = plan(
            &e,
            &Measured {
                vector_len: Some(3200),
                insert_len: Some(762),
                product_len: Some(3950),
                vector_fragments: vec![10, 3198],
                insert_fragments: vec![5, 760, 5],
                vector_circular: Some(true),
                insert_circular: Some(false),
                ..Measured::default()
            },
            &["EcoRI".into(), "BamHI".into()],
            &[],
        );
        let product = p.rows.last().expect("the product row");
        assert_eq!(product.ring, Ring::Closed);
        // The arcs are the backbone and the piece that went in, not the whole constructs.
        let lengths: Vec<Option<u64>> = product.pieces.iter().map(|x| x.length).collect();
        assert_eq!(lengths, vec![Some(3198), Some(760)]);

        // Gateway: the backbone is not the whole vector, so vector + insert overshoots the
        // product badly and the row stays a line rather than lying about the proportions.
        let mut g = experiment(CloningTechnique::GatewayLr);
        g.step = ExperimentStep::Product;
        let p = plan(
            &g,
            &Measured {
                vector_len: Some(5000),
                insert_len: Some(2000),
                product_len: Some(5500),
                vector_circular: Some(true),
                insert_circular: Some(true),
                ..Measured::default()
            },
            &[],
            &[],
        );
        assert_eq!(p.rows.last().map(|r| r.ring), Some(Ring::None));
    }

    /// The schematic draws: circles, the piece in the gap, the line rows and the warning
    /// under them. `SNAPPIESNAP_PLAN_PNG=/path.png` to look at it.
    #[test]
    fn the_schematic_draws() {
        use crate::state::ExperimentStep;
        let mut e = experiment(CloningTechnique::Restriction);
        e.step = ExperimentStep::Product;
        let p = plan(
            &e,
            &Measured {
                vector_len: Some(3200),
                insert_len: Some(762),
                product_len: Some(3950),
                vector_fragments: vec![10, 3198],
                insert_fragments: vec![5, 760, 5],
                vector_circular: Some(true),
                insert_circular: Some(false),
                fits: Some(Fit::EitherWay),
                ..Measured::default()
            },
            &["EcoRI".into(), "BamHI".into()],
            &[],
        );
        for t in [Tokens::light(), Tokens::dark()] {
            let image = crate::offscreen::render(
                660,
                420,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    show(ui, &t, &p, 620.0);
                    legend(ui, &t, &p);
                },
            )
            .unwrap();
            if let Ok(path) = std::env::var("SNAPPIESNAP_PLAN_PNG") {
                let name = if t.dark { "_dark.png" } else { "_light.png" };
                crate::offscreen::write_png(&path.replace(".png", name), &image).unwrap();
            }
        }
    }

    /// A construct that hasn't been read has no length, and is not drawn to a guessed one.
    #[test]
    fn unread_constructs_are_not_given_a_size() {
        let e = experiment(CloningTechnique::Restriction);
        let p = plan(&e, &Measured::default(), &["EcoRI".into()], &[]);
        assert!(
            p.rows
                .iter()
                .all(|r| r.pieces.iter().all(|x| x.length.is_none()))
        );
        // Nothing is claimed about a topology that hasn't been read either.
        assert_eq!(p.rows[0].circular, None);
        assert_eq!(p.rows[1].circular, None);
    }
}
