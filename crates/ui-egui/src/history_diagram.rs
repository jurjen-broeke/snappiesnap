//! The construction history as a diagram: a card with a small map for every step, the product
//! at the top and the molecules it was made from below it, joined by lines labelled with the
//! operation. Branches fold and unfold; steps with a stored molecule open in a tab of their own.

use std::collections::BTreeSet;

use eframe::egui::{
    self, Align2, FontId, Pos2, Rect, RichText, Sense, Stroke, Ui, UiBuilder, vec2,
};
use snappiesnap_engine::doc::{HistoryNode, MAX_HISTORY_DEPTH, Molecule, operation_label};

use crate::theme::Tokens;

const CARD_W: f32 = 210.0;
const THUMB_H: f32 = 150.0;
const CARD_H: f32 = THUMB_H + 74.0;
const GAP_X: f32 = 22.0;
/// Room between a card and the row below it, for the operation label and the lines.
const GAP_Y: f32 = 62.0;
/// Steps below this depth start folded.
const OPEN_DEPTH: usize = 2;

/// What the user asked for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Action {
    Open(u32),
}

struct Placed<'a> {
    node: &'a HistoryNode,
    rect: Rect,
    depth: usize,
    parent: Option<usize>,
    folded: bool,
}

/// Folding is remembered per document while the app runs: the ids whose state differs from
/// the default (open above [`OPEN_DEPTH`], folded below).
fn toggled_id(ui: &Ui, doc: Option<u64>) -> egui::Id {
    ui.scope_id().with(("history-diagram-toggled", doc))
}

fn is_folded(toggled: &BTreeSet<u32>, n: &HistoryNode, depth: usize) -> bool {
    !n.inputs.is_empty()
        && (depth >= MAX_HISTORY_DEPTH || (depth >= OPEN_DEPTH) != toggled.contains(&n.id))
}

/// Width of the subtree under `n` as drawn.
fn subtree_width(n: &HistoryNode, depth: usize, toggled: &BTreeSet<u32>) -> f32 {
    if is_folded(toggled, n, depth) || n.inputs.is_empty() {
        return CARD_W;
    }
    let children: f32 = n
        .inputs
        .iter()
        .map(|i| subtree_width(i, depth + 1, toggled))
        .sum::<f32>()
        + GAP_X * n.inputs.len().saturating_sub(1) as f32;
    children.max(CARD_W)
}

fn place<'a>(
    n: &'a HistoryNode,
    left: f32,
    depth: usize,
    parent: Option<usize>,
    toggled: &BTreeSet<u32>,
    out: &mut Vec<Placed<'a>>,
) {
    let width = subtree_width(n, depth, toggled);
    let top = depth as f32 * (CARD_H + GAP_Y);
    let x = left + (width - CARD_W) / 2.0;
    let folded = is_folded(toggled, n, depth);
    let index = out.len();
    out.push(Placed {
        node: n,
        rect: Rect::from_min_size(Pos2::new(x, top), vec2(CARD_W, CARD_H)),
        depth,
        parent,
        folded,
    });
    if folded || n.inputs.is_empty() {
        return;
    }
    let children: f32 = n
        .inputs
        .iter()
        .map(|i| subtree_width(i, depth + 1, toggled))
        .sum::<f32>()
        + GAP_X * n.inputs.len().saturating_sub(1) as f32;
    let mut x = left + (width - children) / 2.0;
    for i in &n.inputs {
        let w = subtree_width(i, depth + 1, toggled);
        place(i, x, depth + 1, Some(index), toggled, out);
        x += w + GAP_X;
    }
}

/// Draw the diagram of `root` (the history of `current`, the open document).
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    doc_id: Option<u64>,
    current: &Molecule,
    root: &HistoryNode,
) -> Option<Action> {
    let key = toggled_id(ui, doc_id);
    let mut toggled: BTreeSet<u32> = ui.data(|d| d.get_temp(key)).unwrap_or_default();
    let mut placed = Vec::new();
    place(root, 0.0, 0, None, &toggled, &mut placed);
    let width = placed.iter().map(|p| p.rect.max.x).fold(CARD_W, f32::max);
    let height = placed.iter().map(|p| p.rect.max.y).fold(CARD_H, f32::max);
    let (area, _) = ui.allocate_exact_size(vec2(width + 2.0, height + 2.0), Sense::hover());
    let origin = area.min.to_vec2() + vec2(1.0, 1.0);
    let painter = ui.painter_at(area.expand(1.0));
    let line = Stroke::new(1.5, t.border);

    // Lines first, then the operation under each parent, then the cards.
    for (i, p) in placed.iter().enumerate() {
        let children: Vec<&Placed> = placed.iter().filter(|c| c.parent == Some(i)).collect();
        if children.is_empty() {
            continue;
        }
        let from = Pos2::new(p.rect.center().x, p.rect.max.y) + origin;
        let bus = from.y + GAP_Y * 0.62;
        painter.line_segment([from, Pos2::new(from.x, bus)], line);
        let xs = children.iter().map(|c| c.rect.center().x + origin.x);
        let (lo, hi) = xs.fold((from.x, from.x), |(lo, hi), x| (lo.min(x), hi.max(x)));
        painter.line_segment([Pos2::new(lo, bus), Pos2::new(hi, bus)], line);
        for c in &children {
            let x = c.rect.center().x + origin.x;
            painter.line_segment(
                [Pos2::new(x, bus), Pos2::new(x, c.rect.min.y + origin.y)],
                line,
            );
        }
        let label = operation_label(&p.node.operation);
        if !label.is_empty() {
            let galley = painter.layout_no_wrap(label, FontId::proportional(12.0), t.accent);
            let pos = Pos2::new(from.x + 6.0, from.y + 4.0);
            painter.galley(pos, galley, t.accent);
        }
    }

    let mut action = None;
    for (i, p) in placed.iter().enumerate() {
        let rect = p.rect.translate(origin);
        if !ui.is_rect_visible(rect) {
            continue;
        }
        let molecule = if i == 0 {
            Some(current)
        } else {
            p.node.molecule.as_deref()
        };
        let mut ui = ui.new_child(
            UiBuilder::new()
                .max_rect(rect)
                .id_salt(("history-card", p.node.id)),
        );
        if let Some(a) = card(&mut ui, t, p, i == 0, molecule, &mut toggled) {
            action = Some(a);
        }
    }
    ui.data_mut(|d| d.insert_temp(key, toggled));
    action
}

fn card(
    ui: &mut Ui,
    t: &Tokens,
    p: &Placed,
    root: bool,
    molecule: Option<&Molecule>,
    toggled: &mut BTreeSet<u32>,
) -> Option<Action> {
    let mut action = None;
    let n = p.node;
    let rect = ui.max_rect();
    let stroke = if root {
        Stroke::new(1.5, t.accent)
    } else {
        Stroke::new(1.0, t.border)
    };
    ui.painter()
        .rect(rect, t.radius, t.canvas, stroke, egui::StrokeKind::Inside);
    let inner = rect.shrink(6.0);
    let mut ui = ui.new_child(UiBuilder::new().max_rect(inner));
    ui.spacing_mut().item_spacing.y = 2.0;
    let thumb = Rect::from_min_size(inner.min, vec2(inner.width(), THUMB_H));
    match molecule {
        Some(m) => {
            let mut tui = ui.new_child(UiBuilder::new().max_rect(thumb).sense(Sense::hover()));
            crate::map::thumbnail_sized(&mut tui, t, m, thumb.size(), false);
        }
        None => {
            ui.painter().text(
                thumb.center(),
                Align2::CENTER_CENTER,
                "sequence not stored",
                FontId::proportional(12.0),
                t.text_muted,
            );
        }
    }
    ui.advance_cursor_after_rect(thumb);
    let name = n.name.strip_suffix(".dna").unwrap_or(&n.name);
    ui.add(egui::Label::new(RichText::new(name).color(t.text).strong().size(12.0)).truncate())
        .on_hover_text(name);
    // An open branch has its operation on the line below the card.
    let what = if n.operation.is_empty() {
        " · Starting molecule".to_string()
    } else if p.folded || n.inputs.is_empty() {
        format!(" · {}", operation_label(&n.operation))
    } else {
        String::new()
    };
    let topology = if n.topology == snappiesnap_engine::doc::Topology::Circular {
        "circular"
    } else {
        "linear"
    };
    ui.label(
        RichText::new(format!("{} bp, {topology}{what}", n.length))
            .color(t.text_muted)
            .size(11.0),
    );
    ui.horizontal(|ui| {
        if root {
            ui.label(
                RichText::new("this document")
                    .color(t.text_muted)
                    .size(11.0),
            );
        } else if n.molecule.is_some()
            && ui
                .small_button("Open")
                .on_hover_text("Open this step in a new tab")
                .clicked()
        {
            action = Some(Action::Open(n.id));
        }
        if !n.inputs.is_empty() {
            let text = if p.folded {
                format!("Show {} inputs", n.inputs.len())
            } else {
                "Hide inputs".to_string()
            };
            if p.depth < MAX_HISTORY_DEPTH
                && ui.small_button(text).clicked()
                && !toggled.remove(&n.id)
            {
                toggled.insert(n.id);
            }
        }
    });
    if !n.details.is_empty() {
        ui.label(RichText::new(&n.details).color(t.text_muted).size(11.0))
            .on_hover_text(&n.details);
    }
    action
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_engine::doc::{Alphabet, Seq, Topology};

    fn mol(name: &str, n: usize) -> Molecule {
        Molecule::new(
            name,
            Seq::new(Alphabet::Dna, "ACGT".repeat(n)).unwrap(),
            Topology::Circular,
        )
        .unwrap()
    }

    #[test]
    fn layout_centres_parents_over_their_inputs() {
        let (a, b) = (mol("a", 10), mol("b", 20));
        let mut p = mol("p", 30);
        let h = HistoryNode::made(
            &p,
            "gibsonAssembly",
            "",
            vec![HistoryNode::input(&a), HistoryNode::input(&b)],
        );
        p.set_history(Some(h)).unwrap();
        let root = p.history().unwrap();
        let mut placed = Vec::new();
        place(root, 0.0, 0, None, &BTreeSet::new(), &mut placed);
        assert_eq!(placed.len(), 3);
        let [r, c1, c2] = placed.as_slice() else {
            panic!()
        };
        assert!((r.rect.center().x - (c1.rect.center().x + c2.rect.center().x) / 2.0).abs() < 0.01);
        assert!(c1.rect.max.x + GAP_X <= c2.rect.min.x + 0.01);
        assert_eq!(c1.parent, Some(0));
        // Folding the root hides its inputs.
        let mut placed = Vec::new();
        place(root, 0.0, 0, None, &BTreeSet::from([0]), &mut placed);
        assert_eq!(placed.len(), 1);
    }

    #[test]
    fn diagram_draws() {
        let (a, b) = (mol("a", 10), mol("b", 20));
        let mut ab = mol("ab", 30);
        let h = HistoryNode::made(
            &ab,
            "gibsonAssembly",
            "a; b",
            vec![HistoryNode::input(&a), HistoryNode::input(&b)],
        );
        ab.set_history(Some(h)).unwrap();
        let mut p = mol("p", 25);
        let h = HistoryNode::made(&p, "amplifyFragment", "", vec![HistoryNode::input(&ab)]);
        p.set_history(Some(h)).unwrap();
        let t = Tokens::light();
        let image = crate::offscreen::render(
            900,
            900,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let _ = show(ui, &t, None, &p, p.history().unwrap());
            },
        )
        .unwrap();
        assert_eq!(image.size, [900, 900]);
    }
}
