//! The sequence view: numbered lines of residues in groups of ten, the complementary strand,
//! and feature bars underneath, with the translation of CDS features written in their bars.
//! Only visible lines are laid out, so long sequences stay fast.
//!
//! It only reads the document. Clicks, drags and typing become [`SeqAction`]s that the shell
//! runs as engine commands, so every edit is undoable and identical to an agent's.

use eframe::egui::{
    self, Align2, Event, FontId, Key, Pos2, Rect, Sense, Shape, Stroke, Ui, pos2, vec2,
};
use snappiesnap_engine::doc::{Molecule, Range};
use snappiesnap_layout::{self as layout, Arrow};

use crate::state::TailStyle;
use crate::theme::Tokens;

/// Feature lanes shown under each line; more are summarised.
const MAX_LANES: usize = 4;
/// Rows of restriction site names above each line, and their height.
const SITE_ROWS: usize = 2;
const SITE_ROW_H: f32 = 12.0;

#[derive(Debug, Clone, PartialEq)]
pub enum SeqAction {
    Select(Option<Range>),
    Insert {
        at: usize,
        residues: String,
    },
    Replace {
        range: Range,
        residues: String,
    },
    Delete(Range),
    /// Copy the selection to the clipboard: "plain", "fasta" or "reverse_complement".
    Copy(&'static str),
    /// Copy the selection, then delete it.
    Cut,
    /// Paste this text at the cursor or over the selection.
    Paste(String),
    /// BLAST the selection at NCBI against "nucleotide", "protein" or "translated" databases.
    Blast(&'static str),
    /// A restriction site's name was double-clicked: its index in the sites given. The
    /// cloning experiment adds that enzyme to the cut list; the document view ignores it.
    SiteChosen(usize),
}

/// View state that isn't document state.
#[derive(Debug, Default, Clone)]
pub struct SeqViewState {
    drag_anchor: Option<usize>,
    /// Scroll so this position is in view on the next frame (find).
    pub scroll_to: Option<usize>,
    /// Translations of CDS features (see [`cds_codons`]), for this document revision.
    codons: Option<(u64, Codons)>,
    /// The residues drawn last frame, and where the selection started then: a selection that
    /// moves somewhere not on screen (Go to, Find, a command an agent ran) brings itself into
    /// view, while one made with the pointer, which is on screen already, leaves it alone.
    drawn: Option<(usize, usize)>,
    last_selection: Option<usize>,
}

/// Where to scroll so a selection can be seen, if anywhere: only when it has just moved, and
/// only when it landed outside the residues drawn last frame. A selection made with the
/// pointer is on screen already, so the view stays where the reader put it; one from Go to,
/// Find, or a command an agent ran usually isn't.
fn bring_into_view(
    selection: Option<usize>,
    last: Option<usize>,
    drawn: Option<(usize, usize)>,
) -> Option<usize> {
    let p = selection?;
    if selection == last {
        return None;
    }
    drawn.is_none_or(|(a, b)| p < a || p >= b).then_some(p)
}

/// Per feature index: the amino acid of each codon, at the top-strand position of its middle
/// base.
pub type Codons = std::collections::HashMap<usize, Vec<(usize, u8)>>;

/// The translations of the molecule's CDS features: segments joined in order on the
/// feature's strand, from /codon_start, with /transl_table (default 1).
pub fn cds_codons(m: &Molecule) -> Codons {
    use snappiesnap_engine::bio::translate::{code, translate};
    use snappiesnap_engine::doc::Strand;
    let len = m.len();
    let seq = m.seq().as_bytes();
    let mut out = Codons::new();
    if len == 0 || !m.seq().alphabet().is_nucleotide() {
        return out;
    }
    for (index, f) in m.features().iter().enumerate() {
        if f.kind != "CDS" {
            continue;
        }
        let q = |n: &str| {
            f.qualifiers
                .iter()
                .find(|q| q.name == n)
                .and_then(|q| q.values.first())
                .and_then(|v| v.trim().parse::<u8>().ok())
        };
        let Some(code) = code(q("transl_table").unwrap_or(1)) else {
            continue;
        };
        let reverse = f.strand == Strand::Reverse;
        // Top-strand positions of the coding bases, in reading order.
        let mut positions: Vec<usize> = Vec::new();
        for seg in &f.segments {
            let n = seg.range.len_on(len);
            let along = (0..n).map(|k| (seg.range.start + k) % len);
            if reverse {
                positions.extend(along.rev());
            } else {
                positions.extend(along);
            }
            if positions.len() > 30_000_000 {
                break;
            }
        }
        let skip = usize::from(q("codon_start").unwrap_or(1).clamp(1, 3)) - 1;
        let positions = positions.get(skip..).unwrap_or_default();
        let bases: Vec<u8> = positions
            .iter()
            .map(|&p| {
                let b = seq.get(p).copied().unwrap_or(b'N');
                if reverse {
                    match b.to_ascii_uppercase() {
                        b'A' => b'T',
                        b'T' => b'A',
                        b'C' => b'G',
                        b'G' => b'C',
                        o => o,
                    }
                } else {
                    b
                }
            })
            .collect();
        let protein = translate(&bases, code, true);
        let mut marks: Vec<(usize, u8)> = protein
            .iter()
            .enumerate()
            .filter_map(|(k, &aa)| Some((*positions.get(3 * k + 1)?, aa)))
            .collect();
        // By position along the top strand, for drawing line by line.
        marks.sort_unstable_by_key(|&(p, _)| p);
        out.insert(index, marks);
    }
    out
}

/// Height of one translation row.
const AA_ROW_H: f32 = 14.0;
/// Height of one ORF lane, and the most lanes shown.
const ORF_ROW_H: f32 = 7.0;
const MAX_ORF_ROWS: usize = 4;
/// Height of the ruler along the top in single-line mode.
const RULER_H: f32 = 14.0;
/// How far a primer base that doesn't pair steps up out of its row.
const MISMATCH_JOG: f32 = 6.0;

/// The codons of reading `frame` (1..=3 on the top strand from the first base, -1..=-3 on the
/// bottom strand from the last) whose middle base lies in `start..end`: (middle base, amino
/// acid, NCBI code `table`). Codons that would run past either end are left out.
pub fn frame_codons(
    m: &Molecule,
    frame: i8,
    start: usize,
    end: usize,
    table: u8,
) -> Vec<(usize, u8)> {
    use snappiesnap_engine::bio::translate::code;
    let mut out = Vec::new();
    let (Some(code), Some(offset)) = (code(table), Some(usize::from(frame.unsigned_abs()))) else {
        return out;
    };
    let seq = m.seq().as_bytes();
    let len = seq.len();
    if !(1..=3).contains(&offset) || len < offset + 2 {
        return out;
    }
    let comp = |b: u8| match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    };
    let first_mid = if frame > 0 {
        // Codons start at offset-1, offset+2, ...: middles at offset, offset+3, ...
        offset
    } else {
        // Read from the last base down: the first codon covers len-offset-2..=len-offset.
        (len - offset - 1) % 3
    };
    let k0 = start.saturating_sub(first_mid).div_ceil(3);
    let mut mid = first_mid + 3 * k0;
    while mid < end.min(len) {
        let (Some(&a), Some(&b), Some(&c)) = (
            mid.checked_sub(1).and_then(|i| seq.get(i)),
            seq.get(mid),
            seq.get(mid + 1),
        ) else {
            mid += 3;
            continue;
        };
        let in_frame = if frame > 0 {
            mid + 1 < len
        } else {
            mid < len - offset
        };
        if in_frame && mid >= 1 {
            let codon = if frame > 0 {
                [a, b, c]
            } else {
                [comp(c), comp(b), comp(a)]
            };
            out.push((mid, code.amino_acid(codon)));
        }
        mid += 3;
    }
    out
}

struct Metrics {
    /// Height of the restriction site label rows above each line (0 when sites are off).
    sites_h: f32,
    /// Translation rows under the two strands (+1..+3, -1..-3), and their genetic code.
    frames: Vec<i8>,
    table: u8,
    /// Layers shown (as on the map).
    layers: crate::map::Layers,
    /// Height of the translation strip above each feature lane (0 when no CDS is shown).
    /// The amino acids used to sit inside the feature's arrow, where its colour and its name
    /// hid them.
    cds_h: f32,
    /// Feature lanes (0 when features are hidden), primer rows (forward, reverse) and ORF lanes.
    feature_lanes: usize,
    primer_rows: usize,
    orf_rows: usize,
    font: FontId,
    char_w: f32,
    line_h: f32,
    gap: f32,
    lane_h: f32,
    number_w: f32,
    per_line: usize,
    row_h: f32,
    /// The whole sequence on one line, scrolled sideways: no grouping gaps, and a ruler of
    /// positions along the top in place of the number at the start of each line.
    single: bool,
    /// Height of that ruler (0 when wrapped).
    ruler_h: f32,
}

impl Metrics {
    #[allow(clippy::too_many_arguments)]
    fn new(
        ui: &Ui,
        len: usize,
        width: f32,
        lanes: usize,
        sites: bool,
        frames: &[i8],
        layers: crate::map::Layers,
        primer_rows: usize,
        orf_rows: usize,
        translated: bool,
        single: bool,
    ) -> Self {
        let font = FontId::monospace(13.0);
        let char_w = ui.ctx().fonts_mut(|f| f.glyph_width(&font, 'A')).max(1.0);
        let line_h = 17.0;
        // On one line the residues run on without the groups of ten, so a position can be
        // worked out from the ruler by counting characters.
        let gap = if single { 0.0 } else { char_w * 0.8 };
        let digits = len.max(1).to_string().len() as f32;
        let number_w = digits * char_w + 16.0;
        let group_w = 10.0 * char_w + gap;
        let groups = ((width - number_w - 24.0) / group_w).floor().max(1.0) as usize;
        let lane_h = 13.0;
        let sites_h = if sites {
            SITE_ROWS as f32 * SITE_ROW_H + 2.0
        } else {
            0.0
        };
        let aa_h = frames.len() as f32 * AA_ROW_H;
        let feature_lanes = if layers.features {
            lanes.clamp(1, MAX_LANES)
        } else {
            0
        };
        let extra = primer_rows as f32 * line_h + orf_rows as f32 * ORF_ROW_H;
        let cds_h = if layers.features && translated {
            AA_ROW_H
        } else {
            0.0
        };
        let ruler_h = if single && layers.ruler { RULER_H } else { 0.0 };
        Metrics {
            sites_h,
            cds_h,
            frames: frames.to_vec(),
            table: 1,
            layers,
            feature_lanes,
            primer_rows,
            orf_rows,
            font,
            char_w,
            line_h,
            gap,
            lane_h,
            number_w,
            per_line: if single { len.max(1) } else { groups * 10 },
            row_h: ruler_h
                + sites_h
                + 2.0 * line_h
                + aa_h
                + feature_lanes as f32 * (lane_h + cds_h)
                + extra
                + 14.0,
            single,
            ruler_h,
        }
    }

    /// Width of the whole line, gutter included (single-line mode: the scrolled content).
    fn line_width(&self, len: usize) -> f32 {
        self.col_x(0.0, len) + self.char_w + 12.0
    }

    /// x of the left edge of column `col` (0-based within the line).
    fn col_x(&self, left: f32, col: usize) -> f32 {
        left + self.number_w + col as f32 * self.char_w + (col / 10) as f32 * self.gap
    }

    /// The column whose left edge is closest to `x` (0..=per_line).
    fn col_at(&self, left: f32, x: f32) -> usize {
        let rel = (x - left - self.number_w).max(0.0);
        let group_w = 10.0 * self.char_w + self.gap;
        let group = (rel / group_w).floor() as usize;
        let within = ((rel - group as f32 * group_w) / self.char_w).round() as usize;
        (group * 10 + within.min(10)).min(self.per_line)
    }
}

#[allow(clippy::too_many_arguments)]
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    m: &Molecule,
    selection: Option<Range>,
    state: &mut SeqViewState,
    revision: u64,
    marks: &crate::enzymes::Marks,
    layers: crate::map::Layers,
    frames: &[i8],
    table: u8,
    found: &[Range],
    tails: TailStyle,
    // The whole sequence on one line, scrolled sideways, instead of wrapped into rows.
    single_line: bool,
) -> Vec<SeqAction> {
    let sites: &[crate::enzymes::MapSite] = if layers.enzymes { marks.sites } else { &[] };
    let primers: &[crate::enzymes::PrimerMark] = if layers.primers { marks.primers } else { &[] };
    let orfs: &[crate::enzymes::OrfMark] = if layers.orfs { marks.orfs } else { &[] };
    let (primers, primer_rows) = duplexes(
        primers,
        tails,
        &ui.ctx()
            .fonts_mut(|f| f.glyph_width(&FontId::monospace(13.0), 'A'))
            .max(1.0),
    );
    let primers = primers.as_slice();
    let orf_rows = orfs
        .iter()
        .map(|o| o.lane + 1)
        .max()
        .unwrap_or(0)
        .min(MAX_ORF_ROWS);
    let mut actions = Vec::new();
    let sel_start = selection.map(|r| r.start);
    if let Some(p) = bring_into_view(sel_start, state.last_selection, state.drawn) {
        state.scroll_to = Some(p);
    }
    state.last_selection = sel_start;
    if state.codons.as_ref().is_none_or(|(r, _)| *r != revision) {
        state.codons = Some((revision, cds_codons(m)));
    }
    let codons = state.codons.take().map(|(_, c)| c).unwrap_or_default();
    let len = m.len();
    let metrics = Metrics::new(
        ui,
        len,
        ui.available_width(),
        layout::max_lanes(m),
        !sites.is_empty(),
        if m.seq().alphabet().is_nucleotide() {
            frames
        } else {
            &[]
        },
        layers,
        primer_rows,
        orf_rows,
        !codons.is_empty(),
        single_line,
    );
    let metrics = Metrics { table, ..metrics };
    let lines = layout::line_count(len, metrics.per_line);
    let id = ui.scope_id().with("sequence-view");
    let mut rows: Vec<(Rect, usize)> = Vec::new();
    let mut pointer_actions: Vec<(Pos2, bool, bool, bool)> = Vec::new();
    // Where each enzyme's name was drawn on the line being handled, and which was chosen.
    let mut site_labels: Vec<(Rect, usize)> = Vec::new();
    let mut chosen_site: Option<usize> = None;
    let mut menu_actions: Vec<SeqAction> = Vec::new();

    // What each row of the scroll area draws: a wrapped line, or a window of the one long
    // line. The body is the same either way, so the two modes cannot drift apart.
    let draw_row = |ui: &mut Ui,
                    rect: Rect,
                    response: &egui::Response,
                    range: (usize, usize),
                    left: f32,
                    gutter_x: f32,
                    site_labels: &mut Vec<(Rect, usize)>,
                    menu_actions: &mut Vec<SeqAction>,
                    pointer_actions: &mut Vec<(Pos2, bool, bool, bool)>,
                    chosen_site: &mut Option<usize>| {
        let has_selection = selection.is_some_and(|r| !r.is_empty());
        response.context_menu(|ui| {
            let nucleotide = m.seq().alphabet().is_nucleotide();
            let items: [(&str, Option<SeqAction>, bool); 4] = [
                ("Copy", Some(SeqAction::Copy("plain")), true),
                ("Copy as FASTA", Some(SeqAction::Copy("fasta")), true),
                (
                    "Copy reverse complement",
                    Some(SeqAction::Copy("reverse_complement")),
                    nucleotide,
                ),
                ("Cut", Some(SeqAction::Cut), true),
            ];
            for (label, act, ok) in items {
                if ui
                    .add_enabled(has_selection && ok, egui::Button::new(label))
                    .on_disabled_hover_text("Select some residues first")
                    .clicked()
                {
                    if let Some(a) = act {
                        menu_actions.push(a);
                    }
                    ui.close();
                }
            }
            ui.separator();
            ui.add_enabled_ui(has_selection, |ui| {
                ui.menu_button("BLAST selection at NCBI", |ui| {
                    let choices: &[(&str, &'static str)] = if nucleotide {
                        &[
                            ("Nucleotide database (blastn)", "nucleotide"),
                            ("Protein database (blastx)", "protein"),
                            ("Translated nucleotide database (tblastx)", "translated"),
                        ]
                    } else {
                        &[
                            ("Protein database (blastp)", "protein"),
                            ("Translated nucleotide database (tblastn)", "translated"),
                        ]
                    };
                    for (label, against) in choices {
                        if ui
                            .button(*label)
                            .on_hover_text(
                                "Opens NCBI BLAST in your browser with the selection filled in",
                            )
                            .clicked()
                        {
                            menu_actions.push(SeqAction::Blast(against));
                            ui.close();
                        }
                    }
                });
            })
            .response
            .on_disabled_hover_text("Select some residues first");
            ui.separator();
            ui.label(egui::RichText::new("Paste: Cmd+V (Ctrl+V)").color(t.text_muted));
        });
        draw_line(
            ui,
            t,
            m,
            &metrics,
            rect,
            range,
            left,
            gutter_x,
            selection,
            &codons,
            sites,
            primers,
            orfs,
            found,
            response.hover_pos(),
            site_labels,
        );
        // Double-clicking a site's name chooses that enzyme, as on the map.
        if response.double_clicked()
            && let Some(p) = response.interact_pointer_pos()
            && let Some((_, i)) = site_labels.iter().find(|(r, _)| r.contains(p))
        {
            *chosen_site = Some(*i);
        }
        if let Some(p) = response.interact_pointer_pos() {
            pointer_actions.push((
                p,
                response.clicked(),
                response.drag_started(),
                response.dragged(),
            ));
        }
        if response.clicked() || response.drag_started() {
            ui.memory_mut(|mem| mem.request_focus(id));
        }
    };

    if single_line {
        // One line: allocate its whole width and draw only the window in view, so a chromosome
        // costs no more to show than a plasmid.
        let width = metrics.line_width(len).max(ui.available_width());
        let mut area = egui::ScrollArea::horizontal().auto_shrink([false, false]);
        if let Some(pos) = state.scroll_to.take() {
            // A little before it, so it isn't against the left edge.
            let x = metrics.col_x(0.0, pos.min(len)) - metrics.number_w - 80.0;
            area = area.horizontal_scroll_offset(x.max(0.0));
        }
        area.show_viewport(ui, |ui, viewport| {
            let (rect, response) = ui.allocate_exact_size(
                vec2(width, metrics.row_h.max(ui.available_height())),
                Sense::click_and_drag(),
            );
            rows.push((rect, 0));
            if len == 0 {
                ui.label("Empty sequence: type residues to add them.");
                return;
            }
            // The residues the viewport covers, with a column either side.
            let view_left = rect.left() + viewport.min.x;
            let first = ((viewport.min.x - metrics.number_w) / metrics.char_w).floor();
            let first = (first.max(0.0) as usize).min(len);
            let count = (viewport.width() / metrics.char_w).ceil() as usize + 2;
            let last = first.saturating_add(count).min(len);
            state.drawn = Some((first, last));
            // `left` is shifted so that column `first` lands where it belongs on the line.
            let left = rect.left() + first as f32 * metrics.char_w;
            draw_row(
                ui,
                Rect::from_min_max(
                    pos2(view_left, rect.top()),
                    pos2(view_left + viewport.width(), rect.bottom()),
                ),
                &response,
                (first, last),
                left,
                view_left + metrics.number_w - 10.0,
                &mut site_labels,
                &mut menu_actions,
                &mut pointer_actions,
                &mut chosen_site,
            );
        });
    } else {
        let mut area = egui::ScrollArea::vertical().auto_shrink([false, false]);
        if let Some(pos) = state.scroll_to.take() {
            let line = pos / metrics.per_line.max(1);
            // One line of context above.
            area = area.vertical_scroll_offset((line.saturating_sub(1)) as f32 * metrics.row_h);
        }
        area.show_rows(ui, metrics.row_h, lines.max(1), |ui, visible| {
            state.drawn = Some((
                layout::line_range(len, metrics.per_line, visible.start).0,
                layout::line_range(len, metrics.per_line, visible.end.saturating_sub(1)).1,
            ));
            for line in visible {
                let (rect, response) = ui.allocate_exact_size(
                    vec2(ui.available_width(), metrics.row_h),
                    Sense::click_and_drag(),
                );
                rows.push((rect, line));
                draw_row(
                    ui,
                    rect,
                    &response,
                    layout::line_range(len, metrics.per_line, line),
                    rect.left(),
                    rect.left() + metrics.number_w - 10.0,
                    &mut site_labels,
                    &mut menu_actions,
                    &mut pointer_actions,
                    &mut chosen_site,
                );
            }
            if len == 0 {
                ui.label("Empty sequence: type residues to add them.");
            }
        });
    }

    // Pointer → position, using the rows laid out this frame.
    let pos_at = |p: Pos2| -> Option<usize> {
        let (rect, line) = rows
            .iter()
            .find(|(r, _)| p.y >= r.top() && p.y < r.bottom())
            .or_else(|| {
                if rows.first().is_some_and(|(r, _)| p.y < r.top()) {
                    rows.first()
                } else {
                    rows.last()
                }
            })?;
        let (start, end) = layout::line_range(len, metrics.per_line, *line);
        let col = metrics.col_at(rect.left(), p.x);
        Some(start.saturating_add(col).min(end))
    };
    if let Some(i) = chosen_site {
        actions.push(SeqAction::SiteChosen(i));
    }
    for (p, clicked, drag_started, dragged) in pointer_actions {
        let Some(pos) = pos_at(p) else { continue };
        if drag_started {
            state.drag_anchor = Some(pos);
        }
        if dragged && let Some(anchor) = state.drag_anchor {
            let (a, b) = (anchor.min(pos), anchor.max(pos));
            actions.push(SeqAction::Select(Some(Range::new(a, b))));
        } else if clicked {
            state.drag_anchor = None;
            actions.push(SeqAction::Select(Some(Range::at(pos))));
        }
    }

    actions.extend(menu_actions);
    if ui.memory(|mem| mem.has_focus(id)) {
        keyboard(ui, m, selection, &mut actions);
    }
    state.codons = Some((revision, codons));
    actions
}

/// Restriction sites cutting this line: their names in rows above it, and each cut as a stepped
/// line through both strands (down the top strand at the top-strand cut, across, down the
/// bottom strand at the bottom-strand cut), as in the map's hover card.
#[allow(clippy::too_many_arguments)]
fn draw_sites(
    painter: &egui::Painter,
    t: &Tokens,
    mx: &Metrics,
    sites: &[crate::enzymes::MapSite],
    left: f32,
    gutter_x: f32,
    top: f32,
    start: usize,
    end: usize,
    len: usize,
    positions: bool,
    // Where the pointer is: the enzyme whose name it is over has its cut drawn thick, in
    // `t.highlight`, as the map's hover picks out a site.
    hover: Option<Pos2>,
    // Where each name was drawn, so a double-click can tell which enzyme it hit.
    labels: &mut Vec<(Rect, usize)>,
) {
    if sites.is_empty() {
        return;
    }
    // x of the boundary before position `p` of this line (cuts sit between bases).
    let x_at = |p: usize| -> Option<f32> {
        if p < start || p > end {
            return None;
        }
        let col = p - start;
        Some(if col == 0 {
            mx.col_x(left, 0) - 1.0
        } else {
            mx.col_x(left, col - 1)
                + mx.char_w
                + if col.is_multiple_of(10) {
                    mx.gap / 2.0
                } else {
                    0.0
                }
        })
    };
    let font = FontId::proportional(10.0);
    let mut rows: Vec<f32> = vec![f32::NEG_INFINITY; SITE_ROWS];
    let mut hidden = 0;
    let mid = top + mx.line_h;
    for (index, s) in sites.iter().enumerate() {
        // A cut at the origin of a circular molecule is after the last base.
        let cut = if s.mark.pos == 0 && len > 0 {
            len
        } else {
            s.mark.pos
        };
        let Some(xt) = x_at(cut).or_else(|| x_at(s.mark.pos)) else {
            continue;
        };
        let xb = x_at(if s.cut_bottom == 0 { len } else { s.cut_bottom }).unwrap_or(xt);
        // The name in the first row with room; placed before the cut is drawn, because
        // hovering the name changes how the cut looks.
        let label = if positions {
            s.mark.text.clone()
        } else {
            s.detail.label.clone()
        };
        // Laid out without a colour, so the hovered name can be drawn in another one.
        let galley = painter.layout_no_wrap(label, font.clone(), egui::Color32::PLACEHOLDER);
        let w = galley.size().x;
        let h = galley.size().y;
        let lx = xt - w / 2.0;
        let placed = rows.iter().position(|&e| e + 4.0 < lx).map(|r| {
            let y = top - mx.sites_h + 1.0 + r as f32 * SITE_ROW_H;
            (r, y, Rect::from_min_size(pos2(lx, y), vec2(w, h)))
        });
        let hovered = placed.is_some_and(|(_, _, rect)| hover.is_some_and(|p| rect.contains(p)));
        let stroke = if hovered {
            Stroke::new(2.5, t.highlight)
        } else {
            Stroke::new(1.0, t.accent.gamma_multiply(0.8))
        };
        painter.add(Shape::line(
            vec![
                pos2(xt, top - 1.0),
                pos2(xt, mid),
                pos2(xb, mid),
                pos2(xb, top + 2.0 * mx.line_h),
            ],
            stroke,
        ));
        match placed {
            Some((r, y, rect)) => {
                labels.push((rect, index));
                let (text, leader) = if hovered {
                    (t.highlight, Stroke::new(1.5, t.highlight))
                } else {
                    (
                        t.text_muted,
                        Stroke::new(0.6, t.text_muted.gamma_multiply(0.5)),
                    )
                };
                painter.galley(pos2(lx, y), galley, text);
                painter.line_segment(
                    [pos2(xt, y + SITE_ROW_H - 1.0), pos2(xt, top - 1.0)],
                    leader,
                );
                if let Some(e) = rows.get_mut(r) {
                    *e = lx + w;
                }
            }
            None => hidden += 1,
        }
    }
    if hidden > 0 {
        painter.text(
            pos2(gutter_x, top - mx.sites_h + 1.0),
            Align2::RIGHT_TOP,
            format!("+{hidden}"),
            font,
            t.text_muted,
        );
    }
}

/// Typing edits at the cursor or replaces the selection; Backspace/Delete remove residues;
/// arrow keys move the cursor; copy, cut and paste use the clipboard.
fn keyboard(ui: &Ui, m: &Molecule, selection: Option<Range>, actions: &mut Vec<SeqAction>) {
    let len = m.len();
    let sel = selection.unwrap_or(Range::at(len));
    let alphabet = m.seq().alphabet();
    let events = ui.input(|i| i.events.clone());
    for event in events {
        match event {
            Event::Text(text) => {
                let residues: String = text
                    .chars()
                    .filter(|c| c.is_ascii() && alphabet.is_valid(*c as u8))
                    .map(|c| c.to_ascii_uppercase())
                    .collect();
                if residues.is_empty() {
                    continue;
                }
                if sel.is_empty() {
                    actions.push(SeqAction::Insert {
                        at: sel.start.min(len),
                        residues,
                    });
                } else {
                    actions.push(SeqAction::Replace {
                        range: sel,
                        residues,
                    });
                }
            }
            Event::Copy if !sel.is_empty() => actions.push(SeqAction::Copy("plain")),
            Event::Cut if !sel.is_empty() => actions.push(SeqAction::Cut),
            Event::Paste(text) => actions.push(SeqAction::Paste(text)),
            Event::Key {
                key, pressed: true, ..
            } => match key {
                Key::Backspace if !sel.is_empty() => actions.push(SeqAction::Delete(sel)),
                Key::Backspace if sel.start > 0 => {
                    actions.push(SeqAction::Delete(Range::new(sel.start - 1, sel.start)))
                }
                Key::Delete if !sel.is_empty() => actions.push(SeqAction::Delete(sel)),
                Key::Delete if sel.start < len => {
                    actions.push(SeqAction::Delete(Range::new(sel.start, sel.start + 1)))
                }
                Key::ArrowLeft => actions.push(SeqAction::Select(Some(Range::at(
                    sel.start.saturating_sub(1),
                )))),
                Key::ArrowRight => {
                    let to = if sel.is_empty() {
                        sel.start + 1
                    } else {
                        sel.end
                    };
                    actions.push(SeqAction::Select(Some(Range::at(to.min(len)))));
                }
                _ => {}
            },
            _ => {}
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_line(
    ui: &Ui,
    t: &Tokens,
    m: &Molecule,
    mx: &Metrics,
    rect: Rect,
    // The residues to draw, where column `start` goes (`left` is what `col_x` measures from,
    // so a window of a long line can be drawn where it belongs), and the x of the right edge
    // of the row labels, pinned to the viewport in single-line mode so they stay in view.
    (start, end): (usize, usize),
    left: f32,
    gutter_x: f32,
    selection: Option<Range>,
    codons: &Codons,
    sites: &[crate::enzymes::MapSite],
    primers: &[Duplex<'_>],
    orfs: &[crate::enzymes::OrfMark],
    found: &[Range],
    // The pointer, when it is over this line: a hovered enzyme name highlights its cut.
    hover: Option<Pos2>,
    labels: &mut Vec<(Rect, usize)>,
) {
    let painter = ui.painter_at(rect);
    let len = m.len();
    let top = rect.top() + 4.0 + mx.ruler_h + mx.sites_h;
    draw_sites(
        &painter,
        t,
        mx,
        sites,
        left,
        gutter_x,
        top,
        start,
        end,
        len,
        mx.layers.positions,
        hover,
        labels,
    );
    let bottom_y = top + mx.line_h;

    // Find matches, tinted under both strands.
    for r in found {
        for (s, e, _, _) in pieces(*r, len, start, end) {
            let rect = Rect::from_min_max(
                pos2(mx.col_x(left, s - start) - 1.0, top),
                pos2(
                    mx.col_x(left, e - 1 - start) + mx.char_w + 1.0,
                    top + 2.0 * mx.line_h,
                ),
            );
            painter.rect_filled(rect, 2.0, t.accent.gamma_multiply(0.18));
        }
    }

    // Selection and cursor.
    if let Some(sel) = selection {
        if sel.is_empty() {
            if sel.start >= start && sel.start <= end && !(sel.start == end && end < len) {
                let x = mx.col_x(left, sel.start - start) - 1.0;
                painter.line_segment(
                    [pos2(x, top), pos2(x, top + 2.0 * mx.line_h)],
                    Stroke::new(2.0, t.accent),
                );
            }
        } else {
            let (a, b) = sel.segments(len);
            for piece in std::iter::once(a).chain(b) {
                let (s, e) = (piece.start.max(start), piece.end.min(end));
                if s < e {
                    let r = Rect::from_min_max(
                        pos2(mx.col_x(left, s - start) - 1.0, top),
                        pos2(
                            mx.col_x(left, e - 1 - start) + mx.char_w + 1.0,
                            top + 2.0 * mx.line_h,
                        ),
                    );
                    painter.rect_filled(r, 2.0, t.selection);
                }
            }
        }
    }

    // Position number (the ruler), residues in groups of ten, complementary strand.
    if mx.layers.ruler {
        if mx.single {
            // One line: a ruler of positions along the top, every ten residues.
            let first = start.div_ceil(10) * 10;
            let mut p = first.max(10);
            while p <= end {
                let x = mx.col_x(left, p - start) - mx.char_w / 2.0;
                let y = top - mx.sites_h - 3.0;
                painter.line_segment(
                    [pos2(x, y), pos2(x, y - 4.0)],
                    Stroke::new(1.0, t.text_muted),
                );
                painter.text(
                    pos2(x + 2.0, y - RULER_H + 1.0),
                    Align2::LEFT_TOP,
                    p.to_string(),
                    FontId::monospace(10.0),
                    t.text_muted,
                );
                p += 10;
            }
        } else {
            painter.text(
                pos2(gutter_x, top + mx.line_h / 2.0),
                Align2::RIGHT_CENTER,
                (start + 1).to_string(),
                FontId::monospace(11.0),
                t.text_muted,
            );
        }
    }
    let bytes = m.seq().as_bytes().get(start..end).unwrap_or_default();
    let alphabet = m.seq().alphabet();
    for (group, chunk) in bytes.chunks(10).enumerate() {
        let x = mx.col_x(left, group * 10);
        let text = String::from_utf8_lossy(chunk);
        painter.text(
            pos2(x, top),
            Align2::LEFT_TOP,
            text,
            mx.font.clone(),
            t.text,
        );
        if alphabet.is_nucleotide() {
            let comp: String = chunk
                .iter()
                .map(|&b| char::from(alphabet.complement(b).unwrap_or(b'?')))
                .collect();
            painter.text(
                pos2(x, bottom_y),
                Align2::LEFT_TOP,
                comp,
                mx.font.clone(),
                t.text_muted,
            );
        }
    }

    // The features come straight under the two strands, each with its translation in a strip
    // of its own above the bar; the reading frames go under the features (amino acids drawn
    // inside the arrow were hidden by its colour and its name).
    let tracks_top = top + 2.0 * mx.line_h + 2.0;
    let lane_block = mx.lane_h + mx.cds_h;
    let aa_top = tracks_top + mx.feature_lanes as f32 * lane_block + 4.0;
    for (row, &frame) in mx.frames.iter().enumerate() {
        let y = aa_top + row as f32 * AA_ROW_H + AA_ROW_H / 2.0;
        painter.text(
            pos2(gutter_x, y),
            Align2::RIGHT_CENTER,
            if frame > 0 {
                format!("+{frame}")
            } else {
                format!("\u{2212}{}", -frame)
            },
            FontId::monospace(10.0),
            t.text_muted,
        );
        for (mid, aa) in frame_codons(m, frame, start, end, mx.table) {
            let cx = mx.col_x(left, mid - start) + mx.char_w / 2.0;
            let color = match aa {
                b'*' => t.error,
                b'M' => t.accent,
                _ => t.text,
            };
            painter.text(
                pos2(cx, y),
                Align2::CENTER_CENTER,
                (aa as char).to_string(),
                FontId::monospace(12.0),
                color,
            );
        }
    }

    // Feature bars.
    let mut hidden = 0usize;
    let bars = if mx.layers.features {
        layout::line_tracks(m, start, end)
    } else {
        Vec::new()
    };
    for bar in bars {
        if bar.lane >= MAX_LANES {
            hidden += 1;
            continue;
        }
        let color = t.color_at(m, bar.feature);
        // The bar sits below its own translation strip.
        let y = tracks_top + bar.lane as f32 * lane_block + mx.cds_h;
        let x0 = mx.col_x(left, bar.start - start);
        let x1 = mx.col_x(left, bar.end - 1 - start) + mx.char_w;
        let h = mx.lane_h - 4.0;
        let head = (h * 0.8).min((x1 - x0) * 0.5);
        let (b0, b1) = match bar.arrow {
            Arrow::End => (x0, x1 - head),
            Arrow::Start => (x0 + head, x1),
            Arrow::None => (x0, x1),
        };
        if b1 > b0 {
            painter.rect_filled(Rect::from_min_max(pos2(b0, y), pos2(b1, y + h)), 2.0, color);
        }
        let mid = y + h / 2.0;
        let tri = match bar.arrow {
            Arrow::End if head >= 2.0 => Some([
                pos2(x1 - head, y - 1.5),
                pos2(x1, mid),
                pos2(x1 - head, y + h + 1.5),
            ]),
            Arrow::Start if head >= 2.0 => Some([
                pos2(x0 + head, y - 1.5),
                pos2(x0, mid),
                pos2(x0 + head, y + h + 1.5),
            ]),
            _ => None,
        };
        if let Some(tri) = tri {
            painter.add(Shape::convex_polygon(tri.to_vec(), color, Stroke::NONE));
        }
        // The translation: one letter over each codon's middle base, in the strip above the
        // bar, so neither the arrow's colour nor its name can hide it.
        if mx.cds_h > 0.0
            && let Some(list) = codons.get(&bar.feature)
        {
            let lo = list.partition_point(|&(p, _)| p < bar.start);
            let font = FontId::monospace(10.0);
            let aa_y = y - mx.cds_h / 2.0;
            for &(p, aa) in list.iter().skip(lo).take_while(|&&(p, _)| p < bar.end) {
                if p < bar.start {
                    continue;
                }
                let cx = mx.col_x(left, p - start) + mx.char_w / 2.0;
                painter.text(
                    pos2(cx, aa_y),
                    Align2::CENTER_CENTER,
                    (aa as char).to_string(),
                    font.clone(),
                    // Text colour, not the feature's: the strip is bare canvas, where a pale
                    // feature colour would be hard to read.
                    t.text,
                );
            }
        }
        if bar.label
            && let Some(f) = m.features().get(bar.feature)
        {
            let label_w = f.name.chars().count() as f32 * 6.2;
            if label_w + 6.0 < b1 - b0 {
                painter.text(
                    pos2(b0 + 4.0, mid),
                    Align2::LEFT_CENTER,
                    &f.name,
                    FontId::proportional(10.0),
                    t.canvas,
                );
            } else {
                painter.text(
                    pos2(x1 + 4.0, mid),
                    Align2::LEFT_CENTER,
                    &f.name,
                    FontId::proportional(10.0),
                    t.text_muted,
                );
            }
        }
    }
    if hidden > 0 {
        painter.text(
            pos2(gutter_x, tracks_top),
            Align2::RIGHT_TOP,
            format!("+{hidden}"),
            FontId::proportional(10.0),
            t.text_muted,
        );
    }

    // Primers: each in its lane, written base by base against the strand it pairs with.
    let primers_top = aa_top + mx.frames.len() as f32 * AA_ROW_H;
    for d in primers {
        draw_duplex(
            &painter,
            t,
            m,
            mx,
            d,
            left,
            (gutter_x + 6.0, mx.col_x(left, end - start) + 4.0),
            primers_top,
            start,
            end,
        );
    }

    // ORFs: thin arrows, one lane each.
    let orfs_top = primers_top + mx.primer_rows as f32 * mx.line_h + 2.0;
    for o in orfs {
        if o.lane >= mx.orf_rows {
            continue;
        }
        let y = orfs_top + o.lane as f32 * ORF_ROW_H + ORF_ROW_H / 2.0;
        let forward = o.orf.strand == snappiesnap_engine::doc::Strand::Forward;
        for (s, e, has_start, has_end) in pieces(o.orf.range, len, start, end) {
            let x0 = mx.col_x(left, s - start);
            let x1 = mx.col_x(left, e - 1 - start) + mx.char_w;
            painter.line_segment([pos2(x0, y), pos2(x1, y)], Stroke::new(1.6, t.text_muted));
            let tip = if forward && has_end {
                Some((x1, -1.0))
            } else if !forward && has_start {
                Some((x0, 1.0))
            } else {
                None
            };
            if let Some((x, dir)) = tip {
                painter.add(Shape::convex_polygon(
                    vec![
                        pos2(x, y),
                        pos2(x + dir * 5.0, y - 3.0),
                        pos2(x + dir * 5.0, y + 3.0),
                    ],
                    t.text_muted,
                    Stroke::NONE,
                ));
            }
        }
    }
}

/// Most primer lanes per strand in the sequence view; further primers are left out there.
const MAX_PRIMER_LANES: usize = 8;

/// A primer as the sequence view draws it: which position shows which of its bases, the
/// loops standing for parts that aren't written out, and its lane.
pub(crate) struct Duplex<'a> {
    mark: &'a crate::enzymes::PrimerMark,
    /// (position, unwrapped: may be below 0 or past the end; primer base index; paired part)
    cells: Vec<(i64, usize, bool)>,
    /// Loops: (position they attach to, bases they stand for, on its left side).
    loops: Vec<(i64, usize, bool)>,
    /// The position of the primer's 5' base, for its name.
    five: i64,
    /// Row from the top of the primer rows.
    row: usize,
}

/// Lay the primers out: each primer's bases against the positions they pair with (forward
/// primers 5'→3' rightwards; reverse ones written 3'→5', each base over its partner), tails
/// written out or shown as a loop, and lanes so that primers don't overlap. Forward lanes
/// first, then reverse ones. Returns the duplexes and the number of rows.
fn duplexes<'a>(
    marks: &'a [crate::enzymes::PrimerMark],
    tails: TailStyle,
    char_w: &f32,
) -> (Vec<Duplex<'a>>, usize) {
    let full = tails == TailStyle::Full;
    let name_cols = |p: &crate::enzymes::PrimerMark| -> i64 {
        ((p.name.chars().count() as f32 * 5.8) / char_w.max(1.0)).ceil() as i64 + 2
    };
    let mut out: Vec<(Duplex<'a>, (i64, i64))> = Vec::new();
    for p in marks {
        let bases: Vec<u8> = p.seq.bytes().collect();
        let n = bases.len();
        let a = p.annealed.min(n);
        let rs = p.range.start as i64;
        let tail = n - a;
        // Where the primer's base `i` lies, and whether it is in the part that pairs.
        let (five, dir): (i64, i64) = match (p.forward, p.five_prime_only) {
            (true, false) => (rs - tail as i64, 1),
            (true, true) => (rs, 1),
            (false, false) => (rs + a as i64 - 1 + tail as i64, -1),
            (false, true) => (rs + a as i64 - 1, -1),
        };
        let paired = |i: usize| if p.five_prime_only { i < a } else { i >= tail };
        let mut cells = Vec::new();
        for i in 0..n {
            if full || paired(i) {
                cells.push((five + dir * i as i64, i, paired(i)));
            }
        }
        let mut loops = Vec::new();
        if !full && tail > 0 {
            // The part not written out hangs off the end of the paired part it joins.
            let (first_paired, last_paired) = (rs, rs + a as i64 - 1);
            let at_five_end = !p.five_prime_only;
            let left = at_five_end == p.forward;
            let col = if left { first_paired } else { last_paired };
            loops.push((col, tail, left));
        }
        let lo = cells.iter().map(|c| c.0).min().unwrap_or(rs);
        let hi = cells.iter().map(|c| c.0).max().unwrap_or(rs);
        // Room for the name beside the 5' end, and the loop's label.
        let extent = if p.forward {
            (lo - name_cols(p) - 2, hi + 2)
        } else {
            (lo - 2, hi + name_cols(p) + 2)
        };
        out.push((
            Duplex {
                mark: p,
                cells,
                loops,
                five,
                row: 0,
            },
            extent,
        ));
    }
    // Lanes: greedy by start, per strand.
    out.sort_by_key(|(d, e)| (!d.mark.forward, e.0));
    let mut ends: [Vec<i64>; 2] = [Vec::new(), Vec::new()];
    let mut kept = Vec::new();
    for (mut d, (s, e)) in out {
        let k = usize::from(!d.mark.forward);
        let Some(lanes) = ends.get_mut(k) else {
            continue;
        };
        let lane = match lanes.iter().position(|&last| s > last) {
            Some(l) => l,
            None if lanes.len() < MAX_PRIMER_LANES => {
                lanes.push(i64::MIN);
                lanes.len() - 1
            }
            None => continue,
        };
        if let Some(last) = lanes.get_mut(lane) {
            *last = e;
        }
        d.row = lane;
        kept.push((d, k));
    }
    let forward_rows = ends[0].len();
    for (d, k) in &mut kept {
        if *k == 1 {
            d.row += forward_rows;
        }
    }
    let rows = forward_rows + ends[1].len();
    (kept.into_iter().map(|(d, _)| d).collect(), rows)
}

/// The column (0-based on this line) of unwrapped position `pos`, if it is on the line.
fn line_col(pos: i64, len: usize, circular: bool, start: usize, end: usize) -> Option<usize> {
    let l = i64::try_from(len).ok()?;
    let p = if circular {
        pos.rem_euclid(l.max(1))
    } else {
        pos
    };
    let p = usize::try_from(p).ok()?;
    (p >= start && p < end).then(|| p - start)
}

#[allow(clippy::too_many_arguments)]
fn draw_duplex(
    painter: &egui::Painter,
    t: &Tokens,
    m: &Molecule,
    mx: &Metrics,
    d: &Duplex<'_>,
    left: f32,
    // Left and right edge of the residues drawn, to keep a primer's name in view.
    bounds: (f32, f32),
    top: f32,
    start: usize,
    end: usize,
) {
    let len = m.len();
    let circular = m.topology() == snappiesnap_engine::doc::Topology::Circular;
    let seq = m.seq().as_bytes();
    let p = d.mark;
    let y = top + d.row as f32 * mx.line_h + mx.line_h / 2.0;
    let paired_color = if p.primary { t.text } else { t.text_muted };
    let complement = |b: u8| match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    };
    let primer = p.seq.as_bytes();
    let mut last_x: Option<(f32, usize)> = None;
    for &(pos, i, paired) in &d.cells {
        let Some(col) = line_col(pos, len, circular, start, end) else {
            continue;
        };
        let Some(&b) = primer.get(i) else { continue };
        let x = mx.col_x(left, col);
        let up = b.to_ascii_uppercase();
        // The template base this primer base would pair with (bottom strand for reverse).
        let partner = seq.get(col + start).map(|&tb| {
            let tb = tb.to_ascii_uppercase();
            if p.forward { tb } else { complement(tb) }
        });
        let matches = partner.is_some_and(|tb| tb == up || up == b'N');
        let (text, color) = if !paired {
            (
                char::from(b.to_ascii_lowercase()),
                t.text_muted.gamma_multiply(0.8),
            )
        } else if matches {
            (char::from(up), paired_color)
        } else {
            (char::from(up), t.error)
        };
        // A base that doesn't pair steps up out of the row, with the row stepping around it:
        // the primer still reads as one piece, and where it differs is obvious at a glance.
        // A repair or mutagenesis oligo is mostly made of such bases.
        let jog = if paired && !matches {
            MISMATCH_JOG
        } else {
            0.0
        };
        if jog > 0.0 {
            painter.rect_filled(
                Rect::from_min_size(
                    pos2(x - 0.5, y - jog - mx.line_h / 2.0 + 1.0),
                    vec2(mx.char_w + 1.0, mx.line_h - 2.0),
                ),
                2.0,
                t.error.gamma_multiply(0.15),
            );
            let step = Stroke::new(1.0, t.error.gamma_multiply(0.7));
            let (bottom, top) = (y + mx.line_h / 2.0 - 2.0, y - jog + mx.line_h / 2.0 - 2.0);
            for sx in [x - 0.5, x + mx.char_w + 0.5] {
                painter.line_segment([pos2(sx, bottom), pos2(sx, top)], step);
            }
        }
        painter.text(
            pos2(x + mx.char_w / 2.0, y - jog),
            Align2::CENTER_CENTER,
            text,
            mx.font.clone(),
            color,
        );
        if primer.len().checked_sub(1) == Some(i) {
            last_x = Some((x, col));
        }
    }
    // The 3' end: a small arrow pointing the way the primer extends.
    if let Some((x, _)) = last_x {
        let (tip, base) = if p.forward {
            (x + mx.char_w + 6.0, x + mx.char_w + 1.0)
        } else {
            (x - 6.0, x - 1.0)
        };
        painter.add(egui::Shape::convex_polygon(
            vec![pos2(tip, y), pos2(base, y - 4.0), pos2(base, y + 4.0)],
            paired_color,
            Stroke::NONE,
        ));
    }
    // Loops for the parts not written out, with how many bases they stand for.
    for &(pos, count, on_left) in &d.loops {
        let Some(col) = line_col(pos, len, circular, start, end) else {
            continue;
        };
        let x = if on_left {
            mx.col_x(left, col) - 1.0
        } else {
            mx.col_x(left, col) + mx.char_w + 1.0
        };
        let dx = if on_left { -7.0 } else { 7.0 };
        let top_y = y - mx.line_h / 2.0 + 1.0;
        painter.add(egui::Shape::CubicBezier(
            egui::epaint::CubicBezierShape::from_points_stroke(
                [
                    pos2(x, y + 3.0),
                    pos2(x + dx, y + 3.0),
                    pos2(x + dx, top_y),
                    pos2(x, top_y + 2.0),
                ],
                false,
                egui::Color32::TRANSPARENT,
                Stroke::new(1.2, t.text_muted),
            ),
        ));
        painter.text(
            pos2(x + dx + if on_left { -2.0 } else { 2.0 }, y),
            if on_left {
                Align2::RIGHT_CENTER
            } else {
                Align2::LEFT_CENTER
            },
            format!("+{count}"),
            FontId::proportional(9.5),
            t.text_muted,
        );
    }
    // The name beside the 5' end (beyond its loop, when it has one at that end); where the
    // line leaves no room there, beside the primer's other end on this line.
    let name_w = p.name.chars().count() as f32 * 5.8;
    let drawn: Vec<f32> = d
        .cells
        .iter()
        .filter_map(|c| line_col(c.0, len, circular, start, end))
        .map(|col| mx.col_x(left, col))
        .collect();
    let (line_left, line_right) = bounds;
    if let Some(col) = line_col(d.five_name_at(), len, circular, start, end) {
        let loop_here = d.loops.iter().any(|l| l.0 == d.five_name_at());
        let shift = if loop_here { 30.0 } else { 4.0 };
        let (mut x, mut align) = if p.forward {
            (mx.col_x(left, col) - shift, Align2::RIGHT_CENTER)
        } else {
            (mx.col_x(left, col) + mx.char_w + shift, Align2::LEFT_CENTER)
        };
        let fits = if p.forward {
            x - name_w >= line_left
        } else {
            x + name_w <= line_right
        };
        if !fits {
            let lo = drawn.iter().copied().fold(f32::INFINITY, f32::min);
            let hi = drawn.iter().copied().fold(f32::NEG_INFINITY, f32::max);
            (x, align) = if p.forward {
                (hi + mx.char_w + 12.0, Align2::LEFT_CENTER)
            } else {
                (lo - 12.0, Align2::RIGHT_CENTER)
            };
        }
        if x.is_finite() {
            painter.text(
                pos2(x, y),
                align,
                &p.name,
                FontId::proportional(9.5),
                t.text_muted,
            );
        }
    }
}

impl Duplex<'_> {
    /// The position the name sits beside: the outermost drawn base at the 5' end.
    fn five_name_at(&self) -> i64 {
        let drawn_five = self.cells.iter().any(|c| c.0 == self.five);
        if drawn_five {
            self.five
        } else if self.mark.forward {
            self.cells.iter().map(|c| c.0).min().unwrap_or(self.five)
        } else {
            self.cells.iter().map(|c| c.0).max().unwrap_or(self.five)
        }
    }
}

/// The parts of `r` (which may wrap the origin) inside the line `start..end`: (from, to, whether
/// it holds the range's start, whether it holds its end).
fn pieces(r: Range, len: usize, start: usize, end: usize) -> Vec<(usize, usize, bool, bool)> {
    let mut out = Vec::new();
    if len == 0 {
        return out;
    }
    let (a, b) = r.segments(len);
    let parts: Vec<std::ops::Range<usize>> = std::iter::once(a).chain(b).collect();
    let n = parts.len();
    for (k, p) in parts.into_iter().enumerate() {
        let (s, e) = (p.start.max(start), p.end.min(end));
        if s < e {
            out.push((s, e, k == 0 && s == p.start, k + 1 == n && e == p.end));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn metrics() -> Metrics {
        Metrics {
            sites_h: 0.0,
            frames: Vec::new(),
            table: 1,
            layers: crate::map::Layers::default(),
            cds_h: 0.0,
            feature_lanes: 1,
            primer_rows: 0,
            orf_rows: 0,
            font: FontId::monospace(13.0),
            char_w: 8.0,
            line_h: 17.0,
            gap: 6.0,
            lane_h: 13.0,
            number_w: 40.0,
            per_line: 60,
            row_h: 100.0,
            single: false,
            ruler_h: 0.0,
        }
    }

    fn mark(
        name: &str,
        seq: &str,
        at: usize,
        annealed: usize,
        forward: bool,
        five: bool,
    ) -> crate::enzymes::PrimerMark {
        crate::enzymes::PrimerMark {
            primer: 0,
            name: name.into(),
            range: Range::new(at, at + annealed),
            forward,
            primary: true,
            annealed,
            mismatches: 0,
            five_prime_only: five,
            seq: seq.into(),
            lane: 0,
            window: None,
            tm_site: None,
            tm_full: None,
            hairpin_dg: None,
            dimer_dg: None,
            checks: Vec::new(),
        }
    }

    #[test]
    fn duplex_layout() {
        // A forward primer with a 4-base tail annealing at 100..110, a reverse one at 300..308
        // with a 2-base tail, and one whose 5' 6 bases alone pair at 500..506.
        let marks = vec![
            mark("f", "ttttACGTACGTAC", 100, 10, true, false),
            mark("r", "ggCCGGCCGG", 300, 8, false, false),
            mark("p", "AAAAAAcccc", 500, 6, true, true),
        ];
        let (d, rows) = duplexes(&marks, TailStyle::Simplified, &8.0);
        assert_eq!(rows, 2, "one forward lane, one reverse");
        let f = d.iter().find(|d| d.mark.name == "f").unwrap();
        assert_eq!(f.cells.len(), 10);
        assert_eq!(f.cells.first(), Some(&(100, 4, true)));
        assert_eq!(f.cells.last(), Some(&(109, 13, true)));
        assert_eq!(f.loops, vec![(100, 4, true)]);
        assert_eq!(f.five, 96);
        assert_eq!(f.five_name_at(), 100);
        let r = d.iter().find(|d| d.mark.name == "r").unwrap();
        assert_eq!(r.row, 1);
        // Written 3'→5': its first paired base (index 2) over the last template base.
        assert_eq!(r.cells.first(), Some(&(307, 2, true)));
        assert_eq!(r.cells.last(), Some(&(300, 9, true)));
        assert_eq!(r.loops, vec![(307, 2, false)]);
        let p = d.iter().find(|d| d.mark.name == "p").unwrap();
        assert_eq!(p.cells.len(), 6);
        assert_eq!(p.loops, vec![(505, 4, false)]);
        assert_eq!(p.row, 0);

        // Written out, the tails take their own positions, unpaired.
        let (d, _) = duplexes(&marks, TailStyle::Full, &8.0);
        let f = d.iter().find(|d| d.mark.name == "f").unwrap();
        assert_eq!(f.cells.first(), Some(&(96, 0, false)));
        assert!(f.loops.is_empty());
        let r = d.iter().find(|d| d.mark.name == "r").unwrap();
        assert_eq!(r.cells.first(), Some(&(309, 0, false)));
        let p = d.iter().find(|d| d.mark.name == "p").unwrap();
        assert_eq!(p.cells.last(), Some(&(509, 9, false)));

        // Overlapping primers on one strand stack, up to the lane limit.
        let many: Vec<_> = (0..12)
            .map(|i| mark("x", "ACGTACGTACGT", 10 + i, 12, true, false))
            .collect();
        let (d, rows) = duplexes(&many, TailStyle::Simplified, &8.0);
        assert_eq!(rows, MAX_PRIMER_LANES);
        assert_eq!(d.len(), MAX_PRIMER_LANES);
        assert_eq!(duplexes(&[], TailStyle::Full, &8.0).1, 0);
    }

    #[test]
    fn columns_round_trip() {
        let m = metrics();
        for col in [0, 1, 9, 10, 11, 35, 59, 60] {
            let x = m.col_x(0.0, col);
            assert_eq!(m.col_at(0.0, x), col, "col {col}");
            assert_eq!(m.col_at(0.0, x + 3.0), col, "col {col} + 3");
        }
        assert_eq!(m.col_at(0.0, -100.0), 0);
        assert_eq!(m.col_at(0.0, 10_000.0), 60);
    }

    /// Drive `keyboard` with real egui input events.
    fn keys(events: Vec<Event>, selection: Option<Range>) -> Vec<SeqAction> {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, "ACGTACGT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let ctx = egui::Context::default();
        let mut out = Vec::new();
        let input = egui::RawInput {
            events,
            ..Default::default()
        };
        let mut full = ctx.run_ui(input, |ui| keyboard(ui, &m, selection, &mut out));
        full.textures_delta.clear();
        out
    }

    fn key(key: Key) -> Event {
        Event::Key {
            key,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        }
    }

    #[test]
    fn typing_and_deleting() {
        assert_eq!(
            keys(vec![Event::Text("gx t".into())], Some(Range::at(2))),
            vec![SeqAction::Insert {
                at: 2,
                residues: "GT".into()
            }]
        );
        assert_eq!(
            keys(vec![Event::Text("A".into())], Some(Range::new(1, 3))),
            vec![SeqAction::Replace {
                range: Range::new(1, 3),
                residues: "A".into()
            }]
        );
        assert_eq!(
            keys(vec![key(Key::Backspace)], Some(Range::at(3))),
            vec![SeqAction::Delete(Range::new(2, 3))]
        );
        assert!(keys(vec![key(Key::Backspace)], Some(Range::at(0))).is_empty());
        assert_eq!(
            keys(vec![key(Key::Delete)], Some(Range::new(2, 5))),
            vec![SeqAction::Delete(Range::new(2, 5))]
        );
        assert!(keys(vec![key(Key::Delete)], Some(Range::at(8))).is_empty());
        assert_eq!(
            keys(vec![key(Key::ArrowRight)], Some(Range::at(8))),
            vec![SeqAction::Select(Some(Range::at(8)))]
        );
        // Without a selection, typing appends.
        assert_eq!(
            keys(vec![Event::Text("c".into())], None),
            vec![SeqAction::Insert {
                at: 8,
                residues: "C".into()
            }]
        );
    }

    #[test]
    fn cds_translations_on_both_strands() {
        use snappiesnap_engine::doc::{Alphabet, NewFeature, Segment, Seq, Strand, Topology};
        // ATG AAA CCC TAA forward at 2..14; its reverse complement at 16..28 on the bottom.
        let s = "CCATGAAACCCTAACCTTAGGGTTTCATCC";
        let mut m =
            Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Linear).unwrap();
        for (range, strand) in [
            (Range::new(2, 14), Strand::Forward),
            (Range::new(16, 28), Strand::Reverse),
        ] {
            m.add_feature(NewFeature {
                color: None,
                name: "cds".into(),
                kind: "CDS".into(),
                strand,
                segments: vec![Segment { range }],
                qualifiers: vec![],
            })
            .unwrap();
        }
        let c = cds_codons(&m);
        assert_eq!(c[&0], vec![(3, b'M'), (6, b'K'), (9, b'P'), (12, b'*')]);
        // Reverse: the first codon (Met) is at the right end, its middle base at 26.
        assert_eq!(c[&1], vec![(17, b'*'), (20, b'P'), (23, b'K'), (26, b'M')]);
    }

    /// Reading frames by hand: ATGGCCTAA... on the top strand, and the bottom strand read from
    /// the other end.
    #[test]
    fn reading_frames() {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        // 0-based: A0 T1 G2 G3 C4 C5 T6 A7 A8 G9
        let m = Molecule::new(
            "f",
            Seq::new(Alphabet::Dna, "ATGGCCTAAG").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let all = |frame| frame_codons(&m, frame, 0, m.len(), 1);
        // +1: ATG GCC TAA (middles 1, 4, 7).
        assert_eq!(all(1), vec![(1, b'M'), (4, b'A'), (7, b'*')]);
        // +2: TGG CCT AAG.
        assert_eq!(all(2), vec![(2, b'W'), (5, b'P'), (8, b'K')]);
        // +3: GGC CTA.
        assert_eq!(all(3), vec![(3, b'G'), (6, b'L')]);
        // Bottom strand (reverse complement CTTAGGCCAT), read from position 9 down:
        // -1: CTT AGG CCA (top 9..7, 6..4, 3..1) -> L R P, middles 8, 5, 2.
        assert_eq!(all(-1), vec![(2, b'P'), (5, b'R'), (8, b'L')]);
        // -2: TTA GGC CAT (top 8..6, 5..3, 2..0) -> L G H.
        assert_eq!(all(-2), vec![(1, b'H'), (4, b'G'), (7, b'L')]);
        // -3: TAG GCC (top 7..5, 4..2) -> * A.
        assert_eq!(all(-3), vec![(3, b'A'), (6, b'*')]);
        // Only codons whose middle base is in the line.
        assert_eq!(frame_codons(&m, 1, 3, 6, 1), vec![(4, b'A')]);
        // Too short, or not a frame: nothing.
        let short = Molecule::new(
            "s",
            Seq::new(Alphabet::Dna, "AT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert!(frame_codons(&short, 1, 0, 2, 1).is_empty());
        assert!(frame_codons(&short, -1, 0, 2, 1).is_empty());
        assert!(frame_codons(&m, 4, 0, 10, 1).is_empty());
        assert!(frame_codons(&m, 0, 0, 10, 1).is_empty());
    }

    /// Ranges (also across the origin) split over lines, knowing which piece holds each end.
    #[test]
    fn ranges_split_over_lines() {
        // 100 bp; a range 90..10 across the origin; lines of 60.
        let r = Range::new(90, 10);
        assert_eq!(pieces(r, 100, 60, 100), vec![(90, 100, true, false)]);
        assert_eq!(pieces(r, 100, 0, 60), vec![(0, 10, false, true)]);
        // Inside one line: both ends.
        assert_eq!(
            pieces(Range::new(5, 20), 100, 0, 60),
            vec![(5, 20, true, true)]
        );
        // Across a line break.
        assert_eq!(
            pieces(Range::new(50, 70), 100, 0, 60),
            vec![(50, 60, true, false)]
        );
        assert_eq!(
            pieces(Range::new(50, 70), 100, 60, 100),
            vec![(60, 70, false, true)]
        );
        // Elsewhere, or an empty molecule: nothing.
        assert!(pieces(Range::new(5, 20), 100, 60, 100).is_empty());
        assert!(pieces(Range::new(0, 0), 0, 0, 0).is_empty());
    }

    #[test]
    fn copy_cut_and_paste_keys() {
        let sel = Some(Range::new(2, 5));
        assert_eq!(keys(vec![Event::Copy], sel), vec![SeqAction::Copy("plain")]);
        assert_eq!(keys(vec![Event::Cut], sel), vec![SeqAction::Cut]);
        // Nothing selected: nothing to copy or cut.
        assert!(keys(vec![Event::Copy, Event::Cut], Some(Range::at(3))).is_empty());
        assert_eq!(
            keys(vec![Event::Paste(">x\nACGT".into())], Some(Range::at(3))),
            vec![SeqAction::Paste(">x\nACGT".into())]
        );
    }

    /// A primer base that doesn't pair steps up out of its row, so a repair or mutagenesis
    /// oligo still reads as one piece. Set `SNAPPIESNAP_JOG_PNG=/path.png` to look at it.
    #[test]
    fn a_base_that_does_not_pair_steps_out_of_the_row() {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        let template = "GCCATCGTCCCCATTCTGCTGGATGCCAATGTCAGCACTTATGACAAAATCTGCATCATCC";
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, template).expect("a template"),
            Topology::Linear,
        )
        .expect("a molecule");
        // The same oligo as the bio test, cut to this stretch: two changes in it.
        let oligo = "CCATCGTCCCCATTCTGCTGAATGCAAATGTCAGCACTTATGACAAAATCCGCATCATCC";
        let marks = vec![mark("Repair", oligo, 1, oligo.len(), true, false)];
        let t = Tokens::light();
        let mx = metrics();
        let (duplexes, rows) = duplexes(&marks, TailStyle::Simplified, &mx.char_w);
        assert_eq!(rows, 1);
        let d = duplexes.first().expect("the oligo is drawn");
        let image = crate::offscreen::render(
            700,
            120,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let painter = ui.painter().clone();
                draw_duplex(
                    &painter,
                    &t,
                    &m,
                    &mx,
                    d,
                    10.0,
                    (10.0, 690.0),
                    40.0,
                    0,
                    m.len(),
                );
            },
        )
        .expect("the primer draws");
        if let Ok(path) = std::env::var("SNAPPIESNAP_JOG_PNG") {
            crate::offscreen::write_png(&path, &image).expect("writing the png");
        }
        // The two changed bases are drawn above the row the rest of the oligo sits on: count
        // the error-coloured pixels above and below the row's middle.
        let [w, _] = image.size;
        let err = t.error;
        let close = |a: u8, b: u8| a.abs_diff(b) < 40;
        let marked: Vec<usize> = image
            .pixels
            .iter()
            .enumerate()
            .filter(|(_, px)| {
                close(px.r(), err.r()) && close(px.g(), err.g()) && close(px.b(), err.b())
            })
            .map(|(i, _)| i / w.max(1))
            .collect();
        assert!(marked.len() > 20, "the changed bases should be drawn");
        // The row's middle, as `draw_duplex` works it out from the `top` passed in. Without
        // the step up, nothing marked could reach higher than the top of the row itself.
        let row_middle = 40.0 + mx.line_h / 2.0;
        let row_top = row_middle - mx.line_h / 2.0;
        let highest = marked.iter().copied().min().unwrap_or(usize::MAX) as f32;
        assert!(
            highest < row_top,
            "a changed base should reach above the row's top ({row_top}), not stop at {highest}"
        );
        let mean = marked.iter().sum::<usize>() as f32 / marked.len() as f32;
        assert!(
            mean < row_middle,
            "the changed bases sit above the row's middle ({row_middle}), not at {mean}"
        );
    }

    /// One line, scrolled sideways: only the residues in view are drawn, whatever the
    /// molecule's length, and a selection that lands off screen brings itself into view. Set
    /// `SNAPPIESNAP_SINGLE_PNG=/path.png` to look at it.
    #[test]
    fn the_whole_sequence_on_one_line_draws_only_what_is_in_view() {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        let m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, "GAATTCAAATGGTCCCACCA".repeat(2000)).expect("a sequence"),
            Topology::Circular,
        )
        .expect("a molecule");
        assert_eq!(m.len(), 40_000);
        let t = Tokens::dark();
        let marks = crate::enzymes::Marks {
            sites: &[],
            primers: &[],
            orfs: &[],
        };
        let mut state = SeqViewState::default();
        let draw = |state: &mut SeqViewState, selection: Option<Range>| {
            crate::offscreen::render(
                1200,
                300,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    show(
                        ui,
                        &t,
                        &m,
                        selection,
                        state,
                        1,
                        &marks,
                        crate::map::Layers::default(),
                        &[],
                        1,
                        &[],
                        TailStyle::Simplified,
                        true,
                    );
                },
            )
            .expect("the view draws")
        };
        let image = draw(&mut state, None);
        let (first, last) = state.drawn.expect("a window was drawn");
        assert_eq!(first, 0);
        assert!(
            last < 400,
            "a 40 kb molecule should cost a screenful, not {last} residues"
        );
        if let Ok(path) = std::env::var("SNAPPIESNAP_SINGLE_PNG") {
            crate::offscreen::write_png(&path, &image).expect("writing the png");
        }
        // A selection far away brings itself into view, and only then.
        draw(&mut state, Some(Range::new(30_000, 30_020)));
        let (first, last) = state.drawn.expect("a window was drawn");
        assert!(
            (first..last).contains(&30_000),
            "the selection should be in view, not {first}..{last}"
        );
    }

    /// When a selection pulls the view to it, and when it leaves the view alone. (The render
    /// test above cannot check this: each offscreen frame is a context of its own, so a scroll
    /// area never remembers where it was.)
    #[test]
    fn only_a_selection_that_moved_out_of_view_scrolls() {
        let drawn = Some((100, 200));
        // Moved, and off screen: scroll to it.
        assert_eq!(bring_into_view(Some(900), Some(150), drawn), Some(900));
        assert_eq!(bring_into_view(Some(5), Some(150), drawn), Some(5));
        // Moved, but already on screen (a click, a drag): leave the view alone.
        assert_eq!(bring_into_view(Some(150), Some(120), drawn), None);
        assert_eq!(bring_into_view(Some(199), Some(120), drawn), None);
        // The end of the window is past the last residue drawn.
        assert_eq!(bring_into_view(Some(200), Some(120), drawn), Some(200));
        // Unchanged, or nothing selected: nothing to do.
        assert_eq!(bring_into_view(Some(900), Some(900), drawn), None);
        assert_eq!(bring_into_view(None, Some(150), drawn), None);
        // Before anything has been drawn, any selection is worth scrolling to.
        assert_eq!(bring_into_view(Some(7), None, None), Some(7));
    }

    /// Find matches are tinted, and the view scrolls to a requested position. Set
    /// `SNAPPIESNAP_FIND_PNG=/path.png` to look at it.
    #[test]
    fn find_matches_draw() {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        let m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, "GAATTCAAATGGTCCCACCA".repeat(30)).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        let found = [
            Range::new(300, 306),
            Range::new(320, 326),
            Range::new(595, 5),
        ];
        let t = Tokens::light();
        let mut state = SeqViewState {
            scroll_to: Some(300),
            ..Default::default()
        };
        let marks = crate::enzymes::Marks {
            sites: &[],
            primers: &[],
            orfs: &[],
        };
        let image = crate::offscreen::render(
            900,
            400,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let _ = show(
                    ui,
                    &t,
                    &m,
                    Some(Range::new(300, 306)),
                    &mut state,
                    1,
                    &marks,
                    crate::map::Layers::default(),
                    &[],
                    1,
                    &found,
                    TailStyle::Simplified,
                    false,
                );
            },
        )
        .unwrap();
        assert!(state.scroll_to.is_none(), "the scroll request is used once");
        if let Ok(path) = std::env::var("SNAPPIESNAP_FIND_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
    }

    /// A double-click on an enzyme's name is resolved against the rectangles `draw_sites`
    /// records: each drawn name must point back at its own site, sit above the strands, and a
    /// site off this line must record nothing.
    #[test]
    fn site_names_record_where_they_were_drawn() {
        use crate::enzymes::map_sites;
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        use snappiesnap_engine::enzymes::{EnzymeDb, SiteQuery};
        // EcoRI at 10, BamHI at 30, HindIII at 90 (past the 60-base line drawn).
        let s = format!(
            "{}GAATTC{}GGATCC{}AAGCTT{}",
            "T".repeat(10),
            "T".repeat(14),
            "T".repeat(54),
            "T".repeat(20)
        );
        let m = Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Linear).unwrap();
        let query = SiteQuery {
            enzymes: Some(vec!["EcoRI".into(), "BamHI".into(), "HindIII".into()]),
            ..SiteQuery::default()
        };
        let sites = map_sites(&EnzymeDb::default(), &m, &query);
        assert_eq!(sites.len(), 3, "three sites to label");
        let t = Tokens::light();
        let mx = Metrics {
            sites_h: SITE_ROWS as f32 * SITE_ROW_H + 4.0,
            ..metrics()
        };
        let mut labels: Vec<(Rect, usize)> = Vec::new();
        let top = 60.0;
        crate::offscreen::render(
            900,
            200,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                let painter = ui.painter().clone();
                draw_sites(
                    &painter,
                    &t,
                    &mx,
                    &sites,
                    10.0,
                    10.0 + mx.number_w - 10.0,
                    top,
                    0,
                    60,
                    m.len(),
                    false,
                    None,
                    &mut labels,
                );
            },
        )
        .unwrap();
        assert!(!labels.is_empty(), "names were drawn");
        let drawn: Vec<&str> = labels
            .iter()
            .filter_map(|(_, i)| sites.get(*i).map(|s| s.mark.text.as_str()))
            .collect();
        assert!(drawn.iter().any(|n| n.starts_with("EcoRI")), "{drawn:?}");
        assert!(drawn.iter().any(|n| n.starts_with("BamHI")), "{drawn:?}");
        assert!(
            !drawn.iter().any(|n| n.starts_with("HindIII")),
            "a site on another line records nothing: {drawn:?}"
        );
        for (r, _) in &labels {
            assert!(r.bottom() <= top, "names sit above the strands: {r:?}");
        }
        // EcoRI is left of BamHI, so a point inside its rectangle finds EcoRI.
        let eco = labels
            .iter()
            .find(|(_, i)| {
                sites
                    .get(*i)
                    .is_some_and(|s| s.mark.text.starts_with("EcoRI"))
            })
            .unwrap();
        let hit = labels
            .iter()
            .find(|(r, _)| r.contains(eco.0.center()))
            .unwrap();
        assert_eq!(hit.1, eco.1);

        // Hovering EcoRI's name draws its cut (and name) in the highlight colour; with the
        // pointer elsewhere that colour appears nowhere.
        let draw = |hover: Option<Pos2>| {
            let mut scratch = Vec::new();
            crate::offscreen::render(
                900,
                200,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    let painter = ui.painter().clone();
                    draw_sites(
                        &painter,
                        &t,
                        &mx,
                        &sites,
                        10.0,
                        10.0 + mx.number_w - 10.0,
                        top,
                        0,
                        60,
                        m.len(),
                        false,
                        hover,
                        &mut scratch,
                    );
                },
            )
            .unwrap()
        };
        let highlighted =
            |image: &egui::ColorImage| image.pixels.iter().filter(|&&p| p == t.highlight).count();
        let plain = draw(Some(pos2(5.0, 190.0)));
        let hovered = draw(Some(eco.0.center()));
        assert_eq!(highlighted(&plain), 0, "nothing is highlighted unhovered");
        assert!(
            highlighted(&hovered) > 20,
            "the hovered cut is drawn thick in the highlight colour: {} px",
            highlighted(&hovered)
        );
        if let Ok(path) = std::env::var("SNAPPIESNAP_SITE_HOVER_PNG") {
            crate::offscreen::write_png(&path, &hovered).unwrap();
        }
    }
}
