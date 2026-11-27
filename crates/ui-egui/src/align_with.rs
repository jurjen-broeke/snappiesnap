//! "Align with": sequencing reads and other sequences aligned under the open document, in the
//! sequence view. A panel at the top left lists them (ticked ones are shown; add them with a
//! file picker, from the open tabs, or by dropping files on it; the document's stored reads
//! are listed by themselves). A pane under the sequence lays every ticked one out against the
//! document, column by column, with its differences marked and its trace drawn under its
//! bases (fold a trace away to keep one line of bases). The pane draws on the columns of the
//! one-line sequence view above it, so the two scroll as one; the navigator at the bottom of
//! the view shows the whole construct, where each read lies on it and the stretch in view.
//!
//! The alignment is the `check_reads` command's (both strands, poor ends trimmed, the read
//! kept if at least 80% identical), so agents see the same answer.

use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use eframe::egui::{
    self, Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Shape, Stroke, Ui, pos2, vec2,
};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::{Channel, Molecule};

use crate::seq_navigator::ReadSpan;
use crate::seqview::LineGeometry;
use crate::state::{AlignEntry, AlignSource};
use crate::theme::Tokens;

/// The narrowest the names column of the pane gets.
const MIN_NAME_W: f32 = 120.0;
/// Height of a trace under its read, when there is room.
const TRACE_H: f32 = 76.0;
/// The shortest a trace is squeezed to before the reads scroll instead.
const MIN_TRACE_H: f32 = 36.0;
/// Width of the panel listing the sequences.
pub const PANEL_W: f32 = 232.0;
/// A call below this quality is drawn faint.
const FAINT_BELOW: u64 = 20;

/// What the user asked for.
#[derive(Debug, Clone, PartialEq)]
pub enum AlignWithAction {
    /// Pick sequence files to add to the list.
    AddFiles,
    /// Run `check_reads` with these params; give the answer to [`AlignWithView::set_result`].
    Check(Value),
    /// Select this base of the document (and bring it into view above).
    Select(usize),
    /// Close the panel and the pane.
    Close,
}

/// What the pane remembers between frames (not saved: it is all derived again).
#[derive(Default)]
pub struct AlignWithView {
    key: Option<String>,
    result: Option<Result<Value, String>>,
    /// Read files as molecules (for their traces), loaded once.
    molecules: HashMap<String, Option<Molecule>>,
    /// Highest signal of each read's trace, by source key, for scaling.
    peaks: HashMap<String, f32>,
    /// The reads stored with a document, read once per document.
    matched: Option<(u64, Vec<snappiesnap_engine::io::MatchedRead>)>,
    /// The construct's sequence, hashed, at a revision: the alignment depends on the bases
    /// only, and a selection changes the revision too.
    residues_hash: Option<(u64, u64, u64)>,
    /// First construct column in view (fractional, so scrolling is smooth).
    left: f64,
    /// Columns that fitted last frame.
    shown: usize,
    /// The list row picked for Move and Remove.
    picked: Option<usize>,
    /// Where the panel is on screen, so a file dropped on it joins the list.
    pub panel_rect: Option<Rect>,
    /// Where the pane was drawn last frame.
    pub pane_rect: Option<Rect>,
    /// Points the sequence view should scroll sideways: the user swiped in the pane.
    pub scroll_by: f32,
}

/// A stable key for a source: files by path, tabs by id, stored reads by index.
fn key_of(s: &AlignSource) -> String {
    match s {
        AlignSource::File(p) => format!("file:{p}"),
        AlignSource::Document(id) => format!("doc:{id}"),
        AlignSource::Matched(i) => format!("matched:{i}"),
    }
}

/// Which source a read of a `check_reads` answer came from.
fn source_of(read: &Value) -> Option<AlignSource> {
    let s = read.get("source")?;
    if let Some(p) = s.get("file").and_then(Value::as_str) {
        return Some(AlignSource::File(p.to_string()));
    }
    if let Some(id) = s.get("document").and_then(Value::as_u64) {
        return Some(AlignSource::Document(id));
    }
    let i = s.get("matched").and_then(Value::as_u64)?;
    Some(AlignSource::Matched(usize::try_from(i).ok()?))
}

fn count(v: &Value, k: &str) -> u64 {
    v.get(k).and_then(Value::as_u64).unwrap_or(0)
}

fn to_usize(v: Option<&Value>) -> Option<usize> {
    usize::try_from(v?.as_u64()?).ok()
}

/// 4124 as "4,124".
fn grouped(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::new();
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// `name` cut to `n` characters, with an ellipsis when cut.
fn short(name: &str, n: usize) -> String {
    if name.chars().count() <= n {
        return name.to_string();
    }
    let mut s: String = name.chars().take(n.saturating_sub(1)).collect();
    s.push('\u{2026}');
    s
}

impl AlignWithView {
    /// The answer of `check_reads` for the params the pane asked for.
    pub fn set_result(&mut self, answer: Result<Value, String>) {
        self.result = Some(answer);
    }

    /// The reads stored with document `id` (read from its file once).
    fn matched_reads(
        &mut self,
        session: &Session,
        id: u64,
    ) -> &[snappiesnap_engine::io::MatchedRead] {
        if self.matched.as_ref().is_none_or(|(d, _)| *d != id) {
            let reads =
                snappiesnap_engine::matched::matched_reads(session, Some(id)).unwrap_or_default();
            self.matched = Some((id, reads));
        }
        self.matched.as_ref().map_or(&[], |(_, r)| r.as_slice())
    }

    /// Load what drawing `sources` needs: read files and their trace heights.
    fn load(&mut self, session: &Session, id: u64, sources: &[AlignSource]) {
        for s in sources {
            if let AlignSource::File(path) = s
                && !self.molecules.contains_key(path)
            {
                let m = session.platform().read(path).ok().and_then(|bytes| {
                    snappiesnap_engine::io::import(path, &bytes)
                        .ok()
                        .map(|i| i.molecule)
                });
                self.molecules.insert(path.clone(), m);
            }
            let key = key_of(s);
            if !self.peaks.contains_key(&key) {
                let peak = molecule(self, session, id, s).and_then(call_height);
                if let Some(p) = peak {
                    self.peaks.insert(key, p);
                }
            }
        }
    }

    /// A hash of document `id`'s bases, worked out again only when the session changed.
    fn residues_hash(&mut self, session: &Session, id: u64) -> u64 {
        let revision = session.revision();
        if let Some((d, r, h)) = self.residues_hash
            && d == id
            && r == revision
        {
            return h;
        }
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        if let Some(d) = session.document_by_id(id) {
            d.molecule.seq().as_bytes().hash(&mut hasher);
            d.molecule.topology().hash(&mut hasher);
        }
        let h = hasher.finish();
        self.residues_hash = Some((id, revision, h));
        h
    }
}

/// The height a trace is drawn to: its called peaks, not its highest point. Dye blobs and the
/// first bases of a read can stand ten times higher than the rest and would flatten it. The
/// 90th percentile of the called bases' own peak heights, with some headroom.
fn call_height(m: &Molecule) -> Option<f32> {
    let tr = m.trace()?;
    let residues = m.seq().as_bytes();
    let mut heights: Vec<u16> = (0..m.len())
        .filter_map(|i| {
            let channel = match residues.get(i)?.to_ascii_uppercase() {
                b'A' => Channel::A,
                b'C' => Channel::C,
                b'G' => Channel::G,
                b'T' => Channel::T,
                _ => return None,
            };
            tr.value(channel, tr.peak(i)?)
        })
        .filter(|&v| v > 0)
        .collect();
    if heights.is_empty() {
        return Some(f32::from(tr.max_value()));
    }
    heights.sort_unstable();
    let k = (heights.len() * 9 / 10).min(heights.len() - 1);
    heights.get(k).map(|&v| f32::from(v) * 1.15)
}

/// The molecule behind a source (for its trace), if it can be had.
fn molecule<'a>(
    view: &'a AlignWithView,
    session: &'a Session,
    id: u64,
    source: &AlignSource,
) -> Option<&'a Molecule> {
    match source {
        AlignSource::File(p) => view.molecules.get(p)?.as_ref(),
        AlignSource::Document(d) => session.document_by_id(*d).map(|d| &d.molecule),
        AlignSource::Matched(i) => view
            .matched
            .as_ref()
            .filter(|(d, _)| *d == id)?
            .1
            .get(*i)?
            .molecule
            .as_ref()
            .ok(),
    }
}

/// What a list row is called.
fn entry_name(view: &AlignWithView, session: &Session, id: u64, s: &AlignSource) -> String {
    match s {
        AlignSource::File(p) => crate::sidebar::base_name(p).to_string(),
        AlignSource::Document(d) => session.document_by_id(*d).map_or_else(
            || format!("tab {d} (closed)"),
            |d| d.molecule.name().to_string(),
        ),
        AlignSource::Matched(i) => view
            .matched
            .as_ref()
            .filter(|(d, _)| *d == id)
            .and_then(|(_, r)| r.get(*i))
            .map_or_else(|| format!("stored read {}", i + 1), |r| r.name.clone()),
    }
}

/// Whether a source can be aligned now: its file is there, its tab is open.
fn available(session: &Session, id: u64, s: &AlignSource) -> bool {
    match s {
        AlignSource::File(p) => session.platform().exists(p),
        AlignSource::Document(d) => *d != id && session.document_by_id(*d).is_some(),
        AlignSource::Matched(_) => true,
    }
}

/// The `check_reads` params for the ticked sequences of a list; `None` when none is ticked.
pub fn params(session: &Session, id: u64, entries: &[AlignEntry]) -> Option<Value> {
    let on = || {
        entries
            .iter()
            .filter(|e| e.on && available(session, id, &e.source))
    };
    let files: Vec<&str> = on()
        .filter_map(|e| match &e.source {
            AlignSource::File(p) => Some(p.as_str()),
            _ => None,
        })
        .collect();
    let documents: Vec<u64> = on()
        .filter_map(|e| match e.source {
            AlignSource::Document(d) => Some(d),
            _ => None,
        })
        .collect();
    let matched = on().any(|e| matches!(e.source, AlignSource::Matched(_)));
    if files.is_empty() && documents.is_empty() && !matched {
        return None;
    }
    Some(json!({
        "document": id,
        "files": files,
        "documents": documents,
        "matched": matched,
        "rows": true,
    }))
}

/// Add the reads stored with document `id` to its list once, ticked: they were matched to it.
pub fn list_stored_reads(
    view: &mut AlignWithView,
    session: &Session,
    id: u64,
    entries: &mut Vec<AlignEntry>,
) {
    let n = view.matched_reads(session, id).len();
    for i in 0..n {
        let s = AlignSource::Matched(i);
        if entries.len() < crate::state::MAX_ALIGN_WITH && !entries.iter().any(|e| e.source == s) {
            entries.push(AlignEntry {
                source: s,
                on: true,
                collapsed: false,
            });
        }
    }
}

/// The panel: the list, ticks, order, and ways to add more. `dragging` is true while files
/// are dragged over the window, so the panel can say it takes them.
pub fn panel(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    id: u64,
    entries: &mut Vec<AlignEntry>,
    view: &mut AlignWithView,
    dragging: bool,
) -> Option<AlignWithAction> {
    list_stored_reads(view, session, id, entries);
    let mut action = None;
    let stroke = if dragging {
        Stroke::new(2.0, t.accent)
    } else {
        Stroke::new(1.0, t.border)
    };
    let height = ui.available_height();
    let frame = egui::Frame::new()
        .fill(t.panel)
        .stroke(stroke)
        .corner_radius(4.0)
        .inner_margin(6.0);
    // A column of its own: inside a row, its lines would be laid out side by side.
    let response = ui.allocate_ui_with_layout(
        vec2(PANEL_W, height),
        egui::Layout::top_down(egui::Align::Min),
        |ui| frame.show(ui, |ui| {
        ui.set_width(PANEL_W - 14.0);
        ui.set_min_height((height - 14.0).max(80.0));
        ui.horizontal(|ui| {
            ui.label(RichText::new("Align with").strong());
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .small_button("\u{d7}")
                    .on_hover_text("Close the alignment")
                    .clicked()
                {
                    action = Some(AlignWithAction::Close);
                }
            });
        });
        let list_h = (ui.available_height() - 96.0).max(60.0);
        egui::ScrollArea::vertical()
            .id_salt("align-with-list")
            .max_height(list_h)
            .min_scrolled_height(list_h)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                if entries.is_empty() {
                    ui.label(
                        RichText::new(if dragging {
                            "Drop the files here to align them."
                        } else {
                            "Drop reads (.ab1, .scf) or other sequences here, or add them below."
                        })
                        .color(t.text_muted),
                    );
                }
                for (i, e) in entries.iter_mut().enumerate() {
                    let name = entry_name(view, session, id, &e.source);
                    let there = available(session, id, &e.source);
                    ui.horizontal(|ui| {
                        ui.spacing_mut().item_spacing.x = 4.0;
                        ui.checkbox(&mut e.on, "")
                            .on_hover_text("Show it aligned under the sequence");
                        let text = RichText::new(short(&name, 24));
                        let text = if there { text } else { text.color(t.error) };
                        let tip = match &e.source {
                            AlignSource::File(p) if !there => format!("{p}: not found"),
                            AlignSource::File(p) => p.clone(),
                            AlignSource::Document(_) if !there => {
                                "The tab was closed".to_string()
                            }
                            AlignSource::Document(_) => format!("{name}, open in a tab"),
                            AlignSource::Matched(_) => {
                                format!("{name}, stored with this file")
                            }
                        };
                        if ui
                            .selectable_label(view.picked == Some(i), text)
                            .on_hover_text(tip)
                            .clicked()
                        {
                            view.picked = if view.picked == Some(i) { None } else { Some(i) };
                        }
                    });
                }
            });
        if view.picked.is_some_and(|i| i >= entries.len()) {
            view.picked = None;
        }
        ui.horizontal(|ui| {
            ui.label("Move:");
            let picked = view.picked;
            let up = picked.filter(|&i| i > 0);
            let down = picked.filter(|&i| i + 1 < entries.len());
            if ui
                .add_enabled(up.is_some(), egui::Button::new("Up"))
                .on_hover_text("Up the list")
                .clicked()
                && let Some(i) = up
            {
                entries.swap(i, i - 1);
                view.picked = Some(i - 1);
            }
            if ui
                .add_enabled(down.is_some(), egui::Button::new("Down"))
                .on_hover_text("Down the list")
                .clicked()
                && let Some(i) = down
            {
                entries.swap(i, i + 1);
                view.picked = Some(i + 1);
            }
            // A stored read comes back by itself, so it can only be unticked.
            let removable = picked.filter(|&i| {
                entries
                    .get(i)
                    .is_some_and(|e| !matches!(e.source, AlignSource::Matched(_)))
            });
            if ui
                .add_enabled(removable.is_some(), egui::Button::new("Remove"))
                .on_hover_text("Take it off the list (reads stored with the file can only be unticked)")
                .clicked()
                && let Some(i) = removable
                && i < entries.len()
            {
                entries.remove(i);
                view.picked = None;
            }
        });
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Add files\u{2026}")
                .on_hover_text("Reads (.ab1, .scf) or sequences (FASTA, GenBank, .dna)")
                .clicked()
            {
                action = Some(AlignWithAction::AddFiles);
            }
            let tabs: Vec<(u64, String)> = session
                .documents()
                .into_iter()
                .filter(|d| d.id != id)
                .filter(|d| {
                    session
                        .document_by_id(d.id)
                        .is_some_and(|doc| doc.molecule.seq().alphabet().is_nucleotide())
                })
                .filter(|d| {
                    !entries
                        .iter()
                        .any(|e| e.source == AlignSource::Document(d.id))
                })
                .map(|d| (d.id, d.name))
                .collect();
            ui.add_enabled_ui(!tabs.is_empty(), |ui| {
                ui.menu_button("From a tab", |ui| {
                    for (d, name) in &tabs {
                        if ui.button(name).clicked()
                            && entries.len() < crate::state::MAX_ALIGN_WITH
                        {
                            entries.push(AlignEntry {
                                source: AlignSource::Document(*d),
                                on: true,
                                collapsed: false,
                            });
                        }
                    }
                })
                .response
                .on_disabled_hover_text("No other sequence is open in a tab")
                .on_hover_text("A read or sequence already open");
            });
        });
    }).response);
    view.panel_rect = Some(response.inner.rect);
    action
}

/// Where a read's column for construct position `p` is in its row, if it reaches there.
fn offset_in_row(read: &Value, len: usize, p: usize) -> Option<usize> {
    let start = to_usize(read.get("start"))?;
    let row_len = read.get("row")?.as_str()?.len();
    let offset = p.checked_add(len)?.checked_sub(start)? % len.max(1);
    (offset < row_len).then_some(offset)
}

/// One read's base at construct position `p`, with its quality.
fn base_at(read: &Value, len: usize, p: usize) -> Option<(u8, Option<u64>)> {
    let offset = offset_in_row(read, len, p)?;
    let b = read.get("row")?.as_str()?.as_bytes().get(offset).copied()?;
    let q = read
        .get("row_quality")
        .and_then(Value::as_array)
        .and_then(|q| q.get(offset))
        .and_then(Value::as_u64);
    Some((b, q))
}

/// The read base behind construct position `p` (untrimmed read, its own orientation).
fn read_base_at(read: &Value, len: usize, p: usize) -> Option<usize> {
    let offset = offset_in_row(read, len, p)?;
    to_usize(read.get("row_positions")?.as_array()?.get(offset))
}

/// Whether the read has extra bases just before construct position `p`.
fn inserts_before(read: &Value, len: usize, p: usize) -> bool {
    let Some(offset) = offset_in_row(read, len, p) else {
        return false;
    };
    read.get("row_insertions")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .any(|x| x.as_u64() == u64::try_from(offset).ok())
}

/// The read's substitutions, and its insertions and deletions.
fn tally(read: &Value) -> (usize, usize) {
    let diffs = read.get("differences").and_then(Value::as_array);
    let kind = |k: &str| {
        diffs
            .into_iter()
            .flatten()
            .filter(|d| d.get("kind").and_then(Value::as_str) == Some(k))
            .count()
    };
    (kind("substitution"), kind("insertion") + kind("deletion"))
}

/// The stretch of the read that aligns, 1-based and inclusive, in its own coordinates.
fn aligned_part(read: &Value) -> Option<(usize, usize)> {
    let at = read.get("row_positions")?.as_array()?;
    let known = at.iter().filter_map(|v| to_usize(Some(v)));
    let lo = known.clone().min()?;
    let hi = known.max()?;
    Some((lo + 1, hi + 1))
}

/// Column geometry: which construct columns are in view and where each is drawn.
#[derive(Clone, Copy)]
struct Columns {
    /// First column in view (fractional).
    left: f64,
    /// Left edge of the sequence area.
    x0: f32,
    cw: f32,
    from: usize,
    to: usize,
}

impl Columns {
    /// Left edge of column `p`.
    fn x(&self, p: usize) -> f32 {
        self.x0 + ((p as f64 - self.left) * f64::from(self.cw)) as f32
    }

    /// The column under screen x, if it is one of the construct's.
    fn at(&self, x: f32, len: usize) -> Option<usize> {
        if x < self.x0 {
            return None;
        }
        let p = self.left + f64::from((x - self.x0) / self.cw);
        if !p.is_finite() || p < 0.0 {
            return None;
        }
        let p = p.floor() as usize;
        (p < len).then_some(p)
    }
}

/// A read laid out for the pane: its list row, and what the check said of it.
struct Row<'a> {
    entry: usize,
    read: &'a Value,
}

/// The pane under the sequence: the construct, then each ticked sequence with its trace, on
/// the columns of the one-line sequence view above (`line_geometry`). It fills the room it is
/// given and never asks for more: the traces shrink to fit, and only if even short ones don't
/// fit do the reads scroll up and down.
#[allow(clippy::too_many_arguments)]
pub fn pane(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    id: u64,
    entries: &mut [AlignEntry],
    view: &mut AlignWithView,
    line_geometry: Option<LineGeometry>,
) -> Option<AlignWithAction> {
    let doc = session.document_by_id(id)?;
    view.pane_rect = Some(ui.max_rect());
    let mut action = None;
    // Align again when the ticked sequences or the construct's bases change.
    let request = params(session, id, entries);
    let key = format!(
        "{id}|{}|{}",
        view.residues_hash(session, id),
        request.as_ref().map(Value::to_string).unwrap_or_default()
    );
    if view.key.as_deref() != Some(key.as_str()) {
        view.key = Some(key);
        view.result = None;
        if let Some(p) = request.clone() {
            action = Some(AlignWithAction::Check(p));
        }
    }
    if request.is_none() {
        ui.label(
            RichText::new(
                "Tick sequences in the Align with list to see them aligned here, under the construct.",
            )
            .color(t.text_muted),
        );
        return action;
    }
    let r = match &view.result {
        None => {
            ui.label(RichText::new("Aligning\u{2026}").color(t.text_muted));
            return action;
        }
        Some(Err(e)) => {
            ui.label(RichText::new(e).color(t.error));
            return action;
        }
        Some(Ok(r)) => r.clone(),
    };
    let reads: Vec<Value> = r
        .get("reads")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    // Each ticked sequence, in the list's order, with its answer.
    let rows: Vec<Row> = entries
        .iter()
        .enumerate()
        .filter(|(_, e)| e.on)
        .filter_map(|(i, e)| {
            let read = reads
                .iter()
                .find(|x| source_of(x).as_ref() == Some(&e.source))?;
            Some(Row { entry: i, read })
        })
        .collect();
    let sources: Vec<AlignSource> = rows
        .iter()
        .filter_map(|row| entries.get(row.entry).map(|e| e.source.clone()))
        .collect();
    view.load(session, id, &sources);
    let residues = doc.molecule.seq().as_bytes();
    let len = residues.len();
    if len == 0 {
        return action;
    }
    // The reads sit on the sequence view's own columns, so each base is right under the
    // construct's above, and the two move as one.
    let Some(geometry) = line_geometry.filter(|g| {
        g.char_w.is_finite() && g.char_w >= 1.0 && g.origin_x.is_finite() && g.clip_left.is_finite()
    }) else {
        ui.label(
            RichText::new("Show the sequence on one line to see the reads under it.")
                .color(t.text_muted),
        );
        return action;
    };
    let font = FontId::monospace(12.0);
    let cw = geometry.char_w;
    let line = ui.fonts_mut(|f| f.row_height(&font)) + 4.0;
    let full = ui.available_rect_before_wrap();
    // The names column runs up to where the sequence view starts drawing residues.
    let x0 = geometry
        .clip_left
        .max(full.left() + MIN_NAME_W)
        .min((full.right() - cw).max(full.left()));
    let name_w = x0 - full.left();
    let x1 = geometry.clip_right.min(full.right());
    // Negative while the start of the line, with the sequence view's gutter, is in view.
    let left = f64::from((x0 - geometry.origin_x) / cw);
    let left = if left.is_finite() { left } else { 0.0 };
    let shown = (((x1 - x0) / cw).ceil().max(1.0)) as usize;
    view.left = left;
    view.shown = shown;
    // Sideways scrolling, a trackpad swipe or the wheel with Shift held, moves the sequence
    // view, and these columns with it.
    if ui.rect_contains_pointer(full) {
        let d = ui.input(|i| i.smooth_scroll_delta);
        let sideways = if d.x != 0.0 {
            d.x
        } else if ui.input(|i| i.modifiers.shift) {
            d.y
        } else {
            0.0
        };
        if sideways.is_finite() && sideways != 0.0 {
            view.scroll_by -= sideways;
        }
    }
    let from = (left.floor().max(0.0) as usize).min(len);
    let to = from.saturating_add(shown + 2).min(len);
    let cols = Columns {
        left,
        x0,
        cw,
        from,
        to,
    };
    let selection = doc.selection;

    // The construct's positions and bases stay at the top while the reads scroll under them.
    let header_h = line * 2.0;
    let (header, header_resp) =
        ui.allocate_exact_size(vec2(full.width(), header_h), Sense::click());
    {
        let painter = ui.painter_at(header);
        painter.rect_filled(header, 0.0, t.canvas);
        shade_selection(&painter, t, &cols, selection, len, header);
        let mut y = header.top() + 2.0;
        for p in from..to {
            if (p + 1) % 10 == 0 {
                let x = cols.x(p) + cw / 2.0;
                painter.line_segment(
                    [pos2(x, y + line - 6.0), pos2(x, y + line - 2.0)],
                    Stroke::new(1.0, t.text_muted),
                );
                painter.text(
                    pos2(x, y),
                    Align2::CENTER_TOP,
                    grouped((p + 1) as u64),
                    FontId::proportional(10.0),
                    t.text_muted,
                );
            }
        }
        y += line;
        painter.text(
            pos2(header.left() + 4.0, y + 2.0),
            Align2::LEFT_TOP,
            short(doc.molecule.name(), 22),
            FontId::proportional(12.5),
            t.text,
        );
        for p in from..to {
            let c = residues
                .get(p)
                .map_or('?', |b| b.to_ascii_uppercase() as char);
            painter.text(
                pos2(cols.x(p) + cw / 2.0, y + 2.0),
                Align2::CENTER_TOP,
                c,
                font.clone(),
                t.text,
            );
        }
        // Keep the names column clear of residues scrolled under it.
        painter.rect_filled(
            Rect::from_min_max(
                pos2(header.left() + name_w - 4.0, header.top()),
                pos2(header.left() + name_w, header.bottom()),
            ),
            0.0,
            t.canvas,
        );
        if header_resp.clicked()
            && let Some(pos) = header_resp.interact_pointer_pos()
            && let Some(p) = cols.at(pos.x, len)
        {
            action = Some(AlignWithAction::Select(p));
        }
    }
    ui.add_space(2.0);

    // The reads, in the room left under the construct. The traces share it, as tall as
    // TRACE_H and as short as MIN_TRACE_H; only when even short ones don't fit do the reads
    // scroll up and down.
    let reads_h = (full.bottom() - ui.cursor().top()).max(0.0);
    let gap = ui.spacing().item_spacing.y;
    let (mut fixed, mut open_traces) = (0.0_f32, 0.0_f32);
    for row in &rows {
        let open = entries.get(row.entry).is_some_and(|e| !e.collapsed)
            && row.read.get("error").is_none_or(Value::is_null);
        fixed += line + 4.0 + gap;
        if open {
            open_traces += 1.0;
        }
    }
    let trace_h = if open_traces > 0.0 {
        let each = (reads_h - fixed) / open_traces;
        if each.is_finite() {
            each.clamp(MIN_TRACE_H, TRACE_H)
        } else {
            TRACE_H
        }
    } else {
        TRACE_H
    };
    let mut toggle = None;
    egui::ScrollArea::vertical()
        .id_salt("align-with-rows")
        .max_height(reads_h)
        .min_scrolled_height(reads_h)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            if rows.is_empty() {
                ui.label(
                    RichText::new("None of the ticked sequences could be aligned.")
                        .color(t.text_muted),
                );
            }
            for row in &rows {
                let Some(entry) = entries.get(row.entry) else {
                    continue;
                };
                let read = row.read;
                let placed = read.get("error").is_none_or(Value::is_null);
                let open = placed && !entry.collapsed;
                let h = line + if open { trace_h } else { 0.0 } + 4.0;
                let (rect, resp) =
                    ui.allocate_exact_size(vec2(full.width(), h), Sense::click());
                if !ui.is_rect_visible(rect) {
                    continue;
                }
                let painter = ui.painter_at(rect);
                painter.rect_filled(rect, 0.0, t.canvas);
                let name_rect = Rect::from_min_size(rect.min, vec2(name_w, line));
                let name = read.get("name").and_then(Value::as_str).unwrap_or("");
                // The fold triangle, the name and the strand.
                if placed {
                    let c = pos2(name_rect.left() + 8.0, name_rect.center().y);
                    let tri = if open {
                        vec![pos2(c.x - 4.0, c.y - 2.5), pos2(c.x + 4.0, c.y - 2.5), pos2(c.x, c.y + 3.5)]
                    } else {
                        vec![pos2(c.x - 2.5, c.y - 4.0), pos2(c.x + 3.5, c.y), pos2(c.x - 2.5, c.y + 4.0)]
                    };
                    painter.add(Shape::convex_polygon(tri, t.text_muted, Stroke::NONE));
                }
                let forward = read.get("strand").and_then(Value::as_str) == Some("forward");
                painter.text(
                    pos2(name_rect.left() + 18.0, name_rect.center().y),
                    Align2::LEFT_CENTER,
                    short(name, 20),
                    FontId::proportional(12.5),
                    t.text,
                );
                if placed {
                    // Which strand it reads: an arrow drawn, as the fonts may lack one.
                    let (y, x) = (name_rect.center().y, name_rect.right() - 16.0);
                    let d = if forward { 1.0 } else { -1.0 };
                    let stroke = Stroke::new(2.0, t.accent);
                    painter.line_segment([pos2(x - 7.0 * d, y), pos2(x + 6.0 * d, y)], stroke);
                    painter.add(Shape::convex_polygon(
                        vec![pos2(x + 2.0 * d, y - 4.0), pos2(x + 8.0 * d, y), pos2(x + 2.0 * d, y + 4.0)],
                        t.accent,
                        Stroke::NONE,
                    ));
                }
                let seq_rect = Rect::from_min_max(
                    pos2(cols.x0, rect.top()),
                    pos2(rect.right(), rect.bottom()),
                );
                let seq_painter = ui.painter_at(seq_rect);
                if !placed {
                    seq_painter.text(
                        pos2(cols.x0 + 4.0, name_rect.center().y),
                        Align2::LEFT_CENTER,
                        format!(
                            "not aligned: {}",
                            read.get("error").and_then(Value::as_str).unwrap_or("")
                        ),
                        FontId::proportional(12.0),
                        t.error,
                    );
                } else {
                    shade_selection(&seq_painter, t, &cols, selection, len, seq_rect);
                    bases(&seq_painter, t, &cols, residues, read, rect.top(), line, &font);
                }
                if open {
                    let strip = Rect::from_min_max(
                        pos2(rect.left(), rect.top() + line),
                        pos2(rect.right(), rect.bottom() - 4.0),
                    );
                    info(&painter, t, read, strip);
                    let trace_rect =
                        Rect::from_min_max(pos2(cols.x0, strip.top()), strip.max);
                    let source = &entry.source;
                    match molecule(view, session, id, source) {
                        Some(m) if m.trace().is_some() => {
                            let peak = view.peaks.get(&key_of(source)).copied().unwrap_or(1.0);
                            trace(&seq_painter, t, &cols, m, read, len, trace_rect, peak, !forward);
                        }
                        _ => {
                            seq_painter.text(
                                trace_rect.left_center() + vec2(4.0, 0.0),
                                Align2::LEFT_CENTER,
                                "no trace: a sequence only",
                                FontId::proportional(11.0),
                                t.text_muted,
                            );
                        }
                    }
                }
                painter.line_segment(
                    [pos2(rect.left(), rect.bottom() - 0.5), pos2(rect.right(), rect.bottom() - 0.5)],
                    Stroke::new(1.0, t.border),
                );
                if resp.clicked()
                    && let Some(pos) = resp.interact_pointer_pos()
                {
                    if pos.x < cols.x0 {
                        if placed && pos.y < rect.top() + line {
                            toggle = Some(row.entry);
                        }
                    } else if let Some(p) = cols.at(pos.x, len) {
                        action = Some(AlignWithAction::Select(p));
                    }
                }
                if placed {
                    resp.on_hover_text_at_pointer(
                        "Click a column to select that base above; click the name to fold the trace",
                    );
                }
            }
        });
    if let Some(i) = toggle
        && let Some(e) = entries.get_mut(i)
    {
        e.collapsed = !e.collapsed;
    }
    action
}

/// Tint the columns of the document's selection.
fn shade_selection(
    painter: &egui::Painter,
    t: &Tokens,
    cols: &Columns,
    selection: Option<snappiesnap_engine::doc::Range>,
    len: usize,
    rect: Rect,
) {
    let Some(sel) = selection else { return };
    let pieces = if sel.start <= sel.end {
        vec![(sel.start, sel.end.max(sel.start + 1))]
    } else {
        vec![(sel.start, len), (0, sel.end)]
    };
    for (a, b) in pieces {
        let (a, b) = (a.max(cols.from), b.min(cols.to));
        if a >= b {
            continue;
        }
        painter.rect_filled(
            Rect::from_min_max(pos2(cols.x(a), rect.top()), pos2(cols.x(b), rect.bottom())),
            0.0,
            t.selection,
        );
    }
}

/// A read's bases in their columns: a base that differs from the construct (or a gap) is
/// boxed, a poor call is faint, and a mark stands where the read has extra bases.
#[allow(clippy::too_many_arguments)]
fn bases(
    painter: &egui::Painter,
    t: &Tokens,
    cols: &Columns,
    residues: &[u8],
    read: &Value,
    top: f32,
    line: f32,
    font: &FontId,
) {
    let len = residues.len();
    for p in cols.from..cols.to {
        let Some((base, quality)) = base_at(read, len, p) else {
            continue;
        };
        let same = residues
            .get(p)
            .is_some_and(|r| r.eq_ignore_ascii_case(&base));
        let x = cols.x(p);
        if !same {
            let cell = Rect::from_min_size(pos2(x, top + 1.0), vec2(cols.cw, line - 2.0));
            painter.rect_filled(cell, 1.0, t.error.gamma_multiply(0.25));
            painter.rect_stroke(
                cell,
                1.0,
                Stroke::new(1.0, t.error),
                egui::StrokeKind::Inside,
            );
        }
        let color = if same { t.text } else { t.error };
        let color = match quality {
            Some(q) if q < FAINT_BELOW => color.gamma_multiply(0.45),
            _ => color,
        };
        painter.text(
            pos2(x + cols.cw / 2.0, top + 2.0),
            Align2::CENTER_TOP,
            base as char,
            font.clone(),
            color,
        );
        if inserts_before(read, len, p) {
            painter.line_segment(
                [pos2(x, top + 1.0), pos2(x, top + line - 1.0)],
                Stroke::new(2.0, t.accent),
            );
        }
    }
}

/// What the read is, beside its trace: its length, the part that aligns, and its differences.
fn info(painter: &egui::Painter, t: &Tokens, read: &Value, strip: Rect) {
    let (subs, gaps) = tally(read);
    let plural =
        |n: usize, one: &str, many: &str| format!("{n} {}", if n == 1 { one } else { many });
    let mut lines = vec![format!("{} bases", grouped(count(read, "read_length")))];
    if let Some((a, b)) = aligned_part(read) {
        lines.push(format!(
            "{}..{} aligned",
            grouped(a as u64),
            grouped(b as u64)
        ));
    }
    lines.push(plural(subs, "mismatch", "mismatches"));
    lines.push(plural(gaps, "gap or insertion", "gaps or insertions"));
    // A short trace leaves room for fewer lines: run the facts together rather than cut them.
    let room = ((strip.height() - 4.0) / 15.0).floor();
    let room = if room.is_finite() {
        room.max(1.0) as usize
    } else {
        1
    };
    if room < lines.len() {
        let per = lines.len().div_ceil(room);
        lines = lines.chunks(per).map(|c| c.join(" \u{b7} ")).collect();
    }
    for (k, text) in lines.iter().enumerate() {
        painter.text(
            pos2(strip.left() + 18.0, strip.top() + 4.0 + 15.0 * k as f32),
            Align2::LEFT_TOP,
            text,
            FontId::proportional(11.0),
            t.text_muted,
        );
    }
}

/// The read's chromatogram drawn under its bases: the stretch of signal between the peaks of
/// two neighbouring columns is stretched to fit between those columns, so every peak sits
/// under the base it was called as. A read of the other strand is drawn as the construct's top
/// strand reads it (mirrored by the columns, with complementary channel colours).
#[allow(clippy::too_many_arguments)]
fn trace(
    painter: &egui::Painter,
    t: &Tokens,
    cols: &Columns,
    m: &Molecule,
    read: &Value,
    len: usize,
    rect: Rect,
    peak: f32,
    reverse: bool,
) {
    let Some(tr) = m.trace() else { return };
    let xs = crate::trace_view::positions(tr, m.len());
    // Column centres with the sample of the read base behind each, a column either side of the
    // view so the curve runs off its edges.
    let anchors: Vec<(f32, f32)> = (cols.from.saturating_sub(1)
        ..cols.to.saturating_add(1).min(len))
        .filter_map(|p| {
            let i = read_base_at(read, len, p)?;
            Some((cols.x(p) + cols.cw / 2.0, *xs.get(i)?))
        })
        .collect();
    let bottom = rect.bottom() - 2.0;
    let height = (rect.height() - 6.0).max(4.0);
    let peak = peak.max(1.0);
    let complement = |c: Channel| match c {
        Channel::A => Channel::T,
        Channel::C => Channel::G,
        Channel::G => Channel::C,
        Channel::T => Channel::A,
    };
    painter.line_segment(
        [pos2(rect.left(), bottom), pos2(rect.right(), bottom)],
        Stroke::new(1.0, t.border),
    );
    for c in Channel::ALL {
        let shown = if reverse { complement(c) } else { c };
        let k = Channel::ALL.iter().position(|x| *x == shown).unwrap_or(0);
        let color: Color32 = t.trace_channels.get(k).copied().unwrap_or(t.text);
        let mut points: Vec<Pos2> = Vec::new();
        for w in anchors.windows(2) {
            let [(x0, s0), (x1, s1)] = w else { continue };
            let span = s1 - s0;
            if span.abs() < 0.5 || !span.is_finite() {
                continue;
            }
            // Every sample between the two peaks, in the order the columns run.
            let (a, b) = (s0.min(*s1).ceil() as i64, s0.max(*s1).floor() as i64);
            let steps: Box<dyn Iterator<Item = i64>> = if span > 0.0 {
                Box::new(a..=b)
            } else {
                Box::new((a..=b).rev())
            };
            for k in steps {
                let Ok(i) = usize::try_from(k) else { continue };
                let Some(v) = tr.value(c, i) else { continue };
                let f = (k as f32 - s0) / span;
                let x = x0 + (x1 - x0) * f;
                let y = bottom - (f32::from(v) / peak).min(1.0) * height;
                points.push(pos2(x, y));
            }
        }
        if points.len() > 1 {
            painter.add(Shape::line(points, Stroke::new(1.0, color)));
        }
    }
}

/// The ticked reads that placed, for the navigator under the sequence view: where each lies
/// on the construct (`len` residues) and where it differs.
pub fn read_spans(view: &AlignWithView, entries: &[AlignEntry], len: usize) -> Vec<ReadSpan> {
    let Some(Ok(r)) = &view.result else {
        return Vec::new();
    };
    let reads = r.get("reads").and_then(Value::as_array);
    let mut out = Vec::new();
    for entry in entries.iter().filter(|e| e.on) {
        let Some(read) = reads
            .into_iter()
            .flatten()
            .find(|x| source_of(x).as_ref() == Some(&entry.source))
        else {
            continue;
        };
        if !read.get("error").is_none_or(Value::is_null) {
            continue;
        }
        let (Some(a), Some(b), Some(n)) = (
            to_usize(read.get("start")),
            to_usize(read.get("end")),
            to_usize(read.get("length")),
        ) else {
            continue;
        };
        let (a, b) = (a.min(len), b.min(len));
        let pieces = if n >= len {
            vec![(0, len)]
        } else if a <= b {
            vec![(a, b)]
        } else {
            vec![(a, len), (0, b)]
        };
        let differences = read
            .get("differences")
            .and_then(Value::as_array)
            .into_iter()
            .flatten()
            .filter_map(|d| to_usize(d.get("position").or_else(|| d.get("start"))))
            .filter(|&p| p < len)
            .collect();
        out.push(ReadSpan {
            name: read
                .get("name")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string(),
            pieces,
            forward: read.get("strand").and_then(Value::as_str) == Some("forward"),
            differences,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Call, Seq, Signal, Topology, Trace};

    fn construct(n: usize) -> String {
        // A sequence without long repeats, so reads place in one way only.
        let mut x: u32 = 12345;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                b"ACGT"[((x >> 16) % 4) as usize] as char
            })
            .collect()
    }

    /// A read with a trace: one clean peak per base, every 10 samples.
    fn traced(name: &str, residues: &str) -> Molecule {
        let n = residues.len();
        let samples = n * 10 + 10;
        let mut signal = Signal {
            a: vec![0; samples],
            c: vec![0; samples],
            g: vec![0; samples],
            t: vec![0; samples],
        };
        for (i, b) in residues.bytes().enumerate() {
            let ch = match b {
                b'A' => &mut signal.a,
                b'C' => &mut signal.c,
                b'G' => &mut signal.g,
                _ => &mut signal.t,
            };
            for (k, v) in [(-2i64, 300u16), (-1, 700), (0, 1000), (1, 700), (2, 300)] {
                if let Some(s) = usize::try_from(i as i64 * 10 + 5 + k)
                    .ok()
                    .and_then(|j| ch.get_mut(j))
                {
                    *s = v;
                }
            }
        }
        let calls = (0..n)
            .map(|i| Call {
                peak: Some((i * 10 + 5) as u32),
                quality: Some(40),
            })
            .collect();
        let (tr, _) = Trace::new(signal, calls, None).unwrap();
        let mut m = Molecule::new(
            name,
            Seq::new(Alphabet::Dna, residues).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        m.set_trace(Some(tr)).unwrap();
        m
    }

    fn revcomp(s: &str) -> String {
        s.bytes()
            .rev()
            .map(|b| match b {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                _ => 'A',
            })
            .collect()
    }

    /// Two reads open in tabs, one of each strand, aligned under a construct and drawn with
    /// their traces; the list ticks, orders and folds them. Set `SNAPPIESNAP_ALIGN_WITH_PNG`
    /// to a path to look at it.
    #[test]
    fn reads_are_laid_out_under_the_construct_with_their_traces() {
        let c = construct(1500);
        let mut session = Session::new();
        session.open_document(Document::new(
            Molecule::new(
                "pTest",
                Seq::new(Alphabet::Dna, c.as_str()).unwrap(),
                Topology::Circular,
            )
            .unwrap(),
        ));
        let construct_id = session.active_document_id().unwrap();
        // Forward read of 300..800 with a substitution at 500.
        let mut fwd = c.as_bytes()[300..800].to_vec();
        fwd[200] = if fwd[200] == b'A' { b'C' } else { b'A' };
        session.open_document(Document::new(traced(
            "fwd",
            std::str::from_utf8(&fwd).unwrap(),
        )));
        let fwd_id = session.active_document_id().unwrap();
        session.open_document(Document::new(traced("rev", &revcomp(&c[450..950]))));
        let rev_id = session.active_document_id().unwrap();
        session.switch_document(0).unwrap();
        assert_eq!(session.active_document_id(), Some(construct_id));

        let mut entries = vec![
            AlignEntry {
                source: AlignSource::Document(fwd_id),
                on: true,
                collapsed: false,
            },
            AlignEntry {
                source: AlignSource::Document(rev_id),
                on: true,
                collapsed: false,
            },
            AlignEntry {
                source: AlignSource::File("/nowhere/read.ab1".into()),
                on: true,
                collapsed: false,
            },
        ];
        // A missing file is left out rather than failing the whole check.
        let p = params(&session, construct_id, &entries).unwrap();
        assert_eq!(p["documents"], json!([fwd_id, rev_id]));
        assert_eq!(p["files"], json!([]));
        assert_eq!(p["matched"], false);

        let mut view = AlignWithView::default();
        let t = Tokens::light();
        let draw = |session: &Session, entries: &mut Vec<AlignEntry>, view: &mut AlignWithView| {
            let mut asked = None;
            let image = crate::offscreen::render(
                1400,
                520,
                1.0,
                |ctx| t.apply(ctx),
                |ui| {
                    egui::Frame::new().fill(t.panel).show(ui, |ui| {
                        ui.horizontal_top(|ui| {
                            let _ = panel(ui, &t, session, construct_id, entries, view, false);
                            ui.vertical(|ui| {
                                // The sequence view above shows residue 450 at the left of its
                                // residues, 200 points in.
                                let r = ui.max_rect();
                                let line = LineGeometry {
                                    origin_x: r.left() + 200.0 - 450.0 * 8.0,
                                    char_w: 8.0,
                                    clip_left: r.left() + 200.0,
                                    clip_right: r.right(),
                                };
                                // Offscreen frames may run twice: keep what any of them asked.
                                if let Some(a) =
                                    pane(ui, &t, session, construct_id, entries, view, Some(line))
                                {
                                    asked = Some(a);
                                }
                            });
                        });
                    });
                },
            )
            .unwrap();
            (image, asked)
        };
        // The first frame asks for the check; the next draws it.
        let (_, asked) = draw(&session, &mut entries, &mut view);
        let Some(AlignWithAction::Check(p)) = asked else {
            panic!("the pane asks for the check: {asked:?}");
        };
        let answer = session.execute("check_reads", p).map_err(|e| e.message);
        let r = answer.clone().unwrap();
        view.set_result(answer);
        let reads = r["reads"].as_array().unwrap();
        let rev = reads
            .iter()
            .find(|x| source_of(x) == Some(AlignSource::Document(rev_id)))
            .unwrap();
        assert_eq!(rev["strand"], "reverse");
        // The reverse read's base behind construct 450 is its last one.
        assert_eq!(read_base_at(rev, c.len(), 450), Some(499));
        assert_eq!(read_base_at(rev, c.len(), 949), Some(0));
        let fwd_read = reads
            .iter()
            .find(|x| source_of(x) == Some(AlignSource::Document(fwd_id)))
            .unwrap();
        assert_eq!(tally(fwd_read), (1, 0));
        assert_eq!(aligned_part(fwd_read), Some((1, 500)));

        // Selecting 500 above brings it into view.
        session
            .execute("select", json!({"range": {"start": 500, "end": 501}}))
            .unwrap();
        let (image, asked) = draw(&session, &mut entries, &mut view);
        assert_eq!(asked, None, "the same sequences are not aligned again");
        assert!(view.left <= 500.0 && 500.0 < view.left + view.shown as f64);
        if let Ok(path) = std::env::var("SNAPPIESNAP_ALIGN_WITH_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        // Each channel colour of the trace shows up in the pane.
        let has = |c: Color32| image.pixels.contains(&c);
        assert!(
            t.trace_channels.iter().filter(|c| has(**c)).count() >= 3,
            "the traces are drawn"
        );

        // Unticking a read asks for a check without it; folding one does not.
        entries[1].on = false;
        let (_, asked) = draw(&session, &mut entries, &mut view);
        assert!(
            matches!(asked, Some(AlignWithAction::Check(p)) if p["documents"] == json!([fwd_id]))
        );
        entries[0].collapsed = true;
        let (_, asked) = draw(&session, &mut entries, &mut view);
        assert_eq!(asked, None);
    }

    /// The pane draws within the room it is given and never asks for more: the traces shrink
    /// to fit, and the reads scroll only when even short traces don't fit. (Sized by its
    /// contents, it once crept up or down a little with every repaint.) Its columns are the
    /// sequence view's, and a swipe over it is handed to the sequence view.
    #[test]
    fn pane_fits_its_room_on_the_sequence_columns() {
        let c = construct(1200);
        let mut session = Session::new();
        session.open_document(Document::new(
            Molecule::new(
                "pTest",
                Seq::new(Alphabet::Dna, c.as_str()).unwrap(),
                Topology::Circular,
            )
            .unwrap(),
        ));
        let construct_id = session.active_document_id().unwrap();
        let mut reads = Vec::new();
        for (name, from) in [("a", 100), ("b", 300), ("c", 500)] {
            session.open_document(Document::new(traced(name, &c[from..from + 400])));
            reads.push(session.active_document_id().unwrap());
        }
        session.switch_document(0).unwrap();
        let t = Tokens::dark();
        for scale in [1.0_f32, 4.0 / 3.0] {
            for room in [60.0_f32, 140.0, 220.0, 480.0] {
                let mut entries: Vec<AlignEntry> = reads
                    .iter()
                    .map(|&id| AlignEntry {
                        source: AlignSource::Document(id),
                        on: true,
                        collapsed: false,
                    })
                    .collect();
                let mut view = AlignWithView::default();
                let ctx = egui::Context::default();
                t.apply(&ctx);
                let area = Rect::from_min_size(pos2(0.0, 300.0), vec2(1300.0, room));
                let line = LineGeometry {
                    origin_x: 280.0 - 250.0 * 8.5,
                    char_w: 8.5,
                    clip_left: 280.0,
                    clip_right: 1300.0,
                };
                for frame in 0..6u32 {
                    let mut input = egui::RawInput {
                        screen_rect: Some(Rect::from_min_size(Pos2::ZERO, vec2(1300.0, 900.0))),
                        time: Some(f64::from(frame) * 0.1),
                        ..Default::default()
                    };
                    input
                        .viewports
                        .entry(egui::ViewportId::ROOT)
                        .or_default()
                        .native_pixels_per_point = Some(scale);
                    if frame == 5 {
                        input
                            .events
                            .push(egui::Event::PointerMoved(pos2(700.0, 300.0 + room / 2.0)));
                        input.events.push(egui::Event::MouseWheel {
                            unit: egui::MouseWheelUnit::Point,
                            delta: vec2(-40.0, 0.0),
                            modifiers: egui::Modifiers::NONE,
                            phase: egui::TouchPhase::Move,
                            source: egui::MouseWheelSource::Trackpad,
                        });
                    }
                    let mut asked = None;
                    let mut used = Rect::NOTHING;
                    let mut out = ctx.run_ui(input, |ui| {
                        ui.scope_builder(egui::UiBuilder::new().max_rect(area), |ui| {
                            asked = pane(
                                ui,
                                &t,
                                &session,
                                construct_id,
                                &mut entries,
                                &mut view,
                                Some(line),
                            );
                            used = ui.min_rect();
                        });
                    });
                    out.textures_delta.clear();
                    if let Some(AlignWithAction::Check(p)) = asked {
                        view.set_result(session.execute("check_reads", p).map_err(|e| e.message));
                    }
                    if frame >= 2 {
                        assert!(
                            used.bottom() <= area.bottom() + 0.01,
                            "scale {scale}, room {room}: the pane takes {used:?} of {area:?}"
                        );
                    }
                }
                // Residue 250 is at the left of the residues, as in the sequence view.
                assert!((view.left - 250.0).abs() < 1e-3, "{}", view.left);
                // The swipe is passed on (smoothed by egui, so some of it at least).
                assert!(view.scroll_by > 0.0, "{}", view.scroll_by);
            }
        }
    }

    /// Hostile geometry never panics: an empty construct, columns off either end, NaN.
    #[test]
    fn columns_hold_up_at_the_edges() {
        let cols = Columns {
            left: 0.0,
            x0: 100.0,
            cw: 7.0,
            from: 0,
            to: 10,
        };
        assert_eq!(cols.at(50.0, 10), None);
        assert_eq!(cols.at(100.0, 10), Some(0));
        assert_eq!(cols.at(100.0 + 7.0 * 10.5, 10), None);
        assert_eq!(cols.at(f32::NAN, 10), None);
        let read = json!({"start": 8, "row": "ACGT", "row_positions": [0, null, 2, 3]});
        // A row that wraps the origin of a 10 bp circle.
        assert_eq!(base_at(&read, 10, 9).map(|b| b.0), Some(b'C'));
        assert_eq!(read_base_at(&read, 10, 9), None);
        assert_eq!(read_base_at(&read, 10, 1), Some(3));
        assert_eq!(base_at(&read, 10, 2), None);
        assert_eq!(base_at(&read, 0, 0), None);
        assert_eq!(aligned_part(&json!({"row_positions": [null]})), None);
    }
}
