//! Check Reads (an experiment design): a construct checked with sequencing reads. Step 1
//! chooses the construct and the read files; step 2 shows what the reads cover, how each one
//! lies, and every difference pooled over the reads, with the traces at a chosen difference.
//! The check is the `check_reads` command's, so agents see the same results.

use std::collections::HashMap;

use eframe::egui::{self, Align2, FontId, Rect, RichText, Sense, Stroke, Ui, pos2, vec2};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::Molecule;

use crate::state::{Experiment, ExperimentStep};
use crate::theme::Tokens;

/// File extensions offered as reads.
pub const READ_EXTENSIONS: &[&str] = &["ab1", "abi", "scf", "fa", "fasta", "seq", "gb", "gbk"];

/// What the user asked for.
#[derive(Debug, Clone, PartialEq)]
pub enum ReadsAction {
    /// Pick read files to add.
    AddFiles,
    /// Run `check_reads` with these params; give the answer to [`ReadsInputs::set_result`].
    Check(Value),
    /// Run a command (opening a read in a tab of its own).
    Run(&'static str, Value),
}

/// The check's answer, the chosen difference, and the reads loaded for their traces.
#[derive(Default)]
pub struct ReadsInputs {
    key: Option<String>,
    result: Option<Result<Value, String>>,
    selected: Option<usize>,
    molecules: HashMap<String, Option<Molecule>>,
    /// The construct's matched reads, and the document and revision they were read for.
    matched: Option<(u64, u64, Vec<snappiesnap_engine::io::MatchedRead>)>,
    /// Show differences of every quality, not only those called at Q20 or better.
    show_all: bool,
    /// Draw the alignment as one long line scrolled sideways, rather than in blocks.
    single_line: bool,
    /// Move the alignment grid to the chosen difference on the next frame (set when it was
    /// chosen somewhere other than the grid, where it is already in view).
    scroll_to_selected: bool,
}

/// Differences called at this quality or better are always listed.
const GOOD: u64 = 20;

/// Whether a difference is listed: any read shows it at Q20 or better, or its reads have no
/// qualities (plain sequences).
fn visible(f: &Value, show_all: bool) -> bool {
    show_all
        || f["reads"]
            .as_array()
            .into_iter()
            .flatten()
            .any(|x| match x["quality"].as_u64() {
                Some(q) => q >= GOOD,
                None => true,
            })
}

impl ReadsInputs {
    /// The answer of `check_reads`. The first disputed difference (else the first one) is
    /// chosen, so its traces show straight away.
    pub fn set_result(&mut self, answer: Result<Value, String>) {
        let findings = answer
            .as_ref()
            .ok()
            .and_then(|r| r.get("findings"))
            .and_then(Value::as_array);
        let show_all = self.show_all;
        self.selected = findings.and_then(|f| {
            f.iter()
                .position(|x| {
                    visible(x, show_all)
                        && x.get("verdict").and_then(Value::as_str) == Some("disputed")
                })
                .or_else(|| f.iter().position(|x| visible(x, show_all)))
        });
        self.scroll_to_selected = true;
        self.result = Some(answer);
    }

    /// The reads matched to document `id` in its file (read again when it changes).
    fn matched_reads(
        &mut self,
        session: &Session,
        id: u64,
    ) -> &[snappiesnap_engine::io::MatchedRead] {
        let revision = session.revision();
        if self
            .matched
            .as_ref()
            .is_none_or(|(d, r, _)| *d != id || *r != revision)
        {
            let reads =
                snappiesnap_engine::matched::matched_reads(session, Some(id)).unwrap_or_default();
            self.matched = Some((id, revision, reads));
        }
        self.matched.as_ref().map_or(&[], |(_, _, r)| r.as_slice())
    }

    /// A read file as a molecule (with its trace, if it has one); loaded once.
    fn molecule(&mut self, session: &Session, path: &str) -> Option<&Molecule> {
        self.molecules
            .entry(path.to_string())
            .or_insert_with(|| {
                let bytes = session.platform().read(path).ok()?;
                snappiesnap_engine::io::import(path, &bytes)
                    .ok()
                    .map(|i| i.molecule)
            })
            .as_ref()
    }
}

/// A field of a read's `source` in a `check_reads` answer.
fn source<'a>(read: &'a Value, key: &str) -> Option<&'a Value> {
    read.get("source")?.get(key)
}

fn name_of(session: &Session, id: u64) -> Option<String> {
    session
        .document_by_id(id)
        .map(|d| d.molecule.name().to_string())
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

fn steps(ui: &mut Ui, t: &Tokens, current: u8) {
    ui.horizontal(|ui| {
        for (n, label) in [(1u8, "Construct and reads"), (2, "Results")] {
            if n > 1 {
                ui.label(RichText::new("\u{203a}").color(t.text_muted));
            }
            let text = RichText::new(format!("{n}  {label}"));
            ui.label(if n == current {
                text.strong().color(t.text)
            } else {
                text.color(t.text_muted)
            });
        }
    });
    ui.add_space(t.spacing);
}

/// Step 1: the construct and the read files.
pub fn setup(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut ReadsInputs,
    folder: Option<&str>,
) -> Option<ReadsAction> {
    steps(ui, t, 1);
    let mut action = None;
    ui.label(
        RichText::new(
            "Choose the construct the reads should match, and the read files: traces (.ab1, \
             .scf) or sequences (FASTA, GenBank). Each read is trimmed of its poor ends and \
             placed on the construct on either strand.",
        )
        .color(t.text_muted),
    );
    ui.add_space(t.spacing);
    let docs: Vec<(u64, String)> = session
        .documents()
        .into_iter()
        .filter(|d| {
            session
                .document_by_id(d.id)
                .is_some_and(|doc| doc.molecule.seq().alphabet().is_nucleotide())
        })
        .map(|d| (d.id, d.name))
        .collect();
    ui.horizontal(|ui| {
        ui.label("Construct");
        let shown = exp
            .reference
            .and_then(|id| name_of(session, id))
            .unwrap_or_else(|| "Choose\u{2026}".into());
        egui::ComboBox::from_id_salt("reads-construct")
            .selected_text(shown)
            .width(320.0)
            .show_ui(ui, |ui| {
                for (id, name) in &docs {
                    ui.selectable_value(&mut exp.reference, Some(*id), name);
                }
            });
    });
    // Reads already open in tabs: the only way in the browser until a file is loaded, and the
    // quickest way anywhere once a read has been opened to look at its trace.
    let others: Vec<(u64, String)> = docs
        .iter()
        .filter(|(id, _)| Some(*id) != exp.reference)
        .cloned()
        .collect();
    if !others.is_empty() {
        ui.add_space(t.spacing);
        ui.label(RichText::new("Reads open in tabs").strong());
        egui::ScrollArea::vertical()
            .id_salt("reads-tabs")
            .max_height(140.0)
            .show(ui, |ui| {
                for (id, name) in &others {
                    let mut on = exp.read_documents.contains(id);
                    if ui.checkbox(&mut on, name).changed() {
                        if on {
                            exp.read_documents.push(*id);
                        } else {
                            exp.read_documents.retain(|x| x != id);
                        }
                    }
                }
            });
    }
    // Tabs closed since they were ticked, and the construct if it was chosen as a read.
    exp.read_documents
        .retain(|id| others.iter().any(|(i, _)| i == id));
    let matched = exp
        .reference
        .map_or(0, |id| inputs.matched_reads(session, id).len());
    if matched > 0 {
        let name = exp
            .reference
            .and_then(|id| name_of(session, id))
            .unwrap_or_default();
        ui.checkbox(
            &mut exp.matched,
            format!("Include the {matched} reads stored with {name}"),
        )
        .on_hover_text("Reads matched to the construct in its .dna file, kept as trimmed there");
    }
    ui.add_space(t.spacing);
    ui.label(RichText::new(format!("Read files ({})", exp.reads.len())).strong());
    let mut remove = None;
    egui::ScrollArea::vertical()
        .id_salt("reads-files")
        .max_height(260.0)
        .show(ui, |ui| {
            for (i, path) in exp.reads.iter().enumerate() {
                ui.horizontal(|ui| {
                    if ui.small_button("Remove").clicked() {
                        remove = Some(i);
                    }
                    ui.label(crate::sidebar::base_name(path))
                        .on_hover_text(path.as_str());
                    if !session.platform().exists(path) {
                        ui.label(RichText::new("not found").color(t.error));
                    }
                });
            }
        });
    if let Some(i) = remove
        && i < exp.reads.len()
    {
        exp.reads.remove(i);
    }
    ui.horizontal(|ui| {
        if ui.button("Add read files\u{2026}").clicked() {
            action = Some(ReadsAction::AddFiles);
        }
        if let Some(folder) = folder
            && ui
                .button("Add the traces in the folder")
                .on_hover_text(format!("Every .ab1 and .scf file in {folder}"))
                .clicked()
        {
            let mut found: Vec<String> = session
                .platform()
                .list_dir(folder)
                .unwrap_or_default()
                .into_iter()
                .filter(|e| {
                    !e.is_dir
                        && e.name.rsplit_once('.').is_some_and(|(_, x)| {
                            matches!(x.to_ascii_lowercase().as_str(), "ab1" | "abi" | "scf")
                        })
                })
                .map(|e| e.path)
                .collect();
            found.sort();
            for p in found {
                add_read(exp, p);
            }
        }
    });
    ui.add_space(t.spacing * 2.0);
    let ready = exp.reference.is_some()
        && (!exp.reads.is_empty()
            || !exp.read_documents.is_empty()
            || (exp.matched && matched > 0));
    if ui
        .add_enabled(ready, egui::Button::new("Continue: check the reads"))
        .clicked()
    {
        exp.step = ExperimentStep::Results;
    }
    action
}

/// Add a read file once.
pub fn add_read(exp: &mut Experiment, path: String) {
    if !exp.reads.contains(&path) && exp.reads.len() < snappiesnap_engine::bio::reads::MAX_READS {
        exp.reads.push(path);
    }
}

/// Step 2: the results.
pub fn results(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut ReadsInputs,
) -> Option<ReadsAction> {
    steps(ui, t, 2);
    let mut back = false;
    ui.horizontal(|ui| {
        back = ui.button("Back to construct and reads").clicked();
    });
    if back {
        exp.step = ExperimentStep::Constructs;
        return None;
    }
    let Some(reference) = exp
        .reference
        .filter(|id| session.document_by_id(*id).is_some())
    else {
        ui.label(
            RichText::new("The construct isn't open any more: choose it again.").color(t.error),
        );
        return None;
    };
    let key = format!(
        "{}|{reference}|{:?}|{:?}|{}",
        session.revision(),
        exp.reads,
        exp.read_documents,
        exp.matched
    );
    let mut request = None;
    if inputs.key.as_deref() != Some(key.as_str()) {
        inputs.key = Some(key);
        inputs.result = None;
        request = Some(ReadsAction::Check(json!({
            "document": reference,
            "files": exp.reads,
            "documents": exp.read_documents,
            "matched": exp.matched,
            "rows": true,
        })));
    }
    ui.add_space(t.spacing);
    let r = match &inputs.result {
        None => {
            ui.label(RichText::new("Checking\u{2026}").color(t.text_muted));
            return request;
        }
        Some(Err(e)) => {
            ui.label(RichText::new(e).color(t.error));
            return request;
        }
        Some(Ok(r)) => r.clone(),
    };
    let len = r.get("length").and_then(Value::as_u64).unwrap_or(0);
    let reads = r
        .get("reads")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let findings = r
        .get("findings")
        .and_then(Value::as_array)
        .cloned()
        .unwrap_or_default();
    let count = |v: &str| {
        findings
            .iter()
            .filter(|f| f.get("verdict").and_then(Value::as_str) == Some(v))
            .count()
    };
    let placed = r.get("placed").and_then(Value::as_u64).unwrap_or(0);
    let coverage = r.get("coverage").and_then(Value::as_f64).unwrap_or(0.0);
    ui.label(
        RichText::new(format!(
            "{placed} of {} reads placed  \u{b7}  {:.1}% of {} covered ({} bp by two or more)  \u{b7}  differences: {} confirmed, {} disputed, {} in one read",
            reads.len(),
            coverage * 100.0,
            r.get("construct").and_then(Value::as_str).unwrap_or(""),
            grouped(r.get("covered_twice").and_then(Value::as_u64).unwrap_or(0)),
            count("confirmed"),
            count("disputed"),
            count("one_read"),
        ))
        .strong(),
    );
    ui.add_space(t.spacing);
    let residues = session
        .document_by_id(reference)
        .map(|d| d.molecule.seq().as_bytes().to_vec())
        .unwrap_or_default();
    // What was found on the left, the reads themselves on the right: the alignment is wide and
    // tall, and in one column it was squeezed under the summary with the page's width unused.
    let height = ui.available_height();
    let Some(left) = crate::two_columns(ui.available_width()) else {
        egui::ScrollArea::vertical()
            .id_salt("reads-results")
            .show(ui, |ui| {
                if let Some(a) = summary(ui, t, &r, len, &reads, &findings, reference, inputs) {
                    request = Some(a);
                }
                ui.add_space(t.spacing);
                alignment_grid(ui, t, &residues, &reads, &findings, inputs, GRID_HEIGHT);
                if let Some(f) = inputs.selected.and_then(|i| findings.get(i)) {
                    ui.add_space(t.spacing);
                    traces_at(ui, t, session, f, &reads, exp, inputs);
                }
            });
        return request;
    };
    ui.horizontal_top(|ui| {
        ui.allocate_ui_with_layout(
            vec2(left, height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                // Held to its own width, or a long line of the summary would push the
                // alignment off the window instead of wrapping.
                ui.set_max_width(left);
                egui::ScrollArea::vertical()
                    .id_salt("reads-summary")
                    .max_width(left)
                    .show(ui, |ui| {
                        ui.set_max_width(left);
                        if let Some(a) =
                            summary(ui, t, &r, len, &reads, &findings, reference, inputs)
                        {
                            request = Some(a);
                        }
                    });
            },
        );
        ui.separator();
        egui::ScrollArea::vertical()
            .id_salt("reads-alignment-column")
            .show(ui, |ui| {
                // A column of its own: inside a row, the traces would be laid out beside the
                // alignment rather than under it, off the window.
                ui.vertical(|ui| {
                    // With a difference chosen, leave the traces under it room to be seen.
                    let for_grid = if inputs.selected.is_some() {
                        (height * 0.55).clamp(260.0, 700.0)
                    } else {
                        (height - 60.0).clamp(GRID_HEIGHT, 900.0)
                    };
                    alignment_grid(ui, t, &residues, &reads, &findings, inputs, for_grid);
                    if let Some(f) = inputs.selected.and_then(|i| findings.get(i)) {
                        ui.add_space(t.spacing);
                        traces_at(ui, t, session, f, &reads, exp, inputs);
                    }
                });
            });
    });
    request
}

/// What the check found: how the reads lie, and every difference pooled over them.
#[allow(clippy::too_many_arguments)]
fn summary(
    ui: &mut Ui,
    t: &Tokens,
    r: &Value,
    len: u64,
    reads: &[Value],
    findings: &[Value],
    reference: u64,
    inputs: &mut ReadsInputs,
) -> Option<ReadsAction> {
    coverage_map(ui, t, len, reads, findings, r, inputs);
    ui.add_space(t.spacing);
    let action = read_table(ui, t, reads, reference);
    ui.add_space(t.spacing);
    finding_list(ui, t, findings, reads, inputs);
    action
}

/// The construct as a line, each placed read as an arrow above it, stretches no read covers
/// marked under it, and differences as ticks (click one to see its traces).
fn coverage_map(
    ui: &mut Ui,
    t: &Tokens,
    len: u64,
    reads: &[Value],
    findings: &[Value],
    r: &Value,
    inputs: &mut ReadsInputs,
) {
    let placed: Vec<&Value> = reads.iter().filter(|x| x["error"].is_null()).collect();
    let row = 12.0;
    let height = 40.0 + row * placed.len().min(40) as f32;
    let width = ui.available_width().min(1100.0);
    let (rect, response) = ui.allocate_exact_size(vec2(width, height), Sense::click());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, 3.0, t.canvas);
    if len == 0 {
        return;
    }
    let (left, right) = (rect.left() + 10.0, rect.right() - 10.0);
    let x_of = |p: u64| left + (right - left) * p as f32 / len as f32;
    let axis = rect.bottom() - 16.0;
    painter.line_segment(
        [pos2(left, axis), pos2(right, axis)],
        Stroke::new(2.0, t.backbone),
    );
    for p in [0, len] {
        painter.text(
            pos2(x_of(p), axis + 3.0),
            if p == 0 {
                Align2::LEFT_TOP
            } else {
                Align2::RIGHT_TOP
            },
            grouped(p.max(1)),
            FontId::proportional(10.0),
            t.text_muted,
        );
    }
    for gap in r
        .get("uncovered")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
    {
        let (Some(a), Some(b)) = (gap["start"].as_u64(), gap["end"].as_u64()) else {
            continue;
        };
        let pieces = if a <= b {
            vec![(a, b)]
        } else {
            vec![(a, len), (0, b)]
        };
        for (a, b) in pieces {
            painter.line_segment(
                [pos2(x_of(a), axis), pos2(x_of(b), axis)],
                Stroke::new(4.0, t.error),
            );
        }
    }
    for (k, read) in placed.iter().enumerate().take(40) {
        let (Some(a), Some(b), Some(length)) = (
            read["start"].as_u64(),
            read["end"].as_u64(),
            read["length"].as_u64(),
        ) else {
            continue;
        };
        let y = rect.top() + 8.0 + row * k as f32;
        let forward = read["strand"] == "forward";
        let color = if forward {
            t.feature_color(0)
        } else {
            t.feature_color(1)
        };
        let pieces = if length >= len {
            vec![(0, len)]
        } else if a <= b {
            vec![(a, b)]
        } else {
            vec![(a, len), (0, b)]
        };
        for (pa, pb) in pieces {
            painter.line_segment(
                [pos2(x_of(pa), y), pos2(x_of(pb), y)],
                Stroke::new(5.0, color),
            );
        }
        let tip = if forward { x_of(b) } else { x_of(a) };
        let dir = if forward { 1.0 } else { -1.0 };
        painter.add(egui::Shape::convex_polygon(
            vec![
                pos2(tip, y - 4.5),
                pos2(tip + 5.0 * dir, y),
                pos2(tip, y + 4.5),
            ],
            color,
            Stroke::NONE,
        ));
    }
    let mut ticks = Vec::new();
    for (i, f) in findings.iter().enumerate() {
        let p = f
            .get("position")
            .or_else(|| f.get("start"))
            .and_then(Value::as_u64)
            .unwrap_or(0);
        let color = match f.get("verdict").and_then(Value::as_str) {
            Some("confirmed") => t.accent,
            Some("disputed") => t.error,
            _ => t.text_muted,
        };
        let x = x_of(p);
        let w = if inputs.selected == Some(i) { 3.0 } else { 1.5 };
        painter.line_segment(
            [pos2(x, axis - 8.0), pos2(x, axis + 2.0)],
            Stroke::new(w, color),
        );
        ticks.push((x, i));
    }
    if response.clicked()
        && let Some(pos) = response.interact_pointer_pos()
        && let Some((_, i)) = ticks
            .iter()
            .filter(|(x, _)| (x - pos.x).abs() < 6.0)
            .min_by(|a, b| (a.0 - pos.x).abs().total_cmp(&(b.0 - pos.x).abs()))
    {
        inputs.selected = Some(*i);
        inputs.scroll_to_selected = true;
    }
}

/// Width of the name column of the alignment grid.
const GRID_NAME: f32 = 136.0;
/// A call below this quality is drawn faint: near a read's ends such calls are often noise.
const FAINT_BELOW: u64 = 20;
/// The alignment's height when it is stacked under the summary, and its smallest beside it.
const GRID_HEIGHT: f32 = 340.0;

/// One read's base at construct position `p`, with its quality, from the `row` the check laid
/// out against the construct. `None` when the read doesn't reach there.
fn base_at(read: &Value, len: usize, p: usize) -> Option<(u8, Option<u64>)> {
    let row = read.get("row")?.as_str()?.as_bytes();
    let start = usize::try_from(read.get("start")?.as_u64()?).ok()?;
    let offset = p.checked_add(len)?.checked_sub(start)? % len.max(1);
    let b = row.get(offset).copied()?;
    let q = read
        .get("row_quality")
        .and_then(Value::as_array)
        .and_then(|q| q.get(offset))
        .and_then(Value::as_u64);
    Some((b, q))
}

/// Whether the read has extra bases just before construct position `p`.
fn inserts_before(read: &Value, len: usize, p: usize) -> bool {
    let Some(start) = read.get("start").and_then(Value::as_u64) else {
        return false;
    };
    let Ok(start) = usize::try_from(start) else {
        return false;
    };
    let Some(offset) = p.checked_add(len).and_then(|x| x.checked_sub(start)) else {
        return false;
    };
    let offset = offset % len.max(1);
    read.get("row_insertions")
        .and_then(Value::as_array)
        .into_iter()
        .flatten()
        .any(|x| x.as_u64() == u64::try_from(offset).ok())
}

/// The construct position a finding sits at.
fn finding_at(f: &Value) -> Option<u64> {
    f.get("position").or_else(|| f.get("start"))?.as_u64()
}

/// The reads laid out under the construct, one column per base, in blocks. Letters that differ
/// from the construct stand out; a call below Q20 is drawn faint, because the noise near a
/// read's ends is called at low quality. Clicking a column picks the difference there, so its
/// traces show below. Only the stretch the reads lie on is drawn.
#[allow(clippy::too_many_arguments)]
fn alignment_grid(
    ui: &mut Ui,
    t: &Tokens,
    residues: &[u8],
    reads: &[Value],
    findings: &[Value],
    inputs: &mut ReadsInputs,
    height: f32,
) {
    let len = residues.len();
    let placed: Vec<&Value> = reads.iter().filter(|r| r["error"].is_null()).collect();
    if len == 0 || placed.is_empty() {
        return;
    }
    let font = FontId::monospace(12.0);
    let (line, cw) = ui.fonts_mut(|f| (f.row_height(&font), f.glyph_width(&font, 'A').max(1.0)));
    // As many columns per block as the width allows.
    let cols = crate::alignment_columns(ui.available_width() - GRID_NAME - 28.0, cw);
    // The blocks worth drawing: those at least one read reaches.
    let covered: Vec<usize> = (0..len.div_ceil(cols))
        .filter(|b| {
            let from = b * cols;
            (from..(from + cols).min(len))
                .any(|p| placed.iter().any(|r| base_at(r, len, p).is_some()))
        })
        .collect();
    if covered.is_empty() {
        return;
    }
    let header = format!(
        "Alignment ({} of {} reads, {} of the construct)",
        placed.len(),
        reads.len(),
        match (covered.first(), covered.last()) {
            (Some(a), Some(b)) => format!(
                "{}..{}",
                grouped((a * cols + 1) as u64),
                grouped(((b + 1) * cols).min(len) as u64)
            ),
            _ => String::new(),
        }
    );
    egui::CollapsingHeader::new(RichText::new(header).strong())
        .id_salt("reads-alignment")
        .default_open(true)
        .show(ui, |ui| {
            ui.checkbox(&mut inputs.single_line, "One line")
                .on_hover_text("Run the reads along one line and scroll sideways, rather than wrapping them into blocks");
            if inputs.single_line {
                one_line(
                    ui, t, residues, &placed, findings, inputs, &font, line, cw, height,
                );
                return;
            }
            let block_h = line * (placed.len() + 2) as f32 + line;
            let width = GRID_NAME + cw * cols as f32 + 12.0;
            // A difference chosen on the coverage map or in the list brings its block into view.
            let offset = inputs.scroll_to_selected.then(|| {
                let p = inputs
                    .selected
                    .and_then(|i| findings.get(i))
                    .and_then(finding_at)
                    .and_then(|p| usize::try_from(p).ok())?;
                let k = covered.iter().position(|&b| b == p / cols)?;
                Some(block_h * k as f32)
            });
            inputs.scroll_to_selected = false;
            let mut area = egui::ScrollArea::vertical()
                .id_salt("reads-alignment-rows")
                .max_height(height);
            if let Some(Some(y)) = offset {
                area = area.vertical_scroll_offset(y);
            }
            area.show_rows(ui, block_h, covered.len(), |ui, range| {
                for k in range {
                    let Some(&b) = covered.get(k) else { continue };
                    let from = b * cols;
                    let (rect, response) =
                        ui.allocate_exact_size(vec2(width, block_h), Sense::click());
                    block(
                        ui,
                        t,
                        residues,
                        &placed,
                        findings,
                        inputs,
                        (from, (from + cols).min(len)),
                        &font,
                        line,
                        cw,
                        rect,
                        &response,
                        rect.left(),
                        rect.left() + GRID_NAME,
                    );
                }
            });
        });
}

/// The alignment as one line, scrolled sideways: only the residues in view are drawn, so the
/// cost does not grow with the construct. The stretch the reads lie on is what can be reached,
/// since the rest is empty columns.
#[allow(clippy::too_many_arguments)]
fn one_line(
    ui: &mut Ui,
    t: &Tokens,
    residues: &[u8],
    placed: &[&Value],
    findings: &[Value],
    inputs: &mut ReadsInputs,
    font: &FontId,
    line: f32,
    cw: f32,
    height: f32,
) {
    let len = residues.len();
    let (first, last) = match (
        placed
            .iter()
            .filter_map(|r| r.get("start").and_then(Value::as_u64))
            .min(),
        placed
            .iter()
            .filter_map(|r| r.get("end").and_then(Value::as_u64))
            .max(),
    ) {
        (Some(a), Some(b)) if b > a => (
            usize::try_from(a).unwrap_or(0),
            usize::try_from(b).unwrap_or(len).min(len),
        ),
        _ => (0, len),
    };
    let shown = last.saturating_sub(first);
    let rows_h = line * (placed.len() + 2) as f32 + line;
    let width = GRID_NAME + cw * shown as f32 + 12.0;
    let mut area = egui::ScrollArea::horizontal()
        .id_salt("reads-one-line")
        .auto_shrink([false, true])
        .max_height(rows_h.min(height).max(line * 3.0));
    // A difference chosen elsewhere scrolls to it, a little in from the edge.
    if inputs.scroll_to_selected
        && let Some(p) = inputs
            .selected
            .and_then(|i| findings.get(i))
            .and_then(finding_at)
            .and_then(|p| usize::try_from(p).ok())
        && p >= first
    {
        area = area.horizontal_scroll_offset((cw * (p - first) as f32 - 200.0).max(0.0));
    }
    inputs.scroll_to_selected = false;
    area.show_viewport(ui, |ui, viewport| {
        let (rect, response) = ui.allocate_exact_size(vec2(width, rows_h), Sense::click());
        // The residues the viewport covers, with a column either side.
        let from = ((viewport.min.x - GRID_NAME) / cw).floor().max(0.0) as usize;
        let from = first.saturating_add(from).min(last);
        let count = (viewport.width() / cw).ceil() as usize + 2;
        let to = from.saturating_add(count).min(last);
        let view_left = rect.left() + viewport.min.x;
        block(
            ui,
            t,
            residues,
            placed,
            findings,
            inputs,
            (from, to),
            font,
            line,
            cw,
            Rect::from_min_size(pos2(view_left, rect.top()), vec2(viewport.width(), rows_h)),
            &response,
            // The names stay against the left edge while the residues run past them.
            view_left,
            rect.left() + GRID_NAME + cw * (from - first) as f32,
        );
    });
}

/// One block of [`alignment_grid`]: the positions, the construct, then each read.
#[allow(clippy::too_many_arguments)]
fn block(
    ui: &mut Ui,
    t: &Tokens,
    residues: &[u8],
    placed: &[&Value],
    findings: &[Value],
    inputs: &mut ReadsInputs,
    (from, to): (usize, usize),
    font: &FontId,
    line: f32,
    cw: f32,
    // The space to draw in, what the reader clicked on, and where the row names go and column
    // `from` starts (the caller may have scrolled the line sideways under them).
    rect: Rect,
    response: &egui::Response,
    names_x: f32,
    origin: f32,
) {
    let len = residues.len();
    let painter = ui.painter_at(rect);
    let x_of = |p: usize| origin + cw * (p - from) as f32;
    let letter = |p: usize, y: f32, c: char, color: egui::Color32| {
        painter.text(
            pos2(x_of(p) + cw / 2.0, y),
            Align2::CENTER_TOP,
            c,
            font.clone(),
            color,
        );
    };
    // The column of the chosen difference, marked behind everything else.
    if let Some(p) = inputs
        .selected
        .and_then(|i| findings.get(i))
        .and_then(finding_at)
        .and_then(|p| usize::try_from(p).ok())
        && (from..to).contains(&p)
    {
        painter.rect_filled(
            egui::Rect::from_min_size(
                pos2(x_of(p), rect.top() + line),
                vec2(cw, line * (placed.len() + 1) as f32),
            ),
            2.0,
            t.highlight.gamma_multiply(0.25),
        );
    }
    // The first position of the block, and every tenth after it.
    let mut y = rect.top();
    for p in (from..to).step_by(10) {
        painter.text(
            pos2(x_of(p), y),
            Align2::LEFT_TOP,
            grouped((p + 1) as u64),
            FontId::proportional(9.0),
            t.text_muted,
        );
    }
    y += line;
    painter.text(
        pos2(names_x, y),
        Align2::LEFT_TOP,
        "construct",
        font.clone(),
        t.text_muted,
    );
    for p in from..to {
        let c = residues
            .get(p)
            .map_or('?', |b| b.to_ascii_uppercase() as char);
        letter(p, y, c, t.text);
    }
    for read in placed {
        y += line;
        let name: String = read["name"]
            .as_str()
            .unwrap_or("")
            .chars()
            .take(17)
            .collect();
        painter.text(
            pos2(names_x, y),
            Align2::LEFT_TOP,
            &name,
            font.clone(),
            t.text_muted,
        );
        for p in from..to {
            let Some((base, quality)) = base_at(read, len, p) else {
                continue;
            };
            let same = residues
                .get(p)
                .is_some_and(|r| r.eq_ignore_ascii_case(&base));
            let color = if base == b'-' || !same {
                t.error
            } else {
                t.text
            };
            // A poor call is drawn faint, whether or not it differs.
            let color = match quality {
                Some(q) if q < FAINT_BELOW => color.gamma_multiply(0.45),
                _ => color,
            };
            letter(p, y, base as char, color);
            if inserts_before(read, len, p) {
                // Extra bases have no column of their own: mark the join before this one.
                painter.line_segment(
                    [pos2(x_of(p), y + 1.0), pos2(x_of(p), y + line - 1.0)],
                    Stroke::new(1.5, t.accent),
                );
            }
        }
    }
    // Clicking a column picks the difference there.
    if response.clicked()
        && let Some(pos) = response.interact_pointer_pos()
        && pos.x >= origin
    {
        let col = ((pos.x - origin) / cw).floor();
        if col >= 0.0 {
            let p = from.saturating_add(col as usize);
            if let Some(i) = findings
                .iter()
                .position(|f| finding_at(f).and_then(|x| usize::try_from(x).ok()) == Some(p))
            {
                inputs.selected = Some(i);
            }
        }
    }
}

fn read_table(ui: &mut Ui, t: &Tokens, reads: &[Value], reference: u64) -> Option<ReadsAction> {
    let mut action = None;
    ui.label(RichText::new("Reads").strong());
    egui::Grid::new("reads-table")
        .num_columns(6)
        .striped(true)
        .spacing(vec2(16.0, 4.0))
        .show(ui, |ui| {
            for h in ["Read", "Strand", "On the construct", "Identity", "Kept", ""] {
                ui.label(RichText::new(h).color(t.text_muted));
            }
            ui.end_row();
            for read in reads {
                ui.label(read["name"].as_str().unwrap_or(""));
                let open = || -> Option<ReadsAction> {
                    if let Some(i) = source(read, "matched").and_then(Value::as_u64) {
                        return Some(ReadsAction::Run(
                            "open_matched_read",
                            json!({"index": i, "document": reference}),
                        ));
                    }
                    source(read, "file")
                        .and_then(Value::as_str)
                        .map(|p| ReadsAction::Run("open", json!({"path": p})))
                };
                if let Some(e) = read["error"].as_str() {
                    ui.label(RichText::new(format!("not placed: {e}")).color(t.error));
                    ui.label("");
                    ui.label("");
                    ui.label("");
                    if ui.small_button("Open").clicked() {
                        action = open();
                    }
                    ui.end_row();
                    continue;
                }
                ui.label(if read["strand"] == "forward" {
                    "forward"
                } else {
                    "reverse"
                });
                ui.label(format!(
                    "{}..{} ({} bp)",
                    grouped(read["start"].as_u64().unwrap_or(0) + 1),
                    grouped(read["end"].as_u64().unwrap_or(0)),
                    grouped(read["length"].as_u64().unwrap_or(0)),
                ));
                ui.label(format!(
                    "{:.1}%",
                    read["identity"].as_f64().unwrap_or(0.0) * 100.0
                ));
                ui.label(match read["kept"].as_array() {
                    Some(k) => format!(
                        "{}..{} of {}",
                        k.first().and_then(Value::as_u64).unwrap_or(0) + 1,
                        k.get(1).and_then(Value::as_u64).unwrap_or(0),
                        read["read_length"].as_u64().unwrap_or(0)
                    ),
                    None => format!("all {}", read["read_length"].as_u64().unwrap_or(0)),
                });
                if source(read, "document").is_some() {
                    ui.label(RichText::new("open in a tab").color(t.text_muted));
                } else if ui
                    .small_button("Open")
                    .on_hover_text(
                        "Open the read in a tab of its own (traces with their chromatogram)",
                    )
                    .clicked()
                {
                    action = open();
                }
                ui.end_row();
            }
        });
    action
}

/// "1,843  T to A", "400  +GG inserted before it", "600..601  TC missing".
fn describe(f: &Value) -> String {
    match f.get("kind").and_then(Value::as_str) {
        Some("substitution") => format!(
            "{}  {} \u{2192} {}",
            grouped(f["position"].as_u64().unwrap_or(0) + 1),
            f["reference"].as_str().unwrap_or("").to_ascii_uppercase(),
            f["query"].as_str().unwrap_or("").to_ascii_uppercase()
        ),
        Some("insertion") => format!(
            "{}  +{} inserted before it",
            grouped(f["position"].as_u64().unwrap_or(0) + 1),
            f["bases"].as_str().unwrap_or("").to_ascii_uppercase()
        ),
        Some("deletion") => {
            let start = f["start"].as_u64().unwrap_or(0);
            let n = f["length"].as_u64().unwrap_or(0);
            format!(
                "{}..{}  {} missing",
                grouped(start + 1),
                grouped(start + n),
                f["bases"].as_str().unwrap_or("").to_ascii_uppercase()
            )
        }
        _ => String::new(),
    }
}

fn finding_list(
    ui: &mut Ui,
    t: &Tokens,
    findings: &[Value],
    reads: &[Value],
    inputs: &mut ReadsInputs,
) {
    let hidden = findings.iter().filter(|f| !visible(f, false)).count();
    ui.horizontal(|ui| {
        ui.label(RichText::new(format!("Differences ({})", findings.len())).strong());
        if hidden > 0 {
            ui.checkbox(
                &mut inputs.show_all,
                format!("Also the {hidden} seen only in calls below Q{GOOD}"),
            )
            .on_hover_text(
                "Calls below Q20 (more than 1 error in 100) are often noise near a read's ends",
            );
        }
    });
    if findings.is_empty() {
        ui.label(
            RichText::new("The reads match the construct wherever they cover it.")
                .color(t.text_muted),
        );
        return;
    }
    for (i, f) in findings.iter().enumerate() {
        if !visible(f, inputs.show_all) {
            continue;
        }
        let (verdict, color) = match f.get("verdict").and_then(Value::as_str) {
            Some("confirmed") => ("confirmed", t.accent_text),
            Some("disputed") => ("disputed", t.error),
            _ => ("one read", t.text_muted),
        };
        let shown_by: Vec<String> = f["reads"]
            .as_array()
            .into_iter()
            .flatten()
            .map(|x| {
                let name = x["read"]
                    .as_u64()
                    .and_then(|k| reads.get(usize::try_from(k).ok()?))
                    .and_then(|r| r["name"].as_str())
                    .unwrap_or("?");
                match x["quality"].as_u64() {
                    Some(q) => format!("{name} Q{q}"),
                    None => name.to_string(),
                }
            })
            .collect();
        let covering = f["covered_by"].as_array().map_or(0, Vec::len);
        ui.horizontal_wrapped(|ui| {
            let selected = inputs.selected == Some(i);
            if ui
                .selectable_label(selected, RichText::new(describe(f)).monospace())
                .on_hover_text("Show the traces here")
                .clicked()
            {
                inputs.selected = if selected { None } else { Some(i) };
                inputs.scroll_to_selected = true;
            }
            ui.label(RichText::new(verdict).color(color));
            ui.label(
                RichText::new(format!(
                    "shown by {} of {covering}: {}",
                    shown_by.len(),
                    shown_by.join(", ")
                ))
                .color(t.text_muted),
            );
        });
    }
}

/// The chromatogram of every read covering a difference, as the construct's top strand reads.
fn traces_at(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    f: &Value,
    reads: &[Value],
    exp: &Experiment,
    inputs: &mut ReadsInputs,
) {
    ui.label(
        RichText::new(format!("Traces at {}", describe(f)))
            .strong()
            .monospace(),
    );
    let shows: Vec<u64> = f["reads"]
        .as_array()
        .into_iter()
        .flatten()
        .filter_map(|x| x["read"].as_u64())
        .collect();
    for c in f["covered_by"].as_array().into_iter().flatten() {
        let (Some(k), Some(pos)) = (c["read"].as_u64(), c["read_position"].as_u64()) else {
            continue;
        };
        let Some(read) = usize::try_from(k).ok().and_then(|k| reads.get(k)) else {
            continue;
        };
        if let Some(i) = source(read, "matched").and_then(Value::as_u64) {
            let reverse = read["strand"] == "reverse";
            ui.label(
                RichText::new(format!(
                    "{}  {}",
                    read["name"].as_str().unwrap_or(""),
                    if shows.contains(&k) {
                        "shows the difference"
                    } else {
                        "shows the construct's base"
                    }
                ))
                .color(t.text_muted),
            );
            let width = ui.available_width().min(520.0);
            let m = exp.reference.and_then(|id| {
                inputs
                    .matched_reads(session, id)
                    .get(usize::try_from(i).ok()?)
                    .and_then(|r| r.molecule.clone().ok())
            });
            if let Some(m) = m {
                crate::trace_view::snippet(
                    ui,
                    t,
                    &m,
                    usize::try_from(pos).unwrap_or(0),
                    10,
                    reverse,
                    vec2(width, 90.0),
                );
            }
            continue;
        }
        if let Some(id) = source(read, "document").and_then(Value::as_u64) {
            let reverse = read["strand"] == "reverse";
            ui.label(
                RichText::new(format!(
                    "{}  {}",
                    read["name"].as_str().unwrap_or(""),
                    if shows.contains(&k) {
                        "shows the difference"
                    } else {
                        "shows the construct's base"
                    }
                ))
                .color(t.text_muted),
            );
            let width = ui.available_width().min(520.0);
            if let Some(m) = session.document_by_id(id).map(|d| d.molecule.clone()) {
                crate::trace_view::snippet(
                    ui,
                    t,
                    &m,
                    usize::try_from(pos).unwrap_or(0),
                    10,
                    reverse,
                    vec2(width, 90.0),
                );
            }
            continue;
        }
        let Some(path) = read
            .get("source")
            .and_then(|src| src.get("file"))
            .and_then(Value::as_str)
            .map(str::to_string)
            .or_else(|| {
                usize::try_from(k)
                    .ok()
                    .and_then(|k| exp.reads.get(k).cloned())
            })
        else {
            continue;
        };
        let reverse = read["strand"] == "reverse";
        let label = format!(
            "{}  {}",
            read["name"].as_str().unwrap_or(""),
            if shows.contains(&k) {
                "shows the difference"
            } else {
                "shows the construct's base"
            }
        );
        ui.label(RichText::new(label).color(t.text_muted));
        let width = ui.available_width().min(520.0);
        match inputs.molecule(session, &path) {
            Some(m) => {
                let m = m.clone();
                crate::trace_view::snippet(
                    ui,
                    t,
                    &m,
                    usize::try_from(pos).unwrap_or(0),
                    10,
                    reverse,
                    vec2(width, 90.0),
                );
            }
            None => {
                ui.label(RichText::new(format!("{path} can't be read")).color(t.error));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quality_filter() {
        let poor = json!({"reads": [{"quality": 8}, {"quality": 12}]});
        let good = json!({"reads": [{"quality": 8}, {"quality": 31}]});
        let plain = json!({"reads": [{"quality": null}]});
        assert!(!visible(&poor, false) && visible(&poor, true));
        assert!(visible(&good, false));
        assert!(visible(&plain, false));
        // The first difference shown opens by itself: a disputed one that passes the filter.
        let mut inputs = ReadsInputs::default();
        inputs.set_result(Ok(json!({"findings": [
            {"verdict": "disputed", "reads": [{"quality": 5}]},
            {"verdict": "one_read", "reads": [{"quality": 40}]},
            {"verdict": "disputed", "reads": [{"quality": 30}]},
        ]})));
        assert_eq!(inputs.selected, Some(2));
        assert_eq!(grouped(4124), "4,124");
    }

    /// A row laid out against the construct is read column by column, and wraps the origin.
    #[test]
    fn a_row_is_read_in_the_constructs_coordinates() {
        let read = json!({
            "start": 8, "row": "ACG-T", "row_quality": [40, 9, 40, null, 40],
            "row_insertions": [2],
        });
        assert_eq!(base_at(&read, 10, 8), Some((b'A', Some(40))));
        assert_eq!(base_at(&read, 10, 9), Some((b'C', Some(9))));
        // Past the origin of a 10 bp circle.
        assert_eq!(base_at(&read, 10, 0), Some((b'G', Some(40))));
        assert_eq!(base_at(&read, 10, 1), Some((b'-', None)));
        assert_eq!(base_at(&read, 10, 3), None);
        // The extra bases sit before the column at offset 2.
        assert!(inserts_before(&read, 10, 0));
        assert!(!inserts_before(&read, 10, 1));
        // A read that was never placed has no row at all.
        assert_eq!(base_at(&json!({"start": 0}), 10, 0), None);
        assert!(!inserts_before(&json!({"start": 0}), 10, 0));
    }

    /// The alignment takes a column of its own once the page is wide enough for one, and the
    /// summary's column never grows past what its own content needs.
    #[test]
    fn the_alignment_gets_its_own_column_on_a_wide_page() {
        assert_eq!(crate::two_columns(900.0), None, "stacked under the summary");
        assert_eq!(crate::two_columns(999.0), None);
        let left = crate::two_columns(1400.0).expect("two columns");
        assert!((380.0..=620.0).contains(&left), "{left}");
        // Whatever the page's width, the alignment keeps room for a readable block.
        for full in [1000.0, 1400.0, 2400.0, 4000.0] {
            let left = crate::two_columns(full).expect("two columns");
            assert!(
                full - left >= 380.0,
                "{full} wide leaves the alignment only {}",
                full - left
            );
        }
    }

    /// Reads that are already open in tabs are offered as reads, counted towards being able to
    /// go on, and sent to the check. In the browser this is the only way to add a read until a
    /// file has been loaded, and dropping a trace on the page opens it as a tab.
    #[test]
    fn reads_can_be_tabs_as_well_as_files() {
        use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
        use snappiesnap_engine::{Document, Session};

        let construct =
            "ACGTACGTAAAACCCCTTTTGGGGACGTACGTAAAACCCCTTTTGGGGACGTTTTTAAAACCCCGGGGTTTTAAAA";
        let doc = |name: &str, residues: &str| {
            Document::new(
                Molecule::new(
                    name,
                    Seq::new(Alphabet::Dna, residues).expect("test sequence"),
                    Topology::Linear,
                )
                .expect("test molecule"),
            )
        };
        let mut session = Session::new();
        session.open_document(doc("pCheck", construct));
        session.open_document(doc(
            "read from a tab",
            construct.get(10..70).expect("part of the construct"),
        ));
        let mut exp = Experiment {
            reference: Some(1),
            ..Experiment::new(crate::state::ExperimentKind::Reads, "Check Reads", [1u64])
        };
        let mut inputs = ReadsInputs::default();
        let t = Tokens::dark();
        // Step 1 with nothing chosen: the second tab is offered, and there is nothing to run.
        crate::offscreen::render(
            900,
            600,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    setup(ui, &t, &session, &mut exp, &mut inputs, None);
                });
            },
        )
        .expect("step 1 draws");
        assert_eq!(exp.step, ExperimentStep::Constructs, "nothing to check yet");
        // Tick the tab, as the checkbox does, and the check takes it.
        exp.read_documents.push(2);
        exp.step = ExperimentStep::Results;
        let mut asked = None;
        crate::offscreen::render(
            1200,
            700,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    // The check is asked for on the first pass; later passes find it pending.
                    let got = results(ui, &t, &session, &mut exp, &mut inputs);
                    if asked.is_none() {
                        asked = got;
                    }
                });
            },
        )
        .expect("step 2 draws");
        let Some(ReadsAction::Check(params)) = asked else {
            panic!("the check should have been asked for: {asked:?}");
        };
        assert_eq!(params["documents"], json!([2]));
        // Run it as the app does, and the read is placed.
        let answer = session
            .execute("check_reads", params)
            .map_err(|e| e.message);
        let answer = answer.expect("a tab checks like a file");
        assert_eq!(answer["placed"], 1);
        assert_eq!(answer["reads"][0]["name"], "read from a tab");
        assert_eq!(answer["reads"][0]["source"]["document"], 2);
        // A tab closed since it was ticked is dropped rather than checked as a missing id.
        exp.step = ExperimentStep::Constructs;
        let mut session = Session::new();
        session.open_document(doc("pCheck", construct));
        crate::offscreen::render(
            900,
            600,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                egui::Frame::new().fill(t.panel).show(ui, |ui| {
                    setup(ui, &t, &session, &mut exp, &mut inputs, None);
                });
            },
        )
        .expect("step 1 draws again");
        assert!(exp.read_documents.is_empty());
    }

    /// Three reads of one construct, drawn in both themes, and wide enough for the alignment
    /// to get a column of its own beside the summary. Set `SNAPPIESNAP_READS_PNG=/p.png` to
    /// look at it (THEME in the name is replaced).
    #[test]
    fn the_alignment_grid_draws() {
        use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
        use snappiesnap_engine::{Document, Session};

        let mut x: u32 = 99;
        let construct: String = (0..400)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
                ['A', 'C', 'G', 'T'][((x >> 16) & 3) as usize]
            })
            .collect();
        let mut session = Session::new();
        let doc = |name: &str, residues: &str| {
            Document::new(
                Molecule::new(
                    name,
                    Seq::new(Alphabet::Dna, residues).expect("test sequence"),
                    Topology::Linear,
                )
                .expect("test molecule"),
            )
        };
        session.open_document(doc("pCheck", &construct));
        let part = construct.get(40..340).expect("the middle of the construct");
        // Two reads agree on one substitution; the third carries the construct's base there.
        let mut changed = part.to_string();
        changed.replace_range(100..101, "A");
        session.open_document(doc("read A", &changed));
        session.open_document(doc("read B", &changed));
        session.open_document(doc("read C", part));
        let answer = session
            .execute(
                "check_reads",
                json!({"document": 1, "documents": [2, 3, 4], "rows": true}),
            )
            .map_err(|e| e.message);
        let r = answer.as_ref().expect("the check runs");
        assert_eq!(r["placed"], 3);
        let row = r["reads"][0]["row"].as_str().expect("a row per read");
        assert_eq!(row.len(), 300);
        assert_eq!(r["reads"][0]["start"], 40);
        assert_eq!(r["findings"][0]["verdict"], "disputed");
        let mut inputs = ReadsInputs::default();
        inputs.set_result(answer);
        let mut exp = Experiment {
            reference: Some(1),
            step: ExperimentStep::Results,
            ..Experiment::new(crate::state::ExperimentKind::Reads, "Check Reads", [1u64])
        };
        inputs.key = Some(format!(
            "{}|1|{:?}|{}",
            session.revision(),
            exp.reads,
            exp.matched
        ));
        // One line as well as blocks: the same state, so neither can go stale.
        inputs.single_line = true;
        crate::offscreen::render(
            1400,
            760,
            1.0,
            |ctx| Tokens::dark().apply(ctx),
            |ui| {
                egui::Frame::new()
                    .fill(Tokens::dark().panel)
                    .show(ui, |ui| {
                        results(ui, &Tokens::dark(), &session, &mut exp, &mut inputs);
                    });
            },
        )
        .expect("one line draws");
        inputs.single_line = false;
        // Wide: the alignment beside the summary. Narrow: stacked under it. Both are drawn
        // from the same state, so neither layout can go stale.
        for (name, t) in [("dark", Tokens::dark()), ("light", Tokens::light())] {
            for (tag, w) in [("wide", 1400), ("narrow", 860)] {
                let image = crate::offscreen::render(
                    w,
                    760,
                    1.0,
                    |ctx| t.apply(ctx),
                    |ui| {
                        egui::Frame::new().fill(t.panel).show(ui, |ui| {
                            results(ui, &t, &session, &mut exp, &mut inputs);
                        });
                    },
                )
                .expect("the panel draws");
                if let Ok(path) = std::env::var("SNAPPIESNAP_READS_PNG") {
                    crate::offscreen::write_png(
                        &path.replace("THEME", &format!("{name}-{tag}")),
                        &image,
                    )
                    .expect("writing the png");
                }
            }
        }
    }
}
