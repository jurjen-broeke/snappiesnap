//! The folder panel: the contents of a folder as a tree of sequence files, with a search box,
//! a "···" menu, right-click menus, inline rename and new folder, and drag-and-drop moves.
//!
//! Listings go through the engine's [`Platform`], so the panel works on any file system the
//! session uses. It never changes files itself: everything it asks for is a [`SidebarAction`]
//! that the shell runs as an engine command (`open`, `file_rename`, `file_copy`,
//! `folder_create`), exactly as an agent would.

use std::collections::{BTreeSet, HashMap};

use eframe::egui::{self, Align2, FontId, RichText, Sense, Stroke, Ui, pos2, vec2};
use snappiesnap_engine::{DirEntry, Platform};

use crate::theme::Tokens;

/// File extensions the panel shows.
pub const SEQUENCE_EXTENSIONS: &[&str] = &[
    "dna", "gb", "gbk", "genbank", "gbff", "embl", "emb", "gff3", "gff", "sbol", "fa", "fasta",
    "fna", "fas", "ffn", "ab1", "abi", "scf",
];

/// Folders searched at most, and how deep, when filtering or drawing the tree.
const SEARCH_FOLDERS: usize = 2_000;
const MAX_DEPTH: usize = 8;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SidebarAction {
    Open(String),
    ChooseFolder,
    Rename {
        from: String,
        to: String,
    },
    Copy {
        from: String,
        to: String,
    },
    /// Move `from` into the folder `into`.
    Move {
        from: String,
        into: String,
    },
    NewFolder(String),
    Reveal(String),
    CopyPath(String),
    Hide,
}

/// What's being edited inline.
#[derive(Debug, Clone, PartialEq, Eq)]
enum Edit {
    Rename { path: String, text: String },
    NewFolder { parent: String, text: String },
}

/// What's dragged inside the panel.
#[derive(Debug, Clone)]
struct Dragged(String);

/// Panel state: the search text, expanded folders, inline edit and cached listings.
#[derive(Debug, Default)]
pub struct Sidebar {
    pub filter: String,
    pub expanded: BTreeSet<String>,
    edit: Option<Edit>,
    listings: HashMap<String, Result<Vec<DirEntry>, String>>,
}

/// A SnappieSnap experiment file (`.sslab`, or `.ssexp` as 0.6.0 named them).
pub fn is_experiment_file(name: &str) -> bool {
    use snappiesnap_engine::experiment::{EXTENSION, OLD_EXTENSIONS};
    name.rsplit_once('.').is_some_and(|(_, ext)| {
        std::iter::once(EXTENSION)
            .chain(OLD_EXTENSIONS.iter().copied())
            .any(|e| ext.eq_ignore_ascii_case(e))
    })
}

/// Files the panel lists: sequences and experiments.
fn is_listed(name: &str) -> bool {
    is_sequence_file(name) || is_experiment_file(name)
}

pub fn is_sequence_file(name: &str) -> bool {
    name.rsplit_once('.')
        .is_some_and(|(_, ext)| SEQUENCE_EXTENSIONS.contains(&ext.to_ascii_lowercase().as_str()))
}

/// The folder containing `path` (handles `/` and `\`).
pub fn parent_folder(path: &str) -> Option<String> {
    let (dir, _) = path.rsplit_once(['/', '\\'])?;
    Some(if dir.is_empty() {
        "/".to_string()
    } else {
        dir.to_string()
    })
}

/// The last component of a path.
pub fn base_name(path: &str) -> &str {
    path.trim_end_matches(['/', '\\'])
        .rsplit(['/', '\\'])
        .next()
        .filter(|n| !n.is_empty())
        .unwrap_or(path)
}

/// `dir` joined with `name`, using the separator `dir` already uses.
pub fn join(dir: &str, name: &str) -> String {
    let sep = if dir.contains('\\') && !dir.contains('/') {
        '\\'
    } else {
        '/'
    };
    if dir.ends_with(['/', '\\']) {
        format!("{dir}{name}")
    } else {
        format!("{dir}{sep}{name}")
    }
}

/// Whether `path` is `dir` or inside it.
pub fn is_within(path: &str, dir: &str) -> bool {
    let dir = dir.trim_end_matches(['/', '\\']);
    path == dir
        || path
            .strip_prefix(dir)
            .is_some_and(|rest| rest.starts_with(['/', '\\']))
}

/// A name that is usable as a single path component.
fn valid_name(name: &str) -> bool {
    let n = name.trim();
    !n.is_empty()
        && n != "."
        && n != ".."
        && !n.contains(['/', '\\'])
        && n.chars().all(|c| !c.is_control())
}

/// "name copy.ext", or "name copy 2.ext", ... not in `taken`.
pub fn copy_name(name: &str, taken: &[String]) -> String {
    let (stem, ext) = match name.rsplit_once('.') {
        Some((s, e)) if !s.is_empty() => (s, format!(".{e}")),
        _ => (name, String::new()),
    };
    let mut candidate = format!("{stem} copy{ext}");
    let mut n = 2;
    while taken.contains(&candidate) && n < 10_000 {
        candidate = format!("{stem} copy {n}{ext}");
        n += 1;
    }
    candidate
}

impl Sidebar {
    /// Forget cached listings (after file changes, or when the user asks).
    pub fn refresh(&mut self) {
        self.listings.clear();
    }

    fn listing(&mut self, platform: &dyn Platform, path: &str) -> Result<Vec<DirEntry>, String> {
        self.listings
            .entry(path.to_string())
            .or_insert_with(|| platform.list_dir(path).map_err(|e| e.to_string()))
            .clone()
    }

    fn names_in(&mut self, platform: &dyn Platform, dir: &str) -> Vec<String> {
        self.listing(platform, dir)
            .map(|v| v.into_iter().map(|e| e.name).collect())
            .unwrap_or_default()
    }

    /// Start creating a folder inside `parent` (shown as an inline text field).
    pub fn start_new_folder(&mut self, parent: &str) {
        self.expanded.insert(parent.to_string());
        self.edit = Some(Edit::NewFolder {
            parent: parent.to_string(),
            text: "New folder".into(),
        });
    }

    /// Sequence files under `root` whose name contains `filter` (case-insensitive), with paths
    /// relative to `root`. Bounded in folders visited and depth.
    pub fn search(
        &mut self,
        platform: &dyn Platform,
        root: &str,
        filter: &str,
    ) -> Vec<(String, String)> {
        let needle = filter.to_lowercase();
        let mut out = Vec::new();
        let mut queue = vec![(root.to_string(), String::new(), 0usize)];
        let mut visited = 0usize;
        while let Some((dir, rel, depth)) = queue.pop() {
            visited += 1;
            if visited > SEARCH_FOLDERS {
                break;
            }
            let Ok(entries) = self.listing(platform, &dir) else {
                continue;
            };
            for e in entries {
                let rel_name = if rel.is_empty() {
                    e.name.clone()
                } else {
                    format!("{rel}/{}", e.name)
                };
                if e.is_dir {
                    if depth < MAX_DEPTH {
                        queue.push((e.path.clone(), rel_name, depth + 1));
                    }
                } else if is_listed(&e.name) && e.name.to_lowercase().contains(&needle) {
                    out.push((rel_name, e.path));
                }
            }
        }
        out.sort();
        out
    }
}

/// Draw the panel. `current` is the open document's path, highlighted in the tree.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    platform: &dyn Platform,
    folder: Option<&str>,
    current: Option<&str>,
    sb: &mut Sidebar,
) -> Option<SidebarAction> {
    let mut action = None;

    // Header: folder name (also a drop target for its root) and the "···" menu.
    ui.horizontal(|ui| {
        let title = folder.map_or("No folder", base_name);
        let r = ui
            .add(egui::Label::new(RichText::new(title).strong().size(15.0)).sense(Sense::hover()));
        if let Some(f) = folder {
            drop_target(ui, t, &r, f, &mut action);
            r.on_hover_text(f);
        }
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.menu_button("···", |ui| {
                if ui.button("Choose folder…").clicked() {
                    action = Some(SidebarAction::ChooseFolder);
                    ui.close();
                }
                if let Some(f) = folder {
                    if ui.button("New folder").clicked() {
                        sb.start_new_folder(f);
                        ui.close();
                    }
                    if ui.button("Show in file manager").clicked() {
                        action = Some(SidebarAction::Reveal(f.to_string()));
                        ui.close();
                    }
                    if ui.button("Copy folder path").clicked() {
                        action = Some(SidebarAction::CopyPath(f.to_string()));
                        ui.close();
                    }
                }
                ui.separator();
                if ui.button("Refresh").clicked() {
                    sb.refresh();
                    ui.close();
                }
                if ui.button("Collapse all").clicked() {
                    sb.expanded.clear();
                    ui.close();
                }
                if ui.button("Hide panel").clicked() {
                    action = Some(SidebarAction::Hide);
                    ui.close();
                }
            });
        });
    });
    ui.add(
        egui::TextEdit::singleline(&mut sb.filter)
            .hint_text("Search files")
            .desired_width(f32::INFINITY),
    );
    ui.add_space(t.spacing);
    let Some(folder) = folder else {
        ui.label(RichText::new("Open a file or choose a folder.").color(t.text_muted));
        if ui.button("Choose folder…").clicked() {
            action = Some(SidebarAction::ChooseFolder);
        }
        return action;
    };
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            if sb.filter.trim().is_empty() {
                new_folder_field(ui, t, sb, folder, 0, &mut action);
                tree(ui, t, platform, folder, current, sb, 0, &mut action);
            } else {
                let filter = sb.filter.trim().to_string();
                let hits = sb.search(platform, folder, &filter);
                if hits.is_empty() {
                    ui.label(RichText::new("No matching files").color(t.text_muted));
                }
                for (rel, path) in hits {
                    let r = row(
                        ui,
                        t,
                        &rel,
                        0,
                        RowKind::File,
                        current == Some(path.as_str()),
                    );
                    file_interactions(ui, platform, sb, &r, &path, &mut action);
                }
            }
        });
    action
}

#[derive(Clone, Copy, PartialEq)]
enum RowKind {
    File,
    Folder { open: bool },
}

#[allow(clippy::too_many_arguments)]
fn tree(
    ui: &mut Ui,
    t: &Tokens,
    platform: &dyn Platform,
    dir: &str,
    current: Option<&str>,
    sb: &mut Sidebar,
    depth: usize,
    action: &mut Option<SidebarAction>,
) {
    if depth > MAX_DEPTH {
        return;
    }
    let entries = match sb.listing(platform, dir) {
        Ok(e) => e,
        Err(e) => {
            ui.label(RichText::new(format!("Can't read this folder: {e}")).color(t.error));
            return;
        }
    };
    let mut shown = 0usize;
    for e in entries {
        if !e.is_dir && !is_listed(&e.name) {
            continue;
        }
        shown += 1;
        if rename_field(ui, t, sb, &e, depth, action) {
            continue;
        }
        if e.is_dir {
            let open = sb.expanded.contains(&e.path);
            let r = row(ui, t, &e.name, depth, RowKind::Folder { open }, false);
            r.dnd_set_drag_payload(Dragged(e.path.clone()));
            if r.clicked() {
                if open {
                    sb.expanded.remove(&e.path);
                } else {
                    sb.expanded.insert(e.path.clone());
                }
            }
            drop_target(ui, t, &r, &e.path, action);
            r.context_menu(|ui| {
                if ui
                    .button(if open { "Collapse" } else { "Expand" })
                    .clicked()
                {
                    if open {
                        sb.expanded.remove(&e.path);
                    } else {
                        sb.expanded.insert(e.path.clone());
                    }
                    ui.close();
                }
                if ui.button("New folder inside").clicked() {
                    sb.start_new_folder(&e.path);
                    ui.close();
                }
                common_items(ui, sb, &e.path, action);
            });
            if sb.expanded.contains(&e.path) {
                new_folder_field(ui, t, sb, &e.path, depth + 1, action);
                tree(ui, t, platform, &e.path, current, sb, depth + 1, action);
            }
        } else {
            let r = row(
                ui,
                t,
                &e.name,
                depth,
                RowKind::File,
                current == Some(e.path.as_str()),
            );
            file_interactions(ui, platform, sb, &r, &e.path, action);
        }
    }
    if shown == 0 && depth == 0 {
        ui.label(RichText::new("No sequence files here").color(t.text_muted));
    }
}

/// Click to open, drag to move, right-click for the file menu.
fn file_interactions(
    ui: &mut Ui,
    platform: &dyn Platform,
    sb: &mut Sidebar,
    r: &egui::Response,
    path: &str,
    action: &mut Option<SidebarAction>,
) {
    let _ = ui;
    r.dnd_set_drag_payload(Dragged(path.to_string()));
    if r.clicked() {
        *action = Some(SidebarAction::Open(path.to_string()));
    }
    r.context_menu(|ui| {
        if ui.button("Open").clicked() {
            *action = Some(SidebarAction::Open(path.to_string()));
            ui.close();
        }
        if ui.button("Duplicate").clicked() {
            if let Some(dir) = parent_folder(path) {
                let taken = sb.names_in(platform, &dir);
                let to = join(&dir, &copy_name(base_name(path), &taken));
                *action = Some(SidebarAction::Copy {
                    from: path.to_string(),
                    to,
                });
            }
            ui.close();
        }
        common_items(ui, sb, path, action);
    });
}

/// Menu items for both files and folders.
fn common_items(ui: &mut Ui, sb: &mut Sidebar, path: &str, action: &mut Option<SidebarAction>) {
    if ui.button("Rename…").clicked() {
        sb.edit = Some(Edit::Rename {
            path: path.to_string(),
            text: base_name(path).to_string(),
        });
        ui.close();
    }
    ui.separator();
    if ui.button("Show in file manager").clicked() {
        *action = Some(SidebarAction::Reveal(path.to_string()));
        ui.close();
    }
    if ui.button("Copy path").clicked() {
        *action = Some(SidebarAction::CopyPath(path.to_string()));
        ui.close();
    }
}

/// Accept files or folders dragged onto a folder row.
fn drop_target(
    ui: &Ui,
    t: &Tokens,
    r: &egui::Response,
    dir: &str,
    action: &mut Option<SidebarAction>,
) {
    if let Some(dragged) = r.dnd_hover_payload::<Dragged>()
        && can_move(&dragged.0, dir)
    {
        ui.painter().rect_stroke(
            r.rect,
            t.radius,
            Stroke::new(1.5, t.accent),
            egui::StrokeKind::Inside,
        );
    }
    if let Some(dragged) = r.dnd_release_payload::<Dragged>()
        && can_move(&dragged.0, dir)
    {
        *action = Some(SidebarAction::Move {
            from: dragged.0.clone(),
            into: dir.to_string(),
        });
    }
}

/// A move that changes something and doesn't put a folder inside itself.
fn can_move(from: &str, into: &str) -> bool {
    parent_folder(from).as_deref() != Some(into) && !is_within(into, from)
}

/// An inline rename field in place of a row; returns true if it was drawn.
fn rename_field(
    ui: &mut Ui,
    t: &Tokens,
    sb: &mut Sidebar,
    e: &DirEntry,
    depth: usize,
    action: &mut Option<SidebarAction>,
) -> bool {
    let Some(Edit::Rename { path, text }) = sb.edit.as_mut() else {
        return false;
    };
    if *path != e.path {
        return false;
    }
    let (done, cancel) = edit_field(ui, t, text, depth);
    if done {
        let new_name = text.trim().to_string();
        if valid_name(&new_name)
            && new_name != e.name
            && let Some(dir) = parent_folder(&e.path)
        {
            *action = Some(SidebarAction::Rename {
                from: e.path.clone(),
                to: join(&dir, &new_name),
            });
        }
        sb.edit = None;
    } else if cancel {
        sb.edit = None;
    }
    true
}

/// The inline "new folder" field at the top of `parent`'s children, if that's being edited.
fn new_folder_field(
    ui: &mut Ui,
    t: &Tokens,
    sb: &mut Sidebar,
    parent: &str,
    depth: usize,
    action: &mut Option<SidebarAction>,
) {
    let Some(Edit::NewFolder { parent: p, text }) = sb.edit.as_mut() else {
        return;
    };
    if p != parent {
        return;
    }
    let (done, cancel) = edit_field(ui, t, text, depth);
    if done {
        let name = text.trim().to_string();
        if valid_name(&name) {
            *action = Some(SidebarAction::NewFolder(join(parent, &name)));
        }
        sb.edit = None;
    } else if cancel {
        sb.edit = None;
    }
}

/// A focused single-line edit; returns (confirmed with Enter, cancelled).
fn edit_field(ui: &mut Ui, t: &Tokens, text: &mut String, depth: usize) -> (bool, bool) {
    let mut out = (false, false);
    ui.horizontal(|ui| {
        ui.add_space(4.0 + depth as f32 * INDENT + 13.0);
        let r = ui.add(egui::TextEdit::singleline(text).desired_width(f32::INFINITY));
        if !r.has_focus() && !r.lost_focus() {
            r.request_focus();
        }
        let enter = ui.input(|i| i.key_pressed(egui::Key::Enter));
        let escape = ui.input(|i| i.key_pressed(egui::Key::Escape));
        out = (
            r.lost_focus() && enter,
            escape || (r.lost_focus() && !enter),
        );
        let _ = t;
    });
    out
}

const ROW_H: f32 = 22.0;
const INDENT: f32 = 14.0;

/// One row of the tree: disclosure triangle and folder glyph, or a file glyph tinted by format.
fn row(
    ui: &mut Ui,
    t: &Tokens,
    name: &str,
    depth: usize,
    kind: RowKind,
    current: bool,
) -> egui::Response {
    let (rect, r) =
        ui.allocate_exact_size(vec2(ui.available_width(), ROW_H), Sense::click_and_drag());
    if current {
        ui.painter().rect_filled(rect, t.radius, t.accent);
    } else if r.hovered() {
        ui.painter().rect_filled(rect, t.radius, t.border);
    }
    let x = rect.left() + 4.0 + depth as f32 * INDENT;
    let cy = rect.center().y;
    let text_color = if current { t.accent_text } else { t.text };
    match kind {
        RowKind::Folder { open } => {
            let tri = if open {
                vec![
                    pos2(x, cy - 3.0),
                    pos2(x + 8.0, cy - 3.0),
                    pos2(x + 4.0, cy + 3.0),
                ]
            } else {
                vec![
                    pos2(x + 1.0, cy - 4.0),
                    pos2(x + 7.0, cy),
                    pos2(x + 1.0, cy + 4.0),
                ]
            };
            ui.painter()
                .add(egui::Shape::convex_polygon(tri, t.text_muted, Stroke::NONE));
            let f = egui::Rect::from_min_size(pos2(x + 13.0, cy - 5.0), vec2(14.0, 10.0));
            ui.painter().rect_stroke(
                f,
                2.0,
                Stroke::new(1.2, t.text_muted),
                egui::StrokeKind::Inside,
            );
            ui.painter().text(
                pos2(x + 33.0, cy),
                Align2::LEFT_CENTER,
                name,
                FontId::proportional(13.0),
                text_color,
            );
        }
        RowKind::File => {
            let ext = name.rsplit_once('.').map_or("", |(_, e)| e);
            let tint = match ext.to_ascii_lowercase().as_str() {
                "dna" => t.feature_color(5),
                "gb" | "gbk" | "genbank" | "gbff" | "embl" | "emb" => t.feature_color(1),
                "gff3" | "gff" => t.feature_color(3),
                "sbol" => t.feature_color(4),
                "ab1" | "abi" | "scf" => t.feature_color(2),
                "sslab" | "ssexp" => t.feature_color(6),
                _ => t.feature_color(0),
            };
            let page = egui::Rect::from_min_size(pos2(x + 14.0, cy - 6.0), vec2(10.0, 12.0));
            ui.painter()
                .rect_filled(page, 1.5, if current { t.accent_text } else { tint });
            ui.painter().text(
                pos2(x + 33.0, cy),
                Align2::LEFT_CENTER,
                name,
                FontId::proportional(13.0),
                text_color,
            );
        }
    }
    r.on_hover_text(name)
}

/// [`file_picker`] for other views (no "this document" choice).
pub(crate) fn file_picker_pub(
    ui: &mut Ui,
    t: &Tokens,
    id: &str,
    value: &mut String,
    files: &[String],
) {
    file_picker(ui, t, id, value, files, None);
}

/// Pick a sequence file from the open folder (or browse); `this_doc` offers the open document
/// as the empty value.
fn file_picker(
    ui: &mut Ui,
    t: &Tokens,
    id: &str,
    value: &mut String,
    files: &[String],
    this_doc: Option<&str>,
) {
    let shown = if value.is_empty() {
        this_doc.map_or("Choose a file".to_string(), |n| {
            format!("This document ({n})")
        })
    } else {
        base_name(value).to_string()
    };
    egui::ComboBox::from_id_salt(id)
        .width(240.0)
        .selected_text(shown)
        .show_ui(ui, |ui| {
            if let Some(n) = this_doc {
                ui.selectable_value(value, String::new(), format!("This document ({n})"));
            }
            for f in files {
                ui.selectable_value(value, f.clone(), base_name(f));
            }
            if files.is_empty() {
                ui.label(RichText::new("No sequence files in the open folder").color(t.text_muted));
            }
        });
    #[cfg(not(target_arch = "wasm32"))]
    if ui.small_button("Browse\u{2026}").clicked()
        && let Some(p) = rfd::FileDialog::new().pick_file()
    {
        *value = p.to_string_lossy().into_owned();
    }
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::MemoryFs;

    use super::*;

    #[test]
    fn names_and_paths() {
        assert!(is_sequence_file("p.DNA") && is_sequence_file("x.gbk") && is_sequence_file("y.fa"));
        assert!(!is_sequence_file("notes.txt") && !is_sequence_file("dna"));
        assert!(is_experiment_file("run.sslab") && is_experiment_file("old.SSEXP"));
        assert!(!is_experiment_file("run.dna") && !is_sequence_file("run.sslab"));
        assert_eq!(parent_folder("/a/b/c.dna").as_deref(), Some("/a/b"));
        assert_eq!(parent_folder("C:\\x\\y.gb").as_deref(), Some("C:\\x"));
        assert_eq!(parent_folder("/top.dna").as_deref(), Some("/"));
        assert_eq!(parent_folder("plain.dna"), None);
        assert_eq!(base_name("/a/b/"), "b");
        assert_eq!(join("/a", "b.dna"), "/a/b.dna");
        assert_eq!(join("C:\\x", "y.gb"), "C:\\x\\y.gb");
        assert_eq!(join("/", "y.gb"), "/y.gb");
        assert!(is_within("/a/b/c", "/a/b") && is_within("/a/b", "/a/b/"));
        assert!(!is_within("/a/bc", "/a/b"));
    }

    #[test]
    fn moves_and_names() {
        assert!(can_move("/lab/a.dna", "/lab/sub"));
        assert!(!can_move("/lab/a.dna", "/lab"), "already there");
        assert!(
            !can_move("/lab/sub", "/lab/sub/deeper"),
            "folder into itself"
        );
        assert!(!can_move("/lab/sub", "/lab/sub"));
        assert!(valid_name("pUC19.dna"));
        for bad in ["", "  ", ".", "..", "a/b", "a\\b", "x\ny"] {
            assert!(!valid_name(bad), "{bad:?}");
        }
        let taken = vec!["p copy.dna".to_string(), "p copy 2.dna".to_string()];
        assert_eq!(copy_name("p.dna", &[]), "p copy.dna");
        assert_eq!(copy_name("p.dna", &taken), "p copy 3.dna");
        assert_eq!(copy_name("README", &[]), "README copy");
    }

    #[test]
    fn search_finds_sequence_files_in_subfolders() {
        let fs = MemoryFs::new()
            .with_file("lab/pUC19.dna", vec![])
            .with_file("lab/old/pUC18.gb", vec![])
            .with_file("lab/old/notes.txt", vec![])
            .with_file("lab/old/deeper/pucX.fa", vec![]);
        let mut sb = Sidebar::default();
        let hits: Vec<String> = sb
            .search(&fs, "lab", "puc")
            .into_iter()
            .map(|(r, _)| r)
            .collect();
        assert_eq!(
            hits,
            vec!["old/deeper/pucX.fa", "old/pUC18.gb", "pUC19.dna"]
        );
        assert!(sb.search(&fs, "lab", "notes").is_empty());
        assert!(sb.search(&fs, "missing", "x").is_empty());
    }
}
