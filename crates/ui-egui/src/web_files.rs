//! Files in the web build. A browser gives a page no file paths, so files the user picks or
//! drops are copied into the session's in-memory file system ([`snappiesnap_engine::MemoryFs`])
//! and opened from there with the ordinary `open` command. What the app saves there goes back
//! to the browser as a download. Nothing here touches the browser (that is [`crate::web`]), so
//! these rules are tested natively.

#![cfg_attr(not(target_arch = "wasm32"), allow(dead_code))]

use std::sync::{Arc, Mutex, PoisonError};

use eframe::egui::{self, RichText, vec2};
use snappiesnap_engine::{MAX_FILE, Platform};

use crate::theme::Tokens;

/// The longest file name kept, in characters.
const MAX_NAME: usize = 120;

/// What the file picker accepts for Open (the desktop app's Open filter).
pub const OPEN_ACCEPT: &str = ".dna,.gb,.gbk,.genbank,.embl,.emb,.gff3,.gff,.sbol,.xml,.rdf,\
.fa,.fasta,.fna,.fas,.ab1,.abi,.scf,.sslab,.ssexp";

/// A file the browser handed over: its name, and its contents or why they couldn't be read.
#[derive(Debug)]
pub struct Incoming {
    pub name: String,
    pub bytes: Result<Vec<u8>, String>,
}

/// Files read by the browser in the background, waiting for the next frame to open them.
#[derive(Debug, Clone, Default)]
pub struct Inbox(Arc<Mutex<Vec<Incoming>>>);

impl Inbox {
    pub fn push(&self, file: Incoming) {
        self.0
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .push(file);
    }

    pub fn take(&self) -> Vec<Incoming> {
        std::mem::take(&mut *self.0.lock().unwrap_or_else(PoisonError::into_inner))
    }
}

/// Whether a file of `size` bytes may be loaded (the same limit as for files on disk).
pub fn size_allowed(size: f64) -> bool {
    size.is_finite() && size >= 0.0 && size <= MAX_FILE as f64
}

/// A file name the in-memory file system and a download can use: the last path component,
/// without control characters, characters Windows forbids, or leading dots; at most
/// [`MAX_NAME`] characters. Never empty.
pub fn safe_name(name: &str) -> String {
    let base = name.rsplit(['/', '\\']).next().unwrap_or_default();
    let cleaned: String = base
        .chars()
        .filter(|c| !c.is_control() && !matches!(c, '<' | '>' | ':' | '"' | '|' | '?' | '*'))
        .take(MAX_NAME)
        .collect();
    let cleaned = cleaned.trim().trim_start_matches('.').trim();
    if cleaned.is_empty() {
        "Untitled".to_string()
    } else {
        cleaned.to_string()
    }
}

/// `name` with a number added before its extension: `a.dna` → `a (2).dna`.
fn numbered(name: &str, n: u32) -> String {
    match name.rsplit_once('.') {
        Some((stem, ext)) if !stem.is_empty() => format!("{stem} ({n}).{ext}"),
        _ => format!("{name} ({n})"),
    }
}

/// Copy a file the browser handed over into the session's file system, under its own name or,
/// when that is taken (the same file opened twice, maybe changed in between), a numbered one,
/// so an open document is never silently replaced. Returns the path to open.
pub fn store(platform: &dyn Platform, name: &str, bytes: &[u8]) -> Result<String, String> {
    let name = safe_name(name);
    let path = std::iter::once(name.clone())
        .chain((2..1000).map(|n| numbered(&name, n)))
        .find(|p| !platform.exists(p))
        .ok_or_else(|| format!("too many files named {name} are open"))?;
    platform
        .write(&path, bytes)
        .map_err(|e| format!("{name} couldn't be loaded: {e}"))?;
    Ok(path)
}

/// What a Download prompt saves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SaveTarget {
    /// The active document, in the format its extension names.
    Document,
    /// A `.dna` copy of the active document without its construction history.
    DocumentWithoutHistory,
    /// Experiment tab `i`, as `.sslab`.
    Experiment(usize),
    /// The map of the active document, as SVG or PNG.
    Map,
}

impl SaveTarget {
    fn title(self) -> &'static str {
        match self {
            SaveTarget::Document => "Download sequence",
            SaveTarget::DocumentWithoutHistory => "Download copy without history",
            SaveTarget::Experiment(_) => "Download experiment",
            SaveTarget::Map => "Download map",
        }
    }

    fn hint(self) -> &'static str {
        match self {
            SaveTarget::Document => {
                "The extension picks the format: .dna, .gb, .embl, .gff3, .sbol or .fasta."
            }
            SaveTarget::DocumentWithoutHistory => "Saved as a SnapGene .dna file.",
            SaveTarget::Experiment(_) => "Saved as a SnappieSnap experiment (.sslab).",
            SaveTarget::Map => "The extension picks the format: .svg (vector) or .png.",
        }
    }

    /// The file name to save under: `name` made safe, with this target's default extension
    /// added when it doesn't name a format the target can write.
    pub fn file_name(self, name: &str) -> String {
        let name = safe_name(name);
        let lower = name.to_ascii_lowercase();
        let (ok, ext) = match self {
            SaveTarget::Document => (
                snappiesnap_engine::io::format_for_path(&name).is_some_and(|f| {
                    !matches!(
                        f,
                        snappiesnap_engine::io::Format::Abif | snappiesnap_engine::io::Format::Scf
                    )
                }),
                ".dna",
            ),
            SaveTarget::DocumentWithoutHistory => (lower.ends_with(".dna"), ".dna"),
            SaveTarget::Experiment(_) => (lower.ends_with(".sslab"), ".sslab"),
            SaveTarget::Map => (lower.ends_with(".svg") || lower.ends_with(".png"), ".svg"),
        };
        if ok { name } else { format!("{name}{ext}") }
    }
}

/// The web build's stand-in for a Save dialog: a file name, then a download.
#[derive(Debug, Clone)]
pub struct SavePrompt {
    pub target: SaveTarget,
    pub name: String,
    focused: bool,
}

impl SavePrompt {
    pub fn new(target: SaveTarget, name: &str) -> Self {
        SavePrompt {
            target,
            name: target.file_name(name),
            focused: false,
        }
    }
}

/// Show the prompt: `Some(true)` to download, `Some(false)` to cancel, `None` while open.
pub fn prompt(ctx: &egui::Context, t: &Tokens, p: &mut SavePrompt) -> Option<bool> {
    let mut answer = None;
    egui::Window::new(p.target.title())
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label("File name");
            let edit = ui.add(egui::TextEdit::singleline(&mut p.name).desired_width(320.0));
            if !p.focused {
                edit.request_focus();
                p.focused = true;
            }
            ui.label(RichText::new(p.target.hint()).small().color(t.text_muted));
            ui.label(
                RichText::new("Your browser saves it in its downloads folder.")
                    .small()
                    .color(t.text_muted),
            );
            ui.add_space(4.0);
            ui.horizontal(|ui| {
                let enter = edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
                if ui.button("Download").clicked() || enter {
                    answer = Some(true);
                }
                if ui.button("Cancel").clicked() || ui.input(|i| i.key_pressed(egui::Key::Escape)) {
                    answer = Some(false);
                }
            });
        });
    answer
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_engine::MemoryFs;

    #[test]
    fn names_are_made_safe() {
        assert_eq!(safe_name("pUC19.dna"), "pUC19.dna");
        assert_eq!(safe_name("C:\\lab\\run 1\\read.ab1"), "read.ab1");
        assert_eq!(safe_name("../../etc/passwd"), "passwd");
        assert_eq!(safe_name("a<b>:c\"d|e?f*.gb"), "abcdef.gb");
        assert_eq!(safe_name("bad\u{0}\nname.fa"), "badname.fa");
        assert_eq!(safe_name(".hidden.dna"), "hidden.dna");
        for empty in ["", "   ", ".", "..", "/", "dir/", "\u{7}"] {
            assert_eq!(safe_name(empty), "Untitled", "{empty:?}");
        }
        let long = format!("{}.dna", "é".repeat(500));
        assert_eq!(safe_name(&long).chars().count(), MAX_NAME);
    }

    #[test]
    fn stored_files_never_replace_open_ones() {
        let fs = MemoryFs::new();
        assert_eq!(store(&fs, "p.dna", b"one").unwrap(), "p.dna");
        assert_eq!(store(&fs, "p.dna", b"two").unwrap(), "p (2).dna");
        assert_eq!(store(&fs, "dir/p.dna", b"three").unwrap(), "p (3).dna");
        assert_eq!(store(&fs, "README", b"x").unwrap(), "README");
        assert_eq!(store(&fs, "README", b"y").unwrap(), "README (2)");
        assert_eq!(fs.get("p.dna").unwrap(), b"one");
        assert_eq!(fs.get("p (2).dna").unwrap(), b"two");
        assert_eq!(store(&fs, "", b"").unwrap(), "Untitled");
    }

    #[test]
    fn inbox_hands_files_over_once() {
        let inbox = Inbox::default();
        inbox.clone().push(Incoming {
            name: "a.fa".into(),
            bytes: Ok(b">a\nACGT\n".to_vec()),
        });
        inbox.push(Incoming {
            name: "b.fa".into(),
            bytes: Err("unreadable".into()),
        });
        let files = inbox.take();
        assert_eq!(files.len(), 2);
        assert!(inbox.take().is_empty());
    }

    #[test]
    fn size_limit() {
        assert!(size_allowed(0.0) && size_allowed(1e6) && size_allowed(MAX_FILE as f64));
        for bad in [-1.0, f64::NAN, f64::INFINITY, MAX_FILE as f64 + 1.0] {
            assert!(!size_allowed(bad), "{bad}");
        }
    }

    #[test]
    fn save_names_get_a_format() {
        use SaveTarget::*;
        assert_eq!(Document.file_name("pX"), "pX.dna");
        assert_eq!(Document.file_name("pX.gb"), "pX.gb");
        assert_eq!(Document.file_name("pX.FASTA"), "pX.FASTA");
        assert_eq!(Document.file_name("read.ab1"), "read.ab1.dna");
        assert_eq!(Document.file_name("v1.2"), "v1.2.dna");
        assert_eq!(DocumentWithoutHistory.file_name("pX.gb"), "pX.gb.dna");
        assert_eq!(DocumentWithoutHistory.file_name("pX.DNA"), "pX.DNA");
        assert_eq!(Experiment(0).file_name("run"), "run.sslab");
        assert_eq!(Map.file_name("pX map.png"), "pX map.png");
        assert_eq!(Map.file_name("pX map"), "pX map.svg");
        assert_eq!(Map.file_name(""), "Untitled.svg");
    }
}
