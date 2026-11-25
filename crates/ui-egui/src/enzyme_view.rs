//! The Enzymes view: which enzyme database is in use (with the REBASE download), every enzyme
//! (or those of one of the user's enzyme sets, such as the lab's stock) with how often it cuts
//! the open molecule, and the fragments of a digest with the enzymes the user ticks. It only
//! reads; changes go back to the shell as [`EnzymeAction`]s.

use std::collections::BTreeSet;
use std::sync::Arc;

use eframe::egui::{self, RichText, Sense, Ui, vec2};
use snappiesnap_engine::Session;
use snappiesnap_engine::bio::enzymes::{Fragment, digest, find_sites, sites_by_enzyme};
use snappiesnap_engine::doc::{Molecule, Range, Topology};
use snappiesnap_engine::enzymes::{EnzymeDb, REBASE_CREDIT};

use crate::theme::Tokens;

/// What the user asked for in the Enzymes view.
#[derive(Debug, Clone, PartialEq)]
pub enum EnzymeAction {
    /// Download REBASE (in the background) and use it.
    Download,
    /// Switch back to the built-in enzymes.
    UseBuiltin,
    /// Select a range of the molecule.
    Select(Range),
    /// Add a digest lane with these enzymes to the gel and show it.
    ShowOnGel(Vec<String>),
    /// Keep the ticked enzymes as an enzyme set of the user's own, under this name,
    /// replacing a set that already has it.
    SaveSet(String, Vec<String>),
    /// Import an enzyme set from a file (the browser's picker in the web build).
    ImportSet,
}

/// Which enzymes the table lists.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CutFilter {
    #[default]
    Cutters,
    Unique,
    NonCutters,
    All,
}

impl CutFilter {
    const ALL: [CutFilter; 4] = [
        CutFilter::Cutters,
        CutFilter::Unique,
        CutFilter::NonCutters,
        CutFilter::All,
    ];

    fn label(self) -> &'static str {
        match self {
            CutFilter::Cutters => "Cutters",
            CutFilter::Unique => "Unique",
            CutFilter::NonCutters => "Non-cutters",
            CutFilter::All => "All",
        }
    }

    fn keeps(self, cuts: usize) -> bool {
        match self {
            CutFilter::Cutters => cuts > 0,
            CutFilter::Unique => cuts == 1,
            CutFilter::NonCutters => cuts == 0,
            CutFilter::All => true,
        }
    }
}

/// Which enzymes the table lists from (remembered in the settings).
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ListFrom {
    /// The whole database.
    All,
    /// One of the user's enzyme sets, by name.
    Set(String),
}

/// The enzyme set the table lists, if any: the user's pick, else the preferred set (the one
/// the map uses), else the first set. Picking "All" shows the whole database. A set that no
/// longer exists falls back the same way.
fn listed_set<'a>(
    db: &'a EnzymeDb,
    from: Option<&ListFrom>,
    preferred: Option<&str>,
) -> Option<&'a snappiesnap_engine::enzymes::UserSet> {
    match from {
        Some(ListFrom::All) => None,
        Some(ListFrom::Set(name)) if db.set(name).is_some() => db.set(name),
        _ => preferred.and_then(|p| db.set(p)).or(db.sets.first()),
    }
}

/// Per enzyme (database order): the site ranges and top-strand cuts in the molecule.
type Cuts = Vec<Vec<(Range, Option<usize>)>>;

pub struct EnzymeViewState {
    /// What the table lists; `None` until the user picks (see [`listed_set`]). Remembered
    /// between sessions.
    pub from: Option<ListFrom>,
    pub filter: String,
    pub commercial_only: bool,
    pub cuts: CutFilter,
    /// Enzymes ticked for the digest, by name.
    pub chosen: BTreeSet<String>,
    /// The name the ticked enzymes would be saved under.
    pub set_name: String,
    /// A download is running.
    pub downloading: bool,
    cache: Option<(u64, Arc<EnzymeDb>, Arc<Cuts>)>,
    digest_cache: Option<(u64, Vec<String>, Vec<Fragment>)>,
}

impl Default for EnzymeViewState {
    fn default() -> Self {
        EnzymeViewState {
            from: None,
            filter: String::new(),
            commercial_only: true,
            cuts: CutFilter::default(),
            chosen: BTreeSet::new(),
            set_name: String::new(),
            downloading: false,
            cache: None,
            digest_cache: None,
        }
    }
}

impl EnzymeViewState {
    fn cuts(&mut self, session: &Session, m: &Molecule) -> Arc<Cuts> {
        let db = session.enzymes_arc();
        if let Some((rev, cached_db, cuts)) = &self.cache
            && *rev == session.revision()
            && Arc::ptr_eq(cached_db, &db)
        {
            return Arc::clone(cuts);
        }
        let mut cuts: Cuts = vec![Vec::new(); db.enzymes.len()];
        for (i, sites) in sites_by_enzyme(find_sites(m, &db.enzymes)) {
            if let Some(slot) = cuts.get_mut(i) {
                *slot = sites.into_iter().map(|s| (s.range, s.cut)).collect();
                slot.sort_by_key(|(range, cut)| (cut.is_none(), *cut, range.start));
            }
        }
        let cuts = Arc::new(cuts);
        self.cache = Some((session.revision(), db, Arc::clone(&cuts)));
        cuts
    }

    fn fragments(&mut self, session: &Session, m: &Molecule) -> Vec<Fragment> {
        let names: Vec<String> = self.chosen.iter().cloned().collect();
        if let Some((rev, cached_names, frags)) = &self.digest_cache
            && *rev == session.revision()
            && *cached_names == names
        {
            return frags.clone();
        }
        let db = session.enzymes();
        let enzymes: Vec<_> = names
            .iter()
            .filter_map(|n| db.get(n).map(|(_, e)| e.clone()))
            .collect();
        let frags = if enzymes.is_empty() {
            Vec::new()
        } else {
            digest(m, &enzymes)
        };
        self.digest_cache = Some((session.revision(), names, frags.clone()));
        frags
    }
}

/// Most rows shown in the table (the filter narrows the rest).
const MAX_ROWS: usize = 2000;

/// Draw the view. `preferred` is the enzyme set the map draws from, which the table lists by
/// default.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    state: &mut EnzymeViewState,
    preferred: Option<&str>,
) -> Option<EnzymeAction> {
    let mut action = None;
    let doc = session.document()?;
    let m = &doc.molecule;
    let db = session.enzymes();

    // The database in use, and the download.
    ui.horizontal(|ui| {
        let what = match &db.version {
            Some(v) => format!("{} {v}", db.source),
            None => db.source.clone(),
        };
        let commercial = db.enzymes.iter().filter(|e| e.commercial).count();
        // The built-in set is every enzyme a supplier sells, so counting them twice says
        // nothing; REBASE holds research enzymes nobody sells, where the count does.
        let sold = if commercial == db.enzymes.len() {
            ", all sold by a supplier".to_string()
        } else {
            format!(" ({commercial} commercially available)")
        };
        ui.label(
            RichText::new(format!("{what}: {} enzymes{sold}", db.enzymes.len()))
            .strong(),
        );
        if state.downloading {
            ui.spinner();
            ui.label(RichText::new("Downloading REBASE…").color(t.text_muted));
        } else {
            let label = if db.is_builtin() {
                "Download REBASE"
            } else {
                "Update REBASE"
            };
            let available = cfg!(not(target_arch = "wasm32"));
            let r = ui
                .add_enabled(available, egui::Button::new(label))
                .on_hover_text("Download the full enzyme database from rebase.neb.com (about 4 MB) and keep it for later sessions")
                .on_disabled_hover_text("Not available in the browser yet");
            if r.clicked() {
                action = Some(EnzymeAction::Download);
            }
        }
        if !db.is_builtin() && ui.button("Use built-in").clicked() {
            action = Some(EnzymeAction::UseBuiltin);
        }
    });
    if !db.is_builtin() {
        ui.label(RichText::new(REBASE_CREDIT).small().color(t.text_muted));
    }
    ui.add_space(t.spacing);

    // Filters: first which list (when the user has enzyme sets).
    let set = listed_set(db, state.from.as_ref(), preferred);
    ui.horizontal(|ui| {
        if !db.sets.is_empty() {
            ui.label("List");
            let shown = set.map_or("All enzymes".to_string(), |s| s.name.clone());
            egui::ComboBox::from_id_salt("enzyme-list")
                .selected_text(shown)
                .show_ui(ui, |ui| {
                    for user in &db.sets {
                        let on = set.is_some_and(|s| s.name == user.name);
                        if ui.selectable_label(on, &user.name).clicked() {
                            state.from = Some(ListFrom::Set(user.name.clone()));
                        }
                    }
                    if ui.selectable_label(set.is_none(), "All enzymes").clicked() {
                        state.from = Some(ListFrom::All);
                    }
                })
                .response
                .on_hover_text(
                    "List the enzymes of one of your enzyme sets, or every enzyme in the database",
                );
        }
        ui.add(
            egui::TextEdit::singleline(&mut state.filter)
                .hint_text("Filter by name or site")
                .desired_width(200.0),
        );
        for f in CutFilter::ALL {
            ui.selectable_value(&mut state.cuts, f, f.label());
        }
        // A set says what is at hand, so whether an enzyme is sold doesn't matter then.
        let mut commercial_only = state.commercial_only && set.is_none();
        let r = ui.add_enabled(
            set.is_none(),
            egui::Checkbox::new(&mut commercial_only, "Commercial only"),
        );
        if r.changed() {
            state.commercial_only = commercial_only;
        }
    });
    // The set's enzymes that the database lacks, and how to get them.
    if let Some(set) = set {
        let missing = set.missing(db);
        let known = set.enzymes.len().saturating_sub(missing.len());
        let mut text = format!(
            "{}: {known} of {} enzymes in the database",
            set.name,
            set.enzymes.len()
        );
        if !missing.is_empty() {
            let fix = if db.is_builtin() {
                "; download REBASE to include the rest"
            } else {
                ""
            };
            text.push_str(&format!(" (not found: {}{fix})", missing.join(", ")));
        }
        let color = if missing.is_empty() {
            t.text_muted
        } else {
            t.error
        };
        ui.add(egui::Label::new(RichText::new(text).color(color)).wrap());
    }
    ui.add_space(t.spacing);

    let cuts = state.cuts(session, m);
    let needle = state.filter.trim().to_ascii_uppercase();
    let members: Option<BTreeSet<String>> =
        set.map(|s| s.enzymes.iter().map(|n| n.to_ascii_uppercase()).collect());
    let rows: Vec<usize> = db
        .enzymes
        .iter()
        .enumerate()
        .filter(|(_, e)| match &members {
            Some(m) => m.contains(&e.name.to_ascii_uppercase()),
            None => e.commercial || !state.commercial_only,
        })
        .filter(|(_, e)| {
            needle.is_empty()
                || e.name.to_ascii_uppercase().contains(&needle)
                || e.site.contains(&needle)
        })
        .filter(|(i, _)| state.cuts.keeps(cuts.get(*i).map_or(0, Vec::len)))
        .map(|(i, _)| i)
        .take(MAX_ROWS)
        .collect();

    let digest_width = 300.0;
    let table_width = (ui.available_width() - digest_width - t.spacing * 4.0).max(300.0);
    ui.horizontal_top(|ui| {
        // The enzyme table.
        ui.vertical(|ui| {
            ui.set_width(table_width);
            let row_h = ui.text_style_height(&egui::TextStyle::Body) + 4.0;
            let header = |ui: &mut Ui, text: &str, w: f32| {
                cell(ui, w, row_h, egui::Label::new(RichText::new(text).strong()));
            };
            ui.horizontal(|ui| {
                header(ui, "", 20.0);
                header(ui, "Enzyme", 110.0);
                header(ui, "Site", 190.0);
                header(ui, "Cuts", 40.0);
                ui.label(RichText::new("Cut after").strong());
            });
            // Ticking the whole list is how an enzyme set is loaded for editing (List: the
            // set, Tick all, change what is ticked, save it under the same name) and how a
            // set is made from a filtered table.
            ui.horizontal(|ui| {
                let n = rows.len();
                if ui
                    .small_button(format!("Tick all {n}"))
                    .on_hover_text("Tick every enzyme the table is showing")
                    .clicked()
                {
                    for &i in &rows {
                        if let Some(e) = db.enzymes.get(i) {
                            state.chosen.insert(e.name.clone());
                        }
                    }
                }
                if ui
                    .small_button("Untick all")
                    .on_hover_text("Untick every enzyme the table is showing")
                    .clicked()
                {
                    for &i in &rows {
                        if let Some(e) = db.enzymes.get(i) {
                            state.chosen.remove(&e.name);
                        }
                    }
                }
                #[cfg(target_arch = "wasm32")]
                if ui
                    .small_button("Import a set\u{2026}")
                    .on_hover_text("Read an enzyme set from a text file of enzyme names")
                    .clicked()
                {
                    action = Some(EnzymeAction::ImportSet);
                }
            });
            egui::ScrollArea::vertical()
                .id_salt("enzyme-table")
                .auto_shrink([false, false])
                .show_rows(ui, row_h, rows.len(), |ui, visible| {
                    for &i in rows.get(visible).unwrap_or_default() {
                        let Some(e) = db.enzymes.get(i) else { continue };
                        let sites = cuts.get(i).map(Vec::as_slice).unwrap_or_default();
                        ui.horizontal(|ui| {
                            let mut on = state.chosen.contains(&e.name);
                            if ui
                                .add_sized([20.0, row_h], egui::Checkbox::without_text(&mut on))
                                .changed()
                            {
                                if on {
                                    state.chosen.insert(e.name.clone());
                                } else {
                                    state.chosen.remove(&e.name);
                                }
                            }
                            let name = RichText::new(&e.name).color(if e.commercial {
                                t.text
                            } else {
                                t.text_muted
                            });
                            cell(ui, 110.0, row_h, egui::Label::new(name).truncate());
                            cell(
                                ui,
                                190.0,
                                row_h,
                                egui::Label::new(RichText::new(e.notation()).monospace())
                                    .truncate(),
                            );
                            ui.add_sized([40.0, row_h], egui::Label::new(sites.len().to_string()));
                            for (range, cut) in sites.iter().take(8) {
                                let text =
                                    cut.map_or("?".to_string(), |c| shown_cut(m, c).to_string());
                                let r = ui.add(egui::Label::new(text).sense(Sense::click()));
                                if r.on_hover_text("Select this site").clicked() {
                                    action = Some(EnzymeAction::Select(*range));
                                }
                            }
                            if sites.len() > 8 {
                                ui.label(RichText::new("…").color(t.text_muted));
                            }
                        });
                    }
                });
        });

        ui.separator();

        // The digest.
        ui.vertical(|ui| {
            ui.set_width(digest_width);
            ui.label(RichText::new("Digest").strong());
            if state.chosen.is_empty() {
                ui.label(
                    RichText::new("Tick enzymes in the table to digest the molecule with them.")
                        .color(t.text_muted),
                );
                return;
            }
            ui.horizontal_wrapped(|ui| {
                ui.label(state.chosen.iter().cloned().collect::<Vec<_>>().join(", "));
                if ui.small_button("Clear").clicked() {
                    state.chosen.clear();
                }
                if ui.small_button("Show on gel").clicked() {
                    action = Some(EnzymeAction::ShowOnGel(
                        state.chosen.iter().cloned().collect(),
                    ));
                }
            });
            // The ticked enzymes as a set of the user's own: what this lab can run. Saving
            // under the name of a set that exists replaces it, which is how a set is edited:
            // list it, "Tick all", change the ticks, save it again.
            ui.horizontal_wrapped(|ui| {
                ui.add(
                    egui::TextEdit::singleline(&mut state.set_name)
                        .hint_text("Name for an enzyme set")
                        .desired_width(160.0),
                );
                let named = !state.set_name.trim().is_empty();
                let exists = db
                    .sets
                    .iter()
                    .any(|x| x.name.eq_ignore_ascii_case(state.set_name.trim()));
                let label = if exists { "Replace set" } else { "Save as set" };
                let r = ui
                    .add_enabled(named, egui::Button::new(label))
                    .on_disabled_hover_text("Give the set a name first");
                if r.clicked() {
                    action = Some(EnzymeAction::SaveSet(
                        state.set_name.trim().to_string(),
                        state.chosen.iter().cloned().collect(),
                    ));
                }
                if exists {
                    ui.label(
                        RichText::new(format!("replaces \"{}\"", state.set_name.trim()))
                            .color(t.warning),
                    );
                }
            });
            if cfg!(target_arch = "wasm32") {
                ui.label(
                    RichText::new(
                        "In the browser a set is kept for this visit only; the page can't \
                         write files that outlive it.",
                    )
                    .color(t.text_muted),
                );
            }
            let frags = state.fragments(session, m);
            ui.label(
                RichText::new(format!(
                    "{} fragment{}{}",
                    frags.len(),
                    if frags.len() == 1 { "" } else { "s" },
                    if m.topology() == Topology::Circular
                        && frags.len() == 1
                        && frags.first().is_some_and(|f| f.left.is_none())
                    {
                        " (uncut circle)"
                    } else {
                        ""
                    }
                ))
                .color(t.text_muted),
            );
            ui.add_space(t.spacing);
            egui::ScrollArea::vertical()
                .id_salt("digest")
                .show(ui, |ui| {
                    egui::Grid::new("fragments")
                        .striped(true)
                        .spacing(vec2(t.spacing * 2.0, t.spacing * 0.5))
                        .show(ui, |ui| {
                            for h in ["Length", "From", "To", "Ends"] {
                                ui.label(RichText::new(h).strong());
                            }
                            ui.end_row();
                            for f in &frags {
                                let r = ui.add(
                                    egui::Label::new(format!("{} bp", f.length))
                                        .sense(Sense::click()),
                                );
                                if r.on_hover_text("Select this fragment").clicked()
                                    && !f.range.is_empty()
                                {
                                    action = Some(EnzymeAction::Select(f.range));
                                }
                                ui.label((f.range.start + 1).to_string());
                                ui.label(
                                    if f.range.end == 0 {
                                        m.len()
                                    } else {
                                        f.range.end
                                    }
                                    .to_string(),
                                );
                                ui.label(format!(
                                    "{} – {}",
                                    f.left.as_deref().unwrap_or("end"),
                                    f.right.as_deref().unwrap_or("end")
                                ));
                                ui.end_row();
                            }
                        });
                });
        });
    });
    action
}

/// A left-aligned table cell of fixed width.
fn cell(ui: &mut Ui, width: f32, height: f32, widget: impl egui::Widget) -> egui::Response {
    ui.allocate_ui_with_layout(
        vec2(width, height),
        egui::Layout::left_to_right(egui::Align::Center),
        |ui| {
            ui.set_width(width);
            ui.add(widget)
        },
    )
    .inner
}

/// A cut as shown to the user: residues before it, with a cut at the origin of a circular
/// molecule shown after its last residue.
fn shown_cut(m: &Molecule, cut: usize) -> usize {
    if cut == 0 && m.topology() == Topology::Circular {
        m.len()
    } else {
        cut
    }
}

/// Download REBASE on a background thread and save it as the cache file. The receiver gets
/// the cache path to load, or an error message.
#[cfg(not(target_arch = "wasm32"))]
pub fn start_download(session: &Session) -> std::sync::mpsc::Receiver<Result<String, String>> {
    use snappiesnap_engine::enzymes::REBASE_URL;
    let (tx, rx) = std::sync::mpsc::channel();
    let platform = session.platform_arc();
    let cache = session.enzyme_cache_path();
    std::thread::spawn(move || {
        let result = (|| {
            let path = cache.ok_or("this platform has no folder for app data")?;
            let bytes = platform
                .fetch(REBASE_URL)
                .map_err(|e| format!("download failed: {e}"))?;
            // Check before replacing a working cache.
            EnzymeDb::from_rebase(&bytes)?;
            platform
                .write(&path, &bytes)
                .map_err(|e| format!("cannot save {path}: {e}"))?;
            Ok(path)
        })();
        let _ = tx.send(result);
    });
    rx
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::Document;
    use snappiesnap_engine::doc::{Alphabet, Seq};

    use super::*;

    #[test]
    fn draws_table_and_digest() {
        let s = format!(
            "{}GAATTC{}GGATCC{}",
            "A".repeat(10),
            "A".repeat(20),
            "A".repeat(30)
        );
        let m =
            Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Circular).unwrap();
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        session.open_document(Document::new(m));
        let mut state = EnzymeViewState::default();
        state.chosen.insert("EcoRI".into());
        state.chosen.insert("BamHI".into());
        state.chosen.insert("NoSuchI".into());
        let ctx = egui::Context::default();
        for _ in 0..2 {
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                let t = Tokens::dark();
                show(ui, &t, &session, &mut state, None);
            });
            out.textures_delta.clear();
        }
        let frags = state.fragments(&session, &session.document().unwrap().molecule);
        let lens: Vec<usize> = frags.iter().map(|f| f.length).collect();
        assert_eq!(lens, vec![46, 26]);
        let cuts = state.cuts(&session, &session.document().unwrap().molecule);
        let eco = session.enzymes().get("EcoRI").unwrap().0;
        assert_eq!(cuts[eco].len(), 1);
    }

    /// Ticking enzymes and naming them makes a set, with no file and no dialog anywhere, so
    /// it works in the browser too. Saving under the name of a set that exists replaces it,
    /// which is how a set is edited.
    #[test]
    fn ticked_enzymes_become_an_enzyme_set() {
        use snappiesnap_engine::enzymes::UserSet;
        let s = format!(
            "{}GAATTC{}GGATCC{}",
            "A".repeat(10),
            "A".repeat(20),
            "A".repeat(30)
        );
        let m =
            Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Circular).unwrap();
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        session.open_document(Document::new(m));
        session
            .set_enzyme_sets(vec![UserSet {
                name: "Fridge".into(),
                enzymes: vec!["EcoRI".into()],
            }])
            .expect("sets kept");
        let mut state = EnzymeViewState::default();
        state.chosen.insert("EcoRI".into());
        state.chosen.insert("BamHI".into());

        // The button is dead until the set has a name, and no action escapes by itself.
        let run = |state: &mut EnzymeViewState, session: &Session| {
            let ctx = egui::Context::default();
            let mut seen = None;
            for _ in 0..3 {
                let mut out = ctx.run_ui(egui::RawInput::default(), |ui| {
                    let t = Tokens::dark();
                    let a = show(ui, &t, session, state, None);
                    // Three passes run the closure three times: keep the first answer.
                    if seen.is_none() {
                        seen = a;
                    }
                });
                out.textures_delta.clear();
            }
            seen
        };
        assert_eq!(run(&mut state, &session), None);

        // What the button sends when it is pressed: the name, and the ticked enzymes.
        state.set_name = "  My set  ".into();
        let action = EnzymeAction::SaveSet(
            state.set_name.trim().to_string(),
            state.chosen.iter().cloned().collect(),
        );
        assert_eq!(
            action,
            EnzymeAction::SaveSet("My set".into(), vec!["BamHI".into(), "EcoRI".into()])
        );
        // The command behind it is the one an agent would call, and it replaces by name.
        let v = session
            .execute(
                "enzyme_set_create",
                serde_json::json!({"name": "Fridge", "enzymes": ["EcoRI", "BamHI"]}),
            )
            .expect("set saved");
        assert_eq!(v["replaced"], serde_json::json!(true));
        assert_eq!(v["count"], serde_json::json!(2));
        assert_eq!(session.enzymes().sets.len(), 1);
    }

    /// With enzyme sets, the table lists the preferred set (else the first) until the user
    /// picks another set or the whole database. Set
    /// `SNAPPIESNAP_ENZYME_VIEW_PNG=/path.png` to look at it.
    #[test]
    fn lists_an_enzyme_set_by_default() {
        use snappiesnap_engine::enzymes::UserSet;
        let mut db = EnzymeDb::default();
        assert!(listed_set(&db, None, Some("Lab")).is_none());
        let set = |name: &str| UserSet {
            name: name.into(),
            enzymes: vec!["EcoRI".into(), "XcmI".into()],
        };
        db.sets = vec![set("First"), set("Lab")];
        let name = |from: Option<&ListFrom>, preferred| {
            listed_set(&db, from, preferred).map(|s| s.name.clone())
        };
        assert_eq!(name(None, None).as_deref(), Some("First"));
        assert_eq!(name(None, Some("lab")).as_deref(), Some("Lab"));
        assert_eq!(name(None, Some("gone")).as_deref(), Some("First"));
        assert_eq!(name(Some(&ListFrom::All), Some("Lab")), None);
        let gone = ListFrom::Set("gone".into());
        assert_eq!(name(Some(&gone), Some("Lab")).as_deref(), Some("Lab"));

        let s = format!("{0}GAATTC{0}GGATCC{0}", "A".repeat(20));
        let m = Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Linear).unwrap();
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        session
            .execute(
                "enzyme_set_create",
                serde_json::json!({"name": "Lab", "enzymes": ["EcoRI", "SmaI", "XcmI"]}),
            )
            .unwrap();
        session.open_document(Document::new(m));
        let mut state = EnzymeViewState::default();
        // With enzymes ticked and a name typed, the set controls are drawn too.
        state.chosen.insert("EcoRI".into());
        state.chosen.insert("SmaI".into());
        state.set_name = "Lab".into();
        let t = Tokens::dark();
        let image = crate::offscreen::render(
            900,
            420,
            1.0,
            |ctx| t.apply(ctx),
            |ui| {
                show(ui, &t, &session, &mut state, None);
            },
        )
        .unwrap();
        if let Ok(path) = std::env::var("SNAPPIESNAP_ENZYME_VIEW_PNG") {
            crate::offscreen::write_png(&path, &image).unwrap();
        }
        assert_eq!(state.from, None);
    }
}
