//! The Cloning experiment: a construct built in four steps, in a tab of its own, bound to no
//! open document.
//!
//! 1. the technique (every one the engine simulates);
//! 2. the vector, from the open folder or from disk, with its map and how it is prepared;
//! 3. the insert (or, for mutagenesis, the primers), the same way, and whether it fits the
//!    vector: the reaction itself is run as a preview, so the answer is the engine's, not a
//!    guess of the UI's;
//! 4. the product, with its map, and a button that opens it as a new construct.
//!
//! Every step builds the params of the same command an agent would call
//! (`restriction_clone`, `gateway_clone`, `assemble`, `topo_clone`, `mutagenesis`), so what
//! the window does and what the control channel does cannot drift apart.

use std::collections::HashMap;

use eframe::egui::{self, RichText, Ui};
use serde_json::{Value, json};
use snappiesnap_engine::Session;
use snappiesnap_engine::doc::{Molecule, Range, Topology};

use crate::state::{CloningSide, CloningTechnique, Experiment, ExperimentStep, SideUse};
use crate::theme::Tokens;

/// How deep into the open folder constructs are looked for, and how many are listed.
const MAX_DEPTH: usize = 4;
/// Most enzymes offered for one construct, and most sites an offered one may have: past a
/// few cuts the fragments are no longer what you clone with.
const MAX_OFFERED: usize = 200;
const MAX_OFFERED_CUTS: usize = 3;
/// Fewest specific bases in an offered enzyme's site.
const MIN_OFFERED_LENGTH: usize = 6;
/// Most constructs whose sites are kept at once.
const MAX_CACHED: usize = 32;
/// Most enzymes labelled on the construct when none have been typed.
const MAX_MAP_LABELS: usize = 12;
/// Least width the construct's map gets, and the room kept for the plan under it.
const MIN_MAP_W: f32 = 320.0;
/// Room kept under the map for the plan: its caption, three rows and the legend.
const PLAN_H: f32 = 210.0;
/// Room kept at the bottom for the Back and Continue buttons.
const NAV_H: f32 = 48.0;
/// Width of the layer toggles down the left edge of the construct, as in a document tab.
const LAYER_W: f32 = 38.0;
const MAX_FILES: usize = 500;
const MAX_FOLDERS: usize = 200;
/// How many parts an assembly may have in all (the vector, the second part and the rest).
const MAX_PARTS: usize = 10;

/// What a chosen file holds, as `file_info` reports it.
#[derive(Debug, Clone, Default)]
pub struct ConstructInfo {
    pub name: String,
    pub length: u64,
    pub topology: String,
    pub features: u64,
    pub molecule: Option<Molecule>,
}

/// The width the schematic is drawn at where there is a whole window to spare.
const PLAN_W: f32 = 620.0;

/// A piece a digest makes: how long it is and what its two ends are. The ends say whether a
/// piece fits the backbone, and whether its direction is settled.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Frag {
    pub len: u64,
    pub left: snappiesnap_engine::bio::cloning::End,
    pub right: snappiesnap_engine::bio::cloning::End,
}

/// The product of the reaction, as the command's preview reports it.
#[derive(Debug, Clone, Default)]
pub struct Product {
    pub length: u64,
    pub description: String,
    pub molecule: Option<Molecule>,
}

/// An enzyme a construct can be cut with: its name and how many times it cuts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cutter {
    pub name: String,
    pub cuts: usize,
}

/// What this construct can be cut with, for the enzyme fields of one side: the lab's set
/// when there is one, and the enzymes of it that cut this construct.
pub struct Offer<'a> {
    pub pool: Option<&'a str>,
    pub db: &'a snappiesnap_engine::EnzymeDb,
    pub cutters: &'a [Cutter],
}

/// A construct the experiment can be built from: a tab that is already open, or a file of
/// the open folder.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Choice {
    /// An open tab, by its stable id, and its name. In the web build, where the browser gives
    /// the page no file paths, this is the only way to choose a construct.
    Document {
        id: u64,
        name: String,
        path: Option<String>,
    },
    /// A file of the open folder or picked from disk.
    File(String),
}

/// How the construct is shown beside the choices: the same view and the same layers a
/// document tab offers, so what is learned there is true here too. Held by the app (it is UI
/// state, remembered between sessions) and passed in and out of the step.
#[derive(Debug, Clone, Copy)]
pub struct Pane {
    /// The sequence rather than the map.
    pub sequence: bool,
    pub layers: crate::map::Layers,
    /// The rest of what a document tab draws with, so the construct looks the same here.
    pub style: crate::state::MapStyle,
    pub translation: crate::state::Translation,
    pub genetic_code: u8,
    pub tails: crate::state::TailStyle,
    /// The sequence on one line, scrolled sideways (the document views' own setting).
    pub single_line: bool,
}

/// The cloning step's own working state: what it has read and previewed. None of it belongs
/// to the experiment itself, so it is never saved.
pub struct CloningInputs {
    /// The enzyme set the app is drawing from (what the lab holds), set each frame by the
    /// experiment view. When there is one, the enzyme fields offer its enzymes only.
    pub pool: Option<String>,
    /// Which enzymes cut a construct, per construct and set.
    cutters: HashMap<String, Vec<Cutter>>,
    /// The sites drawn on a construct's map, per construct, set and typed enzymes.
    sites: HashMap<String, (String, Vec<crate::enzymes::MapSite>)>,
    /// Its ORFs and primer binding sites, per construct.
    #[allow(clippy::type_complexity)]
    extras: HashMap<
        String,
        (
            Vec<crate::enzymes::OrfMark>,
            Vec<crate::enzymes::PrimerMark>,
        ),
    >,
    /// What each construct holds, or why it couldn't be read, by [`CloningSide::source_key`].
    info: HashMap<String, Result<ConstructInfo, String>>,
    /// The session revision an open tab's details were taken at, so a tab edited while the
    /// experiment is open is seen as it now stands.
    info_rev: HashMap<String, u64>,
    /// The constructs of the open folder (and its subfolders), and the folder they are from.
    files: Vec<String>,
    listed: Option<String>,
    /// The construct's sequence view (scroll position and its cached translations).
    seq: crate::seqview::SeqViewState,
    /// What is selected in it, and in which construct: selecting and copying work here even
    /// though editing doesn't, and a selection means nothing once another construct is shown.
    selection: Option<(String, Range)>,
    /// An edit was attempted in the sequence view, which this pane doesn't apply.
    read_only_note: bool,
    /// What the last copy did, shown under the sequence.
    copy_note: Option<String>,
    /// An enzyme double-clicked in the sequence view, waiting to join the cut list.
    chosen_in_sequence: Option<usize>,
    /// What cutting a construct gives, per construct and enzyme list. The schematic draws
    /// the real pieces rather than a sketch of them.
    fragments: HashMap<String, (String, Vec<Frag>)>,
    /// The pane's own map zoom and pan, as a document tab keeps its own.
    zoom: f32,
    pan: egui::Vec2,
    /// The preview: the call it was made for, and the answer.
    preview_key: Option<String>,
    preview: Option<Result<Product, String>>,
}

impl Default for CloningInputs {
    fn default() -> Self {
        CloningInputs {
            // A map starts whole, as a document tab's does; the rest is empty.
            chosen_in_sequence: None,
            fragments: HashMap::new(),
            zoom: 1.0,
            pan: egui::Vec2::ZERO,
            pool: None,
            cutters: HashMap::new(),
            sites: HashMap::new(),
            extras: HashMap::new(),
            info: HashMap::new(),
            info_rev: HashMap::new(),
            files: Vec::new(),
            listed: None,
            seq: Default::default(),
            selection: None,
            read_only_note: false,
            copy_note: None,
            preview_key: None,
            preview: None,
        }
    }
}

impl CloningInputs {
    /// Keep what `file_info` answered for `path`.
    pub fn set_info(&mut self, path: String, answer: Result<Value, String>) {
        self.info.insert(
            path,
            answer.map(|v| {
                let num = |k: &str| v.get(k).and_then(Value::as_u64).unwrap_or(0);
                ConstructInfo {
                    name: v
                        .get("name")
                        .and_then(Value::as_str)
                        .unwrap_or_default()
                        .to_string(),
                    length: num("length"),
                    topology: v
                        .get("topology")
                        .and_then(Value::as_str)
                        .unwrap_or_default()
                        .to_string(),
                    features: num("features"),
                    molecule: v
                        .get("molecule")
                        .and_then(|m| serde_json::from_value(m.clone()).ok()),
                }
            }),
        );
    }

    /// Keep what the reaction's preview answered.
    pub fn set_preview(&mut self, answer: Result<Value, String>) {
        self.preview = Some(answer.map(product_of));
    }

    /// What is selected in the construct `key` is shown, and a selection made in another
    /// construct is not: they are different molecules, and a range means nothing across them.
    fn selection_in(&self, key: &str) -> Option<Range> {
        self.selection
            .as_ref()
            .filter(|(had, _)| had == key)
            .map(|(_, r)| *r)
    }

    fn select_in(&mut self, key: &str, range: Option<Range>) {
        self.selection = range.map(|r| (key.to_string(), r));
    }

    /// Forget the preview, so the next frame asks again (after an input changed).
    fn forget_preview(&mut self) {
        self.preview = None;
        self.preview_key = None;
    }
}

/// Read a preview answer, whichever command gave it.
fn product_of(v: Value) -> Product {
    let num = |value: &Value, k: &str| value.get(k).and_then(Value::as_u64).unwrap_or(0);
    // Restriction and TOPO list their constructs; the rest report one.
    let first = v
        .get("constructs")
        .and_then(Value::as_array)
        .and_then(|l| l.first())
        .unwrap_or(&v);
    let mut description = first
        .get("description")
        .and_then(Value::as_str)
        .unwrap_or_default()
        .to_string();
    if description.is_empty()
        && let Some(sites) = v.get("sites").and_then(Value::as_array)
    {
        let names: Vec<&str> = sites
            .iter()
            .filter_map(|s| s.get("name")?.as_str())
            .collect();
        description = format!("recombination sites: {}", names.join(" and "));
    }
    if description.is_empty()
        && let Some(junctions) = v.get("junctions").and_then(Value::as_array)
    {
        description = format!("{} junction(s)", junctions.len());
    }
    if description.is_empty() {
        let change = v.get("length_change").and_then(Value::as_i64).unwrap_or(0);
        description = match change {
            0 => format!("{} substitution(s)", num(&v, "substitutions")),
            c if c > 0 => format!("{c} bp inserted"),
            c => format!("{} bp deleted", -c),
        };
    }
    let others = v
        .get("constructs")
        .and_then(Value::as_array)
        .map_or(0, Vec::len);
    if others > 1 {
        description.push_str(&format!(" ({others} constructs; the first is shown)"));
    }
    Product {
        length: num(first, "length"),
        description,
        molecule: first
            .get("molecule")
            .and_then(|m| serde_json::from_value(m.clone()).ok()),
    }
}

/// What the cloning step wants done, handed back to the app.
pub enum CloningAction {
    /// Read a file with `file_info`; the answer goes to [`CloningInputs::set_info`].
    Info(String),
    /// Run the reaction as a preview; the answer goes to [`CloningInputs::set_preview`].
    Preview(&'static str, Value),
    /// Run the reaction for real, opening the product as a new construct.
    Open(&'static str, Value),
    /// Open the sample vector and insert, and set the experiment up to clone them.
    Samples,
}

/// The names in an enzyme or primer field, as typed.
fn names(text: &str) -> Vec<String> {
    text.split([',', ' ', '+'])
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect()
}

/// Add `name` to a field of names, keeping what is already typed there.
fn add_name(text: &mut String, name: &str) {
    if names(text).iter().any(|n| n.eq_ignore_ascii_case(name)) {
        return;
    }
    let kept = text.trim_end_matches(|c: char| c.is_whitespace() || c == ',' || c == '+');
    *text = if kept.is_empty() {
        name.to_string()
    } else {
        format!("{kept}, {name}")
    };
}

/// The constructs the picker offers: the open tabs first (what the user already has in front
/// of them, and the only choice the web build has), then the files of the open folder.
fn choices(session: &Session, files: &[String]) -> Vec<Choice> {
    let mut out: Vec<Choice> = session
        .documents()
        .into_iter()
        .map(|d| Choice::Document {
            id: d.id,
            name: d.name,
            path: d.path,
        })
        .collect();
    out.extend(files.iter().cloned().map(Choice::File));
    out
}

/// What to call a side in a sentence: the name of the open tab it is bound to, the file's
/// name, or (nothing chosen yet) an empty string. A construct picked from the open tabs has
/// no file name to fall back on, which is every construct in the web build.
pub fn side_label(session: &Session, side: &CloningSide) -> String {
    match side.document {
        Some(id) => session
            .document_by_id(id)
            .map(|d| d.molecule.name().to_string())
            .unwrap_or_else(|| "a tab that was closed".to_string()),
        None => side.file_name().to_string(),
    }
}

/// The name of the open tab a side is bound to, when it is still open.
fn chosen_name<'a>(side: &CloningSide, choices: &'a [Choice]) -> Option<&'a str> {
    let id = side.document?;
    choices.iter().find_map(|c| match c {
        Choice::Document {
            id: other, name, ..
        } if *other == id => Some(name.as_str()),
        _ => None,
    })
}

/// Keep an open tab's details in step with the tab itself: they are read again whenever the
/// session has moved on and that tab's molecule really changed, so editing a construct in its
/// own tab is seen here, and the preview is run again rather than showing a stale product.
fn refresh_document(inputs: &mut CloningInputs, session: &Session, side: &CloningSide) {
    let Some(id) = side.document else {
        return;
    };
    let key = side.source_key();
    let revision = session.revision();
    if inputs.info_rev.get(&key) == Some(&revision) {
        return;
    }
    inputs.info_rev.insert(key.clone(), revision);
    let found = session.document_by_id(id).map(|d| &d.molecule);
    // Only a real change costs anything: comparing against what is already shown keeps a
    // preview that is still valid, and can't loop asking for a new one.
    let known = inputs.info.get(&key).map(|i| match i {
        Ok(info) => info.molecule.as_ref(),
        Err(_) => None,
    });
    if let (Some(m), Some(Some(before))) = (found, known)
        && m == before
    {
        return;
    }
    let entry = match found {
        Some(m) => Ok(ConstructInfo {
            name: m.name().to_string(),
            length: m.len() as u64,
            topology: match m.topology() {
                Topology::Circular => "circular".to_string(),
                Topology::Linear => "linear".to_string(),
            },
            features: m.features().len() as u64,
            molecule: Some(m.clone()),
        }),
        None => Err("this tab isn't open any more; choose another construct".to_string()),
    };
    inputs.info.insert(key.clone(), entry);
    inputs
        .cutters
        .retain(|k, _| !k.ends_with(&format!("\u{1}{key}")));
    inputs.forget_preview();
}

/// The enzymes that cut the construct at `path`, from the lab's set when there is one.
/// Computed once per construct and set; empty until the construct has been read.
fn cutters(inputs: &mut CloningInputs, session: &Session, source: &str) -> Vec<Cutter> {
    let source = source.trim();
    if source.is_empty() {
        return Vec::new();
    }
    let key = format!("{}\u{1}{source}", inputs.pool.as_deref().unwrap_or(""));
    if let Some(found) = inputs.cutters.get(&key) {
        return found.clone();
    }
    let Some(Ok(info)) = inputs.info.get(source) else {
        return Vec::new();
    };
    let Some(m) = &info.molecule else {
        return Vec::new();
    };
    let query = snappiesnap_engine::enzymes::SiteQuery {
        enzymes: None,
        set: inputs.pool.clone(),
        min_length: MIN_OFFERED_LENGTH,
        max_cuts: MAX_OFFERED_CUTS,
        commercial_only: true,
    };
    // An unknown set (deleted since) answers with an error: offer nothing rather than the
    // whole database, which would propose enzymes the lab hasn't got.
    let mut found: Vec<Cutter> = snappiesnap_engine::enzymes::sites(session.enzymes(), m, &query)
        .unwrap_or_default()
        .into_iter()
        .filter(|e| !e.sites.is_empty())
        .map(|e| Cutter {
            name: e.name,
            cuts: e.sites.len(),
        })
        .collect();
    found.sort_by(|a, b| a.cuts.cmp(&b.cuts).then_with(|| a.name.cmp(&b.name)));
    found.truncate(MAX_OFFERED);
    inputs.cutters.insert(key, found.clone());
    found
}

/// The restriction sites to draw on a construct's map: the enzymes typed in its fields,
/// which are the ones the choice is about, and otherwise the cutters the lab could use (the
/// chosen enzyme set, filtered as the main map filters it). Cached per construct, set and
/// typed names, because a map is redrawn every frame and a site search is not free.
fn sites_for(
    inputs: &mut CloningInputs,
    session: &Session,
    source: &str,
    typed: &[String],
) -> Vec<crate::enzymes::MapSite> {
    let source = source.trim();
    if source.is_empty() {
        return Vec::new();
    }
    let signature = format!(
        "{}\u{1}{}",
        inputs.pool.as_deref().unwrap_or(""),
        typed.join(",")
    );
    if let Some((had, found)) = inputs.sites.get(source)
        && had == &signature
    {
        return found.clone();
    }
    let Some(Ok(info)) = inputs.info.get(source) else {
        return Vec::new();
    };
    let Some(m) = &info.molecule else {
        return Vec::new();
    };
    let known: Vec<String> = typed
        .iter()
        .filter(|n| session.enzymes().get(n).is_some())
        .cloned()
        .collect();
    // Every enzyme that cuts the construct is drawn, so the choice is made from what is
    // actually there; the ones named in the cut list are marked, and the map draws those in
    // the accent colour and in bold.
    let query = snappiesnap_engine::enzymes::SiteQuery {
        enzymes: None,
        set: inputs.pool.clone(),
        min_length: MIN_OFFERED_LENGTH,
        max_cuts: MAX_OFFERED_CUTS,
        commercial_only: true,
    };
    let mut found = crate::enzymes::map_sites(session.enzymes(), m, &query);
    // Enzymes named in the fields are shown even when the filters above would leave them
    // out, because they are the ones the step is about.
    let extra: Vec<String> = known
        .iter()
        .filter(|n| {
            !found.iter().any(|s| {
                s.detail
                    .names
                    .split(", ")
                    .any(|x| x.eq_ignore_ascii_case(n))
            })
        })
        .cloned()
        .collect();
    if !extra.is_empty() {
        let named = snappiesnap_engine::enzymes::SiteQuery {
            enzymes: Some(extra),
            ..snappiesnap_engine::enzymes::SiteQuery::default()
        };
        found.extend(crate::enzymes::map_sites(session.enzymes(), m, &named));
    }
    for site in &mut found {
        site.chosen = site
            .detail
            .names
            .split(", ")
            .any(|n| known.iter().any(|k| k.eq_ignore_ascii_case(n)));
    }
    found.sort_by_key(|x| x.mark.pos);
    if inputs.sites.len() > MAX_CACHED {
        inputs.sites.clear();
    }
    inputs
        .sites
        .insert(source.to_string(), (signature, found.clone()));
    found
}

/// Add the enzyme of site `index` on this side's construct to its cut list. Returns whether
/// anything changed: double-clicking an enzyme that is already named leaves the field alone.
fn add_site_enzyme(inputs: &CloningInputs, side: &mut CloningSide, index: usize) -> bool {
    let key = side.source_key();
    let Some(name) = inputs
        .sites
        .get(&key)
        .and_then(|(_, found)| found.get(index))
        .map(|s| s.detail.label.clone())
    else {
        return false;
    };
    let before = side.enzymes.clone();
    add_name(&mut side.enzymes, &name);
    before != side.enzymes
}

/// Put the primer at `index` in the forward or reverse field of this side, by the strand it
/// binds: a pair is one primer on each strand. Returns whether anything changed.
fn add_primer(inputs: &CloningInputs, side: &mut CloningSide, key: &str, index: usize) -> bool {
    let Some((_, primers)) = inputs.extras.get(key) else {
        return false;
    };
    let Some(mark) = primers.get(index) else {
        return false;
    };
    let field = if mark.forward {
        &mut side.forward
    } else {
        &mut side.reverse
    };
    if field.trim() == mark.name {
        return false;
    }
    *field = mark.name.clone();
    // A primer is only used when the side is made by PCR, so choosing one says so.
    if side.using != SideUse::Pcr {
        side.using = SideUse::Pcr;
    }
    true
}

/// The ORFs and primer binding sites of a construct, kept per construct so the layers have
/// something to draw without being found again on every frame. The ORF and primer toggles
/// were switched on in 0.26.0 with empty lists behind them, which drew nothing.
fn extras_for(
    inputs: &mut CloningInputs,
    session: &Session,
    key: &str,
    m: &Molecule,
) -> (
    Vec<crate::enzymes::OrfMark>,
    Vec<crate::enzymes::PrimerMark>,
) {
    if let Some(found) = inputs.extras.get(key) {
        return found.clone();
    }
    let orfs = crate::enzymes::orf_marks(m);
    let primers = crate::enzymes::primer_marks(m, session.primer_binding());
    if inputs.extras.len() > MAX_CACHED {
        inputs.extras.clear();
    }
    inputs
        .extras
        .insert(key.to_string(), (orfs.clone(), primers.clone()));
    (orfs, primers)
}

/// The sites to label on the product, cached by its length so the map isn't recomputed on
/// every frame. The product has no enzyme fields of its own, so it shows the cutters the lab
/// could use, the same way a construct does before enzymes are typed.
fn product_sites(
    inputs: &mut CloningInputs,
    session: &Session,
    key: &str,
    m: &Molecule,
) -> Vec<crate::enzymes::MapSite> {
    let signature = inputs.pool.clone().unwrap_or_default();
    if let Some((had, found)) = inputs.sites.get(key)
        && had == &signature
    {
        return found.clone();
    }
    let query = snappiesnap_engine::enzymes::SiteQuery {
        enzymes: None,
        set: inputs.pool.clone(),
        min_length: MIN_OFFERED_LENGTH,
        max_cuts: 1,
        commercial_only: true,
    };
    let mut found = crate::enzymes::map_sites(session.enzymes(), m, &query);
    if found.len() > MAX_MAP_LABELS {
        let db = session.enzymes();
        found.sort_by(|a, b| {
            db.naming_rank(&a.detail.label)
                .cmp(&db.naming_rank(&b.detail.label))
                .then_with(|| a.detail.label.cmp(&b.detail.label))
        });
        found.truncate(MAX_MAP_LABELS);
        found.sort_by_key(|x| x.mark.pos);
    }
    if inputs.sites.len() > MAX_CACHED {
        inputs.sites.clear();
    }
    inputs
        .sites
        .insert(key.to_string(), (signature, found.clone()));
    found
}

/// What is wrong with the enzymes typed in a field, if anything: a name the database hasn't
/// got (the reaction will fail), or one the lab's set hasn't got (it has to be bought or
/// borrowed). Returns the note and whether it is a warning rather than an error.
fn enzyme_note(text: &str, offer: &Offer<'_>) -> Option<(String, bool)> {
    let typed = names(text);
    let unknown: Vec<&str> = typed
        .iter()
        .filter(|n| offer.db.get(n).is_none())
        .map(String::as_str)
        .collect();
    if !unknown.is_empty() {
        return Some((
            format!(
                "Not in the enzyme database, so the reaction can't run: {}. Download REBASE in \
                 the Enzymes tab for the full list.",
                unknown.join(", ")
            ),
            false,
        ));
    }
    let pool = offer.pool?;
    let set = offer.db.set(pool)?;
    let absent: Vec<&str> = typed
        .iter()
        .filter(|n| !set.enzymes.iter().any(|e| e.eq_ignore_ascii_case(n)))
        .map(String::as_str)
        .collect();
    (!absent.is_empty()).then(|| {
        (
            format!(
                "Not in \"{pool}\": you'd have to buy or borrow {}.",
                absent.join(", ")
            ),
            true,
        )
    })
}

/// An enzyme field: what is typed, and a list of the enzymes that cut this construct to add
/// from. With a set chosen, the list holds only what the lab has, so what it offers can be
/// run; what is typed is never restricted.
fn enzyme_field(
    ui: &mut Ui,
    t: &Tokens,
    id: &str,
    text: &mut String,
    hint: &str,
    width: f32,
    offer: &Offer<'_>,
) {
    ui.add(
        egui::TextEdit::singleline(text)
            .hint_text(hint)
            .desired_width(width),
    );
    if offer.cutters.is_empty() {
        return;
    }
    let tip = match offer.pool {
        Some(pool) => format!("The enzymes of \"{pool}\" that cut this construct"),
        None => "The enzymes that cut this construct".to_string(),
    };
    egui::ComboBox::from_id_salt(id)
        .selected_text("Add\u{2026}")
        .width(110.0)
        .show_ui(ui, |ui| {
            for c in offer.cutters {
                let plural = if c.cuts == 1 { "" } else { "s" };
                let label = format!("{} ({} site{plural})", c.name, c.cuts);
                let r = ui.selectable_label(false, RichText::new(label));
                if r.clicked() {
                    add_name(text, &c.name);
                    ui.close();
                }
            }
        })
        .response
        .on_hover_text(tip);
    if offer.pool.is_none() {
        ui.label(RichText::new("every enzyme").color(t.text_muted))
            .on_hover_text(
                "No enzyme set is chosen, so every commercial enzyme is offered. Choose a set \
                 of your own in Preferences to be offered only what the lab holds.",
            );
    }
}

/// Show [`enzyme_note`] under an enzyme field, when there is one.
fn note(ui: &mut Ui, t: &Tokens, text: &str, offer: &Offer<'_>) {
    if let Some((message, warn)) = enzyme_note(text, offer) {
        ui.label(RichText::new(message).color(if warn { t.warning } else { t.error }));
    }
}

/// Name one side of the reaction the way the command takes it: an open document by its tab
/// id, otherwise the file it was chosen from. The document wins, so a construct edited in its
/// tab is cloned as it stands, and so the web build (which has no file paths of its own) can
/// clone from what is open.
fn source(o: &mut serde_json::Map<String, Value>, field: &str, side: &CloningSide) {
    match side.document {
        Some(id) => {
            o.insert(format!("{field}_document"), json!(id));
        }
        None => {
            o.insert(field.to_string(), json!(side.file));
        }
    }
}

/// One part of an assembly, as `assemble` takes it.
fn part(side: &CloningSide) -> Value {
    let mut v = json!({});
    if let Some(o) = v.as_object_mut() {
        match side.document {
            Some(id) => {
                o.insert("document".into(), json!(id));
            }
            None => {
                o.insert("file".into(), json!(side.file));
            }
        }
        match side.using {
            SideUse::Cut => {
                o.insert("enzymes".into(), json!(names(&side.enzymes)));
            }
            SideUse::Pcr => {
                o.insert(
                    "pcr".into(),
                    json!({"forward": side.forward, "reverse": side.reverse}),
                );
            }
            SideUse::Whole => {}
        }
    }
    v
}

/// The command and params for the experiment as it stands, previewing or for real.
/// `None` while something it needs is still missing.
pub fn call(exp: &Experiment, preview: bool) -> Option<(&'static str, Value)> {
    let t = exp.technique;
    if !exp.vector.is_set() {
        return None;
    }
    let command = t.command();
    let mut params = match t {
        CloningTechnique::Restriction => {
            if !exp.insert.is_set() || names(&exp.vector.enzymes).is_empty() {
                return None;
            }
            let mut p = json!({ "vector_enzymes": names(&exp.vector.enzymes) });
            if let Some(o) = p.as_object_mut() {
                source(o, "vector", &exp.vector);
                source(o, "insert", &exp.insert);
                let cuts = names(&exp.insert.enzymes);
                if !cuts.is_empty() {
                    o.insert("insert_enzymes".into(), json!(cuts));
                }
            }
            p
        }
        CloningTechnique::GatewayLr | CloningTechnique::GatewayBp => {
            if !exp.insert.is_set() {
                return None;
            }
            let mut p = json!({
                "reaction": if t == CloningTechnique::GatewayBp { "bp" } else { "lr" },
            });
            if let Some(o) = p.as_object_mut() {
                source(o, "acceptor", &exp.vector);
                source(o, "donor", &exp.insert);
            }
            p
        }
        CloningTechnique::Gibson | CloningTechnique::Infusion | CloningTechnique::GoldenGate => {
            if !exp.insert.is_set() {
                return None;
            }
            let mut parts = vec![part(&exp.vector), part(&exp.insert)];
            parts.extend(exp.extra.iter().filter(|s| s.is_set()).map(part));
            let method = match t {
                CloningTechnique::Gibson => "gibson",
                CloningTechnique::Infusion => "infusion",
                _ => "golden_gate",
            };
            let mut p = json!({ "method": method, "parts": parts });
            if let Some(o) = p.as_object_mut()
                && t == CloningTechnique::GoldenGate
            {
                let enzyme = exp.vector.enzymes.trim();
                o.insert(
                    "enzyme".into(),
                    json!(if enzyme.is_empty() { "BsaI" } else { enzyme }),
                );
            }
            p
        }
        CloningTechnique::TopoTa
        | CloningTechnique::TopoBlunt
        | CloningTechnique::TopoDirectional => {
            let kind = match t {
                CloningTechnique::TopoBlunt => "blunt",
                CloningTechnique::TopoDirectional => "directional",
                _ => "ta",
            };
            let pcr = exp.insert.using == SideUse::Pcr;
            if !exp.insert.is_set() {
                return None;
            }
            let mut p = json!({ "kind": kind });
            if let Some(o) = p.as_object_mut() {
                source(o, "vector", &exp.vector);
                source(o, "insert", &exp.insert);
                if let Ok(at) = exp.vector.position.trim().parse::<usize>() {
                    o.insert("position".into(), json!(at));
                }
                if pcr {
                    o.insert(
                        "pcr".into(),
                        json!({"forward": exp.insert.forward, "reverse": exp.insert.reverse}),
                    );
                }
            }
            p
        }
        CloningTechnique::Mutagenesis => {
            if exp.insert.forward.trim().is_empty() || exp.insert.reverse.trim().is_empty() {
                return None;
            }
            let mut p = json!({
                "forward": exp.insert.forward.trim(),
                "reverse": exp.insert.reverse.trim(),
            });
            if let Some(o) = p.as_object_mut() {
                source(o, "template", &exp.vector);
            }
            p
        }
    };
    if let Some(o) = params.as_object_mut() {
        if preview {
            o.insert("preview".into(), json!(true));
        } else {
            // The experiment is bound to no document, so the product opens in a new tab and
            // nothing of the user's is replaced.
            o.insert("discard_changes".into(), json!(true));
        }
    }
    Some((command, params))
}

/// The constructs of the open folder and its subfolders, cached per folder.
fn folder_files(
    inputs: &mut CloningInputs,
    session: &Session,
    folder: Option<&str>,
) -> Vec<String> {
    let Some(folder) = folder else {
        return Vec::new();
    };
    if inputs.listed.as_deref() == Some(folder) {
        return inputs.files.clone();
    }
    let mut out = Vec::new();
    let mut queue = vec![(folder.to_string(), 0usize)];
    let mut folders = 0usize;
    while let Some((dir, depth)) = queue.pop() {
        folders += 1;
        if folders > MAX_FOLDERS || out.len() >= MAX_FILES {
            break;
        }
        let Ok(entries) = session.platform().list_dir(&dir) else {
            continue;
        };
        for e in entries {
            if e.is_dir {
                if depth < MAX_DEPTH {
                    queue.push((e.path.clone(), depth + 1));
                }
            } else if crate::sidebar::is_sequence_file(&e.name) && out.len() < MAX_FILES {
                out.push(e.path);
            }
        }
    }
    out.sort();
    inputs.files.clone_from(&out);
    inputs.listed = Some(folder.to_string());
    out
}

/// A construct's map, at the size the steps use.
fn map(
    ui: &mut Ui,
    t: &Tokens,
    m: &Molecule,
    it: &Shown<'_>,
    selection: Option<Range>,
    zoom: &mut f32,
    pan: &mut egui::Vec2,
) -> Option<crate::map::MapAction> {
    let (size, pane) = (it.size, it.pane);
    let marks = crate::enzymes::Marks {
        sites: it.sites,
        primers: it.primers,
        orfs: it.orfs,
    };
    let layers = crate::map::Layers {
        // The enzymes are what the choice is about, so they are on unless switched off.
        enzymes: pane.layers.enzymes && !it.sites.is_empty(),
        labels: true,
        ..pane.layers
    };
    let circular = match pane.style {
        crate::state::MapStyle::Auto => m.topology() == Topology::Circular,
        crate::state::MapStyle::Circular => true,
        crate::state::MapStyle::Linear => false,
    };
    let mut action = None;
    ui.allocate_ui(size, |ui| {
        ui.set_min_size(size);
        let area = ui.available_rect_before_wrap();
        // The same zoom and pan a document tab's map has, over its own state.
        action = crate::app::zoomed(ui, "cloning-map", area, zoom, pan, circular, |ui| {
            crate::map::show(ui, t, m, selection, pane.style, &marks, layers)
        });
    });
    action
}

/// The breadcrumb of the four steps.
fn steps(ui: &mut Ui, t: &Tokens, exp: &Experiment, current: u8) {
    let names = [
        "Technique",
        exp.technique.vector_label(),
        exp.technique.insert_label(),
        "Product",
    ];
    ui.horizontal(|ui| {
        for (i, label) in names.iter().enumerate() {
            let n = u8::try_from(i + 1).unwrap_or(u8::MAX);
            if i > 0 {
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

/// Step 1: which technique.
fn technique(
    ui: &mut Ui,
    t: &Tokens,
    exp: &mut Experiment,
    inputs: &mut CloningInputs,
) -> Option<CloningAction> {
    steps(ui, t, exp, 1);
    ui.label(RichText::new("Which cloning would you like to simulate?").strong());
    ui.add_space(t.spacing);
    egui::ScrollArea::vertical()
        .auto_shrink([false, true])
        .max_height(360.0)
        .show(ui, |ui| {
            for technique in CloningTechnique::ALL {
                let chosen = exp.technique == technique;
                if ui
                    .selectable_label(chosen, RichText::new(technique.label()).strong())
                    .clicked()
                    && !chosen
                {
                    exp.technique = technique;
                    inputs.forget_preview();
                }
                ui.label(RichText::new(technique.description()).color(t.text_muted));
                ui.add_space(t.spacing);
            }
        });
    ui.add_space(t.spacing);
    let mut action = None;
    ui.horizontal(|ui| {
        if ui
            .button(format!("Continue: {}", exp.technique.vector_label()))
            .clicked()
        {
            exp.step = ExperimentStep::Vector;
        }
        // Somewhere to start with nothing of your own open: the sample plasmid and a fragment
        // that fits its polylinker, set up as a restriction clone ready to run. In the browser
        // this is the whole of a first try, because there are no files on a web page.
        if ui
            .button("Try it with the sample constructs")
            .on_hover_text(
                "Opens pSample and demoInsert and fills in a restriction clone with EcoRI and                  BamHI, so you can walk the steps through to a product",
            )
            .clicked()
        {
            action = Some(CloningAction::Samples);
        }
    });
    action
}

/// What one side's editor needs besides the side itself: what can be chosen, what the
/// technique does with it, whether it is the vector, and what it can be cut with.
#[derive(Clone, Copy)]
struct Side<'a> {
    choices: &'a [Choice],
    technique: CloningTechnique,
    is_vector: bool,
    offer: &'a Offer<'a>,
}

/// The picker and preparation of one construct. Returns true when something changed.
fn side_editor(ui: &mut Ui, t: &Tokens, id: &str, side: &mut CloningSide, ctx: &Side<'_>) -> bool {
    let Side {
        choices,
        technique,
        is_vector,
        offer,
    } = *ctx;
    let before = side.clone();
    // The controls live in a third of the window now, so every row wraps instead of pushing
    // the construct's map off the side.
    ui.horizontal_wrapped(|ui| {
        let shown = match (&side.document, side.is_set()) {
            (Some(_), _) => match chosen_name(side, choices) {
                Some(name) => format!("{name} (open)"),
                None => "A tab that was closed".to_string(),
            },
            (None, true) => side.file_name().to_string(),
            (None, false) => "Choose a construct".to_string(),
        };
        egui::ComboBox::from_id_salt(id)
            .width(200.0)
            .selected_text(shown)
            .show_ui(ui, |ui| {
                if choices.is_empty() {
                    ui.label(
                        RichText::new("Nothing to choose: open a construct first")
                            .color(t.text_muted),
                    );
                }
                let mut open_shown = false;
                let mut files_shown = false;
                for c in choices {
                    match c {
                        Choice::Document { id, name, path } => {
                            if !open_shown {
                                ui.label(RichText::new("Open now").color(t.text_muted));
                                open_shown = true;
                            }
                            if ui
                                .selectable_label(side.document == Some(*id), name)
                                .clicked()
                            {
                                side.document = Some(*id);
                                // The file is kept beside the tab, so an experiment saved
                                // now and reopened later still finds the construct.
                                side.file = path.clone().unwrap_or_default();
                            }
                        }
                        Choice::File(f) => {
                            if !files_shown {
                                ui.label(RichText::new("In the folder").color(t.text_muted));
                                files_shown = true;
                            }
                            let picked = side.document.is_none() && &side.file == f;
                            if ui
                                .selectable_label(picked, crate::sidebar::base_name(f))
                                .clicked()
                            {
                                side.document = None;
                                side.file = f.clone();
                            }
                        }
                    }
                }
            });
        #[cfg(not(target_arch = "wasm32"))]
        if ui
            .button("Open from disk\u{2026}")
            .on_hover_text("Choose a construct outside the open folder")
            .clicked()
            && let Some(p) = rfd::FileDialog::new()
                .add_filter("Sequence files", crate::sidebar::SEQUENCE_EXTENSIONS)
                .pick_file()
        {
            side.document = None;
            side.file = p.to_string_lossy().into_owned();
        }
    });
    if !side.is_set() {
        return before != *side;
    }
    ui.add_space(t.spacing);
    // How this side is prepared, as far as the technique allows it.
    match technique {
        CloningTechnique::Restriction => {
            ui.horizontal_wrapped(|ui| {
                ui.label(if is_vector {
                    "Cut the vector with"
                } else {
                    "Cut the insert with"
                });
                enzyme_field(
                    ui,
                    t,
                    &format!("{id}-cut"),
                    &mut side.enzymes,
                    if is_vector {
                        "EcoRI, BamHI"
                    } else {
                        "the vector's"
                    },
                    220.0,
                    offer,
                );
            });
            if !is_vector && names(&side.enzymes).is_empty() {
                ui.label(
                    RichText::new("Left empty, the insert is cut with the vector's enzymes.")
                        .color(t.text_muted),
                );
            }
            note(ui, t, &side.enzymes, offer);
        }
        CloningTechnique::Gibson | CloningTechnique::Infusion | CloningTechnique::GoldenGate => {
            ui.horizontal_wrapped(|ui| {
                ui.label("Use it");
                ui.selectable_value(&mut side.using, SideUse::Whole, "whole");
                if technique != CloningTechnique::GoldenGate {
                    ui.selectable_value(&mut side.using, SideUse::Cut, "cut with enzymes");
                    ui.selectable_value(&mut side.using, SideUse::Pcr, "as a PCR product");
                }
            });
            match side.using {
                SideUse::Cut => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Enzymes");
                        enzyme_field(
                            ui,
                            t,
                            &format!("{id}-cut"),
                            &mut side.enzymes,
                            "EcoRI",
                            140.0,
                            offer,
                        );
                        ui.label(RichText::new("the largest fragment is kept").color(t.text_muted));
                    });
                    note(ui, t, &side.enzymes, offer);
                }
                SideUse::Pcr => {
                    primer_fields(ui, side, "Forward", "Reverse");
                }
                SideUse::Whole => {}
            }
            if technique == CloningTechnique::GoldenGate && is_vector {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Type IIS enzyme");
                    enzyme_field(
                        ui,
                        t,
                        &format!("{id}-typeiis"),
                        &mut side.enzymes,
                        "BsaI",
                        120.0,
                        offer,
                    );
                });
                note(ui, t, &side.enzymes, offer);
            }
        }
        CloningTechnique::TopoTa
        | CloningTechnique::TopoBlunt
        | CloningTechnique::TopoDirectional => {
            if is_vector {
                ui.horizontal(|ui| {
                    ui.label("Open a circular vector at");
                    ui.add(
                        egui::TextEdit::singleline(&mut side.position)
                            .hint_text("position")
                            .desired_width(80.0),
                    );
                    ui.label(
                        RichText::new("leave empty for a vector that is already linear")
                            .color(t.text_muted),
                    );
                });
            } else {
                ui.horizontal(|ui| {
                    ui.label("Use it");
                    ui.selectable_value(&mut side.using, SideUse::Whole, "as it is");
                    ui.selectable_value(&mut side.using, SideUse::Pcr, "as a PCR product");
                });
                if side.using == SideUse::Pcr {
                    primer_fields(ui, side, "Forward", "Reverse");
                }
            }
        }
        CloningTechnique::GatewayLr
        | CloningTechnique::GatewayBp
        | CloningTechnique::Mutagenesis => {}
    }
    before != *side
}

/// The two primer fields of a side.
fn primer_fields(ui: &mut Ui, side: &mut CloningSide, forward: &str, reverse: &str) {
    ui.horizontal(|ui| {
        ui.label(forward);
        ui.add(
            egui::TextEdit::singleline(&mut side.forward)
                .hint_text("name or sequence")
                .desired_width(170.0),
        );
    });
    ui.horizontal(|ui| {
        ui.label(reverse);
        ui.add(
            egui::TextEdit::singleline(&mut side.reverse)
                .hint_text("name or sequence")
                .desired_width(170.0),
        );
    });
}

/// The row above the construct: which view, the map's zoom, and (on the sequence) what the
/// ORF layer translates. The same controls a document tab has, in the same order.
fn view_controls(ui: &mut Ui, t: &Tokens, pane: &mut Pane, inputs: &mut CloningInputs) {
    ui.horizontal_wrapped(|ui| {
        if ui
            .selectable_label(!pane.sequence, "Map")
            .on_hover_text("Show the construct as a map")
            .clicked()
        {
            pane.sequence = false;
        }
        if ui
            .selectable_label(pane.sequence, "Sequence")
            .on_hover_text("Show the construct's residues, with its features and sites")
            .clicked()
        {
            pane.sequence = true;
        }
        ui.separator();
        if pane.sequence {
            return;
        }
        // Zoom, as on a document tab.
        if ui
            .small_button("\u{2212}")
            .on_hover_text("Zoom out")
            .clicked()
        {
            inputs.zoom = (inputs.zoom / 1.5).max(1.0);
        }
        ui.label(RichText::new(format!("{:.0}%", inputs.zoom * 100.0)).color(t.text_muted))
            .on_hover_text("Zoom: pinch or Cmd/Ctrl+scroll over the map; drag to move it");
        if ui.small_button("+").on_hover_text("Zoom in").clicked() {
            inputs.zoom = (inputs.zoom * 1.5).min(crate::state::MAX_ZOOM);
        }
        if ui
            .add_enabled(inputs.zoom > 1.0, egui::Button::new("Fit").small())
            .clicked()
        {
            inputs.zoom = 1.0;
            inputs.pan = egui::Vec2::ZERO;
        }
        ui.separator();
        // How the map is drawn, as in the document views.
        for (style, label) in [
            (crate::state::MapStyle::Auto, "Auto"),
            (crate::state::MapStyle::Circular, "Circular"),
            (crate::state::MapStyle::Linear, "Linear"),
        ] {
            if ui.selectable_label(pane.style == style, label).clicked() {
                pane.style = style;
            }
        }
    });
}

/// The layer toggles down the left edge, as the document views have them.
fn layer_column(ui: &mut Ui, t: &Tokens, pane: &mut Pane) {
    use crate::icons::Icon;
    ui.vertical(|ui| {
        ui.spacing_mut().item_spacing.y = 4.0;
        let l = &mut pane.layers;
        for (icon, tip, on) in [
            (Icon::Features, "Features", &mut l.features),
            (Icon::Enzymes, "Restriction sites", &mut l.enzymes),
            (Icon::Primers, "Primer binding sites", &mut l.primers),
            (Icon::Orfs, "Open reading frames", &mut l.orfs),
            (Icon::Ruler, "Ruler", &mut l.ruler),
            (
                Icon::Positions,
                "Cut positions in site labels",
                &mut l.positions,
            ),
        ] {
            if crate::app::layer_toggle(ui, t, icon, tip, *on).clicked() {
                *on = !*on;
            }
        }
    });
}

/// What the construct pane draws: how much room it has, which sites, ORFs and primer sites
/// to show, and the view and layers chosen for it.
struct Shown<'a> {
    size: egui::Vec2,
    sites: &'a [crate::enzymes::MapSite],
    orfs: &'a [crate::enzymes::OrfMark],
    primers: &'a [crate::enzymes::PrimerMark],
    session: &'a Session,
    pane: &'a Pane,
}

/// What is known about a chosen construct, and its map.
fn construct_card(
    ui: &mut Ui,
    t: &Tokens,
    inputs: &mut CloningInputs,
    path: &str,
    it: &Shown<'_>,
) -> Option<crate::map::MapAction> {
    let mut action = None;
    let info = inputs.info.get(path).cloned();
    match &info {
        None => {
            ui.label(RichText::new("Reading the construct\u{2026}").color(t.text_muted));
        }
        Some(Err(e)) => {
            ui.label(RichText::new(format!("This file can't be read: {e}")).color(t.error));
        }
        Some(Ok(info)) => {
            ui.label(
                RichText::new(format!(
                    "{}  \u{b7}  {} bp  \u{b7}  {}  \u{b7}  {} feature(s)",
                    info.name, info.length, info.topology, info.features
                ))
                .color(t.text_muted),
            );
            if let Some(m) = &info.molecule {
                if it.pane.sequence {
                    sequence(ui, t, m, path, inputs, it);
                } else {
                    let selection = inputs.selection_in(path);
                    let (mut zoom, mut pan) = (inputs.zoom, inputs.pan);
                    action = map(ui, t, m, it, selection, &mut zoom, &mut pan);
                    inputs.zoom = zoom;
                    inputs.pan = pan;
                }
            }
        }
    }
    action
}

/// The construct's residues, as the sequence view draws them. The construct belongs to a file
/// or to another tab, so edits typed here are not applied: the note says where to make them.
fn sequence(
    ui: &mut Ui,
    t: &Tokens,
    m: &Molecule,
    key: &str,
    inputs: &mut CloningInputs,
    it: &Shown<'_>,
) {
    use crate::seqview::SeqAction;
    let marks = crate::enzymes::Marks {
        sites: it.sites,
        primers: it.primers,
        orfs: it.orfs,
    };
    // A selection belongs to the construct it was made in.
    let selection = inputs.selection_in(key);
    let mut actions = Vec::new();
    ui.allocate_ui(it.size, |ui| {
        ui.set_min_size(it.size);
        // The ORF layer shows ORFs or translated frames, exactly as in a document tab.
        let mut layers = it.pane.layers;
        layers.orfs = it.pane.layers.orfs && it.pane.translation == crate::state::Translation::Orfs;
        let frames: &[i8] = if it.pane.layers.orfs {
            it.pane.translation.frames()
        } else {
            &[]
        };
        actions = crate::seqview::show(
            ui,
            t,
            m,
            selection,
            &mut inputs.seq,
            it.session.revision(),
            &marks,
            layers,
            frames,
            it.pane.genetic_code,
            &[],
            it.pane.tails,
            it.pane.single_line,
        );
    });
    let mut edited = false;
    for action in actions {
        match action {
            // The same pick the map offers: a double-clicked enzyme joins the cut list.
            SeqAction::SiteChosen(i) => inputs.chosen_in_sequence = Some(i),
            // Reading the construct: selecting and copying work here, from the construct
            // shown rather than from the open document.
            SeqAction::Select(r) => inputs.select_in(key, r),
            SeqAction::Copy(format) => {
                let Some(range) = selection.filter(|r| !r.is_empty()) else {
                    inputs.copy_note = Some("Select some residues first.".to_string());
                    continue;
                };
                match snappiesnap_engine::copy_text(m, range, format) {
                    Ok((text, n)) => {
                        ui.ctx().copy_text(text);
                        inputs.copy_note = Some(format!("Copied {n} residues"));
                    }
                    Err(e) => inputs.copy_note = Some(e.message),
                }
            }
            // Changing it: not here, and said so rather than quietly dropped.
            _ => edited = true,
        }
    }
    if edited {
        inputs.read_only_note = true;
    }
    if let Some(note) = &inputs.copy_note {
        ui.label(RichText::new(note).small().color(t.text_muted));
    }
    if inputs.read_only_note {
        ui.label(
            RichText::new(
                "Selecting and copying work here. To change this construct, open its own tab: \
                 the experiment follows it.",
            )
            .small()
            .color(t.warning),
        );
    }
}

/// The lengths of the fragments `enzymes` cut `source` into, in order along the molecule.
/// Empty when the construct hasn't been read, no enzyme is named, or none of them cuts: the
/// schematic then says nothing about pieces instead of inventing them. Cached per construct
/// and enzyme list, as the digest runs every frame otherwise.
fn fragments_for(
    inputs: &mut CloningInputs,
    session: &Session,
    source: &str,
    enzymes: &[String],
) -> Vec<Frag> {
    let source = source.trim();
    if source.is_empty() || enzymes.is_empty() {
        return Vec::new();
    }
    let signature = enzymes.join(",");
    if let Some((had, found)) = inputs.fragments.get(source)
        && had == &signature
    {
        return found.clone();
    }
    let Some(Ok(info)) = inputs.info.get(source) else {
        return Vec::new();
    };
    let Some(m) = &info.molecule else {
        return Vec::new();
    };
    let chosen: Vec<snappiesnap_engine::bio::enzymes::Enzyme> = enzymes
        .iter()
        .filter_map(|n| session.enzymes().get(n).map(|(_, e)| e.clone()))
        .collect();
    let found: Vec<Frag> = snappiesnap_engine::bio::cloning::digest_fragments(m, &chosen)
        .iter()
        .map(|f| Frag {
            len: f.len() as u64,
            left: f.left.clone(),
            right: f.right.clone(),
        })
        .collect();
    inputs
        .fragments
        .insert(source.to_string(), (signature, found.clone()));
    found
}

/// The schematic of the reaction as it stands, from what has been read so far.
fn plan_of(
    exp: &Experiment,
    inputs: &mut CloningInputs,
    session: &Session,
) -> crate::cloning_plan::Plan {
    let read = |side: &CloningSide| {
        inputs
            .info
            .get(&side.source_key())
            .and_then(|i| i.as_ref().ok())
    };
    let length = |side: &CloningSide| read(side).map(|i| i.length);
    // `file_info` reports the topology as a word; anything else is left unsaid.
    let circular = |side: &CloningSide| match read(side).map(|i| i.topology.as_str()) {
        Some("circular") => Some(true),
        Some("linear") => Some(false),
        _ => None,
    };
    let product = inputs
        .preview
        .as_ref()
        .and_then(|p| p.as_ref().ok())
        .map(|p| p.length);
    let measured = crate::cloning_plan::Measured {
        vector_len: length(&exp.vector),
        insert_len: length(&exp.insert),
        product_len: product,
        vector_circular: circular(&exp.vector),
        insert_circular: circular(&exp.insert),
        vector_fragments: Vec::new(),
        insert_fragments: Vec::new(),
        vector_name: Some(side_label(session, &exp.vector)),
        insert_name: Some(side_label(session, &exp.insert)),
        fits: None,
    };
    let vector_enzymes = names(&exp.vector.enzymes);
    // An insert with no enzymes of its own is cut with the vector's, as the reaction does it.
    let insert_enzymes = {
        let own = names(&exp.insert.enzymes);
        if own.is_empty() {
            vector_enzymes.clone()
        } else {
            own
        }
    };
    // Only a side that is actually cut has fragments: one used whole or as a PCR product
    // isn't digested, whatever enzymes the other side names.
    let cut_with = |side: &CloningSide, enzymes: &[String]| -> Vec<String> {
        match side.using {
            SideUse::Cut => enzymes.to_vec(),
            SideUse::Pcr => Vec::new(),
            // Restriction cloning digests both sides with the enzymes named for them and has
            // no "use it whole" choice; the assembly and TOPO techniques do, and there only
            // `cut with enzymes` means a digest.
            SideUse::Whole if exp.technique == CloningTechnique::Restriction => enzymes.to_vec(),
            SideUse::Whole => Vec::new(),
        }
    };
    let (vk, ik) = (exp.vector.source_key(), exp.insert.source_key());
    let vf = cut_with(&exp.vector, &vector_enzymes);
    let if_ = cut_with(&exp.insert, &insert_enzymes);
    let vector_cut = fragments_for(inputs, session, &vk, &vf);
    let insert_cut = fragments_for(inputs, session, &ik, &if_);
    // Whether the piece goes into the backbone, and whether the way round is settled. The
    // reaction keeps the largest fragment of each side, so those are the ends that meet.
    let largest = |f: &[Frag]| f.iter().max_by_key(|x| x.len).cloned();
    let fits = match (largest(&vector_cut), largest(&insert_cut)) {
        (Some(b), Some(i)) => Some(snappiesnap_engine::bio::cloning::how_it_fits(
            (&b.left, &b.right),
            (&i.left, &i.right),
        )),
        _ => None,
    };
    let measured = crate::cloning_plan::Measured {
        vector_fragments: vector_cut.iter().map(|f| f.len).collect(),
        insert_fragments: insert_cut.iter().map(|f| f.len).collect(),
        fits,
        ..measured
    };
    crate::cloning_plan::plan(
        exp,
        &measured,
        &names(&exp.vector.enzymes),
        &names(&exp.insert.enzymes),
    )
}

/// What a step works with besides the experiment itself: the session it reads constructs
/// from, the folder it offers files from, and how the construct is shown.
struct Step<'a> {
    session: &'a Session,
    folder: Option<&'a str>,
    pane: &'a mut Pane,
}

/// Step 2 and step 3 share their shape: pick a construct, prepare it, see its map.
fn construct_step(
    ui: &mut Ui,
    t: &Tokens,
    exp: &mut Experiment,
    inputs: &mut CloningInputs,
    ctx: &mut Step<'_>,
    is_vector: bool,
) -> Option<CloningAction> {
    let Step {
        session,
        folder,
        pane,
    } = ctx;
    let (session, folder) = (*session, *folder);
    let technique = exp.technique;
    steps(ui, t, exp, if is_vector { 2 } else { 3 });
    let files = folder_files(inputs, session, folder);
    // Open tabs count as constructs: in the web build the browser gives the page no file
    // paths at all, so what is open is the only thing there is to clone from.
    let offered = choices(session, &files);
    let label = if is_vector {
        technique.vector_label()
    } else {
        technique.insert_label()
    };
    ui.label(RichText::new(label).strong());
    // A column for the choices, the rest of the window for the construct: the map is what
    // makes a choice of enzymes obvious, so it gets the room (as the gel experiment does).
    let full = ui.available_width();
    // The columns take what is left after the row of step buttons: filling the whole height
    // pushed Back and Continue off the bottom of a tall window.
    let height = (ui.available_height() - NAV_H).max(320.0);
    // A third of the window for the choices, but never so narrow that the pickers can't be
    // read; the rest for the construct, which is what makes a choice of enzymes obvious.
    let controls_w = (full / 3.0)
        .clamp(280.0, 460.0)
        .min(full - MIN_MAP_W)
        .max(240.0);
    let map_w = (full - controls_w - t.spacing * 2.0).max(MIN_MAP_W);
    let mut pending: Option<CloningAction> = None;
    let mut request = None;
    // What the map was asked for: an enzyme or primer to put in a field, or a range to
    // select (a feature click, which selects here as it does in a document tab).
    let mut picked: Option<crate::map::MapAction> = None;
    ui.horizontal_top(|ui| {
        // `allocate_ui_with_layout`, not `vertical` + `set_width`: a plain vertical lets a
        // wide child push the next column off the window, which is what it did.
        ui.allocate_ui_with_layout(
            egui::vec2(controls_w, height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                ui.set_width(controls_w);
                // Mutagenesis has no second construct: step 3 takes the two back-to-back primers.
                if !is_vector && technique.insert_is_primers() {
                    ui.label(
            RichText::new(
                "Two primers that point away from each other and meet back to back around the \
                 plasmid. Changes in their 5' parts become the mutation, a gap between them a \
                 deletion, 5' tails an insertion.",
            )
            .color(t.text_muted),
        );
                    ui.add_space(t.spacing);
                    let before = exp.insert.clone();
                    primer_fields(ui, &mut exp.insert, "Forward primer", "Reverse primer");
                    if before != exp.insert {
                        inputs.forget_preview();
                    }
                } else {
                    // What this construct can be cut with, from the lab's set when one is chosen.
                    {
                        let side = if is_vector { &exp.vector } else { &exp.insert };
                        refresh_document(inputs, session, side);
                    }
                    let chosen = if is_vector {
                        exp.vector.source_key()
                    } else {
                        exp.insert.source_key()
                    };
                    let found = cutters(inputs, session, &chosen);
                    let pool = inputs.pool.clone();
                    let offer = Offer {
                        pool: pool.as_deref(),
                        db: session.enzymes(),
                        cutters: &found,
                    };
                    let changed = {
                        let side = if is_vector {
                            &mut exp.vector
                        } else {
                            &mut exp.insert
                        };
                        side_editor(
                            ui,
                            t,
                            if is_vector {
                                "clone-vector"
                            } else {
                                "clone-insert"
                            },
                            side,
                            &Side {
                                choices: &offered,
                                technique,
                                is_vector,
                                offer: &offer,
                            },
                        )
                    };
                    if changed {
                        inputs.forget_preview();
                    }
                    let side = if is_vector { &exp.vector } else { &exp.insert };
                    let key = side.source_key();
                    if side.is_set() {
                        refresh_document(inputs, session, side);
                        // Only a file has to be read; an open tab is already in hand.
                        if side.document.is_none() && !inputs.info.contains_key(&key) {
                            pending = Some(CloningAction::Info(key));
                        }
                    }
                }
                // Assemblies may have more than two parts: the rest are added here, in order.
                if !is_vector && technique.is_assembly() {
                    ui.add_space(t.spacing);
                    ui.separator();
                    let mut remove = None;
                    let mut changed = false;
                    let pool = inputs.pool.clone();
                    for (i, side) in exp.extra.iter_mut().enumerate() {
                        ui.horizontal(|ui| {
                            ui.label(RichText::new(format!("Part {}", i + 3)).strong());
                            if ui.small_button("Remove").clicked() {
                                remove = Some(i);
                            }
                        });
                        refresh_document(inputs, session, side);
                        let found = cutters(inputs, session, &side.source_key());
                        let offer = Offer {
                            pool: pool.as_deref(),
                            db: session.enzymes(),
                            cutters: &found,
                        };
                        changed |= side_editor(
                            ui,
                            t,
                            &format!("clone-extra-{i}"),
                            side,
                            &Side {
                                choices: &offered,
                                technique,
                                is_vector: false,
                                offer: &offer,
                            },
                        );
                        ui.add_space(t.spacing);
                    }
                    if let Some(i) = remove
                        && i < exp.extra.len()
                    {
                        exp.extra.remove(i);
                        changed = true;
                    }
                    if exp.extra.len() < MAX_PARTS
                        && ui
                            .button("Add another part")
                            .on_hover_text("Assemble more than two pieces in one reaction")
                            .clicked()
                    {
                        exp.extra.push(CloningSide::default());
                        changed = true;
                    }
                    if changed {
                        inputs.forget_preview();
                    }
                    // The extra parts have no card, but their enzymes are offered from what they hold, so
                    // they are read too.
                    if let Some(key) = exp
                        .extra
                        .iter()
                        .filter(|s| s.is_set() && s.document.is_none())
                        .map(CloningSide::source_key)
                        .find(|k| !inputs.info.contains_key(k))
                    {
                        pending = Some(CloningAction::Info(key));
                    }
                }
                // In step 3 the reaction itself says whether the two fit.
                if !is_vector {
                    ui.add_space(t.spacing);
                    request = fit_note(ui, t, exp, inputs);
                }
            },
        );
        // The construct itself, with the enzymes that cut it, and the plan underneath.
        ui.allocate_ui_with_layout(
            egui::vec2(map_w, height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                let w = map_w;
                ui.set_width(w);
                // Scrolls when the window is too short for the map and the plan together,
                // so neither is ever cut off.
                egui::ScrollArea::vertical()
                    .id_salt("cloning-construct")
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        // Step 3 of mutagenesis has no construct of its own, so it keeps
                        // showing the vector: the primers are what is being prepared.
                        let side = if is_vector || technique.insert_is_primers() {
                            &exp.vector
                        } else {
                            &exp.insert
                        };
                        let key = side.source_key();
                        if side.is_set() {
                            let typed = names(&side.enzymes);
                            let sites = sites_for(inputs, session, &key, &typed);
                            let molecule = inputs
                                .info
                                .get(&key)
                                .and_then(|i| i.as_ref().ok())
                                .and_then(|i| i.molecule.clone());
                            view_controls(ui, t, pane, inputs);
                            // Leave the plan its room; the column scrolls when the window
                            // is too short for both.
                            let height =
                                (ui.available_height() - PLAN_H - 24.0).clamp(160.0, 420.0);
                            let (orfs, primers) = match &molecule {
                                Some(m) => extras_for(inputs, session, &key, m),
                                None => (Vec::new(), Vec::new()),
                            };
                            // The layer toggles down the left edge, then the construct:
                            // the shape a document tab has.
                            ui.horizontal_top(|ui| {
                                layer_column(ui, t, pane);
                                // The construct's own column: its line of detail above the
                                // canvas, not beside it.
                                ui.vertical(|ui| {
                                    picked = construct_card(
                                        ui,
                                        t,
                                        inputs,
                                        &key,
                                        &Shown {
                                            size: egui::vec2(
                                                (w - LAYER_W).max(MIN_MAP_W - LAYER_W),
                                                height,
                                            ),
                                            sites: &sites,
                                            orfs: &orfs,
                                            primers: &primers,
                                            session,
                                            pane,
                                        },
                                    );
                                });
                            });
                            ui.label(
                                RichText::new(
                                    "Read-only: double-click an enzyme to cut with it, or a \
                                     primer to use it; click a feature to select it.",
                                )
                                .small()
                                .color(t.text_muted),
                            );
                        }
                        ui.add_space(t.spacing);
                        ui.separator();
                        ui.label(RichText::new("The plan").strong());
                        let plan = plan_of(exp, inputs, session);
                        crate::cloning_plan::show(ui, t, &plan, w);
                        crate::cloning_plan::legend(ui, t, &plan);
                    });
            },
        );
    });
    // Double-clicking an enzyme or a primer on the map puts it in a field of this side,
    // which is the quick way to pick from what the construct actually offers.
    // A pick made in the sequence view is the same thing as one made on the map.
    if let Some(i) = inputs.chosen_in_sequence.take() {
        picked = Some(crate::map::MapAction::SiteChosen(i));
    }
    if let Some(action) = picked {
        let key = if is_vector {
            exp.vector.source_key()
        } else {
            exp.insert.source_key()
        };
        let side = if is_vector {
            &mut exp.vector
        } else {
            &mut exp.insert
        };
        match action {
            crate::map::MapAction::SiteChosen(i) => {
                if add_site_enzyme(inputs, side, i) {
                    inputs.forget_preview();
                }
            }
            crate::map::MapAction::PrimerChosen(i) => {
                if add_primer(inputs, side, &key, i) {
                    inputs.forget_preview();
                }
            }
            // A feature or the backbone: selected, as in a document tab.
            crate::map::MapAction::Select(range) => inputs.select_in(&key, range),
        }
    }
    ui.add_space(t.spacing);
    let ready = if is_vector {
        exp.vector.is_set()
    } else {
        call(exp, true).is_some()
    };
    let action = request
        .map(|(c, p)| CloningAction::Preview(c, p))
        .or(pending);
    ui.horizontal(|ui| {
        if ui.button("Back").clicked() {
            exp.step = if is_vector {
                ExperimentStep::Technique
            } else {
                ExperimentStep::Vector
            };
        }
        let next = if is_vector {
            format!("Continue: {}", technique.insert_label())
        } else {
            "Continue: the product".to_string()
        };
        if ui.add_enabled(ready, egui::Button::new(next)).clicked() {
            exp.step = if is_vector {
                ExperimentStep::Insert
            } else {
                ExperimentStep::Product
            };
        }
        if !ready {
            ui.label(
                RichText::new(if is_vector {
                    "Choose a construct to go on."
                } else {
                    "Fill in what the reaction needs to go on."
                })
                .color(t.text_muted),
            );
        }
    });
    action
}

/// Whether the two sides fit, answered by running the reaction as a preview. Returns the
/// call to make when the answer isn't in yet.
fn fit_note(
    ui: &mut Ui,
    t: &Tokens,
    exp: &Experiment,
    inputs: &mut CloningInputs,
) -> Option<(&'static str, Value)> {
    let (command, params) = call(exp, true)?;
    let key = format!("{command}{params}");
    if inputs.preview_key.as_deref() != Some(key.as_str()) {
        inputs.preview_key = Some(key);
        inputs.preview = None;
        return Some((command, params));
    }
    match &inputs.preview {
        None => {
            ui.label(RichText::new("Checking that they fit\u{2026}").color(t.text_muted));
        }
        Some(Ok(p)) => {
            ui.label(
                RichText::new(format!(
                    "These fit: the product is {} bp \u{2014} {}",
                    p.length, p.description
                ))
                .color(t.accent_text),
            );
        }
        Some(Err(e)) => {
            ui.label(RichText::new(format!("These don't fit yet: {e}")).color(t.error));
        }
    }
    None
}

/// Step 4: the product, and opening it as a construct of its own.
fn product_step(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut CloningInputs,
    pane: &mut Pane,
) -> Option<CloningAction> {
    steps(ui, t, exp, 4);
    let Some((command, params)) = call(exp, true) else {
        ui.label(
            RichText::new("Something the reaction needs is still missing; go back a step.")
                .color(t.text_muted),
        );
        if ui.button("Back").clicked() {
            exp.step = ExperimentStep::Insert;
        }
        return None;
    };
    let key = format!("{command}{params}");
    if inputs.preview_key.as_deref() != Some(key.as_str()) {
        inputs.preview_key = Some(key);
        inputs.preview = None;
        return Some(CloningAction::Preview(command, params));
    }
    let mut action = None;
    // The schematic is worked out before the preview is read, because the map has to leave
    // room for it (without that the step's buttons land off the bottom) and because reading
    // the preview borrows `inputs` for the rest of the step.
    let plan = plan_of(exp, inputs, session);
    let plan_h = crate::cloning_plan::height(&plan) + 70.0;
    match &inputs.preview {
        None => {
            ui.label(RichText::new("Running the reaction\u{2026}").color(t.text_muted));
        }
        Some(Err(e)) => {
            ui.label(RichText::new(format!("This reaction doesn't work: {e}")).color(t.error));
        }
        Some(Ok(p)) => {
            ui.label(
                RichText::new(format!(
                    "{}: {} bp \u{b7} {}",
                    exp.technique.label(),
                    p.length,
                    p.description
                ))
                .strong(),
            );
            ui.label(
                RichText::new(format!(
                    "{} and {}",
                    side_label(session, &exp.vector),
                    if exp.technique.insert_is_primers() {
                        "its primers".to_string()
                    } else {
                        side_label(session, &exp.insert)
                    }
                ))
                .color(t.text_muted),
            );
            if let Some(m) = &p.molecule {
                ui.add_space(t.spacing);
                // The product gets the room and the choices the other steps give a
                // construct: it is the thing you are about to open, after all.
                let w = ui.available_width().max(MIN_MAP_W);
                let h = (ui.available_height() - 80.0 - plan_h).clamp(220.0, 460.0);
                let m = m.clone();
                let key = format!("product\u{1}{}", m.len());
                let sites = product_sites(inputs, session, &key, &m);
                view_controls(ui, t, pane, inputs);
                let (orfs, primers) = extras_for(inputs, session, &key, &m);
                // The toggles down the left edge, as on the other steps and in a tab.
                ui.horizontal_top(|ui| {
                    layer_column(ui, t, pane);
                    ui.vertical(|ui| {
                        let it = Shown {
                            size: egui::vec2((w - LAYER_W).max(MIN_MAP_W - LAYER_W), h),
                            sites: &sites,
                            orfs: &orfs,
                            primers: &primers,
                            session,
                            pane,
                        };
                        if pane.sequence {
                            sequence(ui, t, &m, &key, inputs, &it);
                        } else {
                            let selection = inputs.selection_in(&key);
                            let (mut zoom, mut pan) = (inputs.zoom, inputs.pan);
                            let picked = map(ui, t, &m, &it, selection, &mut zoom, &mut pan);
                            inputs.zoom = zoom;
                            inputs.pan = pan;
                            // The product is finished: a click selects, nothing fills a field.
                            if let Some(crate::map::MapAction::Select(range)) = picked {
                                inputs.select_in(&key, range);
                            }
                        }
                    });
                });
            }
            // The schematic here too, so the same picture of the reaction runs through
            // every step and the product can be read against what went in.
            ui.add_space(t.spacing);
            ui.separator();
            ui.label(RichText::new("The plan").strong());
            // Narrower than the map: rows spread over a whole wide window are harder to
            // compare than rows of a readable width.
            let w = ui
                .available_width()
                .clamp(crate::cloning_plan::MIN_WIDTH, PLAN_W);
            crate::cloning_plan::show(ui, t, &plan, w);
            crate::cloning_plan::legend(ui, t, &plan);
            ui.add_space(t.spacing);
            if ui
                .button("Open as a new construct")
                .on_hover_text(
                    "Opens the product in a tab of its own, with its construction history",
                )
                .clicked()
                && let Some((c, p)) = call(exp, false)
            {
                action = Some(CloningAction::Open(c, p));
            }
        }
    }
    ui.add_space(t.spacing);
    if ui.button("Back").clicked() {
        exp.step = ExperimentStep::Insert;
    }
    action
}

/// The cloning experiment, step by step.
pub fn show(
    ui: &mut Ui,
    t: &Tokens,
    session: &Session,
    exp: &mut Experiment,
    inputs: &mut CloningInputs,
    folder: Option<&str>,
    pane: &mut Pane,
) -> Option<CloningAction> {
    let mut ctx = Step {
        session,
        folder,
        pane,
    };
    match exp.step {
        ExperimentStep::Vector => construct_step(ui, t, exp, inputs, &mut ctx, true),
        ExperimentStep::Insert => construct_step(ui, t, exp, inputs, &mut ctx, false),
        ExperimentStep::Product => product_step(ui, t, session, exp, inputs, ctx.pane),
        _ => technique(ui, t, exp, inputs),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::{ExperimentKind, ExperimentStep};

    fn experiment(technique: CloningTechnique) -> Experiment {
        let mut e = Experiment::new(ExperimentKind::Cloning, "Cloning", []);
        e.technique = technique;
        e.vector.file = "vector.gb".into();
        e.insert.file = "insert.fa".into();
        e
    }

    /// A cloning experiment opens on its first step and builds no params until it has what
    /// the reaction needs.
    #[test]
    fn nothing_runs_until_the_reaction_is_complete() {
        let mut e = Experiment::new(ExperimentKind::Cloning, "Cloning", []);
        assert_eq!(e.step, ExperimentStep::Technique);
        // No vector yet.
        assert!(call(&e, true).is_none());
        e.vector.file = "vector.gb".into();
        // Restriction needs the vector's enzymes and an insert.
        assert!(call(&e, true).is_none());
        e.insert.file = "insert.fa".into();
        assert!(call(&e, true).is_none());
        e.vector.enzymes = "HindIII".into();
        assert!(call(&e, true).is_some());
        // Mutagenesis needs no second construct, but it does need both primers.
        let mut m = Experiment::new(ExperimentKind::Cloning, "Cloning", []);
        m.technique = CloningTechnique::Mutagenesis;
        m.vector.file = "plasmid.dna".into();
        assert!(call(&m, true).is_none());
        m.insert.forward = "ACGT".into();
        assert!(call(&m, true).is_none());
        m.insert.reverse = "TTTT".into();
        let (command, params) = call(&m, true).expect("ready");
        assert_eq!(command, "mutagenesis");
        assert_eq!(params["template"], "plasmid.dna");
        assert_eq!(params["preview"], true);
    }

    /// Each technique builds the params of the command an agent would call, with both
    /// constructs given as files: an experiment is bound to no open document.
    #[test]
    fn every_technique_builds_its_command() {
        let mut restriction = experiment(CloningTechnique::Restriction);
        restriction.vector.enzymes = "EcoRI, BamHI".into();
        restriction.insert.enzymes = "EcoRI".into();
        let (c, p) = call(&restriction, true).expect("ready");
        assert_eq!(c, "restriction_clone");
        assert_eq!(p["vector"], "vector.gb");
        assert_eq!(p["insert"], "insert.fa");
        assert_eq!(p["vector_enzymes"], json!(["EcoRI", "BamHI"]));
        assert_eq!(p["insert_enzymes"], json!(["EcoRI"]));

        let lr = experiment(CloningTechnique::GatewayLr);
        let (c, p) = call(&lr, true).expect("ready");
        assert_eq!((c, &p["reaction"]), ("gateway_clone", &json!("lr")));
        assert_eq!(p["acceptor"], "vector.gb");
        assert_eq!(p["donor"], "insert.fa");
        let bp = experiment(CloningTechnique::GatewayBp);
        assert_eq!(call(&bp, true).expect("ready").1["reaction"], "bp");

        // An assembly takes its parts in order, with how each is prepared.
        let mut gibson = experiment(CloningTechnique::Gibson);
        gibson.vector.using = SideUse::Cut;
        gibson.vector.enzymes = "EcoRI".into();
        gibson.insert.using = SideUse::Pcr;
        gibson.insert.forward = "f".into();
        gibson.insert.reverse = "r".into();
        gibson.extra.push(CloningSide {
            file: "third.fa".into(),
            ..CloningSide::default()
        });
        let (c, p) = call(&gibson, true).expect("ready");
        assert_eq!((c, &p["method"]), ("assemble", &json!("gibson")));
        let parts = p["parts"].as_array().expect("parts");
        assert_eq!(parts.len(), 3);
        assert_eq!(parts[0]["enzymes"], json!(["EcoRI"]));
        assert_eq!(parts[1]["pcr"], json!({"forward": "f", "reverse": "r"}));
        assert_eq!(parts[2]["file"], "third.fa");
        // A part left empty is left out.
        gibson.extra.push(CloningSide::default());
        assert_eq!(
            call(&gibson, true).expect("ready").1["parts"]
                .as_array()
                .map(Vec::len),
            Some(3)
        );

        // Golden Gate names its Type IIS enzyme, BsaI unless another is given.
        let mut gg = experiment(CloningTechnique::GoldenGate);
        assert_eq!(call(&gg, true).expect("ready").1["enzyme"], "BsaI");
        gg.vector.enzymes = "BsmBI".into();
        assert_eq!(call(&gg, true).expect("ready").1["enzyme"], "BsmBI");

        let mut topo = experiment(CloningTechnique::TopoDirectional);
        topo.vector.position = "10".into();
        topo.insert.using = SideUse::Pcr;
        topo.insert.forward = "f".into();
        topo.insert.reverse = "r".into();
        let (c, p) = call(&topo, true).expect("ready");
        assert_eq!((c, &p["kind"]), ("topo_clone", &json!("directional")));
        assert_eq!(p["position"], 10);
        assert_eq!(p["pcr"]["forward"], "f");
        // A position that isn't a number is simply left out, not a crash.
        topo.vector.position = "soon".into();
        assert!(
            call(&topo, true)
                .expect("ready")
                .1
                .get("position")
                .is_none()
        );
    }

    /// A construct and a session whose database holds `set` as the user's own enzyme set.
    fn lab(set: &[&str]) -> (Session, CloningInputs) {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        // One EcoRI site (GAATTC), two BamHI (GGATCC), one HindIII (AAGCTT).
        let residues = format!(
            "{0}GAATTC{0}GGATCC{0}AAGCTT{0}GGATCC{0}",
            "ACGTTACGATCGTAGC".repeat(20)
        );
        let m = Molecule::new(
            "pTest",
            Seq::new(Alphabet::Dna, residues).expect("residues"),
            Topology::Circular,
        )
        .expect("molecule");
        // A memory platform: the test never touches the user's own data folder.
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        if !set.is_empty() {
            session
                .set_enzyme_sets(vec![snappiesnap_engine::enzymes::UserSet {
                    name: "Fridge".into(),
                    enzymes: set.iter().map(|s| (*s).to_string()).collect(),
                }])
                .expect("sets kept");
        }
        let mut inputs = CloningInputs::default();
        inputs.info.insert(
            "v.gb".into(),
            Ok(ConstructInfo {
                name: "pTest".into(),
                length: m.len() as u64,
                topology: "circular".into(),
                features: 0,
                molecule: Some(m.clone()),
            }),
        );
        (session, inputs)
    }

    /// The schematic's pieces come from the digest itself, so what it draws is what the
    /// reaction will really produce.
    #[test]
    fn the_schematic_gets_the_real_fragments() {
        let (session, mut inputs) = lab(&[]);
        // One EcoRI site in a circle: one fragment, the plasmid opened (overhang at both
        // ends, so it counts four bases more than the plasmid's length).
        let one = fragments_for(&mut inputs, &session, "v.gb", &["EcoRI".into()]);
        assert_eq!(one.len(), 1, "a single cut opens the circle: {one:?}");
        // The plasmid is 1624 bp; an opened circle carries EcoRI's 4-base overhang at both
        // ends, so the fragment counts four bases more.
        assert_eq!(one[0].len, 1628);
        // Both ends are EcoRI's, so an opened plasmid closes on itself again.
        assert_eq!(one[0].left, one[0].right);
        // EcoRI and the two BamHI sites: three pieces, and together they cover the plasmid.
        let three = fragments_for(
            &mut inputs,
            &session,
            "v.gb",
            &["EcoRI".into(), "BamHI".into()],
        );
        assert_eq!(three.len(), 3, "{three:?}");
        // The cached answer follows the enzymes, not the construct alone.
        assert_eq!(
            fragments_for(&mut inputs, &session, "v.gb", &["EcoRI".into()]),
            one,
            "a different enzyme list is digested again, not served from the cache"
        );
        // Nothing to say without enzymes, for an unread construct, or for a name the
        // database doesn't know.
        assert!(fragments_for(&mut inputs, &session, "v.gb", &[]).is_empty());
        assert!(fragments_for(&mut inputs, &session, "unread.gb", &["EcoRI".into()]).is_empty());
        assert!(fragments_for(&mut inputs, &session, "v.gb", &["NotAnEnzyme".into()]).is_empty());
    }

    /// The schematic warns when the piece has the same cut at both ends: it then goes into
    /// the backbone either way round, and half the colonies carry it backwards. With two
    /// different enzymes the direction is settled and nothing is said.
    #[test]
    fn the_plan_warns_when_the_direction_is_open() {
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        let (session, mut inputs) = lab(&[]);
        // A vector with one EcoRI and one BamHI site, a few bases apart: cut with both, its
        // backbone carries one of each, which is what makes the clone directional. (lab()'s
        // own plasmid has two BamHI sites, so its largest fragment is BamHI at both ends and
        // an EcoRI-BamHI piece genuinely doesn't fit it.)
        let vector = Molecule::new(
            "pDir",
            Seq::new(
                Alphabet::Dna,
                format!(
                    "GAATTC{}GGATCC{}",
                    "ACGTACGTAC".repeat(2),
                    "ACGTACGTAC".repeat(60)
                ),
            )
            .expect("residues"),
            Topology::Circular,
        )
        .expect("molecule");
        inputs.info.insert(
            "pDir.gb".into(),
            Ok(ConstructInfo {
                name: "pDir".into(),
                length: vector.len() as u64,
                topology: "circular".into(),
                features: 0,
                molecule: Some(vector),
            }),
        );
        // An insert with an EcoRI site at each end of its middle, and a BamHI site inside it.
        let residues = format!(
            "TTTT{0}GAATTC{1}GGATCC{1}GAATTC{0}TTTT",
            "ACGTAGCT".repeat(4),
            "GTCAGTCAGTCA".repeat(6)
        );
        let insert = Molecule::new(
            "demo",
            Seq::new(Alphabet::Dna, residues).expect("residues"),
            Topology::Linear,
        )
        .expect("molecule");
        inputs.info.insert(
            "i.fa".into(),
            Ok(ConstructInfo {
                name: "demo".into(),
                length: insert.len() as u64,
                topology: "linear".into(),
                features: 0,
                molecule: Some(insert),
            }),
        );
        let mut exp = Experiment::new(ExperimentKind::Cloning, "Cloning", []);
        exp.technique = CloningTechnique::Restriction;
        exp.step = ExperimentStep::Insert;
        exp.vector = CloningSide {
            file: "pDir.gb".into(),
            enzymes: "EcoRI".into(),
            ..Default::default()
        };
        exp.insert = CloningSide {
            file: "i.fa".into(),
            ..Default::default()
        };
        // EcoRI alone: the same cut at both ends of the piece and of the opened vector.
        let p = plan_of(&exp, &mut inputs, &session);
        assert!(
            p.warnings.iter().any(|w| w.contains("either way round")),
            "{:?}",
            p.warnings
        );
        // EcoRI and BamHI: the ends differ, so the way round is settled.
        exp.vector.enzymes = "EcoRI, BamHI".into();
        let p = plan_of(&exp, &mut inputs, &session);
        assert!(p.warnings.is_empty(), "{:?}", p.warnings);
        // The vector's own step is about the backbone, so it says nothing about the piece.
        exp.vector.enzymes = "EcoRI".into();
        exp.step = ExperimentStep::Vector;
        assert!(plan_of(&exp, &mut inputs, &session).warnings.is_empty());
    }

    /// The enzyme fields offer the enzymes that cut *this* construct, fewest cuts first; with
    /// an enzyme set chosen, only the ones the lab holds.
    #[test]
    fn the_fields_offer_what_cuts_this_construct() {
        let (session, mut inputs) = lab(&[]);
        let all = cutters(&mut inputs, &session, "v.gb");
        let names: Vec<&str> = all.iter().map(|c| c.name.as_str()).collect();
        assert!(names.contains(&"EcoRI"), "{names:?}");
        assert!(names.contains(&"BamHI"), "{names:?}");
        assert!(names.contains(&"HindIII"), "{names:?}");
        // Fewest cuts first: the single cutters come before BamHI's two.
        let bam = all.iter().position(|c| c.name == "BamHI").expect("BamHI");
        let eco = all.iter().position(|c| c.name == "EcoRI").expect("EcoRI");
        assert!(eco < bam, "{names:?}");
        assert_eq!(all[bam].cuts, 2);
        assert_eq!(all[eco].cuts, 1);
        // A construct that hasn't been read yet offers nothing rather than guessing.
        assert!(cutters(&mut inputs, &session, "unread.gb").is_empty());

        // With the lab's set chosen, only what the lab holds is offered, and only if it cuts.
        let (session, mut inputs) = lab(&["BamHI", "SpeI"]);
        inputs.pool = Some("Fridge".into());
        let held = cutters(&mut inputs, &session, "v.gb");
        assert_eq!(
            held.iter().map(|c| c.name.as_str()).collect::<Vec<_>>(),
            vec!["BamHI"]
        );
        // A set that was deleted since offers nothing, never the whole database.
        inputs.pool = Some("Gone".into());
        assert!(cutters(&mut inputs, &session, "v.gb").is_empty());
    }

    /// What is typed is never restricted, but a name the database or the lab hasn't got is
    /// said so before the reaction is run.
    #[test]
    fn typing_an_enzyme_the_lab_lacks_is_told() {
        let (session, _) = lab(&["BamHI", "SpeI"]);
        let db = session.enzymes();
        let offer = |pool| Offer {
            pool,
            db,
            cutters: &[],
        };
        // Nothing typed, nothing said.
        assert!(enzyme_note("", &offer(Some("Fridge"))).is_none());
        // In the set: nothing to say.
        assert!(enzyme_note("BamHI, SpeI", &offer(Some("Fridge"))).is_none());
        // Case and separators don't matter.
        assert!(enzyme_note("bamhi + spei", &offer(Some("Fridge"))).is_none());
        // Known to the database, but not in the lab's set: a warning, naming it.
        let (message, warn) = enzyme_note("BamHI, EcoRI", &offer(Some("Fridge"))).expect("note");
        assert!(warn);
        assert!(message.contains("EcoRI"), "{message}");
        assert!(!message.contains("BamHI"), "{message}");
        // Unknown to the database: an error, whether or not a set is chosen.
        for pool in [Some("Fridge"), None] {
            let (message, warn) = enzyme_note("NotAnEnzyme", &offer(pool)).expect("note");
            assert!(!warn);
            assert!(message.contains("NotAnEnzyme"), "{message}");
        }
        // With no set chosen, an enzyme the lab may not hold is not second-guessed.
        assert!(enzyme_note("EcoRI", &offer(None)).is_none());
    }

    /// Choosing from the list adds to what is typed instead of replacing it.
    #[test]
    fn choosing_an_enzyme_keeps_what_is_typed() {
        let mut text = String::new();
        add_name(&mut text, "EcoRI");
        assert_eq!(text, "EcoRI");
        add_name(&mut text, "BamHI");
        assert_eq!(text, "EcoRI, BamHI");
        // Already there (in any case): added once.
        add_name(&mut text, "ecori");
        assert_eq!(text, "EcoRI, BamHI");
        // A half-typed list doesn't grow a second separator.
        let mut text = "EcoRI, ".to_string();
        add_name(&mut text, "BamHI");
        assert_eq!(text, "EcoRI, BamHI");
        let mut text = "EcoRI + ".to_string();
        add_name(&mut text, "BamHI");
        assert_eq!(text, "EcoRI, BamHI");
    }

    /// The pane reads a construct: a selection belongs to the construct it was made in, and
    /// copying it gives the same text the Copy command gives.
    #[test]
    fn the_pane_selects_and_copies_but_does_not_edit() {
        use snappiesnap_engine::doc::{Alphabet, Seq};
        let mut inputs = CloningInputs::default();
        let range = Range::new(2, 8);
        inputs.select_in("doc:1", Some(range));
        assert_eq!(inputs.selection_in("doc:1"), Some(range));
        // Another construct is another molecule: the range would mean something else there.
        assert_eq!(inputs.selection_in("doc:2"), None);
        assert_eq!(inputs.selection_in("lab/vector.gb"), None);
        inputs.select_in("doc:1", None);
        assert_eq!(inputs.selection_in("doc:1"), None);

        // The text a copy puts on the clipboard comes from the engine, so the pane and the
        // Copy command cannot drift apart.
        let m = Molecule::new(
            "pTest",
            Seq::new(Alphabet::Dna, "AACCGGTTAA".to_string()).expect("residues"),
            Topology::Linear,
        )
        .expect("molecule");
        let (text, n) = snappiesnap_engine::copy_text(&m, Range::new(0, 4), "plain").expect("copy");
        assert_eq!((text.as_str(), n), ("AACC", 4));
        let (rc, _) = snappiesnap_engine::copy_text(&m, Range::new(0, 4), "reverse_complement")
            .expect("copy");
        assert_eq!(rc, "GGTT");
        let (fasta, _) =
            snappiesnap_engine::copy_text(&m, Range::new(0, 4), "fasta").expect("copy");
        assert!(fasta.starts_with(">pTest 1..4\n"), "{fasta}");
    }

    /// Double-clicking an enzyme on the map adds it to the cut list of the side being
    /// prepared, and naming it twice doesn't double it up.
    #[test]
    fn double_clicking_an_enzyme_adds_it_to_the_cut_list() {
        let (session, mut inputs) = lab(&[]);
        let sites = sites_for(&mut inputs, &session, "v.gb", &[]);
        let eco = sites
            .iter()
            .position(|s| s.detail.label == "EcoRI")
            .expect("EcoRI cuts the test construct");
        let mut side = CloningSide {
            file: "v.gb".into(),
            ..CloningSide::default()
        };
        assert!(add_site_enzyme(&inputs, &mut side, eco));
        assert_eq!(side.enzymes, "EcoRI");
        // Again: nothing to add, and the field is left exactly as it was.
        assert!(!add_site_enzyme(&inputs, &mut side, eco));
        assert_eq!(side.enzymes, "EcoRI");
        // A second enzyme joins the first, in the form the fields use.
        let bam = sites
            .iter()
            .position(|s| s.detail.label == "BamHI")
            .expect("BamHI cuts it twice");
        assert!(add_site_enzyme(&inputs, &mut side, bam));
        assert_eq!(side.enzymes, "EcoRI, BamHI");
        // An index that isn't there changes nothing rather than panicking.
        assert!(!add_site_enzyme(&inputs, &mut side, 9_999));
    }

    /// Double-clicking a primer puts it in the forward or reverse field by the strand it
    /// binds, and says the side is made by PCR, which is the only way a primer is used.
    #[test]
    fn double_clicking_a_primer_fills_the_right_field() {
        use crate::enzymes::PrimerMark;
        use snappiesnap_engine::doc::Range;
        let mut inputs = CloningInputs::default();
        let mark = |name: &str, forward: bool| PrimerMark {
            primer: 0,
            name: name.to_string(),
            range: Range::new(0, 20),
            forward,
            primary: true,
            annealed: 20,
            mismatches: 0,
            five_prime_only: false,
            seq: "ACGTACGTACGTACGTACGT".to_string(),
            lane: 0,
            window: None,
            tm_site: None,
            tm_full: None,
            hairpin_dg: None,
            dimer_dg: None,
            checks: Vec::new(),
        };
        inputs.extras.insert(
            "doc:1".to_string(),
            (
                Vec::new(),
                vec![mark("IS_2667", true), mark("IS_2669", false)],
            ),
        );
        let mut side = CloningSide {
            document: Some(1),
            ..CloningSide::default()
        };
        assert!(add_primer(&inputs, &mut side, "doc:1", 0));
        assert_eq!(side.forward, "IS_2667");
        assert_eq!(side.using, SideUse::Pcr, "a primer means a PCR product");
        // The one on the other strand goes in the other field, not over the first.
        assert!(add_primer(&inputs, &mut side, "doc:1", 1));
        assert_eq!(
            (side.forward.as_str(), side.reverse.as_str()),
            ("IS_2667", "IS_2669")
        );
        // The same primer again changes nothing.
        assert!(!add_primer(&inputs, &mut side, "doc:1", 1));
        // An index or a construct that isn't there is ignored, not a panic.
        assert!(!add_primer(&inputs, &mut side, "doc:1", 7));
        assert!(!add_primer(&inputs, &mut side, "doc:9", 0));
    }

    /// The construct's map shows the enzymes: the ones typed in the fields when there are
    /// any, otherwise a readable number of the cutters the user is most likely to own.
    #[test]
    fn the_map_labels_the_enzymes_the_choice_is_about() {
        let (session, mut inputs) = lab(&[]);
        // Every enzyme that cuts the construct is offered: the choice is made from what is
        // there, not from a shortlist.
        let shown = sites_for(&mut inputs, &session, "v.gb", &[]);
        assert!(shown.len() > 3, "{} sites", shown.len());
        let names: Vec<&str> = shown.iter().map(|s| s.detail.label.as_str()).collect();
        for expected in ["HindIII", "EcoRI", "BamHI"] {
            assert!(
                names.contains(&expected),
                "{expected} missing from {names:?}"
            );
        }
        assert!(shown.iter().all(|s| !s.chosen), "nothing is named yet");
        // In order along the molecule, so the map reads left to right.
        let positions: Vec<usize> = shown.iter().map(|s| s.mark.pos).collect();
        let sorted = {
            let mut p = positions.clone();
            p.sort_unstable();
            p
        };
        assert_eq!(positions, sorted);

        // Naming an enzyme marks it; the rest stay on the map, so another choice is one
        // double-click away.
        let typed = vec!["BamHI".to_string()];
        let marked = sites_for(&mut inputs, &session, "v.gb", &typed);
        assert_eq!(marked.len(), shown.len(), "the others are still shown");
        let bam: Vec<&crate::enzymes::MapSite> = marked.iter().filter(|s| s.chosen).collect();
        assert_eq!(bam.len(), 2, "both BamHI sites are marked");
        assert!(bam.iter().all(|s| s.detail.names.contains("BamHI")));
        // Case doesn't matter: what the user typed is what they meant.
        let lower = sites_for(&mut inputs, &session, "v.gb", &["bamhi".to_string()]);
        assert_eq!(lower.iter().filter(|s| s.chosen).count(), 2);
        // A name the database hasn't got marks nothing, and is not an error.
        let unknown = sites_for(&mut inputs, &session, "v.gb", &["Nonesuch".to_string()]);
        assert_eq!(unknown.len(), shown.len());
        assert!(unknown.iter().all(|s| !s.chosen));
    }

    /// Constructs that are open in tabs can be cloned: the picker offers them, the params
    /// name them by id, and (the web build's case) no file is needed anywhere.
    #[test]
    fn open_tabs_can_be_cloned_without_any_file() {
        use snappiesnap_engine::Document;
        use snappiesnap_engine::doc::{Alphabet, Seq, Topology};
        let mut session =
            Session::new().with_platform(std::sync::Arc::new(snappiesnap_engine::MemoryFs::new()));
        for (name, residues) in [("pVec", "GGGGGAAGCTTGGGGGGGGGG"), ("pIns", "AAGCTTCCGG")] {
            let m = Molecule::new(
                name,
                Seq::new(Alphabet::Dna, residues.to_string()).expect("residues"),
                Topology::Linear,
            )
            .expect("molecule");
            session.open_document(Document::new(m));
        }
        let open = session.documents();
        let (vector_id, insert_id) = (open[0].id, open[1].id);
        // The picker offers the open tabs first, then the folder's files.
        let offered = choices(&session, &["lab/other.gb".to_string()]);
        assert_eq!(
            offered,
            vec![
                Choice::Document {
                    id: vector_id,
                    name: "pVec".into(),
                    path: None
                },
                Choice::Document {
                    id: insert_id,
                    name: "pIns".into(),
                    path: None
                },
                Choice::File("lab/other.gb".into()),
            ]
        );

        let mut e = Experiment::new(ExperimentKind::Cloning, "Cloning", []);
        e.technique = CloningTechnique::Restriction;
        e.vector.document = Some(vector_id);
        e.insert.document = Some(insert_id);
        e.vector.enzymes = "HindIII".into();
        // Both sides are set although neither has a file: this is the web build's case.
        assert!(e.vector.is_set() && e.insert.is_set());
        let (command, p) = call(&e, true).expect("ready");
        assert_eq!(command, "restriction_clone");
        assert_eq!(p["vector_document"], json!(vector_id));
        assert_eq!(p["insert_document"], json!(insert_id));
        assert!(p.get("vector").is_none() && p.get("insert").is_none());

        // An open tab wins over the file kept beside it, so the construct as it stands in
        // its tab is what is cloned.
        e.vector.file = "lab/vector.gb".into();
        let (_, p) = call(&e, true).expect("ready");
        assert_eq!(p["vector_document"], json!(vector_id));
        assert!(p.get("vector").is_none());
        // Dropping the tab choice falls back to that file, which is what a reopened .sslab
        // experiment does (the tab id is never saved).
        e.vector.document = None;
        let (_, p) = call(&e, true).expect("ready");
        assert_eq!(p["vector"], json!("lab/vector.gb"));
        assert!(p.get("vector_document").is_none());

        // Assemblies name a part the same way.
        e.technique = CloningTechnique::Gibson;
        e.vector.document = Some(vector_id);
        let (command, p) = call(&e, true).expect("ready");
        assert_eq!(command, "assemble");
        assert_eq!(p["parts"][0]["document"], json!(vector_id));
        assert_eq!(p["parts"][1]["document"], json!(insert_id));

        // The tab's details are read from the session, with no file_info call to make.
        let mut inputs = CloningInputs::default();
        refresh_document(&mut inputs, &session, &e.vector);
        let key = e.vector.source_key();
        assert_eq!(key, format!("doc:{vector_id}"));
        let info = inputs.info.get(&key).expect("read").as_ref().expect("ok");
        assert_eq!((info.name.as_str(), info.length), ("pVec", 21));
        assert!(info.molecule.is_some());
        // A tab that has been closed says so instead of showing the construct it held.
        let closed = CloningSide {
            document: Some(99_999),
            ..CloningSide::default()
        };
        refresh_document(&mut inputs, &session, &closed);
        assert!(inputs.info.get("doc:99999").expect("read").is_err());
    }

    /// Running it for real opens the product in a tab of its own and touches nothing else.
    #[test]
    fn running_for_real_discards_no_work_of_the_user() {
        let mut e = experiment(CloningTechnique::GatewayLr);
        e.step = ExperimentStep::Product;
        let (_, preview) = call(&e, true).expect("ready");
        assert_eq!(preview["preview"], true);
        assert!(preview.get("discard_changes").is_none());
        let (_, real) = call(&e, false).expect("ready");
        assert!(real.get("preview").is_none());
        assert_eq!(real["discard_changes"], true);
    }
}
