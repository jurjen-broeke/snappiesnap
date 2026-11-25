//! The session: open document, undo/redo, log, and the guarded command runner.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::Arc;

use serde_json::Value;
use snappiesnap_doc::{Detached, Molecule, Range, Seq, SeqDelta};
use snappiesnap_io::Source;

use crate::enzymes::{EnzymeDb, MAX_SETS, REBASE_CACHE, SETS_FILE, UserSet};
use crate::log::{Level, Log, LogEntry};
use crate::platform::Platform;
use crate::registry::Run;
use crate::{CmdError, CommandSpec, Effect, ErrorCode, Registry};

/// Undo snapshots kept per document; older ones are dropped.
pub const MAX_UNDO: usize = 100;

/// An open document.
#[derive(Debug, Clone, PartialEq)]
pub struct Document {
    pub molecule: Molecule,
    /// Where it was opened from or last saved to.
    pub path: Option<String>,
    /// Changed since it was opened or saved.
    pub dirty: bool,
    /// The current selection (view state, not undoable).
    pub selection: Option<Range>,
    /// The file the document was read from, kept so saving can be faithful.
    pub source: Option<Source>,
}

impl Document {
    pub fn new(molecule: Molecule) -> Self {
        Document {
            molecule,
            path: None,
            dirty: false,
            selection: None,
            source: None,
        }
    }
}

/// One undoable (or redoable) edit: the molecule to go back to, kept without its residues plus
/// what turns the current residues into its residues. After an edit that is a few bases, not a
/// copy of the whole sequence.
#[derive(Debug, Clone)]
struct Step {
    target: Detached,
    delta: SeqDelta,
    label: String,
}

impl Step {
    /// A step back to `target` from a molecule whose residues are `base`.
    fn new(target: Molecule, base: &Seq, label: String) -> Step {
        let delta = base.delta_to(target.seq());
        let (target, _) = target.detach();
        Step {
            target,
            delta,
            label,
        }
    }

    /// The molecule this step goes back to, from the current residues `base`.
    fn restore(&self, base: &Seq) -> Result<Molecule, CmdError> {
        let seq = base.apply(&self.delta).map_err(|e| {
            CmdError::new(ErrorCode::Internal, format!("undo record doesn't fit: {e}"))
        })?;
        self.target.clone().attach(seq).map_err(|e| {
            CmdError::new(ErrorCode::Internal, format!("undo record doesn't fit: {e}"))
        })
    }
}

/// A document that is open but not active, with its undo history.
#[derive(Debug, Clone)]
struct Parked {
    doc: Document,
    undo: Vec<Step>,
    redo: Vec<Step>,
}

/// One open document as listed for tabs.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct OpenDocument {
    pub index: usize,
    /// Stays the same while the document is open (tab positions shift when others close).
    pub id: u64,
    pub name: String,
    pub path: Option<String>,
    pub dirty: bool,
    pub active: bool,
}

/// The edit history, for display: labels of undoable edits (oldest first) and of undone edits
/// that can be redone (next first).
#[derive(Debug, Clone, Default, PartialEq, Eq, serde::Serialize)]
pub struct History {
    pub undo: Vec<String>,
    pub redo: Vec<String>,
}

/// A successful mutating command, as a script step would replay it.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct JournalEntry {
    pub command: String,
    pub params: Value,
}

/// Cloning is cheap apart from the document and its undo history: the command registry is
/// shared.
#[derive(Clone)]
pub struct Session {
    registry: Arc<Registry>,
    platform: Arc<dyn Platform>,
    /// The active document and its history. Other open documents wait in `tabs`.
    doc: Option<Document>,
    undo: Vec<Step>,
    redo: Vec<Step>,
    /// Every open document in tab order; the active one's slot is empty (it lives in `doc`).
    /// Invariant: `doc` is set exactly when `tabs` is non-empty, and `tabs[active]` is `None`.
    tabs: Vec<Option<Parked>>,
    /// A stable id per tab, parallel to `tabs`.
    tab_ids: Vec<u64>,
    next_tab_id: u64,
    active: usize,
    /// Successful mutating commands, when journaling is on (for `--record`).
    journal: Option<Vec<JournalEntry>>,
    log: Log,
    enzymes: Arc<EnzymeDb>,
    /// Which primer binding sites count (Preferences; `primer_binding`).
    primer_binding: snappiesnap_bio::primers::BindingPrefs,
    /// Changes whenever the document or the enzymes may have changed (see [`Session::revision`]).
    revision: u64,
}

impl Default for Session {
    fn default() -> Self {
        Session::new()
    }
}

impl Session {
    /// A session with the standard commands.
    pub fn new() -> Self {
        Session::with_commands(crate::cmd::all())
    }

    /// A session with a custom set of commands (used by tests).
    pub fn with_commands(specs: Vec<CommandSpec>) -> Self {
        let (registry, problems) = Registry::new(specs);
        let mut session = Session {
            registry: Arc::new(registry),
            platform: crate::platform::default_platform(),
            doc: None,
            undo: Vec::new(),
            redo: Vec::new(),
            tabs: Vec::new(),
            tab_ids: Vec::new(),
            next_tab_id: 1,
            active: 0,
            journal: None,
            log: Log::default(),
            enzymes: Arc::new(EnzymeDb::default()),
            primer_binding: Default::default(),
            revision: 0,
        };
        for p in problems {
            session.log.push(Level::Error, "registry", p);
        }
        session
    }

    /// Use a different file system (for example [`crate::MemoryFs`] in tests).
    pub fn with_platform(mut self, platform: Arc<dyn Platform>) -> Self {
        self.platform = platform;
        self
    }

    pub fn platform(&self) -> &dyn Platform {
        self.platform.as_ref()
    }

    /// A shared handle to the platform (for UI code that holds it while editing the session).
    pub fn platform_arc(&self) -> Arc<dyn Platform> {
        Arc::clone(&self.platform)
    }

    pub fn registry(&self) -> &Registry {
        &self.registry
    }

    /// The enzymes used for sites and digests.
    pub fn enzymes(&self) -> &EnzymeDb {
        &self.enzymes
    }

    /// A shared handle to the enzyme database (cheap to hold across edits).
    pub fn enzymes_arc(&self) -> Arc<EnzymeDb> {
        Arc::clone(&self.enzymes)
    }

    /// Use another enzyme database. The user's enzyme sets stay.
    pub fn set_enzymes(&mut self, mut db: EnzymeDb) {
        db.sets = self.enzymes.sets.clone();
        self.enzymes = Arc::new(db);
        self.touch();
    }

    /// Replace the user's enzyme sets and keep them in the data folder. The new sets are in
    /// use even when saving fails, which the error reports. Returns where they were saved.
    pub fn set_enzyme_sets(&mut self, sets: Vec<UserSet>) -> Result<Option<String>, String> {
        let mut db = (*self.enzymes).clone();
        db.sets = sets;
        let json = serde_json::to_vec_pretty(&serde_json::json!({ "sets": db.sets }))
            .map_err(|e| e.to_string());
        self.enzymes = Arc::new(db);
        self.touch();
        let Some(path) = self.enzyme_sets_path() else {
            return Ok(None);
        };
        self.platform
            .write(&path, &json?)
            .map_err(|e| format!("enzyme sets not saved to {path}: {e}"))?;
        Ok(Some(path))
    }

    /// Path of the user's enzyme sets, if the platform has a data folder.
    pub fn enzyme_sets_path(&self) -> Option<String> {
        let dir = self.platform.data_dir()?;
        Some(format!("{}/{SETS_FILE}", dir.trim_end_matches(['/', '\\'])))
    }

    /// Load the user's enzyme sets from the data folder. A missing file means none; a broken
    /// one is logged and left alone. Returns how many sets were loaded.
    pub fn load_enzyme_sets(&mut self) -> usize {
        let Some(path) = self.enzyme_sets_path() else {
            return 0;
        };
        if !self.platform.exists(&path) {
            return 0;
        }
        #[derive(serde::Deserialize)]
        struct File {
            sets: Vec<UserSet>,
        }
        let loaded = self
            .platform
            .read(&path)
            .map_err(|e| e.to_string())
            .and_then(|bytes| serde_json::from_slice::<File>(&bytes).map_err(|e| e.to_string()));
        match loaded {
            Ok(file) => {
                // Checked again: the file may have been edited by hand.
                let sets: Vec<UserSet> = file
                    .sets
                    .iter()
                    .filter_map(|s| UserSet::new(&s.name, &s.enzymes).ok().map(|(s, _)| s))
                    .take(MAX_SETS)
                    .collect();
                let n = sets.len();
                let mut db = (*self.enzymes).clone();
                db.sets = sets;
                self.enzymes = Arc::new(db);
                self.touch();
                n
            }
            Err(e) => {
                self.log(
                    Level::Warn,
                    "enzymes",
                    format!("enzyme sets in {path} not used: {e}"),
                );
                0
            }
        }
    }

    /// Which primer binding sites count: the rules `primer_sites`, the oligo search, PCR's
    /// weaker sites and the map use.
    pub fn primer_binding(&self) -> &snappiesnap_bio::primers::BindingPrefs {
        &self.primer_binding
    }

    /// Change the binding rules (kept in range); views showing sites update.
    pub fn set_primer_binding(&mut self, prefs: snappiesnap_bio::primers::BindingPrefs) {
        let prefs = prefs.clamped();
        if prefs != self.primer_binding {
            self.primer_binding = prefs;
            self.touch();
        }
    }

    /// A number that changes whenever the document, its molecule, the enzyme database or the
    /// primer binding rules may have changed, so views can cache what they derive from them.
    pub fn revision(&self) -> u64 {
        self.revision
    }

    fn touch(&mut self) {
        self.revision = self.revision.wrapping_add(1);
    }

    /// Path of the cached REBASE download, if the platform has a data folder.
    pub fn enzyme_cache_path(&self) -> Option<String> {
        let dir = self.platform.data_dir()?;
        Some(format!(
            "{}/{REBASE_CACHE}",
            dir.trim_end_matches(['/', '\\'])
        ))
    }

    /// Use the cached REBASE download if there is one. A missing cache is normal; a broken one
    /// is logged and the built-in enzymes stay in use. Returns whether REBASE was loaded.
    pub fn load_enzyme_cache(&mut self) -> bool {
        let Some(path) = self.enzyme_cache_path() else {
            return false;
        };
        if !self.platform.exists(&path) {
            return false;
        }
        let loaded = self
            .platform
            .read(&path)
            .map_err(|e| e.to_string())
            .and_then(|bytes| EnzymeDb::from_rebase(&bytes));
        match loaded {
            Ok((db, _)) => {
                self.set_enzymes(db);
                true
            }
            Err(e) => {
                self.log(
                    Level::Warn,
                    "enzymes",
                    format!("cached REBASE file {path} not used: {e}"),
                );
                false
            }
        }
    }

    pub fn document(&self) -> Option<&Document> {
        self.doc.as_ref()
    }

    pub fn document_mut(&mut self) -> Option<&mut Document> {
        self.touch();
        self.doc.as_mut()
    }

    pub fn molecule(&self) -> Result<&Molecule, CmdError> {
        self.doc
            .as_ref()
            .map(|d| &d.molecule)
            .ok_or_else(CmdError::no_document)
    }

    pub fn molecule_mut(&mut self) -> Result<&mut Molecule, CmdError> {
        self.touch();
        self.doc
            .as_mut()
            .map(|d| &mut d.molecule)
            .ok_or_else(CmdError::no_document)
    }

    /// Open a document in a new tab and make it the active one.
    pub fn open_document(&mut self, doc: Document) {
        self.touch();
        self.park();
        self.tabs.push(None);
        self.tab_ids.push(self.next_tab_id);
        self.next_tab_id = self.next_tab_id.wrapping_add(1);
        self.active = self.tabs.len() - 1;
        self.doc = Some(doc);
        self.undo.clear();
        self.redo.clear();
    }

    /// Close the active document; its neighbour (if any) becomes active.
    pub fn close_document(&mut self) {
        self.touch();
        self.doc = None;
        self.undo.clear();
        self.redo.clear();
        if self.active < self.tabs.len() {
            self.tabs.remove(self.active);
        }
        if self.active < self.tab_ids.len() {
            self.tab_ids.remove(self.active);
        }
        if !self.tabs.is_empty() {
            self.active = self.active.min(self.tabs.len() - 1);
            self.unpark();
        } else {
            self.active = 0;
        }
    }

    /// Move the active document into its tab slot.
    fn park(&mut self) {
        if let Some(doc) = self.doc.take()
            && let Some(slot) = self.tabs.get_mut(self.active)
        {
            *slot = Some(Parked {
                doc,
                undo: std::mem::take(&mut self.undo),
                redo: std::mem::take(&mut self.redo),
            });
        }
    }

    /// Make the document in the active slot the active document.
    fn unpark(&mut self) {
        if let Some(p) = self.tabs.get_mut(self.active).and_then(Option::take) {
            self.doc = Some(p.doc);
            self.undo = p.undo;
            self.redo = p.redo;
        }
    }

    /// Make tab `index` the active document.
    pub fn switch_document(&mut self, index: usize) -> Result<(), CmdError> {
        if index >= self.tabs.len() {
            return Err(CmdError::invalid_params(format!(
                "no document {index} ({} open)",
                self.tabs.len()
            )));
        }
        if index != self.active {
            self.touch();
            self.park();
            self.active = index;
            self.unpark();
        }
        Ok(())
    }

    /// The open documents in tab order.
    pub fn documents(&self) -> Vec<OpenDocument> {
        self.tabs
            .iter()
            .enumerate()
            .filter_map(|(index, slot)| {
                let d = match slot {
                    Some(p) => &p.doc,
                    None => self.doc.as_ref()?,
                };
                Some(OpenDocument {
                    index,
                    id: self.tab_ids.get(index).copied().unwrap_or(0),
                    name: d.molecule.name().to_string(),
                    path: d.path.clone(),
                    dirty: d.dirty,
                    active: index == self.active && slot.is_none(),
                })
            })
            .collect()
    }

    /// Any open document, by its tab's stable id (see [`OpenDocument::id`]): to show it beside
    /// the active one. Read-only; edits go to the active document.
    pub fn document_by_id(&self, id: u64) -> Option<&Document> {
        let index = self.tab_ids.iter().position(|&t| t == id)?;
        match self.tabs.get(index)? {
            Some(p) => Some(&p.doc),
            None => self.doc.as_ref(),
        }
    }

    /// The stable id of the active document's tab (see [`OpenDocument::id`]).
    pub fn active_document_id(&self) -> Option<u64> {
        self.doc.as_ref()?;
        self.tab_ids.get(self.active).copied()
    }

    /// The tab holding the document read from `path`, if it is open.
    pub fn document_with_path(&self, path: &str) -> Option<usize> {
        self.documents()
            .into_iter()
            .find(|d| d.path.as_deref() == Some(path))
            .map(|d| d.index)
    }

    pub fn can_undo(&self) -> bool {
        self.doc.is_some() && !self.undo.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        self.doc.is_some() && !self.redo.is_empty()
    }

    /// Undo the last edit. Checks everything before changing anything.
    pub fn undo(&mut self) -> Result<(), CmdError> {
        self.touch();
        let doc = self.doc.as_mut().ok_or_else(CmdError::no_document)?;
        let step = self
            .undo
            .last()
            .ok_or_else(|| CmdError::new(ErrorCode::NotEnabled, "nothing to undo"))?;
        // Rebuild first: if the record doesn't fit, nothing changes.
        let previous = step.restore(doc.molecule.seq())?;
        let label = step.label.clone();
        self.undo.pop();
        let current = std::mem::replace(&mut doc.molecule, previous);
        self.redo
            .push(Step::new(current, doc.molecule.seq(), label));
        doc.dirty = true;
        doc.selection = None;
        Ok(())
    }

    /// Redo the last undone edit.
    pub fn redo(&mut self) -> Result<(), CmdError> {
        self.touch();
        let doc = self.doc.as_mut().ok_or_else(CmdError::no_document)?;
        let step = self
            .redo
            .last()
            .ok_or_else(|| CmdError::new(ErrorCode::NotEnabled, "nothing to redo"))?;
        let next = step.restore(doc.molecule.seq())?;
        let label = step.label.clone();
        self.redo.pop();
        let current = std::mem::replace(&mut doc.molecule, next);
        self.undo
            .push(Step::new(current, doc.molecule.seq(), label));
        doc.dirty = true;
        doc.selection = None;
        Ok(())
    }

    /// Labels of the edits that can be undone and redone.
    pub fn history(&self) -> History {
        History {
            undo: self.undo.iter().map(|s| s.label.clone()).collect(),
            redo: self.redo.iter().rev().map(|s| s.label.clone()).collect(),
        }
    }

    /// Start or stop recording successful mutating commands (see [`Session::take_journal`]).
    pub fn set_journaling(&mut self, on: bool) {
        self.journal = on.then(Vec::new);
    }

    /// The commands recorded since the last call, oldest first.
    pub fn take_journal(&mut self) -> Vec<JournalEntry> {
        self.journal
            .as_mut()
            .map(std::mem::take)
            .unwrap_or_default()
    }

    pub fn log(&mut self, level: Level, source: &str, message: impl Into<String>) {
        self.log.push(level, source, message);
    }

    /// The last `n` log entries, oldest first.
    pub fn log_tail(&self, n: usize) -> Vec<LogEntry> {
        self.log.tail(n)
    }

    /// Run a command by id with JSON params.
    ///
    /// Never panics. A failing or panicking command leaves the document exactly as it was;
    /// a successful edit that changed the molecule is recorded for undo.
    pub fn execute(&mut self, id: &str, params: Value) -> Result<Value, CmdError> {
        let result = self.execute_inner(id, params);
        if let Err(e) = &result {
            let level = if e.code == ErrorCode::Internal {
                Level::Error
            } else {
                Level::Warn
            };
            self.log.push(level, id, e.to_string());
        }
        result
    }

    fn execute_inner(&mut self, id: &str, params: Value) -> Result<Value, CmdError> {
        let registry = Arc::clone(&self.registry);
        let spec = registry.get(id).ok_or_else(|| {
            CmdError::new(ErrorCode::UnknownCommand, format!("unknown command {id:?}"))
                .with_hint("`commands.list` lists every command")
        })?;
        if spec.needs_document && self.doc.is_none() {
            return Err(CmdError::no_document());
        }
        match &spec.run {
            Run::Query(run) => {
                let this: &Session = self;
                catch_unwind(AssertUnwindSafe(|| run(this, params)))
                    .unwrap_or_else(|payload| Err(internal(id, payload.as_ref(), false)))
            }
            Run::Mutate(run) => {
                self.touch();
                let saved = self.doc.clone();
                let journal_params = self.journal.is_some().then(|| params.clone());
                match catch_unwind(AssertUnwindSafe(|| run(self, params))) {
                    Ok(Ok(value)) => {
                        if spec.effect == Effect::Edit {
                            self.record_edit(saved, spec.label);
                        }
                        if let (Some(journal), Some(params)) =
                            (self.journal.as_mut(), journal_params)
                        {
                            journal.push(JournalEntry {
                                command: id.to_string(),
                                params,
                            });
                        }
                        Ok(value)
                    }
                    Ok(Err(e)) => {
                        self.doc = saved;
                        Err(e)
                    }
                    Err(payload) => {
                        self.doc = saved;
                        Err(internal(id, payload.as_ref(), true))
                    }
                }
            }
        }
    }

    /// After a successful edit: if the molecule changed, push the previous version for undo.
    fn record_edit(&mut self, before: Option<Document>, label: &str) {
        let (Some(before), Some(doc)) = (before, self.doc.as_mut()) else {
            return;
        };
        if before.molecule == doc.molecule {
            return;
        }
        doc.dirty = true;
        let step = Step::new(before.molecule, doc.molecule.seq(), label.to_string());
        self.undo.push(step);
        if self.undo.len() > MAX_UNDO {
            self.undo.remove(0);
        }
        self.redo.clear();
    }
}

fn internal(id: &str, payload: &(dyn std::any::Any + Send), restored: bool) -> CmdError {
    let detail = payload
        .downcast_ref::<&str>()
        .map(|s| s.to_string())
        .or_else(|| payload.downcast_ref::<String>().cloned())
        .unwrap_or_else(|| "unknown panic".to_string());
    let mut message = format!("internal error in command {id:?}: {detail}");
    if restored {
        message.push_str("; the document was restored");
    }
    CmdError::new(ErrorCode::Internal, message)
        .with_hint("this is a bug in SnappieSnap; please report it with the steps that led here")
}

#[cfg(test)]
mod tests {
    use serde::Deserialize;
    use serde_json::json;
    use snappiesnap_doc::{Alphabet, Seq, Topology};

    use super::*;
    use crate::NoParams;

    #[derive(Deserialize, schemars::JsonSchema)]
    #[serde(deny_unknown_fields)]
    struct At {
        at: usize,
    }

    fn commands() -> Vec<CommandSpec> {
        vec![
            CommandSpec::edit("add_a", "Add A", "insert one A", |s, p: At| {
                s.molecule_mut()?.insert(p.at, b"A")?;
                Ok(json!({}))
            }),
            CommandSpec::edit(
                "half_then_fail",
                "",
                "edits, then fails",
                |s, _: NoParams| {
                    s.molecule_mut()?.insert(0, b"C")?;
                    Err(CmdError::new(ErrorCode::InvalidInput, "deliberate"))
                },
            ),
            CommandSpec::edit(
                "half_then_panic",
                "",
                "edits, then panics",
                |s, _: NoParams| {
                    s.molecule_mut()?.insert(0, b"G")?;
                    #[allow(clippy::panic)]
                    {
                        panic!("deliberate panic")
                    }
                },
            ),
            CommandSpec::query("peek_panic", "", "panics", |_, _: NoParams| {
                #[allow(clippy::panic)]
                {
                    panic!("deliberate")
                }
            }),
            CommandSpec::edit("noop", "", "changes nothing", |_, _: NoParams| {
                Ok(json!({}))
            }),
            CommandSpec::session("undo", "", "", |s, _: NoParams| {
                s.undo()?;
                Ok(json!({}))
            }),
            CommandSpec::session("redo", "", "", |s, _: NoParams| {
                s.redo()?;
                Ok(json!({}))
            }),
        ]
    }

    fn session() -> Session {
        let mut s = Session::with_commands(commands());
        let m = Molecule::new(
            "t",
            Seq::new(Alphabet::Dna, "TTTT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        s.open_document(Document::new(m));
        s
    }

    fn residues(s: &Session) -> String {
        s.molecule().unwrap().seq().as_str().to_string()
    }

    #[test]
    fn failing_and_panicking_commands_leave_the_document_untouched() {
        let mut s = session();
        let err = s.execute("half_then_fail", json!({})).unwrap_err();
        assert_eq!(err.code, ErrorCode::InvalidInput);
        assert_eq!(residues(&s), "TTTT");

        let err = s.execute("half_then_panic", Value::Null).unwrap_err();
        assert_eq!(err.code, ErrorCode::Internal);
        assert!(err.message.contains("deliberate panic"), "{}", err.message);
        assert_eq!(residues(&s), "TTTT");
        assert!(!s.can_undo());

        let err = s.execute("peek_panic", json!({})).unwrap_err();
        assert_eq!(err.code, ErrorCode::Internal);
        assert!(
            s.log_tail(10)
                .iter()
                .any(|e| e.level == Level::Error && e.source == "half_then_panic")
        );
    }

    #[test]
    fn undo_and_redo() {
        let mut s = session();
        s.execute("add_a", json!({"at": 0})).unwrap();
        s.execute("add_a", json!({"at": 5})).unwrap();
        s.execute("noop", json!({})).unwrap();
        assert_eq!(residues(&s), "ATTTTA");
        s.execute("undo", json!({})).unwrap();
        assert_eq!(residues(&s), "ATTTT");
        s.execute("undo", json!({})).unwrap();
        assert_eq!(residues(&s), "TTTT");
        assert_eq!(
            s.execute("undo", json!({})).unwrap_err().code,
            ErrorCode::NotEnabled
        );
        s.execute("redo", json!({})).unwrap();
        assert_eq!(residues(&s), "ATTTT");
        // A new edit clears the redo stack.
        s.execute("add_a", json!({"at": 0})).unwrap();
        assert!(!s.can_redo());
        assert!(s.document().unwrap().dirty);
    }

    #[test]
    fn history_labels_follow_undo_and_redo() {
        let mut s = session();
        s.execute("add_a", json!({"at": 0})).unwrap();
        s.execute("add_a", json!({"at": 0})).unwrap();
        assert_eq!(s.history().undo, vec!["Add A", "Add A"]);
        s.execute("undo", json!({})).unwrap();
        assert_eq!(
            s.history(),
            History {
                undo: vec!["Add A".into()],
                redo: vec!["Add A".into()]
            }
        );
    }

    #[test]
    fn journal_records_successful_mutations_only() {
        let mut s = session();
        assert!(s.take_journal().is_empty());
        s.set_journaling(true);
        s.execute("add_a", json!({"at": 0})).unwrap();
        let _ = s.execute("add_a", json!({"at": 99}));
        let _ = s.execute("half_then_panic", json!({}));
        s.execute("undo", json!({})).unwrap();
        let j = s.take_journal();
        assert_eq!(j.len(), 2);
        assert_eq!(
            j[0],
            JournalEntry {
                command: "add_a".into(),
                params: json!({"at": 0})
            }
        );
        assert_eq!(j[1].command, "undo");
        assert!(s.take_journal().is_empty());
    }

    #[test]
    fn undo_history_is_bounded() {
        let mut s = session();
        for _ in 0..MAX_UNDO + 10 {
            s.execute("add_a", json!({"at": 0})).unwrap();
        }
        let mut n = 0;
        while s.execute("undo", json!({})).is_ok() {
            n += 1;
        }
        assert_eq!(n, MAX_UNDO);
    }

    #[test]
    fn bad_calls_are_errors() {
        let mut s = session();
        assert_eq!(
            s.execute("nope", json!({})).unwrap_err().code,
            ErrorCode::UnknownCommand
        );
        for bad in [
            json!({}),
            json!({"at": -1}),
            json!({"at": "0"}),
            json!({"at": 0, "x": 1}),
            json!([]),
        ] {
            assert_eq!(
                s.execute("add_a", bad).unwrap_err().code,
                ErrorCode::InvalidParams
            );
        }
        assert_eq!(
            s.execute("add_a", json!({"at": 99})).unwrap_err().code,
            ErrorCode::InvalidInput
        );
        s.close_document();
        assert_eq!(
            s.execute("add_a", json!({"at": 0})).unwrap_err().code,
            ErrorCode::NoDocument
        );
    }

    #[test]
    fn duplicate_and_malformed_ids_are_reported() {
        let mut specs = commands();
        specs.push(CommandSpec::edit("noop", "", "", |_, _: NoParams| {
            Ok(json!({}))
        }));
        specs.push(CommandSpec::edit("Bad-Id", "", "", |_, _: NoParams| {
            Ok(json!({}))
        }));
        let s = Session::with_commands(specs);
        assert_eq!(s.registry().len(), commands().len());
        assert_eq!(s.log_tail(10).len(), 2);
    }

    #[test]
    fn schemas_describe_params() {
        let s = session();
        let schema = s.registry().get("add_a").unwrap().params_schema();
        assert_eq!(schema["properties"]["at"]["type"], "integer");
        assert_eq!(schema["required"], json!(["at"]));
    }
}
