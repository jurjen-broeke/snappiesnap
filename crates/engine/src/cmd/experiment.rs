//! Gel experiments: one gel across constructs, each lane a construct (an open document) and a
//! method (uncut, digest, PCR), run exactly as on that construct's own gel. Experiments are
//! saved to and opened from `.sslab` files ([`crate::experiment`]).

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use std::collections::BTreeMap;

use snappiesnap_doc::{Molecule, Topology};
use snappiesnap_io::Format;

use super::pcr::{MAX_LANES, bands_json, one};
use crate::experiment::{
    self as file, Construct, EXTENSION, ExperimentFile, FORMAT, GelSettings, Lane, VERSION,
    fingerprint, folder_of, relative_path, resolve,
};
use crate::gel::{ExperimentLane, experiment};
use crate::{CmdError, CommandSpec, Document, ErrorCode, Level, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "gel_experiment",
            "Simulate Experiment Gel",
            "One agarose gel across several constructs: the ladder in lane 1, then a lane per entry of `lanes` in order, each an open document (its `id` from `documents`) and a method (`uncut`, `digest`, `pcr`, as in `gel`; none: an empty lane). Bands as in `gel`; a lane that can't be run reports `error`.",
            gel_experiment,
        )
        .needs_document(false)
        .example(json!({"lanes": [{"document": 1, "method": {"kind": "uncut"}}]}))
        .example(json!({"agarose": 0.8, "ladder": "1 kb", "lanes": [
            {"document": 1, "method": {"kind": "digest", "enzymes": ["EcoRI"]}},
            {"document": 1, "method": {"kind": "pcr", "forward": "f", "reverse": "r"}},
            {"document": 1}
        ]})),
        CommandSpec::query(
            "diagnostic_digest",
            "Design Diagnostic Digest",
            "Find restriction enzymes, alone or in pairs, whose digests tell open constructs apart on an agarose gel (a clone from its empty vector, one orientation from the other), or for one construct give a pattern that confirms it. `documents`: open documents by id (1-12). Every construct must be cut into at most `max_bands` fragments; bands of one lane should run apart and every two constructs must differ by a clearly visible band, judged with the gel model at `agarose` %. Enzymes that cut the same way are listed once with `alternatives`. Results are best first, with each lane's fragments and notes (double digests need a shared buffer). `enzymes` limits the search to those names. `set` searches only an enzyme set of the user's own, so every suggestion is one the lab can run; `include_out_of_stock` widens that search and marks the suggestions that need an enzyme the set hasn't got. Default without either: every enzyme sold by a supplier.",
            diagnostic_digest,
        )
        .needs_document(false)
        .example(json!({"documents": [1]}))
        .example(json!({"documents": [1], "agarose": 0.8, "max_enzymes": 1, "limit": 5}))
        .example(json!({"documents": [1], "set": "Lab freezer"}))
        .example(json!({"documents": [1], "set": "Lab freezer", "include_out_of_stock": true})),
        CommandSpec::session(
            "experiment_save",
            "Save Experiment",
            "Save a gel experiment as a SnappieSnap experiment file (.sslab; the format is described in docs/experiment-format.md). `lanes` as in `gel_experiment` (open documents by id, each with a method or none). Each construct is stored as a reference to its file (relative to the experiment file, and as it was opened) with a fingerprint of its sequence, and as an embedded GenBank copy of what is open now, unsaved edits included.",
            save,
        )
        .example(json!({"path": "reads.sslab", "name": "Sequencing check", "kind": "reads", "reference": 1, "reads": ["read.ab1"]}))
        .example(json!({"path": "run.sslab", "name": "Digest check", "lanes": [
            {"document": 1, "method": {"kind": "digest", "enzymes": ["EcoRI"]}}
        ]})),
        CommandSpec::session(
            "experiment_open",
            "Open Experiment",
            "Open a SnappieSnap experiment file (.sslab). Each construct opens from its file when the file is found and its sequence is unchanged (an open tab of that file is reused), otherwise from the copy saved in the experiment, in a new unsaved tab; `warnings` says which and why. The active document stays active. Returns the experiment with the open documents' ids in `lanes`, ready for `gel_experiment`.",
            open,
        )
        .example(json!({"path": "experiment.sslab"})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Params {
    /// Agarose, % (default 1).
    #[serde(default = "one")]
    agarose: f64,
    /// The ladder in lane 1 (default "1 kb").
    #[serde(default = "one_kb")]
    ladder: String,
    lanes: Vec<ExperimentLane>,
}

fn one_kb() -> String {
    "1 kb".to_string()
}

fn gel_experiment(s: &Session, p: Params) -> Result<Value, CmdError> {
    if !p.agarose.is_finite() || !(0.1..=5.0).contains(&p.agarose) {
        return Err(CmdError::invalid_params("agarose must be 0.1-5 %"));
    }
    if p.lanes.is_empty() || p.lanes.len() >= MAX_LANES {
        return Err(CmdError::invalid_params(format!(
            "give 1 to {} lanes besides the ladder",
            MAX_LANES - 1
        )));
    }
    if let Some(l) = p
        .lanes
        .iter()
        .find(|l| s.document_by_id(l.document).is_none())
    {
        return Err(CmdError::new(
            ErrorCode::NotFound,
            format!("document {} is not open", l.document),
        )
        .with_hint("the `documents` command lists the open documents and their ids"));
    }
    let run = experiment(s, &p.ladder, p.agarose, &p.lanes);
    let mut lanes = Vec::new();
    for (i, (title, result)) in run.iter().enumerate() {
        let mut lane = json!({ "lane": i + 1, "title": title });
        if let Some(obj) = lane.as_object_mut() {
            if let Some(l) = i.checked_sub(1).and_then(|k| p.lanes.get(k)) {
                obj.insert("document".into(), json!(l.document));
            }
            match result {
                Ok(b) => obj.insert("bands".into(), json!(bands_json(b))),
                Err(e) if i == 0 => {
                    return Err(CmdError::new(ErrorCode::InvalidInput, e.clone()));
                }
                Err(e) => obj.insert("error".into(), json!(e)),
            };
        }
        lanes.push(lane);
    }
    Ok(json!({ "agarose": p.agarose, "lanes": lanes }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct SaveParams {
    /// Where to save; the name must end in .sslab.
    path: String,
    /// The experiment's name (default "Gel experiment").
    #[serde(default = "default_name")]
    name: String,
    /// Agarose, % (default 1).
    #[serde(default = "one")]
    agarose: f64,
    /// The ladder in lane 1 (default "1 kb").
    #[serde(default = "one_kb")]
    ladder: String,
    /// The lanes (kind "gel").
    #[serde(default)]
    lanes: Vec<ExperimentLane>,
    /// How the lanes were chosen: "diagnostic_digest" for a Diagnostic Digest design.
    design: Option<String>,
    /// "gel" (default) or "reads" (a construct checked with sequencing reads).
    #[serde(default = "gel_kind")]
    kind: String,
    /// Kind "reads": the construct, an open document by id.
    reference: Option<u64>,
    /// Kind "reads": the read files.
    #[serde(default)]
    reads: Vec<String>,
    /// Kind "reads": also check the reads stored with the construct in its .dna file.
    #[serde(default)]
    matched: bool,
}

fn gel_kind() -> String {
    "gel".to_string()
}

fn default_name() -> String {
    "Gel experiment".to_string()
}

fn topology_name(t: Topology) -> &'static str {
    match t {
        Topology::Circular => "circular",
        Topology::Linear => "linear",
    }
}

fn save(s: &mut Session, p: SaveParams) -> Result<Value, CmdError> {
    let lower = p.path.to_ascii_lowercase();
    if !lower.ends_with(&format!(".{EXTENSION}")) {
        return Err(CmdError::invalid_params(format!(
            "experiment files end in .{EXTENSION}"
        )));
    }
    if p.name.trim().is_empty() || p.name.chars().count() > 200 {
        return Err(CmdError::invalid_params(
            "give the experiment a name of up to 200 characters",
        ));
    }
    if p.design
        .as_deref()
        .is_some_and(|d| d != "diagnostic_digest")
    {
        return Err(CmdError::invalid_params(
            "`design` is \"diagnostic_digest\" or absent",
        ));
    }
    if !snappiesnap_bio::gel::LADDERS
        .iter()
        .any(|(n, _)| *n == p.ladder)
    {
        return Err(CmdError::invalid_params(format!(
            "unknown ladder {:?}",
            p.ladder
        )));
    }
    if p.lanes.len() > file::MAX_LANES {
        return Err(CmdError::invalid_params(format!(
            "at most {} lanes besides the ladder",
            file::MAX_LANES
        )));
    }
    let folder = folder_of(&p.path).to_string();
    // One construct per document, in order of first use.
    let mut order: Vec<u64> = Vec::new();
    let mut read_files = Vec::new();
    let mut warnings = Vec::new();
    match p.kind.as_str() {
        "gel" => {
            for l in &p.lanes {
                if !order.contains(&l.document) {
                    order.push(l.document);
                }
            }
        }
        "reads" => {
            let reference = p.reference.ok_or_else(|| {
                CmdError::invalid_params("a reads experiment needs its construct as `reference`")
            })?;
            if !p.lanes.is_empty() {
                return Err(CmdError::invalid_params("a reads experiment has no lanes"));
            }
            if p.reads.len() > file::MAX_READS {
                return Err(CmdError::invalid_params(format!(
                    "at most {} reads",
                    file::MAX_READS
                )));
            }
            order.push(reference);
            for path in &p.reads {
                if !s.platform().exists(path) {
                    warnings.push(format!("{path} wasn't found; it is saved as given"));
                }
                read_files.push(file::ReadFile {
                    name: snappiesnap_io::file_stem(path),
                    file: relative_path(&folder, path),
                    original_path: Some(path.clone()),
                    other: BTreeMap::new(),
                });
            }
        }
        other => {
            return Err(CmdError::invalid_params(format!(
                "kind {other:?}: use \"gel\" or \"reads\""
            )));
        }
    }
    let mut constructs = Vec::new();
    for id in &order {
        let doc = s.document_by_id(*id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
                .with_hint("the `documents` command lists the open documents and their ids")
        })?;
        let m = &doc.molecule;
        let genbank = snappiesnap_io::export(m, None, Format::GenBank)?;
        let genbank = String::from_utf8(genbank.bytes)
            .map_err(|_| CmdError::new(ErrorCode::Internal, "the GenBank copy isn't UTF-8 text"))?;
        if doc.path.is_none() {
            warnings.push(format!(
                "{} has never been saved: the experiment keeps a copy of it",
                m.name()
            ));
        } else if doc.dirty {
            warnings.push(format!(
                "{} has unsaved changes: the experiment keeps them in its copy, not in the file",
                m.name()
            ));
        }
        constructs.push(Construct {
            name: m.name().to_string(),
            file: doc.path.as_deref().and_then(|d| relative_path(&folder, d)),
            original_path: doc.path.clone(),
            length: m.len(),
            topology: topology_name(m.topology()).to_string(),
            fingerprint: fingerprint(m.seq().as_bytes()),
            genbank,
            other: BTreeMap::new(),
        });
    }
    let lanes = p
        .lanes
        .iter()
        .map(|l| Lane {
            construct: order.iter().position(|d| *d == l.document).unwrap_or(0),
            method: l.method.clone(),
            other: BTreeMap::new(),
        })
        .collect();
    let f = ExperimentFile {
        format: FORMAT.to_string(),
        version: VERSION,
        generator: Some(format!("SnappieSnap {}", env!("CARGO_PKG_VERSION"))),
        kind: p.kind.clone(),
        design: p.design.clone(),
        name: p.name.trim().to_string(),
        gel: GelSettings {
            agarose: p.agarose,
            ladder: p.ladder,
            other: BTreeMap::new(),
        },
        constructs,
        lanes,
        reads: read_files,
        matched_reads: p.matched && p.kind == "reads",
        other: BTreeMap::new(),
    };
    let bytes = f.to_bytes().map_err(CmdError::invalid_params)?;
    s.platform()
        .write(&p.path, &bytes)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {}: {e}", p.path)))?;
    Ok(json!({
        "path": p.path,
        "constructs": f.constructs.len(),
        "lanes": f.lanes.len(),
        "reads": f.reads.len(),
        "bytes": bytes.len(),
        "warnings": warnings,
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct OpenParams {
    /// The .sslab file.
    path: String,
}

/// How a construct will be opened.
enum Plan {
    /// Already open, at this tab.
    Tab(usize),
    /// From its file.
    File(String, Box<snappiesnap_io::Imported>),
    /// From the copy in the experiment, and why.
    Copy(Box<Molecule>, Option<String>),
}

fn open(s: &mut Session, p: OpenParams) -> Result<Value, CmdError> {
    let bytes = s
        .platform()
        .read(&p.path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {}: {e}", p.path)))?;
    let f = ExperimentFile::parse(&bytes)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, format!("{}: {e}", p.path)))?;
    let folder = folder_of(&p.path).to_string();

    // First decide for every construct, so a damaged copy fails before anything opens.
    let mut plans = Vec::new();
    let mut warnings = Vec::new();
    for c in &f.constructs {
        let mut candidates: Vec<String> = Vec::new();
        for path in [
            c.file.as_deref().and_then(|r| resolve(&folder, r)),
            c.original_path.clone(),
        ]
        .into_iter()
        .flatten()
        {
            if !candidates.contains(&path) {
                candidates.push(path);
            }
        }
        let mut plan = None;
        let mut why = None;
        for path in &candidates {
            if let Some(index) = s.document_with_path(path) {
                let same = s
                    .documents()
                    .get(index)
                    .and_then(|d| s.document_by_id(d.id))
                    .is_some_and(|d| fingerprint(d.molecule.seq().as_bytes()) == c.fingerprint);
                if same {
                    plan = Some(Plan::Tab(index));
                    break;
                }
                why = Some(format!(
                    "the open {path} differs from the experiment's copy"
                ));
                continue;
            }
            let Ok(bytes) = s.platform().read(path) else {
                continue;
            };
            match snappiesnap_io::import(path, &bytes) {
                Ok(i) if fingerprint(i.molecule.seq().as_bytes()) == c.fingerprint => {
                    plan = Some(Plan::File(path.clone(), Box::new(i)));
                    break;
                }
                Ok(_) => why = Some(format!("{path} has changed since the experiment was saved")),
                Err(e) => why = Some(format!("{path} can't be read: {e}")),
            }
        }
        // A construct opened earlier from the copy: the unsaved tab with its name and sequence.
        if plan.is_none() {
            plan = s
                .documents()
                .into_iter()
                .find(|d| {
                    d.path.is_none()
                        && d.name == c.name
                        && s.document_by_id(d.id).is_some_and(|doc| {
                            fingerprint(doc.molecule.seq().as_bytes()) == c.fingerprint
                        })
                })
                .map(|d| Plan::Tab(d.index));
        }
        let plan = match plan {
            Some(plan) => plan,
            None => {
                let mut m = snappiesnap_io::import("copy.gb", c.genbank.as_bytes())
                    .map_err(|e| {
                        CmdError::new(
                            ErrorCode::InvalidInput,
                            format!("{}: the saved copy of {} is damaged: {e}", p.path, c.name),
                        )
                    })?
                    .molecule;
                m.rename(&c.name).ok();
                let why = why.or_else(|| {
                    (!candidates.is_empty())
                        .then(|| format!("{} wasn't found", candidates.join(" or ")))
                });
                Plan::Copy(Box::new(m), why)
            }
        };
        plans.push(plan);
    }

    // Then open, keeping the active document active.
    let active = s
        .documents()
        .into_iter()
        .find(|d| d.active)
        .map(|d| d.index);
    let mut ids = Vec::new();
    let mut constructs = Vec::new();
    for (c, plan) in f.constructs.iter().zip(plans) {
        let (source, path) = match plan {
            Plan::Tab(index) => {
                let id = s.documents().get(index).map_or(0, |d| d.id);
                ids.push(id);
                constructs.push(json!({"name": c.name, "document": id, "source": "open"}));
                continue;
            }
            Plan::File(path, i) => {
                for w in &i.warnings {
                    s.log(Level::Warn, "open", format!("{path}: {w}"));
                }
                s.open_document(Document {
                    molecule: i.molecule,
                    path: Some(path.clone()),
                    dirty: false,
                    selection: None,
                    source: Some(i.source),
                });
                ("file", Some(path))
            }
            Plan::Copy(m, why) => {
                if let Some(why) = why {
                    warnings.push(format!(
                        "{}: {why}; opened the copy saved in the experiment",
                        c.name
                    ));
                }
                let mut doc = Document::new(*m);
                doc.dirty = true;
                s.open_document(doc);
                ("copy", None)
            }
        };
        let id = s.active_document_id().unwrap_or(0);
        ids.push(id);
        constructs.push(json!({"name": c.name, "document": id, "source": source, "path": path}));
    }
    if let Some(index) = active {
        s.switch_document(index)?;
    }
    // Reads are found by their relative path, else as originally saved.
    let mut reads = Vec::new();
    for r in &f.reads {
        let candidates = [
            r.file.as_deref().and_then(|rel| resolve(&folder, rel)),
            r.original_path.clone(),
        ];
        match candidates
            .into_iter()
            .flatten()
            .find(|c| s.platform().exists(c))
        {
            Some(path) => reads.push(path),
            None => warnings.push(format!("the read {} wasn't found", r.name)),
        }
    }
    for w in &warnings {
        s.log(Level::Warn, "experiment", format!("{}: {w}", p.path));
    }
    let lanes: Vec<Value> = f
        .lanes
        .iter()
        .map(|l| json!({"document": ids.get(l.construct), "method": l.method}))
        .collect();
    Ok(json!({
        "path": p.path,
        "name": f.name,
        "kind": f.kind,
        "design": f.design,
        "agarose": f.gel.agarose,
        "ladder": f.gel.ladder,
        "lanes": lanes,
        "constructs": constructs,
        "reference": if f.kind == "reads" { ids.first().copied() } else { None },
        "reads": reads,
        "matched": f.matched_reads,
        "warnings": warnings,
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DiagnosticParams {
    /// The constructs to tell apart: open documents by id (see `documents`).
    documents: Vec<u64>,
    /// Agarose, % (default 1).
    #[serde(default = "one")]
    agarose: f64,
    /// 1: single enzymes only; 2 (default): also double digests.
    #[serde(default = "two")]
    max_enzymes: usize,
    /// Most fragments per construct (default 6).
    #[serde(default = "six")]
    max_bands: usize,
    /// Most suggestions (default 10, at most 50).
    #[serde(default = "ten")]
    limit: usize,
    /// Only these enzymes (names); default: the enzymes of `set`, or every commercially
    /// available enzyme when no set is given.
    enzymes: Option<Vec<String>>,
    /// Search only an enzyme set of the user's own (see `enzyme_sets_list`): the enzymes the
    /// lab actually holds, so no suggestion needs one it hasn't got. `commercial_only` doesn't
    /// apply; set names the database doesn't know are reported in `missing_from_database`.
    set: Option<String>,
    /// With a `set`, also search the enzymes outside it. Those suggestions are marked
    /// `in_stock: false` and come after the ones the set can run (default false).
    #[serde(default)]
    include_out_of_stock: bool,
}

fn two() -> usize {
    2
}

fn six() -> usize {
    6
}

fn ten() -> usize {
    10
}

fn diagnostic_digest(s: &Session, p: DiagnosticParams) -> Result<Value, CmdError> {
    use snappiesnap_bio::diagnostic::{MAX_CONSTRUCTS, Options, design};
    if !p.agarose.is_finite() || !(0.1..=5.0).contains(&p.agarose) {
        return Err(CmdError::invalid_params("agarose must be 0.1-5 %"));
    }
    if !(1..=50).contains(&p.limit) {
        return Err(CmdError::invalid_params("`limit` must be 1-50"));
    }
    let mut ids: Vec<u64> = Vec::new();
    for id in &p.documents {
        if !ids.contains(id) {
            ids.push(*id);
        }
    }
    if ids.is_empty() || ids.len() > MAX_CONSTRUCTS {
        return Err(CmdError::invalid_params(format!(
            "give 1 to {MAX_CONSTRUCTS} documents"
        )));
    }
    let mut molecules = Vec::new();
    for id in &ids {
        let doc = s.document_by_id(*id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
                .with_hint("the `documents` command lists the open documents and their ids")
        })?;
        molecules.push(&doc.molecule);
    }
    let db = s.enzymes();
    // The enzymes the set holds, by their database name, and the set names the database
    // doesn't know (an enzyme the lab owns under a name we lack: worth saying so).
    let mut in_stock: Vec<String> = Vec::new();
    let mut missing: Vec<String> = Vec::new();
    if let Some(name) = &p.set {
        let set = db.set(name).ok_or_else(|| {
            CmdError::new(
                ErrorCode::NotFound,
                format!("no enzyme set called {name:?}"),
            )
            .with_hint("`enzyme_sets_list` lists the sets")
        })?;
        for n in &set.enzymes {
            match db.get(n) {
                Some((_, e)) => {
                    if !in_stock.iter().any(|x| x.eq_ignore_ascii_case(&e.name)) {
                        in_stock.push(e.name.clone());
                    }
                }
                None => missing.push(n.clone()),
            }
        }
        if in_stock.is_empty() {
            return Err(CmdError::new(
                ErrorCode::InvalidInput,
                format!(
                    "none of the {} enzymes of set {name:?} is in the enzyme database",
                    set.enzymes.len()
                ),
            )
            .with_hint("download REBASE (enzymes_download) to use the whole set"));
        }
    }
    let held = |e: &snappiesnap_bio::enzymes::Enzyme| {
        in_stock.iter().any(|x| x.eq_ignore_ascii_case(&e.name))
    };
    let pool: Vec<_> = match (&p.enzymes, &p.set) {
        (Some(names), _) => {
            let mut pool = Vec::new();
            for n in names {
                let (_, e) = db.get(n).ok_or_else(|| {
                    CmdError::new(ErrorCode::NotFound, format!("unknown enzyme {n:?}"))
                })?;
                pool.push(e.clone());
            }
            pool
        }
        // A set says what is at hand, so being sold commercially doesn't come into it.
        (None, Some(_)) if !p.include_out_of_stock => {
            db.enzymes.iter().filter(|e| held(e)).cloned().collect()
        }
        (None, Some(_)) => db
            .enzymes
            .iter()
            .filter(|e| held(e) || e.commercial)
            .cloned()
            .collect(),
        (None, None) => db
            .enzymes
            .iter()
            .filter(|e| e.commercial)
            .cloned()
            .collect(),
    };
    // Enzymes that cut the same way are shown under the first name: put the ones the lab uses
    // and the widely sold ones first (see `EnzymeDb::naming_rank`), so NotI isn't named CciNI.
    let mut pool: Vec<_> = pool
        .into_iter()
        .map(|e| (db.naming_rank(&e.name), e))
        .collect();
    pool.sort_by_key(|a| a.0);
    let pool: Vec<_> = pool.into_iter().map(|(_, e)| e).collect();
    let opts = Options {
        agarose: p.agarose,
        max_enzymes: p.max_enzymes,
        max_bands: p.max_bands,
        limit: p.limit,
    };
    let mut found = design(&molecules, &pool, &opts).map_err(CmdError::invalid_params)?;
    // With a set in play, a suggestion the lab can run comes before a better-scoring one it
    // can't. `design` already ordered them best first, and the sort is stable.
    let runnable = |f: &snappiesnap_bio::diagnostic::Suggestion| {
        f.enzymes
            .iter()
            .all(|n| in_stock.iter().any(|x| x.eq_ignore_ascii_case(n)))
    };
    if p.set.is_some() && p.include_out_of_stock {
        found.sort_by_key(|f| !runnable(f));
    }
    let round = |x: f64| (x * 1000.0).round() / 1000.0;
    let suggestions: Vec<Value> = found
        .iter()
        .map(|f| {
            let lanes: Vec<Value> = ids
                .iter()
                .zip(&molecules)
                .zip(&f.fragments)
                .map(|((id, m), frags)| {
                    json!({"document": id, "name": m.name(), "fragments": frags})
                })
                .collect();
            let mut v = json!({
                "enzymes": f.enzymes,
                "alternatives": f.alternatives,
                "lanes": lanes,
                "separation": f.separation.map(round),
                "closest_bands": f.closest_bands.map(round),
                "score": round(f.score),
                "notes": f.notes,
            });
            if let (Some(o), true) = (v.as_object_mut(), p.set.is_some()) {
                o.insert("in_stock".into(), json!(runnable(f)));
            }
            v
        })
        .collect();
    let mut out = json!({
        "constructs": ids.len(),
        "enzymes_searched": pool.len(),
        "agarose": p.agarose,
        "suggestions": suggestions,
    });
    if let (Some(o), Some(name)) = (out.as_object_mut(), &p.set) {
        o.insert("set".into(), json!(name));
        o.insert("in_stock_searched".into(), json!(in_stock.len()));
        if !missing.is_empty() {
            o.insert("missing_from_database".into(), json!(missing));
        }
    }
    Ok(out)
}
