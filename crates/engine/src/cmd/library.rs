//! A feature library learned from the user's own annotated files, and finding its features in
//! the open document ("detect common features" without anyone else's database).
//!
//! Building reads every sequence file in a folder (and its subfolders), keeps each single-part
//! feature of at least `min_length` bp (default 15) as its own 5'->3' sequence with its name, type,
//! direction, colour and qualifiers, and drops duplicates. The library is a JSON file in the
//! app's data folder; it never leaves the user's machine. Detection finds copies on both strands,
//! across the origin of circular molecules, allowing a few substitutions (seeded search). Built-in features from public records
//! (assets/data/features/) are searched too.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use snappiesnap_bio::translate::strand_residues;
use snappiesnap_doc::{FeatureId, NewFeature, Qualifier, Range, Segment, Strand, Topology};

use crate::{CmdError, CommandSpec, ErrorCode, NoParams, Session};

mod proteins;

/// The built-in features: cut from public INSDC/RefSeq records by `cargo xtask features`
/// (assets/data/features/; sources and terms in SOURCES.md there).
const BUILTIN: &str = include_str!("../../../../assets/data/features/insdc_features.json");

/// The library file, inside the data folder.
const LIBRARY_FILE: &str = "feature_library.json";
/// Features the user added one by one (from a document, or of their own design); kept apart
/// so that rebuilding from a folder never loses them.
const MY_FILE: &str = "my_features.json";
/// The library built from downloaded vector records (feature_library_download).
const DOWNLOADED_FILE: &str = "downloaded_features.json";
/// The vector records offered for download (accessions only; never the records themselves).
const VECTOR_RECORDS: &str = include_str!("../../../../assets/data/features/vector_records.json");
const EFETCH: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id=";
/// What the user must accept before downloading (shown in the app before it asks).
const DOWNLOAD_TERMS: &str = "This downloads vector records from NCBI GenBank, many of them submitted by the companies that sell those vectors, and builds a feature library from their annotations on this computer only. Companies may hold rights in these records, so SnappieSnap never ships them: the downloaded library is for your own use on this computer and may not be shared or redistributed. Its annotations are as the submitters deposited them, years ago in some cases; SnappieSnap does not check or curate them, so names and boundaries may differ from current products. Features added from it name their source record.";
const MAX_FILES: usize = 5_000;
const MAX_DEPTH: usize = 6;
const MAX_FEATURES: usize = 20_000;
const MAX_FILE_BYTES: usize = 64 << 20;
const MAX_FEATURE_LEN: usize = 50_000;

/// Sequence files the library reads.
const EXTENSIONS: &[&str] = &[
    "dna", "gb", "gbk", "genbank", "gbff", "embl", "emb", "gff3", "gff", "sbol", "fa", "fasta",
    "fna",
];

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "feature_library_build",
            "Build Feature Library",
            "Learn a feature library from the annotated files in `folder` (and its subfolders): every single-part feature of at least `min_length` bp (default 15) becomes an entry with its sequence (5' to 3' on its strand), name, type, direction, colour and qualifiers; duplicates are dropped. Replaces the previous library; kept in the app's data folder, on this machine only.",
            build,
        )
        .needs_document(false)
        .menu("Features/Build Feature Library…")
        .example(json!({"folder": "lab", "min_length": 10})),
        CommandSpec::query(
            "feature_library",
            "Feature Library",
            "What your feature library holds: how many features were learned from which folder and how many files, how many you added (`mine`), how many are built in, the data folder they're kept in (`data_folder`: always one of the user's own, never the application's), and the features (name, type, length, source, mine; yours first, at most 1000 listed).",
            library_command,
        )
        .needs_document(false)
        .example(json!({})),
        CommandSpec::session(
            "feature_library_add",
            "Add to My Feature Library",
            "Keep a feature in your own library, so detect_features finds it in other sequences: a feature of the open document (`feature`: its id; its sequence, name, type, direction, colour and qualifiers are kept, `name` renames it), or one of your own design (`name`, `residues` in IUPAC DNA letters, optional `kind` (default misc_feature), `directional` (default true) and `note`). An entry with the same name and sequence is replaced. Kept in the app's data folder apart from the folder-built library, so rebuilding never loses it.",
            add_entry,
        )
        .needs_document(false)
        .menu("Features/Add to My Feature Library")
        .example(json!({"name": "my linker", "residues": "GGTGGCGGTGGCTCGGGCGGTGGTGGGTCG", "kind": "misc_feature", "note": "(G4S)x2"}))
        .example(json!({"feature": 1})),
        CommandSpec::session(
            "feature_library_remove",
            "Remove from My Feature Library",
            "Remove every entry called `name` from your library (the features you added and those learned from a folder). The built-in features can't be removed; leave them out of detection with `builtin: false`.",
            remove_entry,
        )
        .needs_document(false)
        .example(json!({"name": "start block"})),
        CommandSpec::session(
            "feature_library_download",
            "Download Vector Feature Library",
            "Download common vector records from NCBI GenBank (many submitted by vector suppliers) and build a feature library from their annotations, on this computer only. Requires `accept: true`: the records may carry their submitters' rights, so the library is for personal use and may not be shared (without it, the error message gives the full terms). Replaces a previous download. `url` overrides the source (testing).",
            download,
        )
        .needs_document(false)
        .menu("Features/Download Vector Feature Library…")
        .example(json!({"accept": true, "url": "http://fixtures.test/vectors.gb"})),
        CommandSpec::session(
            "feature_library_delete_download",
            "Delete Downloaded Vector Library",
            "Delete the library built by feature_library_download from this computer.",
            delete_download,
        )
        .needs_document(false)
        .example(json!({})),
        CommandSpec::query(
            "detect_features",
            "Detect Features",
            "Find known features in the open document: the built-in features (elements from public INSDC/RefSeq records, and fluorescent proteins from FPbase found by six-frame translation, the closest protein per place) and the library learned from your own files (`builtin`, `own`; both by default). Copies on either strand, across the origin of a circular molecule, with up to `max_mismatch_percent` differing bases (default 2; 0 for exact): substituted, and with `indels` (default true) also inserted or missing bases, leaving out ones already annotated at the same place. Lists {name, kind, range, strand, length, compared (dna or protein), mismatches (all differences), substitutions, insertions, deletions, identity, source}; add_detected_features adds them (a near match gets a note saying how it differs).",
            detect_command,
        )
        .menu("Features/Detect Features")
        .example(json!({}))
        .example(json!({"own": false}))
        .example(json!({"max_mismatch_percent": 0})),
        CommandSpec::edit(
            "add_detected_features",
            "Add Detected Features",
            "Add the features detect_features finds (all, or only those whose names are in `names`) as one undoable edit.",
            add_command,
        )
        .example(json!({}))
        .example(json!({"names": ["AmpR"]}))
        .example(json!({"builtin": true, "own": false})),
    ]
}

/// One feature of the library.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
struct Entry {
    name: String,
    kind: String,
    /// The feature's own sequence, 5' to 3' on its strand, upper case.
    residues: String,
    /// Whether it has a direction (forward or reverse where it was found).
    directional: bool,
    #[serde(default)]
    color: Option<String>,
    #[serde(default)]
    qualifiers: Vec<Qualifier>,
    /// The file it was first seen in.
    source: String,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
struct Library {
    folder: String,
    files: usize,
    features: Vec<Entry>,
}

/// The built-in features as a library.
fn builtin() -> &'static Library {
    static CELL: std::sync::OnceLock<Library> = std::sync::OnceLock::new();
    CELL.get_or_init(|| {
        let v: Value = serde_json::from_str(BUILTIN).unwrap_or(Value::Null);
        let features = v
            .get("features")
            .and_then(Value::as_array)
            .map(|list| {
                list.iter()
                    .filter_map(|f| {
                        let text = |k: &str| f.get(k).and_then(Value::as_str).map(str::to_string);
                        let note = text("note").unwrap_or_default();
                        Some(Entry {
                            name: text("name")?,
                            kind: text("kind")?,
                            residues: text("residues")?.to_ascii_uppercase(),
                            directional: f
                                .get("directional")
                                .and_then(Value::as_bool)
                                .unwrap_or(true),
                            color: None,
                            qualifiers: {
                                let citation = text("citation").unwrap_or_default();
                                let mut q = Vec::new();
                                if !note.is_empty() {
                                    q.push(Qualifier {
                                        name: "note".into(),
                                        values: vec![note],
                                    });
                                }
                                // Boundaries taken from a publication: say which.
                                if !citation.is_empty() {
                                    q.push(Qualifier {
                                        name: "note".into(),
                                        values: vec![format!("boundaries from {citation}")],
                                    });
                                }
                                q
                            },
                            source: format!("{} (built-in)", text("accession").unwrap_or_default()),
                        })
                    })
                    .collect()
            })
            .unwrap_or_default();
        Library {
            folder: String::new(),
            files: 0,
            features,
        }
    })
}

/// The libraries to search: the built-in features and/or the user's own (if built).
fn sources(s: &Session, which: Which) -> Result<Vec<Library>, CmdError> {
    let mut out = Vec::new();
    if which.builtin {
        out.push(builtin().clone());
    }
    if which.own {
        let mut any = false;
        for file in [LIBRARY_FILE, MY_FILE] {
            if let Some(l) = load_file(s, file)? {
                out.push(l);
                any = true;
            }
        }
        if !any && !which.builtin && !which.downloaded {
            return Err(CmdError::new(
                ErrorCode::NotFound,
                "there is no feature library of yours yet",
            )
            .with_hint(
                "build one with feature_library_build, or add features with feature_library_add",
            ));
        }
    }
    if which.downloaded
        && let Some(l) = load_file(s, DOWNLOADED_FILE)?
    {
        out.push(l);
    }
    if out.is_empty() {
        return Err(CmdError::invalid_params(
            "choose at least one source: `builtin`, `own` or `downloaded` (when there is one)",
        ));
    }
    Ok(out)
}

/// Which feature sets to search.
#[derive(Deserialize, JsonSchema, Clone, Copy)]
#[serde(deny_unknown_fields)]
struct Which {
    /// The built-in features: public INSDC/RefSeq elements and FPbase proteins (default true).
    #[serde(default = "yes")]
    builtin: bool,
    /// The library learned from your own files (default true, when there is one).
    #[serde(default = "yes")]
    own: bool,
    /// The downloaded vector library (default true, when there is one).
    #[serde(default = "yes")]
    downloaded: bool,
    /// Differences allowed (substituted, inserted or missing bases), as a percentage of a
    /// feature's length (default 2, at most 10; 0 for exact copies only). Fewer than one per
    /// 12 bp, so short sites stay exact.
    #[serde(default = "default_percent")]
    max_mismatch_percent: f64,
    /// Allow inserted and missing bases, not only substitutions (default true). Proteins found
    /// by translation allow substitutions only.
    #[serde(default = "yes")]
    indels: bool,
}

fn yes() -> bool {
    true
}

impl Default for Which {
    fn default() -> Self {
        Which {
            builtin: true,
            own: true,
            downloaded: true,
            max_mismatch_percent: default_percent(),
            indels: true,
        }
    }
}

fn library_path(s: &Session) -> Result<String, CmdError> {
    data_path(s, LIBRARY_FILE)
}

fn data_path(s: &Session, file: &str) -> Result<String, CmdError> {
    let dir = s.platform().data_dir().ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotEnabled,
            "this platform has no data folder for the feature library",
        )
    })?;
    Ok(format!("{}/{file}", dir.trim_end_matches(['/', '\\'])))
}

/// A library file in the data folder, or None when there is none yet.
fn load_file(s: &Session, file: &str) -> Result<Option<Library>, CmdError> {
    let path = data_path(s, file)?;
    if !s.platform().exists(&path) {
        return Ok(None);
    }
    let bytes = s
        .platform()
        .read(&path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {path}: {e}")))?;
    serde_json::from_slice(&bytes).map(Some).map_err(|e| {
        CmdError::new(
            ErrorCode::InvalidInput,
            format!("the library file {file} is damaged ({e})"),
        )
    })
}

fn save_file(s: &Session, file: &str, library: &Library) -> Result<String, CmdError> {
    let path = data_path(s, file)?;
    if let Some(dir) = s.platform().data_dir() {
        let _ = s.platform().create_dir(&dir);
    }
    let bytes = serde_json::to_vec(library)
        .map_err(|e| CmdError::new(ErrorCode::Internal, e.to_string()))?;
    s.platform()
        .write(&path, &bytes)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {path}: {e}")))?;
    Ok(path)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct BuildParams {
    folder: String,
    /// Shortest feature kept (default 15 bp; an exact 15 bp match is about a 1-in-20,000 chance in a plasmid).
    min_length: Option<usize>,
}

/// Sequence files under `folder`, depth first, at most [`MAX_FILES`].
fn files_under(s: &Session, folder: &str) -> Result<Vec<String>, CmdError> {
    let mut out = Vec::new();
    let mut stack = vec![(folder.to_string(), 0usize)];
    let mut first = true;
    while let Some((dir, depth)) = stack.pop() {
        let entries = match s.platform().list_dir(&dir) {
            Ok(e) => e,
            Err(e) if first => {
                return Err(CmdError::new(
                    ErrorCode::Io,
                    format!("cannot read the folder {folder}: {e}"),
                ));
            }
            Err(_) => continue,
        };
        first = false;
        for e in entries {
            if e.is_dir {
                if depth < MAX_DEPTH && !e.name.starts_with('.') {
                    stack.push((e.path, depth + 1));
                }
            } else if e
                .name
                .rsplit_once('.')
                .is_some_and(|(_, x)| EXTENSIONS.contains(&x.to_ascii_lowercase().as_str()))
            {
                out.push(e.path);
                if out.len() >= MAX_FILES {
                    return Ok(out);
                }
            }
        }
    }
    Ok(out)
}

fn build(s: &mut Session, p: BuildParams) -> Result<Value, CmdError> {
    let min = p.min_length.unwrap_or(15).max(1);
    let files = files_under(s, &p.folder)?;
    let mut library = Library {
        folder: p.folder.clone(),
        files: 0,
        features: Vec::new(),
    };
    let mut skipped = Vec::new();
    for path in &files {
        let bytes = match s.platform().read(path) {
            Ok(b) if b.len() <= MAX_FILE_BYTES => b,
            Ok(_) => {
                skipped.push(format!("{path}: too large"));
                continue;
            }
            Err(e) => {
                skipped.push(format!("{path}: {e}"));
                continue;
            }
        };
        let m = match snappiesnap_io::import(path, &bytes) {
            Ok(i) => i.molecule,
            Err(e) => {
                skipped.push(format!("{path}: {e}"));
                continue;
            }
        };
        library.files += 1;
        let source = path.rsplit(['/', '\\']).next().unwrap_or(path).to_string();
        learn(&mut library, &m, &source, min);
    }
    let path = library_path(s)?;
    if let Some(dir) = s.platform().data_dir() {
        let _ = s.platform().create_dir(&dir);
    }
    let bytes = serde_json::to_vec(&library)
        .map_err(|e| CmdError::new(ErrorCode::Internal, e.to_string()))?;
    s.platform()
        .write(&path, &bytes)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {path}: {e}")))?;
    skipped.truncate(50);
    Ok(json!({
        "features": library.features.len(),
        "files": library.files,
        "skipped": skipped,
        "path": path,
    }))
}

/// A feature's name for the library: its own, or, when that's only its type (as for many
/// GenBank features with just a `/note`), the first of its label, standard name, gene, product
/// or note, up to the first `;` or `,` and at most 60 characters. None when there is nothing
/// better than the type.
fn entry_name(f: &snappiesnap_doc::Feature) -> Option<String> {
    let own = f.name.trim();
    if !own.is_empty() && !own.eq_ignore_ascii_case(&f.kind) {
        return Some(own.to_string());
    }
    ["label", "standard_name", "gene", "product", "note"]
        .iter()
        .find_map(|key| {
            f.qualifiers
                .iter()
                .find(|q| q.name == *key)
                .and_then(|q| q.values.first())
        })
        .map(|v| {
            v.split([';', ','])
                .next()
                .unwrap_or(v)
                .split_whitespace()
                .collect::<Vec<_>>()
                .join(" ")
                .chars()
                .take(60)
                .collect::<String>()
        })
        .filter(|n| !n.is_empty() && !n.eq_ignore_ascii_case(&f.kind))
}

/// Add the single-part features of `m` (at least `min` bp) to the library, once each.
fn learn(library: &mut Library, m: &snappiesnap_doc::Molecule, source: &str, min: usize) {
    if !m.seq().alphabet().is_nucleotide() {
        return;
    }
    for f in m.features() {
        let [segment] = f.segments.as_slice() else {
            continue;
        };
        if f.kind == "source" {
            continue;
        }
        let length = segment.range.len_on(m.len());
        if length < min || length > MAX_FEATURE_LEN {
            continue;
        }
        let Some(residues) = strand_residues(m, segment.range, f.strand) else {
            continue;
        };
        let residues = String::from_utf8_lossy(&residues).to_ascii_uppercase();
        let Some(name) = entry_name(f) else {
            continue;
        };
        if library
            .features
            .iter()
            .any(|e| e.residues == residues && e.name == name && e.kind == f.kind)
        {
            continue;
        }
        if library.features.len() >= MAX_FEATURES {
            return;
        }
        library.features.push(Entry {
            name,
            kind: f.kind.clone(),
            residues,
            directional: matches!(f.strand, Strand::Forward | Strand::Reverse),
            color: f.color.clone(),
            qualifiers: f
                .qualifiers
                .iter()
                .filter(|q| q.name != "translation")
                .cloned()
                .collect(),
            source: source.to_string(),
        });
    }
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DownloadLibraryParams {
    /// That you accept the terms (personal use on this computer; no sharing). Required.
    #[serde(default)]
    accept: bool,
    /// Where to download from instead of NCBI (testing).
    url: Option<String>,
}

/// The accessions offered for download.
fn vector_accessions() -> Vec<String> {
    let v: Value = serde_json::from_str(VECTOR_RECORDS).unwrap_or(Value::Null);
    v.get("records")
        .and_then(Value::as_array)
        .map(|list| {
            list.iter()
                .filter_map(|r| r.get("accession").and_then(Value::as_str))
                .filter(|a| {
                    a.chars()
                        .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.')
                })
                .map(str::to_string)
                .collect()
        })
        .unwrap_or_default()
}

fn download(s: &mut Session, p: DownloadLibraryParams) -> Result<Value, CmdError> {
    if !p.accept {
        return Err(CmdError::invalid_params(DOWNLOAD_TERMS)
            .with_hint("run again with `accept: true` to accept these terms and download"));
    }
    let accessions = vector_accessions();
    let url = p
        .url
        .clone()
        .unwrap_or_else(|| format!("{EFETCH}{}", accessions.join(",")));
    let bytes = s.platform().fetch(&url).map_err(|e| {
        CmdError::new(ErrorCode::Io, format!("download failed: {e}"))
            .with_hint("check the network connection and try again")
    })?;
    let text = String::from_utf8_lossy(&bytes);
    let mut library = Library {
        folder: "NCBI GenBank (downloaded)".into(),
        files: 0,
        features: Vec::new(),
    };
    let mut skipped = Vec::new();
    // The records, one after another, each ending in `//`.
    for record in text.split("\n//").filter(|r| r.contains("LOCUS")) {
        let record = record.trim_start_matches(['\r', '\n']);
        let accession = record
            .lines()
            .find_map(|l| l.strip_prefix("VERSION"))
            .and_then(|v| v.split_whitespace().next())
            .unwrap_or("record")
            .to_string();
        let gb = format!("{record}\n//\n");
        match snappiesnap_io::import("record.gb", gb.as_bytes()) {
            Ok(i) => {
                library.files += 1;
                learn(
                    &mut library,
                    &i.molecule,
                    &format!("{accession} (downloaded)"),
                    15,
                );
            }
            Err(e) => skipped.push(format!("{accession}: {e}")),
        }
    }
    if library.files == 0 {
        return Err(CmdError::new(
            ErrorCode::InvalidInput,
            "the download held no GenBank records",
        ));
    }
    let path = save_file(s, DOWNLOADED_FILE, &library)?;
    skipped.truncate(20);
    Ok(json!({
        "records": library.files,
        "features": library.features.len(),
        "skipped": skipped,
        "path": path,
    }))
}

fn delete_download(s: &mut Session, _: NoParams) -> Result<Value, CmdError> {
    let path = data_path(s, DOWNLOADED_FILE)?;
    if !s.platform().exists(&path) {
        return Ok(json!({ "deleted": false }));
    }
    // An empty library in its place: the platform has no delete, and an empty file is the same
    // as none to detection.
    save_file(s, DOWNLOADED_FILE, &Library::default())?;
    Ok(json!({ "deleted": true }))
}

fn library_command(s: &Session, _: NoParams) -> Result<Value, CmdError> {
    let library = load_file(s, LIBRARY_FILE)?.unwrap_or_default();
    let mine = load_file(s, MY_FILE)?.unwrap_or_default();
    let listed: Vec<Value> = mine
        .features
        .iter()
        .map(|e| (e, true))
        .chain(library.features.iter().map(|e| (e, false)))
        .take(1000)
        .map(|(e, own)| json!({"name": e.name, "kind": e.kind, "length": e.residues.len(), "source": e.source, "mine": own}))
        .collect();
    Ok(json!({
        "folder": library.folder,
        "files": library.files,
        "count": library.features.len(),
        "mine": mine.features.len(),
        "features": listed,
        "builtin": builtin().features.len(),
        "proteins": proteins::entries().len(),
        "downloaded": load_file(s, DOWNLOADED_FILE)?.map_or(0, |l| l.features.len()),
        "data_folder": s.platform().data_dir(),
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct AddEntryParams {
    /// A feature of the open document to keep.
    feature: Option<FeatureId>,
    /// The entry's name (required for a feature of your own design; renames a kept feature).
    name: Option<String>,
    /// The feature type of your own design (default misc_feature).
    kind: Option<String>,
    /// The sequence of your own design, 5' to 3' (IUPAC DNA letters; U is read as T).
    residues: Option<String>,
    /// Whether your design has a direction (default true).
    directional: Option<bool>,
    /// A note kept with your design.
    note: Option<String>,
}

const IUPAC_DNA: &[u8] = b"ACGTRYKMSWBDHVN";

fn add_entry(s: &mut Session, p: AddEntryParams) -> Result<Value, CmdError> {
    let name = p.name.as_deref().map(str::trim).filter(|n| !n.is_empty());
    if name.is_some_and(|n| n.chars().count() > 200) {
        return Err(CmdError::invalid_params(
            "`name` is longer than 200 characters",
        ));
    }
    let entry = if let Some(id) = p.feature {
        let m = s.molecule()?;
        let f = m.feature(id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("there is no feature {}", id.0))
        })?;
        let [segment] = f.segments.as_slice() else {
            return Err(CmdError::invalid_params(
                "only single-part features can be kept in the library",
            ));
        };
        if !m.seq().alphabet().is_nucleotide() {
            return Err(CmdError::invalid_params(
                "only DNA or RNA features can be kept",
            ));
        }
        let length = segment.range.len_on(m.len());
        if length == 0 || length > MAX_FEATURE_LEN {
            return Err(CmdError::invalid_params(format!(
                "the feature must be 1 to {MAX_FEATURE_LEN} bases long"
            )));
        }
        let residues = strand_residues(m, segment.range, f.strand).ok_or_else(|| {
            CmdError::invalid_params("the feature's range is outside the sequence")
        })?;
        let doc = s
            .document()
            .map_or(String::new(), |d| d.molecule.name().to_string());
        Entry {
            name: name.map_or_else(|| f.name.clone(), str::to_string),
            kind: f.kind.clone(),
            residues: String::from_utf8_lossy(&residues)
                .to_ascii_uppercase()
                .replace('U', "T"),
            directional: matches!(f.strand, Strand::Forward | Strand::Reverse),
            color: f.color.clone(),
            qualifiers: f
                .qualifiers
                .iter()
                .filter(|q| q.name != "translation")
                .cloned()
                .collect(),
            source: if doc.is_empty() {
                "added by you".into()
            } else {
                format!("added by you from {doc}")
            },
        }
    } else {
        let (Some(name), Some(residues)) = (name, p.residues.as_deref()) else {
            return Err(CmdError::invalid_params(
                "give `feature` (a feature of the open document), or `name` and `residues` for one of your own design",
            ));
        };
        let residues: String = residues
            .chars()
            .filter(|c| !c.is_whitespace())
            .map(|c| match c.to_ascii_uppercase() {
                'U' => 'T',
                other => other,
            })
            .collect();
        if residues.is_empty() || residues.len() > MAX_FEATURE_LEN {
            return Err(CmdError::invalid_params(format!(
                "`residues` must be 1 to {MAX_FEATURE_LEN} bases"
            )));
        }
        if let Some(bad) = residues.bytes().find(|b| !IUPAC_DNA.contains(b)) {
            return Err(CmdError::invalid_params(format!(
                "`residues` holds {:?}, which isn't an IUPAC DNA letter",
                bad as char
            )));
        }
        let kind = p
            .kind
            .as_deref()
            .map(str::trim)
            .filter(|k| !k.is_empty())
            .unwrap_or("misc_feature");
        Entry {
            name: name.to_string(),
            kind: kind.chars().take(60).collect(),
            residues,
            directional: p.directional.unwrap_or(true),
            color: None,
            qualifiers: p
                .note
                .as_deref()
                .map(str::trim)
                .filter(|n| !n.is_empty())
                .map(|n| {
                    vec![Qualifier {
                        name: "note".into(),
                        values: vec![n.chars().take(2000).collect()],
                    }]
                })
                .unwrap_or_default(),
            source: "designed by you".into(),
        }
    };
    let mut mine = load_file(s, MY_FILE)?.unwrap_or_default();
    let replaced = mine
        .features
        .iter()
        .position(|e| e.name == entry.name && e.residues == entry.residues);
    let reply = json!({
        "name": entry.name,
        "length": entry.residues.len(),
        "replaced": replaced.is_some(),
    });
    let full = mine.features.len() >= MAX_FEATURES;
    match replaced.and_then(|i| mine.features.get_mut(i)) {
        Some(slot) => *slot = entry,
        None if full => {
            return Err(CmdError::invalid_params(format!(
                "your library holds the most features it can ({MAX_FEATURES})"
            )));
        }
        None => mine.features.push(entry),
    }
    save_file(s, MY_FILE, &mine)?;
    Ok(reply)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct RemoveParams {
    /// The name of the entries to remove.
    name: String,
}

fn remove_entry(s: &mut Session, p: RemoveParams) -> Result<Value, CmdError> {
    let mut removed = 0usize;
    for file in [MY_FILE, LIBRARY_FILE, DOWNLOADED_FILE] {
        if let Some(mut l) = load_file(s, file)? {
            let before = l.features.len();
            l.features.retain(|e| e.name != p.name);
            if l.features.len() != before {
                removed += before - l.features.len();
                save_file(s, file, &l)?;
            }
        }
    }
    if removed == 0 {
        return Err(CmdError::new(
            ErrorCode::NotFound,
            format!("your library has no feature called {:?}", p.name),
        ));
    }
    Ok(json!({ "removed": removed }))
}

/// A library feature found in the document.
struct Hit<'a> {
    entry: &'a Entry,
    range: Range,
    strand: Strand,
    /// Bases (amino acids, for a protein) that differ from the entry; 0 for an exact copy.
    mismatches: usize,
    /// Of those: substituted, inserted (in the document, not in the entry) and missing bases.
    edits: Edits,
    /// Found by translation: the entry's residues are amino acids.
    protein: bool,
    /// The position of its library among those searched (lower is preferred).
    rank: usize,
    /// Names of other finds of (nearly) the same stretch that this one stands for.
    also: Vec<String>,
}

/// Start and length of a range on a molecule of `len` (a range past the origin wraps).
fn extent(r: Range, len: usize) -> (usize, usize) {
    if r.end > r.start || (r.start == 0 && r.end == len) {
        (r.start, r.end - r.start)
    } else {
        (r.start, len - r.start + r.end)
    }
}

/// Whether two ranges of a molecule of `len` share a base (also across the origin).
fn overlaps(a: Range, b: Range, len: usize) -> bool {
    let ((a0, al), (b0, bl)) = (extent(a, len), extent(b, len));
    let meet = |x: usize, xl: usize, y: usize, yl: usize| x < y + yl && y < x + xl;
    meet(a0, al, b0, bl) || meet(a0 + len, al, b0, bl) || meet(a0, al, b0 + len, bl)
}

fn revcomp(s: &str) -> String {
    s.bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'T' | b'U' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'N',
        })
        .collect()
}

/// Seed length of the near-match search: each entry is cut into non-overlapping seeds of this
/// length, and a copy with fewer differences than seeds shares at least one seed exactly
/// (pigeonhole; the filter of Baeza-Yates and Perleberg 1992, "Fast and practical approximate
/// string matching"). The sequence is then scanned once for all seeds of all entries.
const SEED: usize = 12;
/// Upper bound for `max_mismatch_percent`.
const MAX_PERCENT: f64 = 10.0;
/// Candidate positions checked at most (low-complexity sequences can seed everywhere).
const MAX_CANDIDATES: usize = 2_000_000;
/// Hits returned at most.
const MAX_HITS: usize = 10_000;

fn default_percent() -> f64 {
    2.0
}

/// The 2-bit code of a seed, or None when it holds anything but A, C, G, T.
fn seed_code(bases: &[u8]) -> Option<u32> {
    bases.iter().try_fold(0u32, |acc, b| {
        let v = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        Some((acc << 2) | v)
    })
}

/// Differences between `pattern` (DNA, IUPAC codes allowed) and `text`, or None once they
/// exceed `max`.
fn mismatches(pattern: &[u8], text: &[u8], max: usize) -> Option<usize> {
    differences(pattern, text, max, base_matches)
}

/// Differences under `same`, or None once they exceed `max`.
fn differences(pattern: &[u8], text: &[u8], max: usize, same: fn(u8, u8) -> bool) -> Option<usize> {
    let mut n = 0usize;
    for (a, b) in pattern.iter().zip(text) {
        if !same(*a, *b) {
            n += 1;
            if n > max {
                return None;
            }
        }
    }
    Some(n)
}

/// Whether a library base (an IUPAC code in a designed feature) matches a sequence base.
fn base_matches(pattern: u8, base: u8) -> bool {
    let set: &[u8] = match pattern {
        b'R' => b"AG",
        b'Y' => b"CT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'S' => b"CG",
        b'W' => b"AT",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => return pattern == base,
    };
    set.contains(&base)
}

/// How a find differs from its entry.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Edits {
    substitutions: usize,
    insertions: usize,
    deletions: usize,
}

/// Where `pattern` lies in `text` near `start`, with at most `max` edits (substitutions,
/// inserted text bases, missing pattern bases): the whole pattern aligned to a stretch of
/// text that may begin and end anywhere within `max` of the seeded place. Semi-global edit
/// distance (Sellers P.H., "The theory and computation of evolutionary distances: pattern
/// recognition", J. Algorithms 1:359–373, 1980) in a diagonal band of width 2·max+1
/// (Ukkonen E., "Algorithms for approximate string matching", Inform. Control 64:100–118,
/// 1985). Returns the text stretch `start..end`, the edit count and its kinds; ties prefer
/// the seeded diagonal and substitutions.
fn fit(pattern: &[u8], text: &[u8], start: usize, max: usize) -> Option<(usize, usize, Edits)> {
    const INF: u32 = u32::MAX / 4;
    let k = pattern.len();
    let ws = start.saturating_sub(max);
    let off = start - ws;
    let window = text.get(ws..)?;
    let w = window.len().min(k + off + max);
    let band = 2 * max + 1;
    // Text column of band cell `c` on pattern row `i`: j = i + off + c - max.
    let col = |i: usize, c: usize| -> Option<usize> {
        (i + off + c).checked_sub(max).filter(|&j| j <= w)
    };
    // Steps for the traceback: 0 match, 1 substitution, 2 missing base, 3 inserted base.
    let mut steps = vec![0u8; (k + 1) * band];
    let mut prev = vec![INF; band];
    for (c, cell) in prev.iter_mut().enumerate() {
        if col(0, c).is_some() {
            *cell = 0;
        }
    }
    let mut cur = vec![INF; band];
    let at = |v: &[u32], c: usize| v.get(c).copied().unwrap_or(INF);
    for i in 1..=k {
        let mut best = INF;
        let p = *pattern.get(i - 1)?;
        for c in 0..band {
            let mut v = INF;
            let mut step = 0u8;
            if let Some(j) = col(i, c) {
                if let Some(&t) = j.checked_sub(1).and_then(|j1| window.get(j1)) {
                    let same = base_matches(p, t);
                    v = at(&prev, c).saturating_add(u32::from(!same));
                    step = u8::from(!same);
                }
                if at(&prev, c + 1).saturating_add(1) < v {
                    v = at(&prev, c + 1) + 1;
                    step = 2;
                }
                if let Some(left) = c.checked_sub(1)
                    && at(&cur, left).saturating_add(1) < v
                {
                    v = at(&cur, left) + 1;
                    step = 3;
                }
            }
            if let Some(cell) = cur.get_mut(c) {
                *cell = v;
            }
            if v < INF
                && let Some(slot) = steps.get_mut(i * band + c)
            {
                *slot = step;
            }
            best = best.min(v);
        }
        if best > max as u32 {
            return None;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    // The best end on the last row; on a tie the one nearest the seeded diagonal.
    let (end_c, edits) = (0..band)
        .filter(|&c| col(k, c).is_some())
        .map(|c| (c, at(&prev, c)))
        .filter(|&(_, v)| v <= max as u32)
        .min_by_key(|&(c, v)| (v, c.abs_diff(max)))?;
    let end_j = col(k, end_c)?;
    let (mut i, mut c) = (k, end_c);
    let mut e = Edits::default();
    while i > 0 {
        match steps.get(i * band + c).copied()? {
            0 => i -= 1,
            1 => {
                e.substitutions += 1;
                i -= 1;
            }
            2 => {
                e.deletions += 1;
                i -= 1;
                c += 1;
            }
            _ => {
                e.insertions += 1;
                c = c.checked_sub(1)?;
            }
        }
    }
    let start_j = col(0, c)?;
    debug_assert_eq!(e.substitutions + e.insertions + e.deletions, edits as usize);
    Some((ws + start_j, ws + end_j, e))
}

/// One entry on one strand, as searched.
struct Probe<'a> {
    entry: &'a Entry,
    pattern: Vec<u8>,
    forward: bool,
    max_mismatches: usize,
    /// The position of its library among those searched.
    rank: usize,
}

fn detect<'a>(
    s: &Session,
    libraries: &'a [Library],
    percent: f64,
    with_proteins: bool,
    indels: bool,
) -> Result<Vec<Hit<'a>>, CmdError> {
    if !percent.is_finite() || !(0.0..=MAX_PERCENT).contains(&percent) {
        return Err(CmdError::invalid_params(format!(
            "`max_mismatch_percent` must be between 0 and {MAX_PERCENT}"
        )));
    }
    let m = s.molecule()?;
    if !m.seq().alphabet().is_nucleotide() {
        return Err(CmdError::invalid_params(
            "features are detected in DNA or RNA only",
        ));
    }
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let mut text: Vec<u8> = m
        .seq()
        .as_bytes()
        .iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect();
    let mut probes = Vec::new();
    // Each entry with the position of its library: earlier libraries win ties (built-in first).
    for (rank, e) in libraries
        .iter()
        .enumerate()
        .flat_map(|(r, l)| l.features.iter().map(move |e| (r, e)))
    {
        let k = e.residues.len();
        if k == 0 || k > len {
            continue;
        }
        // Substitutions allowed: the percentage, rounded down, and fewer than the seeds.
        let allowed = (k as f64 * percent / 100.0).floor() as usize;
        let max_mismatches = allowed.min((k / SEED).saturating_sub(1));
        let rc = revcomp(&e.residues);
        let both = rc != e.residues;
        probes.push(Probe {
            entry: e,
            pattern: e.residues.as_bytes().to_vec(),
            forward: true,
            max_mismatches,
            rank,
        });
        if both {
            probes.push(Probe {
                entry: e,
                pattern: rc.into_bytes(),
                forward: false,
                max_mismatches,
                rank,
            });
        }
    }
    // A circular molecule continues past its origin (as far as a copy with insertions may).
    let longest = probes
        .iter()
        .map(|p| p.pattern.len() + if indels { p.max_mismatches } else { 0 })
        .max()
        .unwrap_or(0);
    if circular && longest > 1 {
        let wrap: Vec<u8> = text
            .get(..(longest - 1).min(len))
            .unwrap_or_default()
            .to_vec();
        text.extend(wrap);
    }
    // Candidate starts per probe: from the seeds, or (entries shorter than a seed) exact finds.
    let mut seeds: std::collections::HashMap<u32, Vec<(usize, usize)>> =
        std::collections::HashMap::new();
    let mut candidates: std::collections::HashSet<(usize, usize)> =
        std::collections::HashSet::new();
    for (i, p) in probes.iter_mut().enumerate() {
        let codes: Vec<Option<u32>> = p
            .pattern
            .as_chunks::<SEED>()
            .0
            .iter()
            .map(|bases| seed_code(bases))
            .collect();
        let usable = codes.iter().flatten().count();
        if usable == 0 {
            // Shorter than a seed, or every seed holds an IUPAC code: checked at every
            // position, exact copies only.
            p.max_mismatches = 0;
            let mut from = 0usize;
            while let Some(found) = text.get(from..).and_then(|t| {
                t.windows(p.pattern.len())
                    .position(|w| mismatches(&p.pattern, w, 0).is_some())
            }) {
                let start = from + found;
                from = start + 1;
                if start >= len || candidates.len() >= MAX_CANDIDATES {
                    break;
                }
                candidates.insert((i, start));
            }
            continue;
        }
        // Fewer differences than usable seeds, so a copy always shares one exactly.
        p.max_mismatches = p.max_mismatches.min(usable - 1);
        for (piece, code) in codes.into_iter().enumerate() {
            if let Some(code) = code {
                seeds.entry(code).or_default().push((i, piece * SEED));
            }
        }
    }
    if !seeds.is_empty() {
        let mask = (1u32 << (2 * SEED)) - 1;
        let mut code = 0u32;
        let mut run = 0usize;
        'scan: for (end, b) in text.iter().enumerate() {
            let v = match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => {
                    run = 0;
                    continue;
                }
            };
            code = ((code << 2) | v) & mask;
            run += 1;
            if run < SEED {
                continue;
            }
            let Some(at) = (end + 1).checked_sub(SEED) else {
                continue;
            };
            for &(i, offset) in seeds.get(&code).map(Vec::as_slice).unwrap_or_default() {
                let Some(start) = at.checked_sub(offset) else {
                    continue;
                };
                if start < len {
                    candidates.insert((i, start));
                    if candidates.len() >= MAX_CANDIDATES {
                        break 'scan;
                    }
                }
            }
        }
    }
    // Check each candidate (allowing for edits around it); per probe, of overlapping finds
    // keep the closest.
    let mut found: Vec<(usize, usize, usize, usize, Edits)> = candidates
        .into_iter()
        .filter_map(|(i, start)| {
            let p = probes.get(i)?;
            if indels && p.max_mismatches > 0 {
                let (from, to, e) = fit(&p.pattern, &text, start, p.max_mismatches)?;
                let n = e.substitutions + e.insertions + e.deletions;
                Some((i, from, to, n, e))
            } else {
                let end = start.checked_add(p.pattern.len())?;
                let window = text.get(start..end)?;
                mismatches(&p.pattern, window, p.max_mismatches).map(|n| {
                    let e = Edits {
                        substitutions: n,
                        ..Edits::default()
                    };
                    (i, start, end, n, e)
                })
            }
        })
        // A copy starting past the origin is the same copy found from the start.
        .filter(|&(_, from, to, _, _)| from < len && to > from && to - from <= len)
        .collect();
    found.sort_unstable_by_key(|&(i, from, to, n, _)| (i, from, to, n));
    let mut kept: Vec<(usize, usize, usize, usize, Edits)> = Vec::new();
    for f in found {
        match kept.last_mut() {
            Some(last) if last.0 == f.0 && f.1 < last.2 => {
                if f.3 < last.3 {
                    *last = f;
                }
            }
            _ => kept.push(f),
        }
    }
    let mut hits = Vec::new();
    for (i, start, to, n, edits) in kept {
        let Some(p) = probes.get(i) else { continue };
        let end = if to <= len { to } else { to - len };
        let range = Range::new(start, if end == len && start == 0 { len } else { end });
        let strand = match (p.entry.directional, p.forward) {
            (false, _) => Strand::None,
            (true, true) => Strand::Forward,
            (true, false) => Strand::Reverse,
        };
        // Already annotated here under this name: not a new find.
        let present = m.features().iter().any(|f| {
            f.name == p.entry.name
                && f.segments.len() == 1
                && f.segments.first().is_some_and(|sg| sg.range == range)
        });
        if !present {
            hits.push(Hit {
                entry: p.entry,
                range,
                strand,
                mismatches: n,
                edits,
                protein: false,
                rank: p.rank,
                also: Vec::new(),
            });
        }
        if hits.len() >= MAX_HITS {
            break;
        }
    }
    // Fluorescent proteins, by translation; not where the same name is found as DNA or is
    // already annotated.
    if with_proteins {
        let dna = text.get(..len).unwrap_or_default();
        for h in proteins::find(dna, circular, percent) {
            let same = |name: &str| name.eq_ignore_ascii_case(&h.entry.name);
            let known = hits
                .iter()
                .any(|d: &Hit| same(&d.entry.name) && overlaps(d.range, h.range, len))
                || m.features().iter().any(|f| {
                    same(&f.name) && f.segments.iter().any(|sg| overlaps(sg.range, h.range, len))
                });
            if !known && hits.len() < MAX_HITS {
                hits.push(Hit {
                    entry: h.entry,
                    range: h.range,
                    strand: h.strand,
                    mismatches: h.mismatches,
                    edits: Edits {
                        substitutions: h.mismatches,
                        ..Edits::default()
                    },
                    protein: true,
                    rank: 0,
                    also: Vec::new(),
                });
            }
        }
    }
    hits.sort_by_key(|h| (h.range.start, h.entry.name.clone(), h.mismatches));
    // The same stretch found under one name twice (two library entries, or a palindrome on
    // both strands): once, the closest.
    hits.dedup_by(|a, b| a.range == b.range && a.entry.name == b.entry.name);
    Ok(merge_same_stretch(hits, len))
}

/// Finds that cover (nearly) the same stretch, such as one gene known under several names
/// from several records: kept once, the best (closest match, then the earliest library, then
/// the longest), with the other names in `also`. "Nearly": each covers at least 90% of the
/// other, on compatible strands.
fn merge_same_stretch(hits: Vec<Hit<'_>>, len: usize) -> Vec<Hit<'_>> {
    // Quadratic, so only for lists of a sensible size.
    if hits.len() > 2_000 {
        return hits;
    }
    let shared = |a: &Hit, b: &Hit| -> bool {
        let compatible =
            a.strand == b.strand || a.strand == Strand::None || b.strand == Strand::None;
        if !compatible || !overlaps(a.range, b.range, len) {
            return false;
        }
        let ((a0, al), (b0, bl)) = (extent(a.range, len), extent(b.range, len));
        // Common bases, also across the origin.
        let common = |x: usize, xl: usize, y: usize, yl: usize| {
            (x + xl).min(y + yl).saturating_sub(x.max(y))
        };
        let both = common(a0, al, b0, bl)
            .max(common(a0 + len, al, b0, bl))
            .max(common(a0, al, b0 + len, bl));
        both * 10 >= al * 9 && both * 10 >= bl * 9
    };
    let score = |h: &Hit| {
        let k = h.entry.residues.len().max(1);
        // Fewer differences per base first (in millionths), then library, then length.
        (h.mismatches * 1_000_000 / k, h.rank, std::cmp::Reverse(k))
    };
    let mut kept: Vec<Hit> = Vec::new();
    for h in hits {
        match kept.iter_mut().find(|k| shared(k, &h)) {
            Some(k) => {
                if score(&h) < score(k) {
                    let mut also = std::mem::take(&mut k.also);
                    also.push(k.entry.name.clone());
                    let mut h = h;
                    h.also = also;
                    *k = h;
                } else if !k.also.contains(&h.entry.name) && h.entry.name != k.entry.name {
                    k.also.push(h.entry.name.clone());
                }
            }
            None => kept.push(h),
        }
    }
    for k in &mut kept {
        let name = k.entry.name.clone();
        k.also.retain(|n| *n != name);
        k.also.sort();
        k.also.dedup();
    }
    kept
}

/// " (2 substituted, 1 inserted, 3 missing)" when a find has inserted or missing bases;
/// nothing for substitutions alone.
fn describe_edits(e: &Edits) -> String {
    if e.insertions == 0 && e.deletions == 0 {
        return String::new();
    }
    let mut parts = Vec::new();
    if e.substitutions > 0 {
        parts.push(format!("{} substituted", e.substitutions));
    }
    if e.insertions > 0 {
        parts.push(format!("{} inserted", e.insertions));
    }
    if e.deletions > 0 {
        parts.push(format!("{} missing", e.deletions));
    }
    format!(" ({})", parts.join(", "))
}

/// Percent identity of a hit, to one decimal.
fn identity(h: &Hit) -> f64 {
    let k = h.entry.residues.len().max(1) as f64;
    ((k - h.mismatches as f64) / k * 1000.0).round() / 10.0
}

fn detect_command(s: &Session, which: Which) -> Result<Value, CmdError> {
    let libraries = sources(s, which)?;
    let hits = detect(
        s,
        &libraries,
        which.max_mismatch_percent,
        which.builtin,
        which.indels,
    )?;
    let list: Vec<Value> = hits
        .iter()
        .map(|h| {
            json!({
                "name": h.entry.name,
                "kind": h.entry.kind,
                "range": h.range,
                "strand": h.strand,
                "length": h.entry.residues.len() * if h.protein { 3 } else { 1 },
                "compared": if h.protein { "protein" } else { "dna" },
                "mismatches": h.mismatches,
                "substitutions": h.edits.substitutions,
                "insertions": h.edits.insertions,
                "deletions": h.edits.deletions,
                "identity": identity(h),
                "also": h.also,
                "source": h.entry.source,
            })
        })
        .collect();
    Ok(json!({ "features": list }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct AddParams {
    /// Only these (by name); all when not given.
    names: Option<Vec<String>>,
    #[serde(default = "yes")]
    builtin: bool,
    #[serde(default = "yes")]
    own: bool,
    #[serde(default = "yes")]
    downloaded: bool,
    /// As for detect_features (default 2).
    #[serde(default = "default_percent")]
    max_mismatch_percent: f64,
    /// As for detect_features (default true).
    #[serde(default = "yes")]
    indels: bool,
}

fn add_command(s: &mut Session, p: AddParams) -> Result<Value, CmdError> {
    let libraries = sources(
        s,
        Which {
            builtin: p.builtin,
            own: p.own,
            downloaded: p.downloaded,
            max_mismatch_percent: p.max_mismatch_percent,
            indels: p.indels,
        },
    )?;
    let new: Vec<NewFeature> = detect(s, &libraries, p.max_mismatch_percent, p.builtin, p.indels)?
        .into_iter()
        .filter(|h| {
            p.names
                .as_ref()
                .is_none_or(|n| n.iter().any(|x| x == &h.entry.name))
        })
        .map(|h| NewFeature {
            color: h.entry.color.clone(),
            name: h.entry.name.clone(),
            kind: h.entry.kind.clone(),
            strand: h.strand,
            segments: vec![Segment { range: h.range }],
            qualifiers: {
                let mut q = h.entry.qualifiers.clone();
                // A near match says so, so nobody takes it for the reference sequence.
                if h.mismatches > 0 {
                    q.push(Qualifier {
                        name: "note".into(),
                        values: vec![format!(
                            "detected with {} of {} {} differing{} ({}% identity) from {}",
                            h.mismatches,
                            h.entry.residues.len(),
                            if h.protein { "amino acids" } else { "bases" },
                            describe_edits(&h.edits),
                            identity(&h),
                            h.entry.source
                        )],
                    });
                }
                q
            },
        })
        .collect();
    let m = s.molecule_mut()?;
    let mut added = Vec::new();
    for f in new {
        added.push(m.add_feature(f)?);
    }
    Ok(json!({ "added": added.len(), "ids": added }))
}

#[cfg(test)]
mod tests {
    use super::*;

    const P: &[u8] = b"GATTACAGCCTAGGCATTCC";

    fn edits(s: usize, i: usize, d: usize) -> Edits {
        Edits {
            substitutions: s,
            insertions: i,
            deletions: d,
        }
    }

    #[test]
    fn fits_by_hand() {
        let text = |middle: &str| format!("TTTTTTTTTT{middle}TTTTTTTTTT").into_bytes();
        // Exact, seeded at 10.
        assert_eq!(
            fit(P, &text("GATTACAGCCTAGGCATTCC"), 10, 2),
            Some((10, 30, edits(0, 0, 0)))
        );
        // One substitution (C>G at pattern 8).
        assert_eq!(
            fit(P, &text("GATTACAGGCTAGGCATTCC"), 10, 2),
            Some((10, 30, edits(1, 0, 0)))
        );
        // One extra base in the text: the copy is 21 long.
        assert_eq!(
            fit(P, &text("GATTACAGCCTAAGGCATTCC"), 10, 2),
            Some((10, 31, edits(0, 1, 0)))
        );
        // One base missing from the text: 19 long.
        assert_eq!(
            fit(P, &text("GATTACAGCCTGGCATTCC"), 10, 2),
            Some((10, 29, edits(0, 0, 1)))
        );
        // Seeded a few bases off (an indel before the seed): still found at its true place.
        assert_eq!(
            fit(P, &text("GATTACAGCCTAGGCATTCC"), 12, 2),
            Some((10, 30, edits(0, 0, 0)))
        );
        // Too many differences: not found.
        assert_eq!(fit(P, &text("GATTACAGGGGGGGCATTCC"), 10, 2), None);
        // At the very start of the text, and a pattern longer than the text.
        assert_eq!(
            fit(P, b"GATTACAGCCTAGGCATTCC", 0, 2),
            Some((0, 20, edits(0, 0, 0)))
        );
        assert_eq!(fit(P, b"GATTACA", 0, 2), None);
        // No edits allowed: exact or nothing.
        assert_eq!(
            fit(P, &text("GATTACAGCCTAGGCATTCC"), 10, 0),
            Some((10, 30, edits(0, 0, 0)))
        );
        assert_eq!(fit(P, &text("GATTACAGGCTAGGCATTCC"), 10, 0), None);
    }

    #[test]
    fn edit_notes() {
        assert_eq!(describe_edits(&edits(2, 0, 0)), "");
        assert_eq!(
            describe_edits(&edits(1, 1, 3)),
            " (1 substituted, 1 inserted, 3 missing)"
        );
    }
}
