//! The lab's oligo library: import a list, look at it, find the oligos that bind the open
//! molecule, and add those as its primers.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::tm::{TmConditions, tm};
use snappiesnap_doc::{Alphabet, Seq};

use crate::oligos::{LIBRARY_FILE, OligoLibrary, parse_oligo_text};
use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "oligo_library_import",
            "Import Oligo List",
            "Add oligos from a list file to the lab's oligo library, kept between sessions: one oligo per line as name, sequence and optional description separated by tabs (as SnapGene and spreadsheets export, UTF-8 or UTF-16), or FASTA. A file holding nothing but one sequence is named after the file, which is how a lab keeps one primer per file. `path` may also be a folder, and then every list, FASTA and single-primer file in it is read (its subfolders are not). An oligo already in the library under the same name is kept as it is. With `replace`, the library is emptied first.",
            import,
        )
        .example(json!({"path": "oligos.txt"}))
        .example(json!({"path": "primers/"}))
        .example(json!({"path": "oligos.txt", "replace": true})),
        CommandSpec::query(
            "oligo_library_list",
            "Oligo Library",
            "The lab's oligo library, filtered by name, sequence or description (case-insensitive substring).",
            list,
        )
        .needs_document(false)
        .example(json!({}))
        .example(json!({"filter": "IS_1", "limit": 20})),
        CommandSpec::query(
            "oligo_library_find",
            "Find Lab Oligos",
            "The library oligos that bind the open molecule, with their sites, ordered by position. Sites follow primer_sites' rules (the session's binding rules, or explicit `window`, `mismatches`, `clamp`); by default only primary sites count (exact pairing over at least 15 bases, or the whole oligo); `primary_only: false` adds the weaker ones. `in_molecule` says whether the molecule already has a primer with that sequence.",
            find,
        )
        .example(json!({}))
        .example(json!({"primary_only": false, "mismatches": 0, "window": 15})),
        CommandSpec::edit(
            "oligo_library_add",
            "Add Lab Oligos",
            "Add library oligos (by name) to the open molecule as primers, in one undoable step. An oligo the molecule already has (same name and sequence) is passed over.",
            add,
        )
        .example(json!({"names": ["M13F"]})),
        CommandSpec::session(
            "oligo_library_remove",
            "Remove Lab Oligos",
            "Remove oligos from the lab's library by name, or all of them with `all`.",
            remove,
        )
        .example(json!({"names": ["M13F"]}))
        .example(json!({"all": true})),
    ]
}

fn library_path(s: &Session) -> Result<String, CmdError> {
    let dir = s.platform().data_dir().ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotEnabled,
            "this platform has no data folder for the oligo library",
        )
    })?;
    Ok(format!(
        "{}/{LIBRARY_FILE}",
        dir.trim_end_matches(['/', '\\'])
    ))
}

/// The library in the data folder (empty when there is none yet). A damaged file is an error,
/// never silently replaced.
fn load(s: &Session) -> Result<OligoLibrary, CmdError> {
    let path = library_path(s)?;
    if !s.platform().exists(&path) {
        return Ok(OligoLibrary::default());
    }
    let bytes = s
        .platform()
        .read(&path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {path}: {e}")))?;
    serde_json::from_slice::<OligoLibrary>(&bytes)
        .map(OligoLibrary::checked)
        .map_err(|e| {
            CmdError::new(
                ErrorCode::InvalidInput,
                format!("the oligo library {path} is damaged ({e})"),
            )
            .with_hint("import a list again with \"replace\": true to start a new library")
        })
}

fn save(s: &Session, library: &OligoLibrary) -> Result<String, CmdError> {
    let path = library_path(s)?;
    let bytes = serde_json::to_vec_pretty(library)
        .map_err(|e| CmdError::new(ErrorCode::Internal, e.to_string()))?;
    s.platform()
        .write(&path, &bytes)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot write {path}: {e}")))?;
    Ok(path)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ImportParams {
    path: String,
    /// Empty the library first (default: add to it).
    #[serde(default)]
    replace: bool,
}

/// The file's own name without its folders or extension, which names an oligo the file doesn't
/// name itself.
fn file_stem(path: &str) -> &str {
    let name = path
        .rsplit_once(['/', '\\'])
        .map_or(path, |(_, name)| name)
        .trim();
    // Only a real extension is dropped, never a leading dot or a name that is all dots.
    match name.rsplit_once('.') {
        Some((stem, _)) if !stem.is_empty() => stem,
        _ => name,
    }
}

/// Extensions read from a folder. Anything else there is left alone, so a folder of primers
/// beside its maps and notes still imports.
const LIST_EXTENSIONS: [&str; 6] = ["txt", "tsv", "seq", "fa", "fas", "fasta"];
/// Most files read from one folder.
const MAX_FILES: usize = 1_000;

/// The files to read for `path`: the file itself, or the lists in a folder.
fn files_to_read(s: &Session, path: &str) -> Result<Vec<String>, CmdError> {
    let read_error = match s.platform().read(path) {
        Ok(_) => return Ok(vec![path.to_string()]),
        Err(e) => e,
    };
    let Ok(entries) = s.platform().list_dir(path) else {
        return Err(CmdError::new(
            ErrorCode::Io,
            format!("cannot read {path}: {read_error}"),
        ));
    };
    let files: Vec<String> = entries
        .into_iter()
        .filter(|e| {
            !e.is_dir
                && LIST_EXTENSIONS.contains(&file_extension(&e.name).to_ascii_lowercase().as_str())
        })
        .map(|e| e.path)
        .collect();
    if files.is_empty() {
        return Err(CmdError::new(
            ErrorCode::InvalidInput,
            format!("the folder {path} holds no oligo lists"),
        )
        .with_hint(format!(".{} are read", LIST_EXTENSIONS.join(", ."))));
    }
    if files.len() > MAX_FILES {
        return Err(CmdError::new(
            ErrorCode::InvalidInput,
            format!(
                "the folder {path} holds {} files; at most {MAX_FILES} are read at once",
                files.len()
            ),
        ));
    }
    Ok(files)
}

fn file_extension(name: &str) -> &str {
    name.rsplit_once('.').map_or("", |(_, ext)| ext)
}

fn import(s: &mut Session, p: ImportParams) -> Result<Value, CmdError> {
    let files = files_to_read(s, &p.path)?;
    let one_file = files.len() == 1 && files.first().map(String::as_str) == Some(p.path.as_str());
    let mut oligos = Vec::new();
    let mut skipped = Vec::new();
    // Files that held nothing readable. One file on its own is an error; in a folder, the
    // others are still imported and the failures reported.
    let mut unreadable: Vec<Value> = Vec::new();
    let mut read = 0usize;
    for file in &files {
        let note = |e: String| json!({"file": file_name(file), "reason": e});
        let bytes = match s.platform().read(file) {
            Ok(b) => b,
            Err(e) => {
                unreadable.push(note(e.to_string()));
                continue;
            }
        };
        match parse_oligo_text(&bytes, Some(file_stem(file))) {
            Ok((found, lines)) => {
                read = read.saturating_add(1);
                oligos.extend(found);
                skipped.extend(if one_file {
                    lines
                } else {
                    lines
                        .into_iter()
                        .map(|l| l.in_file(file_name(file)))
                        .collect()
                });
            }
            Err(e) => unreadable.push(note(e)),
        }
    }
    if read == 0 {
        let why = unreadable
            .first()
            .and_then(|v| v.get("reason"))
            .and_then(Value::as_str)
            .unwrap_or("nothing could be read");
        return Err(CmdError::new(
            ErrorCode::InvalidInput,
            format!("{}: {why}", p.path),
        ));
    }
    let mut library = if p.replace {
        OligoLibrary::default()
    } else {
        load(s)?
    };
    let report = library
        .add(oligos)
        .map_err(CmdError::invalid_params)?
        .with_skipped(skipped);
    let saved = save(s, &library)?;
    let mut v = json!(report);
    if let Some(o) = v.as_object_mut() {
        o.insert("total".into(), json!(library.oligos.len()));
        o.insert("saved".into(), json!(saved));
        if !one_file {
            o.insert("files".into(), json!(read));
        }
        if !unreadable.is_empty() {
            o.insert("unreadable".into(), json!(unreadable));
        }
    }
    Ok(v)
}

fn file_name(path: &str) -> &str {
    path.rsplit_once(['/', '\\']).map_or(path, |(_, n)| n)
}

#[derive(Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct ListParams {
    filter: String,
    limit: usize,
}

impl Default for ListParams {
    fn default() -> Self {
        ListParams {
            filter: String::new(),
            limit: 500,
        }
    }
}

fn list(s: &Session, p: ListParams) -> Result<Value, CmdError> {
    let library = load(s)?;
    let needle = p.filter.trim().to_lowercase();
    let matching: Vec<Value> = library
        .oligos
        .iter()
        .filter(|o| {
            needle.is_empty()
                || o.name.to_lowercase().contains(&needle)
                || o.residues.to_lowercase().contains(&needle)
                || o.description.to_lowercase().contains(&needle)
        })
        .map(|o| {
            json!({
                "name": o.name,
                "residues": o.residues,
                "length": o.residues.len(),
                "description": o.description,
            })
        })
        .collect();
    Ok(json!({
        "total": library.oligos.len(),
        "matching": matching.len(),
        "oligos": matching.into_iter().take(p.limit.min(10_000)).collect::<Vec<_>>(),
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct FindParams {
    /// 3'-terminal bases that must pair (at least 10).
    window: Option<usize>,
    /// Mismatches allowed in the window.
    mismatches: Option<usize>,
    /// 3'-terminal bases that must pair exactly.
    clamp: Option<usize>,
    /// Only primary sites: exact pairing over at least 15 bases, or the whole oligo (default).
    primary_only: bool,
}

impl Default for FindParams {
    fn default() -> Self {
        FindParams {
            window: None,
            mismatches: None,
            clamp: None,
            primary_only: true,
        }
    }
}

fn find(s: &Session, p: FindParams) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let library = load(s)?;
    let seqs: Vec<&[u8]> = library
        .oligos
        .iter()
        .map(|o| o.residues.as_bytes())
        .collect();
    let conditions = TmConditions::default();
    let sites =
        super::primer::find_sites(s, m, &seqs, [p.window, p.mismatches, p.clamp], &conditions);
    let mut found: Vec<(usize, Value)> = Vec::new();
    for (i, o) in library.oligos.iter().enumerate() {
        let mut mine: Vec<_> = sites
            .iter()
            .filter(|b| b.primer == i)
            .filter(|b| !p.primary_only || b.is_primary(o.residues.len()))
            .collect();
        mine.sort_by_key(|b| b.range.start);
        let Some(first) = mine.iter().map(|b| b.range.start).min() else {
            continue;
        };
        let upper = o.residues.to_ascii_uppercase();
        let in_molecule = m
            .primers()
            .iter()
            .any(|x| x.seq.as_bytes().eq_ignore_ascii_case(o.residues.as_bytes()));
        let round = |v: Option<f64>| v.map(|t| (t * 10.0).round() / 10.0);
        found.push((
            first,
            json!({
                "name": o.name,
                "residues": o.residues,
                "length": o.residues.len(),
                "description": o.description,
                "tm": round(tm(upper.as_bytes(), &conditions)),
                "in_molecule": in_molecule,
                "sites": mine.iter().map(|b| json!({
                    "strand": b.strand,
                    "range": b.range,
                    "annealed": b.annealed,
                    "mismatches": b.mismatches,
                    "primary": b.is_primary(o.residues.len()),
                })).collect::<Vec<_>>(),
            }),
        ));
    }
    found.sort_by_key(|(first, _)| *first);
    Ok(json!({
        "library": library.oligos.len(),
        "oligos": found.into_iter().map(|(_, v)| v).collect::<Vec<_>>(),
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct NamesParams {
    names: Vec<String>,
}

/// Most oligos added in one step.
const MAX_ADD: usize = 1_000;

fn add(s: &mut Session, p: NamesParams) -> Result<Value, CmdError> {
    if p.names.is_empty() || p.names.len() > MAX_ADD {
        return Err(CmdError::invalid_params(format!(
            "name 1 to {MAX_ADD} oligos"
        )));
    }
    let library = load(s)?;
    let mut chosen = Vec::new();
    let mut unknown = Vec::new();
    for name in &p.names {
        match library
            .oligos
            .iter()
            .find(|o| o.name.eq_ignore_ascii_case(name.trim()))
        {
            Some(o) if !chosen.contains(&o) => chosen.push(o),
            Some(_) => {}
            None => unknown.push(name.as_str()),
        }
    }
    if !unknown.is_empty() {
        return Err(CmdError::new(
            ErrorCode::NotFound,
            format!("not in the oligo library: {}", unknown.join(", ")),
        )
        .with_hint("oligo_library_list lists the library"));
    }
    let m = s.molecule()?;
    let (new, present): (Vec<_>, Vec<_>) = chosen.into_iter().partition(|o| {
        !m.primers().iter().any(|x| {
            x.name.eq_ignore_ascii_case(&o.name)
                && x.seq.as_bytes().eq_ignore_ascii_case(o.residues.as_bytes())
        })
    });
    if new.is_empty() {
        return Err(CmdError::new(
            ErrorCode::NotEnabled,
            "the molecule already has these primers",
        ));
    }
    let mut ids = Vec::new();
    let m = s.molecule_mut()?;
    for o in &new {
        let seq = Seq::new(Alphabet::Dna, o.residues.as_bytes())?;
        ids.push(m.add_primer(&o.name, seq, &o.description)?);
    }
    Ok(json!({
        "added": new.iter().zip(&ids).map(|(o, id)| json!({"id": id, "name": o.name})).collect::<Vec<_>>(),
        "already_there": present.iter().map(|o| o.name.as_str()).collect::<Vec<_>>(),
    }))
}

#[derive(Default, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct RemoveParams {
    names: Vec<String>,
    /// Remove every oligo.
    all: bool,
}

fn remove(s: &mut Session, p: RemoveParams) -> Result<Value, CmdError> {
    if !p.all && p.names.is_empty() {
        return Err(CmdError::invalid_params(
            "name the oligos to remove, or pass \"all\": true",
        ));
    }
    let mut library = load(s)?;
    let before = library.oligos.len();
    if p.all {
        library.oligos.clear();
    } else {
        library.oligos.retain(|o| {
            !p.names
                .iter()
                .any(|n| n.trim().eq_ignore_ascii_case(&o.name))
        });
    }
    let removed = before - library.oligos.len();
    save(s, &library)?;
    Ok(json!({ "removed": removed, "total": library.oligos.len() }))
}
