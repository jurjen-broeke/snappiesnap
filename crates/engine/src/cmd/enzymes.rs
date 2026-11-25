//! Restriction enzymes: sites, digests and the enzyme database (built-in or REBASE).

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::enzymes::digest as digest_molecule;

use crate::enzymes::{EnzymeDb, REBASE_CREDIT, REBASE_URL, SiteQuery, sites};
use crate::{CmdError, CommandSpec, ErrorCode, Level, NoParams, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "enzyme_sites",
            "Restriction Sites",
            "Where enzymes cut. Without `enzymes`: commercial enzymes with at least `min_length` specific bases that cut at most `max_cuts` times (default: unique 6+ cutters). `cut` is the top-strand cut as the number of residues before it.",
            enzyme_sites,
        )
        .example(json!({}))
        .example(json!({"enzymes": ["EcoRI", "BamHI"]}))
        .example(json!({"min_length": 4, "max_cuts": 2, "commercial_only": false})),
        CommandSpec::query(
            "digest",
            "Digest",
            "Cut the molecule with enzymes; fragments longest first, with the enzyme at each end.",
            digest,
        )
        .example(json!({"enzymes": ["EcoRI", "BamHI"]})),
        CommandSpec::query(
            "enzymes_info",
            "Enzyme Database",
            "Which enzyme database is in use (built-in or REBASE), its version and size.",
            enzymes_info,
        )
        .needs_document(false)
        .example(json!({})),
        CommandSpec::query(
            "enzymes_list",
            "List Enzymes",
            "Enzymes in the database, filtered by name or site (case-insensitive substring).",
            enzymes_list,
        )
        .needs_document(false)
        .example(json!({"filter": "GGATCC"}))
        .example(json!({"filter": "bsa", "commercial_only": true, "limit": 10})),
        CommandSpec::session(
            "enzymes_download",
            "Download REBASE",
            "Download the REBASE enzyme database, use it, and keep a copy for later sessions. REBASE is © Dr. Richard J. Roberts and is never shipped with SnappieSnap.",
            enzymes_download,
        )
        .example(json!({})),
        CommandSpec::session(
            "enzymes_load",
            "Load Enzyme File",
            "Use enzymes from a REBASE withrefm file on disk.",
            enzymes_load,
        )
        .example(json!({"path": "withrefm.txt"})),
        CommandSpec::session(
            "enzymes_reset",
            "Use Built-in Enzymes",
            "Switch back to the built-in enzyme set (the downloaded copy is kept).",
            enzymes_reset,
        )
        .example(json!({})),
    ]
}

fn info(db: &EnzymeDb) -> Value {
    json!({
        "source": db.source,
        "version": db.version,
        "count": db.enzymes.len(),
        "commercial": db.enzymes.iter().filter(|e| e.commercial).count(),
        "credit": (!db.is_builtin()).then_some(REBASE_CREDIT),
    })
}

fn enzyme_sites(s: &Session, p: SiteQuery) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let found = sites(s.enzymes(), m, &p).map_err(|e| CmdError::new(ErrorCode::NotFound, e))?;
    Ok(json!({ "length": m.len(), "enzymes": found }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DigestParams {
    enzymes: Vec<String>,
}

fn digest(s: &Session, p: DigestParams) -> Result<Value, CmdError> {
    if p.enzymes.is_empty() {
        return Err(CmdError::invalid_params("name at least one enzyme"));
    }
    let q = SiteQuery {
        enzymes: Some(p.enzymes),
        ..SiteQuery::default()
    };
    let m = s.molecule()?;
    let chosen = sites(s.enzymes(), m, &q).map_err(|e| CmdError::new(ErrorCode::NotFound, e))?;
    let enzymes: Vec<_> = chosen
        .iter()
        .filter_map(|c| s.enzymes().get(&c.name).map(|(_, e)| e.clone()))
        .collect();
    let fragments = digest_molecule(m, &enzymes);
    Ok(json!({
        "length": m.len(),
        "topology": m.topology(),
        "fragments": fragments,
    }))
}

fn enzymes_info(s: &Session, _: NoParams) -> Result<Value, CmdError> {
    let mut v = info(s.enzymes());
    if let Some(obj) = v.as_object_mut() {
        let cache = s.enzyme_cache_path().filter(|p| s.platform().exists(p));
        obj.insert("cache".into(), json!(cache));
    }
    Ok(v)
}

#[derive(Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct ListParams {
    filter: String,
    commercial_only: bool,
    limit: usize,
}

impl Default for ListParams {
    fn default() -> Self {
        ListParams {
            filter: String::new(),
            commercial_only: false,
            limit: 100,
        }
    }
}

fn enzymes_list(s: &Session, p: ListParams) -> Result<Value, CmdError> {
    let needle = p.filter.trim().to_ascii_uppercase();
    let matching: Vec<Value> = s
        .enzymes()
        .enzymes
        .iter()
        .filter(|e| e.commercial || !p.commercial_only)
        .filter(|e| {
            needle.is_empty()
                || e.name.to_ascii_uppercase().contains(&needle)
                || e.site.contains(&needle)
        })
        .map(|e| {
            json!({
                "name": e.name,
                "notation": e.notation(),
                "length": e.specific_length(),
                "overhang": e.overhang(),
                "commercial": e.commercial,
            })
        })
        .collect();
    let total = matching.len();
    let limit = p.limit.min(10_000);
    Ok(json!({
        "total": total,
        "enzymes": matching.into_iter().take(limit).collect::<Vec<_>>(),
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct DownloadParams {
    /// Defaults to REBASE's withrefm file.
    url: String,
}

impl Default for DownloadParams {
    fn default() -> Self {
        DownloadParams {
            url: REBASE_URL.to_string(),
        }
    }
}

fn enzymes_download(s: &mut Session, p: DownloadParams) -> Result<Value, CmdError> {
    let bytes = s.platform().fetch(&p.url).map_err(|e| {
        CmdError::new(ErrorCode::Io, format!("download failed: {e}"))
            .with_hint("check the network connection, or load a saved file with enzymes_load")
    })?;
    let (db, skipped) = EnzymeDb::from_rebase(&bytes)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, format!("{}: {e}", p.url)))?;
    // Keep a copy so later sessions start with it. Failing to cache isn't fatal.
    let cached = match s.enzyme_cache_path() {
        Some(path) => match s.platform().write(&path, &bytes) {
            Ok(()) => Some(path),
            Err(e) => {
                s.log(Level::Warn, "enzymes", format!("REBASE not cached: {e}"));
                None
            }
        },
        None => None,
    };
    let mut v = info(&db);
    s.set_enzymes(db);
    if let Some(obj) = v.as_object_mut() {
        obj.insert("skipped".into(), json!(skipped));
        obj.insert("cache".into(), json!(cached));
    }
    Ok(v)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct LoadParams {
    path: String,
}

fn enzymes_load(s: &mut Session, p: LoadParams) -> Result<Value, CmdError> {
    let bytes = s
        .platform()
        .read(&p.path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {}: {e}", p.path)))?;
    let (db, skipped) = EnzymeDb::from_rebase(&bytes)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, format!("{}: {e}", p.path)))?;
    let mut v = info(&db);
    s.set_enzymes(db);
    if let Some(obj) = v.as_object_mut() {
        obj.insert("skipped".into(), json!(skipped));
    }
    Ok(v)
}

fn enzymes_reset(s: &mut Session, _: NoParams) -> Result<Value, CmdError> {
    s.set_enzymes(EnzymeDb::default());
    Ok(info(s.enzymes()))
}
