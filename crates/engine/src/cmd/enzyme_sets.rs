//! The user's own enzyme sets (such as the enzymes a lab has in stock): import one from a text
//! file, create, list and delete them. A set is chosen for the map with `enzyme_sites`' `set`.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::enzymes::{EnzymeDb, MAX_SETS, UserSet, parse_set_text};
use crate::{CmdError, CommandSpec, ErrorCode, Level, NoParams, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "enzyme_sets_list",
            "Enzyme Sets",
            "The user's own enzyme sets, with the enzymes each holds and those the enzyme database in use doesn't have.",
            list,
        )
        .needs_document(false)
        .example(json!({})),
        CommandSpec::session(
            "enzyme_set_import",
            "Import Enzyme Set",
            "Add an enzyme set from a text file of enzyme names (an optional title line, a blank line, then names separated by commas, semicolons or line breaks, as SnapGene's \"Export Enzyme Set\" writes). It is named `name`, else the file's title, else the file name, and replaces a set of that name. Kept for later sessions.",
            import,
        )
        .example(json!({"path": "enzymes.txt"}))
        .example(json!({"path": "enzymes.txt", "name": "Lab freezer"})),
        CommandSpec::session(
            "enzyme_set_create",
            "New Enzyme Set",
            "Add an enzyme set from a list of names, replacing a set of the same name. Kept for later sessions.",
            create,
        )
        .example(json!({"name": "Lab freezer", "enzymes": ["EcoRI", "BamHI", "NotI"]})),
        CommandSpec::session(
            "enzyme_set_delete",
            "Delete Enzyme Set",
            "Remove one of the user's enzyme sets.",
            delete,
        )
        .example(json!({"name": "Lab freezer"})),
    ]
}

/// A set as reported: its enzymes and those the database lacks.
fn describe(set: &UserSet, db: &EnzymeDb) -> Value {
    json!({
        "name": set.name,
        "count": set.enzymes.len(),
        "enzymes": set.enzymes,
        "missing": set.missing(db),
    })
}

fn list(s: &Session, _: NoParams) -> Result<Value, CmdError> {
    let db = s.enzymes();
    Ok(json!({
        "database": db.source,
        "sets": db.sets.iter().map(|set| describe(set, db)).collect::<Vec<_>>(),
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ImportParams {
    path: String,
    /// The set's name (default: the file's title line, else the file name).
    #[serde(default)]
    name: Option<String>,
}

fn import(s: &mut Session, p: ImportParams) -> Result<Value, CmdError> {
    let bytes = s
        .platform()
        .read(&p.path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {}: {e}", p.path)))?;
    let (title, names) = parse_set_text(&bytes)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, format!("{}: {e}", p.path)))?;
    let name = p
        .name
        .filter(|n| !n.trim().is_empty())
        .or(title)
        .unwrap_or_else(|| file_stem(&p.path));
    add(s, &name, &names, ErrorCode::InvalidInput)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct CreateParams {
    name: String,
    enzymes: Vec<String>,
}

fn create(s: &mut Session, p: CreateParams) -> Result<Value, CmdError> {
    add(s, &p.name, &p.enzymes, ErrorCode::InvalidParams)
}

/// Add (or replace) a set and save the sets. Names that can't be enzyme names are refused
/// with `bad` (invalid params, or invalid input for a file).
fn add(s: &mut Session, name: &str, names: &[String], bad: ErrorCode) -> Result<Value, CmdError> {
    const SHOWN: usize = 10;
    let (set, rejected) = UserSet::new(name, names).map_err(|e| CmdError::new(bad, e))?;
    if !rejected.is_empty() {
        let mut list = rejected
            .iter()
            .take(SHOWN)
            .cloned()
            .collect::<Vec<_>>()
            .join(", ");
        if rejected.len() > SHOWN {
            list.push_str(&format!(" and {} more", rejected.len() - SHOWN));
        }
        return Err(CmdError::new(bad, format!("not enzyme names: {list}"))
            .with_hint("enzyme names are letters and digits, such as EcoRI or Hpy188III"));
    }
    let mut sets = s.enzymes().sets.clone();
    let replaced = match sets
        .iter()
        .position(|x| x.name.eq_ignore_ascii_case(&set.name))
    {
        Some(i) => {
            sets.remove(i);
            true
        }
        None => false,
    };
    if sets.len() >= MAX_SETS {
        return Err(CmdError::invalid_params(format!(
            "at most {MAX_SETS} enzyme sets; delete one first"
        )));
    }
    sets.push(set.clone());
    let saved = save(s, sets);
    let db = s.enzymes();
    let mut v = describe(&set, db);
    if let Some(obj) = v.as_object_mut() {
        obj.insert("replaced".into(), json!(replaced));
        obj.insert("saved".into(), json!(saved));
        if !set.missing(db).is_empty() {
            let hint = if db.is_builtin() {
                "the built-in enzymes lack these; download REBASE (enzymes_download) to use them"
            } else {
                "the enzyme database in use lacks these; check their spelling"
            };
            obj.insert("hint".into(), json!(hint));
        }
    }
    Ok(v)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DeleteParams {
    name: String,
}

fn delete(s: &mut Session, p: DeleteParams) -> Result<Value, CmdError> {
    let mut sets = s.enzymes().sets.clone();
    let i = sets
        .iter()
        .position(|x| x.name.eq_ignore_ascii_case(p.name.trim()))
        .ok_or_else(|| {
            CmdError::new(
                ErrorCode::NotFound,
                format!("there is no enzyme set called {:?}", p.name),
            )
            .with_hint("enzyme_sets_list lists them")
        })?;
    let removed = sets.remove(i);
    save(s, sets);
    Ok(json!({ "deleted": removed.name }))
}

/// Use and save `sets`. A failed save is logged: the sets are still in use this session.
fn save(s: &mut Session, sets: Vec<UserSet>) -> Option<String> {
    match s.set_enzyme_sets(sets) {
        Ok(path) => path,
        Err(e) => {
            s.log(Level::Warn, "enzymes", e);
            None
        }
    }
}

/// The file name without folders or extension.
fn file_stem(path: &str) -> String {
    let name = path.rsplit(['/', '\\']).next().unwrap_or(path);
    let stem = match name.rsplit_once('.') {
        Some((stem, _)) if !stem.is_empty() => stem,
        _ => name,
    };
    if stem.trim().is_empty() {
        "Imported enzymes".to_string()
    } else {
        stem.trim().to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::file_stem;

    #[test]
    fn file_stems() {
        assert_eq!(file_stem("a/b/Lab set.txt"), "Lab set");
        assert_eq!(file_stem("C:\\x\\fridge"), "fridge");
        assert_eq!(file_stem(".txt"), ".txt");
        assert_eq!(file_stem("dir/"), "Imported enzymes");
    }
}
