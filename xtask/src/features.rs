//! `cargo xtask features`: build the built-in feature dataset from public INSDC/RefSeq records.
//!
//! Reads `assets/data/features/manifest.json` (our curation: which record, which annotated
//! feature, and what we call it), downloads each record once from NCBI E-utilities with `curl`
//! and a generic User-Agent (cached in `target/feature-records/`), checks that the record has a
//! feature with exactly that key and location whose qualifiers contain `expect`, and writes
//! `assets/data/features/insdc_features.json`: each element's sequence (5' to 3' on its own
//! strand), name, type, whether it has a direction, a note, and its provenance (accession.version,
//! location, record title). Any mismatch stops the build.

use std::fs;
use std::path::Path;
use std::process::Command;

use serde_json::{Value, json};

use crate::util::{self, Result};

const EFETCH: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gbwithparts&retmode=text&id=";

pub fn run() -> Result {
    let root = util::root()?;
    let dir = root.join("assets/data/features");
    let manifest: Value = serde_json::from_str(
        &fs::read_to_string(dir.join("manifest.json"))
            .map_err(|e| format!("cannot read the manifest: {e}"))?,
    )
    .map_err(|e| format!("manifest: {e}"))?;
    let cache = root.join("target/feature-records");
    fs::create_dir_all(&cache).map_err(|e| format!("cannot create {}: {e}", cache.display()))?;
    let records = manifest
        .get("records")
        .and_then(Value::as_object)
        .ok_or("manifest: no `records`")?;
    let features = manifest
        .get("features")
        .and_then(Value::as_array)
        .ok_or("manifest: no `features`")?;
    let mut out = Vec::new();
    let mut problems = Vec::new();
    for f in features {
        let get = |k: &str| f.get(k).and_then(Value::as_str).unwrap_or("").to_string();
        let accession = get("accession");
        if !records.contains_key(&accession) {
            problems.push(format!(
                "{}: record {accession} isn't listed (with why it may be used) under `records`",
                get("name")
            ));
            continue;
        }
        let text = match record(&cache, &accession) {
            Ok(t) => t,
            Err(e) => {
                problems.push(e);
                continue;
            }
        };
        let result = if get("key").is_empty() {
            literature(&text, &get("location"), &get("citation"), &get("check"))
        } else {
            extract(&text, &get("key"), &get("location"), &get("expect"))
        };
        match result {
            Ok((residues, version, title)) => {
                let directional = f
                    .get("directional")
                    .and_then(Value::as_bool)
                    .unwrap_or(true);
                out.push(json!({
                    "name": get("name"),
                    "kind": get("kind"),
                    "residues": residues,
                    "directional": directional,
                    "note": get("note"),
                    "accession": version,
                    "location": get("location"),
                    "record": title,
                    "citation": get("citation"),
                }));
            }
            Err(e) => problems.push(format!("{} ({accession}): {e}", get("name"))),
        }
    }
    if !problems.is_empty() {
        return util::verdict("features", problems);
    }
    let dataset = json!({
        "about": "Built-in features, cut from public INSDC/RefSeq records by `cargo xtask features` from manifest.json. See SOURCES.md.",
        "features": out,
    });
    let path = dir.join("insdc_features.json");
    let text = serde_json::to_string_pretty(&dataset).map_err(|e| e.to_string())?;
    fs::write(&path, format!("{text}\n"))
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    println!("wrote {} features to {}", features.len(), path.display());
    util::verdict("features", Vec::new())
}

/// The record's GenBank text, downloaded once.
fn record(cache: &Path, accession: &str) -> Result<String> {
    if !accession
        .chars()
        .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '.')
    {
        return Err(format!("{accession}: not an accession"));
    }
    let path = cache.join(format!("{accession}.gb"));
    // A cached record without its sequence (an assembled CON record fetched as plain `gb`)
    // is fetched again.
    let complete = fs::read_to_string(&path).is_ok_and(|t| t.contains("\nORIGIN"));
    if !complete {
        let status = Command::new("curl")
            .args(["-sSf", "-A", "SnappieSnap-dev", "-o"])
            .arg(&path)
            .arg(format!("{EFETCH}{accession}"))
            .status()
            .map_err(|e| format!("cannot run curl: {e}"))?;
        // NCBI asks for at most three requests a second.
        std::thread::sleep(std::time::Duration::from_millis(400));
        if !status.success() {
            let _ = fs::remove_file(&path);
            return Err(format!("{accession}: download failed ({status})"));
        }
    }
    let text = fs::read_to_string(&path).map_err(|e| format!("{accession}: {e}"))?;
    if !text.starts_with("LOCUS") {
        let _ = fs::remove_file(&path);
        return Err(format!("{accession}: not a GenBank record"));
    }
    Ok(text)
}

/// An element whose coordinates come from a publication, not the record's feature table: the
/// residues at `location`, which must be the `check` sequence the publication gives (so the
/// coordinates are right), and the record's accession.version and title.
fn literature(
    text: &str,
    location: &str,
    citation: &str,
    check: &str,
) -> Result<(String, String, String)> {
    if citation.trim().is_empty() || check.trim().is_empty() {
        return Err("an element without `key` needs a `citation` and a `check` sequence".into());
    }
    let (version, title) = header(text)?;
    let residues = cut(&record_sequence(text)?, location)?;
    let got = String::from_utf8_lossy(&residues).into_owned();
    if !got.eq_ignore_ascii_case(check.trim()) {
        return Err(format!("{location} reads {got}, not the published {check}"));
    }
    Ok((got, version, title))
}

/// The record's accession.version and title.
fn header(text: &str) -> Result<(String, String)> {
    let version = text
        .lines()
        .find_map(|l| l.strip_prefix("VERSION"))
        .and_then(|v| v.split_whitespace().next())
        .ok_or("no VERSION line")?
        .to_string();
    let mut title = String::new();
    let mut in_definition = false;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("DEFINITION") {
            title.push_str(rest.trim());
            in_definition = true;
        } else if in_definition && line.starts_with("            ") {
            title.push(' ');
            title.push_str(line.trim());
        } else if in_definition {
            break;
        }
    }
    Ok((version, title))
}

/// The sequence lines after ORIGIN (the ORIGIN line itself may carry a comment in old
/// records, such as "1 bp upstream of lacI", which isn't sequence).
fn record_sequence(text: &str) -> Result<Vec<u8>> {
    Ok(text
        .split("\nORIGIN")
        .nth(1)
        .ok_or("no ORIGIN")?
        .split_once('\n')
        .map_or("", |(_, rest)| rest)
        .bytes()
        .filter(u8::is_ascii_alphabetic)
        .map(|b| b.to_ascii_uppercase())
        .collect())
}

/// The residues of the feature with `key` at `location` (checked to mention `expect`), the
/// record's accession.version and its title.
fn extract(
    text: &str,
    key: &str,
    location: &str,
    expect: &str,
) -> Result<(String, String, String)> {
    let (version, title) = header(text)?;
    let sequence = record_sequence(text)?;
    // The feature table: (key, location, qualifier text).
    let table = text
        .split("\nFEATURES")
        .nth(1)
        .and_then(|t| t.split("\nORIGIN").next())
        .ok_or("no FEATURES")?;
    let mut found: Option<String> = None;
    let mut current: Option<(String, String, String)> = None;
    let mut check = |entry: &Option<(String, String, String)>| {
        if let Some((k, loc, quals)) = entry
            && k == key
            && loc == location
            && quals.contains(expect)
        {
            found = Some(quals.clone());
        }
    };
    for line in table.lines().skip(1) {
        let (head, rest) = (line.get(..21).unwrap_or(""), line.get(21..).unwrap_or(""));
        if !head.trim().is_empty() && line.starts_with("     ") {
            check(&current);
            current = Some((
                head.trim().to_string(),
                rest.trim().to_string(),
                String::new(),
            ));
        } else if let Some((_, loc, quals)) = current.as_mut() {
            if quals.is_empty() && !rest.starts_with('/') {
                loc.push_str(rest.trim());
            } else {
                quals.push_str(rest.trim());
                quals.push(' ');
            }
        }
    }
    check(&current);
    let Some(qualifiers) = found else {
        return Err(format!(
            "no {key} at {location} mentioning {expect:?} in the record"
        ));
    };
    let residues = cut(&sequence, location)?;
    // A CDS must translate to the protein its record gives.
    if let Some(given) = qualifiers
        .split("/translation=\"")
        .nth(1)
        .and_then(|t| t.split('"').next())
    {
        let given: String = given.chars().filter(|c| !c.is_whitespace()).collect();
        let ours = translate(&residues);
        let ours = ours.strip_suffix('*').unwrap_or(&ours);
        // The start codon is read as M whatever it is (GTG, TTG).
        if ours.get(1..) != given.get(1..) {
            return Err(format!(
                "the cut CDS doesn't translate to the record's /translation (ours starts {}, the record's {})",
                ours.get(..12).unwrap_or(ours),
                given.get(..12).unwrap_or(&given)
            ));
        }
    }
    Ok((
        String::from_utf8_lossy(&residues).into_owned(),
        version,
        title,
    ))
}

/// Translate with the standard code (also bacterial table 11 for these genes).
fn translate(seq: &[u8]) -> String {
    const BASES: &[u8; 4] = b"TCAG";
    const AA: &[u8; 64] = b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
    seq.as_chunks::<3>()
        .0
        .iter()
        .map(|c| {
            let index = c.iter().try_fold(0usize, |acc, b| {
                BASES.iter().position(|x| x == b).map(|i| acc * 4 + i)
            });
            index.and_then(|i| AA.get(i)).map_or('X', |&a| a as char)
        })
        .collect()
}

/// The bases a location covers, 5' to 3' on the feature's strand: `a..b`, `complement(...)`,
/// `join(..., ...)` (1-based, inclusive).
fn cut(seq: &[u8], location: &str) -> Result<Vec<u8>> {
    if let Some(inner) = location
        .strip_prefix("complement(")
        .and_then(|l| l.strip_suffix(')'))
    {
        let forward = cut(seq, inner)?;
        return Ok(forward
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                other => other,
            })
            .collect());
    }
    if let Some(inner) = location
        .strip_prefix("join(")
        .and_then(|l| l.strip_suffix(')'))
    {
        let mut out = Vec::new();
        for part in inner.split(',') {
            out.extend(cut(seq, part.trim())?);
        }
        return Ok(out);
    }
    let (a, b) = location
        .split_once("..")
        .ok_or_else(|| format!("unsupported location {location}"))?;
    let parse = |s: &str| {
        s.trim()
            .parse::<usize>()
            .map_err(|_| format!("bad position in {location}"))
    };
    let (a, b) = (parse(a)?, parse(b)?);
    if a == 0 || b < a || b > seq.len() {
        return Err(format!(
            "{location} is outside the record ({} bp)",
            seq.len()
        ));
    }
    seq.get(a - 1..b)
        .map(<[u8]>::to_vec)
        .ok_or_else(|| format!("{location} is outside the record"))
}
