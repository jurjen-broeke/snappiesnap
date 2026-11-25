//! Copy and paste: the selection as text, and pasted text (FASTA, GenBank-style numbered
//! blocks, alignment rows) cleaned into residues.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_doc::{Alphabet, Range};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "copy_selection",
            "Copy",
            "The residues of the selection (or `range`) as text, for the clipboard: `format` \"plain\" (default), \"fasta\" (with the document's name as header) or \"reverse_complement\" (nucleotides). A selection across the origin of a circular molecule is copied in order.",
            copy_selection,
        )
        .menu("Edit/Copy")
        .shortcut("Cmd+C")
        .example(json!({"range": {"start": 0, "end": 1}}))
        .example(json!({"range": {"start": 0, "end": 1}, "format": "fasta"})),
        CommandSpec::edit(
            "paste",
            "Paste",
            "Paste `text` at the cursor, or over the selection (or at `at`, or over `range`). The text may be plain residues, FASTA (the first record is used), or copied with numbers, spaces and alignment gaps, which are dropped; anything else that isn't a residue of this molecule is refused. Returns where the pasted residues end.",
            paste,
        )
        .menu("Edit/Paste")
        .shortcut("Cmd+V")
        .example(json!({"text": ">x\nACGT", "at": 0}))
        .example(json!({"text": "1 acgtacgt 9", "at": 0})),
    ]
}

#[derive(Deserialize, JsonSchema, Clone, Copy, Default)]
#[serde(rename_all = "snake_case")]
enum Format {
    #[default]
    Plain,
    Fasta,
    ReverseComplement,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct CopyParams {
    /// What to copy (default: the selection).
    range: Option<Range>,
    #[serde(default)]
    format: Format,
}

/// The text a copy puts on the clipboard, and how many residues it holds. Shared with the
/// window, which copies from constructs that are not the open document (the cloning
/// experiment shows one beside its choices), so both produce exactly the same text.
pub fn copy_text(
    m: &snappiesnap_doc::Molecule,
    range: Range,
    format: &str,
) -> Result<(String, usize), CmdError> {
    let format = match format {
        "fasta" => Format::Fasta,
        "reverse_complement" => Format::ReverseComplement,
        _ => Format::Plain,
    };
    let residues = m.residues(range)?;
    let text = formatted(m, range, &residues, format)?;
    let n = residues.len();
    Ok((text, n))
}

fn copy_selection(s: &Session, p: CopyParams) -> Result<Value, CmdError> {
    let doc = s.document().ok_or_else(CmdError::no_document)?;
    let range = p
        .range
        .or(doc.selection)
        .filter(|r| !r.is_empty())
        .ok_or_else(|| CmdError::invalid_params("select something to copy"))?;
    let m = &doc.molecule;
    let residues = m.residues(range)?;
    let text = formatted(m, range, &residues, p.format)?;
    Ok(json!({ "text": text, "length": residues.len() }))
}

/// `residues` of `m` written in `format`.
fn formatted(
    m: &snappiesnap_doc::Molecule,
    range: Range,
    residues: &[u8],
    format: Format,
) -> Result<String, CmdError> {
    let text = match format {
        Format::Plain => String::from_utf8_lossy(residues).into_owned(),
        Format::ReverseComplement => {
            let alphabet = m.seq().alphabet();
            if !alphabet.is_nucleotide() {
                return Err(CmdError::invalid_params(
                    "only nucleotides have a reverse complement",
                ));
            }
            residues
                .iter()
                .rev()
                .map(|&b| char::from(alphabet.complement(b).unwrap_or(b'N')))
                .collect()
        }
        Format::Fasta => {
            let body = String::from_utf8_lossy(residues);
            let mut out = format!(
                ">{} {}..{}\n",
                m.name().replace(char::is_whitespace, "_"),
                range.start + 1,
                if range.end == 0 { m.len() } else { range.end }
            );
            for chunk in body.as_bytes().chunks(80) {
                out.push_str(&String::from_utf8_lossy(chunk));
                out.push('\n');
            }
            out
        }
    };
    Ok(text)
}

/// Residues in pasted text, and notes on what was dropped. FASTA: the first record (headers
/// and `;` comments dropped). Digits, whitespace and alignment gaps (`-`, `.`) are dropped;
/// any other character must be a residue of `alphabet`.
pub(crate) fn clean_pasted(
    text: &str,
    alphabet: Alphabet,
) -> Result<(Vec<u8>, Vec<String>), String> {
    let text = text.trim_start_matches('\u{feff}');
    let mut notes = Vec::new();
    let mut body = String::new();
    if text.trim_start().starts_with('>') {
        let records = text.lines().filter(|l| l.starts_with('>')).count();
        if records > 1 {
            notes.push(format!(
                "only the first of {records} FASTA records was pasted"
            ));
        }
        for line in text.lines().skip_while(|l| !l.starts_with('>')).skip(1) {
            if line.starts_with('>') {
                break;
            }
            if !line.starts_with(';') {
                body.push_str(line);
                body.push('\n');
            }
        }
    } else {
        body.push_str(text);
    }
    let mut out = Vec::with_capacity(body.len());
    let mut gaps = false;
    for (i, c) in body.chars().enumerate() {
        if c.is_whitespace() || c.is_ascii_digit() || c == '/' {
            continue;
        }
        if c == '-' || c == '.' {
            gaps = true;
            continue;
        }
        if c.is_ascii() && alphabet.is_valid(c as u8) {
            out.push(c as u8);
        } else {
            return Err(format!(
                "{c:?} (character {}) isn't a {} residue",
                i + 1,
                if alphabet.is_nucleotide() {
                    "nucleotide"
                } else {
                    "protein"
                }
            ));
        }
    }
    if gaps {
        notes.push("alignment gaps (- and .) were left out".to_string());
    }
    if out.is_empty() {
        return Err("there are no residues to paste".into());
    }
    Ok((out, notes))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PasteParams {
    text: String,
    /// Insert before this position instead of at the cursor.
    at: Option<usize>,
    /// Replace this range instead of the selection.
    range: Option<Range>,
}

fn paste(s: &mut Session, p: PasteParams) -> Result<Value, CmdError> {
    let doc = s.document().ok_or_else(CmdError::no_document)?;
    let len = doc.molecule.len();
    let (residues, notes) = clean_pasted(&p.text, doc.molecule.seq().alphabet())
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
    let target = match (p.range, p.at) {
        (Some(r), None) => r,
        (None, Some(at)) => Range::at(at),
        (None, None) => doc.selection.unwrap_or(Range::at(len)),
        (Some(_), Some(_)) => {
            return Err(CmdError::invalid_params("give `at` or `range`, not both"));
        }
    };
    let n = residues.len();
    // Where the pasted residues start: a selection across the origin is replaced by residues
    // at the end of the shortened molecule (see `Molecule::replace`).
    let start = if target.is_wrapping() {
        len.saturating_sub(target.len_on(len))
    } else {
        target.start
    };
    let report = if target.is_empty() {
        s.molecule_mut()?.insert(target.start, &residues)?
    } else {
        s.molecule_mut()?.replace(target, &residues)?
    };
    let new_len = s.molecule()?.len();
    // The cursor goes after the pasted residues.
    let end = (start + n).min(new_len);
    if let Some(d) = s.document_mut() {
        d.selection = Some(Range::at(end));
    }
    Ok(json!({
        "pasted": n,
        "end": end,
        "length": new_len,
        "removed_features": report.removed_features,
        "notes": notes,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn clean(t: &str) -> Result<String, String> {
        clean_pasted(t, Alphabet::Dna).map(|(r, _)| String::from_utf8(r).unwrap_or_default())
    }

    #[test]
    fn pasted_text_is_cleaned() {
        assert_eq!(clean("ACGT").unwrap(), "ACGT");
        assert_eq!(clean(">x desc\nACGT\nTT\n").unwrap(), "ACGTTT");
        assert_eq!(clean(">a\nAC\n>b\nGG").unwrap(), "AC");
        // GenBank ORIGIN style and alignment rows.
        assert_eq!(
            clean("        1 acgtacgtac gtacg\n       16 tt\n//").unwrap(),
            "acgtacgtacgtacgtt"
        );
        assert_eq!(clean("AC--GT..A").unwrap(), "ACGTA");
        let (_, notes) = clean_pasted(">a\nAC-G\n>b\nT", Alphabet::Dna).unwrap();
        assert_eq!(notes.len(), 2);
        // Not residues: refused, with where.
        let e = clean("ACGT!").unwrap_err();
        assert!(e.contains("character 5"), "{e}");
        assert!(clean("ACGTE").is_err());
        assert!(clean("  12 \n").is_err());
        assert!(clean("").is_err());
        assert!(clean("ACGTé").is_err());
        // Proteins accept amino acids.
        assert!(clean_pasted("MKVLE*", Alphabet::Protein).is_ok());
    }
}
