//! Going to a position or range by its numbers, and what a stretch of sequence is made of.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::composition::dna_weight;
use snappiesnap_bio::tm::{TmConditions, gc_fraction, tm};
use snappiesnap_bio::translate::protein_properties;
use snappiesnap_doc::{Alphabet, Range, Topology};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

/// Tm is given for stretches of this many bases (the nearest-neighbour model is for oligos).
const TM_LENGTHS: std::ops::RangeInclusive<usize> = 8..=60;

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "go_to",
            "Go To",
            "Select a position or range given by its numbers, as people write them: 1-based and inclusive, like \"1200\" (that base), \"100-250\", \"100..250\" or \"100 to 250\"; thousands may be grouped with commas or spaces (\"1,200\"). On a circular molecule a range whose start is past its end runs across the origin (\"4000-200\"). Returns the selection (0-based, end exclusive).",
            go_to,
        )
        .shortcut("Cmd+L")
        .needs_document(true)
        .example(json!({"location": "3-8"}))
        .example(json!({"location": "12"})),
        CommandSpec::query(
            "sequence_info",
            "Sequence Info",
            "What a stretch of the open document is made of: `range` (0-based, end exclusive, across the origin if it wraps), else the selection, else the whole molecule. DNA: length, GC %, molecular weight of one strand and of both (OligoCalc residue weights; a whole circular molecule has no ends), and the Tm of 8-60 nt stretches (nearest-neighbour, 100 nM, 100 mM Na+, as for primers). Protein: length and average mass.",
            sequence_info,
        )
        .example(json!({}))
        .example(json!({"range": {"start": 0, "end": 2}})),
    ]
}

/// A position or range as people write it (see `go_to`), as a 0-based, end-exclusive range.
pub(crate) fn parse_location(text: &str, len: usize, topology: Topology) -> Result<Range, String> {
    if len == 0 {
        return Err("the sequence is empty".into());
    }
    let cleaned: String = text
        .trim()
        .to_ascii_lowercase()
        .replace(['\u{2013}', '\u{2014}'], "-")
        .replace("..", "-")
        .replace(" to ", "-");
    let number = |s: &str| -> Result<usize, String> {
        let digits: String = s
            .chars()
            .filter(|c| !matches!(c, ',' | ' ' | '_'))
            .collect();
        if digits.is_empty() || !digits.chars().all(|c| c.is_ascii_digit()) {
            return Err(format!(
                "{:?} isn't a position: give a number like 1200, or a range like 100-250",
                text.trim()
            ));
        }
        let n: usize = digits
            .parse()
            .map_err(|_| format!("{digits} is too large"))?;
        if n == 0 || n > len {
            return Err(format!("positions run from 1 to {len}, not {n}"));
        }
        Ok(n)
    };
    let (start, end) = match cleaned.split_once('-') {
        Some((a, b)) => (number(a)?, number(b)?),
        None => {
            let n = number(&cleaned)?;
            (n, n)
        }
    };
    if start > end && topology == Topology::Linear {
        return Err(format!(
            "{start} comes after {end}; only a circular molecule has ranges across the origin"
        ));
    }
    // 1-based inclusive to 0-based end-exclusive; a wrapping range keeps start > end.
    Range::new(start - 1, end)
        .validate(len, topology)
        .map_err(|e| e.to_string())
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct GoToParams {
    /// "1200", "100-250", "100..250", "100 to 250" (1-based, inclusive).
    location: String,
}

fn go_to(s: &mut Session, p: GoToParams) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let range = parse_location(&p.location, m.len(), m.topology())
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
    let length = range.len_on(m.len());
    if let Some(doc) = s.document_mut() {
        doc.selection = Some(range);
    }
    Ok(json!({ "selection": range, "length": length }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct InfoParams {
    /// The stretch to describe; default the selection, else the whole molecule.
    range: Option<Range>,
}

fn round(x: f64, places: i32) -> f64 {
    let f = 10f64.powi(places);
    (x * f).round() / f
}

fn sequence_info(s: &Session, p: InfoParams) -> Result<Value, CmdError> {
    let doc = s.document().ok_or_else(CmdError::no_document)?;
    let m = &doc.molecule;
    let range = match p.range.or(doc.selection.filter(|r| r.len_on(m.len()) > 0)) {
        Some(r) => Some(r.validate(m.len(), m.topology())?),
        None => None,
    };
    let whole = range.is_none_or(|r| r.len_on(m.len()) == m.len());
    let residues = match range {
        Some(r) => m.residues(r)?,
        None => m.seq().as_bytes().to_vec(),
    };
    let mut out = json!({
        "range": range,
        "whole": whole,
        "length": residues.len(),
        "alphabet": m.seq().alphabet(),
    });
    let Some(o) = out.as_object_mut() else {
        return Ok(out);
    };
    match m.seq().alphabet() {
        Alphabet::Protein => {
            let props = protein_properties(&residues);
            o.insert("mass".into(), json!(round(props.mass, 2)));
        }
        alphabet => {
            o.insert(
                "gc_percent".into(),
                json!(gc_fraction(&residues).map(|f| round(f * 100.0, 1))),
            );
            if alphabet == Alphabet::Dna {
                let circular = whole && m.topology() == Topology::Circular;
                if let Some(w) = dna_weight(&residues, circular) {
                    o.insert(
                        "weight".into(),
                        json!({
                            "single_strand": round(w.single, 2),
                            "double_strand": round(w.double, 2),
                            "approximate": w.approximate,
                            "circular": circular,
                        }),
                    );
                }
                if TM_LENGTHS.contains(&residues.len()) {
                    let upper: Vec<u8> = residues.iter().map(u8::to_ascii_uppercase).collect();
                    o.insert(
                        "tm".into(),
                        json!(tm(&upper, &TmConditions::default()).map(|t| round(t, 1))),
                    );
                }
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locations() {
        let lin = |t: &str| parse_location(t, 5000, Topology::Linear);
        let circ = |t: &str| parse_location(t, 5000, Topology::Circular);
        assert_eq!(lin("1200"), Ok(Range::new(1199, 1200)));
        assert_eq!(lin(" 100-250 "), Ok(Range::new(99, 250)));
        assert_eq!(lin("100..250"), Ok(Range::new(99, 250)));
        assert_eq!(lin("100 to 250"), Ok(Range::new(99, 250)));
        assert_eq!(lin("1,000 \u{2013} 2 000"), Ok(Range::new(999, 2000)));
        assert_eq!(lin("1-5000"), Ok(Range::new(0, 5000)));
        assert_eq!(circ("4000-200"), Ok(Range::new(3999, 200)));
        assert_eq!(circ("5000-1"), Ok(Range::new(4999, 1)));
        assert!(lin("4000-200").unwrap_err().contains("circular"));
        assert!(lin("0").unwrap_err().contains("1 to 5000"));
        assert!(lin("5001").is_err());
        assert!(lin("99999999999999999999999").is_err());
        assert!(lin("abc").is_err());
        assert!(lin("").is_err());
        assert!(lin("-").is_err());
        assert!(lin("1-2-3").is_err());
        assert!(lin("1.5").is_err());
        assert_eq!(
            parse_location("1", 1, Topology::Linear),
            Ok(Range::new(0, 1))
        );
        assert_eq!(
            parse_location("1", 1, Topology::Circular),
            Ok(Range::new(0, 1))
        );
        assert!(parse_location("1", 0, Topology::Linear).is_err());
    }
}
