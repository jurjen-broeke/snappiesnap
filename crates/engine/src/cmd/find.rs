//! Find: a sequence (both strands, IUPAC codes, mismatches, across the origin) and features or
//! primers by name.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::primers::{BindingRules, PrimerStrand, binding_sites};
use snappiesnap_bio::search::{HitStrand, MAX_HITS, MAX_MISMATCHES, find, is_nucleotide_pattern};
use snappiesnap_doc::{Range, Strand};

use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "find",
            "Find",
            "Find `query` in the open document: as a sequence (when typed in one case) on both strands (IUPAC codes match the bases they stand for; `mismatches` up to 3; across the origin of a circular molecule), and as part of a feature's or primer's name (case-insensitive; primers at their binding sites). Matches come in position order: {kind: sequence | feature | primer, range, strand, name, mismatches}.",
            find_command,
        )
        .menu("Edit/Find…")
        .shortcut("Cmd+F")
        .example(json!({"query": "ACGT"}))
        .example(json!({"query": "GAATTC", "mismatches": 1, "names": false})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct FindParams {
    query: String,
    /// Differences allowed in a sequence match (0-3; default 0).
    #[serde(default)]
    mismatches: usize,
    /// Also search feature and primer names (default true).
    names: Option<bool>,
}

/// The stretch a feature covers (from the start of its first piece along the top strand to the
/// end of its last).
fn span(segments: &[snappiesnap_doc::Segment], strand: Strand) -> Option<Range> {
    let (first, last) = (segments.first()?, segments.last()?);
    Some(if strand == Strand::Reverse {
        Range::new(last.range.start, first.range.end)
    } else {
        Range::new(first.range.start, last.range.end)
    })
}

fn find_command(s: &Session, p: FindParams) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let query: String = p.query.split_whitespace().collect();
    if query.is_empty() {
        return Err(CmdError::invalid_params("give something to find"));
    }
    if query.len() > 10_000 {
        return Err(CmdError::invalid_params("the query is too long"));
    }
    if p.mismatches > MAX_MISMATCHES {
        return Err(CmdError::invalid_params(format!(
            "at most {MAX_MISMATCHES} mismatches"
        )));
    }
    let mut out: Vec<(usize, Value)> = Vec::new();
    let mut truncated = false;
    let bytes = query.as_bytes();
    // Sequences are typed in one case (ACGT or acgt); names mix them (attB, EcoRI), even when
    // all their letters happen to be IUPAC codes.
    let one_case =
        !(bytes.iter().any(u8::is_ascii_lowercase) && bytes.iter().any(u8::is_ascii_uppercase));
    let sequence_like = one_case
        && if m.seq().alphabet().is_nucleotide() {
            is_nucleotide_pattern(bytes)
        } else {
            bytes.iter().all(u8::is_ascii_alphabetic)
        };
    if sequence_like {
        let (hits, cut) = find(m, bytes, p.mismatches);
        truncated = cut;
        out.extend(hits.into_iter().map(|h| {
            (
                h.range.start,
                json!({
                    "kind": "sequence",
                    "range": h.range,
                    "strand": match h.strand {
                        HitStrand::Forward => "forward",
                        HitStrand::Reverse => "reverse",
                    },
                    "mismatches": h.mismatches,
                }),
            )
        }));
    }
    if p.names.unwrap_or(true) {
        let needle = query.to_lowercase();
        for f in m.features() {
            if f.name.to_lowercase().contains(&needle)
                && let Some(r) = span(&f.segments, f.strand)
            {
                out.push((
                    r.start,
                    json!({"kind": "feature", "id": f.id, "name": f.name, "range": r}),
                ));
            }
        }
        let named: Vec<_> = m
            .primers()
            .iter()
            .filter(|pr| pr.name.to_lowercase().contains(&needle))
            .collect();
        if !named.is_empty() {
            let seqs: Vec<&[u8]> = named.iter().map(|pr| pr.seq.as_bytes()).collect();
            for site in binding_sites(m, &seqs, BindingRules::default()) {
                let Some(pr) = named.get(site.primer) else {
                    continue;
                };
                if !site.is_primary(pr.seq.len()) {
                    continue;
                }
                out.push((
                    site.range.start,
                    json!({
                        "kind": "primer",
                        "name": pr.name,
                        "range": site.range,
                        "strand": if site.strand == PrimerStrand::Forward { "forward" } else { "reverse" },
                    }),
                ));
            }
        }
    }
    out.sort_by_key(|(start, _)| *start);
    if out.len() > MAX_HITS {
        out.truncate(MAX_HITS);
        truncated = true;
    }
    Ok(json!({
        "matches": out.into_iter().map(|(_, v)| v).collect::<Vec<_>>(),
        "truncated": truncated,
    }))
}
