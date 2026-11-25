//! Built-in fluorescent proteins from FPbase, found by translated search.
//!
//! The sequence is translated in all six frames (standard code, NCBI table 1) and searched for
//! every protein with the same pigeonhole filter as the DNA search (Baeza-Yates and Perleberg
//! 1992), on seeds of six amino acids; each candidate is checked for substitutions. Of finds
//! that overlap on one strand only the closest protein is kept, since fluorescent-protein
//! variants differ by a residue or two. Data: FPbase (https://www.fpbase.org), CC BY-SA 4.0;
//! each protein's primary reference is credited on the features it adds.

use std::collections::{HashMap, HashSet};

use serde_json::Value;
use snappiesnap_bio::translate;
use snappiesnap_doc::{Qualifier, Range, Strand};

use super::{Entry, differences};

const FPBASE: &str = include_str!("../../../../../assets/data/features/fpbase_proteins.json");
/// Seed length, in amino acids.
const SEED: usize = 6;
/// Candidates checked at most per frame.
const MAX_CANDIDATES: usize = 1_000_000;

/// The proteins, as library entries whose residues are amino acids.
pub(super) fn entries() -> &'static [Entry] {
    static CELL: std::sync::OnceLock<Vec<Entry>> = std::sync::OnceLock::new();
    CELL.get_or_init(|| {
        let v: Value = serde_json::from_str(FPBASE).unwrap_or(Value::Null);
        v.get("proteins")
            .and_then(Value::as_array)
            .map(|list| {
                list.iter()
                    .filter_map(|p| {
                        let text = |k: &str| p.get(k).and_then(Value::as_str).unwrap_or("");
                        let (name, protein) = (text("name"), text("protein"));
                        if name.is_empty() || protein.is_empty() {
                            return None;
                        }
                        let doi = text("doi");
                        let reference = if doi.is_empty() {
                            String::new()
                        } else {
                            format!("; primary reference doi:{doi}")
                        };
                        Some(Entry {
                            name: name.to_string(),
                            kind: "CDS".into(),
                            residues: protein.to_ascii_uppercase(),
                            directional: true,
                            color: None,
                            qualifiers: vec![
                                Qualifier {
                                    name: "product".into(),
                                    values: vec![name.to_string()],
                                },
                                Qualifier {
                                    name: "note".into(),
                                    values: vec![format!(
                                        "fluorescent protein; FPbase {} (CC BY-SA 4.0){reference}",
                                        text("page")
                                    )],
                                },
                            ],
                            source: if doi.is_empty() {
                                "FPbase".into()
                            } else {
                                format!("FPbase, doi:{doi}")
                            },
                        })
                    })
                    .collect()
            })
            .unwrap_or_default()
    })
}

/// A protein found in the sequence.
pub(super) struct ProteinHit {
    pub entry: &'static Entry,
    pub range: Range,
    pub strand: Strand,
    pub mismatches: usize,
}

/// The 5-bit-per-letter code of a seed, or None for anything but A-Z.
fn seed_key(window: &[u8]) -> Option<u32> {
    window.iter().try_fold(0u32, |acc, &b| {
        b.is_ascii_uppercase()
            .then(|| (acc << 5) | u32::from(b - b'A'))
    })
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// The proteins in `dna` (uppercase, T for U) with at most `percent` of residues substituted.
pub(super) fn find(dna: &[u8], circular: bool, percent: f64) -> Vec<ProteinHit> {
    let Some(code) = translate::code(1) else {
        return Vec::new();
    };
    let len = dna.len();
    let entries = entries();
    // Per entry: the substitutions allowed (None: longer than the molecule).
    let allowed: Vec<Option<usize>> = entries
        .iter()
        .map(|e| {
            let k = e.residues.len();
            (k >= SEED && k.saturating_mul(3) <= len)
                .then(|| ((k as f64 * percent / 100.0).floor() as usize).min(k / SEED - 1))
        })
        .collect();
    let longest = entries
        .iter()
        .zip(&allowed)
        .filter(|(_, a)| a.is_some())
        .map(|(e, _)| e.residues.len() * 3)
        .max()
        .unwrap_or(0);
    if longest == 0 {
        return Vec::new();
    }
    let mut seeds: HashMap<u32, Vec<(usize, usize)>> = HashMap::new();
    for (i, e) in entries.iter().enumerate() {
        if allowed.get(i).copied().flatten().is_none() {
            continue;
        }
        for (piece, window) in e
            .residues
            .as_bytes()
            .as_chunks::<SEED>()
            .0
            .iter()
            .enumerate()
        {
            if let Some(key) = seed_key(window) {
                seeds.entry(key).or_default().push((i, piece * SEED));
            }
        }
    }
    // A circular molecule continues past its origin.
    let mut text = dna.to_vec();
    if circular {
        text.extend_from_slice(dna.get(..longest - 1).unwrap_or_default());
    }
    let total = text.len();
    let reverse: Vec<u8> = text.iter().rev().map(|&b| complement(b)).collect();
    // (forward, start in the extended text, entry, mismatches)
    let mut found: Vec<(bool, usize, usize, usize)> = Vec::new();
    for (forward, strand_text) in [(true, &text), (false, &reverse)] {
        for frame in 0..3 {
            let aa =
                translate::translate(strand_text.get(frame..).unwrap_or_default(), code, false);
            let mut seen: HashSet<(usize, usize)> = HashSet::new();
            for (at, window) in aa.windows(SEED).enumerate() {
                let Some(list) = seed_key(window).and_then(|key| seeds.get(&key)) else {
                    continue;
                };
                for &(i, offset) in list {
                    let Some(start) = at.checked_sub(offset) else {
                        continue;
                    };
                    if seen.len() >= MAX_CANDIDATES || !seen.insert((i, start)) {
                        continue;
                    }
                    let (Some(e), Some(Some(max))) = (entries.get(i), allowed.get(i)) else {
                        continue;
                    };
                    let k = e.residues.len();
                    let Some(stretch) = aa.get(start..start + k) else {
                        continue;
                    };
                    let Some(n) = differences(e.residues.as_bytes(), stretch, *max, |a, b| a == b)
                    else {
                        continue;
                    };
                    // Back to forward-strand coordinates of the (extended) sequence.
                    let p = frame + 3 * start;
                    let begin = if forward {
                        p
                    } else {
                        total.saturating_sub(p + 3 * k)
                    };
                    if begin < len {
                        found.push((forward, begin, i, n));
                    }
                }
            }
        }
    }
    // Overlapping finds on one strand: the closest protein (fewest differences, then the
    // longest, then by name).
    found.sort_unstable_by_key(|&(forward, begin, _, _)| (forward, begin));
    let mut kept: Vec<(bool, usize, usize, usize)> = Vec::new();
    let span = |i: usize| entries.get(i).map_or(0, |e| e.residues.len() * 3);
    let better = |a: (usize, usize), b: (usize, usize)| {
        let (ea, eb) = (entries.get(a.0), entries.get(b.0));
        (a.1, std::cmp::Reverse(span(a.0)), ea.map(|e| &e.name))
            < (b.1, std::cmp::Reverse(span(b.0)), eb.map(|e| &e.name))
    };
    for hit in found {
        match kept.last_mut() {
            Some(last) if last.0 == hit.0 && hit.1 < last.1 + span(last.2) => {
                if better((hit.2, hit.3), (last.2, last.3)) {
                    *last = hit;
                }
            }
            _ => kept.push(hit),
        }
    }
    kept.into_iter()
        .filter_map(|(forward, begin, i, n)| {
            let entry = entries.get(i)?;
            let bp = entry.residues.len() * 3;
            let end = if begin + bp <= len {
                begin + bp
            } else {
                begin + bp - len
            };
            Some(ProteinHit {
                entry,
                range: Range::new(begin, if end == len && begin == 0 { len } else { end }),
                strand: if forward {
                    Strand::Forward
                } else {
                    Strand::Reverse
                },
                mismatches: n,
            })
        })
        .collect()
}
