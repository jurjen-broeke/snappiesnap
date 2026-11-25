//! Finding a short sequence in a molecule: on both strands of DNA (IUPAC codes match the bases
//! they stand for, on either side), across the origin of a circular molecule, with up to a few
//! mismatches; proteins are searched on their one strand (`X` matches anything).

use serde::Serialize;
use snappiesnap_doc::{Molecule, Range, Topology};

/// Most hits reported.
pub const MAX_HITS: usize = 10_000;
/// Most mismatches allowed.
pub const MAX_MISMATCHES: usize = 3;

/// Which strand a hit is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum HitStrand {
    Forward,
    Reverse,
}

/// One place the query was found (top-strand coordinates; may wrap the origin).
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Hit {
    pub range: Range,
    pub strand: HitStrand,
    pub mismatches: usize,
}

/// IUPAC nucleotide code as a set of bases (A=1, C=2, G=4, T/U=8).
fn bases(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

fn complement(b: u8) -> u8 {
    let s = bases(b);
    // Swap A<->T and C<->G bits.
    ((s & 1) << 3) | ((s & 8) >> 3) | ((s & 2) << 1) | ((s & 4) >> 1)
}

/// Whether `query` is a nucleotide pattern (IUPAC codes only).
pub fn is_nucleotide_pattern(query: &[u8]) -> bool {
    !query.is_empty() && query.iter().all(|&b| bases(b) != 0)
}

/// Places `query` occurs in `m` with at most `mismatches` differences, in position order
/// (one hit per place, forward first), at most [`MAX_HITS`]. Returns whether the list
/// was cut short.
pub fn find(m: &Molecule, query: &[u8], mismatches: usize) -> (Vec<Hit>, bool) {
    let seq = m.seq().as_bytes();
    let (len, k) = (seq.len(), query.len());
    let mut hits = Vec::new();
    if k == 0 || k > len {
        return (hits, false);
    }
    let mismatches = mismatches.min(MAX_MISMATCHES).min(k.saturating_sub(1));
    let circular = m.topology() == Topology::Circular;
    let starts = if circular { len } else { len - k + 1 };
    let nucleotide = m.seq().alphabet().is_nucleotide();
    // Patterns as base sets: the query, and (DNA) its reverse complement read on the top strand.
    let forward: Vec<u8> = if nucleotide {
        query.iter().map(|&b| bases(b)).collect()
    } else {
        query.iter().map(u8::to_ascii_uppercase).collect()
    };
    let reverse: Option<Vec<u8>> = if nucleotide {
        let r: Vec<u8> = query.iter().rev().map(|&b| complement(b)).collect();
        // A palindrome is found once.
        (r != forward).then_some(r)
    } else {
        None
    };
    let matches = |pattern: &[u8], start: usize| -> Option<usize> {
        let mut diff = 0usize;
        for (i, &p) in pattern.iter().enumerate() {
            let b = *seq.get((start + i) % len)?;
            let same = if nucleotide {
                bases(b) & p != 0
            } else {
                let b = b.to_ascii_uppercase();
                b == p || p == b'X' || b == b'X'
            };
            if !same {
                diff += 1;
                if diff > mismatches {
                    return None;
                }
            }
        }
        Some(diff)
    };
    for start in 0..starts {
        // Where a match from here ends (past the origin it wraps; a match ending at the very
        // end ends at `len`).
        let end = if start + k <= len {
            start + k
        } else {
            start + k - len
        };
        let mut here = false;
        for (pattern, strand) in [
            (Some(forward.as_slice()), HitStrand::Forward),
            (reverse.as_deref(), HitStrand::Reverse),
        ] {
            let Some(pattern) = pattern else { continue };
            // The same stretch found on both strands is reported once (forward).
            if here {
                continue;
            }
            if let Some(d) = matches(pattern, start) {
                if hits.len() >= MAX_HITS {
                    return (hits, true);
                }
                here = true;
                hits.push(Hit {
                    range: Range::new(start, end),
                    strand,
                    mismatches: d,
                });
            }
        }
    }
    (hits, false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_doc::{Alphabet, Seq};

    fn mol(s: &str, topology: Topology) -> Molecule {
        Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), topology).unwrap()
    }

    fn at(hits: &[Hit]) -> Vec<(usize, usize, HitStrand)> {
        hits.iter()
            .map(|h| (h.range.start, h.range.end, h.strand))
            .collect()
    }

    #[test]
    fn both_strands_and_the_origin() {
        // GAATTC is a palindrome (found once); GGATCC too. ACCA on the bottom strand is TGGT.
        let m = mol("GAATTCAAATGGTCCCACCA", Topology::Circular);
        assert_eq!(
            at(&find(&m, b"GAATTC", 0).0),
            vec![(0, 6, HitStrand::Forward)]
        );
        assert_eq!(
            at(&find(&m, b"ACCA", 0).0),
            vec![(9, 13, HitStrand::Reverse), (16, 20, HitStrand::Forward)]
        );
        // Across the origin: ...CACCA|GAATT...
        assert_eq!(
            at(&find(&m, b"CCAGAA", 0).0),
            vec![(17, 3, HitStrand::Forward)]
        );
        // Linear: not across the end.
        let lin = mol("GAATTCAAATGGTCCCACCA", Topology::Linear);
        assert!(find(&lin, b"CCAGAA", 0).0.is_empty());
    }

    #[test]
    fn iupac_codes_and_mismatches() {
        let m = mol("AAAAGATCAAAAGGTCAAAA", Topology::Linear);
        // R = A or G: GRTC matches GATC and GGTC.
        assert_eq!(at(&find(&m, b"GRTC", 0).0).len(), 2);
        // One mismatch: GATC also finds GGTC; the count is reported.
        let (hits, _) = find(&m, b"GATC", 1);
        let exact: Vec<usize> = hits.iter().map(|h| h.mismatches).collect();
        assert!(exact.contains(&0) && exact.contains(&1));
        // Never more mismatches than the query has bases minus one.
        assert!(find(&m, b"G", 3).0.iter().all(|h| h.mismatches == 0));
    }

    #[test]
    fn edge_cases() {
        let m = mol("ACGT", Topology::Circular);
        assert!(find(&m, b"", 0).0.is_empty());
        assert!(find(&m, b"ACGTA", 0).0.is_empty());
        // The whole circle from every start: the query equals a rotation.
        assert_eq!(
            at(&find(&m, b"GTAC", 0).0),
            vec![(2, 2, HitStrand::Forward)]
        );
        let one = mol("A", Topology::Linear);
        assert_eq!(find(&one, b"A", 0).0.len(), 1);
        assert!(is_nucleotide_pattern(b"acgtn") && !is_nucleotide_pattern(b"lacZ"));
        // Too many hits: cut short and said so.
        let poly = mol(&"A".repeat(20_000), Topology::Linear);
        let (hits, cut) = find(&poly, b"A", 0);
        assert_eq!((hits.len(), cut), (MAX_HITS, true));
    }

    #[test]
    fn proteins_on_one_strand() {
        let p = Molecule::new(
            "p",
            Seq::new(Alphabet::Protein, "MKTAYIAKQR").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert_eq!(at(&find(&p, b"tay", 0).0), vec![(2, 5, HitStrand::Forward)]);
        assert_eq!(find(&p, b"TXY", 0).0.len(), 1);
    }
}
