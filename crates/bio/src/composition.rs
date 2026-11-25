//! Composition of a stretch of DNA: GC content and molecular weight.
//!
//! Molecular weights follow OligoCalc (Kibbe W.A., "OligoCalc: an online oligonucleotide
//! properties calculator", Nucleic Acids Research 35:W43–W46, 2007): the anhydrous weight of
//! a single strand is `313.21 A + 304.2 T + 289.18 C + 329.21 G − 61.96` g/mol, the last term
//! taking off the 5' phosphate and adding the water of the free ends. A circular strand has no
//! ends, so nothing is taken off. An ambiguity code weighs the mean of the bases it stands for.

use serde::Serialize;

/// Residue weights of deoxynucleotides in a strand (OligoCalc), g/mol.
const A: f64 = 313.21;
const C: f64 = 289.18;
const G: f64 = 329.21;
const T: f64 = 304.2;
/// The 5' phosphate taken off (HPO3, 79.98) less the water of the free ends (18.02).
const ENDS: f64 = 61.96;

/// The mean weight of the bases an IUPAC code stands for; `None` for other bytes.
fn base_weight(b: u8) -> Option<f64> {
    let set: &[f64] = match b.to_ascii_uppercase() {
        b'A' => &[A],
        b'C' => &[C],
        b'G' => &[G],
        b'T' | b'U' => &[T],
        b'R' => &[A, G],
        b'Y' => &[C, T],
        b'S' => &[C, G],
        b'W' => &[A, T],
        b'K' => &[G, T],
        b'M' => &[A, C],
        b'B' => &[C, G, T],
        b'D' => &[A, G, T],
        b'H' => &[A, C, T],
        b'V' => &[A, C, G],
        b'N' => &[A, C, G, T],
        _ => return None,
    };
    Some(set.iter().sum::<f64>() / set.len() as f64)
}

/// The weight of the complementary base (mean for ambiguity codes).
fn complement_weight(b: u8) -> Option<f64> {
    let c = match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other, // S, W and N are their own complements.
    };
    base_weight(c)
}

/// Molecular weights of a stretch of DNA, g/mol (Da).
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct DnaWeight {
    /// This strand alone.
    pub single: f64,
    /// Both strands.
    pub double: f64,
    /// Some bases were ambiguity codes, weighed as the mean of what they stand for.
    pub approximate: bool,
}

/// The weight of `residues` (DNA, IUPAC codes), as a linear piece or a whole circle. `None`
/// when it is empty or holds a byte that isn't a nucleotide code.
pub fn dna_weight(residues: &[u8], circular: bool) -> Option<DnaWeight> {
    if residues.is_empty() {
        return None;
    }
    let (mut single, mut other) = (0.0, 0.0);
    let mut approximate = false;
    for &b in residues {
        single += base_weight(b)?;
        other += complement_weight(b)?;
        approximate |= !matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
    }
    if !circular {
        single -= ENDS;
        other -= ENDS;
    }
    Some(DnaWeight {
        single,
        double: single + other,
        approximate,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn worked_by_hand() {
        // ACGT: 313.21 + 289.18 + 329.21 + 304.2 − 61.96 = 1173.84; its complement ACGT too.
        let w = dna_weight(b"ACGT", false).unwrap();
        assert!((w.single - 1173.84).abs() < 1e-9);
        assert!((w.double - 2347.68).abs() < 1e-9);
        assert!(!w.approximate);
        // AAAA: 4 × 313.21 − 61.96 = 1190.88; TTTT: 4 × 304.2 − 61.96 = 1154.84.
        let w = dna_weight(b"aaaa", false).unwrap();
        assert!((w.single - 1190.88).abs() < 1e-9);
        assert!((w.double - (1190.88 + 1154.84)).abs() < 1e-9);
        // A circle has no ends.
        let w = dna_weight(b"ACGT", true).unwrap();
        assert!((w.single - 1235.8).abs() < 1e-9);
        // N is the mean of the four.
        let w = dna_weight(b"N", false).unwrap();
        assert!((w.single - (308.95 - 61.96)).abs() < 1e-9);
        assert!(w.approximate);
    }

    #[test]
    fn not_dna() {
        assert_eq!(dna_weight(b"", false), None);
        assert_eq!(dna_weight(b"ACGX", false), None);
        assert_eq!(dna_weight(b"AC-G", false), None);
    }
}
