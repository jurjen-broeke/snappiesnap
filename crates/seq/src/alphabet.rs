use std::fmt;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

/// The kind of residues a sequence holds.
///
/// Nucleotide alphabets accept the IUPAC codes (NC-IUB 1984, *Nomenclature for incompletely
/// specified bases in nucleic acid sequences*); proteins accept the one-letter amino-acid codes
/// including B, J, O, U, X, Z and `*` for a stop. Case is accepted and preserved.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "lowercase")]
pub enum Alphabet {
    Dna,
    Rna,
    Protein,
}

impl Alphabet {
    /// Whether `byte` is a valid residue in this alphabet.
    pub fn is_valid(self, byte: u8) -> bool {
        let up = byte.to_ascii_uppercase();
        match self {
            Alphabet::Dna => is_iupac_ambiguity(up) || matches!(up, b'A' | b'C' | b'G' | b'T'),
            Alphabet::Rna => is_iupac_ambiguity(up) || matches!(up, b'A' | b'C' | b'G' | b'U'),
            Alphabet::Protein => up.is_ascii_uppercase() || byte == b'*',
        }
    }

    pub fn is_nucleotide(self) -> bool {
        matches!(self, Alphabet::Dna | Alphabet::Rna)
    }

    /// The IUPAC complement of a nucleotide, preserving case; `None` for proteins or bytes
    /// that aren't valid in this alphabet.
    pub fn complement(self, byte: u8) -> Option<u8> {
        if !self.is_nucleotide() || !self.is_valid(byte) {
            return None;
        }
        let up = byte.to_ascii_uppercase();
        let comp = match up {
            b'A' if self == Alphabet::Rna => b'U',
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
            // S, W and N are their own complements.
            other => other,
        };
        Some(if byte.is_ascii_lowercase() {
            comp.to_ascii_lowercase()
        } else {
            comp
        })
    }
}

fn is_iupac_ambiguity(up: u8) -> bool {
    matches!(
        up,
        b'R' | b'Y' | b'S' | b'W' | b'K' | b'M' | b'B' | b'D' | b'H' | b'V' | b'N'
    )
}

impl fmt::Display for Alphabet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Alphabet::Dna => "DNA",
            Alphabet::Rna => "RNA",
            Alphabet::Protein => "protein",
        })
    }
}

#[cfg(test)]
mod tests {
    use super::Alphabet::{self, *};

    #[test]
    fn validity() {
        for b in b"ACGTRYSWKMBDHVNacgtn" {
            assert!(Dna.is_valid(*b), "{}", *b as char);
        }
        assert!(!Dna.is_valid(b'U'));
        assert!(Rna.is_valid(b'U') && !Rna.is_valid(b'T'));
        assert!(Protein.is_valid(b'*') && Protein.is_valid(b'w'));
        for bad in [b'-', b' ', b'1', 0, 255] {
            for a in [Dna, Rna, Protein] {
                assert!(!a.is_valid(bad));
            }
        }
    }

    #[test]
    fn complement_is_an_involution_on_every_valid_byte() {
        for a in [Dna, Rna] {
            for b in 0..=255u8 {
                if let Some(c) = a.complement(b) {
                    assert_eq!(a.complement(c), Some(b), "{a} {}", b as char);
                }
            }
        }
        assert_eq!(Dna.complement(b'a'), Some(b't'));
        assert_eq!(Rna.complement(b'A'), Some(b'U'));
        assert_eq!(Alphabet::Protein.complement(b'A'), None);
    }
}
