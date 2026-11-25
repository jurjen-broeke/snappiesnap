//! Translation, open reading frames and protein properties.
//!
//! Codes are NCBI's genetic code tables ([`crate::codes`]). An ambiguous codon (IUPAC codes)
//! translates to an amino acid only if every codon it stands for gives that amino acid, and to
//! `X` otherwise.
//!
//! Protein properties:
//! - Molecular weight: sum of average residue masses plus one water (average isotopic masses as
//!   tabulated by ExPASy ProtParam / IUPAC).
//! - Isoelectric point: the pH where the net charge (Henderson-Hasselbalch) is zero, with the pK
//!   values of Bjellqvist B. et al. (1993) Electrophoresis 14:1023-1031 for the termini and the
//!   ionisable side chains (without their position-specific corrections).
//! - Extinction coefficient at 280 nm in water: 5500 per Trp, 1490 per Tyr and 125 per cystine
//!   (Pace C. N. et al. (1995) Protein Sci. 4:2411-2423).

use serde::Serialize;
use snappiesnap_doc::{Molecule, Range, Strand, Topology};

use crate::codes::{CODES, GeneticCode};

/// The standard code (NCBI table 1).
pub const STANDARD: u8 = 1;

/// An NCBI genetic code by table number.
pub fn code(id: u8) -> Option<&'static GeneticCode> {
    CODES.iter().find(|c| c.id == id)
}

/// The bases an IUPAC nucleotide code stands for, in T, C, A, G order (as NCBI indexes codons).
fn bases(b: u8) -> &'static [usize] {
    match b.to_ascii_uppercase() {
        b'T' | b'U' => &[0],
        b'C' => &[1],
        b'A' => &[2],
        b'G' => &[3],
        b'Y' => &[0, 1],
        b'R' => &[2, 3],
        b'W' => &[0, 2],
        b'S' => &[1, 3],
        b'K' => &[0, 3],
        b'M' => &[1, 2],
        b'B' => &[0, 1, 3],
        b'D' => &[0, 2, 3],
        b'H' => &[0, 1, 2],
        b'V' => &[1, 2, 3],
        b'N' => &[0, 1, 2, 3],
        _ => &[],
    }
}

/// Look a codon up in one of the code's 64-entry strings; `None` unless all readings agree.
fn lookup(codon: [u8; 3], table: &[u8; 64]) -> Option<u8> {
    let mut found = None;
    for &a in bases(codon[0]) {
        for &b in bases(codon[1]) {
            for &c in bases(codon[2]) {
                let v = *table.get(a * 16 + b * 4 + c)?;
                match found {
                    None => found = Some(v),
                    Some(f) if f != v => return None,
                    _ => {}
                }
            }
        }
    }
    found
}

impl GeneticCode {
    /// The amino acid of a codon (`*` for stop, `X` if ambiguous or not a codon).
    pub fn amino_acid(&self, codon: [u8; 3]) -> u8 {
        lookup(codon, self.amino_acids).unwrap_or(b'X')
    }

    pub fn is_stop(&self, codon: [u8; 3]) -> bool {
        lookup(codon, self.amino_acids) == Some(b'*')
    }

    /// Whether the codon can start translation in this code. With `only_atg`, just ATG.
    pub fn is_start(&self, codon: [u8; 3], only_atg: bool) -> bool {
        if only_atg {
            return codon.iter().map(u8::to_ascii_uppercase).eq(*b"ATG");
        }
        lookup(codon, self.starts) == Some(b'M')
    }
}

/// Translate `seq` (5'->3') codon by codon; a trailing partial codon is left out. With
/// `as_start`, a first codon that can initiate is read as Met (as at a real start).
pub fn translate(seq: &[u8], code: &GeneticCode, as_start: bool) -> Vec<u8> {
    let mut out: Vec<u8> = seq
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| code.amino_acid(*c))
        .collect();
    if as_start
        && let (Some(first), Some(a), Some(b), Some(c)) =
            (out.first_mut(), seq.first(), seq.get(1), seq.get(2))
        && code.is_start([*a, *b, *c], false)
    {
        *first = b'M';
    }
    out
}

fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
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
        o => o,
    }
}

/// The residues of `range` on `m`, read 5'->3' on `strand` (reverse: the reverse complement).
pub fn strand_residues(m: &Molecule, range: Range, strand: Strand) -> Option<Vec<u8>> {
    let r = m.residues(range).ok()?;
    Some(if strand == Strand::Reverse {
        r.iter().rev().map(|&b| complement(b)).collect()
    } else {
        r
    })
}

/// What to look for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, serde::Deserialize)]
#[serde(default)]
pub struct OrfOptions {
    /// Shortest ORF, in amino acids (not counting the stop).
    pub min_aa: usize,
    /// NCBI genetic code.
    pub code: u8,
    /// Only ATG starts (otherwise every start codon of the code).
    pub only_atg: bool,
    /// Search the reverse strand too.
    pub both_strands: bool,
}

impl Default for OrfOptions {
    fn default() -> Self {
        OrfOptions {
            min_aa: 75,
            code: STANDARD,
            only_atg: true,
            both_strands: true,
        }
    }
}

/// An open reading frame: from a start codon to (and including) a stop codon.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Orf {
    /// On the top strand; wraps on circular molecules.
    pub range: Range,
    pub strand: Strand,
    /// Reading frame on its strand, 1-3.
    pub frame: u8,
    /// Amino acids (not counting the stop).
    pub length: usize,
}

/// Most ORFs reported (an absurdly low minimum on a big genome would make millions).
pub const MAX_ORFS: usize = 100_000;

/// The ORFs of `m`, longest first. For each stop codon only the longest ORF (from the first
/// start after the previous stop in that frame) is reported. On circular molecules ORFs may
/// run across the origin; a frame without any stop has no ORF.
pub fn find_orfs(m: &Molecule, opts: &OrfOptions) -> Vec<Orf> {
    let Some(code) = code(opts.code) else {
        return Vec::new();
    };
    let len = m.len();
    let circular = m.topology() == Topology::Circular;
    let top: Vec<u8> = m.seq().as_bytes().to_vec();
    let mut out = Vec::new();
    if len < 3 || !m.seq().alphabet().is_nucleotide() {
        return out;
    }
    let rev: Vec<u8> = top.iter().rev().map(|&b| complement(b)).collect();
    let mut strands = vec![(Strand::Forward, &top)];
    if opts.both_strands {
        strands.push((Strand::Reverse, &rev));
    }
    let min_aa = opts.min_aa.max(1);
    // On a circle read twice round, so ORFs across the origin are found; on a line, once.
    let reach = if circular { 2 * len } else { len };
    for (strand, s) in strands {
        let codon_at = |p: usize| -> Option<[u8; 3]> {
            let get = |i: usize| s.get(if circular { i % len } else { i }).copied();
            Some([get(p)?, get(p + 1)?, get(p + 2)?])
        };
        for frame in 0..3 {
            let mut open: Option<usize> = None;
            let mut p = frame;
            while p + 3 <= reach {
                let Some(codon) = codon_at(p) else { break };
                if code.is_stop(codon) {
                    if let Some(start) = open.take() {
                        let aa = (p - start) / 3;
                        // Starts in the second lap repeat ORFs already found.
                        if aa >= min_aa && start < len && p + 3 - start <= len {
                            let (a, b) = (start, p + 3);
                            let range = match strand {
                                Strand::Reverse => {
                                    Range::new((2 * len - b) % len, (2 * len - a) % len)
                                }
                                _ => Range::new(a % len, b % len),
                            };
                            let range = if !circular && range.end == 0 {
                                Range::new(range.start, len)
                            } else {
                                range
                            };
                            out.push(Orf {
                                range,
                                strand,
                                frame: u8::try_from(start % 3 + 1).unwrap_or(1),
                                length: aa,
                            });
                            if out.len() >= MAX_ORFS {
                                return out;
                            }
                        }
                    }
                } else if open.is_none() && code.is_start(codon, opts.only_atg) {
                    open = Some(p);
                }
                p += 3;
                if open.is_none() && p >= len && circular {
                    // Past the origin with nothing open: anything new repeats the first lap.
                    break;
                }
            }
        }
    }
    out.sort_by(|a, b| {
        b.length
            .cmp(&a.length)
            .then(a.range.start.cmp(&b.range.start))
    });
    out
}

/// Average residue masses (Da) of the 20 amino acids, plus selenocysteine and pyrrolysine.
fn residue_mass(aa: u8) -> Option<f64> {
    Some(match aa.to_ascii_uppercase() {
        b'A' => 71.0788,
        b'R' => 156.1875,
        b'N' => 114.1038,
        b'D' => 115.0886,
        b'C' => 103.1388,
        b'E' => 129.1155,
        b'Q' => 128.1307,
        b'G' => 57.0519,
        b'H' => 137.1411,
        b'I' => 113.1594,
        b'L' => 113.1594,
        b'K' => 128.1741,
        b'M' => 131.1926,
        b'F' => 147.1766,
        b'P' => 97.1167,
        b'S' => 87.0782,
        b'T' => 101.1051,
        b'W' => 186.2132,
        b'Y' => 163.1760,
        b'V' => 99.1326,
        b'U' => 150.0388,
        b'O' => 237.3018,
        _ => return None,
    })
}

const WATER: f64 = 18.01524;

/// Properties of a protein sequence (one-letter codes; a trailing `*` is ignored).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ProteinProperties {
    pub length: usize,
    /// Average mass, Da. Residues without a mass (X, B, Z) are left out and counted below.
    pub mass: f64,
    pub unknown_residues: usize,
    pub isoelectric_point: f64,
    /// ε280 (M⁻¹ cm⁻¹) assuming all Cys form cystines, and with all Cys reduced.
    pub extinction_oxidised: u32,
    pub extinction_reduced: u32,
}

fn charge(protein: &[u8], ph: f64) -> f64 {
    // Bjellqvist et al. 1993 pK values.
    let pos = |pk: f64| 1.0 / (1.0 + 10f64.powf(ph - pk));
    let neg = |pk: f64| -1.0 / (1.0 + 10f64.powf(pk - ph));
    let mut q = pos(7.5) + neg(3.55);
    for &aa in protein {
        q += match aa.to_ascii_uppercase() {
            b'K' => pos(10.0),
            b'R' => pos(12.0),
            b'H' => pos(5.98),
            b'D' => neg(4.05),
            b'E' => neg(4.45),
            b'C' => neg(9.0),
            b'Y' => neg(10.0),
            _ => 0.0,
        };
    }
    q
}

pub fn protein_properties(protein: &[u8]) -> ProteinProperties {
    let p: Vec<u8> = protein.iter().copied().filter(|&b| b != b'*').collect();
    let mut mass = if p.is_empty() { 0.0 } else { WATER };
    let mut unknown = 0;
    for &aa in &p {
        match residue_mass(aa) {
            Some(m) => mass += m,
            None => unknown += 1,
        }
    }
    // Net charge falls with pH: bisect for zero.
    let (mut lo, mut hi) = (0.0f64, 14.0f64);
    for _ in 0..60 {
        let mid = (lo + hi) / 2.0;
        if charge(&p, mid) > 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let count = |c: u8| p.iter().filter(|&&x| x.to_ascii_uppercase() == c).count() as u32;
    let (w, y, c) = (count(b'W'), count(b'Y'), count(b'C'));
    let reduced = 5500 * w + 1490 * y;
    ProteinProperties {
        length: p.len(),
        mass,
        unknown_residues: unknown,
        isoelectric_point: (lo + hi) / 2.0,
        extinction_oxidised: reduced + 125 * (c / 2),
        extinction_reduced: reduced,
    }
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, Seq};

    use super::*;

    fn mol(s: &str, t: Topology) -> Molecule {
        Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), t).unwrap()
    }

    #[test]
    fn codes_and_codons() {
        assert_eq!(CODES.len(), 27);
        let std = code(1).unwrap();
        assert_eq!(translate(b"ATGGCCTAA", std, false), b"MA*");
        assert_eq!(
            translate(b"atggcctaag", std, false),
            b"MA*",
            "partial codon left out"
        );
        // Ambiguity: GCN is always Ala; RAY (AAT/AAC/GAT/GAC) isn't one amino acid.
        assert_eq!(translate(b"GCNRAY", std, false), b"AX");
        // Mitochondrial code: TGA is Trp, AGA a stop.
        let mito = code(2).unwrap();
        assert_eq!(translate(b"TGAAGA", mito, false), b"W*");
        // Bacterial code: GTG starts, read as Met at a start.
        let bact = code(11).unwrap();
        assert!(bact.is_start(*b"GTG", false) && !bact.is_start(*b"GTG", true));
        assert_eq!(translate(b"GTGAAA", bact, true), b"MK");
        assert_eq!(translate(b"GTGAAA", bact, false), b"VK");
        assert!(code(7).is_none());
    }

    #[test]
    fn orfs_on_both_strands_and_across_the_origin() {
        // ATG + 4 codons + TAA forward; its reverse complement placed later reads on the bottom.
        let orf = "ATGAAACCCGGGTTTTAA";
        let rc: String = orf.bytes().rev().map(|b| complement(b) as char).collect();
        let s = format!("CC{orf}CCCC{rc}CC");
        let m = mol(&s, Topology::Linear);
        let opts = OrfOptions {
            min_aa: 5,
            ..OrfOptions::default()
        };
        let found = find_orfs(&m, &opts);
        assert_eq!(found.len(), 2, "{found:?}");
        let f = found.iter().find(|o| o.strand == Strand::Forward).unwrap();
        assert_eq!((f.range, f.length, f.frame), (Range::new(2, 20), 5, 3));
        let r = found.iter().find(|o| o.strand == Strand::Reverse).unwrap();
        assert_eq!((r.range, r.length), (Range::new(24, 42), 5));
        // Rotate so the forward ORF spans the origin.
        let rotated = format!("{}{}", &s[10..], &s[..10]);
        let m = mol(&rotated, Topology::Circular);
        let found = find_orfs(&m, &opts);
        let f = found.iter().find(|o| o.strand == Strand::Forward).unwrap();
        assert_eq!(f.range, Range::new(rotated.len() - 8, 10));
        assert_eq!(found.len(), 2);
        // Linear: the split ORF is gone.
        assert_eq!(find_orfs(&mol(&rotated, Topology::Linear), &opts).len(), 1);
        // The minimum length applies.
        assert!(
            find_orfs(
                &mol(&s, Topology::Linear),
                &OrfOptions { min_aa: 6, ..opts }
            )
            .is_empty()
        );
    }

    #[test]
    fn longest_orf_per_stop() {
        // Two starts in frame before one stop: one ORF from the first.
        let m = mol("ATGAAAATGAAAAAATAA", Topology::Linear);
        let o = find_orfs(
            &m,
            &OrfOptions {
                min_aa: 1,
                both_strands: false,
                ..OrfOptions::default()
            },
        );
        assert_eq!(o.len(), 1);
        assert_eq!(o[0].length, 5);
    }

    /// Reference values from an independent implementation (Biopython 1.88 ProteinAnalysis:
    /// molecular_weight, and isoelectric_point for a sequence whose ends get no special pK).
    #[test]
    fn protein_properties_match_reference_values() {
        let p = protein_properties(b"GSHMWKYCLLCRDEAAAG*");
        assert_eq!(p.length, 18);
        assert!((p.mass - REF_MASS).abs() < 0.05, "{}", p.mass);
        assert!(
            (p.isoelectric_point - REF_PI).abs() < 0.01,
            "{}",
            p.isoelectric_point
        );
        assert_eq!(
            (p.extinction_reduced, p.extinction_oxidised),
            (5500 + 1490, 5500 + 1490 + 125)
        );
        let empty = protein_properties(b"");
        assert_eq!((empty.length, empty.mass), (0, 0.0));
        assert_eq!(protein_properties(b"AXB").unknown_residues, 2);
    }

    const REF_MASS: f64 = 2011.3077;
    const REF_PI: f64 = 6.7332;

    #[test]
    fn hostile() {
        for s in ["", "A", "NN", "NNNNNNNNNNNN", "ATG"] {
            for t in [Topology::Linear, Topology::Circular] {
                let m = mol(s, t);
                for min_aa in [0, 1, 1000] {
                    let _ = find_orfs(
                        &m,
                        &OrfOptions {
                            min_aa,
                            ..OrfOptions::default()
                        },
                    );
                }
            }
        }
        assert!(
            find_orfs(
                &mol("ATGAAATAA", Topology::Linear),
                &OrfOptions {
                    code: 99,
                    ..OrfOptions::default()
                }
            )
            .is_empty()
        );
        let _ = protein_properties(&[0xff, b'*', b'*']);
    }
}
