//! Melting temperature of a DNA duplex by the nearest-neighbour method.
//!
//! Thermodynamics: the unified nearest-neighbour parameters, initiation terms and symmetry
//! correction of SantaLucia J. (1998) "A unified view of polymer, dumbbell, and oligonucleotide
//! DNA nearest-neighbor thermodynamics", PNAS 95:1460-1465 (Table 1). Salt: the same paper's
//! entropy correction, ΔS(Na+) = ΔS(1 M) + 0.368 · (N − 1) · ln[Na+], with Mg2+ (less what dNTPs
//! bind) converted to a sodium equivalent of 120·√([Mg2+] − [dNTP]) mM (von Ahsen N., Wittwer
//! C. T., Schütz E. (2001) Clin. Chem. 47:1956-1961). Two-state melting:
//! Tm = ΔH / (ΔS + R ln(CT/x)) with x = 4 for two different strands at equal concentration and
//! x = 1 for a self-complementary strand.

use serde::{Deserialize, Serialize};

/// Reaction conditions for a Tm.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct TmConditions {
    /// Oligo (strand) concentration, nM.
    pub oligo_nm: f64,
    /// Monovalent cations (Na+, K+), mM.
    pub na_mm: f64,
    /// Mg2+, mM.
    pub mg_mm: f64,
    /// dNTPs, mM (they bind Mg2+).
    pub dntp_mm: f64,
}

impl Default for TmConditions {
    /// 100 nM oligo in 100 mM monovalent salt, no Mg2+.
    fn default() -> Self {
        TmConditions {
            oligo_nm: 100.0,
            na_mm: 100.0,
            mg_mm: 0.0,
            dntp_mm: 0.0,
        }
    }
}

/// Gas constant, cal/(K·mol).
const R: f64 = 1.987;

/// ΔH (kcal/mol) and ΔS (cal/(K·mol)) of a nearest-neighbour pair, 5'->3' on one strand
/// (SantaLucia 1998, Table 1; each entry also covers its reverse complement).
pub(crate) fn nn_pair(a: u8, b: u8) -> Option<(f64, f64)> {
    Some(match (a, b) {
        (b'A', b'A') | (b'T', b'T') => (-7.9, -22.2),
        (b'A', b'T') => (-7.2, -20.4),
        (b'T', b'A') => (-7.2, -21.3),
        (b'C', b'A') | (b'T', b'G') => (-8.5, -22.7),
        (b'G', b'T') | (b'A', b'C') => (-8.4, -22.4),
        (b'C', b'T') | (b'A', b'G') => (-7.8, -21.0),
        (b'G', b'A') | (b'T', b'C') => (-8.2, -22.2),
        (b'C', b'G') => (-10.6, -27.2),
        (b'G', b'C') => (-9.8, -24.4),
        (b'G', b'G') | (b'C', b'C') => (-8.0, -19.9),
        _ => return None,
    })
}

/// Initiation for a terminal pair (SantaLucia 1998): G·C or A·T.
pub(crate) fn nn_init(b: u8) -> (f64, f64) {
    if matches!(b, b'G' | b'C') {
        (0.1, -2.8)
    } else {
        (2.3, 4.1)
    }
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        _ => b'C',
    }
}

/// Longest sequence we compute a Tm for (the two-state model is meant for oligos).
pub const MAX_TM_LEN: usize = 1000;

/// ΔH (kcal/mol) and ΔS at 1 M Na+ (cal/(K·mol)) of `seq` paired with its exact complement,
/// and whether it is self-complementary. `None` unless 2..=1000 unambiguous bases (A, C, G, T).
pub fn thermodynamics(seq: &[u8]) -> Option<(f64, f64, bool)> {
    let s: Vec<u8> = seq
        .iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'U' => b'T',
            other => other,
        })
        .collect();
    if s.len() < 2 || s.len() > MAX_TM_LEN || !s.iter().all(|b| b"ACGT".contains(b)) {
        return None;
    }
    let (mut dh, mut ds) = (0.0, 0.0);
    for w in s.windows(2) {
        let (h, e) = match w {
            [a, b] => nn_pair(*a, *b)?,
            _ => return None,
        };
        dh += h;
        ds += e;
    }
    for end in [s.first()?, s.last()?] {
        let (h, e) = nn_init(*end);
        dh += h;
        ds += e;
    }
    let self_comp = s.iter().rev().map(|&b| complement(b)).eq(s.iter().copied());
    if self_comp {
        ds += -1.4;
    }
    Some((dh, ds, self_comp))
}

/// The sodium-equivalent concentration (M) of the conditions.
pub(crate) fn sodium_equivalent(c: &TmConditions) -> f64 {
    let free_mg = (c.mg_mm - c.dntp_mm).max(0.0);
    (c.na_mm.max(0.0) + 120.0 * free_mg.sqrt()) / 1000.0
}

/// Melting temperature (°C) of `seq` with its exact complement, or `None` for sequences the
/// model doesn't cover (see [`thermodynamics`]) or conditions without salt or oligo.
pub fn tm(seq: &[u8], c: &TmConditions) -> Option<f64> {
    let (dh, ds, self_comp) = thermodynamics(seq)?;
    let na = sodium_equivalent(c);
    let ct = c.oligo_nm * 1e-9;
    if !(na > 0.0 && ct > 0.0 && na.is_finite() && ct.is_finite()) {
        return None;
    }
    let n = seq.len() as f64;
    let ds = ds + 0.368 * (n - 1.0) * na.ln();
    let x = if self_comp { 1.0 } else { 4.0 };
    let denominator = ds + R * (ct / x).ln();
    if denominator >= 0.0 {
        return None;
    }
    let t = dh * 1000.0 / denominator - 273.15;
    t.is_finite().then_some(t)
}

/// GC content as a fraction (0..=1) of the unambiguous bases; `None` if there are none.
pub fn gc_fraction(seq: &[u8]) -> Option<f64> {
    let (mut gc, mut n) = (0usize, 0usize);
    for b in seq {
        match b.to_ascii_uppercase() {
            b'G' | b'C' | b'S' => {
                gc += 1;
                n += 1;
            }
            b'A' | b'T' | b'U' | b'W' => n += 1,
            _ => {}
        }
    }
    (n > 0).then(|| gc as f64 / n as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: Option<f64>, b: f64) {
        let a = a.unwrap();
        assert!((a - b).abs() < 0.01, "{a} vs {b}");
    }

    /// Reference values from an independent implementation of the same model (Biopython 1.88
    /// `Tm_NN`, SantaLucia 1998 table, salt correction 5, equal strand concentrations).
    #[test]
    fn matches_reference_values() {
        let c = TmConditions::default();
        close(tm(b"CGTTGCATGCAAGCTT", &c), 52.2172);
        close(tm(b"GTAAAACGACGGCCAGTGCCA", &c), 62.2225);
        close(tm(b"ACGTTTGCA", &c), 24.7184);
        close(tm(b"agcgtaatctggaacgtc", &c), 53.1591);
        close(tm(b"ATATATATCCGG", &c), 25.2261);
        close(
            tm(b"TTGGCTTTTTTGTTAGACGAACCATGGCCGCCTCTAAACCTGTGGAGG", &c),
            74.6504,
        );
        let mg = TmConditions {
            mg_mm: 1.5,
            dntp_mm: 0.2,
            ..c
        };
        close(tm(b"CGTTGCATGCAAGCTT", &mg), 56.2492);
        let low = TmConditions { na_mm: 50.0, ..c };
        close(tm(b"CGTTGCATGCAAGCTT", &low), 49.0464);
    }

    #[test]
    fn by_hand() {
        // CG/GC alone: ΔH = -10.6 + 0.1 + 0.1, ΔS = -27.2 - 2.8 - 2.8 - 1.4 (self-complementary).
        let (dh, ds, sc) = thermodynamics(b"CG").unwrap();
        assert!(sc);
        assert!((dh - -10.4).abs() < 1e-9 && (ds - -34.2).abs() < 1e-9);
        // A strand and its reverse complement have the same thermodynamics.
        let a = thermodynamics(b"ACGGTTCAG").unwrap();
        let b = thermodynamics(b"CTGAACCGT").unwrap();
        assert!((a.0 - b.0).abs() < 1e-9 && (a.1 - b.1).abs() < 1e-9);
    }

    #[test]
    fn hostile_inputs() {
        let c = TmConditions::default();
        for s in [&b""[..], b"A", b"ACGN", b"AC-GT", &[0xff; 4]] {
            assert_eq!(tm(s, &c), None);
        }
        assert_eq!(tm(&vec![b'A'; MAX_TM_LEN + 1], &c), None);
        for bad in [
            TmConditions { na_mm: 0.0, ..c },
            TmConditions { oligo_nm: 0.0, ..c },
            TmConditions {
                oligo_nm: f64::NAN,
                ..c
            },
            TmConditions {
                na_mm: -5.0,
                mg_mm: -1.0,
                ..c
            },
            TmConditions {
                na_mm: f64::INFINITY,
                ..c
            },
        ] {
            assert!(tm(b"ACGTACGTAC", &bad).is_none_or(f64::is_finite));
        }
        assert_eq!(gc_fraction(b"GGCCAT"), Some(4.0 / 6.0));
        assert_eq!(gc_fraction(b"NNN"), None);
    }
}
