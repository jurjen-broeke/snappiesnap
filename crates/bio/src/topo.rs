//! TA and TOPO cloning: a PCR product put into a vector that is opened at one place.
//!
//! - **TA** (Marchuk D. et al. (1991) Nucleic Acids Res. 19:1154; Holton T. A., Graham M. W.
//!   (1991) Nucleic Acids Res. 19:1156): *Taq* polymerase leaves a single 3' A on each end of a
//!   PCR product, and a T-vector has a single 3' T on each end. They anneal and ligate in either
//!   orientation; read along the top strand the construct is the vector's left arm, `T`, the
//!   insert, `A`, the right arm. TOPO TA vectors (vaccinia topoisomerase I-charged; Shuman S.
//!   (1994) J. Biol. Chem. 269:32678-32684) have the same ends and give the same construct.
//! - **Blunt TOPO**: a blunt PCR product goes in as it is, in either orientation.
//! - **Directional TOPO**: the vector's 5' `GTGG` overhang anneals to `CACC` at the start of the
//!   insert (put there by the forward primer), so only an orientation in which the insert starts
//!   with `CACC` is cloned.
//!
//! A circular vector is opened at a given position (between two bases); a linear vector is
//! taken as already opened, its two ends joining the insert. The construct keeps the vector's
//! numbering; vector features cut by the opening are clipped and noted as partial.

use serde::Serialize;
use snappiesnap_doc::{Molecule, Topology};

use crate::cloning::{Part, assemble, extract};

/// Which kind of TA/TOPO cloning.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum TopoKind {
    /// T-vector (or TOPO TA vector) and a PCR product with 3' A overhangs.
    Ta,
    /// Blunt TOPO: a blunt PCR product.
    Blunt,
    /// Directional TOPO: the insert must start with `CACC`.
    Directional,
}

impl TopoKind {
    pub fn label(self) -> &'static str {
        match self {
            TopoKind::Ta => "TA cloning",
            TopoKind::Blunt => "blunt TOPO cloning",
            TopoKind::Directional => "directional TOPO cloning",
        }
    }
}

/// One possible construct.
#[derive(Debug, Clone)]
pub struct TopoClone {
    pub molecule: Molecule,
    /// Whether the insert went in turned around (reverse complement).
    pub reversed: bool,
    pub description: String,
}

/// The overhang the directional vector anneals to, at the start of the insert.
pub const DIRECTIONAL_START: &[u8] = b"CACC";

/// The constructs `kind` cloning of `insert` into `vector` gives (forward orientation first).
/// `position` is where a circular vector is opened (0..=length, between bases); a linear vector
/// takes none.
pub fn topo_clone(
    vector: &Molecule,
    position: Option<usize>,
    insert: &Molecule,
    kind: TopoKind,
    name: &str,
) -> Result<Vec<TopoClone>, String> {
    let len = vector.len();
    if len == 0 {
        return Err("the vector is empty".into());
    }
    if insert.is_empty() {
        return Err("the insert is empty".into());
    }
    if insert.topology() == Topology::Circular {
        return Err("the insert must be a linear fragment, such as a PCR product".into());
    }
    let circular = vector.topology() == Topology::Circular;
    // The opened vector as one linear arm, from the right of the opening round to its left.
    let (arm, origin) = match (circular, position) {
        (true, None) => {
            return Err("give the position where the circular vector is opened".into());
        }
        (true, Some(p)) if p > len => {
            return Err(format!("position {p} is outside the vector (0 to {len})"));
        }
        (true, Some(p)) => (extract(vector, p % len, len), (len - p % len) % len),
        (false, None) => (extract(vector, 0, len), 0),
        (false, Some(p)) if p == 0 || p == len => (extract(vector, 0, len), 0),
        (false, Some(_)) => {
            return Err("a linear vector is taken as already opened: give no position".into());
        }
    };
    let whole = extract(insert, 0, insert.len());
    let base = |b: &[u8]| Part {
        seq: b.to_vec(),
        features: Vec::new(),
    };
    let mut out = Vec::new();
    for (reversed, piece) in [(false, whole.clone()), (true, whole.reversed())] {
        let starts = piece
            .seq
            .get(..DIRECTIONAL_START.len())
            .is_some_and(|s| s.eq_ignore_ascii_case(DIRECTIONAL_START));
        if kind == TopoKind::Directional && !starts {
            continue;
        }
        let parts = match kind {
            TopoKind::Ta => vec![arm.clone(), base(b"T"), piece, base(b"A")],
            TopoKind::Blunt | TopoKind::Directional => vec![arm.clone(), piece],
        };
        let mut m = assemble(&parts, Topology::Circular, name).map_err(|e| e.to_string())?;
        if origin != 0 {
            m.rotate_origin(origin).map_err(|e| e.to_string())?;
        }
        let orientation = if reversed { "reversed" } else { "forward" };
        let at = match position {
            Some(p) if circular => format!(" at {p}"),
            _ => String::new(),
        };
        out.push(TopoClone {
            molecule: m,
            reversed,
            description: format!("{}{at}, insert {orientation}", kind.label()),
        });
    }
    if out.is_empty() {
        return Err(format!(
            "for directional TOPO cloning the insert must start with {} (from the forward primer) on one strand",
            String::from_utf8_lossy(DIRECTIONAL_START)
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_doc::{Alphabet, NewFeature, Range, Segment, Seq, Strand};

    fn mol(name: &str, s: &str, topology: Topology) -> Molecule {
        Molecule::new(name, Seq::new(Alphabet::Dna, s).unwrap(), topology).unwrap()
    }

    fn residues(m: &Molecule) -> String {
        String::from_utf8_lossy(m.seq().as_bytes()).into_owned()
    }

    fn feature(name: &str, start: usize, end: usize) -> NewFeature {
        NewFeature {
            color: Some("#123456".into()),
            name: name.into(),
            kind: "misc_feature".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(start, end),
            }],
            qualifiers: Vec::new(),
        }
    }

    #[test]
    fn ta_cloning_adds_t_and_a_in_both_orientations() {
        let mut vector = mol("v", "GGGGGCCCCC", Topology::Circular);
        vector.add_feature(feature("left", 0, 3)).unwrap();
        vector.add_feature(feature("right", 6, 9)).unwrap();
        let mut insert = mol("i", "ACGTTG", Topology::Linear);
        insert.add_feature(feature("gene", 0, 4)).unwrap();
        let c = topo_clone(&vector, Some(5), &insert, TopoKind::Ta, "c").unwrap();
        assert_eq!(c.len(), 2);
        // Vector numbering kept: left arm, T, insert, A, right arm.
        assert_eq!(residues(&c[0].molecule), "GGGGGTACGTTGACCCCC");
        assert_eq!(residues(&c[1].molecule), "GGGGGTCAACGTACCCCC");
        assert!(!c[0].reversed && c[1].reversed);
        assert_eq!(c[0].molecule.topology(), Topology::Circular);
        let at = |m: &Molecule, n: &str| {
            m.features()
                .iter()
                .find(|f| f.name == n)
                .map(|f| (f.segments[0].range, f.strand, f.color.clone()))
        };
        // Features follow: the vector's stay put or shift by the insert; the insert's sit after
        // the T, turned with it.
        assert_eq!(at(&c[0].molecule, "left").unwrap().0, Range::new(0, 3));
        assert_eq!(at(&c[0].molecule, "right").unwrap().0, Range::new(14, 17));
        assert_eq!(at(&c[0].molecule, "gene").unwrap().0, Range::new(6, 10));
        let (r, strand, color) = at(&c[1].molecule, "gene").unwrap();
        assert_eq!((r, strand), (Range::new(8, 12), Strand::Reverse));
        assert_eq!(color.as_deref(), Some("#123456"));
    }

    #[test]
    fn opening_across_the_origin_and_at_its_ends() {
        let vector = mol("v", "AAAACCCC", Topology::Circular);
        let insert = mol("i", "GG", Topology::Linear);
        for p in [0, 8] {
            let c = topo_clone(&vector, Some(p), &insert, TopoKind::Blunt, "c").unwrap();
            assert_eq!(residues(&c[0].molecule), "AAAACCCCGG", "position {p}");
        }
        let c = topo_clone(&vector, Some(1), &insert, TopoKind::Blunt, "c").unwrap();
        assert_eq!(residues(&c[0].molecule), "AGGAAACCCC");
        assert_eq!(residues(&c[1].molecule), "ACCAAACCCC");
    }

    #[test]
    fn a_linear_vector_is_already_open() {
        let vector = mol("v", "AAAACCCC", Topology::Linear);
        let insert = mol("i", "GTT", Topology::Linear);
        let c = topo_clone(&vector, None, &insert, TopoKind::Ta, "c").unwrap();
        assert_eq!(residues(&c[0].molecule), "AAAACCCCTGTTA");
        assert_eq!(c[0].molecule.topology(), Topology::Circular);
        assert!(topo_clone(&vector, Some(3), &insert, TopoKind::Ta, "c").is_err());
        assert!(topo_clone(&vector, Some(8), &insert, TopoKind::Ta, "c").is_ok());
    }

    #[test]
    fn directional_topo_needs_cacc() {
        let vector = mol("v", "AAAATTTT", Topology::Circular);
        let insert = mol("i", "CACCATGTAA", Topology::Linear);
        let c = topo_clone(&vector, Some(4), &insert, TopoKind::Directional, "c").unwrap();
        assert_eq!(c.len(), 1);
        assert_eq!(residues(&c[0].molecule), "AAAACACCATGTAATTTT");
        // An insert that also ends with GGTG fits both ways.
        let both = mol("i", "CACCATGGGTG", Topology::Linear);
        let c = topo_clone(&vector, Some(4), &both, TopoKind::Directional, "c").unwrap();
        assert_eq!(c.len(), 2);
        let none = mol("i", "ATGTAA", Topology::Linear);
        let e = topo_clone(&vector, Some(4), &none, TopoKind::Directional, "c").unwrap_err();
        assert!(e.contains("CACC"));
    }

    #[test]
    fn bad_input_is_an_error() {
        let v = mol("v", "AAAA", Topology::Circular);
        let i = mol("i", "GG", Topology::Linear);
        let empty =
            Molecule::new("e", Seq::new(Alphabet::Dna, "").unwrap(), Topology::Linear).unwrap();
        let ring = mol("r", "GG", Topology::Circular);
        assert!(topo_clone(&v, None, &i, TopoKind::Ta, "c").is_err());
        assert!(topo_clone(&v, Some(5), &i, TopoKind::Ta, "c").is_err());
        assert!(topo_clone(&v, Some(2), &empty, TopoKind::Ta, "c").is_err());
        assert!(topo_clone(&empty, None, &i, TopoKind::Ta, "c").is_err());
        assert!(topo_clone(&v, Some(2), &ring, TopoKind::Ta, "c").is_err());
        // One base each.
        let one = mol("o", "A", Topology::Circular);
        let c = topo_clone(
            &one,
            Some(1),
            &mol("i", "G", Topology::Linear),
            TopoKind::Ta,
            "c",
        )
        .unwrap();
        assert_eq!(residues(&c[0].molecule), "ATGA");
    }
}
