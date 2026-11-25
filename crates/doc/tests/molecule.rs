//! Molecule edits keep annotations attached to their residues and the document valid.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

use proptest::prelude::*;
use snappiesnap_doc::*;

fn plasmid() -> Molecule {
    // 20 bp circular, one feature spanning the origin and one in the middle.
    let seq = Seq::new(Alphabet::Dna, "AAAACCCCGGGGTTTTACGT").unwrap();
    let mut m = Molecule::new("pTest", seq, Topology::Circular).unwrap();
    m.add_feature(feature("ori-span", 16, 4)).unwrap();
    m.add_feature(feature("mid", 8, 12)).unwrap();
    m
}

fn feature(name: &str, start: usize, end: usize) -> NewFeature {
    NewFeature {
        color: None,
        name: name.into(),
        kind: "misc_feature".into(),
        strand: Strand::Forward,
        segments: vec![Segment {
            range: Range::new(start, end),
        }],
        qualifiers: vec![],
    }
}

fn residues_of(m: &Molecule, name: &str) -> String {
    let f = m.features().iter().find(|f| f.name == name).unwrap();
    let r = f.segments[0].range;
    String::from_utf8(m.residues(r).unwrap()).unwrap()
}

#[test]
fn features_follow_their_residues_through_edits() {
    let mut m = plasmid();
    assert_eq!(residues_of(&m, "ori-span"), "ACGTAAAA");
    assert_eq!(residues_of(&m, "mid"), "GGGG");

    m.insert(10, b"tt").unwrap();
    assert_eq!(residues_of(&m, "mid"), "GGttGG");

    m.rotate_origin(9).unwrap();
    assert_eq!(residues_of(&m, "ori-span"), "ACGTAAAA");
    assert_eq!(residues_of(&m, "mid"), "GGttGG");

    m.reverse_complement().unwrap();
    assert_eq!(residues_of(&m, "mid"), "CCaaCC");
    assert_eq!(residues_of(&m, "ori-span"), "TTTTACGT");
    assert!(m.features().iter().all(|f| f.strand == Strand::Reverse));
    m.validate().unwrap();
}

#[test]
fn deleting_everything_a_feature_covers_removes_it() {
    let mut m = plasmid();
    let report = m.delete(Range::new(7, 13)).unwrap();
    assert_eq!(report.removed_features, vec![FeatureId(2)]);
    assert_eq!(m.features().len(), 1);
    // A deletion across the origin trims the wrapping feature.
    m.delete(Range::new(12, 2)).unwrap();
    // Before: ori-span = ACGT|AAAA on 14 bp; [12,2) removes "GT" and "AA".
    assert_eq!(residues_of(&m, "ori-span"), "ACAA");
    m.validate().unwrap();
}

#[test]
fn linearising_needs_no_feature_across_the_origin() {
    let mut m = plasmid();
    assert!(matches!(
        m.set_topology(Topology::Linear),
        Err(DocError::FeatureSpansOrigin { .. })
    ));
    m.rotate_origin(16).unwrap();
    m.set_topology(Topology::Linear).unwrap();
    assert!(m.rotate_origin(1).is_err());
}

#[test]
fn failed_replace_changes_nothing() {
    let mut m = plasmid();
    let before = m.clone();
    assert!(m.replace(Range::new(2, 6), b"AXC").is_err());
    assert!(m.replace(Range::new(2, 60), b"A").is_err());
    assert_eq!(m, before);
    m.replace(Range::new(18, 2), b"gg").unwrap();
    assert_eq!(m.len(), 18);
    m.validate().unwrap();
}

#[test]
fn invalid_features_are_rejected() {
    let mut m = plasmid();
    // Unnamed features are allowed: many files have them.
    m.add_feature(feature("", 0, 2)).unwrap();
    for bad in [
        feature("empty", 3, 3),
        feature("oob", 0, 21),
        NewFeature {
            segments: vec![],
            ..feature("none", 0, 1)
        },
        NewFeature {
            kind: " ".into(),
            ..feature("kind", 0, 1)
        },
    ] {
        let name = bad.name.clone();
        assert!(m.add_feature(bad).is_err(), "{name}");
    }
    let mut linear = Molecule::new("lin", m.seq().clone(), Topology::Linear).unwrap();
    assert!(linear.add_feature(feature("wraps", 16, 4)).is_err());
}

#[test]
fn serde_round_trip_and_validation_on_load() {
    let m = plasmid();
    let json = serde_json::to_value(&m).unwrap();
    assert_eq!(json["length"], 20);
    let back: Molecule = serde_json::from_value(json.clone()).unwrap();
    assert_eq!(back, m);

    let mut broken = json;
    broken["features"][0]["segments"][0]["range"]["end"] = 99.into();
    assert!(serde_json::from_value::<Molecule>(broken).is_err());
}

#[derive(Debug, Clone)]
enum Op {
    Insert(usize, String),
    Delete(usize, usize),
    Replace(usize, usize, String),
    Rotate(usize),
    RevComp,
    Linear,
    Circular,
    AddFeature(usize, usize),
}

fn op() -> impl Strategy<Value = Op> {
    let residues = "[ACGTNacgt]{0,6}";
    prop_oneof![
        (0usize..40, residues).prop_map(|(a, s)| Op::Insert(a, s)),
        (0usize..40, 0usize..40).prop_map(|(a, b)| Op::Delete(a, b)),
        (0usize..40, 0usize..40, residues).prop_map(|(a, b, s)| Op::Replace(a, b, s)),
        (0usize..40).prop_map(Op::Rotate),
        Just(Op::RevComp),
        Just(Op::Linear),
        Just(Op::Circular),
        (0usize..40, 0usize..40).prop_map(|(a, b)| Op::AddFeature(a, b)),
    ]
}

proptest! {
    /// Any sequence of edits, valid or not, keeps the document valid and never panics.
    #[test]
    fn random_edits_keep_the_document_valid(ops in proptest::collection::vec(op(), 0..30)) {
        let mut m = plasmid();
        for op in ops {
            let before = m.clone();
            let ok = match op {
                Op::Insert(a, s) => m.insert(a, s.as_bytes()).is_ok(),
                Op::Delete(a, b) => m.delete(Range::new(a, b)).is_ok(),
                Op::Replace(a, b, s) => m.replace(Range::new(a, b), s.as_bytes()).is_ok(),
                Op::Rotate(o) => m.rotate_origin(o).is_ok(),
                Op::RevComp => m.reverse_complement().is_ok(),
                Op::Linear => m.set_topology(Topology::Linear).is_ok(),
                Op::Circular => m.set_topology(Topology::Circular).is_ok(),
                Op::AddFeature(a, b) => m.add_feature(feature("f", a, b)).is_ok(),
            };
            if !ok {
                prop_assert_eq!(&m, &before, "a failed edit changed the document");
            }
            prop_assert!(m.validate().is_ok(), "{:?}", m.validate());
        }
    }
}
