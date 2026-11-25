//! Synthetic sample constructs for demos, screenshots and a first cloning experiment. The
//! sequences are made here from a fixed pseudo-random generator, so they are the same on every
//! platform and every run, and no third-party sequence data is involved.
//!
//! `pSample` is a 3.2 kb circular plasmid with made-up features and a polylinker, and
//! `demoInsert` is a linear fragment with ends that fit it: cutting both with EcoRI and BamHI
//! and ligating gives one directional product, which is what the Cloning experiment's sample
//! setup walks through.

use serde_json::json;

use crate::{CmdError, Session};

/// Length of the sample plasmid.
const LENGTH: usize = 3200;
/// Where the polylinker replaces plasmid residues. The stretch between the terminator (ends at
/// 1640) and the origin (starts at 1900) is free, so no feature moves and the plasmid keeps its
/// length: every coordinate below stays as it was before the polylinker existed.
const MCS_AT: usize = 1700;
/// EcoRI, BamHI, SalI, XhoI, NotI and PstI, each cutting the finished plasmid exactly once.
const MCS: &str = "GAATTCGGATCCGTCGACCTCGAGGCGGCCGCTGCAGT";
/// The demo insert's middle, between its EcoRI and BamHI ends.
const INSERT_CORE: usize = 750;

/// Residues from a xorshift32 generator: deterministic, so screenshots and tests match.
fn residues(seed: u32, n: usize) -> String {
    let mut x = seed;
    (0..n)
        .map(|_| {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            match x % 4 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                _ => 'T',
            }
        })
        .collect()
}

/// The sample plasmid's residues: the pseudo-random backbone with the polylinker written over
/// the free stretch at [`MCS_AT`].
fn plasmid_residues() -> String {
    let mut s = residues(0x5eed_1234, LENGTH);
    let end = MCS_AT + MCS.len();
    // Both bounds are inside the string and on char boundaries (ASCII), but splice by
    // rebuilding rather than indexing, so a changed constant can never slice mid-character.
    let before: String = s.chars().take(MCS_AT).collect();
    let after: String = s.chars().skip(end).collect();
    s = format!("{before}{MCS}{after}");
    s
}

/// The demo insert: EcoRI at the 5' end, BamHI at the 3' end, and a middle that holds neither
/// (nor any other polylinker site), so one cut with each enzyme gives one fragment that fits
/// the opened plasmid in one orientation.
fn insert_residues() -> String {
    format!("GAATTC{}GGATCC", residues(0x779b_1000, INSERT_CORE))
}

/// Open the sample plasmid as the active document. `discard_changes` is passed to `new`, so
/// unsaved work is refused unless the caller has asked to replace it; the apps pass `true`
/// only at start-up, where the session is empty.
pub fn load(session: &mut Session, discard_changes: bool) -> Result<(), CmdError> {
    session.execute(
        "new",
        json!({
            "name": "pSample",
            "residues": plasmid_residues(),
            "topology": "circular",
            "discard_changes": discard_changes,
        }),
    )?;
    let features = [
        ("demo promoter", "promoter", "forward", 120, 330),
        ("demo gene", "CDS", "forward", 360, 1420),
        ("tag", "CDS", "forward", 1420, 1480),
        ("terminator", "terminator", "none", 1520, 1640),
        (
            "polylinker",
            "misc_feature",
            "none",
            MCS_AT,
            MCS_AT + MCS.len(),
        ),
        ("origin", "rep_origin", "reverse", 1900, 2490),
        ("marker", "CDS", "reverse", 2700, 60),
    ];
    for (name, kind, strand, start, end) in features {
        session.execute(
            "add_feature",
            json!({"feature": {
                "name": name, "kind": kind, "strand": strand,
                "segments": [{"range": {"start": start, "end": end}}],
            }}),
        )?;
    }
    Ok(())
}

/// Open the demo insert as a new document.
pub fn load_insert(session: &mut Session, discard_changes: bool) -> Result<(), CmdError> {
    session.execute(
        "new",
        json!({
            "name": "demoInsert",
            "residues": insert_residues(),
            "topology": "linear",
            "discard_changes": discard_changes,
        }),
    )?;
    session.execute(
        "add_feature",
        json!({"feature": {
            "name": "demo insert", "kind": "misc_feature", "strand": "forward",
            "segments": [{"range": {"start": 6, "end": 6 + INSERT_CORE}}],
        }}),
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn count(session: &mut Session, enzyme: &str) -> usize {
        let v = session
            .execute("enzyme_sites", json!({"enzymes": [enzyme]}))
            .expect("sites");
        v["enzymes"][0]["sites"]
            .as_array()
            .map_or(0, |s: &Vec<serde_json::Value>| s.len())
    }

    /// The polylinker gives the plasmid six enzymes that cut it exactly once, and writing it
    /// over a free stretch leaves the length and every feature where they were.
    #[test]
    fn the_sample_plasmid_has_a_usable_polylinker() {
        let mut s = Session::new();
        load(&mut s, true).expect("sample");
        let m = &s.document().expect("open").molecule;
        assert_eq!(m.len(), LENGTH, "the polylinker replaced, not inserted");
        assert_eq!(m.features().len(), 7);
        for enzyme in ["EcoRI", "BamHI", "SalI", "XhoI", "NotI", "PstI"] {
            assert_eq!(count(&mut s, enzyme), 1, "{enzyme} is not a unique cutter");
        }
    }

    /// The insert carries one site of each cloning enzyme, at its ends and nowhere else.
    #[test]
    fn the_demo_insert_fits_the_polylinker() {
        let mut s = Session::new();
        load_insert(&mut s, true).expect("insert");
        let m = &s.document().expect("open").molecule;
        assert_eq!(m.len(), INSERT_CORE + 12);
        for enzyme in ["EcoRI", "BamHI"] {
            assert_eq!(count(&mut s, enzyme), 1, "{enzyme} cuts the insert twice");
        }
        for enzyme in ["SalI", "XhoI", "NotI", "PstI", "HindIII"] {
            assert_eq!(count(&mut s, enzyme), 0, "{enzyme} cuts the insert");
        }
    }

    /// The two together make one product, which is what the sample cloning experiment shows.
    #[test]
    fn the_samples_clone_into_one_product() {
        let mut s = Session::new();
        load(&mut s, true).expect("sample");
        load_insert(&mut s, true).expect("insert");
        let ids: Vec<u64> = s.documents().iter().map(|d| d.id).collect();
        let r = s
            .execute(
                "restriction_clone",
                json!({
                    "vector_document": ids[0],
                    "insert_document": ids[1],
                    "vector_enzymes": ["EcoRI", "BamHI"],
                    "preview": true,
                }),
            )
            .expect("the samples clone");
        let constructs = r["constructs"].as_array().expect("constructs");
        assert_eq!(constructs.len(), 1, "one orientation fits: {r}");
        // The plasmid loses the 6 bp between its two cuts, and the insert brings 6 bp of its
        // own site halves back, so the product is the plasmid plus the insert's core exactly.
        assert_eq!(constructs[0]["length"], LENGTH + INSERT_CORE);
    }
}
