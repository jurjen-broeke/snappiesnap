//! The standard commands, driven through `Session::execute` exactly as agents and the UI do.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::panic)]

mod common;

use serde_json::{Value, json};
use snappiesnap_engine::{ErrorCode, Session};

fn run(s: &mut Session, id: &str, params: Value) -> Value {
    s.execute(id, params)
        .unwrap_or_else(|e| panic!("{id} failed: {e}"))
}

fn code(s: &mut Session, id: &str, params: Value) -> ErrorCode {
    s.execute(id, params).unwrap_err().code
}

fn residues(s: &mut Session, start: usize, end: usize) -> String {
    let v = run(
        s,
        "get_residues",
        json!({"range": {"start": start, "end": end}}),
    );
    v["residues"].as_str().unwrap().to_string()
}

fn plasmid() -> Session {
    let mut s = Session::new();
    run(
        &mut s,
        "new",
        json!({"name": "pDemo", "residues": "AAAACCCC GGGGTTTT\nACGT", "topology": "circular"}),
    );
    s
}

#[test]
fn standard_registry_is_clean() {
    let s = Session::new();
    assert!(
        s.log_tail(100).is_empty(),
        "registry problems: {:?}",
        s.log_tail(100)
    );
    assert!(s.registry().len() >= 15);
    for c in s.registry().iter() {
        assert_eq!(c.params_schema()["type"], "object", "{}", c.id);
        assert!(!c.label.is_empty() && !c.description.is_empty(), "{}", c.id);
    }
}

#[test]
fn a_cloning_style_session() {
    let mut s = plasmid();
    assert_eq!(s.molecule().unwrap().len(), 20);

    let id = run(
        &mut s,
        "add_feature",
        json!({"feature": {
            "name": "ori-span", "kind": "rep_origin", "strand": "forward",
            "segments": [{"range": {"start": 16, "end": 4}}]
        }}),
    )["id"]
        .clone();
    assert_eq!(id, 1);

    run(&mut s, "insert", json!({"at": 0, "residues": "tt"}));
    let f = &s.molecule().unwrap().features()[0];
    assert_eq!(f.segments[0].range.start, 18);
    assert_eq!(f.segments[0].range.end, 6);

    // Linearising is refused while the feature spans the origin, with a hint.
    let err = s
        .execute("set_topology", json!({"topology": "linear"}))
        .unwrap_err();
    assert_eq!(err.code, ErrorCode::InvalidInput);
    assert!(err.hint.unwrap().contains("rotate_origin"));

    run(&mut s, "rotate_origin", json!({"origin": 18}));
    run(&mut s, "set_topology", json!({"topology": "linear"}));
    assert_eq!(residues(&mut s, 0, 4), "ACGT");

    let out = run(&mut s, "delete", json!({"range": {"start": 0, "end": 10}}));
    assert_eq!(out["removed_features"], json!([1]));

    // Every edit is undoable, back to the freshly created molecule.
    for _ in 0..5 {
        run(&mut s, "undo", json!({}));
    }
    assert_eq!(residues(&mut s, 0, 20), "AAAACCCCGGGGTTTTACGT");
    assert!(s.molecule().unwrap().features().is_empty());
    assert_eq!(code(&mut s, "undo", json!({})), ErrorCode::NotEnabled);
}

#[test]
fn unsaved_work_is_protected_and_documents_open_in_tabs() {
    let mut s = plasmid();
    run(&mut s, "insert", json!({"at": 0, "residues": "A"}));
    // A new document opens next to the unsaved one, with its own undo history.
    run(&mut s, "new", json!({"name": "second", "residues": "GGGG"}));
    let docs = run(&mut s, "documents", json!({}));
    let names: Vec<(&str, bool, bool)> = docs["documents"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| {
            (
                d["name"].as_str().unwrap(),
                d["dirty"].as_bool().unwrap(),
                d["active"].as_bool().unwrap(),
            )
        })
        .collect();
    assert_eq!(names, vec![("pDemo", true, false), ("second", false, true)]);
    assert!(!s.can_undo());
    run(&mut s, "switch_document", json!({"index": 0}));
    assert!(s.can_undo());
    assert_eq!(s.molecule().unwrap().len(), 21);
    // Closing unsaved work needs consent; then the other tab becomes active.
    assert_eq!(code(&mut s, "close", json!({})), ErrorCode::NotEnabled);
    run(&mut s, "close", json!({"discard_changes": true}));
    assert_eq!(s.molecule().unwrap().name(), "second");
    run(&mut s, "close", json!({}));
    assert_eq!(
        code(&mut s, "insert", json!({"at": 0, "residues": "A"})),
        ErrorCode::NoDocument
    );
    assert_eq!(
        code(&mut s, "switch_document", json!({"index": 0})),
        ErrorCode::InvalidParams
    );
}

#[test]
fn selection_is_validated_and_not_undoable() {
    let mut s = plasmid();
    let out = run(&mut s, "select", json!({"range": {"start": 18, "end": 2}}));
    assert_eq!(out["selection"], json!({"start": 18, "end": 2}));
    assert_eq!(
        code(&mut s, "select", json!({"range": {"start": 0, "end": 21}})),
        ErrorCode::InvalidInput
    );
    run(&mut s, "select", json!({"range": null}));
    assert!(s.document().unwrap().selection.is_none());
    assert!(!s.can_undo());
}

#[test]
fn features_and_primers() {
    let mut s = plasmid();
    run(
        &mut s,
        "add_feature",
        json!({"feature": {"name": "a", "segments": [{"range": {"start": 0, "end": 4}}]}}),
    );
    let mut f = serde_json::to_value(&s.molecule().unwrap().features()[0]).unwrap();
    assert_eq!(f["kind"], "misc_feature");
    f["name"] = json!("renamed");
    run(&mut s, "update_feature", json!({"feature": f}));
    assert_eq!(s.molecule().unwrap().features()[0].name, "renamed");
    assert_eq!(
        code(&mut s, "remove_feature", json!({"id": 99})),
        ErrorCode::NotFound
    );
    run(&mut s, "remove_feature", json!({"id": 1}));

    run(
        &mut s,
        "add_primer",
        json!({"name": "fwd", "residues": "ACGTACGTAC"}),
    );
    assert_eq!(
        code(
            &mut s,
            "add_primer",
            json!({"name": "bad", "residues": "ACGU"})
        ),
        ErrorCode::InvalidInput
    );
    run(&mut s, "remove_primer", json!({"id": 1}));
    assert!(s.molecule().unwrap().primers().is_empty());
}

#[test]
fn restriction_enzymes() {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    let mut s = common::session();
    let m = Molecule::new(
        "m",
        Seq::new(
            Alphabet::Dna,
            format!(
                "{}GAATTC{}GGATCC{}",
                "A".repeat(10),
                "A".repeat(20),
                "A".repeat(30)
            ),
        )
        .unwrap(),
        Topology::Circular,
    )
    .unwrap();
    s.open_document(snappiesnap_engine::Document::new(m));

    let info = s.execute("enzymes_info", json!({})).unwrap();
    assert_eq!(info["source"], "built-in");
    let sites = s.execute("enzyme_sites", json!({})).unwrap();
    let names: Vec<&str> = sites["enzymes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["name"].as_str().unwrap())
        .collect();
    assert!(
        names.contains(&"EcoRI") && names.contains(&"BamHI"),
        "{names:?}"
    );
    let eco = sites["enzymes"]
        .as_array()
        .unwrap()
        .iter()
        .find(|e| e["name"] == "EcoRI")
        .unwrap();
    assert_eq!(eco["notation"], "G^AATTC");
    assert_eq!(eco["sites"][0]["cut"], 11);

    let d = s
        .execute("digest", json!({"enzymes": ["ecori", "BamHI"]}))
        .unwrap();
    let lens: Vec<u64> = d["fragments"]
        .as_array()
        .unwrap()
        .iter()
        .map(|f| f["length"].as_u64().unwrap())
        .collect();
    assert_eq!(lens, vec![46, 26]);
    let e = s
        .execute("digest", json!({"enzymes": ["NoSuchI"]}))
        .unwrap_err();
    assert!(e.message.contains("NoSuchI"));

    let list = s
        .execute("enzymes_list", json!({"filter": "ggatcc"}))
        .unwrap();
    assert!(
        list["enzymes"]
            .as_array()
            .unwrap()
            .iter()
            .any(|e| e["name"] == "BamHI")
    );

    // Download (served by the in-memory platform), cached for the next session.
    let r = s.execute("enzymes_download", json!({})).unwrap();
    assert_eq!(
        (
            r["source"].as_str(),
            r["count"].as_u64(),
            r["version"].as_str()
        ),
        (Some("REBASE"), Some(2), Some("999"))
    );
    assert_eq!(r["cache"], "app-data/rebase_withrefm.txt");
    assert!(
        s.execute("enzymes_info", json!({})).unwrap()["credit"]
            .as_str()
            .unwrap()
            .contains("REBASE")
    );
    let sites = s
        .execute("enzyme_sites", json!({"commercial_only": false}))
        .unwrap();
    let names: Vec<&str> = sites["enzymes"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, vec!["AaaI"]);

    s.execute("enzymes_reset", json!({})).unwrap();
    assert_eq!(s.enzymes().source, "built-in");
    assert!(s.load_enzyme_cache());
    assert_eq!(s.enzymes().source, "REBASE");

    // Broken files and failed downloads leave the database alone.
    s.execute("enzymes_reset", json!({})).unwrap();
    assert!(
        s.execute("enzymes_load", json!({"path": "example.dna"}))
            .is_err()
    );
    assert!(
        s.execute(
            "enzymes_download",
            json!({"url": "http://nowhere.invalid/x"})
        )
        .is_err()
    );
    assert_eq!(s.enzymes().source, "built-in");
}

#[test]
fn enzyme_sets() {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    let mut s = common::session();
    // The set saved by an earlier session is loaded.
    let list = run(&mut s, "enzyme_sets_list", json!({}));
    assert_eq!(list["sets"][0]["name"], "Lab freezer");

    // An exported set: named after its title. BlnI is not a built-in enzyme: REBASE 610
    // lists no supplier for that name (AvrII and XmaJI are the sold names for its site), and
    // the built-in set is what is sold.
    let r = run(&mut s, "enzyme_set_import", json!({"path": "enzymes.txt"}));
    assert_eq!(r["name"], "LAB ENZYMES");
    assert_eq!(r["count"], 4);
    assert_eq!(r["missing"], json!(["BlnI"]));
    assert!(r["hint"].as_str().unwrap().contains("REBASE"));
    assert_eq!(r["saved"], "app-data/enzyme_sets.json");
    // Importing again replaces it; another name adds a set.
    let r = run(&mut s, "enzyme_set_import", json!({"path": "enzymes.txt"}));
    assert_eq!(r["replaced"], true);
    run(
        &mut s,
        "enzyme_set_import",
        json!({"path": "enzymes.txt", "name": "Copy"}),
    );
    assert_eq!(s.enzymes().sets.len(), 3);

    // Sites from a set: the filters apply within it, commercial_only doesn't.
    let m = Molecule::new(
        "m",
        Seq::new(
            Alphabet::Dna,
            format!("GAATTC{0}GGATCC{0}CCCGGG{0}", "A".repeat(10)),
        )
        .unwrap(),
        Topology::Circular,
    )
    .unwrap();
    s.open_document(snappiesnap_engine::Document::new(m));
    let names = |s: &mut Session, set: &str| -> Vec<String> {
        run(s, "enzyme_sites", json!({"set": set}))["enzymes"]
            .as_array()
            .unwrap()
            .iter()
            .map(|e| e["name"].as_str().unwrap().to_string())
            .collect()
    };
    assert_eq!(names(&mut s, "lab enzymes"), vec!["BamHI", "EcoRI"]);
    assert_eq!(names(&mut s, "Lab freezer"), vec!["EcoRI", "SmaI"]);
    assert_eq!(
        code(&mut s, "enzyme_sites", json!({"set": "nope"})),
        ErrorCode::NotFound
    );

    // Sets survive switching the database, and a new session loads them from the data folder.
    run(&mut s, "enzymes_load", json!({"path": "withrefm.txt"}));
    assert_eq!(s.enzymes().sets.len(), 3);
    run(&mut s, "enzymes_reset", json!({}));
    let mut again = Session::new().with_platform(s.platform_arc());
    assert_eq!(again.load_enzyme_sets(), 3);

    // Bad input is refused and changes nothing.
    for bad in [
        json!({"name": "x", "enzymes": []}),
        json!({"name": " ", "enzymes": ["EcoRI"]}),
        json!({"name": "x", "enzymes": ["EcoRI", "not an enzyme"]}),
    ] {
        assert_eq!(
            code(&mut s, "enzyme_set_create", bad),
            ErrorCode::InvalidParams
        );
    }
    assert_eq!(
        code(&mut s, "enzyme_set_import", json!({"path": "example.dna"})),
        ErrorCode::InvalidInput
    );
    assert_eq!(s.enzymes().sets.len(), 3);

    run(&mut s, "enzyme_set_delete", json!({"name": "copy"}));
    assert_eq!(
        code(&mut s, "enzyme_set_delete", json!({"name": "copy"})),
        ErrorCode::NotFound
    );
    assert_eq!(again.load_enzyme_sets(), 2);
}

#[test]
fn gel_experiment_across_constructs() {
    let mut s = common::session();
    let a = format!("{0}GAATTC{1}GAATTC{0}", "A".repeat(500), "C".repeat(1000));
    run(
        &mut s,
        "new",
        json!({"name": "pA", "residues": a, "topology": "circular"}),
    );
    run(
        &mut s,
        "new",
        json!({"name": "pB", "residues": "ACGT".repeat(300), "topology": "linear"}),
    );
    let docs = run(&mut s, "documents", json!({}));
    let id = |name: &str| -> u64 {
        docs["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["name"] == name)
            .unwrap()["id"]
            .as_u64()
            .unwrap()
    };
    let (a, b) = (id("pA"), id("pB"));
    let gel = run(
        &mut s,
        "gel_experiment",
        json!({"lanes": [
            {"document": a, "method": {"kind": "digest", "enzymes": ["EcoRI"]}},
            {"document": b, "method": {"kind": "uncut"}},
            {"document": a},
            {"document": b, "method": {"kind": "pcr", "forward": "nope", "reverse": "nope"}},
        ]}),
    );
    let lanes = gel["lanes"].as_array().unwrap();
    assert_eq!(lanes.len(), 5);
    assert_eq!(lanes[0]["title"], "1 kb");
    // The digest of pA (not the active document) runs as on pA's own gel.
    assert_eq!(lanes[1]["title"], "pA: EcoRI");
    assert_eq!(lanes[1]["document"], a);
    // Two 1006 bp fragments co-migrate as one band.
    assert_eq!(lanes[1]["bands"].as_array().unwrap().len(), 1);
    assert_eq!(lanes[1]["bands"][0]["bp"], 1006);
    assert_eq!(
        lanes[1]["bands"][0]["fragments"].as_array().unwrap().len(),
        2
    );
    assert_eq!(lanes[2]["bands"][0]["bp"], 1200);
    assert_eq!(lanes[3]["bands"], json!([]));
    assert!(lanes[4]["error"].as_str().unwrap().contains("nope"));

    assert_eq!(
        code(
            &mut s,
            "gel_experiment",
            json!({"lanes": [{"document": 999}]})
        ),
        ErrorCode::NotFound
    );
    assert_eq!(
        code(
            &mut s,
            "gel_experiment",
            json!({"ladder": "nope", "lanes": [{"document": a}]})
        ),
        ErrorCode::InvalidInput
    );
    assert_eq!(
        code(&mut s, "gel_experiment", json!({"lanes": []})),
        ErrorCode::InvalidParams
    );
}

#[test]
fn lab_oligo_library() {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    let mut s = common::session();
    // The library saved by an earlier session.
    assert_eq!(run(&mut s, "oligo_library_list", json!({}))["total"], 2);

    // A tab-separated list with a header: T7 is new, M13F is already there.
    let r = run(
        &mut s,
        "oligo_library_import",
        json!({"path": "oligos.txt"}),
    );
    assert_eq!(
        (r["added"].as_u64(), r["duplicates"].as_u64()),
        (Some(1), Some(1))
    );
    assert_eq!(
        (r["total"].as_u64(), r["skipped_count"].as_u64()),
        (Some(3), Some(0))
    );
    assert_eq!(r["saved"], "app-data/oligo_library.json");
    let found = run(&mut s, "oligo_library_list", json!({"filter": "promoter"}));
    assert_eq!(found["oligos"][0]["name"], "T7");

    // A construct with the tail oligo's 3' part near its start and M13F's site after it
    // (reverse-complemented, so M13F binds the top strand as a reverse primer).
    let m13f_rc = "ACTGGCCGTCGTTTTAC";
    let residues = format!(
        "{0}ACGTTTTTACGTCAGT{0}{m13f_rc}{0}TTACGTCAGT{0}",
        "C".repeat(30)
    );
    let m = Molecule::new(
        "pTest",
        Seq::new(Alphabet::Dna, residues).unwrap(),
        Topology::Circular,
    )
    .unwrap();
    s.open_document(snappiesnap_engine::Document::new(m));
    let found = run(&mut s, "oligo_library_find", json!({}));
    assert_eq!(found["library"], 3);
    let names: Vec<&str> = found["oligos"]
        .as_array()
        .unwrap()
        .iter()
        .map(|o| o["name"].as_str().unwrap())
        .collect();
    assert_eq!(names, vec!["tail", "M13F"], "by position; T7 doesn't bind");
    assert_eq!(found["oligos"][1]["sites"][0]["strand"], "reverse");
    assert_eq!(found["oligos"][0]["in_molecule"], false);
    // Only primary sites by default: "tail"'s 3' end also pairs with a 10-base stretch.
    assert_eq!(found["oligos"][0]["sites"].as_array().unwrap().len(), 1);
    // That 10-base site melts below the default minimum Tm (40 °C), so it only counts once
    // the minimum is off.
    let weak = run(&mut s, "oligo_library_find", json!({"primary_only": false}));
    assert_eq!(weak["oligos"][0]["sites"].as_array().unwrap().len(), 1);
    // The construct is short, so the short sequences' rules apply.
    let rules = run(
        &mut s,
        "primer_binding",
        json!({"short": {"min_tm_on": false}}),
    );
    assert_eq!(
        (
            rules["short"]["min_tm_on"].as_bool(),
            rules["short"]["min_tm"].as_f64()
        ),
        (Some(false), Some(40.0))
    );
    let weak = run(&mut s, "oligo_library_find", json!({"primary_only": false}));
    assert_eq!(weak["oligos"][0]["sites"].as_array().unwrap().len(), 2);
    // Explicit rules still work as before, whatever the session's.
    let exact = run(
        &mut s,
        "oligo_library_find",
        json!({"primary_only": false, "window": 15, "mismatches": 0, "clamp": 15}),
    );
    assert_eq!(exact["oligos"][0]["sites"].as_array().unwrap().len(), 1);
    run(
        &mut s,
        "primer_binding",
        json!({"short": {"min_tm_on": true}}),
    );

    // Added as primers in one undoable step, with their case kept.
    let r = run(
        &mut s,
        "oligo_library_add",
        json!({"names": ["m13f", "tail"]}),
    );
    assert_eq!(r["added"].as_array().unwrap().len(), 2);
    let primers = s.document().unwrap().molecule.primers().to_vec();
    assert_eq!(primers.len(), 2);
    assert_eq!(primers[1].seq.as_bytes(), b"gaattcACGTTTTTACGTCAGT");
    let found = run(&mut s, "oligo_library_find", json!({}));
    assert_eq!(found["oligos"][0]["in_molecule"], true);
    assert_eq!(
        code(&mut s, "oligo_library_add", json!({"names": ["M13F"]})),
        ErrorCode::NotEnabled
    );
    assert_eq!(
        code(&mut s, "oligo_library_add", json!({"names": ["nope"]})),
        ErrorCode::NotFound
    );
    run(&mut s, "undo", json!({}));
    assert!(s.document().unwrap().molecule.primers().is_empty());

    // Removing, then replacing the whole library.
    let r = run(&mut s, "oligo_library_remove", json!({"names": ["TAIL"]}));
    assert_eq!(
        (r["removed"].as_u64(), r["total"].as_u64()),
        (Some(1), Some(2))
    );
    let r = run(
        &mut s,
        "oligo_library_import",
        json!({"path": "oligos.txt", "replace": true}),
    );
    assert_eq!(
        (r["added"].as_u64(), r["total"].as_u64()),
        (Some(2), Some(2))
    );
    run(&mut s, "oligo_library_remove", json!({"all": true}));
    assert_eq!(run(&mut s, "oligo_library_list", json!({}))["total"], 0);
    assert_eq!(
        code(&mut s, "oligo_library_remove", json!({})),
        ErrorCode::InvalidParams
    );
    // A file that isn't an oligo list changes nothing.
    assert_eq!(
        code(
            &mut s,
            "oligo_library_import",
            json!({"path": "example.dna"})
        ),
        ErrorCode::InvalidInput
    );
}

/// A lab keeps its cloning primers one per file, the name only in the file name.
#[test]
fn oligo_folder_import() {
    let mut s = common::session();
    // Each file takes its name from its own name, and the .dna beside them is left alone.
    // Case is kept, so a 5' tail stays readable.
    let r = run(&mut s, "oligo_library_import", json!({"path": "primers"}));
    assert_eq!(
        (r["added"].as_u64(), r["files"].as_u64()),
        (Some(2), Some(2))
    );
    assert_eq!(r["skipped_count"].as_u64(), Some(0));
    assert!(r["unreadable"].is_null(), "{r}");
    let listed = run(&mut s, "oligo_library_list", json!({"filter": "SP_20"}));
    assert_eq!(listed["oligos"].as_array().map(Vec::len), Some(2));
    assert_eq!(listed["oligos"][0]["name"], "SP_201");
    assert_eq!(
        listed["oligos"][0]["residues"],
        "ggatccagatctacgtacgtacgttgcatgcatgcaag"
    );
    assert_eq!(listed["oligos"][1]["length"], 43);

    // Importing the same folder again changes nothing.
    let again = run(&mut s, "oligo_library_import", json!({"path": "primers"}));
    assert_eq!(
        (again["added"].as_u64(), again["duplicates"].as_u64()),
        (Some(0), Some(2))
    );

    // One such file on its own is read the same way.
    let one = run(
        &mut s,
        "oligo_library_import",
        json!({"path": "primers/SP_201.txt", "replace": true}),
    );
    assert_eq!(one["added"].as_u64(), Some(1));
    assert!(one["files"].is_null(), "a single file needs no file count");

    // A folder with no oligo files, and a file that holds no oligos, both say so.
    assert_eq!(
        code(&mut s, "oligo_library_import", json!({"path": "lab"})),
        ErrorCode::InvalidInput
    );
    assert_eq!(
        code(
            &mut s,
            "oligo_library_import",
            json!({"path": "primers/notes.dna"})
        ),
        ErrorCode::InvalidInput
    );
}

#[test]
fn primer_binding_rules() {
    let mut s = Session::new();
    let rules = run(&mut s, "primer_binding", json!({}));
    let class = |matched: u64, isolated: bool| json!({"matched": matched, "isolated_mismatches": isolated, "min_tm_on": true, "min_tm": 40.0});
    assert_eq!(
        rules,
        json!({"long_from_bp": 1_000_000, "short": class(10, true), "long": class(15, false),
               "five_prime_on": true, "five_prime_matched": 15})
    );
    let rules = run(
        &mut s,
        "primer_binding",
        json!({"five_prime_on": false, "five_prime_matched": 20}),
    );
    assert_eq!(
        (
            rules["five_prime_on"].clone(),
            rules["five_prime_matched"].clone()
        ),
        (json!(false), json!(20))
    );
    let before = s.revision();
    let rules = run(
        &mut s,
        "primer_binding",
        json!({"long_from_bp": 50_000, "long": {"matched": 20, "min_tm_on": false}}),
    );
    assert_eq!(rules["long_from_bp"], 50_000);
    assert_eq!(
        rules["long"],
        json!({"matched": 20, "isolated_mismatches": false, "min_tm_on": false, "min_tm": 40.0})
    );
    assert_eq!(
        rules["short"],
        class(10, true),
        "the other class is left alone"
    );
    assert_ne!(s.revision(), before, "views showing sites update");
    for bad in [
        json!({"short": {"matched": 9}}),
        json!({"long": {"matched": 61}}),
        json!({"long_from_bp": 5}),
        json!({"long_from": 25}),
        json!({"short": {"min_tm": -1}}),
        json!({"long": {"min_tm": 120}}),
        json!({"short_match": 12}),
        json!({"short": {"isolated": true}}),
        json!({"five_prime_matched": 9}),
    ] {
        assert_eq!(
            code(&mut s, "primer_binding", bad.clone()),
            ErrorCode::InvalidParams,
            "{bad}"
        );
    }
    assert_eq!(
        run(&mut s, "primer_binding", json!({}))["long"]["matched"],
        20
    );
}

#[test]
fn primer_binding_sites() {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    let mut s = common::session();
    let m = Molecule::new(
        "m",
        Seq::new(
            Alphabet::Dna,
            "AAAAAAAAAAACGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAAAAAAAAAA",
        )
        .unwrap(),
        Topology::Circular,
    )
    .unwrap();
    s.open_document(snappiesnap_engine::Document::new(m));
    s.execute(
        "add_primer",
        json!({"name": "f", "residues": "gaattcCGTTGCATGCAAGCTT"}),
    )
    .unwrap();
    s.execute(
        "add_primer",
        json!({"name": "r", "residues": "GTAAAACGACGGCCAGTGCCA"}),
    )
    .unwrap();
    let r = s.execute("primer_sites", json!({})).unwrap();
    let f = &r["primers"][0]["sites"][0];
    assert_eq!(
        (
            f["strand"].as_str(),
            f["annealed"].as_u64(),
            f["primary"].as_bool()
        ),
        (Some("forward"), Some(16), Some(true))
    );
    assert_eq!(f["range"], json!({"start": 11, "end": 27}));
    let rv = &r["primers"][1]["sites"][0];
    assert_eq!(
        (rv["strand"].as_str(), rv["range"]["start"].as_u64()),
        (Some("reverse"), Some(26))
    );
}

#[test]
fn melting_temperatures() {
    let mut s = common::session();
    let r = s
        .execute("tm", json!({"residues": "CGTTGCATGCAAGCTT"}))
        .unwrap();
    assert_eq!(
        (r["tm"].as_f64(), r["gc_percent"].as_f64()),
        (Some(52.2), Some(50.0))
    );
    let r = s
        .execute(
            "tm",
            json!({"residues": "CGTTGCATGCAAGCTT", "conditions": {"mg_mm": 1.5, "dntp_mm": 0.2}}),
        )
        .unwrap();
    assert_eq!(r["tm"].as_f64(), Some(56.2));
    assert!(s.execute("tm", json!({"residues": "ACGN"})).is_err());
    assert!(
        s.execute(
            "tm",
            json!({"residues": "ACGTACGT", "conditions": {"na_mm": -1}})
        )
        .is_err()
    );
}

#[test]
fn primer_design() {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    let mut x: u32 = 777;
    let s: String = (0..400)
        .map(|_| {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12345);
            b"ACGT"[((x >> 16) % 4) as usize] as char
        })
        .collect();
    let mut session = common::session();
    let m = Molecule::new(
        "m",
        Seq::new(Alphabet::Dna, s.as_str()).unwrap(),
        Topology::Linear,
    )
    .unwrap();
    session.open_document(snappiesnap_engine::Document::new(m));
    let r = session
        .execute(
            "design_primers",
            json!({"range": {"start": 100, "end": 300}, "count": 3}),
        )
        .unwrap();
    let pairs = r["pairs"].as_array().unwrap();
    assert_eq!(pairs.len(), 3);
    let best = &pairs[0];
    assert_eq!(best["product"], 200);
    let fwd = best["forward"]["residues"].as_str().unwrap();
    assert!(s[100..].starts_with(fwd));
    // Keeping it: the designed primers bind where they were designed.
    session
        .execute("add_primer", json!({"name": "F", "residues": fwd}))
        .unwrap();
    let sites = session.execute("primer_sites", json!({})).unwrap();
    let starts: Vec<u64> = sites["primers"][0]["sites"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|s| s["primary"] == true)
        .map(|s| s["range"]["start"].as_u64().unwrap())
        .collect();
    assert_eq!(starts, vec![100]);
    assert!(
        session
            .execute("design_primers", json!({"range": {"start": 0, "end": 10}}))
            .is_err()
    );
}

#[test]
fn primer_structure_checks() {
    let mut s = common::session();
    let ok = s
        .execute("primer_checks", json!({"residues": "GTAAAACGACGGCCAGTG"}))
        .unwrap();
    assert_eq!(ok["warnings"], json!([]));
    let bad = s
        .execute(
            "primer_checks",
            json!({"residues": "ACGTACGTGCGGCCGC", "partner": "TTTTTTTTGCGGCCGC"}),
        )
        .unwrap();
    assert_eq!(bad["self_dimer_3prime"]["three_prime"], true);
    assert!(bad["cross_dimer"]["dg"].as_f64().unwrap() < -5.0);
    assert!(bad["warnings"].as_array().unwrap().len() >= 2, "{bad}");
    let hp = s
        .execute("primer_checks", json!({"residues": "GGGGAAAACCCC"}))
        .unwrap();
    assert_eq!(
        (
            hp["hairpin"]["stem"].as_u64(),
            hp["hairpin"]["loop"].as_u64()
        ),
        (Some(4), Some(4))
    );
    assert!(
        s.execute("primer_checks", json!({"residues": "ACGN"}))
            .is_err()
    );
}

#[test]
fn pcr_and_gel() {
    use snappiesnap_engine::doc::{Alphabet, Molecule, Seq, Topology};
    let mut s = common::session();
    let template = "GATTACAGATTACACCGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCG";
    let mut m = Molecule::new(
        "t",
        Seq::new(Alphabet::Dna, template).unwrap(),
        Topology::Circular,
    )
    .unwrap();
    m.add_primer(
        "f",
        Seq::new(Alphabet::Dna, "gaattcCACCGTTGCATGCAAGCTTG").unwrap(),
        "",
    )
    .unwrap();
    m.add_primer(
        "r",
        Seq::new(Alphabet::Dna, "CAGGGTTTTCCCAGTCACGA").unwrap(),
        "",
    )
    .unwrap();
    s.open_document(snappiesnap_engine::Document::new(m));
    let r = s
        .execute("pcr", json!({"forward": "f", "reverse": "r"}))
        .unwrap();
    let products = r["products"].as_array().unwrap();
    assert_eq!(products.len(), 1);
    assert_eq!(products[0]["length"], 6 + 63);
    assert!(
        products[0]["residues"]
            .as_str()
            .unwrap()
            .starts_with("GAATTCCACCG")
    );
    let gel = s
        .execute("gel", json!({"lanes": [{"kind": "ladder", "name": "100 bp"}, {"kind": "pcr", "forward": "f", "reverse": "r"}, {"kind": "uncut"}]}))
        .unwrap();
    assert_eq!(gel["lanes"].as_array().unwrap().len(), 3);
    assert_eq!(gel["lanes"][1]["bands"][0]["bp"], 69);
    assert!(
        s.execute(
            "gel",
            json!({"lanes": [{"kind": "ladder", "name": "nope"}]})
        )
        .is_err()
    );
    let opened = s
        .execute("pcr_open", json!({"forward": "f", "reverse": "r"}))
        .unwrap();
    assert_eq!(opened["length"], 69);
    let doc = s.document().unwrap();
    assert!(doc.dirty);
    let kinds: Vec<&str> = doc
        .molecule
        .features()
        .iter()
        .map(|f| f.kind.as_str())
        .collect();
    assert_eq!(kinds, vec!["primer_bind", "primer_bind"]);
    // Each product opens in its own tab; the template stays open next to it.
    s.execute("switch_document", json!({"index": 0})).unwrap();
    s.execute("pcr_open", json!({"forward": "f", "reverse": "r"}))
        .unwrap();
    assert_eq!(s.documents().len(), 3);
}

#[test]
fn translation_orfs_and_proteins() {
    use snappiesnap_engine::doc::{
        Alphabet, Molecule, NewFeature, Range, Segment, Seq, Strand, Topology,
    };
    let mut s = common::session();
    let r = s
        .execute("translate", json!({"residues": "ATGGCCTGA"}))
        .unwrap();
    assert_eq!(r["protein"], "MA*");
    let r = s
        .execute("translate", json!({"residues": "ATGGCCTGA", "table": 2}))
        .unwrap();
    assert_eq!(r["protein"], "MAW");
    assert!(
        s.execute("translate", json!({"residues": "ATG", "table": 7}))
            .is_err()
    );
    let orf = "ATGAAACCCGGGTTTTAA";
    let rc: String = orf
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect();
    let seq = format!("CC{orf}CCCC{rc}CC");
    let mut m = Molecule::new(
        "m",
        Seq::new(Alphabet::Dna, seq.as_str()).unwrap(),
        Topology::Linear,
    )
    .unwrap();
    let id = m
        .add_feature(NewFeature {
            color: None,
            name: "rev cds".into(),
            kind: "CDS".into(),
            strand: Strand::Reverse,
            segments: vec![Segment {
                range: Range::new(24, 42),
            }],
            qualifiers: vec![],
        })
        .unwrap();
    s.open_document(snappiesnap_engine::Document::new(m));
    let o = s.execute("orfs", json!({"min_length": 5})).unwrap();
    assert_eq!(o["total"], 2);
    let r = s
        .execute(
            "translate",
            json!({"range": {"start": 24, "end": 42}, "strand": "reverse"}),
        )
        .unwrap();
    assert_eq!(r["protein"], "MKPGF*");
    let p = s
        .execute("protein_properties", json!({"feature": id}))
        .unwrap();
    assert_eq!(
        (p["protein"].as_str(), p["length"].as_u64()),
        (Some("MKPGF"), Some(5))
    );
    assert!(
        s.execute("protein_properties", json!({"feature": 999}))
            .is_err()
    );
    let codes = s.execute("genetic_codes", json!({})).unwrap();
    assert_eq!(codes["codes"].as_array().unwrap().len(), 27);
}

#[test]
fn cloning_commands() {
    let mut s = common::session();
    // Gateway between two files: preview, then open.
    let p = s
        .execute(
            "gateway_clone",
            json!({"acceptor": "dest.gb", "donor": "entry.gb", "preview": true}),
        )
        .unwrap();
    let expected =
        "GGGGGGGGGGACAAGTTTGTACAAAAAAGCAGGCTATGAAAGGGCCCATGAAAACCCAGCTTTCTTGTACAAAGTGGTCCCCCCCCCC";
    assert_eq!(p["length"], expected.len());
    assert_eq!(p["sites"].as_array().unwrap().len(), 2);
    s.execute(
        "gateway_clone",
        json!({"acceptor": "dest.gb", "donor": "entry.gb"}),
    )
    .unwrap();
    let doc = s.document().unwrap();
    assert_eq!(
        String::from_utf8_lossy(doc.molecule.seq().as_bytes()),
        expected
    );
    assert!(doc.dirty);
    // Restriction cloning into the open document (the Gateway product has no HindIII site,
    // so use a template with one).
    s.execute("new", json!({"name": "vec", "residues": "GGGGGAAGCTTGGGGGGGGGGGGGGGGGGG", "topology": "circular", "discard_changes": true}))
        .unwrap();
    let r = s
        .execute(
            "restriction_clone",
            json!({"insert": "insert.fa", "vector_enzymes": ["HindIII"], "preview": true}),
        )
        .unwrap();
    let constructs = r["constructs"].as_array().unwrap();
    assert_eq!(constructs.len(), 2, "{r}");
    assert_eq!(constructs[0]["length"], 30 + 19);
    assert!(
        s.execute(
            "restriction_clone",
            json!({"insert": "insert.fa", "vector_enzymes": ["EcoRI"], "preview": true})
        )
        .is_err()
    );
    assert!(
        s.execute(
            "restriction_clone",
            json!({"insert": "missing.fa", "vector_enzymes": ["HindIII"]})
        )
        .is_err()
    );
}

/// A construct that is open in a tab can be cloned without having a file at all: this is how
/// the web build works, where the browser gives the page no file paths, and how an edited tab
/// is cloned as it now stands.
#[test]
fn cloning_from_open_documents() {
    let mut s = common::session();
    // Two tabs, neither saved anywhere.
    s.execute(
        "new",
        json!({"name": "vec", "residues": "GGGGGAAGCTTGGGGGGGGGGGGGGGGGGG", "topology": "circular"}),
    )
    .unwrap();
    s.execute(
        "new",
        json!({"name": "ins", "residues": "AAGCTTCCGGAATTCGGCCAAGCTT", "topology": "linear"}),
    )
    .unwrap();
    let docs = s.execute("documents", json!({})).unwrap();
    let ids: Vec<u64> = docs["documents"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|d| d["id"].as_u64())
        .collect();
    let (vector, insert) = (ids[ids.len() - 2], ids[ids.len() - 1]);
    let r = s
        .execute(
            "restriction_clone",
            json!({
                "vector_document": vector,
                "insert_document": insert,
                "vector_enzymes": ["HindIII"],
                "preview": true,
            }),
        )
        .unwrap();
    assert_eq!(r["constructs"].as_array().unwrap().len(), 2, "{r}");
    assert_eq!(r["constructs"][0]["length"], 30 + 19);
    // The tab is read as it now stands: shortening the insert between its two HindIII sites
    // shortens the product by as much.
    let index = docs["documents"].as_array().unwrap().len() - 1;
    s.execute("switch_document", json!({"index": index}))
        .unwrap();
    s.execute("delete", json!({"range": {"start": 6, "end": 12}}))
        .unwrap();
    let r = s
        .execute(
            "restriction_clone",
            json!({
                "vector_document": vector,
                "insert_document": insert,
                "vector_enzymes": ["HindIII"],
                "preview": true,
            }),
        )
        .unwrap();
    assert_eq!(r["constructs"][0]["length"], 30 + 13, "{r}");
    // Every other reaction takes a document the same way.
    let a = s
        .execute(
            "assemble",
            json!({"method": "gibson", "parts": [{"file": "gibson_a.fa"}, {"document": vector}], "preview": true}),
        )
        .unwrap_err();
    assert_eq!(a.code, snappiesnap_engine::ErrorCode::InvalidInput, "{a:?}");
    let g = s
        .execute(
            "gateway_clone",
            json!({"acceptor_document": vector, "donor": "entry.gb", "preview": true}),
        )
        .unwrap_err();
    assert!(!g.message.is_empty());
    // A tab that isn't open, and naming both a file and a tab, are reported, not guessed at.
    let gone = s
        .execute(
            "restriction_clone",
            json!({"vector_document": 9999, "insert_document": insert, "vector_enzymes": ["HindIII"]}),
        )
        .unwrap_err();
    assert_eq!(gone.code, snappiesnap_engine::ErrorCode::NotFound);
    assert!(gone.message.contains("9999"), "{}", gone.message);
    let both = s
        .execute(
            "restriction_clone",
            json!({"vector": "dest.gb", "vector_document": vector, "vector_enzymes": ["HindIII"]}),
        )
        .unwrap_err();
    assert_eq!(both.code, snappiesnap_engine::ErrorCode::InvalidParams);
}

#[test]
fn assembly_commands() {
    let mut s = common::session();
    let g = s
        .execute("assemble", json!({"method": "gibson", "parts": [{"file": "gibson_a.fa"}, {"file": "gibson_b.fa"}], "preview": true}))
        .unwrap();
    assert_eq!(g["length"], common::GIBSON_A.len() + common::GIBSON_B.len());
    assert_eq!(g["junctions"].as_array().unwrap().len(), 2);
    let gg = s
        .execute("assemble", json!({"method": "golden_gate", "enzyme": "BsaI", "parts": [{"file": "gg_vector.gb"}, {"file": "gg_part.fa"}]}))
        .unwrap();
    let seq = String::from_utf8_lossy(s.document().unwrap().molecule.seq().as_bytes()).to_string();
    assert!(!seq.contains("GGTCTC") && !seq.contains("GAGACC"), "{gg}");
    assert!(seq.contains("ATGCCCGGGAAATTT"));
    // Mutagenesis on a template with back-to-back primers.
    let t = "ACGTTGCAACGTAACCGGTTGGCATGCATTTCGAGCTAGCATCGGATCGATCGTAGCTAGGCATCGATCGACGACTAGCAGCTAG";
    s.execute(
        "new",
        json!({"name": "t", "residues": t, "topology": "circular", "discard_changes": true}),
    )
    .unwrap();
    let fwd = format!("A{}", &t[41..60]);
    let rev: String = t[20..40]
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect();
    let r = s
        .execute(
            "mutagenesis",
            json!({"forward": fwd, "reverse": rev, "preview": true}),
        )
        .unwrap();
    // One base (40) deleted, one inserted (the A tail): same length, a substitution or so.
    assert_eq!(r["length_change"], 0);
    assert!(
        s.execute("mutagenesis", json!({"forward": "ACGT", "reverse": "ACGT"}))
            .is_err()
    );
}

#[test]
fn construction_history() {
    let mut s = common::session();
    s.execute("open", json!({"path": "history.dna"})).unwrap();
    let t = s.execute("history_tree", json!({})).unwrap();
    let root = &t["steps"];
    assert_eq!(
        (
            root["id"].as_u64(),
            root["operation"].as_str(),
            root["length"].as_u64()
        ),
        (Some(3), Some("Insertion"), Some(18))
    );
    let start = &root["inputs"][0];
    assert_eq!(
        (start["operation"].as_str(), start["available"].as_bool()),
        (Some(""), Some(true))
    );
    let opened = s.execute("history_open", json!({"id": 1})).unwrap();
    assert_eq!(
        (
            opened["name"].as_str(),
            opened["length"].as_u64(),
            opened["features"].as_u64()
        ),
        (Some("start"), Some(16), Some(1))
    );
    assert_eq!(
        String::from_utf8_lossy(s.molecule().unwrap().seq().as_bytes()),
        "ACGTACGTAAAACCCC"
    );
    assert_eq!(s.documents().len(), 2);
    // No history: null, not an error; unknown steps are errors.
    s.execute("new", json!({"residues": "ACGT"})).unwrap();
    assert!(s.execute("history_tree", json!({})).unwrap()["steps"].is_null());
    s.execute("switch_document", json!({"index": 0})).unwrap();
    assert!(s.execute("history_open", json!({"id": 99})).is_err());
    assert!(s.execute("history_open", json!({"id": 3})).is_err());
    // doc.get summarises the history instead of repeating every ancestor.
    let d = snappiesnap_engine::api::call(&mut s, "doc.get", json!({})).unwrap();
    assert_eq!(d["construction"]["id"], 3);
    assert!(d["molecule"].get("history").is_none());

    // The parts of a construct open again with their own history.
    // Cloning records its inputs.
    s.execute("open", json!({"path": "dest.gb"})).unwrap();
    s.execute("gateway_clone", json!({"donor": "entry.gb"}))
        .unwrap();
    let t = s.execute("history_tree", json!({})).unwrap();
    let root = &t["steps"];
    assert_eq!(root["operation"], "Gateway LR reaction");
    assert_eq!(root["inputs"].as_array().map(Vec::len), Some(2));
    assert_eq!(root["inputs"][0]["available"], true);

    // Products record how they were made, on top of their inputs' history.
    let clone = s.molecule().unwrap().seq().as_bytes().to_vec();
    let forward = String::from_utf8(clone[0..20].to_vec()).unwrap();
    let reverse: String = clone[60..80]
        .iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    s.execute("pcr_open", json!({"forward": forward, "reverse": reverse}))
        .unwrap();
    let t = s.execute("history_tree", json!({})).unwrap();
    let root = &t["steps"];
    assert_eq!(root["operation"], "PCR");
    assert_eq!(root["details"], "primers 20 nt and 20 nt");
    let parent = &root["inputs"][0];
    assert_eq!(parent["operation"], "Gateway LR reaction");
    let ids: Vec<u64> = [
        &root["id"],
        &parent["id"],
        &parent["inputs"][0]["id"],
        &parent["inputs"][1]["id"],
    ]
    .iter()
    .filter_map(|v| v.as_u64())
    .collect();
    assert_eq!(ids, vec![0, 1, 2, 3]);
    // An ancestor opens with the history that led to it.
    s.execute("history_open", json!({"id": 1})).unwrap();
    let t = s.execute("history_tree", json!({})).unwrap();
    assert_eq!(t["steps"]["operation"], "Gateway LR reaction");
    assert_eq!(t["steps"]["inputs"].as_array().map(Vec::len), Some(2));
    // GenBank can't hold the history and says so; .dna keeps it.
    let saved = s.execute("save", json!({ "path": "lr.gb" })).unwrap();
    assert!(
        saved["warnings"]
            .as_array()
            .unwrap()
            .iter()
            .any(|w| w == "GenBank files don't hold the construction history")
    );
    let before = s.execute("history_tree", json!({})).unwrap()["steps"].clone();
    let saved = s.execute("save", json!({ "path": "lr.dna" })).unwrap();
    assert_eq!(saved["warnings"], json!([]));
    // Read back from the file (a copy: lr.dna is the open tab).
    let bytes = s.platform().read("lr.dna").unwrap();
    s.platform().write("lr_copy.dna", &bytes).unwrap();
    s.execute("open", json!({ "path": "lr_copy.dna" })).unwrap();
    let after = s.execute("history_tree", json!({})).unwrap()["steps"].clone();
    // The same steps, operations, lengths and stored molecules (ids are renumbered; the
    // details are SnappieSnap's own notes and aren't part of the file).
    fn shape(v: &Value) -> Value {
        json!({
            "name": v["name"], "operation": v["operation"], "length": v["length"],
            "circular": v["circular"], "available": v["available"],
            "inputs": v["inputs"].as_array().map(|i| i.iter().map(shape).collect::<Vec<_>>()),
        })
    }
    assert_eq!(shape(&after), shape(&before));
    // Numbered as the reference application does: inputs first, the construct last.
    assert_eq!(after["id"], 2);
    // An ancestor from the saved file opens with its features.
    s.execute("history_open", json!({ "id": 0 })).unwrap();
    let opened = s.molecule().unwrap();
    assert_eq!(opened.name(), before["inputs"][0]["name"].as_str().unwrap());
    // Saving the reopened file unchanged gives the same bytes.
    let index = s.documents().len() - 2;
    s.execute("switch_document", json!({ "index": index }))
        .unwrap();
    s.execute("save", json!({ "path": "lr2.dna" })).unwrap();
    assert_eq!(s.platform().read("lr2.dna").unwrap(), bytes);
}

#[test]
fn ta_and_topo_cloning() {
    let mut s = common::session();
    s.execute("open", json!({"path": "dest.gb"})).unwrap();
    let insert = "AAGCTTCCGGAATTCGGCCAAGCTT";
    let rc: String = insert
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let p = s
        .execute(
            "topo_clone",
            json!({"insert": "insert.fa", "position": 10, "preview": true}),
        )
        .unwrap();
    let list = p["constructs"].as_array().unwrap();
    assert_eq!(list.len(), 2);
    assert_eq!(list[0]["length"], common::DEST.len() + insert.len() + 2);
    assert_eq!(list[1]["reversed"], true);

    // Open the reversed one: vector numbering kept, T + insert + A at the opening.
    s.execute(
        "topo_clone",
        json!({"insert": "insert.fa", "position": 10, "product": 1}),
    )
    .unwrap();
    let got = String::from_utf8_lossy(s.molecule().unwrap().seq().as_bytes()).to_uppercase();
    let (left, right) = common::DEST.split_at(10);
    assert_eq!(got, format!("{left}T{rc}A{right}"));
    let t = s.execute("history_tree", json!({})).unwrap();
    assert_eq!(t["steps"]["operation"], "TA cloning");
    assert_eq!(t["steps"]["inputs"].as_array().map(Vec::len), Some(2));

    // Blunt TOPO of a PCR product of the open document's primers.
    s.execute("switch_document", json!({"index": 0})).unwrap();
    let blunt = s
        .execute(
            "topo_clone",
            json!({"vector": "dest.gb", "pcr": {"forward": "AAGCTTCCGG", "reverse": "AAGCTTGGCC"}, "position": 0, "kind": "blunt", "preview": true}),
        );
    // dest.gb has no such primers' sites: a clear error, not a crash.
    assert!(blunt.is_err());
    // Directional TOPO needs CACC at the start of the insert.
    let e = s
        .execute(
            "topo_clone",
            json!({"insert": "insert.fa", "position": 10, "kind": "directional"}),
        )
        .unwrap_err();
    assert!(e.message.contains("CACC"), "{}", e.message);
    // A circular vector needs a position.
    assert!(
        s.execute("topo_clone", json!({"insert": "insert.fa"}))
            .is_err()
    );
}

#[test]
fn copy_and_paste() {
    let mut s = common::session();
    s.execute(
        "new",
        json!({"name": "p", "residues": "AAAACCCCGG", "topology": "circular"}),
    )
    .unwrap();
    // Copy across the origin, in order; FASTA and reverse complement.
    s.execute("select", json!({"range": {"start": 8, "end": 2}}))
        .unwrap();
    let c = s.execute("copy_selection", json!({})).unwrap();
    assert_eq!(
        (c["text"].as_str(), c["length"].as_u64()),
        (Some("GGAA"), Some(4))
    );
    let f = s
        .execute("copy_selection", json!({"format": "fasta"}))
        .unwrap();
    assert_eq!(f["text"], ">p 9..2\nGGAA\n");
    let r = s
        .execute("copy_selection", json!({"format": "reverse_complement"}))
        .unwrap();
    assert_eq!(r["text"], "TTCC");

    // Paste over the selection across the origin: the residues go at the end, the cursor after.
    let p = s.execute("paste", json!({"text": "ttt"})).unwrap();
    assert_eq!(
        (p["pasted"].as_u64(), p["length"].as_u64()),
        (Some(3), Some(9))
    );
    assert_eq!(
        String::from_utf8_lossy(s.molecule().unwrap().seq().as_bytes()),
        "AACCCCttt"
    );
    assert_eq!(
        s.document().unwrap().selection,
        Some(snappiesnap_engine::doc::Range::at(9))
    );

    // Paste FASTA at the cursor, with numbers and gaps dropped and noted.
    s.execute("select", json!({"range": {"start": 2, "end": 2}}))
        .unwrap();
    let p = s
        .execute("paste", json!({"text": ">x\n1 GG-A\n>y\nCC"}))
        .unwrap();
    assert_eq!(p["notes"].as_array().map(Vec::len), Some(2));
    assert_eq!(
        String::from_utf8_lossy(s.molecule().unwrap().seq().as_bytes()),
        "AAGGACCCCttt"
    );
    assert_eq!(p["end"], 5);

    // Not residues: refused, nothing changed; undo takes back one paste at a time.
    assert!(s.execute("paste", json!({"text": "AC!T"})).is_err());
    assert_eq!(s.molecule().unwrap().len(), 12);
    s.execute("undo", json!({})).unwrap();
    assert_eq!(
        String::from_utf8_lossy(s.molecule().unwrap().seq().as_bytes()),
        "AACCCCttt"
    );
    // Nothing selected to copy.
    s.execute("select", json!({"range": {"start": 3, "end": 3}}))
        .unwrap();
    assert!(s.execute("copy_selection", json!({})).is_err());
}

#[test]
fn find_sequences_and_names() {
    let mut s = common::session();
    s.execute(
        "new",
        json!({"name": "p", "residues": "GAATTCAAATGGTCCCACCAGGGGGGGGGG", "topology": "circular"}),
    )
    .unwrap();
    s.execute(
        "add_feature",
        json!({"feature": {"name": "lacZ alpha", "kind": "CDS", "strand": "forward", "segments": [{"range": {"start": 6, "end": 14}}]}}),
    )
    .unwrap();
    s.execute(
        "add_primer",
        json!({"name": "lac fwd", "residues": "CAAATGGTCCCACCAG"}),
    )
    .unwrap();
    let kinds = |v: &Value| -> Vec<String> {
        v["matches"]
            .as_array()
            .unwrap()
            .iter()
            .map(|m| m["kind"].as_str().unwrap().to_string())
            .collect()
    };
    // A sequence: once (palindrome), at the start.
    let r = s.execute("find", json!({"query": "gaattc"})).unwrap();
    assert_eq!(kinds(&r), vec!["sequence"]);
    assert_eq!(r["matches"][0]["range"], json!({"start": 0, "end": 6}));
    // On the bottom strand, and across the origin.
    let r = s.execute("find", json!({"query": "TGGT"})).unwrap();
    assert_eq!(r["matches"][0]["strand"], "forward");
    let r = s.execute("find", json!({"query": "GGGGAATT"})).unwrap();
    // Position 0 is a G too: GGG (27..30) + G (0) + AATT (1..5).
    assert_eq!(r["matches"][0]["range"], json!({"start": 27, "end": 5}));
    // Names: the feature and the primer (at its binding site); case doesn't matter.
    let r = s.execute("find", json!({"query": "LAC"})).unwrap();
    // In position order: the primer binds at 5, the feature starts at 6.
    assert_eq!(kinds(&r), vec!["primer", "feature"]);
    assert_eq!(r["matches"][0]["range"], json!({"start": 5, "end": 21}));
    // Mismatches.
    let r = s
        .execute(
            "find",
            json!({"query": "GAATTA", "mismatches": 1, "names": false}),
        )
        .unwrap();
    assert_eq!(r["matches"][0]["mismatches"], 1);
    // Nothing, or nonsense: no matches, or a clear error.
    let r = s.execute("find", json!({"query": "TTTTTTTT"})).unwrap();
    assert_eq!(kinds(&r), Vec::<String>::new());
    assert!(s.execute("find", json!({"query": "  "})).is_err());
    assert!(
        s.execute("find", json!({"query": "A", "mismatches": 9}))
            .is_err()
    );
}

#[test]
fn saving_without_history() {
    let mut s = common::session();
    s.execute("open", json!({"path": "history.dna"})).unwrap();
    let full = s.platform().read("history.dna").unwrap().len();
    // A copy without history: smaller, and the open document is untouched.
    let r = s
        .execute("save", json!({"path": "lean.dna", "without_history": true}))
        .unwrap();
    assert_eq!(r["copy"], true);
    assert!(r["bytes"].as_u64().unwrap() < full as u64);
    assert_eq!(s.document().unwrap().path.as_deref(), Some("history.dna"));
    assert!(!s.document().unwrap().dirty);
    assert!(s.molecule().unwrap().history().is_some());
    // Its own file isn't a copy.
    assert!(
        s.execute(
            "save",
            json!({"path": "history.dna", "without_history": true})
        )
        .is_err()
    );
    s.execute("open", json!({"path": "lean.dna"})).unwrap();
    assert!(s.execute("history_tree", json!({})).unwrap()["steps"].is_null());
    assert_eq!(s.molecule().unwrap().len(), 18);

    // Removing it from the document is an undoable edit.
    s.execute("switch_document", json!({"index": 0})).unwrap();
    let r = s.execute("remove_history", json!({})).unwrap();
    assert_eq!(r["removed_steps"], 2);
    assert!(s.molecule().unwrap().history().is_none());
    assert!(s.execute("remove_history", json!({})).is_err());
    s.execute("undo", json!({})).unwrap();
    assert!(s.molecule().unwrap().history().is_some());
}

/// Undo keeps small records (the changed bases and the rest of the document): after any mix of
/// edits, undoing everything gives back exactly the original, and redoing gives back the end.
#[test]
fn undo_and_redo_restore_exactly() {
    let mut s = common::session();
    s.execute(
        "new",
        json!({"name": "u", "residues": "ACGTTGCAAGGCTTACGATCGGATCCTTAGGCATGCAAGCTT", "topology": "circular"}),
    )
    .unwrap();
    let original = s.molecule().unwrap().clone();
    let mut seed = 0x0bad_5eed_u64;
    let mut next = |k: u64| {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) % k.max(1)
    };
    let mut done = 0;
    // Fewer edits than the undo limit, so everything can be undone.
    for _ in 0..snappiesnap_engine::MAX_UNDO {
        let before = s.molecule().unwrap().clone();
        let len = before.len() as u64;
        let a = next(len + 1) as usize;
        let b = next(len + 1) as usize;
        let ok = match next(7) {
            0 => s.execute("insert", json!({"at": a, "residues": "GAT"})),
            1 if len > 4 => s.execute("delete", json!({"range": {"start": a.min(b), "end": a.max(b).max(a.min(b) + 1)}})),
            2 if len > 4 => s.execute("replace", json!({"range": {"start": a, "end": b}, "residues": "CC"})),
            3 if len > 4 => s.execute(
                "add_feature",
                json!({"feature": {"name": format!("f{a}"), "kind": "misc_feature", "strand": "forward", "segments": [{"range": {"start": a.min(b), "end": a.max(b).max(a.min(b) + 1)}}]}}),
            ),
            4 => s.execute("rename", json!({"name": format!("u{a}")})),
            5 if len > 4 => s.execute("rotate_origin", json!({"origin": a})),
            _ => s.execute("reverse_complement", json!({})),
        };
        // An edit that changes nothing isn't recorded.
        if ok.is_ok() && s.molecule().unwrap() != &before {
            done += 1;
        }
    }
    assert!(done > 60, "only {done} edits applied");
    let end = s.molecule().unwrap().clone();
    let mut undone = 0;
    while s.execute("undo", json!({})).is_ok() {
        undone += 1;
    }
    assert_eq!(undone, done);
    assert_eq!(s.molecule().unwrap(), &original);
    while s.execute("redo", json!({})).is_ok() {}
    assert_eq!(s.molecule().unwrap(), &end);
}

#[test]
fn builtin_features_need_no_library() {
    let mut s = common::session();
    // loxP (34 bp, from X03453) on the bottom strand, and the lac operator (J01636).
    let loxp = "ATAACTTCGTATAGCATACATTATACGAAGTTAT";
    let rc: String = loxp
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let seq = format!("GGGGG{rc}CCCCCAATTGTGAGCGGATAACAATTGGGGGTAATACGACTCACTATAGGGAGAC");
    s.execute(
        "new",
        json!({"name": "n", "residues": seq, "topology": "linear"}),
    )
    .unwrap();
    // No library of one's own yet: the built-in features are still found.
    let d = s.execute("detect_features", json!({})).unwrap();
    let names: Vec<&str> = d["features"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|f| f["name"].as_str())
        .collect();
    assert!(names.contains(&"loxP"), "{d}");
    assert!(names.contains(&"lac operator"), "{d}");
    // An element with boundaries from a publication.
    assert!(names.contains(&"T7 promoter"), "{d}");
    let lox = d["features"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == "loxP")
        .unwrap();
    assert_eq!(lox["strand"], "reverse");
    assert_eq!(lox["range"], json!({"start": 5, "end": 39}));
    assert!(lox["source"].as_str().unwrap().starts_with("X03453"));
    // Only one's own library (the fixture's): no built-in loxP; neither source: invalid.
    let own = s
        .execute("detect_features", json!({"builtin": false}))
        .unwrap();
    assert!(
        !own["features"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["name"] == "loxP")
    );
    assert!(
        s.execute("detect_features", json!({"builtin": false, "own": false}))
            .is_err()
    );
    let a = s
        .execute("add_detected_features", json!({"names": ["loxP"]}))
        .unwrap();
    assert_eq!(a["added"], 1);
    assert_eq!(s.molecule().unwrap().features()[0].name, "loxP");
}

#[test]
fn near_matches_are_found_and_noted() {
    let mut s = common::session();
    // The built-in AmpR with two substitutions, on the bottom strand, across the origin.
    let builtin: serde_json::Value = serde_json::from_str(include_str!(
        "../../../assets/data/features/insdc_features.json"
    ))
    .unwrap();
    let amp = builtin["features"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == "AmpR")
        .unwrap()["residues"]
        .as_str()
        .unwrap()
        .to_string();
    let mut variant = amp.clone().into_bytes();
    for at in [100, 500] {
        variant[at] = if variant[at] == b'A' { b'C' } else { b'A' };
    }
    let rc: String = variant
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let seq = format!("{}{}{}", &rc[300..], "GCTAGC".repeat(50), &rc[..300]);
    let len = seq.len();
    s.execute(
        "new",
        json!({"name": "n", "residues": seq, "topology": "circular"}),
    )
    .unwrap();
    let found = |s: &mut snappiesnap_engine::Session, p: serde_json::Value| {
        s.execute("detect_features", p).unwrap()["features"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|f| f["name"] == "AmpR")
            .cloned()
            .collect::<Vec<_>>()
    };
    let hits = found(&mut s, json!({}));
    assert_eq!(hits.len(), 1, "{hits:?}");
    let h = &hits[0];
    assert_eq!(h["mismatches"], 2);
    assert_eq!(h["strand"], "reverse");
    assert_eq!(
        h["range"],
        json!({"start": len - 300, "end": amp.len() - 300})
    );
    let k = amp.len() as f64;
    assert_eq!(
        h["identity"].as_f64().unwrap(),
        ((k - 2.0) / k * 1000.0).round() / 10.0
    );
    // Exact only: not found; beyond the limit: refused.
    assert!(found(&mut s, json!({"max_mismatch_percent": 0})).is_empty());
    assert!(
        s.execute("detect_features", json!({"max_mismatch_percent": 11}))
            .is_err()
    );
    // Added with a note saying how it differs.
    s.execute("add_detected_features", json!({"names": ["AmpR"]}))
        .unwrap();
    let m = s.molecule().unwrap();
    let f = &m.features()[0];
    assert!(
        f.qualifiers.iter().any(|q| q
            .values
            .iter()
            .any(|v| v.contains("2 of 861 bases differing"))),
        "{:?}",
        f.qualifiers
    );
}

#[test]
fn fluorescent_proteins_are_found_by_translation() {
    let mut s = common::session();
    let fpbase: serde_json::Value = serde_json::from_str(include_str!(
        "../../../assets/data/features/fpbase_proteins.json"
    ))
    .unwrap();
    let egfp = fpbase["proteins"]
        .as_array()
        .unwrap()
        .iter()
        .find(|p| p["name"] == "EGFP")
        .unwrap()["protein"]
        .as_str()
        .unwrap()
        .to_string();
    // Our own codons, so only the protein can match.
    let codon = |a: char| match a {
        'A' => "GCG",
        'C' => "TGC",
        'D' => "GAT",
        'E' => "GAA",
        'F' => "TTT",
        'G' => "GGC",
        'H' => "CAT",
        'I' => "ATT",
        'K' => "AAA",
        'L' => "CTG",
        'M' => "ATG",
        'N' => "AAC",
        'P' => "CCG",
        'Q' => "CAG",
        'R' => "CGT",
        'S' => "AGC",
        'T' => "ACC",
        'V' => "GTG",
        'W' => "TGG",
        _ => "TAT",
    };
    let cds: String = egfp.chars().map(codon).collect();
    let rc: String = cds
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let seq = format!("{}{}{}", &rc[100..], "ACGTTGCA".repeat(40), &rc[..100]);
    let len = seq.len();
    s.execute(
        "new",
        json!({"name": "n", "residues": seq, "topology": "circular"}),
    )
    .unwrap();
    let d = s.execute("detect_features", json!({})).unwrap();
    let hits: Vec<_> = d["features"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|f| f["compared"] == "protein")
        .collect();
    // One place, one protein: the closest.
    assert_eq!(hits.len(), 1, "{d}");
    let h = hits[0];
    assert_eq!(h["name"], "EGFP");
    assert_eq!(h["mismatches"], 0);
    assert_eq!(h["strand"], "reverse");
    assert_eq!(
        h["range"],
        json!({"start": len - 100, "end": cds.len() - 100})
    );
    assert!(h["source"].as_str().unwrap().starts_with("FPbase"));
    // Added with its FPbase credit.
    s.execute("add_detected_features", json!({"names": ["EGFP"]}))
        .unwrap();
    let m = s.molecule().unwrap();
    let notes: Vec<&String> = m.features()[0]
        .qualifiers
        .iter()
        .flat_map(|q| q.values.iter())
        .collect();
    assert!(
        notes
            .iter()
            .any(|v| v.contains("FPbase") && v.contains("CC BY-SA")),
        "{notes:?}"
    );
    // Annotated now: not offered again.
    let again = s.execute("detect_features", json!({})).unwrap();
    assert!(
        !again["features"]
            .as_array()
            .unwrap()
            .iter()
            .any(|f| f["compared"] == "protein")
    );
}

#[test]
fn own_and_designed_features_join_the_library() {
    let mut s = common::session();
    // A design with IUPAC codes: N matches any base.
    let r = s
        .execute(
            "feature_library_add",
            json!({"name": "my tag", "residues": "GATTACA NNN CCGGAA", "kind": "misc_feature", "note": "lab design"}),
        )
        .unwrap();
    assert_eq!(
        (r["length"].as_u64(), r["replaced"].as_bool()),
        (Some(16), Some(false))
    );
    assert_eq!(s.execute("feature_library", json!({})).unwrap()["mine"], 1);
    s.execute(
        "new",
        json!({"name": "n", "residues": "CCCCCGATTACAGTCCCGGAACCCCC", "topology": "linear"}),
    )
    .unwrap();
    let only_mine = json!({"builtin": false, "max_mismatch_percent": 0});
    let d = s.execute("detect_features", only_mine.clone()).unwrap();
    let tag = d["features"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == "my tag")
        .unwrap();
    assert_eq!(tag["range"], json!({"start": 5, "end": 21}));
    assert_eq!(tag["source"], "designed by you");
    // A feature of the document joins it too, under a new name; adding it again replaces it.
    s.execute(
        "add_feature",
        json!({"feature": {
            "name": "spacer", "kind": "misc_feature", "strand": "forward",
            "segments": [{"range": {"start": 0, "end": 5}}]
        }}),
    )
    .unwrap();
    let id = s.molecule().unwrap().features()[0].id;
    let r = s
        .execute("feature_library_add", json!({"feature": id, "name": "C5"}))
        .unwrap();
    assert_eq!(r["replaced"], false);
    let r = s
        .execute("feature_library_add", json!({"feature": id, "name": "C5"}))
        .unwrap();
    assert_eq!(r["replaced"], true);
    // Rebuilding from a folder keeps them.
    s.execute(
        "feature_library_build",
        json!({"folder": "lab", "min_length": 10}),
    )
    .unwrap();
    assert_eq!(s.execute("feature_library", json!({})).unwrap()["mine"], 2);
    // Removing; bad input is a clear error.
    assert_eq!(
        s.execute("feature_library_remove", json!({"name": "my tag"}))
            .unwrap()["removed"],
        1
    );
    assert!(
        s.execute("feature_library_remove", json!({"name": "my tag"}))
            .is_err()
    );
    for bad in [
        json!({"name": "x", "residues": "ACGTZ"}),
        json!({"name": "", "residues": "ACGT"}),
        json!({"residues": "ACGT"}),
        json!({"feature": 9999}),
    ] {
        assert!(
            s.execute("feature_library_add", bad.clone()).is_err(),
            "{bad}"
        );
    }
}

#[test]
fn downloaded_vector_library_needs_consent_and_stays_separate() {
    let mut s = common::session();
    // Without consent: the terms, and nothing downloaded.
    let e = s
        .execute(
            "feature_library_download",
            json!({"url": common::VECTORS_URL}),
        )
        .unwrap_err();
    assert!(e.message.contains("may not be shared"), "{}", e.message);
    let r = s
        .execute(
            "feature_library_download",
            json!({"accept": true, "url": common::VECTORS_URL}),
        )
        .unwrap();
    assert_eq!(
        (r["records"].as_u64(), r["features"].as_u64()),
        (Some(2), Some(2))
    );
    assert_eq!(
        s.execute("feature_library", json!({})).unwrap()["downloaded"],
        2
    );
    s.execute(
        "new",
        json!({"name": "n", "residues": "AAAAAGACTACAAGGACGACGATGACAAGAAAAA", "topology": "linear"}),
    )
    .unwrap();
    let flag = |s: &mut snappiesnap_engine::Session, p: serde_json::Value| {
        s.execute("detect_features", p).unwrap()["features"]
            .as_array()
            .unwrap()
            .iter()
            .find(|f| f["name"] == "FLAG tag")
            .cloned()
    };
    let f = flag(&mut s, json!({})).unwrap();
    assert!(
        f["source"].as_str().unwrap().ends_with("(downloaded)"),
        "{f}"
    );
    // Left out on request, and gone once deleted.
    assert!(flag(&mut s, json!({"downloaded": false})).is_none());
    assert_eq!(
        s.execute("feature_library_delete_download", json!({}))
            .unwrap()["deleted"],
        true
    );
    assert!(flag(&mut s, json!({})).is_none());
    assert_eq!(
        s.execute("feature_library", json!({})).unwrap()["downloaded"],
        0
    );
}

#[test]
fn embl_and_gff3_save_open_and_annotate() {
    let mut s = common::session();
    s.execute(
        "new",
        json!({"name": "demo", "residues": "ACGTACGTACGTACGTACGT", "topology": "circular"}),
    )
    .unwrap();
    // Annotations without a sequence: one undoable edit, no duplicates on a second import.
    let r = s
        .execute("import_annotations", json!({"path": "annotations.gff3"}))
        .unwrap();
    assert_eq!(
        (r["seqid"].as_str(), r["added"].as_u64()),
        (Some("demo"), Some(2))
    );
    let cds = s
        .molecule()
        .unwrap()
        .features()
        .iter()
        .find(|f| f.name == "demo cds")
        .cloned()
        .unwrap();
    assert_eq!(cds.strand, snappiesnap_doc::Strand::Reverse);
    assert_eq!(cds.segments[0].range, snappiesnap_doc::Range::new(2, 9));
    let again = s
        .execute("import_annotations", json!({"path": "annotations.gff3"}))
        .unwrap();
    assert_eq!(
        (again["added"].as_u64(), again["skipped"].as_u64()),
        (Some(0), Some(2))
    );
    // The second import changed nothing, so one undo removes the annotations.
    s.execute("undo", json!({})).unwrap();
    assert!(s.molecule().unwrap().features().is_empty());
    s.execute("redo", json!({})).unwrap();
    // A wrong seqid or a file that isn't GFF3: clear errors.
    assert!(
        s.execute(
            "import_annotations",
            json!({"path": "annotations.gff3", "seqid": "other"})
        )
        .is_err()
    );
    assert!(
        s.execute("import_annotations", json!({"path": "withrefm.txt"}))
            .is_err()
    );
    // Saved as EMBL and as GFF3, and opened again: the same sequence and features.
    let before = s.molecule().unwrap().clone();
    for path in ["demo.embl", "demo.gff3", "demo.sbol"] {
        s.execute("save", json!({"path": path})).unwrap();
        s.execute("close", json!({})).ok();
        let r = s.execute("open", json!({"path": path})).unwrap();
        let m = s.molecule().unwrap();
        assert_eq!(m.seq(), before.seq(), "{path}: {r}");
        assert_eq!(m.topology(), before.topology(), "{path}");
        let names = |m: &snappiesnap_doc::Molecule| {
            let mut v: Vec<(String, Vec<snappiesnap_doc::Range>, snappiesnap_doc::Strand)> = m
                .features()
                .iter()
                .map(|f| {
                    (
                        f.name.clone(),
                        f.segments.iter().map(|x| x.range).collect(),
                        f.strand,
                    )
                })
                .collect();
            v.sort_by(|a, b| a.0.cmp(&b.0));
            v
        };
        assert_eq!(names(m), names(&before), "{path}");
    }
}

#[test]
fn blast_addresses_for_dna_and_protein() {
    let mut s = common::session();
    s.execute(
        "new",
        json!({"name": "my seq", "residues": "ATGGCCAAGTTTCCCGGGTAA", "topology": "circular"}),
    )
    .unwrap();
    let url = |s: &mut snappiesnap_engine::Session, p: serde_json::Value| {
        s.execute("blast_search", p).unwrap()
    };
    // The whole sequence, nucleotides by default.
    let r = url(&mut s, json!({}));
    assert_eq!(r["program"], "blastn");
    let u = r["url"].as_str().unwrap();
    assert!(u.starts_with("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROGRAM=blastn&QUERY=ATGGCCAAGTTTCCCGGGTAA&JOB_TITLE=my%20seq%201..21"), "{u}");
    // The selection, across the origin, against proteins and translated nucleotides.
    s.execute("select", json!({"range": {"start": 18, "end": 3}}))
        .unwrap();
    let r = url(&mut s, json!({"against": "protein"}));
    assert_eq!(
        (r["program"].as_str(), r["query_length"].as_u64()),
        (Some("blastx"), Some(6))
    );
    assert!(r["url"].as_str().unwrap().contains("QUERY=TAAATG&"));
    assert_eq!(
        url(&mut s, json!({"against": "translated"}))["program"],
        "tblastx"
    );
    let r = url(
        &mut s,
        json!({"database": "refseq_rna", "range": {"start": 0, "end": 9}}),
    );
    assert!(r["url"].as_str().unwrap().ends_with("&DATABASE=refseq_rna"));
    assert!(
        s.execute("blast_search", json!({"database": "nt&x=1"}))
            .is_err()
    );
    // A protein: blastp or tblastn, never straight against nucleotides.
    s.execute(
        "new",
        json!({"name": "p", "residues": "MAKFPG", "alphabet": "protein"}),
    )
    .unwrap();
    assert_eq!(url(&mut s, json!({}))["program"], "blastp");
    assert_eq!(
        url(&mut s, json!({"against": "translated"}))["program"],
        "tblastn"
    );
    assert!(
        s.execute("blast_search", json!({"against": "nucleotide"}))
            .is_err()
    );
    // Too long for a web address: a clear error.
    s.execute(
        "new",
        json!({"name": "long", "residues": "ACGT".repeat(2000)}),
    )
    .unwrap();
    let e = s.execute("blast_search", json!({})).unwrap_err();
    assert!(e.message.contains("8000 residues"), "{}", e.message);
}

#[test]
fn feature_library_from_own_files() {
    let mut s = common::session();
    // Learn from the folder (one file, one feature).
    let r = s
        .execute(
            "feature_library_build",
            json!({"folder": "lab", "min_length": 10}),
        )
        .unwrap();
    assert_eq!(
        (r["features"].as_u64(), r["files"].as_u64()),
        (Some(1), Some(1))
    );
    let lib = s.execute("feature_library", json!({})).unwrap();
    assert_eq!(lib["features"][0]["name"], "start block");
    assert_eq!(lib["features"][0]["length"], 16);
    // Its feature in a new sequence: on the bottom strand, across the origin.
    let block = "ACGTACGTAAAACCCC";
    let rc: String = block
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    let seq = format!("{}GGGGGGGGGG{}", &rc[6..], &rc[..6]);
    s.execute(
        "new",
        json!({"name": "n", "residues": seq, "topology": "circular"}),
    )
    .unwrap();
    let d = s.execute("detect_features", json!({})).unwrap();
    let f = &d["features"][0];
    // A plain GenBank location has no direction, so neither does the find.
    assert_eq!(f["strand"], "none");
    assert_eq!(f["range"], json!({"start": 20, "end": 10}));
    // Added as one edit; found again, it's already there.
    let a = s.execute("add_detected_features", json!({})).unwrap();
    assert_eq!(a["added"], 1);
    assert_eq!(
        s.execute("detect_features", json!({})).unwrap()["features"],
        json!([])
    );
    s.execute("undo", json!({})).unwrap();
    assert!(s.molecule().unwrap().features().is_empty());
    // Names that aren't found add nothing; no folder is a clear error.
    let a = s
        .execute("add_detected_features", json!({"names": ["nothing"]}))
        .unwrap();
    assert_eq!(a["added"], 0);
    assert!(
        s.execute("feature_library_build", json!({"folder": "nowhere"}))
            .is_err()
    );
}

#[test]
fn sequencing_reads_open_trim_and_align() {
    let mut s = common::session();
    let opened = run(&mut s, "open", json!({"path": "read.ab1"}));
    assert_eq!(opened["length"], 48);
    let doc = s.document().unwrap();
    assert_eq!(doc.molecule.seq().as_bytes()[..4], *b"NTGA");
    let got = snappiesnap_engine::api::call(&mut s, "doc.get", json!({})).unwrap();
    assert_eq!(got["trace"]["sample_name"], "read 1 fwd");
    assert!(got["molecule"].get("trace").is_none());

    let info = run(
        &mut s,
        "trace_info",
        json!({"from": 4, "to": 6, "signal": true}),
    );
    assert_eq!(info["trim"], json!({"start": 4, "end": 44}));
    assert_eq!(info["quality"]["q30"], 40);
    assert_eq!(
        info["calls"][0],
        json!({"position": 4, "base": "G", "quality": 38, "peak": 54})
    );
    // Signal from 6 samples before the first peak to 6 after the last (54, 66).
    assert_eq!(info["signal"]["start"], 48);
    assert_eq!(info["signal"]["g"].as_array().unwrap().len(), 66 + 7 - 48);
    assert_eq!(info["signal"]["g"][6], 900);
    assert_eq!(
        code(&mut s, "trace_info", json!({"from": 40, "to": 49})),
        ErrorCode::InvalidInput
    );
    assert_eq!(
        code(&mut s, "trace_info", json!({"cutoff": 1.5})),
        ErrorCode::InvalidParams
    );

    // Reverse-complementing mirrors the trace: residue 0 is now the old last base.
    run(&mut s, "reverse_complement", json!({}));
    let info = run(&mut s, "trace_info", json!({"from": 0, "to": 1}));
    assert_eq!(info["reversed"], true);
    assert_eq!(info["calls"][0]["peak"], 12 * 48 - 1 - (12 * 47 + 6));
    assert_eq!(info["trim"], json!({"start": 4, "end": 44}));
    run(&mut s, "undo", json!({}));

    let trimmed = run(&mut s, "trace_trim", json!({}));
    assert_eq!(
        trimmed,
        json!({"removed_start": 4, "removed_end": 4, "length": 40})
    );
    assert_eq!(residues(&mut s, 0, 7), "GATTACA");
    let info = run(&mut s, "trace_info", json!({"from": 0, "to": 1}));
    assert_eq!(info["calls"][0]["peak"], 54);
    run(&mut s, "undo", json!({}));
    assert_eq!(s.document().unwrap().molecule.len(), 48);

    // An edited read can't be saved over its .ab1 file; .dna and FASTA work, with a warning.
    run(&mut s, "trace_trim", json!({}));
    assert_eq!(code(&mut s, "save", json!({})), ErrorCode::Unsupported);
    let saved = run(&mut s, "save", json!({"path": "read.fa"}));
    assert!(saved.to_string().contains("chromatogram"), "{saved}");

    // Aligned to the template, the read is trimmed first and matches exactly.
    run(
        &mut s,
        "new",
        json!({"name": "t", "residues": "GATTACAGATTACACCGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACG", "discard_changes": true}),
    );
    let a = run(&mut s, "align_to_reference", json!({"file": "read.ab1"}));
    assert_eq!(
        a["trimmed"],
        json!({"start": 4, "end": 44, "read_length": 48})
    );
    assert_eq!(a["identity"], 1.0);
    let whole = run(
        &mut s,
        "align_to_reference",
        json!({"file": "read.ab1", "trim": false}),
    );
    assert!(whole["identity"].as_f64().unwrap() < 1.0);
    assert!(whole.get("trimmed").is_none());
    assert_eq!(code(&mut s, "trace_trim", json!({})), ErrorCode::NotEnabled);
}

/// `tests/scenarios/data/read1.ab1` is the synthetic read made by `common::read_ab1` (rewrite
/// it with `SNAPPIESNAP_WRITE_FIXTURES=1`).
#[test]
fn scenario_read_is_the_synthetic_fixture() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/scenarios/data/read1.ab1");
    if std::env::var_os("SNAPPIESNAP_WRITE_FIXTURES").is_some() {
        std::fs::write(&path, common::read_ab1()).unwrap();
    }
    assert_eq!(std::fs::read(&path).unwrap(), common::read_ab1());
}

#[test]
fn experiments_save_and_reopen() {
    let mut s = common::session();
    // The fixture: example.dna is found next to the runs/ folder; pCopy only exists as a copy.
    let opened = run(
        &mut s,
        "experiment_open",
        json!({"path": "runs/experiment.sslab"}),
    );
    assert_eq!(opened["name"], "Fixture run");
    assert_eq!(opened["agarose"], 1.2);
    assert_eq!(opened["constructs"][0]["source"], "file");
    assert_eq!(opened["constructs"][0]["path"], "example.dna");
    assert_eq!(opened["constructs"][1]["source"], "copy");
    assert_eq!(
        opened["warnings"],
        json!([]),
        "a construct without a file needs no warning"
    );
    let docs = run(&mut s, "documents", json!({}))["documents"].clone();
    assert_eq!(docs.as_array().unwrap().len(), 2);
    let (a, b) = (
        opened["constructs"][0]["document"].clone(),
        opened["constructs"][1]["document"].clone(),
    );
    assert_eq!(
        opened["lanes"],
        json!([
            {"document": a, "method": {"kind": "uncut"}},
            {"document": b, "method": {"kind": "digest", "enzymes": ["EcoRI"]}},
            {"document": a, "method": null}
        ])
    );
    // The lanes run as they are.
    let gel = run(&mut s, "gel_experiment", json!({"lanes": opened["lanes"]}));
    assert_eq!(
        gel["lanes"][2]["bands"].as_array().unwrap().len(),
        1,
        "{gel}"
    );

    // Opening it again reuses both tabs (example.dna's, and the unsaved copy is a new tab).
    let again = run(
        &mut s,
        "experiment_open",
        json!({"path": "runs/experiment.sslab"}),
    );
    assert_eq!(again["constructs"][0]["source"], "open");
    assert_eq!(again["constructs"][0]["document"], a);

    // Save: the copy-only construct and an edited file are kept as copies, with warnings.
    let a = a.as_u64().unwrap();
    let index = |s: &mut Session, id: u64| {
        run(s, "documents", json!({}))["documents"]
            .as_array()
            .unwrap()
            .iter()
            .find(|d| d["id"] == id)
            .unwrap()["index"]
            .clone()
    };
    let i = index(&mut s, a);
    run(&mut s, "switch_document", json!({"index": i}));
    run(&mut s, "insert", json!({"at": 0, "residues": "GG"}));
    let saved = run(
        &mut s,
        "experiment_save",
        json!({"path": "out/run2.sslab", "name": "Run 2", "ladder": "100 bp",
               "lanes": [{"document": a, "method": {"kind": "uncut"}}, {"document": b}]}),
    );
    assert_eq!(saved["constructs"], 2);
    assert_eq!(saved["warnings"].as_array().unwrap().len(), 2, "{saved}");
    let bytes = s.platform().read("out/run2.sslab").unwrap();
    let f: Value = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(f["format"], "snappiesnap-experiment");
    assert_eq!(f["constructs"][0]["file"], "../example.dna");
    assert_eq!(f["constructs"][0]["length"], 18);
    assert!(f["constructs"][1].get("file").is_none());
    assert!(
        f["constructs"][0]["genbank"]
            .as_str()
            .unwrap()
            .starts_with("LOCUS")
    );

    // In a fresh session the file on disk lacks the edit: the copy (18 bp) opens instead.
    let mut fresh = common::session();
    fresh.platform().write("out/run2.sslab", &bytes).unwrap();
    let r = run(
        &mut fresh,
        "experiment_open",
        json!({"path": "out/run2.sslab"}),
    );
    assert_eq!(r["constructs"][0]["source"], "copy");
    assert!(
        r["warnings"][0].as_str().unwrap().contains("has changed"),
        "{r}"
    );
    assert_eq!(r["ladder"], "100 bp");
    let id = r["constructs"][0]["document"].as_u64().unwrap();
    let copy = fresh.document_by_id(id).unwrap();
    assert_eq!(copy.molecule.len(), 18);
    assert!(copy.path.is_none() && copy.dirty);

    // Bad input changes nothing.
    let mut t = common::session();
    assert_eq!(
        code(&mut t, "experiment_open", json!({"path": "example.dna"})),
        ErrorCode::InvalidInput
    );
    assert_eq!(
        code(&mut t, "experiment_open", json!({"path": "missing.sslab"})),
        ErrorCode::Io
    );
    assert!(t.documents().is_empty());
    assert_eq!(
        code(
            &mut s,
            "experiment_save",
            json!({"path": "x.dna", "lanes": []})
        ),
        ErrorCode::InvalidParams
    );
    assert_eq!(
        code(
            &mut s,
            "experiment_save",
            json!({"path": "x.sslab", "lanes": [{"document": 999}]})
        ),
        ErrorCode::NotFound
    );
}

#[test]
fn go_to_and_sequence_info() {
    let mut s = common::session();
    run(
        &mut s,
        "new",
        json!({"name": "p", "residues": "GGGGCCCCAAAATTTTACGTACGT", "topology": "circular"}),
    );
    let r = run(&mut s, "go_to", json!({"location": "5-8"}));
    assert_eq!(r, json!({"selection": {"start": 4, "end": 8}, "length": 4}));
    let info = run(&mut s, "sequence_info", json!({}));
    assert_eq!(info["length"], 4);
    assert_eq!(info["gc_percent"], 100.0);
    // CCCC: 4 × 289.18 − 61.96 = 1094.76; GGGG: 4 × 329.21 − 61.96 = 1254.88.
    assert_eq!(info["weight"]["single_strand"], 1094.76);
    assert_eq!(info["weight"]["double_strand"], 2349.64);
    assert!(info.get("tm").is_none(), "4 nt: no Tm");
    // Across the origin.
    let r = run(&mut s, "go_to", json!({"location": "21 to 4"}));
    assert_eq!(r["length"], 8);
    let info = run(&mut s, "sequence_info", json!({}));
    assert_eq!(info["gc_percent"], 75.0);
    assert!(info["tm"].as_f64().is_some(), "{info}");
    // The whole circle: no ends taken off.
    let whole = run(
        &mut s,
        "sequence_info",
        json!({"range": {"start": 0, "end": 24}}),
    );
    assert_eq!(whole["whole"], true);
    assert_eq!(whole["weight"]["circular"], true);
    assert_eq!(
        code(&mut s, "go_to", json!({"location": "25"})),
        ErrorCode::InvalidInput
    );
    assert_eq!(
        code(&mut s, "go_to", json!({"location": "x"})),
        ErrorCode::InvalidInput
    );
    run(
        &mut s,
        "new",
        json!({"name": "prot", "residues": "MKV", "alphabet": "protein"}),
    );
    let info = run(&mut s, "sequence_info", json!({}));
    assert!(info["mass"].as_f64().unwrap() > 300.0);
    assert!(info.get("gc_percent").is_none());
}

#[test]
fn diagnostic_digest_tells_a_clone_from_its_vector() {
    let mut s = common::session();
    let filler = |n: usize| "ACGTTGCA".repeat(n / 8 + 1)[..n].to_string();
    let vector = format!("GAATTC{}AAGCTT{}", filler(2994), filler(2988));
    let clone = format!(
        "GAATTC{}GGATCC{}AAGCTT{}",
        filler(2994),
        filler(1494),
        filler(2988)
    );
    run(
        &mut s,
        "new",
        json!({"name": "vector", "residues": vector, "topology": "circular"}),
    );
    run(
        &mut s,
        "new",
        json!({"name": "clone", "residues": clone, "topology": "circular"}),
    );
    let ids: Vec<u64> = run(&mut s, "documents", json!({}))["documents"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| d["id"].as_u64().unwrap())
        .collect();
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "limit": 50}),
    );
    let suggestions = r["suggestions"].as_array().unwrap();
    assert!(!suggestions.is_empty(), "{r}");
    let pair = suggestions
        .iter()
        .find(|x| x["enzymes"] == json!(["EcoRI", "HindIII"]))
        .unwrap_or_else(|| panic!("EcoRI + HindIII missing: {r}"));
    assert_eq!(pair["lanes"][0]["fragments"], json!([3000, 2994]));
    assert_eq!(pair["lanes"][1]["fragments"], json!([4500, 2994]));
    assert_eq!(pair["lanes"][1]["name"], "clone");
    // Single enzymes only.
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "max_enzymes": 1}),
    );
    assert!(
        r["suggestions"]
            .as_array()
            .unwrap()
            .iter()
            .all(|x| x["enzymes"].as_array().unwrap().len() == 1)
    );
    // EcoRI alone linearises both: 5,994 against 7,494 bp, far enough apart.
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "enzymes": ["EcoRI"]}),
    );
    assert_eq!(r["suggestions"][0]["lanes"][1]["fragments"], json!([7494]));
    // BamHI leaves the vector uncut: nothing.
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "enzymes": ["BamHI"]}),
    );
    assert_eq!(r["suggestions"], json!([]));
    // An enzyme set is what the lab holds, so a search of it only suggests runnable digests.
    // "Lab freezer" is EcoRI + SmaI; SmaI doesn't cut here, so EcoRI alone is all there is,
    // and the EcoRI + HindIII pair above is gone.
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "set": "Lab freezer", "limit": 50}),
    );
    assert_eq!(r["set"], "Lab freezer");
    assert_eq!(r["in_stock_searched"], 2);
    assert_eq!(r["enzymes_searched"], 2);
    let suggestions = r["suggestions"].as_array().unwrap();
    assert!(!suggestions.is_empty(), "{r}");
    assert!(
        suggestions
            .iter()
            .all(|x| x["in_stock"] == json!(true) && x["enzymes"] == json!(["EcoRI"])),
        "{r}"
    );

    // Widened on request: what the lab can't run comes back, marked, and after what it can.
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "set": "Lab freezer", "include_out_of_stock": true, "limit": 50}),
    );
    let suggestions = r["suggestions"].as_array().unwrap();
    assert!(r["enzymes_searched"].as_u64().unwrap() > 2, "{r}");
    let first_out = suggestions
        .iter()
        .position(|x| x["in_stock"] == json!(false))
        .unwrap_or(suggestions.len());
    let last_in = suggestions
        .iter()
        .rposition(|x| x["in_stock"] == json!(true))
        .unwrap_or(0);
    assert!(last_in < first_out, "in-stock suggestions come first: {r}");
    let pair = suggestions
        .iter()
        .find(|x| x["enzymes"] == json!(["EcoRI", "HindIII"]))
        .unwrap_or_else(|| panic!("EcoRI + HindIII missing: {r}"));
    assert_eq!(pair["in_stock"], json!(false));

    // An enzyme the lab owns under a name the database lacks is reported, not passed over.
    run(
        &mut s,
        "enzyme_set_create",
        json!({"name": "Fridge", "enzymes": ["EcoRI", "NotAnEnzymeI"]}),
    );
    let r = run(
        &mut s,
        "diagnostic_digest",
        json!({"documents": ids, "set": "Fridge"}),
    );
    assert_eq!(r["missing_from_database"], json!(["NotAnEnzymeI"]));
    assert_eq!(r["in_stock_searched"], 1);
    assert_eq!(
        code(
            &mut s,
            "diagnostic_digest",
            json!({"documents": ids, "set": "no such set"})
        ),
        ErrorCode::NotFound
    );
    assert_eq!(
        code(&mut s, "diagnostic_digest", json!({"documents": [999]})),
        ErrorCode::NotFound
    );
    assert_eq!(
        code(&mut s, "diagnostic_digest", json!({"documents": []})),
        ErrorCode::InvalidParams
    );
    assert_eq!(
        code(
            &mut s,
            "diagnostic_digest",
            json!({"documents": ids, "enzymes": ["NoSuchI"]})
        ),
        ErrorCode::NotFound
    );
}

#[test]
fn check_reads_against_a_construct() {
    let mut s = common::session();
    let template = "GATTACAGATTACACCGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCG";
    // A second read of the template's first 60 bases with a substitution at position 20 (A > C).
    let mut second: Vec<u8> = template.as_bytes()[..60].to_vec();
    assert_eq!(second[20], b'C');
    second[20] = b'A';
    run(
        &mut s,
        "new",
        json!({"name": "read2", "residues": String::from_utf8(second).unwrap()}),
    );
    let read2 = s.active_document_id().unwrap();
    run(
        &mut s,
        "new",
        json!({"name": "pT", "residues": template, "topology": "circular"}),
    );
    let r = run(
        &mut s,
        "check_reads",
        json!({"files": ["read.ab1"], "documents": [read2]}),
    );
    assert_eq!(r["placed"], 2, "{r}");
    // read.ab1 is trimmed to the template's first 40 bases.
    assert_eq!(r["reads"][0]["kept"], json!([4, 44]));
    assert_eq!(
        r["reads"][0]["source"],
        json!({"file": "read.ab1", "trace": true})
    );
    assert_eq!(
        (r["reads"][0]["start"].clone(), r["reads"][0]["end"].clone()),
        (json!(0), json!(40))
    );
    assert_eq!(r["covered"], 60);
    assert_eq!(r["uncovered"], json!([{"start": 60, "end": 78}]));
    // Base 20 is covered by both reads; only read2 shows the change: disputed.
    let f = &r["findings"][0];
    assert_eq!(f["position"], 20);
    assert_eq!(f["verdict"], "disputed");
    assert_eq!(
        f["reads"],
        json!([{"read": 1, "read_position": 20, "quality": null}])
    );
    // read.ab1 (trimmed by 4) has construct base 20 at its base 24.
    assert_eq!(
        f["covered_by"],
        json!([{"read": 0, "read_position": 24}, {"read": 1, "read_position": 20}])
    );
    assert_eq!(
        code(&mut s, "check_reads", json!({})),
        ErrorCode::InvalidParams
    );
    assert_eq!(
        code(&mut s, "check_reads", json!({"files": ["missing.ab1"]})),
        ErrorCode::Io
    );
    assert_eq!(
        code(
            &mut s,
            "check_reads",
            json!({"document": 999, "files": ["read.ab1"]})
        ),
        ErrorCode::NotFound
    );
}

#[test]
fn reads_experiments_save_and_reopen() {
    let mut s = common::session();
    run(&mut s, "open", json!({"path": "example.dna"}));
    let id = s.active_document_id().unwrap();
    let saved = run(
        &mut s,
        "experiment_save",
        json!({"path": "runs/check.sslab", "name": "Check", "kind": "reads", "reference": id,
               "reads": ["read.ab1", "gone.ab1"]}),
    );
    assert_eq!(saved["reads"], 2);
    assert!(
        saved["warnings"][0].as_str().unwrap().contains("gone.ab1"),
        "{saved}"
    );
    let f: Value = serde_json::from_slice(&s.platform().read("runs/check.sslab").unwrap()).unwrap();
    assert_eq!(f["kind"], "reads");
    assert_eq!(
        f["reads"][0],
        json!({"name": "read", "file": "../read.ab1", "original_path": "read.ab1"})
    );
    assert_eq!(f["lanes"], json!([]));
    let opened = run(
        &mut s,
        "experiment_open",
        json!({"path": "runs/check.sslab"}),
    );
    assert_eq!(opened["kind"], "reads");
    assert_eq!(opened["reference"], id);
    assert_eq!(opened["reads"], json!(["read.ab1"]));
    assert!(
        opened["warnings"][0].as_str().unwrap().contains("gone"),
        "{opened}"
    );
    assert_eq!(
        code(
            &mut s,
            "experiment_save",
            json!({"path": "x.sslab", "kind": "reads", "reads": []})
        ),
        ErrorCode::InvalidParams
    );
    assert_eq!(
        code(
            &mut s,
            "experiment_save",
            json!({"path": "x.sslab", "kind": "western", "lanes": []})
        ),
        ErrorCode::InvalidParams
    );
}

#[test]
fn matched_reads_list_open_and_check() {
    let mut s = common::session();
    run(&mut s, "open", json!({"path": "matched.dna"}));
    let construct = s.active_document_id().unwrap();
    let list = run(&mut s, "matched_reads", json!({}));
    assert_eq!(
        list["reads"],
        json!([{"index": 0, "name": "read 1", "length": 30, "trace": true, "trimmed": [0, 30]}])
    );
    // Checked against the construct: placed at 10..40, the A at construct base 20.
    let r = run(&mut s, "check_reads", json!({"matched": true}));
    assert_eq!(
        r["reads"][0]["source"],
        json!({"matched": 0, "trace": true})
    );
    assert_eq!(
        (r["reads"][0]["start"].clone(), r["reads"][0]["end"].clone()),
        (json!(10), json!(40))
    );
    assert_eq!(r["findings"][0]["position"], 20);
    assert_eq!(r["findings"][0]["reads"][0]["quality"], 40);
    // Opened in a tab of its own, with its trace.
    let opened = run(&mut s, "open_matched_read", json!({"index": 0}));
    assert_eq!(opened["trace"], true);
    let info = run(&mut s, "trace_info", json!({"from": 0, "to": 1}));
    assert_eq!(info["calls"][0]["peak"], 5);
    assert_eq!(
        code(
            &mut s,
            "open_matched_read",
            json!({"index": 3, "document": construct})
        ),
        ErrorCode::NotFound
    );
    // A document without matched reads lists none, and can't be checked with them alone.
    run(
        &mut s,
        "new",
        json!({"name": "plain", "residues": "ACGTACGTACGT"}),
    );
    assert_eq!(run(&mut s, "matched_reads", json!({}))["reads"], json!([]));
    assert_eq!(
        code(&mut s, "check_reads", json!({"matched": true})),
        ErrorCode::InvalidParams
    );
}

/// `tests/scenarios/data/matched.dna` is the synthetic file made by `common::matched_dna`
/// (rewrite it with `SNAPPIESNAP_WRITE_FIXTURES=1`).
#[test]
fn scenario_matched_dna_is_the_synthetic_fixture() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/scenarios/data/matched.dna");
    if std::env::var_os("SNAPPIESNAP_WRITE_FIXTURES").is_some() {
        std::fs::write(&path, common::matched_dna()).unwrap();
    }
    assert_eq!(std::fs::read(&path).unwrap(), common::matched_dna());
}

#[test]
fn near_matches_with_inserted_and_missing_bases() {
    let mut s = common::session();
    let builtin: serde_json::Value = serde_json::from_str(include_str!(
        "../../../assets/data/features/insdc_features.json"
    ))
    .unwrap();
    let amp = builtin["features"]
        .as_array()
        .unwrap()
        .iter()
        .find(|f| f["name"] == "AmpR")
        .unwrap()["residues"]
        .as_str()
        .unwrap()
        .to_string();
    // Two bases inserted after 200 and one base missing at 600: 862 + 2 - 1 bases.
    let variant = format!("{}GG{}{}", &amp[..200], &amp[200..600], &amp[601..]);
    let rc: String = variant
        .bytes()
        .rev()
        .map(|b| match b {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            _ => 'A',
        })
        .collect();
    // On the bottom strand, across the origin of a circle.
    let seq = format!("{}{}{}", &rc[300..], "GCTAGC".repeat(50), &rc[..300]);
    let len = seq.len();
    s.execute(
        "new",
        json!({"name": "n", "residues": seq, "topology": "circular"}),
    )
    .unwrap();
    let found = |s: &mut snappiesnap_engine::Session, p: serde_json::Value| {
        s.execute("detect_features", p).unwrap()["features"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|f| f["name"] == "AmpR")
            .cloned()
            .collect::<Vec<_>>()
    };
    let hits = found(&mut s, json!({}));
    assert_eq!(hits.len(), 1, "{hits:?}");
    let h = &hits[0];
    assert_eq!(h["mismatches"], 3, "{h}");
    assert_eq!(h["insertions"], 2, "{h}");
    assert_eq!(h["deletions"], 1, "{h}");
    assert_eq!(h["strand"], "reverse");
    // The copy (863 bases) wraps the origin: 300 bases at the end, the rest at the start.
    assert_eq!(
        h["range"],
        json!({"start": len - 300, "end": variant.len() - 300})
    );
    // Substitutions only: not found.
    assert!(found(&mut s, json!({"indels": false})).is_empty());
    let a = s
        .execute("add_detected_features", json!({"names": ["AmpR"]}))
        .unwrap();
    assert_eq!(a["added"], 1);
    let f = &s.molecule().unwrap().features()[0];
    let notes: Vec<&String> = f
        .qualifiers
        .iter()
        .filter(|q| q.name == "note")
        .flat_map(|q| &q.values)
        .collect();
    assert!(
        notes
            .iter()
            .any(|n| n.contains("3 of 861 bases differing (2 inserted, 1 missing)")),
        "{notes:?}"
    );
}

/// The eight Gateway att sites are built in and found where they sit, each in the orientation
/// the dataset holds, with the core its number implies (Hartley et al. 2000).
#[test]
fn gateway_att_sites_are_built_in() {
    const DATASET: &str = include_str!("../../../assets/data/features/insdc_features.json");
    let data: Value = serde_json::from_str(DATASET).unwrap();
    let entry = |name: &str| -> (String, String) {
        let f = data["features"]
            .as_array()
            .unwrap()
            .iter()
            .find(|f| f["name"] == name)
            .unwrap_or_else(|| panic!("{name} is not a built-in feature"));
        (
            f["residues"].as_str().unwrap().to_string(),
            f["accession"].as_str().unwrap().to_string(),
        )
    };
    let sites = [
        "attB1", "attB2", "attL1", "attL2", "attR1", "attR2", "attP1", "attP2",
    ];
    // Each site carries the core of its number, and comes to its expected length.
    for name in sites {
        let (residues, _) = entry(name);
        let core = if name.ends_with('1') {
            "TTTGTACAAAAAAGC"
        } else {
            "CTTTCTTGTACAAAG"
        };
        assert!(residues.contains(core), "{name} has no {core} core");
        let expected = match &name[..4] {
            "attB" => 25,
            "attL" => 100,
            "attR" => 125,
            _ => 232,
        };
        assert_eq!(residues.len(), expected, "{name} is the wrong length");
    }
    // A synthetic construct holding every site in turn, 10 bases apart.
    let spacer = "GCTAGCTAGC";
    let mut seq = String::from(spacer);
    let mut placed = Vec::new();
    for name in sites {
        let (residues, accession) = entry(name);
        placed.push((name, seq.len(), seq.len() + residues.len(), accession));
        seq.push_str(&residues);
        seq.push_str(spacer);
    }
    let mut s = common::session();
    run(
        &mut s,
        "new",
        json!({"name": "att", "residues": seq, "topology": "linear"}),
    );
    let d = run(&mut s, "detect_features", json!({"own": false}));
    let found = d["features"].as_array().unwrap();
    for (name, start, end, accession) in placed {
        let range = json!({"start": start, "end": end});
        let hit = found
            .iter()
            .find(|f| f["name"] == name && f["range"] == range)
            .unwrap_or_else(|| panic!("{name} was not detected at {range}: {d}"));
        assert_eq!(hit["strand"], "forward", "{name}");
        assert_eq!(hit["mismatches"], 0, "{name}");
        assert!(
            hit["source"]
                .as_str()
                .unwrap()
                .starts_with(accession.split('.').next().unwrap()),
            "{name}: {hit}"
        );
    }
}

/// Mutagenesis takes its template from a file, so a cloning experiment that is bound to no
/// document can run it like the other cloning commands.
#[test]
fn mutagenesis_from_a_file() {
    let mut s = common::session();
    let t = common::MUT_TEMPLATE;
    let fwd = format!("A{}", &t[41..60]);
    let rev: String = t[20..40]
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect();
    // No document is open: the template is the file.
    assert!(s.document().is_none());
    let r = run(
        &mut s,
        "mutagenesis",
        json!({"template": "mut_template.gb", "forward": fwd, "reverse": rev, "preview": true}),
    );
    assert_eq!(r["length_change"], 0);
    // Without a template and without a document, it says so instead of panicking.
    let e = s
        .execute("mutagenesis", json!({"forward": "ACGT", "reverse": "TTTT"}))
        .unwrap_err();
    assert_eq!(e.code, ErrorCode::NoDocument);
}

/// `file_info` describes a sequence file without opening it, and hands over the molecule when
/// asked, so a map can be drawn of a construct before it is chosen.
#[test]
fn file_info_reads_without_opening() {
    let mut s = common::session();
    let info = run(&mut s, "file_info", json!({"path": "example.dna"}));
    assert!(info["length"].as_u64().unwrap_or(0) > 0, "{info}");
    assert!(info["name"].as_str().is_some(), "{info}");
    assert!(matches!(
        info["topology"].as_str(),
        Some("circular" | "linear")
    ));
    // The molecule is left out unless it is asked for, and no document was opened.
    assert!(info.get("molecule").is_none(), "{info}");
    assert!(s.document().is_none());
    let full = run(
        &mut s,
        "file_info",
        json!({"path": "example.dna", "molecule": true}),
    );
    assert_eq!(full["molecule"]["name"], full["name"]);
    assert!(s.document().is_none());
    // A missing file is an error, not a panic.
    assert_eq!(
        s.execute("file_info", json!({"path": "nowhere.dna"}))
            .unwrap_err()
            .code,
        ErrorCode::Io
    );
    assert_eq!(
        s.execute("file_info", json!({"path": "  "}))
            .unwrap_err()
            .code,
        ErrorCode::InvalidParams
    );
}

/// Gateway refuses pairs that cannot recombine. The att cores are the same in attB, attL,
/// attR and attP, so without checking the arms an entry clone would happily "recombine" with
/// an expression clone and give a confident, wrong answer.
#[test]
fn gateway_refuses_sites_that_cannot_recombine() {
    const ATTL1: &str = "CAAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAAATTGATGAGCAATGCTTTTTTATAATGCCAACTTTGTACAAAAAAGCAGGCT";
    const ATTL2: &str = "ACCCAGCTTTCTTGTACAAAGTTGGCATTATAAAAAAGCATTGCTTATCAATTTGTTGCAACGAACAGGTCACTATCAGTCAAAATAAAATCATTATTTG";
    const ATTB1: &str = "ACAAGTTTGTACAAAAAAGCAGGCT";
    const ATTB2: &str = "ACCCAGCTTTCTTGTACAAAGTGGT";
    let cassette = "ATGAAACCCGGGTTTAAACCCGGGATGAAACCCGGG";
    let pad = "G".repeat(60);

    let mut s = common::session();
    // An entry clone, saved so it can be named as a file.
    run(
        &mut s,
        "new",
        json!({"name": "entry", "residues": format!("{pad}{ATTL1}{cassette}{ATTL2}{pad}"), "topology": "circular"}),
    );
    run(&mut s, "save", json!({"path": "entry_real.gb"}));
    // An expression clone as the open document: attB, the product of an LR, not a substrate.
    run(
        &mut s,
        "new",
        json!({"name": "expression", "residues": format!("{pad}{ATTB1}{cassette}{ATTB2}{pad}"), "topology": "circular", "discard_changes": true}),
    );
    let e = s
        .execute(
            "gateway_clone",
            json!({"donor": "entry_real.gb", "reaction": "lr", "preview": true}),
        )
        .unwrap_err();
    assert_eq!(e.code, ErrorCode::InvalidInput);
    assert!(
        e.message.contains("attB") && e.message.contains("attR"),
        "the message says which sites were found and which are needed: {}",
        e.message
    );
    assert!(
        e.hint.unwrap_or_default().contains("entry clone"),
        "the hint explains what LR joins"
    );
    // The same two the other way round is still no reaction, and says so plainly.
    let e = s
        .execute(
            "gateway_clone",
            json!({"acceptor": "entry_real.gb", "reaction": "lr", "preview": true}),
        )
        .unwrap_err();
    assert!(e.message.contains("attL"), "{}", e.message);
}
