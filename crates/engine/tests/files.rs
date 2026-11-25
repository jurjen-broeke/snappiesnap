//! Opening files through the platform layer (in memory, never the real file system).

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(clippy::unwrap_used, clippy::indexing_slicing)]

mod common;

use std::sync::Arc;

use serde_json::json;
use snappiesnap_engine::doc::{Range, Strand, Topology};
use snappiesnap_engine::{ErrorCode, MemoryFs, Session};

#[test]
fn open_a_dna_file() {
    let mut s = common::session();
    let out = s.execute("open", json!({"path": "example.dna"})).unwrap();
    assert_eq!(out["name"], "example");
    assert_eq!(out["length"], 16);
    assert_eq!(out["features"], 1);
    assert_eq!(out["warnings"], json!([]));
    let doc = s.document().unwrap();
    assert!(!doc.dirty && doc.source.is_some());
    assert_eq!(doc.path.as_deref(), Some("example.dna"));
    let m = &doc.molecule;
    assert_eq!(m.topology(), Topology::Circular);
    let f = &m.features()[0];
    assert_eq!(f.segments[0].range, Range::new(14, 2));
    assert_eq!(f.strand, Strand::Forward);
    assert_eq!(m.residues(f.segments[0].range).unwrap(), b"CCAC");
}

#[test]
fn open_failures_are_clear_and_change_nothing() {
    let fs = MemoryFs::new()
        .with_file("junk.dna", b"not a dna file".to_vec())
        .with_file("cut.dna", common::example_dna()[..30].to_vec());
    let mut s = Session::new().with_platform(Arc::new(fs));
    s.execute("new", json!({"residues": "ACGT"})).unwrap();
    let before = s.document().cloned();
    for (path, code) in [
        ("missing.dna", ErrorCode::Io),
        ("junk.dna", ErrorCode::InvalidInput),
        ("cut.dna", ErrorCode::InvalidInput),
    ] {
        let e = s.execute("open", json!({"path": path})).unwrap_err();
        assert_eq!(e.code, code, "{path}: {e}");
    }
    assert_eq!(s.document().cloned(), before);
}

#[test]
fn save_round_trips_and_keeps_edits() {
    let fs = Arc::new(MemoryFs::new().with_file("example.dna", common::example_dna()));
    let mut s = Session::new().with_platform(fs.clone());
    s.execute("open", json!({"path": "example.dna"})).unwrap();

    // Unchanged: identical bytes.
    s.execute("save", json!({"path": "same.dna"})).unwrap();
    assert_eq!(fs.get("same.dna").unwrap(), common::example_dna());
    assert_eq!(s.document().unwrap().path.as_deref(), Some("same.dna"));

    // Edit, save to the document's own path, reopen.
    s.execute("insert", json!({"at": 0, "residues": "GGG"}))
        .unwrap();
    s.execute(
        "add_feature",
        json!({"feature": {"name": "new", "segments": [{"range": {"start": 0, "end": 3}}]}}),
    )
    .unwrap();
    assert!(s.document().unwrap().dirty);
    let out = s.execute("save", json!({})).unwrap();
    assert_eq!(out["path"], "same.dna");
    assert!(!s.document().unwrap().dirty);
    s.execute("open", json!({"path": "same.dna"})).unwrap();
    let m = &s.document().unwrap().molecule;
    assert_eq!(m.len(), 19);
    assert_eq!(m.features().len(), 2);
    assert_eq!(m.residues(Range::new(0, 3)).unwrap(), b"GGG");
}

#[test]
fn save_failures_are_clear_and_change_nothing() {
    let mut s = common::session();
    s.execute("new", json!({"residues": "ACGT"})).unwrap();
    let e = s.execute("save", json!({})).unwrap_err();
    assert_eq!(e.code, ErrorCode::InvalidParams);
    let e = s.execute("save", json!({"path": "x.txt"})).unwrap_err();
    assert_eq!(e.code, ErrorCode::Unsupported);
    // A protein can't be a .dna file: refused, and the document is as it was.
    s.execute("new", json!({"residues": "MKV", "alphabet": "protein"}))
        .unwrap();
    let before = s.document().cloned();
    let e = s.execute("save", json!({"path": "x.dna"})).unwrap_err();
    assert_eq!(e.code, ErrorCode::Unsupported, "{e}");
    assert_eq!(s.document().cloned(), before);
}

#[test]
fn file_management_never_overwrites_and_follows_the_document() {
    let fs = Arc::new(
        MemoryFs::new()
            .with_file("lab/a.dna", common::example_dna())
            .with_file("lab/b.dna", common::example_dna()),
    );
    let mut s = Session::new().with_platform(fs.clone());
    s.execute("open", json!({"path": "lab/a.dna"})).unwrap();

    let list = s.execute("list_folder", json!({"path": "lab"})).unwrap();
    assert_eq!(list["entries"].as_array().unwrap().len(), 2);

    // Rename the open document: it follows.
    s.execute(
        "file_rename",
        json!({"from": "lab/a.dna", "to": "lab/c.dna"}),
    )
    .unwrap();
    assert_eq!(s.document().unwrap().path.as_deref(), Some("lab/c.dna"));
    // Move its folder: it still follows.
    s.execute("file_rename", json!({"from": "lab", "to": "lab2"}))
        .unwrap();
    assert_eq!(s.document().unwrap().path.as_deref(), Some("lab2/c.dna"));

    // Never overwrite.
    let e = s
        .execute(
            "file_rename",
            json!({"from": "lab2/c.dna", "to": "lab2/b.dna"}),
        )
        .unwrap_err();
    assert_eq!(e.code, ErrorCode::InvalidInput);
    let e = s
        .execute(
            "file_copy",
            json!({"from": "lab2/c.dna", "to": "lab2/b.dna"}),
        )
        .unwrap_err();
    assert_eq!(e.code, ErrorCode::InvalidInput);
    assert_eq!(fs.get("lab2/b.dna").unwrap(), common::example_dna());

    s.execute(
        "file_copy",
        json!({"from": "lab2/c.dna", "to": "lab2/c copy.dna"}),
    )
    .unwrap();
    s.execute("folder_create", json!({"path": "lab2/new"}))
        .unwrap();
    let names: Vec<String> = s.execute("list_folder", json!({"path": "lab2"})).unwrap()["entries"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["name"].as_str().unwrap().to_string())
        .collect();
    assert_eq!(names, vec!["new", "b.dna", "c copy.dna", "c.dna"]);
    assert_eq!(
        s.execute("list_folder", json!({"path": "missing"}))
            .unwrap_err()
            .code,
        ErrorCode::Io
    );
}
