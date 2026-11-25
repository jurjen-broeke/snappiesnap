//! Never crash (AGENTS.md): every registered command, in every kind of session state, with
//! adversarial params. Commands are discovered from the registry, so new commands are covered
//! automatically.
//!
//! For every call: no panic (no `internal` error), the document stays valid, and a failed call
//! leaves the document exactly as it was.

// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::indexing_slicing,
    clippy::panic
)]

mod common;

use serde_json::{Value, json};
use snappiesnap_engine::doc::{Alphabet, Molecule, NewFeature, Range, Segment, Seq, Topology};
use snappiesnap_engine::{Document, ErrorCode, Session};

/// Session states every command is tried in.
fn states() -> Vec<(&'static str, Session)> {
    let with = |m: Molecule| {
        let mut s = common::session();
        s.open_document(Document::new(m));
        s
    };
    let mol = |res: &str, alphabet, topo| {
        Molecule::new("m", Seq::new(alphabet, res).unwrap(), topo).unwrap()
    };
    let mut plasmid = mol("AAAACCCCGGGG", Alphabet::Dna, Topology::Circular);
    plasmid
        .add_feature(NewFeature {
            color: None,
            name: "span".into(),
            kind: "misc_feature".into(),
            strand: Default::default(),
            segments: vec![Segment {
                range: Range::new(10, 2),
            }],
            qualifiers: vec![],
        })
        .unwrap();
    plasmid
        .add_primer("p", Seq::new(Alphabet::Dna, "ACGT").unwrap(), "")
        .unwrap();
    vec![
        ("no document", common::session()),
        (
            "empty linear",
            with(mol("", Alphabet::Dna, Topology::Linear)),
        ),
        ("1 bp", with(mol("A", Alphabet::Dna, Topology::Linear))),
        ("plasmid", {
            // With history, so undo and redo are both possible.
            let mut s = with(plasmid);
            s.execute("insert", json!({"at": 0, "residues": "A"}))
                .unwrap();
            s.execute("insert", json!({"at": 0, "residues": "C"}))
                .unwrap();
            s.execute("undo", json!({})).unwrap();
            s
        }),
        (
            "protein",
            with(mol("MKV*", Alphabet::Protein, Topology::Linear)),
        ),
        ("with history", {
            let mut s = common::session();
            s.execute("open", json!({"path": "history.dna"})).unwrap();
            s
        }),
        ("matched reads", {
            let mut s = common::session();
            s.execute("open", json!({"path": "matched.dna"})).unwrap();
            s
        }),
        ("read", {
            // A sequencing read with its trace, reverse-complemented so mapping is exercised.
            let mut s = common::session();
            s.execute("open", json!({"path": "read.ab1"})).unwrap();
            s.execute("reverse_complement", json!({})).unwrap();
            s
        }),
        ("template", {
            // A made-up template with a primer pair that makes a 63 bp product.
            let mut m = mol(
                "GATTACAGATTACACCGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCG",
                Alphabet::Dna,
                Topology::Circular,
            );
            m.add_primer(
                "f",
                Seq::new(Alphabet::Dna, "CACCGTTGCATGCAAGCTTG").unwrap(),
                "",
            )
            .unwrap();
            m.add_primer(
                "r",
                Seq::new(Alphabet::Dna, "CAGGGTTTTCCCAGTCACGA").unwrap(),
                "",
            )
            .unwrap();
            with(m)
        }),
        ("two tabs", {
            // Two documents open, for commands that take another tab by its id.
            let mut s = with(mol("ACGTACGTAAAACCCCTTTT", Alphabet::Dna, Topology::Linear));
            s.open_document(Document::new(mol(
                "ACGTACGTAAAACCCCTT",
                Alphabet::Dna,
                Topology::Linear,
            )));
            s
        }),
    ]
}

/// Values that have broken parsers before: wrong types, extremes, odd text, deep nesting.
fn hostile_values() -> Vec<Value> {
    let mut deep = json!(0);
    for _ in 0..100 {
        deep = json!([deep]);
    }
    vec![
        Value::Null,
        json!(true),
        json!(-1),
        json!(0),
        json!(1),
        json!(12),
        json!(13),
        json!(u64::MAX),
        json!(i64::MIN),
        json!(1e308),
        json!(0.5),
        json!(""),
        json!(" "),
        json!("\u{0}"),
        json!("💥🧬"),
        json!("ACGT"),
        json!("acgtNNNN"),
        json!("A".repeat(100_000)),
        json!("linear"),
        json!("circular"),
        json!([]),
        json!([null]),
        json!({}),
        json!({"start": u64::MAX, "end": 0}),
        json!({"start": 10, "end": 2}),
        json!({"start": 2, "end": 2}),
        json!({"start": 0, "end": 1e20}),
        json!({"start": 0, "end": 4, "extra": 1}),
        json!({"name": "x", "segments": [{"range": {"start": 0, "end": 1}}]}),
        json!({"name": "x", "segments": []}),
        json!({"id": 1, "name": "x", "kind": "", "strand": "forward", "segments": [{"range": {"start": 11, "end": 1}}]}),
        deep,
    ]
}

/// Params to try for one command: the hostile values as whole params, then each example with
/// one value at a time (at any depth) replaced by every hostile value or removed. Everything
/// else stays valid, so the command's logic is reached, not just its parser.
fn params_for(examples: &[Value]) -> Vec<Value> {
    let mut out = hostile_values();
    for example in examples {
        out.push(example.clone());
        for path in paths(example, &mut Vec::new()) {
            for v in hostile_values() {
                out.push(replaced(example, &path, Some(v)));
            }
            out.push(replaced(example, &path, None));
        }
    }
    out
}

#[derive(Clone)]
enum Step {
    Key(String),
    Index(usize),
}

/// Every path to a value inside `v` (not including the root), depth-limited.
fn paths(v: &Value, prefix: &mut Vec<Step>) -> Vec<Vec<Step>> {
    let mut out = Vec::new();
    if prefix.len() >= 6 {
        return out;
    }
    let children: Vec<(Step, &Value)> = match v {
        Value::Object(m) => m.iter().map(|(k, v)| (Step::Key(k.clone()), v)).collect(),
        Value::Array(a) => a
            .iter()
            .enumerate()
            .map(|(i, v)| (Step::Index(i), v))
            .collect(),
        _ => Vec::new(),
    };
    for (step, child) in children {
        prefix.push(step);
        out.push(prefix.clone());
        out.extend(paths(child, prefix));
        prefix.pop();
    }
    out
}

/// A copy of `v` with the value at `path` replaced (or removed when `with` is `None`).
fn replaced(v: &Value, path: &[Step], with: Option<Value>) -> Value {
    let mut out = v.clone();
    let Some((last, parents)) = path.split_last() else {
        return out;
    };
    let mut cur = &mut out;
    for step in parents {
        cur = match step {
            Step::Key(k) => &mut cur[k.as_str()],
            Step::Index(i) => &mut cur[*i],
        };
    }
    match (last, with) {
        (Step::Key(k), Some(v)) => cur[k.as_str()] = v,
        (Step::Index(i), Some(v)) => cur[*i] = v,
        (Step::Key(k), None) => {
            if let Some(m) = cur.as_object_mut() {
                m.remove(k);
            }
        }
        (Step::Index(i), None) => {
            if let Some(a) = cur.as_array_mut() {
                a.remove(*i);
            }
        }
    }
    out
}

#[test]
fn every_command_survives_hostile_params_in_every_state() {
    let ids: Vec<(String, Vec<Value>)> = Session::new()
        .registry()
        .iter()
        .map(|c| (c.id.to_string(), c.examples.clone()))
        .collect();
    assert!(ids.len() >= 15);
    let fresh = states();
    let mut calls = 0usize;

    // Pass 1: every call starts from a fresh copy of its state, so no earlier call (such as
    // `close`) can hide later commands from their intended state.
    for (state_name, state) in &fresh {
        for (id, examples) in &ids {
            for params in params_for(examples) {
                let mut s = state.clone();
                check_call(&mut s, state_name, id, &params);
                calls += 1;
            }
        }
    }

    // Pass 2: the same calls one after another in each state, so damage could accumulate.
    for (state_name, mut s) in fresh {
        for (id, examples) in &ids {
            for params in params_for(examples) {
                check_call(&mut s, state_name, id, &params);
                calls += 1;
            }
        }
        // Whatever happened, undo/redo still work.
        while s.can_undo() {
            s.execute("undo", json!({})).unwrap();
        }
        while s.can_redo() {
            s.execute("redo", json!({})).unwrap();
        }
    }
    assert!(calls > 20_000, "only {calls} calls");
}

/// One call: no panic, a failure changes nothing, and the document stays valid.
fn check_call(s: &mut Session, state_name: &str, id: &str, params: &Value) {
    let before = s.document().cloned();
    if let Err(e) = s.execute(id, params.clone()) {
        assert_ne!(
            e.code,
            ErrorCode::Internal,
            "{id} panicked in state {state_name:?} with {params}: {e}"
        );
        assert_eq!(
            s.document(),
            before.as_ref(),
            "failed {id} changed the document ({state_name}, {params})"
        );
    }
    if let Some(doc) = s.document() {
        doc.molecule.validate().unwrap_or_else(|e| {
            panic!("{id} left an invalid document ({state_name}, {params}): {e}")
        });
    }
}

/// Every command documents at least one example, and every example works in at least one of:
/// no document, a 1 bp molecule, the plasmid (which has unsaved edits, undo/redo history and a
/// primer), the template with a working primer pair, the history file, the read, or the two
/// open tabs.
#[test]
fn every_command_has_working_examples() {
    let ids: Vec<(String, Vec<Value>)> = Session::new()
        .registry()
        .iter()
        .map(|c| (c.id.to_string(), c.examples.clone()))
        .collect();
    for (id, examples) in ids {
        assert!(!examples.is_empty(), "command {id} has no .example(...)");
        for example in examples {
            let tried: Vec<(&'static str, Session)> = states()
                .into_iter()
                .filter(|(n, _)| {
                    matches!(
                        *n,
                        "no document"
                            | "1 bp"
                            | "plasmid"
                            | "template"
                            | "with history"
                            | "read"
                            | "matched reads"
                            | "two tabs"
                    )
                })
                .collect();
            let count = tried.len();
            let errors: Vec<String> = tried
                .into_iter()
                .filter_map(|(n, mut s)| {
                    s.execute(&id, example.clone())
                        .err()
                        .map(|e| format!("{n}: {e}"))
                })
                .collect();
            assert!(
                errors.len() < count,
                "example for {id} fails everywhere: {example} -> {errors:?}"
            );
        }
    }
}
