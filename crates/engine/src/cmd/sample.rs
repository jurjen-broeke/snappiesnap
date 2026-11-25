//! Open the synthetic sample constructs (see [`crate::sample`]): a demo plasmid, a fragment
//! that clones into it, or both ready for a cloning experiment.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "sample_open",
            "Open Sample Constructs",
            "Open SnappieSnap's synthetic sample constructs, for a demo, a screenshot or a \
             first experiment. `which`: \"plasmid\" (default) is pSample, a 3.2 kb circular \
             plasmid with made-up features and a polylinker; \"insert\" is demoInsert, a linear \
             fragment with EcoRI and BamHI ends; \"cloning\" opens both, which cut with EcoRI \
             and BamHI give one product. The sequences are generated, not taken from any \
             database, and are the same on every platform.",
            sample_open,
        )
        .needs_document(false)
        .example(json!({}))
        .example(json!({ "which": "cloning" }))
        .example(json!({ "which": "insert", "discard_changes": true })),
    ]
}

#[derive(Deserialize, JsonSchema, Default, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum Which {
    #[default]
    Plasmid,
    Insert,
    /// Both, in the order a cloning experiment takes them: the vector, then the insert.
    Cloning,
}

#[derive(Deserialize, JsonSchema, Default)]
#[serde(default, deny_unknown_fields)]
struct Params {
    /// Which sample to open (default the plasmid).
    which: Which,
    /// Replace unsaved work in the active document instead of refusing. A sample is a demo,
    /// so it never takes the user's work away without being told to.
    discard_changes: bool,
}

fn sample_open(s: &mut Session, p: Params) -> Result<Value, CmdError> {
    // A sample is a demo, not the user's work, so it is never dirty on arrival; `new` already
    // refuses to discard unsaved work unless asked.
    if p.which != Which::Insert {
        crate::sample::load(s, p.discard_changes)?;
        clean(s);
    }
    if p.which != Which::Plasmid {
        // The second one opens in a tab of its own, so it discards nothing either way.
        crate::sample::load_insert(s, p.discard_changes)?;
        clean(s);
    }
    let documents: Vec<Value> = s
        .documents()
        .iter()
        .map(|d| json!({"id": d.id, "name": d.name}))
        .collect();
    let opened = match p.which {
        Which::Plasmid => json!(["pSample"]),
        Which::Insert => json!(["demoInsert"]),
        Which::Cloning => json!(["pSample", "demoInsert"]),
    };
    Ok(json!({ "opened": opened, "documents": documents }))
}

/// A freshly opened sample has nothing worth saving, so closing the page or the window needn't
/// ask about it.
fn clean(s: &mut Session) {
    if let Some(doc) = s.document_mut() {
        doc.dirty = false;
    }
}
