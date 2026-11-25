//! Adding, changing and removing features.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_doc::{Feature, FeatureId, NewFeature};

use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
            CommandSpec::edit(
                "add_feature",
                "Add Feature",
                "Annotate a region (one or more segments). Returns the new feature's id.",
                add_feature,
            )
            .menu("Features/Add Feature…")
    .example(json!({"feature": {"name": "site", "kind": "misc_feature", "strand": "forward", "segments": [{"range": {"start": 1, "end": 3}}], "qualifiers": [{"name": "note", "values": ["example"]}]}})),
            CommandSpec::edit(
                "update_feature",
                "Edit Feature",
                "Replace a feature (matched by id) with new values.",
                update_feature,
            )
            .menu("Features/Edit Feature…")
    .example(json!({"feature": {"id": 1, "name": "renamed", "kind": "misc_feature", "strand": "reverse", "segments": [{"range": {"start": 10, "end": 2}}], "qualifiers": []}})),
            CommandSpec::edit(
                "remove_feature",
                "Delete Feature",
                "Remove a feature by id.",
                remove_feature,
            )
            .menu("Features/Delete Feature")
    .example(json!({"id": 1})),
        ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct AddParams {
    feature: NewFeature,
}

fn add_feature(s: &mut Session, p: AddParams) -> Result<Value, CmdError> {
    let id = s.molecule_mut()?.add_feature(p.feature)?;
    Ok(json!({ "id": id }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct UpdateParams {
    feature: Feature,
}

fn update_feature(s: &mut Session, p: UpdateParams) -> Result<Value, CmdError> {
    s.molecule_mut()?.update_feature(p.feature)?;
    Ok(json!({}))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct IdParams {
    id: FeatureId,
}

fn remove_feature(s: &mut Session, p: IdParams) -> Result<Value, CmdError> {
    let removed = s.molecule_mut()?.remove_feature(p.id)?;
    Ok(json!({ "removed": removed }))
}
