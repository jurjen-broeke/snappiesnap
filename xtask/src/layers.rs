//! Crate layering (AGENTS.md › Architecture rules).
//!
//! A crate may depend only on workspace crates in strictly lower layers, and nothing below
//! [`UI_LAYER`] may pull in a UI or windowing crate, directly or transitively. Every workspace
//! crate must be registered in [`LAYERS`].

use std::collections::{BTreeSet, HashMap, VecDeque};

use serde_json::Value;

use crate::util::{self, Result};

/// Layer used for developer tooling that sits outside the product graph.
const TOOL: u8 = u8::MAX;

/// Every workspace crate and its layer. Register new crates here.
pub const LAYERS: &[(&str, u8)] = &[
    ("snappiesnap-seq", 0),
    ("snappiesnap-formats", 0),
    ("snappiesnap-doc", 1),
    ("snappiesnap-bio", 2),
    ("snappiesnap-layout", 2),
    ("snappiesnap-render", 3),
    ("snappiesnap-io", 3),
    ("snappiesnap-engine", 4),
    ("snappiesnap-ui-egui", 5),
    ("snappiesnap-mcp", 5),
    ("snappiesnap", 6),
    ("snappiesnap-cli", 6),
    ("snappiesnap-web", 6),
    ("xtask", TOOL),
];

/// L0 crates that must not depend on anything in the workspace (publishable on their own).
const STANDALONE: &[&str] = &["snappiesnap-seq", "snappiesnap-formats"];

/// The first layer allowed to use UI crates.
pub const UI_LAYER: u8 = 5;

/// UI and windowing crates that are banned below [`UI_LAYER`].
const UI_CRATES: &[&str] = &["egui", "eframe", "winit", "wgpu", "rfd"];

pub fn layer_of(name: &str) -> Option<u8> {
    LAYERS.iter().find(|(n, _)| *n == name).map(|(_, l)| *l)
}

/// Product crates in layers L0..=L5 (the ones that must build for wasm).
pub fn core_crates() -> Vec<&'static str> {
    LAYERS
        .iter()
        .filter(|(_, l)| *l <= UI_LAYER)
        .map(|(n, _)| *n)
        .collect()
}

pub fn run() -> Result {
    let meta = util::metadata()?;
    util::verdict("layers", check(&meta))
}

fn check(meta: &Value) -> Vec<String> {
    let mut problems = Vec::new();
    let packages = util::packages(meta);
    let members: BTreeSet<&str> = util::workspace_members(meta).into_iter().collect();
    let names: HashMap<&str, &str> = packages.iter().map(|p| (p.id, p.name)).collect();

    for pkg in packages.iter().filter(|p| members.contains(p.id)) {
        let Some(layer) = layer_of(pkg.name) else {
            problems.push(format!(
                "`{}` is not registered in xtask/src/layers.rs",
                pkg.name
            ));
            continue;
        };
        if layer == TOOL {
            continue;
        }
        for dep in direct_workspace_deps(pkg.raw, &packages, &members) {
            match layer_of(dep) {
                Some(_) if STANDALONE.contains(&pkg.name) => problems.push(format!(
                    "`{}` is a standalone L0 crate but depends on `{dep}`",
                    pkg.name
                )),
                Some(dep_layer) if dep_layer >= layer => problems.push(format!(
                    "`{}` (L{layer}) depends on `{dep}` (L{dep_layer}); only lower layers are allowed",
                    pkg.name
                )),
                _ => {}
            }
        }
        if layer < UI_LAYER {
            for ui in reachable_ui_crates(meta, pkg.id, &names) {
                problems.push(format!(
                    "`{}` (L{layer}) pulls in UI crate `{ui}`; UI crates are allowed from L{UI_LAYER} only",
                    pkg.name
                ));
            }
        }
    }
    problems
}

/// Names of workspace crates this package depends on (normal and build dependencies).
fn direct_workspace_deps<'a>(
    raw: &'a Value,
    packages: &[util::Package<'a>],
    members: &BTreeSet<&str>,
) -> Vec<&'a str> {
    let member_names: BTreeSet<&str> = packages
        .iter()
        .filter(|p| members.contains(p.id))
        .map(|p| p.name)
        .collect();
    raw.get("dependencies")
        .and_then(Value::as_array)
        .map(|deps| {
            deps.iter()
                .filter(|d| d.get("kind").and_then(Value::as_str) != Some("dev"))
                .filter_map(|d| d.get("name").and_then(Value::as_str))
                .filter(|n| member_names.contains(n))
                .collect()
        })
        .unwrap_or_default()
}

/// UI crates reachable from `root_id` through non-dev edges of the resolve graph.
fn reachable_ui_crates<'a>(
    meta: &'a Value,
    root_id: &'a str,
    names: &HashMap<&str, &'a str>,
) -> BTreeSet<&'a str> {
    let mut found = BTreeSet::new();
    let Some(nodes) = meta
        .get("resolve")
        .and_then(|r| r.get("nodes"))
        .and_then(Value::as_array)
    else {
        return found;
    };
    let edges: HashMap<&str, Vec<&str>> = nodes
        .iter()
        .filter_map(|n| {
            let id = n.get("id")?.as_str()?;
            let deps = n
                .get("deps")?
                .as_array()?
                .iter()
                .filter(|d| !is_dev_only(d))
                .filter_map(|d| d.get("pkg").and_then(Value::as_str))
                .collect();
            Some((id, deps))
        })
        .collect();

    // Breadth-first with a seen-set: the graph is finite, so this always terminates.
    let mut seen = BTreeSet::from([root_id]);
    let mut queue = VecDeque::from([root_id]);
    while let Some(id) = queue.pop_front() {
        for &next in edges.get(id).map(Vec::as_slice).unwrap_or_default() {
            if !seen.insert(next) {
                continue;
            }
            if let Some(name) = names.get(next)
                && UI_CRATES.contains(name)
            {
                found.insert(*name);
            }
            queue.push_back(next);
        }
    }
    found
}

fn is_dev_only(dep: &Value) -> bool {
    dep.get("dep_kinds")
        .and_then(Value::as_array)
        .is_some_and(|kinds| {
            !kinds.is_empty()
                && kinds
                    .iter()
                    .all(|k| k.get("kind").and_then(Value::as_str) == Some("dev"))
        })
}
