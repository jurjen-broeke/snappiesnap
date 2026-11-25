//! Construction history: how a molecule was made, as a tree of steps.
//!
//! The root step is the molecule itself; its inputs are the molecules it was made from, each
//! with the steps that made it. A step keeps a copy of its molecule (without history: the
//! steps below carry that), so any ancestor can be opened again.

use serde::{Deserialize, Serialize};
use snappiesnap_seq::Topology;

use crate::{DocError, Molecule};

/// Deepest history we keep; deeper ancestry is cut off (it stays in the files it came from).
pub const MAX_DEPTH: usize = 32;
/// Most steps in one history.
pub const MAX_NODES: usize = 10_000;

/// One step of a construction history.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HistoryNode {
    /// Unique within the history.
    pub id: u32,
    pub name: String,
    pub length: usize,
    pub topology: Topology,
    /// What made this molecule, as an id (`insertFragment`, `amplifyFragment`, ...); empty for
    /// a starting molecule. See [`operation_label`].
    #[serde(default)]
    pub operation: String,
    /// The particulars, in words (enzymes, primers, ...).
    #[serde(default)]
    pub details: String,
    /// The molecule as it was at this step. `None` for the root (the document itself) and for
    /// steps whose sequence wasn't kept.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub molecule: Option<Box<Molecule>>,
    /// What it was made from.
    #[serde(default)]
    pub inputs: Vec<HistoryNode>,
    /// What each input contributed, one per input, as attribute name/value pairs in the
    /// `.dna` style (`manipulation` first: `amplify`, `insert`, `replace`, ...).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub summaries: Vec<Vec<(String, String)>>,
    /// Primers the step used (name, sequence), as for a PCR.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub oligos: Vec<(String, String)>,
    /// Settings of the step (name, value), such as the polymerase of a PCR.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub parameters: Vec<(String, String)>,
}

/// Most attributes in one input summary, and longest value.
const MAX_SUMMARY_ATTRS: usize = 64;
const MAX_SUMMARY_TEXT: usize = 10_000;

/// A step for display: no molecules, but whether each can be opened.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct HistoryStep {
    pub id: u32,
    pub name: String,
    /// The operation in words ("" for a starting molecule).
    pub operation: String,
    /// The operation id as stored.
    pub operation_id: String,
    pub details: String,
    pub length: usize,
    pub circular: bool,
    /// Whether the step's molecule is kept (and can be opened).
    pub available: bool,
    pub inputs: Vec<HistoryStep>,
}

impl HistoryNode {
    /// A molecule as an input to a new step: its own history (if any) becomes the step, and the
    /// step keeps a copy of it.
    pub fn input(m: &Molecule) -> HistoryNode {
        let mut snapshot = m.clone();
        let made = snapshot.take_history();
        let made = made.unwrap_or_else(|| HistoryNode {
            id: 0,
            name: String::new(),
            length: 0,
            topology: m.topology(),
            operation: String::new(),
            details: String::new(),
            molecule: None,
            inputs: Vec::new(),
            summaries: Vec::new(),
            oligos: Vec::new(),
            parameters: Vec::new(),
        });
        HistoryNode {
            id: 0,
            name: m.name().to_string(),
            length: m.len(),
            topology: m.topology(),
            operation: made.operation,
            details: made.details,
            molecule: Some(Box::new(snapshot)),
            inputs: made.inputs,
            summaries: made.summaries,
            oligos: made.oligos,
            parameters: made.parameters,
        }
    }

    /// The step that made `product` from `inputs` (ids renumbered, ancestry deeper than
    /// [`MAX_DEPTH`] or beyond [`MAX_NODES`] steps cut off).
    pub fn made(
        product: &Molecule,
        operation: &str,
        details: &str,
        inputs: Vec<HistoryNode>,
    ) -> HistoryNode {
        let mut root = HistoryNode {
            id: 0,
            name: product.name().to_string(),
            length: product.len(),
            topology: product.topology(),
            operation: operation.to_string(),
            details: details.chars().take(2000).collect(),
            molecule: None,
            inputs,
            summaries: Vec::new(),
            oligos: Vec::new(),
            parameters: Vec::new(),
        };
        root.prune();
        root.renumber();
        root
    }

    /// Visit every step (depth first, parents before their inputs).
    pub fn walk(&self, f: &mut dyn FnMut(&HistoryNode)) {
        let mut stack = vec![self];
        while let Some(n) = stack.pop() {
            f(n);
            stack.extend(n.inputs.iter().rev());
        }
    }

    /// The step with this id.
    pub fn find(&self, id: u32) -> Option<&HistoryNode> {
        let mut stack = vec![self];
        while let Some(n) = stack.pop() {
            if n.id == id {
                return Some(n);
            }
            stack.extend(n.inputs.iter());
        }
        None
    }

    /// The tree for display.
    pub fn steps(&self) -> HistoryStep {
        fn convert(n: &HistoryNode, root: bool, depth: usize) -> HistoryStep {
            HistoryStep {
                id: n.id,
                name: n.name.clone(),
                operation: operation_label(&n.operation),
                operation_id: n.operation.clone(),
                details: n.details.clone(),
                length: n.length,
                circular: n.topology == Topology::Circular,
                available: root || n.molecule.is_some(),
                inputs: if depth < MAX_DEPTH {
                    n.inputs
                        .iter()
                        .map(|i| convert(i, false, depth + 1))
                        .collect()
                } else {
                    Vec::new()
                },
            }
        }
        convert(self, true, 0)
    }

    /// Check the limits: depth, number of steps, unique ids, copies without history.
    pub fn validate(&self) -> Result<(), DocError> {
        let mut ids = std::collections::BTreeSet::new();
        let mut stack = vec![(self, 0usize)];
        while let Some((n, depth)) = stack.pop() {
            if depth > MAX_DEPTH {
                return Err(DocError::Inconsistent(format!(
                    "history deeper than {MAX_DEPTH} steps"
                )));
            }
            if !ids.insert(n.id) {
                return Err(DocError::Inconsistent(format!(
                    "history step id {} is duplicated",
                    n.id
                )));
            }
            if ids.len() > MAX_NODES {
                return Err(DocError::Inconsistent(format!(
                    "history has more than {MAX_NODES} steps"
                )));
            }
            if n.summaries.len() > n.inputs.len().max(1).saturating_mul(4)
                || n.summaries.iter().any(|s| {
                    s.len() > MAX_SUMMARY_ATTRS
                        || s.iter()
                            .any(|(k, v)| k.len() > MAX_SUMMARY_TEXT || v.len() > MAX_SUMMARY_TEXT)
                })
            {
                return Err(DocError::Inconsistent(format!(
                    "history step {} has too many or too long input summaries",
                    n.id
                )));
            }
            if n.molecule.as_ref().is_some_and(|m| m.history().is_some()) {
                return Err(DocError::Inconsistent(format!(
                    "history step {} keeps a molecule with a history of its own",
                    n.id
                )));
            }
            stack.extend(n.inputs.iter().map(|i| (i, depth + 1)));
        }
        Ok(())
    }

    /// Cut off steps deeper than [`MAX_DEPTH`] and beyond [`MAX_NODES`] (breadth first, so
    /// the nearest ancestors stay).
    pub fn prune(&mut self) {
        // Steps kept so far: those visited plus those waiting.
        let mut total = 1usize;
        let mut level: Vec<&mut HistoryNode> = vec![self];
        let mut depth = 0usize;
        while !level.is_empty() {
            let mut next = Vec::new();
            for n in level {
                if depth >= MAX_DEPTH {
                    n.inputs.clear();
                }
                n.inputs.truncate(MAX_NODES.saturating_sub(total));
                total += n.inputs.len();
                next.extend(n.inputs.iter_mut());
            }
            level = next;
            depth += 1;
        }
    }

    /// Number the steps 0, 1, 2, ... depth first.
    pub fn renumber(&mut self) {
        let mut id = 0u32;
        let mut stack: Vec<&mut HistoryNode> = vec![self];
        while let Some(n) = stack.pop() {
            n.id = id;
            id = id.saturating_add(1);
            stack.extend(n.inputs.iter_mut().rev());
        }
    }
}

/// The operation ids in words.
pub fn operation_label(op: &str) -> String {
    let fixed = match op {
        "" => "",
        "amplifyFragment" => "PCR",
        "changeMethylation" => "Methylation changed",
        "gatewayLRCloning" => "Gateway LR reaction",
        "gatewayBPCloning" => "Gateway BP reaction",
        "gibsonAssembly" => "Gibson assembly",
        "hifiAssembly" => "HiFi DNA assembly",
        "inFusionCloning" => "In-Fusion cloning",
        "goldenGateAssembly" => "Golden Gate assembly",
        "insert" => "Insertion",
        "remove" => "Deletion",
        "replace" => "Replacement",
        "insertFragment" => "Restriction cloning",
        "removeRestrictionFragment" => "Restriction fragment removed",
        "primerDirectedMutagenesis" => "Site-directed mutagenesis",
        "setOrigin" => "Origin changed",
        "reverseComplement" => "Reverse complemented",
        "vnti" => "Vector NTI operation",
        "taCloning" => "TA cloning",
        "topoCloning" => "TOPO cloning",
        other => {
            // camelCase in words.
            let mut s = String::new();
            for (i, c) in other.chars().take(200).enumerate() {
                if c.is_uppercase() && i > 0 {
                    s.push(' ');
                    s.extend(c.to_lowercase());
                } else if i == 0 {
                    s.extend(c.to_uppercase());
                } else {
                    s.push(c);
                }
            }
            return s;
        }
    };
    fixed.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Seq;

    fn mol(name: &str, residues: &str, topology: Topology) -> Molecule {
        Molecule::new(
            name,
            Seq::new(crate::Alphabet::Dna, residues).unwrap(),
            topology,
        )
        .unwrap()
    }

    #[test]
    fn products_carry_the_history_of_their_inputs() {
        let a = mol("a", "ACGT", Topology::Linear);
        let b = mol("b", "GGCC", Topology::Circular);
        let mut ab = mol("ab", "ACGTGGCC", Topology::Circular);
        let h = HistoryNode::made(
            &ab,
            "gibsonAssembly",
            "",
            vec![HistoryNode::input(&a), HistoryNode::input(&b)],
        );
        ab.set_history(Some(h)).unwrap();
        let mut abc = mol("abc", "ACGTGGCCTT", Topology::Linear);
        let h = HistoryNode::made(
            &abc,
            "amplifyFragment",
            "F/R",
            vec![HistoryNode::input(&ab)],
        );
        abc.set_history(Some(h)).unwrap();

        let root = abc.history().unwrap();
        assert_eq!(root.name, "abc");
        assert!(root.molecule.is_none());
        let [first] = root.inputs.as_slice() else {
            panic!("one input")
        };
        assert_eq!(first.operation, "gibsonAssembly");
        assert_eq!(first.length, 8);
        assert_eq!(first.molecule.as_ref().unwrap().history(), None);
        assert_eq!(first.inputs.len(), 2);
        let ids: Vec<u32> = {
            let mut v = Vec::new();
            root.walk(&mut |n| v.push(n.id));
            v
        };
        assert_eq!(ids, vec![0, 1, 2, 3]);
        assert_eq!(root.find(3).unwrap().name, "b");
        let steps = root.steps();
        assert_eq!(steps.operation, "PCR");
        assert!(steps.inputs[0].inputs[1].circular);
        assert!(steps.inputs[0].inputs[1].available);

        // The native form keeps it all.
        let json = serde_json::to_string(&abc).unwrap();
        let back: Molecule = serde_json::from_str(&json).unwrap();
        assert_eq!(back, abc);
        // Older files without history still load.
        let plain = serde_json::to_string(&a).unwrap();
        assert!(!plain.contains("history"));
    }

    #[test]
    fn deep_ancestry_is_cut_off() {
        let mut m = mol("m", "A", Topology::Linear);
        for i in 0..(MAX_DEPTH + 10) {
            let mut next = mol(&format!("m{i}"), "A", Topology::Linear);
            let h = HistoryNode::made(&next, "insert", "", vec![HistoryNode::input(&m)]);
            next.set_history(Some(h)).unwrap();
            m = next;
        }
        let mut depth = 0;
        let mut n = m.history().unwrap();
        while let Some(i) = n.inputs.first() {
            n = i;
            depth += 1;
        }
        assert_eq!(depth, MAX_DEPTH);
        m.history().unwrap().validate().unwrap();
    }

    #[test]
    fn invalid_histories_are_refused() {
        let a = mol("a", "ACGT", Topology::Linear);
        let mut b = mol("b", "ACGT", Topology::Linear);
        let mut h = HistoryNode::made(&b, "x", "", vec![HistoryNode::input(&a)]);
        if let Some(i) = h.inputs.first_mut() {
            i.id = 0;
        }
        assert!(b.set_history(Some(h)).is_err());
        let mut with = a.clone();
        let made = HistoryNode::made(&a, "x", "", vec![]);
        with.set_history(Some(made)).unwrap();
        let mut node = HistoryNode::input(&a);
        node.molecule = Some(Box::new(with));
        let h = HistoryNode::made(&b, "x", "", vec![node]);
        assert!(b.set_history(Some(h)).is_err());
    }

    #[test]
    fn operations_read_as_words() {
        assert_eq!(operation_label("insertFragment"), "Restriction cloning");
        assert_eq!(operation_label("someNewThing"), "Some new thing");
        assert_eq!(operation_label(""), "");
    }
}
