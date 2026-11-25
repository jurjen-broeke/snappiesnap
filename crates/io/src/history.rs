//! Construction history of `.dna` files: read into the document model ([`HistoryNode`]) and
//! written back in the layout of the reference application's files.

use snappiesnap_doc::{HistoryNode, MAX_HISTORY_DEPTH, Molecule, Topology};
use snappiesnap_formats::dna::history::{
    self as dna_history, History, NodeOut, encode_sequence, features_element, node_record, tree_xml,
};
use snappiesnap_formats::dna::{DnaFile, Packet, tag};

use crate::IoError;

/// Most steps read from one file.
const MAX_STEPS: usize = 10_000;

/// The construction history of a `.dna` file, with a copy of every step whose sequence the file
/// stores or lets us rebuild (`Ok(None)` if it has no history). The root is the file itself.
pub(crate) fn read(file: &DnaFile) -> Result<Option<HistoryNode>, IoError> {
    let Some(h) = file.history()? else {
        return Ok(None);
    };
    let root = file.sequence()?.residues;
    let mut count = 0usize;
    let mut node = convert(file, &h, &h.tree, &root, 0, &mut count);
    node.molecule = None;
    if node.validate().is_err() {
        // Ids that aren't unique (or similar): number the steps afresh.
        node.prune();
        node.renumber();
        node.validate()?;
    }
    Ok(Some(node))
}

fn convert(
    file: &DnaFile,
    h: &History,
    n: &dna_history::HistoryNode,
    root: &[u8],
    depth: usize,
    count: &mut usize,
) -> HistoryNode {
    *count += 1;
    let inputs = if depth < MAX_HISTORY_DEPTH {
        let mut v = Vec::new();
        for i in &n.inputs {
            if *count >= MAX_STEPS {
                break;
            }
            v.push(convert(file, h, i, root, depth + 1, count));
        }
        v
    } else {
        Vec::new()
    };
    let name = clean_name(&n.name, n.id);
    let molecule = if depth == 0 {
        None
    } else {
        step_molecule(file, h, n, root, &name).ok().map(Box::new)
    };
    HistoryNode {
        id: n.id,
        name,
        length: n.length,
        topology: if n.circular {
            Topology::Circular
        } else {
            Topology::Linear
        },
        // "invalid" marks a starting molecule (nothing made it).
        operation: if n.operation == "invalid" {
            String::new()
        } else {
            n.operation.clone()
        },
        details: String::new(),
        molecule,
        inputs,
        summaries: n.summaries.clone(),
        oligos: n.oligos.clone(),
        parameters: n.parameters.clone(),
    }
}

fn clean_name(name: &str, id: u32) -> String {
    let name = name.strip_suffix(".dna").unwrap_or(name);
    let name: String = name.chars().filter(|c| !c.is_control()).take(200).collect();
    if name.trim().is_empty() {
        format!("step {id}")
    } else {
        name
    }
}

/// A step as a molecule: its residues (stored or rebuilt), topology, and the features, primers
/// and notes stored with it.
fn step_molecule(
    file: &DnaFile,
    h: &History,
    node: &dna_history::HistoryNode,
    root: &[u8],
    name: &str,
) -> Result<Molecule, IoError> {
    let residues = h
        .residues(node.id, root)
        .ok_or_else(|| IoError::Format(format!("step {} isn't stored", node.id)))?;
    // A small .dna file of its own: the main file's cookie, the step's sequence and its
    // stored packets.
    let flags = file.sequence()?.flags;
    let flags = if node.circular {
        flags | 0x01
    } else {
        flags & !0x01
    };
    let mut seq = vec![flags];
    seq.extend_from_slice(&residues);
    let mut packets = vec![
        file.packet(tag::COOKIE)
            .cloned()
            .ok_or_else(|| IoError::Format("no cookie".into()))?,
        Packet {
            tag: tag::SEQUENCE,
            payload: seq,
        },
    ];
    if let Some(record) = h.record(node.id) {
        packets.extend(
            record
                .packets
                .iter()
                .filter(|p| {
                    matches!(
                        p.tag,
                        tag::FEATURES | tag::PRIMERS | tag::NOTES | tag::PROPERTIES
                    )
                })
                .cloned(),
        );
    }
    let bytes = DnaFile::from_packets(packets)?.to_bytes()?;
    let mut m = crate::dna::import_plain(name, &bytes)?;
    m.set_history(None)?;
    Ok(m)
}

/// Write `molecule`'s construction history into `file` (replacing any there; none removes it).
/// Steps whose sequence the format can't hold (ambiguous bases) are listed but can't be opened;
/// `warnings` says what was left out.
pub(crate) fn write(
    file: &mut DnaFile,
    molecule: &Molecule,
    warnings: &mut Vec<String>,
) -> Result<(), IoError> {
    let Some(root) = molecule.history() else {
        file.set_history(None, Vec::new())?;
        return Ok(());
    };
    let mut records = Vec::new();
    let mut next_id = 0u32;
    let mut left_out = Left::default();
    let mut tree = out_node(root, 0, &mut next_id, &mut records, &mut left_out)?;
    // The root is the file itself: its sequence is the file's.
    tree.name = file_name(molecule.name());
    tree.length = molecule.len();
    tree.circular = molecule.topology() == Topology::Circular;
    tree.stored = false;
    tree.features = None;
    if left_out.ambiguous > 0 {
        warnings.push(format!(
            "{} earlier step(s) have bases other than A, C, G and T, which the history can't store: they are listed but can't be opened",
            left_out.ambiguous
        ));
    }
    if left_out.primers > 0 {
        warnings.push(format!(
            "the primers of {} earlier step(s) aren't written into the history",
            left_out.primers
        ));
    }
    file.set_history(Some(&tree_xml(&tree)?), records)?;
    Ok(())
}

#[derive(Default)]
struct Left {
    ambiguous: usize,
    primers: usize,
}

fn file_name(name: &str) -> String {
    if name.ends_with(".dna") {
        name.to_string()
    } else {
        format!("{name}.dna")
    }
}

/// A step and (first) its inputs, numbered inputs first so the root gets the highest ID, as in
/// the reference application's files.
fn out_node(
    n: &HistoryNode,
    depth: usize,
    next_id: &mut u32,
    records: &mut Vec<Vec<u8>>,
    left: &mut Left,
) -> Result<NodeOut, IoError> {
    let mut inputs = Vec::new();
    if depth < MAX_HISTORY_DEPTH {
        for i in &n.inputs {
            inputs.push(out_node(i, depth + 1, next_id, records, left)?);
        }
    }
    let id = *next_id;
    *next_id = next_id.saturating_add(1);
    let mut node = NodeOut {
        id,
        name: file_name(&n.name),
        length: n.length,
        circular: n.topology == Topology::Circular,
        operation: n.operation.clone(),
        stored: false,
        features: None,
        oligos: n.oligos.clone(),
        parameters: n.parameters.clone(),
        summaries: Vec::new(),
        colors: !n.operation.is_empty(),
        inputs,
    };
    // One summary per input, as the reference application writes them; where none was
    // recorded, the input taken whole.
    if !n.inputs.is_empty() {
        node.summaries = if n.summaries.len() == n.inputs.len() {
            n.summaries.clone()
        } else {
            n.inputs
                .iter()
                .map(|i| {
                    vec![
                        ("manipulation".to_string(), "overlapAndInsert".to_string()),
                        ("val1".to_string(), "0".to_string()),
                        ("val2".to_string(), i.length.saturating_sub(1).to_string()),
                    ]
                })
                .collect()
        };
    }
    if let Some(m) = n.molecule.as_deref() {
        let flags = if m.topology() == Topology::Circular {
            0x03
        } else {
            0x02
        };
        match encode_sequence(m.seq().as_bytes(), flags)? {
            None => left.ambiguous += 1,
            Some(seq) => {
                let mut plain = m.clone();
                plain.take_history();
                if !plain.primers().is_empty() {
                    left.primers += 1;
                    plain = without_primers(&plain)?;
                }
                let (bytes, _, _) = crate::dna::export(&plain, None)?;
                let mini = DnaFile::parse(&bytes)?;
                let features: Vec<Packet> = mini
                    .packets()
                    .iter()
                    .filter(|p| p.tag == tag::FEATURES)
                    .cloned()
                    .collect();
                node.features = features
                    .first()
                    .and_then(|p| features_element(&String::from_utf8_lossy(&p.payload)));
                node.length = m.len();
                node.circular = m.topology() == Topology::Circular;
                node.stored = true;
                records.push(node_record(id, &seq, &features)?);
            }
        }
    }
    Ok(node)
}

fn without_primers(m: &Molecule) -> Result<Molecule, IoError> {
    let mut out = m.clone();
    let ids: Vec<_> = out.primers().iter().map(|p| p.id).collect();
    for id in ids {
        out.remove_primer(id)?;
    }
    Ok(out)
}
