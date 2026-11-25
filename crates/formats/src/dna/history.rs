//! The construction history of a `.dna` file: how the construct was made, step by step, with
//! the sequences of the molecules it was made from (learned from public descriptions of the
//! format and from `.dna` files we authored ourselves, never from the reference app's code).
//!
//! - Packet 0x07 holds the tree: xz-compressed XML, `<HistoryTree>` with nested `<Node>`s. The
//!   root is the file's own construct; a node's child nodes are the inputs of the step that
//!   made it.
//! - Packets 0x0B hold one record per node: a node ID, then a packet stream with the node's
//!   sequence (0x01, two bits per base), its other packets (0x1e: features, primers, notes,
//!   xz-compressed) and, for nodes stored as an edit, the edit (0x1d: xz-compressed
//!   `<Manipulation>`).
//!
//! Everything is read with bounds and size limits: a damaged history is an error, never a
//! panic, and the rest of the file is unaffected (its packets are kept as read).

use std::io::Write;

use super::{DnaError, DnaFile, Packet};
use crate::xml::{Token, Tokenizer, decode, encode_attr};

/// Packet tags of the history.
pub mod tag {
    /// The history tree (xz-compressed XML).
    pub const TREE: u8 = 0x07;
    /// One node record.
    pub const NODE: u8 = 0x0B;
    /// In a node record: its sequence.
    pub const NODE_SEQUENCE: u8 = 0x01;
    /// In a node record: an edit relative to the step built from it (xz-compressed XML).
    pub const NODE_EDIT: u8 = 0x1d;
    /// In a node record: its other packets (xz-compressed packet stream).
    pub const NODE_PACKETS: u8 = 0x1e;
}

/// Largest decompressed history part we accept.
const MAX_DECOMPRESSED: usize = 512 << 20;
/// Deepest tree we read.
const MAX_DEPTH: usize = 200;
/// Most nodes we read.
const MAX_NODES: usize = 100_000;

fn bad(tag: u8, reason: impl Into<String>) -> DnaError {
    DnaError::InvalidPacket {
        tag,
        reason: reason.into(),
    }
}

/// A `Write` that refuses to grow past a limit (so a hostile stream can't exhaust memory).
struct Limited {
    out: Vec<u8>,
    limit: usize,
}

impl Write for Limited {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if self.out.len().saturating_add(buf.len()) > self.limit {
            return Err(std::io::Error::other("decompressed data too large"));
        }
        self.out.extend_from_slice(buf);
        Ok(buf.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Decompress an xz stream.
pub fn xz_decompress(data: &[u8], tag: u8) -> Result<Vec<u8>, DnaError> {
    let mut out = Limited {
        out: Vec::new(),
        limit: MAX_DECOMPRESSED,
    };
    let mut input = std::io::BufReader::new(data);
    lzma_rs::xz_decompress(&mut input, &mut out).map_err(|e| bad(tag, format!("xz: {e}")))?;
    Ok(out.out)
}

/// A packet stream without a cookie: tag, u32 length, payload, repeated.
pub fn packet_stream(bytes: &[u8], tag: u8) -> Result<Vec<Packet>, DnaError> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < bytes.len() {
        let header = bytes
            .get(i..i + 5)
            .ok_or_else(|| bad(tag, "truncated packet header in a record"))?;
        let t = header.first().copied().unwrap_or(0);
        let len = header
            .get(1..5)
            .and_then(|b| b.try_into().ok())
            .map(u32::from_be_bytes)
            .ok_or_else(|| bad(tag, "bad packet length"))? as usize;
        let start = i + 5;
        let payload = start
            .checked_add(len)
            .and_then(|end| bytes.get(start..end))
            .ok_or_else(|| bad(tag, "packet runs past the end of its record"))?;
        out.push(Packet {
            tag: t,
            payload: payload.to_vec(),
        });
        i = start + len;
    }
    Ok(out)
}

/// One step of the history.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct HistoryNode {
    pub id: u32,
    pub name: String,
    /// "DNA", "RNA", "Protein".
    pub kind: String,
    pub length: usize,
    pub circular: bool,
    /// The operation that made it ("" for a starting molecule).
    pub operation: String,
    /// Inputs of the operation, in file order.
    pub inputs: Vec<HistoryNode>,
    /// What each input contributed (`manipulation`, then any names and values), in file order.
    pub summaries: Vec<Vec<(String, String)>>,
    /// Primers of the step (`<Oligo>`: name, sequence).
    pub oligos: Vec<(String, String)>,
    /// Settings of the step (`<Parameter>`: name, value).
    pub parameters: Vec<(String, String)>,
    /// A step the reference app hides (`hideHistory="1"`); such steps may have no record.
    pub hidden: bool,
}

impl HistoryNode {
    /// Every node of the tree, depth first (this one first).
    pub fn walk(&self) -> Vec<&HistoryNode> {
        let mut out = vec![self];
        let mut i = 0;
        while let Some(n) = out.get(i) {
            let n: &HistoryNode = n;
            out.extend(n.inputs.iter());
            i += 1;
            if out.len() > MAX_NODES {
                break;
            }
        }
        out
    }

    /// The node with `id` and the node it is an input of (`None` for the root).
    pub fn find(&self, id: u32) -> Option<(&HistoryNode, Option<&HistoryNode>)> {
        if self.id == id {
            return Some((self, None));
        }
        let mut stack: Vec<&HistoryNode> = vec![self];
        while let Some(n) = stack.pop() {
            for input in &n.inputs {
                if input.id == id {
                    return Some((input, Some(n)));
                }
                stack.push(input);
            }
        }
        None
    }
}

/// Parse the history tree XML.
pub fn parse_tree(xml: &str) -> Result<HistoryNode, DnaError> {
    let t = tag::TREE;
    let mut tokens = Tokenizer::new(xml);
    // Open nodes, innermost last.
    let mut stack: Vec<HistoryNode> = Vec::new();
    let mut root: Option<HistoryNode> = None;
    let mut count = 0usize;
    while let Some(next) = tokens.next_token() {
        let (token, _) =
            next.map_err(|e| bad(t, format!("XML at byte {}: {}", e.offset, e.reason)))?;
        match token {
            Token::Start {
                name: "Node",
                attrs,
                self_closing,
            } => {
                count += 1;
                if count > MAX_NODES || stack.len() >= MAX_DEPTH {
                    return Err(bad(t, "history tree too large or too deep"));
                }
                let get = |n: &str| attrs.iter().find(|a| a.name == n).map(|a| decode(a.raw));
                let node = HistoryNode {
                    id: get("ID")
                        .and_then(|v| v.parse().ok())
                        .ok_or_else(|| bad(t, "node without ID"))?,
                    name: get("name").unwrap_or_default(),
                    kind: get("type").unwrap_or_else(|| "DNA".to_string()),
                    length: get("seqLen").and_then(|v| v.parse().ok()).unwrap_or(0),
                    circular: get("circular").as_deref() == Some("1"),
                    operation: match get("operation").unwrap_or_default().as_str() {
                        // A starting molecule: no operation made it.
                        "invalid" => String::new(),
                        other => other.to_string(),
                    },
                    inputs: Vec::new(),
                    summaries: Vec::new(),
                    oligos: Vec::new(),
                    parameters: Vec::new(),
                    hidden: get("hideHistory").as_deref() == Some("1"),
                };
                if self_closing {
                    match stack.last_mut() {
                        Some(parent) => parent.inputs.push(node),
                        None => root = Some(node),
                    }
                } else {
                    stack.push(node);
                }
            }
            Token::Start {
                name: "InputSummary",
                attrs,
                ..
            } => {
                if let Some(n) = stack.last_mut() {
                    n.summaries.push(
                        attrs
                            .iter()
                            .map(|a| (a.name.to_string(), decode(a.raw)))
                            .collect(),
                    );
                }
            }
            Token::Start {
                name: tag_name @ ("Oligo" | "Parameter"),
                attrs,
                ..
            } => {
                let get = |n: &str| attrs.iter().find(|a| a.name == n).map(|a| decode(a.raw));
                if let Some(n) = stack.last_mut() {
                    if tag_name == "Oligo" {
                        n.oligos.push((
                            get("name").unwrap_or_default(),
                            get("sequence").unwrap_or_default(),
                        ));
                    } else {
                        n.parameters.push((
                            get("name").unwrap_or_default(),
                            get("val").unwrap_or_default(),
                        ));
                    }
                }
            }
            Token::End { name: "Node" } => {
                let node = stack.pop().ok_or_else(|| bad(t, "unbalanced </Node>"))?;
                match stack.last_mut() {
                    Some(parent) => parent.inputs.push(node),
                    None => {
                        if root.is_none() {
                            root = Some(node);
                        }
                    }
                }
            }
            _ => {}
        }
    }
    root.ok_or_else(|| bad(t, "no <Node> in the history tree"))
}

/// A node's own record (packet 0x0B).
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct NodeRecord {
    pub id: u32,
    /// Its sequence, if stored (with lowercase where the original was).
    pub residues: Option<Vec<u8>>,
    /// Its other packets (features, primers, notes, ...), if stored.
    pub packets: Vec<Packet>,
    /// The edit from this node to the step built from it, if stored that way.
    pub edit: Option<Edit>,
}

/// An edit between two history nodes.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Edit {
    pub redo: Vec<Action>,
    pub undo: Vec<Action>,
}

/// One action of an edit. Positions are 0-based (as in the primers packet): `Insert` goes before
/// base `position`, `Remove` covers `first..=last`. A `Shift` in a Redo moves the origin: the
/// product is the node rotated left by `-by`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Action {
    Insert {
        position: usize,
        residues: Vec<u8>,
    },
    /// Remove `first..=last`.
    Remove {
        first: usize,
        last: usize,
    },
    /// Move the origin (see above).
    Shift {
        by: i64,
    },
    /// Properties (topology, ...) changed; the residues didn't.
    Properties,
    Other(String),
}

/// Decode a sequence record (0x01): u32 length, u8 flags (as packet 0x00), u32 segment count,
/// u32 count of lowercase ranges, then the segments (u8 kind, u32 bases, data), then the
/// lowercase ranges (u32 start, u32 end, inclusive). Kind 1 holds two bits per base (G=00,
/// A=01, T=10, C=11, high bits first, a final partial byte right-aligned); kind 2 one byte per
/// base, where 4 is N (the only value seen; others are refused rather than guessed); kind 3
/// a run of N with no data bytes. (Kinds 2 and 3 were confirmed against the reference app.)
pub fn decode_sequence(d: &[u8]) -> Result<Vec<u8>, DnaError> {
    let t = tag::NODE_SEQUENCE;
    let u32_at = |i: usize| -> Result<usize, DnaError> {
        i.checked_add(4)
            .and_then(|e| d.get(i..e))
            .and_then(|b| b.try_into().ok())
            .map(|b: [u8; 4]| u32::from_be_bytes(b) as usize)
            .ok_or_else(|| bad(t, "truncated sequence header"))
    };
    let length = u32_at(0)?;
    let segments = u32_at(5)?;
    let n_lower = u32_at(9)?;
    // Every base takes at least two bits, every segment at least five bytes: bound the
    // header's numbers by the record's size before allocating.
    if length > d.len().saturating_mul(4) || segments > d.len() / 5 {
        return Err(bad(t, "sequence longer than its record"));
    }
    let base = |code: u8| match code & 3 {
        0 => b'G',
        1 => b'A',
        2 => b'T',
        _ => b'C',
    };
    let mut out = Vec::with_capacity(length);
    let mut at = 13usize;
    for _ in 0..segments {
        let kind = *d
            .get(at)
            .ok_or_else(|| bad(t, "truncated sequence segment"))?;
        let n = u32_at(at + 1)?;
        at += 5;
        if n > length.saturating_sub(out.len()) {
            return Err(bad(t, "segments longer than the sequence"));
        }
        match kind {
            1 => {
                let bytes = n.div_ceil(4);
                let body = d
                    .get(at..at.saturating_add(bytes))
                    .ok_or_else(|| bad(t, "sequence shorter than its length"))?;
                let full = n / 4;
                for &byte in body.iter().take(full) {
                    for shift in [6, 4, 2, 0] {
                        out.push(base(byte >> shift));
                    }
                }
                let rem = n % 4;
                if rem > 0
                    && let Some(&last) = body.get(full)
                {
                    // The last partial byte is right-aligned.
                    for k in 0..rem {
                        out.push(base(last >> (2 * (rem - 1 - k))));
                    }
                }
                at += bytes;
            }
            2 => {
                let body = d
                    .get(at..at.saturating_add(n))
                    .ok_or_else(|| bad(t, "sequence shorter than its length"))?;
                for &code in body {
                    match code {
                        4 => out.push(b'N'),
                        other => {
                            return Err(bad(
                                t,
                                format!("unknown base code {other} in a one-byte segment"),
                            ));
                        }
                    }
                }
                at += n;
            }
            3 => out.extend(std::iter::repeat_n(b'N', n)),
            other => {
                return Err(bad(t, format!("unknown sequence segment kind {other}")));
            }
        }
    }
    if out.len() != length {
        return Err(bad(t, "segments don't add up to the sequence length"));
    }
    let trailer = d.get(at..).unwrap_or_default();
    if trailer.len() != n_lower.saturating_mul(8) {
        return Err(bad(t, "lowercase ranges don't match their count"));
    }
    for pair in trailer.as_chunks::<8>().0 {
        let [a0, a1, a2, a3, b0, b1, b2, b3] = *pair;
        let a = u32::from_be_bytes([a0, a1, a2, a3]) as usize;
        let b = u32::from_be_bytes([b0, b1, b2, b3]) as usize;
        if a > b || b >= length {
            return Err(bad(t, "lowercase range outside the sequence"));
        }
        for base in out.get_mut(a..=b).unwrap_or_default() {
            base.make_ascii_lowercase();
        }
    }
    Ok(out)
}

/// Parse an edit (0x1d, decompressed).
pub fn parse_edit(xml: &str) -> Result<Edit, DnaError> {
    let t = tag::NODE_EDIT;
    let mut edit = Edit::default();
    let mut tokens = Tokenizer::new(xml);
    // Which list the current actions go to, and the open INSERT waiting for its residues.
    let mut redo = true;
    let mut pending: Option<usize> = None;
    while let Some(next) = tokens.next_token() {
        let (token, _) =
            next.map_err(|e| bad(t, format!("XML at byte {}: {}", e.offset, e.reason)))?;
        match token {
            Token::Start { name: "Redo", .. } => redo = true,
            Token::Start { name: "Undo", .. } => redo = false,
            Token::Start {
                name: "Action",
                attrs,
                self_closing,
            } => {
                let get = |n: &str| attrs.iter().find(|a| a.name == n).map(|a| decode(a.raw));
                let kind = get("type").unwrap_or_default();
                let action = match kind.as_str() {
                    "INSERT" => {
                        let position = get("position")
                            .and_then(|v| v.parse().ok())
                            .ok_or_else(|| bad(t, "INSERT without position"))?;
                        if self_closing {
                            Some(Action::Insert {
                                position,
                                residues: Vec::new(),
                            })
                        } else {
                            pending = Some(position);
                            None
                        }
                    }
                    "REMOVE" => {
                        let range = get("range").unwrap_or_default();
                        let (a, b) = range
                            .split_once('-')
                            .ok_or_else(|| bad(t, "REMOVE without a range"))?;
                        Some(Action::Remove {
                            first: a.trim().parse().map_err(|_| bad(t, "bad REMOVE range"))?,
                            last: b.trim().parse().map_err(|_| bad(t, "bad REMOVE range"))?,
                        })
                    }
                    "SHIFT" => Some(Action::Shift {
                        by: get("by")
                            .and_then(|v| v.parse().ok())
                            .ok_or_else(|| bad(t, "SHIFT without by"))?,
                    }),
                    "PROPERTIES" => Some(Action::Properties),
                    other => Some(Action::Other(other.to_string())),
                };
                if let Some(a) = action {
                    if redo {
                        edit.redo.push(a)
                    } else {
                        edit.undo.push(a)
                    }
                }
            }
            Token::Start {
                name: "Residues",
                attrs,
                ..
            } => {
                if let Some(position) = pending.take() {
                    let residues = attrs
                        .iter()
                        .find(|a| a.name == "residues")
                        .map(|a| decode(a.raw).into_bytes())
                        .unwrap_or_default();
                    let a = Action::Insert { position, residues };
                    if redo {
                        edit.redo.push(a)
                    } else {
                        edit.undo.push(a)
                    }
                }
            }
            _ => {}
        }
    }
    Ok(edit)
}

/// Read one node record (packet 0x0B payload).
pub fn parse_node_record(payload: &[u8]) -> Result<NodeRecord, DnaError> {
    let t = tag::NODE;
    let id = payload
        .get(..4)
        .and_then(|b| b.try_into().ok())
        .map(u32::from_be_bytes)
        .ok_or_else(|| bad(t, "node record without an ID"))?;
    let mut record = NodeRecord {
        id,
        ..NodeRecord::default()
    };
    for p in packet_stream(payload.get(4..).unwrap_or_default(), t)? {
        match p.tag {
            tag::NODE_SEQUENCE => record.residues = Some(decode_sequence(&p.payload)?),
            tag::NODE_PACKETS => {
                record.packets = packet_stream(&xz_decompress(&p.payload, p.tag)?, p.tag)?;
            }
            tag::NODE_EDIT => {
                let xml = xz_decompress(&p.payload, p.tag)?;
                record.edit = Some(parse_edit(&String::from_utf8_lossy(&xml))?);
            }
            // Features, primers, notes, ... stored directly (not compressed).
            _ => record.packets.push(p),
        }
    }
    Ok(record)
}

/// The history of a file: the tree, and the node records it has.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct History {
    pub tree: HistoryNode,
    pub records: Vec<NodeRecord>,
}

impl History {
    pub fn record(&self, id: u32) -> Option<&NodeRecord> {
        self.records.iter().find(|r| r.id == id)
    }

    /// The residues of node `id`: stored, or rebuilt by undoing its edit on the step built from
    /// it (and so on up to the root, whose residues are `root_residues`). `None` if they can't
    /// be recovered.
    pub fn residues(&self, id: u32, root_residues: &[u8]) -> Option<Vec<u8>> {
        let mut chain = Vec::new();
        let mut current = id;
        // Walk up until a node with stored residues (or the root).
        let base = loop {
            if chain.len() > MAX_DEPTH {
                return None;
            }
            let (_, product) = self.tree.find(current)?;
            if let Some(r) = self.record(current).and_then(|r| r.residues.clone()) {
                break r;
            }
            match product {
                None => break root_residues.to_vec(),
                Some(p) => {
                    chain.push(current);
                    current = p.id;
                }
            }
        };
        // Back down: undo each edit on its product.
        let mut residues = base;
        for node in chain.iter().rev() {
            let edit = self.record(*node)?.edit.as_ref()?;
            residues = undo(residues, edit)?;
        }
        let expected = self.tree.find(id)?.0.length;
        (residues.len() == expected).then_some(residues)
    }
}

/// Rotate left by `k` (any integer, taken modulo the length).
fn rotate(s: &mut [u8], k: i64) -> Option<()> {
    let len = i64::try_from(s.len()).ok()?;
    if len > 0 {
        s.rotate_left(usize::try_from(k.rem_euclid(len)).ok()?);
    }
    Some(())
}

/// Apply an edit backwards: from the product's residues to the node's. The Undo list is used
/// when there is one (its actions lead from product to node); otherwise the Redo list is undone,
/// which works for origin shifts and property changes.
fn undo(mut s: Vec<u8>, edit: &Edit) -> Option<Vec<u8>> {
    if edit.undo.is_empty() {
        for a in edit.redo.iter().rev() {
            match a {
                // The product is the node rotated left by -by: rotate back by `by`.
                Action::Shift { by } => rotate(&mut s, *by)?,
                Action::Properties | Action::Other(_) => {}
                Action::Insert { .. } | Action::Remove { .. } => return None,
            }
        }
        return Some(s);
    }
    for a in &edit.undo {
        match a {
            Action::Insert { position, residues } => {
                if *position > s.len() {
                    return None;
                }
                s.splice(*position..*position, residues.iter().copied());
            }
            Action::Remove { first, last } => {
                if first > last || *last >= s.len() {
                    return None;
                }
                s.drain(*first..=*last);
            }
            Action::Shift { by } => rotate(&mut s, -*by)?,
            Action::Properties | Action::Other(_) => {}
        }
    }
    Some(s)
}

impl DnaFile {
    /// The construction history, if the file has one.
    pub fn history(&self) -> Result<Option<History>, DnaError> {
        let Some(tree) = self.packet(tag::TREE) else {
            return Ok(None);
        };
        let xml = xz_decompress(&tree.payload, tag::TREE)?;
        let tree = parse_tree(&String::from_utf8_lossy(&xml))?;
        let mut records = Vec::new();
        for p in self.packets().iter().filter(|p| p.tag == tag::NODE) {
            records.push(parse_node_record(&p.payload)?);
        }
        Ok(Some(History { tree, records }))
    }
}

// --- Writing -------------------------------------------------------------------------------

/// A step to write into a history tree.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct NodeOut {
    pub id: u32,
    /// File name style (`name.dna`), as the reference application writes it.
    pub name: String,
    pub length: usize,
    pub circular: bool,
    /// Operation id; empty for a starting molecule (written as `invalid`).
    pub operation: String,
    /// Whether a record with the step's sequence is written (it can be opened).
    pub stored: bool,
    /// The step's `<Features>` element (from its features packet), if any.
    pub features: Option<String>,
    /// PCR primers: (name, sequence).
    pub oligos: Vec<(String, String)>,
    /// Settings (`<Parameter name val>`), such as the polymerase of a PCR.
    pub parameters: Vec<(String, String)>,
    /// Whether to mark the whole molecule as made by the step (`<HistoryColors>`).
    pub colors: bool,
    /// One `<InputSummary>` per input: attribute name/value pairs.
    pub summaries: Vec<Vec<(String, String)>>,
    pub inputs: Vec<NodeOut>,
}

/// The history tree as XML (packet 0x07 before compression).
pub fn tree_xml(root: &NodeOut) -> Result<String, DnaError> {
    fn node(out: &mut String, n: &NodeOut, depth: usize) -> Result<(), DnaError> {
        if depth > MAX_DEPTH {
            return Err(bad(tag::TREE, "history tree too deep to write"));
        }
        let op = if n.operation.is_empty() {
            "invalid"
        } else {
            n.operation.as_str()
        };
        out.push_str(&format!(
            "<Node name=\"{}\" type=\"DNA\" seqLen=\"{}\" strandedness=\"double\" ID=\"{}\" circular=\"{}\"",
            encode_attr(&n.name),
            n.length,
            n.id,
            if n.circular { 1 } else { 0 }
        ));
        if !n.circular {
            out.push_str(
                " upstreamModification=\"Unmodified\" downstreamModification=\"Unmodified\"",
            );
        }
        if n.stored {
            out.push_str(" resurrectable=\"1\"");
        }
        out.push_str(&format!(" operation=\"{}\">", encode_attr(op)));
        if let Some(f) = &n.features {
            out.push_str(f);
        }
        for (name, seq) in &n.oligos {
            out.push_str(&format!(
                "<Oligo name=\"{}\" sequence=\"{}\"/>",
                encode_attr(name),
                encode_attr(seq)
            ));
        }
        for (name, val) in &n.parameters {
            out.push_str(&format!(
                "<Parameter name=\"{}\" val=\"{}\"/>",
                encode_attr(name),
                encode_attr(val)
            ));
        }
        if n.colors && n.length > 0 {
            let range = format!("<ColorRange range=\"0..{}\" colors=\"red\"/>", n.length - 1);
            out.push_str(&format!(
                "<HistoryColors><StrandColors type=\"Product\"><TopStrand>{range}</TopStrand><BottomStrand>{range}</BottomStrand></StrandColors></HistoryColors>"
            ));
        }
        for summary in &n.summaries {
            out.push_str("<InputSummary");
            for (k, v) in summary {
                // Attribute names come from our own code or a parsed file: keep only safe ones.
                if !k.is_empty() && k.chars().all(|c| c.is_ascii_alphanumeric()) {
                    out.push_str(&format!(" {k}=\"{}\"", encode_attr(v)));
                }
            }
            out.push_str("/>");
        }
        for i in &n.inputs {
            node(out, i, depth + 1)?;
        }
        out.push_str("</Node>");
        Ok(())
    }
    let mut out = String::from(r#"<?xml version="1.0" encoding="UTF-8"?><HistoryTree>"#);
    node(&mut out, root, 0)?;
    out.push_str("</HistoryTree>");
    Ok(out)
}

/// The `<Features>` element of a features packet, for a tree node (without the XML
/// declaration).
pub fn features_element(features_packet: &str) -> Option<String> {
    let start = features_packet.find("<Features")?;
    Some(features_packet.get(start..)?.trim_end().to_string())
}

/// Compress with xz.
pub fn xz_compress(data: &[u8], tag: u8) -> Result<Vec<u8>, DnaError> {
    let mut out = Vec::new();
    let mut input = std::io::BufReader::new(data);
    lzma_rs::xz_compress(&mut input, &mut out).map_err(|e| bad(tag, format!("xz: {e}")))?;
    Ok(out)
}

fn u32_of(n: usize, tag: u8) -> Result<[u8; 4], DnaError> {
    u32::try_from(n)
        .map(u32::to_be_bytes)
        .map_err(|_| bad(tag, "too long for the file format"))
}

/// Encode a node's sequence (0x01; layout in [`decode_sequence`]): runs of A, C, G and T as
/// two-bit segments, a single N as a one-byte segment and longer runs of N as data-free runs
/// (as reference files have them), lowercase stretches listed after them.
/// `None` if a residue is another ambiguity code (the format's code for it isn't known).
pub fn encode_sequence(residues: &[u8], flags: u8) -> Result<Option<Vec<u8>>, DnaError> {
    let t = tag::NODE_SEQUENCE;
    let len = u32_of(residues.len(), t)?;
    let mut lower: Vec<(usize, usize)> = Vec::new();
    // (kind, bases, data)
    let mut segments: Vec<(u8, usize, Vec<u8>)> = Vec::new();
    let mut byte = 0u8;
    let mut in_byte = 0usize;
    for (i, &b) in residues.iter().enumerate() {
        let up = b.to_ascii_uppercase();
        let code = match up {
            b'G' => Some(0),
            b'A' => Some(1),
            b'T' => Some(2),
            b'C' => Some(3),
            b'N' => None,
            _ => return Ok(None),
        };
        if b.is_ascii_lowercase() {
            match lower.last_mut() {
                Some((_, end)) if *end + 1 == i => *end = i,
                _ => lower.push((i, i)),
            }
        }
        let kind = if code.is_some() { 1 } else { 2 };
        if segments.last().is_none_or(|s| s.0 != kind) {
            // Close a two-bit segment's last, partial byte (right-aligned).
            if let Some(last) = segments.last_mut()
                && last.0 == 1
                && in_byte > 0
            {
                last.2.push(byte);
            }
            segments.push((kind, 0, Vec::new()));
            byte = 0;
            in_byte = 0;
        }
        let Some(seg) = segments.last_mut() else {
            continue;
        };
        seg.1 += 1;
        match code {
            Some(c) => {
                byte = (byte << 2) | c;
                in_byte += 1;
                if in_byte == 4 {
                    seg.2.push(byte);
                    byte = 0;
                    in_byte = 0;
                }
            }
            None => seg.2.push(4),
        }
    }
    if let Some(last) = segments.last_mut()
        && last.0 == 1
        && in_byte > 0
    {
        last.2.push(byte);
    }
    // Runs of more than one N need no data.
    for seg in &mut segments {
        if seg.0 == 2 && seg.1 > 1 {
            seg.0 = 3;
            seg.2.clear();
        }
    }
    let data: usize = segments.iter().map(|s| 5 + s.2.len()).sum();
    let mut out = Vec::with_capacity(13 + data + lower.len() * 8);
    out.extend_from_slice(&len);
    out.push(flags);
    out.extend_from_slice(&u32_of(segments.len(), t)?);
    out.extend_from_slice(&u32_of(lower.len(), t)?);
    for (kind, n, bytes) in &segments {
        out.push(*kind);
        out.extend_from_slice(&u32_of(*n, t)?);
        out.extend_from_slice(bytes);
    }
    for (a, b) in lower {
        out.extend_from_slice(&u32_of(a, t)?);
        out.extend_from_slice(&u32_of(b, t)?);
    }
    Ok(Some(out))
}

/// A node record (0x0B payload): the ID, the sequence, then the other packets (features,
/// primers, notes), xz-compressed together in one 0x1e packet as in recent reference files.
pub fn node_record(id: u32, sequence: &[u8], packets: &[Packet]) -> Result<Vec<u8>, DnaError> {
    fn push(out: &mut Vec<u8>, t: u8, payload: &[u8]) -> Result<(), DnaError> {
        out.push(t);
        out.extend_from_slice(&u32_of(payload.len(), tag::NODE)?);
        out.extend_from_slice(payload);
        Ok(())
    }
    let mut out = id.to_be_bytes().to_vec();
    push(&mut out, tag::NODE_SEQUENCE, sequence)?;
    if !packets.is_empty() {
        let mut stream = Vec::new();
        for p in packets {
            push(&mut stream, p.tag, &p.payload)?;
        }
        push(
            &mut out,
            tag::NODE_PACKETS,
            &xz_compress(&stream, tag::NODE_PACKETS)?,
        )?;
    }
    Ok(out)
}

impl DnaFile {
    /// Replace the history: `tree` is the uncompressed tree XML (`None` removes the history),
    /// `records` the node records. Records go after the sequence packets and the tree after
    /// them, as in the reference application's files; every other packet keeps its place.
    pub fn set_history(
        &mut self,
        tree: Option<&str>,
        records: Vec<Vec<u8>>,
    ) -> Result<(), DnaError> {
        let mut packets: Vec<Packet> = self
            .packets()
            .iter()
            .filter(|p| p.tag != tag::TREE && p.tag != tag::NODE)
            .cloned()
            .collect();
        if let Some(xml) = tree {
            let at = packets
                .iter()
                .rposition(|p| matches!(p.tag, 0x09 | 0x00 | 0x03))
                .map_or(packets.len(), |i| i + 1);
            let mut new: Vec<Packet> = records
                .into_iter()
                .map(|payload| Packet {
                    tag: tag::NODE,
                    payload,
                })
                .collect();
            new.push(Packet {
                tag: tag::TREE,
                payload: xz_compress(xml.as_bytes(), tag::TREE)?,
            });
            packets.splice(at..at, new);
        }
        *self = DnaFile::from_packets(packets)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2-bit encoding of `s` in the record layout (no lowercase ranges unless given).
    fn encode(s: &[u8], lower: &[(u32, u32)]) -> Vec<u8> {
        let len = s.len() as u32;
        let mut d = len.to_be_bytes().to_vec();
        d.push(0x1f);
        d.extend(1u32.to_be_bytes());
        d.extend((lower.len() as u32).to_be_bytes());
        d.push(1);
        d.extend(len.to_be_bytes());
        let code = |b: u8| match b.to_ascii_uppercase() {
            b'G' => 0,
            b'A' => 1,
            b'T' => 2,
            _ => 3,
        };
        let chunks = s.chunks(4);
        for c in chunks {
            let mut byte = 0u8;
            for &b in c {
                byte = (byte << 2) | code(b);
            }
            d.push(byte);
        }
        for (a, b) in lower {
            d.extend(a.to_be_bytes());
            d.extend(b.to_be_bytes());
        }
        d
    }

    #[test]
    fn segmented_sequences() {
        // As in reference files: GATCA, one N, then CC; one lowercase range (the N).
        let mut d = Vec::new();
        d.extend(8u32.to_be_bytes());
        d.push(0x1e);
        d.extend(3u32.to_be_bytes());
        d.extend(1u32.to_be_bytes());
        d.extend([1, 0, 0, 0, 5, 0b0001_1011, 0b01]); // GATC, then A right-aligned
        d.extend([2, 0, 0, 0, 1, 4]);
        d.extend([1, 0, 0, 0, 2, 0b1111]);
        d.extend(5u32.to_be_bytes());
        d.extend(5u32.to_be_bytes());
        assert_eq!(decode_sequence(&d).unwrap(), b"GATCAnCC");
        // Unknown codes and kinds are refused, not guessed.
        let mut odd = d.clone();
        odd[25] = 9;
        assert!(
            decode_sequence(&odd)
                .unwrap_err()
                .to_string()
                .contains("base code 9")
        );
        let mut kind = d.clone();
        kind[13] = 5;
        assert!(
            decode_sequence(&kind)
                .unwrap_err()
                .to_string()
                .contains("kind 5")
        );
        // A run of N carries no data: 3 bases, then a lowercase range over the last two.
        let mut run = Vec::new();
        run.extend(3u32.to_be_bytes());
        run.push(0x02);
        run.extend(1u32.to_be_bytes());
        run.extend(1u32.to_be_bytes());
        run.extend([3, 0, 0, 0, 3]);
        run.extend(1u32.to_be_bytes());
        run.extend(2u32.to_be_bytes());
        assert_eq!(decode_sequence(&run).unwrap(), b"Nnn");
        // A run can't claim more bases than the sequence has.
        let mut huge = run.clone();
        huge[14..18].copy_from_slice(&u32::MAX.to_be_bytes());
        assert!(decode_sequence(&huge).is_err());
        // Segments that don't add up, and hostile counts.
        let mut short = d.clone();
        short[3] = 9;
        assert!(decode_sequence(&short).is_err());
        let mut many = d.clone();
        many[5..9].copy_from_slice(&u32::MAX.to_be_bytes());
        assert!(decode_sequence(&many).is_err());
        let mut long = d.clone();
        long[0..4].copy_from_slice(&u32::MAX.to_be_bytes());
        assert!(decode_sequence(&long).is_err());
        for cut in 0..d.len() {
            assert!(decode_sequence(&d[..cut]).is_err(), "cut at {cut}");
        }
    }

    #[test]
    fn two_bit_sequences() {
        for s in [&b"GATC"[..], b"GATCA", b"ACGTACGTT", b"C", b""] {
            assert_eq!(decode_sequence(&encode(s, &[])).unwrap(), s.to_vec());
        }
        assert_eq!(
            decode_sequence(&encode(b"GATCAA", &[(1, 3)])).unwrap(),
            b"GatcAA"
        );
        // Damaged records are errors.
        let good = encode(b"GATCAA", &[(1, 3)]);
        for cut in 0..good.len() {
            assert!(decode_sequence(&good[..cut]).is_err(), "cut at {cut}");
        }
        assert!(decode_sequence(&encode(b"GATC", &[(2, 9)])).is_err());
    }

    #[test]
    fn tree_and_edits() {
        let xml = r#"<?xml version="1.0"?><HistoryTree><Node name="c" type="DNA" seqLen="10" ID="3" circular="1" operation="insert"><InputSummary manipulation="insert" val1="2"/><Node name="b" seqLen="8" ID="2" operation="setOrigin"><Node name="a" seqLen="8" ID="1" operation="invalid"/></Node></Node></HistoryTree>"#;
        let tree = parse_tree(xml).unwrap();
        assert_eq!((tree.id, tree.inputs.len(), tree.circular), (3, 1, true));
        let a = &tree.inputs[0].inputs[0];
        assert_eq!((a.id, a.operation.as_str()), (1, ""));
        assert_eq!(
            tree.summaries[0][0],
            ("manipulation".to_string(), "insert".to_string())
        );
        assert_eq!(tree.walk().len(), 3);
        assert_eq!(tree.find(1).unwrap().1.unwrap().id, 2);
        // Root "ACGTTTACGT" was made from b by inserting TT before base 4; b from a by a SHIFT
        // of -3.
        let ins = parse_edit(r#"<Manipulation><Redo><Action type="INSERT" position="4"><Residues type="GENERIC" residues="TT"/></Action></Redo><Undo><Action type="REMOVE" range="4-5"/></Undo></Manipulation>"#).unwrap();
        assert_eq!(ins.undo, vec![Action::Remove { first: 4, last: 5 }]);
        let shift = parse_edit(
            r#"<Manipulation><Redo><Action type="SHIFT" by="-3"/></Redo></Manipulation>"#,
        )
        .unwrap();
        let h = History {
            tree,
            records: vec![
                NodeRecord {
                    id: 2,
                    edit: Some(ins),
                    ..NodeRecord::default()
                },
                NodeRecord {
                    id: 1,
                    edit: Some(shift),
                    ..NodeRecord::default()
                },
            ],
        };
        let root = b"ACGTTTACGT";
        assert_eq!(h.residues(2, root).unwrap(), b"ACGTACGT");
        // b = a shifted by -3 (a's residue at -3 became first); so a = b shifted back by 3.
        // b "ACGTACGT" came from a by moving the origin by -3, so a is b shifted by +3.
        // b "ACGTACGT" is a rotated left by 3 (SHIFT -3), so a is b rotated left by -3.
        assert_eq!(h.residues(1, root).unwrap(), b"CGTACGTA".to_vec());
        assert_eq!(h.residues(99, root), None);
    }

    #[test]
    fn hostile_input() {
        for xml in [
            "",
            "<HistoryTree>",
            "<Node>",
            "<Node ID=\"x\"/>",
            "</Node>",
            &"<Node ID=\"1\">".repeat(300),
        ] {
            let _ = parse_tree(xml);
        }
        let _ = parse_edit(r#"<Action type="REMOVE" range="9-"/><Action type="SHIFT" by="a"/>"#);
        assert!(xz_decompress(b"not xz", 7).is_err());
        assert!(packet_stream(&[1, 0, 0, 0, 9, 1], 0x0b).is_err());
        assert!(parse_node_record(&[0, 0]).is_err());
        let _ = parse_node_record(&[0, 0, 0, 1, 1, 0, 0, 0, 2, 0xff, 0xff]);
        // A REMOVE outside the sequence can't be undone.
        let e = Edit {
            redo: vec![],
            undo: vec![Action::Remove { first: 5, last: 99 }],
        };
        assert_eq!(undo(b"ACGT".to_vec(), &e), None);
    }

    #[test]
    fn written_sequences_read_back() {
        for s in [
            &b""[..],
            b"G",
            b"ACGT",
            b"ACGTA",
            b"acgTTGca",
            b"GATTACAGATTACAgattaca",
        ] {
            let enc = encode_sequence(s, 0x03).unwrap().unwrap();
            assert_eq!(
                decode_sequence(&enc).unwrap(),
                s,
                "{}",
                String::from_utf8_lossy(s)
            );
        }
        // Matches the layout of the stored records.
        assert_eq!(
            encode_sequence(b"ACGTAcg", 0x1f).unwrap().unwrap(),
            encode(b"ACGTACG", &[(5, 6)])
                .into_iter()
                .collect::<Vec<u8>>()
        );
        // N is written as a one-byte segment (code 4); other ambiguity codes can't be.
        let n = encode_sequence(b"ACNNTgn", 0x03).unwrap().unwrap();
        assert_eq!(decode_sequence(&n).unwrap(), b"ACNNTgn");
        // Segments: AC (2-bit), NN (a run), Tg (2-bit), n (one byte).
        assert_eq!(&n[5..9], &4u32.to_be_bytes());
        assert_eq!(&n[19..24], &[3, 0, 0, 0, 2]);
        assert_eq!(encode_sequence(b"ACRT", 0x03).unwrap(), None);
        assert_eq!(
            decode_sequence(&encode_sequence(b"NNNN", 0x03).unwrap().unwrap()).unwrap(),
            b"NNNN"
        );
    }

    #[test]
    fn written_histories_read_back() {
        let leaf = |id, name: &str, circular| {
            NodeOut {
            id,
            name: name.into(),
            length: 8,
            circular,
            stored: true,
            features: Some(r#"<Features nextValidID="1"><Feature name="f &amp; g" type="misc_feature"><Segment range="1-4" type="standard"/></Feature></Features>"#.into()),
            ..NodeOut::default()
        }
        };
        let mut pcr = leaf(2, "pcr.dna", false);
        pcr.operation = "amplifyFragment".into();
        pcr.oligos = vec![("F".into(), "ACGT".into()), ("R".into(), "TTGG".into())];
        pcr.inputs = vec![leaf(1, "template \"x\".dna", true)];
        let root = NodeOut {
            id: 3,
            name: "product.dna".into(),
            length: 16,
            circular: true,
            operation: "gibsonAssembly".into(),
            inputs: vec![leaf(0, "vector.dna", true), pcr],
            ..NodeOut::default()
        };
        let xml = tree_xml(&root).unwrap();
        let tree = parse_tree(&xml).unwrap();
        assert_eq!(tree.id, 3);
        assert_eq!(tree.operation, "gibsonAssembly");
        let ids: Vec<u32> = tree.walk().iter().map(|n| n.id).collect();
        assert_eq!(ids, vec![3, 0, 2, 1]);
        let (template, _) = tree.find(1).unwrap();
        assert_eq!(
            (
                template.name.as_str(),
                template.circular,
                template.operation.as_str()
            ),
            ("template \"x\".dna", true, "")
        );

        // Into a file: records before the tree, after the sequence.
        let mut file = DnaFile::from_packets(vec![
            Packet {
                tag: 0x09,
                payload: b"SnapGene\x00\x01\x00\x0f\x00\x13".to_vec(),
            },
            Packet {
                tag: 0x00,
                payload: b"\x03ACGTACGTACGTACGT".to_vec(),
            },
            Packet {
                tag: 0x0a,
                payload: b"<Features/>".to_vec(),
            },
        ])
        .unwrap();
        let seq = encode_sequence(b"ACGTacgt", 0x03).unwrap().unwrap();
        let feats = Packet {
            tag: 0x0a,
            payload: b"<Features/>".to_vec(),
        };
        let records = vec![
            node_record(0, &seq, std::slice::from_ref(&feats)).unwrap(),
            node_record(1, &seq, &[]).unwrap(),
        ];
        file.set_history(Some(&xml), records).unwrap();
        let tags: Vec<u8> = file.packets().iter().map(|p| p.tag).collect();
        assert_eq!(tags, vec![0x09, 0x00, 0x0b, 0x0b, 0x07, 0x0a]);
        let file = DnaFile::parse(&file.to_bytes().unwrap()).unwrap();
        let h = file.history().unwrap().unwrap();
        assert_eq!(h.tree, tree);
        assert_eq!(h.residues(0, b"").unwrap(), b"ACGTacgt");
        assert_eq!(h.record(0).unwrap().packets, vec![feats]);
        // Removing it leaves the rest alone.
        let mut file = file;
        file.set_history(None, Vec::new()).unwrap();
        let tags: Vec<u8> = file.packets().iter().map(|p| p.tag).collect();
        assert_eq!(tags, vec![0x09, 0x00, 0x0a]);
    }
}
