//! Cloning simulation: restriction-ligation and Gateway. Each input is the open document or a
//! file; the result is previewed or opened as a new, unsaved document.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::assembly::{golden_gate, homology_assembly, mutagenesis, to_molecule};
use snappiesnap_bio::cloning::{
    GatewayReaction, Part, assemble, att_kinds, digest_fragments, extract, gateway,
    restriction_clone,
};
use snappiesnap_bio::pcr::{DEFAULT_MAX_PRODUCT, pcr, product_molecule};
use snappiesnap_doc::{HistoryNode, Molecule, Topology};

use crate::enzymes::{SiteQuery, sites};
use crate::{CmdError, CommandSpec, Document, ErrorCode, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "restriction_clone",
            "Restriction Cloning",
            "Cut a vector and an insert with restriction enzymes and ligate the vector's largest fragment with each insert fragment whose ends fit (both orientations where they do; directional clones first). Each of `vector` and `insert` is a file path, `vector_document`/`insert_document` an open document's id (see `documents`); omit both for the active document. Without insert enzymes a linear insert goes in blunt. With `preview`, lists the constructs (each with its `molecule`); otherwise opens construct `product` (default 0) as a new document.",
            restriction,
        )
        .menu("Actions/Restriction Cloning…")
        .example(json!({"insert": "insert.fa", "vector_enzymes": ["HindIII"], "insert_enzymes": ["HindIII"], "preview": true}))
        .example(json!({"insert": "insert.fa", "vector_enzymes": ["HindIII"], "insert_enzymes": ["HindIII"], "product": 1, "discard_changes": true})),
        CommandSpec::session(
            "gateway_clone",
            "Gateway Cloning",
            "Gateway LR (entry clone into destination vector, giving attB sites) or BP (attB insert into donor vector, giving attL sites): the stretch between the donor's att1 and att2 cores replaces the acceptor's. `acceptor` (destination or donor vector) and `donor` (entry clone or attB insert) are file paths, or `acceptor_document`/`donor_document` an open document's id (see `documents`); omit both for the active document. With `preview`, reports the product; otherwise opens it as a new document.",
            gateway_command,
        )
        .menu("Actions/Gateway Cloning…")
        .example(json!({"acceptor": "dest.gb", "donor": "entry.gb", "preview": true}))
        .example(json!({"acceptor": "dest.gb", "donor": "entry.gb", "reaction": "lr", "discard_changes": true})),
        CommandSpec::session(
            "assemble",
            "Assemble",
            "Gibson/HiFi or In-Fusion (join pieces where the end of one matches the start of the next over at least `min_overlap` bp; default 20 for Gibson, 15 for In-Fusion; a piece is turned around if only that fits) or Golden Gate (cut every part with a Type IIS `enzyme` and chain the fragments by overhang, the first part being the vector). Each part is a file, an open document (`document`: its id, see `documents`) or the active document, used whole, cut with `enzymes` (largest fragment kept), or the PCR product of two of its primers. With `preview`, reports the construct and its junctions; otherwise opens it.",
            assemble_command,
        )
        .menu("Actions/Assemble…")
        .example(json!({"method": "gibson", "parts": [{"file": "gibson_a.fa"}, {"file": "gibson_b.fa"}], "preview": true}))
        .example(json!({"method": "golden_gate", "enzyme": "BsaI", "parts": [{"file": "gg_vector.gb"}, {"file": "gg_part.fa"}], "discard_changes": true})),
        CommandSpec::session(
            "mutagenesis",
            "Site-Directed Mutagenesis",
            "Copy a circular plasmid (`template`: a file, or `template_document` an open document's id; default the active document) with back-to-back primers (a primer's name or a sequence; changes in their 5' parts become the mutation, a gap between them a deletion, 5' tails an insertion) and ligate the product closed. With `preview`, reports the result; otherwise opens it.",
            mutagenesis_command,
        )
        .menu("Actions/Site-Directed Mutagenesis…")
        .example(json!({"forward": "f", "reverse": "r", "preview": true}))
        .example(json!({"template_document": 1, "forward": "f", "reverse": "r", "preview": true}))
        .example(json!({"template": "mut_template.gb", "forward": "ATCGGATCGATCGTAGCTAG", "reverse": "GCTAGCTCGAAATGCATGCC", "preview": true})),
    ]
}

/// A molecule from a file, or the open document's.
pub(super) fn load(s: &Session, path: Option<&str>) -> Result<Molecule, CmdError> {
    match path {
        None => Ok(s.molecule()?.clone()),
        Some(p) => {
            let bytes = s
                .platform()
                .read(p)
                .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {p}: {e}")))?;
            Ok(snappiesnap_io::import(p, &bytes)?.molecule)
        }
    }
}

/// A molecule from a file, from an open document by its tab id, or (both left out) from the
/// active document. Naming both a file and a document is a mistake worth reporting rather
/// than silently preferring one: the two can hold different sequences.
pub(super) fn load_source(
    s: &Session,
    path: Option<&str>,
    document: Option<u64>,
    what: &str,
) -> Result<Molecule, CmdError> {
    match (path, document) {
        (Some(_), Some(_)) => Err(CmdError::invalid_params(format!(
            "give {what} as a file or as a document, not both"
        ))),
        (None, Some(id)) => s
            .document_by_id(id)
            .map(|d| d.molecule.clone())
            .ok_or_else(|| {
                CmdError::new(ErrorCode::NotFound, format!("no open document {id}"))
                    .with_hint("`documents` lists the open documents and their ids")
            }),
        (path, None) => load(s, path),
    }
}

/// Enzymes by name, from the session's database.
fn enzymes(
    s: &Session,
    m: &Molecule,
    names: &[String],
) -> Result<Vec<snappiesnap_bio::enzymes::Enzyme>, CmdError> {
    if names.is_empty() {
        return Ok(Vec::new());
    }
    let q = SiteQuery {
        enzymes: Some(names.to_vec()),
        ..SiteQuery::default()
    };
    let chosen = sites(s.enzymes(), m, &q).map_err(|e| CmdError::new(ErrorCode::NotFound, e))?;
    Ok(chosen
        .iter()
        .filter_map(|c| s.enzymes().get(&c.name).map(|(_, e)| e.clone()))
        .collect())
}

fn product_name(a: &str, b: &str) -> String {
    let n = format!("{a} + {b}");
    n.chars().take(200).collect()
}

/// Record in `m`'s history that `operation` made it from `inputs`.
/// `summaries` says what each input contributed (one per input, in the `.dna` style; see
/// [`summary`]).
pub(super) fn record(
    m: &mut Molecule,
    operation: &str,
    details: &str,
    inputs: Vec<HistoryNode>,
    summaries: Vec<Vec<(String, String)>>,
) -> Result<(), CmdError> {
    let mut h = HistoryNode::made(m, operation, details, inputs);
    if summaries.len() == h.inputs.len() {
        h.summaries = summaries;
    }
    m.set_history(Some(h))?;
    Ok(())
}

/// An input summary: the manipulation, then attribute pairs.
pub(super) fn summary(manipulation: &str, attrs: &[(&str, String)]) -> Vec<(String, String)> {
    let mut v = vec![("manipulation".to_string(), manipulation.to_string())];
    v.extend(attrs.iter().map(|(k, val)| (k.to_string(), val.clone())));
    v
}

/// The amplified stretch of a PCR template (0-based, inclusive; start > end when it wraps).
pub(super) fn amplified(
    product: &snappiesnap_bio::pcr::Product,
    len: usize,
) -> Vec<(String, String)> {
    let start = product.template.start;
    let last = (product.template.end + len.max(1) - 1) % len.max(1);
    summary(
        "amplify",
        &[("val1", start.to_string()), ("val2", last.to_string())],
    )
}

/// Note the primers of the step that made `m` (and for a PCR, its polymerase, as the
/// reference application records it).
pub(super) fn set_primers(
    m: &mut Molecule,
    primers: &[(&str, &[u8])],
    pcr: bool,
) -> Result<(), CmdError> {
    if let Some(mut h) = m.take_history() {
        h.oligos = primers
            .iter()
            .map(|(name, seq)| (name.to_string(), String::from_utf8_lossy(seq).into_owned()))
            .collect();
        if pcr {
            h.parameters = vec![(
                "polymerase".to_string(),
                "polymerase that creates blunt ends".to_string(),
            )];
        }
        m.set_history(Some(h))?;
    }
    Ok(())
}

/// An input used whole (assemblies).
fn whole(node: &HistoryNode) -> Vec<(String, String)> {
    summary(
        "overlapAndInsert",
        &[
            ("val1", "0".to_string()),
            ("val2", node.length.saturating_sub(1).to_string()),
        ],
    )
}

/// The first and last enzyme of a list, as `name1`/`name2`.
fn enzyme_names(names: &[String]) -> Vec<(&'static str, String)> {
    let mut v = Vec::new();
    if let Some(first) = names.first() {
        v.push(("name1", first.clone()));
        v.push(("name2", names.last().unwrap_or(first).clone()));
    }
    v
}

/// Open `m` as a new, unsaved document.
pub(super) fn open_product(s: &mut Session, m: Molecule, discard: bool) -> Result<Value, CmdError> {
    let _ = discard; // Constructs open in a new tab; nothing is replaced.
    let result = json!({
        "name": m.name(),
        "length": m.len(),
        "topology": m.topology(),
        "features": m.features().len(),
    });
    let mut doc = Document::new(m);
    doc.dirty = true;
    s.open_document(doc);
    Ok(result)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct RestrictionParams {
    /// Vector file (default: the active document).
    vector: Option<String>,
    /// Vector as an open document's id instead of a file (see `documents`).
    vector_document: Option<u64>,
    /// Insert file (default: the active document).
    insert: Option<String>,
    /// Insert as an open document's id instead of a file.
    insert_document: Option<u64>,
    /// Enzymes that cut the vector.
    vector_enzymes: Vec<String>,
    /// Enzymes that cut the insert (default: the vector's).
    insert_enzymes: Option<Vec<String>>,
    /// Only list the constructs.
    #[serde(default)]
    preview: bool,
    /// Which construct to open (default 0, the first).
    #[serde(default)]
    product: usize,
    #[serde(default)]
    discard_changes: bool,
}

fn restriction(s: &mut Session, p: RestrictionParams) -> Result<Value, CmdError> {
    if p.vector.is_none()
        && p.insert.is_none()
        && p.vector_document.is_none()
        && p.insert_document.is_none()
    {
        return Err(CmdError::invalid_params(
            "give `vector` or `insert` as a file or an open document's id (the other is the \
             active document)",
        ));
    }
    if p.vector_enzymes.is_empty() {
        return Err(CmdError::invalid_params(
            "name the enzymes that cut the vector",
        ));
    }
    let vector = load_source(s, p.vector.as_deref(), p.vector_document, "the vector")?;
    let insert = load_source(s, p.insert.as_deref(), p.insert_document, "the insert")?;
    let venz = enzymes(s, &vector, &p.vector_enzymes)?;
    let inames = p
        .insert_enzymes
        .clone()
        .unwrap_or_else(|| p.vector_enzymes.clone());
    let ienz = enzymes(s, &insert, &inames)?;
    let constructs = restriction_clone(&vector, &venz, &insert, &ienz).map_err(|e| {
        CmdError::new(ErrorCode::InvalidInput, e)
            .with_hint("check the enzymes: the ends of the vector and insert fragments must match")
    })?;
    if p.preview {
        let name = product_name(vector.name(), insert.name());
        let list: Vec<Value> = constructs
            .iter()
            .map(|c| {
                let molecule = assemble(std::slice::from_ref(&c.part), Topology::Circular, &name)
                    .ok()
                    .and_then(|m| serde_json::to_value(m).ok());
                json!({"length": c.part.seq.len(), "description": c.description, "molecule": molecule})
            })
            .collect();
        return Ok(json!({ "constructs": list }));
    }
    let c = constructs.get(p.product).ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotFound,
            format!("no construct {} ({} possible)", p.product, constructs.len()),
        )
    })?;
    let mut m = assemble(
        std::slice::from_ref(&c.part),
        Topology::Circular,
        &product_name(vector.name(), insert.name()),
    )?;
    let details = format!(
        "vector cut with {}, insert cut with {}; {}",
        p.vector_enzymes.join(" + "),
        inames.join(" + "),
        c.description
    );
    let summaries = vec![
        summary("replace", &enzyme_names(&p.vector_enzymes)),
        summary("insert", &enzyme_names(&inames)),
    ];
    record(
        &mut m,
        "insertFragment",
        &details,
        vec![HistoryNode::input(&vector), HistoryNode::input(&insert)],
        summaries,
    )?;
    let mut out = open_product(s, m, p.discard_changes)?;
    if let Some(o) = out.as_object_mut() {
        o.insert("description".into(), json!(c.description));
    }
    Ok(out)
}

/// Which Gateway reaction.
#[derive(Deserialize, JsonSchema, Clone, Copy)]
#[serde(rename_all = "lowercase")]
enum Reaction {
    /// Entry clone × destination vector → expression clone (attB sites).
    Lr,
    /// attB insert × donor vector → entry clone (attL sites).
    Bp,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct GatewayParams {
    /// Destination (LR) or donor (BP) vector file (default: the active document).
    acceptor: Option<String>,
    /// The acceptor as an open document's id instead of a file (see `documents`).
    acceptor_document: Option<u64>,
    /// Entry clone (LR) or attB insert (BP) file (default: the active document).
    donor: Option<String>,
    /// The donor as an open document's id instead of a file.
    donor_document: Option<u64>,
    /// "lr" (default) or "bp".
    reaction: Option<Reaction>,
    #[serde(default)]
    preview: bool,
    #[serde(default)]
    discard_changes: bool,
}

fn gateway_command(s: &mut Session, p: GatewayParams) -> Result<Value, CmdError> {
    if p.acceptor.is_none()
        && p.donor.is_none()
        && p.acceptor_document.is_none()
        && p.donor_document.is_none()
    {
        return Err(CmdError::invalid_params(
            "give `acceptor` or `donor` as a file or an open document's id (the other is the \
             active document)",
        ));
    }
    let acceptor = load_source(
        s,
        p.acceptor.as_deref(),
        p.acceptor_document,
        "the acceptor",
    )?;
    let donor = load_source(s, p.donor.as_deref(), p.donor_document, "the donor")?;
    let reaction = match p.reaction {
        Some(Reaction::Bp) => GatewayReaction::Bp,
        _ => GatewayReaction::Lr,
    };
    let name = product_name(acceptor.name(), donor.name());
    let mut m = gateway(&acceptor, &donor, reaction, &name).map_err(|e| {
        CmdError::new(ErrorCode::InvalidInput, e).with_hint(match reaction {
            GatewayReaction::Lr => "LR joins an entry clone (attL sites) to a destination vector (attR sites); both need one att1 and one att2 site (cores TTTGTACAAAAAAGC and CTTTCTTGTACAAAG)",
            GatewayReaction::Bp => "BP joins an attB insert to a donor vector (attP sites); both need one att1 and one att2 site (cores TTTGTACAAAAAAGC and CTTTCTTGTACAAAG)",
        })
    })?;
    // A site whose arms are too unusual to name is let through, but said out loud: it is the
    // one case where an impossible reaction could still be simulated.
    let unnamed: Vec<String> = [(&acceptor, "acceptor"), (&donor, "donor")]
        .into_iter()
        .filter(|(m, _)| att_kinds(m).is_ok_and(|k| k.iter().any(Option::is_none)))
        .map(|(m, role)| {
            format!(
                "the {role} {}: its att sites aren't the usual attB, attL, attR or attP, so whether this reaction can happen wasn't checked",
                m.name()
            )
        })
        .collect();
    if p.preview {
        let sites: Vec<Value> = m
            .features()
            .iter()
            .filter(|f| f.name.starts_with("att"))
            .map(|f| json!({"name": f.name, "range": f.segments.first().map(|x| x.range)}))
            .collect();
        return Ok(json!({
            "name": m.name(),
            "length": m.len(),
            "sites": sites,
            "warnings": unnamed,
            "molecule": m,
        }));
    }
    let (operation, sites) = match reaction {
        GatewayReaction::Lr => ("gatewayLRCloning", ["attR", "attL"]),
        GatewayReaction::Bp => ("gatewayBPCloning", ["attP", "attB"]),
    };
    let names = |prefix: &str| {
        [
            ("name1", format!("{prefix}1")),
            ("name2", format!("{prefix}2")),
        ]
    };
    record(
        &mut m,
        operation,
        "",
        vec![HistoryNode::input(&acceptor), HistoryNode::input(&donor)],
        vec![
            summary("replace", &names(sites[0])),
            summary("insert", &names(sites[1])),
        ],
    )?;
    let mut out = open_product(s, m, p.discard_changes)?;
    if !unnamed.is_empty()
        && let Some(o) = out.as_object_mut()
    {
        o.insert("warnings".into(), json!(unnamed));
    }
    Ok(out)
}

/// How to assemble.
#[derive(Deserialize, JsonSchema, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
enum Method {
    Gibson,
    Infusion,
    GoldenGate,
}

/// One part of an assembly.
#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PartSpec {
    /// File (default: the active document).
    file: Option<String>,
    /// This part as an open document's id instead of a file (see `documents`).
    document: Option<u64>,
    /// Cut with these enzymes and keep the largest fragment (to linearise a vector).
    #[serde(default)]
    enzymes: Vec<String>,
    /// Use the PCR product of these two primers instead (a primer's name or a sequence).
    pcr: Option<PcrSpec>,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct PcrSpec {
    forward: String,
    reverse: String,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct AssembleParams {
    method: Method,
    parts: Vec<PartSpec>,
    /// Golden Gate: the Type IIS enzyme (default BsaI).
    enzyme: Option<String>,
    /// Homology assembly: shortest overlap (default 20 for Gibson, 15 for In-Fusion).
    min_overlap: Option<usize>,
    /// Close the construct into a circle (default true).
    circular: Option<bool>,
    #[serde(default)]
    preview: bool,
    #[serde(default)]
    discard_changes: bool,
}

/// The (shortest) PCR product of two primers on `m` (each a primer name or a sequence), with
/// the PCR recorded in its history.
pub(super) fn pcr_molecule(
    m: &Molecule,
    forward: &str,
    reverse: &str,
) -> Result<Molecule, CmdError> {
    let name = m.name().to_string();
    let (fname, f) = crate::gel::primer_residues(m, forward)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
    let (rname, r) = crate::gel::primer_residues(m, reverse)
        .map_err(|e| CmdError::new(ErrorCode::InvalidInput, e))?;
    let products = pcr(m, [&f, &r], None, DEFAULT_MAX_PRODUCT);
    let product = products.first().ok_or_else(|| {
        CmdError::new(
            ErrorCode::InvalidInput,
            format!("{fname} and {rname} make no product on {name}"),
        )
    })?;
    let mut pm = product_molecule(m, product, [&fname, &rname], &name)?;
    record(
        &mut pm,
        "amplifyFragment",
        &format!("primers {fname} and {rname}"),
        vec![HistoryNode::input(m)],
        vec![amplified(product, m.len())],
    )?;
    set_primers(&mut pm, &[(&fname, &f), (&rname, &r)], true)?;
    Ok(pm)
}

/// A part as a linear piece (for homology assembly).
fn piece(s: &Session, spec: &PartSpec) -> Result<(Part, String, HistoryNode), CmdError> {
    let m = load_source(s, spec.file.as_deref(), spec.document, "a part")?;
    let name = m.name().to_string();
    if let Some(p) = &spec.pcr {
        let pm = pcr_molecule(&m, &p.forward, &p.reverse)?;
        let label = format!(
            "PCR {} on {name}",
            pm.history()
                .map(|h| {
                    h.oligos
                        .iter()
                        .map(|(n, _)| n.as_str())
                        .collect::<Vec<_>>()
                        .join("/")
                })
                .unwrap_or_default()
        );
        return Ok((extract(&pm, 0, pm.len()), label, HistoryNode::input(&pm)));
    }
    if !spec.enzymes.is_empty() {
        let enz = enzymes(s, &m, &spec.enzymes)?;
        let f = digest_fragments(&m, &enz)
            .into_iter()
            .max_by_key(|f| f.len())
            .ok_or_else(|| {
                CmdError::new(
                    ErrorCode::InvalidInput,
                    format!("the enzymes don't cut {name}"),
                )
            })?;
        return Ok((
            f.part,
            format!("{name} cut with {}", spec.enzymes.join(" + ")),
            HistoryNode::input(&m),
        ));
    }
    if m.topology() == Topology::Circular {
        return Err(CmdError::new(
            ErrorCode::InvalidInput,
            format!("{name} is circular: name enzymes to linearise it, or a PCR"),
        ));
    }
    Ok((extract(&m, 0, m.len()), name, HistoryNode::input(&m)))
}

fn assemble_command(s: &mut Session, p: AssembleParams) -> Result<Value, CmdError> {
    if p.parts.is_empty() || p.parts.len() > 30 {
        return Err(CmdError::invalid_params("give 1 to 30 parts"));
    }
    let circular = p.circular.unwrap_or(true);
    let bad = |e: String| {
        CmdError::new(ErrorCode::InvalidInput, e)
            .with_hint("check the parts' order, their overlaps or overhangs, and the enzyme")
    };
    let mut inputs = Vec::new();
    let (assembly, labels) = match p.method {
        Method::Gibson | Method::Infusion => {
            let min = p
                .min_overlap
                .unwrap_or(if p.method == Method::Gibson { 20 } else { 15 });
            let mut pieces = Vec::new();
            let mut labels = Vec::new();
            for spec in &p.parts {
                let (part, label, node) = piece(s, spec)?;
                pieces.push(part);
                labels.push(label);
                inputs.push(node);
            }
            (
                homology_assembly(&pieces, min, circular).map_err(bad)?,
                labels,
            )
        }
        Method::GoldenGate => {
            let name = p.enzyme.clone().unwrap_or_else(|| "BsaI".to_string());
            let (_, enzyme) = s.enzymes().get(&name).ok_or_else(|| {
                CmdError::new(ErrorCode::NotFound, format!("unknown enzyme {name}"))
            })?;
            let enzyme = enzyme.clone();
            let mut molecules = Vec::new();
            for spec in &p.parts {
                molecules.push(load_source(
                    s,
                    spec.file.as_deref(),
                    spec.document,
                    "a part",
                )?);
            }
            let labels = molecules.iter().map(|m| m.name().to_string()).collect();
            inputs.extend(molecules.iter().map(HistoryNode::input));
            (golden_gate(&molecules, &enzyme).map_err(bad)?, labels)
        }
    };
    let junctions: Vec<Value> = assembly
        .junctions
        .iter()
        .map(|j| json!({"after": j.after, "length": j.length, "bases": j.bases, "tm": j.tm.map(|t| (t * 10.0).round() / 10.0)}))
        .collect();
    let labels: Vec<String> = labels;
    let name: String = labels.join(" + ").chars().take(200).collect();
    let mut m = to_molecule(&assembly, &name, circular).map_err(bad)?;
    if p.preview {
        return Ok(json!({
            "length": assembly.part.seq.len(),
            "parts": labels,
            "reversed": assembly.reversed,
            "junctions": junctions,
            "molecule": m,
        }));
    }
    let (operation, details) = match p.method {
        Method::Gibson => ("gibsonAssembly", labels.join("; ")),
        Method::Infusion => ("inFusionCloning", labels.join("; ")),
        Method::GoldenGate => (
            "goldenGateAssembly",
            format!("cut with {}", p.enzyme.as_deref().unwrap_or("BsaI")),
        ),
    };
    let summaries = inputs.iter().map(whole).collect();
    record(&mut m, operation, &details, inputs, summaries)?;
    let mut out = open_product(s, m, p.discard_changes)?;
    if let Some(o) = out.as_object_mut() {
        o.insert("junctions".into(), json!(junctions));
    }
    Ok(out)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct MutagenesisParams {
    /// The plasmid to copy, as a file (default: the active document).
    template: Option<String>,
    /// The plasmid as an open document's id instead of a file (see `documents`).
    template_document: Option<u64>,
    /// Forward primer (a primer of the template, or a sequence).
    forward: String,
    /// Reverse primer, likewise; it binds just upstream of the forward one.
    reverse: String,
    #[serde(default)]
    preview: bool,
    #[serde(default)]
    discard_changes: bool,
}

fn mutagenesis_command(s: &mut Session, p: MutagenesisParams) -> Result<Value, CmdError> {
    let template = load_source(
        s,
        p.template.as_deref(),
        p.template_document,
        "the template",
    )?;
    let bad = |e: String| CmdError::new(ErrorCode::InvalidInput, e);
    let (fname, f) = crate::gel::primer_residues(&template, &p.forward).map_err(bad)?;
    let (rname, r) = crate::gel::primer_residues(&template, &p.reverse).map_err(bad)?;
    let name: String = format!("{} (mutant)", template.name())
        .chars()
        .take(200)
        .collect();
    let mut m = mutagenesis(&template, &f, &r, [&fname, &rname], &name).map_err(|e| {
        bad(e).with_hint("the primers must bind back to back around the whole plasmid")
    })?;
    let change = m.len() as i64 - template.len() as i64;
    let differences = if change == 0 {
        m.seq()
            .as_bytes()
            .iter()
            .zip(template.seq().as_bytes())
            .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
            .count()
    } else {
        0
    };
    if p.preview {
        return Ok(json!({
            "length": m.len(),
            "length_change": change,
            "substitutions": differences,
            "molecule": m,
        }));
    }
    // The changed stretch of the template: between the longest common start and end.
    let (a, b) = (template.seq().as_bytes(), m.seq().as_bytes());
    let same = |x: &u8, y: &u8| x.eq_ignore_ascii_case(y);
    let start = a.iter().zip(b).take_while(|(x, y)| same(x, y)).count();
    let end = a
        .iter()
        .rev()
        .zip(b.iter().rev())
        .take_while(|(x, y)| same(x, y))
        .count()
        .min(a.len().saturating_sub(start));
    let last = a.len().saturating_sub(end + 1).max(start);
    record(
        &mut m,
        "primerDirectedMutagenesis",
        &format!("primers {fname} and {rname}"),
        vec![HistoryNode::input(&template)],
        vec![summary(
            "replace",
            &[("val1", start.to_string()), ("val2", last.to_string())],
        )],
    )?;
    set_primers(&mut m, &[(&fname, &f)], false)?;
    let mut out = open_product(s, m, p.discard_changes)?;
    if let Some(o) = out.as_object_mut() {
        o.insert("length_change".into(), json!(change));
        o.insert("substitutions".into(), json!(differences));
    }
    Ok(out)
}
