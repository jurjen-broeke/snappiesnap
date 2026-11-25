//! SBOL 2 ↔ document.
//!
//! The document shows one `ComponentDefinition`: the top of the design that has a sequence
//! (see [`SbolDocument::main_definition`]); a file with sequences but no definitions shows its
//! first sequence. Each `SequenceAnnotation` with `Range` locations becomes a feature: its
//! name from its title, else the title or displayId of the component it annotates, else its
//! displayId; its type from its Sequence Ontology role (or the component's), shown with the
//! INSDC name; `reverseComplement` ranges are on the reverse strand, `inline` ones forward,
//! others have no direction; two ranges meeting at the origin of a circular design are one
//! feature across it. `Cut` and `GenericLocation` annotations have no stretch to show and are
//! kept in the file.
//!
//! Saving splices only what changed into the original text (the sequence, the topology type,
//! the description, edited, removed and new annotations), so everything SnappieSnap doesn't
//! model (provenance, sub-components, interactions, other designs) is written back as it was.

use std::ops::Range as Span;

use snappiesnap_doc::{
    Alphabet, DocError, Feature, FeatureId, Molecule, NewFeature, Notes, Qualifier, Range, Segment,
    Seq, Strand, Topology,
};
use snappiesnap_formats::sbol::{
    self, ComponentDefinition, DCTERMS, DNA_REGION, INLINE, IUPAC_DNA, IUPAC_PROTEIN, Location,
    PROTEIN, REVERSE_COMPLEMENT, RNA_REGION, SBOL2, SO_CIRCULAR, SO_LINEAR, SbolDocument,
    escape_attr, escape_text,
};

use crate::genbank::join_across_origin;
use crate::{Exported, Format, Imported, IoError, Source};

#[derive(Debug, Clone, PartialEq)]
pub struct SbolSource {
    original: Vec<u8>,
    doc: SbolDocument,
    /// The shown definition's index, if the document shows one.
    definition: Option<usize>,
    /// The shown sequence's index.
    sequence: usize,
    imported: Molecule,
    /// Document feature id → annotation index in the shown definition.
    feature_map: Vec<(FeatureId, usize)>,
}

pub(crate) fn looks_like_sbol(bytes: &[u8]) -> bool {
    sbol::looks_like_sbol(bytes)
}

/// Sequence Ontology roles and the feature types they're shown as (INSDC names where there is
/// one). Roles are matched by their `SO:nnnnnnn` accession, whatever the URI prefix.
const ROLES: &[(&str, &str)] = &[
    ("SO:0000167", "promoter"),
    ("SO:0000139", "RBS"),
    ("SO:0000316", "CDS"),
    ("SO:0000141", "terminator"),
    ("SO:0000296", "rep_origin"),
    ("SO:0000057", "protein_bind"),
    ("SO:0000704", "gene"),
    ("SO:0000234", "mRNA"),
    ("SO:0000165", "enhancer"),
    ("SO:0000551", "polyA_signal"),
    ("SO:0000553", "polyA_site"),
    ("SO:0005850", "primer_bind"),
    ("SO:0000409", "misc_binding"),
    ("SO:0000418", "sig_peptide"),
    ("SO:0000188", "intron"),
    ("SO:0000147", "exon"),
    ("SO:0000204", "5'UTR"),
    ("SO:0000205", "3'UTR"),
    ("SO:0000252", "rRNA"),
    ("SO:0000253", "tRNA"),
    ("SO:0000673", "misc_RNA"),
    ("SO:0005836", "regulatory"),
    ("SO:0000657", "repeat_region"),
    ("SO:0001037", "mobile_element"),
    ("SO:0000413", "misc_difference"),
    ("SO:0000298", "misc_recomb"),
    ("SO:0000313", "stem_loop"),
    ("SO:0000419", "mat_peptide"),
    ("SO:0000725", "transit_peptide"),
    ("SO:0000286", "LTR"),
    ("SO:0001060", "variation"),
    ("SO:0000110", "misc_feature"),
];

fn accession(role: &str) -> Option<&str> {
    let at = role.rfind("SO:")?;
    role.get(at..)
}

fn kind_for(roles: &[String]) -> Option<String> {
    roles.iter().find_map(|r| {
        let acc = accession(r)?;
        ROLES
            .iter()
            .find(|(so, _)| *so == acc)
            .map(|(_, kind)| (*kind).to_string())
    })
}

fn role_for(kind: &str) -> String {
    let acc = ROLES
        .iter()
        .find(|(_, k)| *k == kind)
        .map_or("SO:0000110", |(so, _)| so);
    format!("http://identifiers.org/so/{acc}")
}

/// The alphabet of a sequence, from its encoding and the definition's types.
fn alphabet(encoding: Option<&str>, types: &[String]) -> Alphabet {
    if encoding == Some(IUPAC_PROTEIN) || types.iter().any(|t| t == PROTEIN) {
        Alphabet::Protein
    } else if types.iter().any(|t| t == RNA_REGION) {
        Alphabet::Rna
    } else {
        Alphabet::Dna
    }
}

/// A shown definition's annotations as document features (with their annotation index).
fn features(
    doc: &SbolDocument,
    cd: &ComponentDefinition,
    len: usize,
    topology: Topology,
    warnings: &mut Vec<String>,
) -> Vec<(usize, NewFeature)> {
    let mut out = Vec::new();
    for (index, a) in cd.annotations.iter().enumerate() {
        let label = a.display_id.clone().unwrap_or_else(|| a.uri.clone());
        let ranges: Vec<(u64, u64, Option<&str>)> = a
            .locations
            .iter()
            .filter_map(|l| match l {
                Location::Range {
                    start,
                    end,
                    orientation,
                } => Some((*start, *end, orientation.as_deref())),
                _ => None,
            })
            .collect();
        if ranges.is_empty() {
            warnings.push(format!(
                "annotation {label} marks a cut site or has no position; it is kept in the file"
            ));
            continue;
        }
        if ranges.len() < a.locations.len() {
            warnings.push(format!(
                "annotation {label}: its cut-site or unpositioned parts are kept in the file only"
            ));
        }
        let mut segments = Vec::new();
        let mut bad = None;
        for (start, end, _) in &ranges {
            if *start == 0 || end < start || *end > len as u64 {
                bad = Some(format!("{start}..{end} is outside the {len} bp sequence"));
                break;
            }
            segments.push(Range::new(*start as usize - 1, *end as usize));
        }
        if let Some(reason) = bad {
            warnings.push(format!(
                "annotation {label} not shown ({reason}); it is kept in the file"
            ));
            continue;
        }
        let orientations: Vec<Option<&str>> = ranges.iter().map(|(_, _, o)| *o).collect();
        let strand = if orientations.iter().all(|o| *o == Some(REVERSE_COMPLEMENT)) {
            Strand::Reverse
        } else if orientations.iter().all(|o| *o == Some(INLINE)) {
            Strand::Forward
        } else {
            Strand::None
        };
        segments.sort_by_key(|r| r.start);
        if topology == Topology::Circular {
            segments = join_across_origin(segments, len);
        }
        // Biological (5'→3') order: descending on the reverse strand.
        if strand == Strand::Reverse {
            segments.reverse();
        }
        let component_definition = a
            .component
            .as_deref()
            .and_then(|c| cd.components.iter().find(|x| x.uri == c))
            .and_then(|c| c.definition.as_deref())
            .and_then(|d| doc.definition(d));
        let roles = if a.roles.is_empty() {
            component_definition
                .map(|d| d.roles.clone())
                .unwrap_or_default()
        } else {
            a.roles.clone()
        };
        let kind = kind_for(&roles).unwrap_or_else(|| "misc_feature".into());
        let name = a
            .title
            .clone()
            .or_else(|| {
                component_definition.and_then(|d| d.title.clone().or_else(|| d.display_id.clone()))
            })
            .or_else(|| a.display_id.clone())
            .unwrap_or_else(|| kind.clone());
        let qualifiers = a
            .description
            .iter()
            .map(|d| Qualifier {
                name: "note".into(),
                values: vec![d.clone()],
            })
            .collect();
        out.push((
            index,
            NewFeature {
                color: None,
                name,
                kind,
                strand,
                segments: segments
                    .into_iter()
                    .map(|range| Segment { range })
                    .collect(),
                qualifiers,
            },
        ));
    }
    out
}

pub(crate) fn import(name: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    let doc = SbolDocument::parse(bytes).map_err(|e| IoError::Format(e.to_string()))?;
    let mut warnings = Vec::new();
    let definition = doc.main_definition().and_then(|cd| {
        doc.component_definitions
            .iter()
            .position(|c| c.uri == cd.uri)
    });
    let cd = definition.and_then(|i| doc.component_definitions.get(i));
    let sequence = match cd {
        Some(cd) => cd.sequences.iter().find_map(|s| {
            doc.sequences
                .iter()
                .position(|x| Some(x) == doc.sequence(s))
        }),
        None => {
            if !doc.sequences.is_empty() {
                warnings.push("no design uses the sequence; it is shown without features".into());
            }
            (!doc.sequences.is_empty()).then_some(0)
        }
    }
    .ok_or_else(|| {
        IoError::Format(format!(
            "this SBOL file has no sequence to show ({} designs, none with a sequence)",
            doc.component_definitions.len()
        ))
    })?;
    let seq_obj = doc
        .sequences
        .get(sequence)
        .ok_or_else(|| IoError::Format("missing sequence".into()))?;
    let residues: Vec<u8> = seq_obj
        .elements
        .bytes()
        .filter(|b| !b.is_ascii_whitespace())
        .map(|b| b.to_ascii_uppercase())
        .collect();
    // Without a design to give the type, the residues tell (as for FASTA).
    let alphabet = match cd {
        Some(c) => alphabet(seq_obj.encoding.as_deref(), &c.types),
        None if seq_obj.encoding.as_deref() == Some(IUPAC_PROTEIN) => Alphabet::Protein,
        // The nucleic acid encoding: RNA when it has U and no T, else DNA.
        None if seq_obj.encoding.as_deref() == Some(IUPAC_DNA) => {
            if residues.contains(&b'U') && !residues.contains(&b'T') {
                Alphabet::Rna
            } else {
                Alphabet::Dna
            }
        }
        None => crate::fasta::guess_alphabet(&residues),
    };
    let seq = Seq::new(alphabet, residues).map_err(DocError::from)?;
    let len = seq.len();
    let topology = if cd.is_some_and(ComponentDefinition::circular) {
        Topology::Circular
    } else {
        Topology::Linear
    };
    let mut m = Molecule::new(name, seq, topology)?;
    if let Some(description) = cd.and_then(|c| c.description.clone()) {
        let notes = Notes {
            description,
            ..Notes::default()
        };
        if let Err(e) = m.set_notes(notes) {
            warnings.push(format!("description not shown ({e})"));
        }
    }
    let mut feature_map = Vec::new();
    if let Some(cd) = cd {
        for (index, nf) in features(&doc, cd, len, topology, &mut warnings) {
            match m.add_feature(nf) {
                Ok(id) => feature_map.push((id, index)),
                Err(e) => warnings.push(format!(
                    "an annotation isn't shown ({e}); it is kept in the file"
                )),
            }
        }
    }
    let others = doc.component_definitions.len().saturating_sub(1);
    if others > 0 {
        warnings.push(format!(
            "{others} more design(s) in the file are kept but not shown"
        ));
    }
    Ok(Imported {
        source: Source::Sbol(Box::new(SbolSource {
            original: bytes.to_vec(),
            doc,
            definition,
            sequence,
            imported: m.clone(),
            feature_map,
        })),
        molecule: m,
        format: Format::Sbol,
        warnings,
    })
}

/// Qualified names in the document's own prefixes (or with an inline declaration where the
/// document element doesn't declare the namespace).
struct Names {
    sbol: String,
    rdf: String,
    dcterms: String,
    /// Extra declarations for the outermost new element.
    declare: String,
}

impl Names {
    fn for_document(doc: Option<&SbolDocument>) -> Names {
        let mut declare = String::new();
        let mut pick = |ns: &str, fallback: &str| match doc.and_then(|d| d.prefix(ns)) {
            Some(p) if !p.is_empty() => p.to_string(),
            _ => {
                if doc.is_some() {
                    declare.push_str(&format!(" xmlns:{fallback}=\"{ns}\""));
                }
                fallback.to_string()
            }
        };
        let sbol = pick(SBOL2, "sbol");
        let rdf = pick(sbol::RDF, "rdf");
        let dcterms = pick(DCTERMS, "dcterms");
        Names {
            sbol,
            rdf,
            dcterms,
            declare,
        }
    }
}

/// A displayId: letters, digits and underscores, not starting with a digit.
fn display_id(name: &str) -> String {
    let mut s: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if s.is_empty() || s.starts_with(|c: char| c.is_ascii_digit()) {
        s.insert(0, '_');
    }
    s
}

/// One annotation as an `sbol:sequenceAnnotation` property, indented by `indent`.
#[allow(clippy::too_many_arguments)]
fn annotation_xml(
    n: &Names,
    f: &Feature,
    len: usize,
    uri: &str,
    roles: &[String],
    component: Option<&str>,
    indent: &str,
    declare: &str,
) -> String {
    let (s, r, d) = (&n.sbol, &n.rdf, &n.dcterms);
    let step = "  ";
    let i1 = format!("{indent}{step}");
    let i2 = format!("{i1}{step}");
    let i3 = format!("{i2}{step}");
    let i4 = format!("{i3}{step}");
    let id = uri
        .rsplit('/')
        .next()
        .map_or_else(|| display_id(&f.name), display_id);
    let mut out = format!("{indent}<{s}:sequenceAnnotation{declare}>\n");
    out.push_str(&format!(
        "{i1}<{s}:SequenceAnnotation {r}:about=\"{}\">\n",
        escape_attr(uri)
    ));
    out.push_str(&format!("{i2}<{s}:displayId>{id}</{s}:displayId>\n"));
    if !f.name.is_empty() {
        out.push_str(&format!(
            "{i2}<{d}:title>{}</{d}:title>\n",
            escape_text(&f.name)
        ));
    }
    if let Some(note) = f
        .qualifiers
        .iter()
        .find(|q| q.name == "note")
        .and_then(|q| q.values.first())
    {
        out.push_str(&format!(
            "{i2}<{d}:description>{}</{d}:description>\n",
            escape_text(note)
        ));
    }
    // Each part as a Range (a part across the origin as two), ascending.
    let mut pieces: Vec<(usize, usize)> = Vec::new();
    for seg in &f.segments {
        let (first, second) = seg.range.segments(len);
        pieces.extend(
            std::iter::once(first)
                .chain(second)
                .map(|p| (p.start + 1, p.end)),
        );
    }
    pieces.sort_unstable();
    let orientation = match f.strand {
        Strand::Forward => Some(INLINE),
        Strand::Reverse => Some(REVERSE_COMPLEMENT),
        _ => None,
    };
    for (k, (start, end)) in pieces.iter().enumerate() {
        let rid = format!("range{}", k + 1);
        out.push_str(&format!("{i2}<{s}:location>\n"));
        out.push_str(&format!(
            "{i3}<{s}:Range {r}:about=\"{}/{rid}\">\n",
            escape_attr(uri)
        ));
        out.push_str(&format!("{i4}<{s}:displayId>{rid}</{s}:displayId>\n"));
        out.push_str(&format!("{i4}<{s}:start>{start}</{s}:start>\n"));
        out.push_str(&format!("{i4}<{s}:end>{end}</{s}:end>\n"));
        if let Some(o) = orientation {
            out.push_str(&format!("{i4}<{s}:orientation {r}:resource=\"{o}\"/>\n"));
        }
        out.push_str(&format!("{i3}</{s}:Range>\n"));
        out.push_str(&format!("{i2}</{s}:location>\n"));
    }
    for role in roles {
        out.push_str(&format!(
            "{i2}<{s}:role {r}:resource=\"{}\"/>\n",
            escape_attr(role)
        ));
    }
    if let Some(c) = component {
        out.push_str(&format!(
            "{i2}<{s}:component {r}:resource=\"{}\"/>\n",
            escape_attr(c)
        ));
    }
    out.push_str(&format!("{i1}</{s}:SequenceAnnotation>\n"));
    out.push_str(&format!("{indent}</{s}:sequenceAnnotation>\n"));
    out
}

/// The whitespace before `at` on its line.
fn indent_at(text: &str, at: usize) -> String {
    let line_start = text
        .get(..at)
        .and_then(|t| t.rfind('\n'))
        .map_or(0, |i| i + 1);
    text.get(line_start..at)
        .unwrap_or("")
        .chars()
        .take_while(|c| *c == ' ' || *c == '\t')
        .collect()
}

/// A span widened to whole lines when it sits alone on its lines, so removing it leaves no
/// blank line.
fn whole_lines(text: &str, span: &Span<usize>) -> Span<usize> {
    let line_start = text
        .get(..span.start)
        .and_then(|t| t.rfind('\n'))
        .map_or(0, |i| i + 1);
    let before_blank = text
        .get(line_start..span.start)
        .is_some_and(|t| t.chars().all(|c| c == ' ' || c == '\t'));
    let after_newline = text.get(span.end..).is_some_and(|t| t.starts_with('\n'));
    if before_blank && after_newline {
        line_start..span.end + 1
    } else {
        span.clone()
    }
}

pub(crate) fn export(m: &Molecule, source: Option<&SbolSource>) -> Result<Exported, IoError> {
    let mut warnings = Vec::new();
    if !m.primers().is_empty() {
        warnings.push("primers aren't saved in SBOL files".to_string());
    }
    let Some(src) = source else {
        let bytes = new_document(m).into_bytes();
        return finish(m, bytes, warnings);
    };
    if *m == src.imported {
        return Ok(Exported {
            bytes: src.original.clone(),
            source: Source::Sbol(Box::new(src.clone())),
            warnings,
        });
    }
    let text = std::str::from_utf8(&src.original)
        .map_err(|_| IoError::Format("the original file isn't text".into()))?;
    let Some(cd) = src
        .definition
        .and_then(|i| src.doc.component_definitions.get(i))
    else {
        // A bare sequence without a design: write a new document for it.
        warnings.push(
            "the file had no design for this sequence; it is saved as a new SBOL document".into(),
        );
        let bytes = new_document(m).into_bytes();
        return finish(m, bytes, warnings);
    };
    let names = Names::for_document(Some(&src.doc));
    let len = m.len();
    let mut edits: Vec<(Span<usize>, String)> = Vec::new();

    // The sequence.
    if m.seq() != src.imported.seq()
        && let Some(span) = src
            .doc
            .sequences
            .get(src.sequence)
            .and_then(|s| s.elements_span.clone())
    {
        let residues = String::from_utf8_lossy(m.seq().as_bytes());
        let residues = if m.seq().alphabet() == Alphabet::Protein {
            residues.into_owned()
        } else {
            residues.to_ascii_lowercase()
        };
        edits.push((span, residues));
    }
    let child_indent = cd
        .title_span
        .as_ref()
        .map(|s| indent_at(text, s.start))
        .unwrap_or_else(|| format!("{}  ", indent_at(text, cd.span.start)));
    // New properties go on their own lines before the end tag's line (or right before the end
    // tag when it shares a line with other content).
    let line_start = text
        .get(..cd.inner.end)
        .and_then(|t| t.rfind('\n'))
        .map_or(0, |i| i + 1);
    let insert_at = if text
        .get(line_start..cd.inner.end)
        .is_some_and(|t| t.chars().all(|c| c == ' ' || c == '\t'))
    {
        line_start
    } else {
        cd.inner.end
    };
    let mut inserts = String::new();

    // Topology.
    if m.topology() != src.imported.topology() {
        let uri = if m.topology() == Topology::Circular {
            SO_CIRCULAR
        } else {
            SO_LINEAR
        };
        let element = format!("<{}:type {}:resource=\"{uri}\"/>", names.sbol, names.rdf);
        match &cd.topology_span {
            Some(span) => edits.push((span.clone(), element)),
            None => inserts.push_str(&format!("{child_indent}{element}\n")),
        }
    }
    // Description.
    if m.notes().description != src.imported.notes().description {
        let d = &m.notes().description;
        let element = format!(
            "<{0}:description>{1}</{0}:description>",
            names.dcterms,
            escape_text(d)
        );
        match (&cd.description_span, d.is_empty()) {
            (Some(span), true) => edits.push((whole_lines(text, span), String::new())),
            (Some(span), false) => edits.push((span.clone(), element)),
            (None, false) => {
                inserts.push_str(&format!("{child_indent}{element}{}\n", names.declare))
            }
            (None, true) => {}
        }
    }
    // Annotations.
    for (id, index) in &src.feature_map {
        let Some(a) = cd.annotations.get(*index) else {
            continue;
        };
        let then = src.imported.feature(*id);
        match (m.feature(*id), then) {
            (None, _) => edits.push((whole_lines(text, &a.span), String::new())),
            (Some(now), Some(then)) if now != then => {
                let roles = if now.kind == then.kind && !a.roles.is_empty() {
                    a.roles.clone()
                } else if now.kind == then.kind {
                    Vec::new()
                } else {
                    vec![role_for(&now.kind)]
                };
                let indent = indent_at(text, a.span.start);
                let xml = annotation_xml(
                    &names,
                    now,
                    len,
                    &a.uri,
                    &roles,
                    a.component.as_deref(),
                    &indent,
                    &names.declare,
                );
                edits.push((whole_lines(text, &a.span), xml));
            }
            _ => {}
        }
    }
    let mapped: Vec<FeatureId> = src.feature_map.iter().map(|(id, _)| *id).collect();
    let base = cd
        .persistent_identity
        .clone()
        .unwrap_or_else(|| cd.uri.clone());
    let taken: Vec<&str> = cd.annotations.iter().map(|a| a.uri.as_str()).collect();
    for f in m.features().iter().filter(|f| !mapped.contains(&f.id)) {
        let mut uri = format!("{base}/{}", display_id(&f.name));
        let mut k = 2;
        while taken.contains(&uri.as_str()) {
            uri = format!("{base}/{}_{k}", display_id(&f.name));
            k += 1;
        }
        inserts.push_str(&annotation_xml(
            &names,
            f,
            len,
            &uri,
            &[role_for(&f.kind)],
            None,
            &child_indent,
            &names.declare,
        ));
    }
    if !inserts.is_empty() {
        edits.push((insert_at..insert_at, inserts));
    }
    // Apply from the end, so earlier spans stay valid.
    edits.sort_by_key(|(span, _)| std::cmp::Reverse((span.start, span.end)));
    let mut out = text.to_string();
    let mut limit = out.len();
    for (span, replacement) in edits {
        if span.end > limit || span.start > span.end {
            return Err(IoError::Format("overlapping edits in the SBOL file".into()));
        }
        out.replace_range(span.clone(), &replacement);
        limit = span.start;
    }
    finish(m, out.into_bytes(), warnings)
}

/// Read back what was written: the next save starts from it.
fn finish(m: &Molecule, bytes: Vec<u8>, warnings: Vec<String>) -> Result<Exported, IoError> {
    let again = import(m.name(), &bytes)?;
    let mut source = match again.source {
        Source::Sbol(s) => s,
        _ => {
            return Err(IoError::Format(
                "the written SBOL file didn't read back".into(),
            ));
        }
    };
    // The written file shows the same features: map them to the document's ids by order.
    let ids: Vec<FeatureId> = m.features().iter().map(|f| f.id).collect();
    let written: Vec<FeatureId> = again.molecule.features().iter().map(|f| f.id).collect();
    if ids.len() == written.len() {
        source.feature_map = source
            .feature_map
            .iter()
            .filter_map(|(wid, index)| {
                let pos = written.iter().position(|w| w == wid)?;
                Some((*ids.get(pos)?, *index))
            })
            .collect();
        source.imported = m.clone();
    }
    Ok(Exported {
        bytes,
        source: Source::Sbol(source),
        warnings,
    })
}

/// A new SBOL 2 document: one ComponentDefinition with its annotations, and its Sequence.
/// URIs use the `urn:snappiesnap:` prefix.
fn new_document(m: &Molecule) -> String {
    let names = Names::for_document(None);
    let id = display_id(m.name());
    let base = format!("urn:snappiesnap:{id}");
    let len = m.len();
    let region = match m.seq().alphabet() {
        Alphabet::Protein => PROTEIN,
        Alphabet::Rna => RNA_REGION,
        Alphabet::Dna => DNA_REGION,
    };
    let encoding = if m.seq().alphabet() == Alphabet::Protein {
        IUPAC_PROTEIN
    } else {
        IUPAC_DNA
    };
    let topology = if m.topology() == Topology::Circular {
        SO_CIRCULAR
    } else {
        SO_LINEAR
    };
    let mut out = String::from("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    out.push_str(&format!(
        "<rdf:RDF xmlns:rdf=\"{}\" xmlns:sbol=\"{SBOL2}\" xmlns:dcterms=\"{DCTERMS}\" xmlns:prov=\"{}\">\n",
        sbol::RDF,
        sbol::PROV
    ));
    out.push_str(&format!(
        "  <sbol:ComponentDefinition rdf:about=\"{}\">\n",
        escape_attr(&base)
    ));
    out.push_str(&format!("    <sbol:displayId>{id}</sbol:displayId>\n"));
    out.push_str(&format!(
        "    <dcterms:title>{}</dcterms:title>\n",
        escape_text(m.name())
    ));
    if !m.notes().description.is_empty() {
        out.push_str(&format!(
            "    <dcterms:description>{}</dcterms:description>\n",
            escape_text(&m.notes().description)
        ));
    }
    out.push_str(&format!("    <sbol:type rdf:resource=\"{region}\"/>\n"));
    if m.seq().alphabet() != Alphabet::Protein {
        out.push_str(&format!("    <sbol:type rdf:resource=\"{topology}\"/>\n"));
    }
    out.push_str("    <sbol:role rdf:resource=\"http://identifiers.org/so/SO:0000804\"/>\n");
    let mut used: Vec<String> = Vec::new();
    for f in m.features() {
        let mut fid = display_id(&f.name);
        let mut k = 2;
        while used.contains(&fid) {
            fid = format!("{}_{k}", display_id(&f.name));
            k += 1;
        }
        used.push(fid.clone());
        out.push_str(&annotation_xml(
            &names,
            f,
            len,
            &format!("{base}/{fid}"),
            &[role_for(&f.kind)],
            None,
            "    ",
            "",
        ));
    }
    out.push_str(&format!(
        "    <sbol:sequence rdf:resource=\"{}/sequence\"/>\n",
        escape_attr(&base)
    ));
    out.push_str("  </sbol:ComponentDefinition>\n");
    out.push_str(&format!(
        "  <sbol:Sequence rdf:about=\"{}/sequence\">\n",
        escape_attr(&base)
    ));
    out.push_str("    <sbol:displayId>sequence</sbol:displayId>\n");
    let residues = String::from_utf8_lossy(m.seq().as_bytes());
    let residues = if m.seq().alphabet() == Alphabet::Protein {
        residues.into_owned()
    } else {
        residues.to_ascii_lowercase()
    };
    out.push_str(&format!("    <sbol:elements>{residues}</sbol:elements>\n"));
    out.push_str(&format!(
        "    <sbol:encoding rdf:resource=\"{encoding}\"/>\n"
    ));
    out.push_str("  </sbol:Sequence>\n");
    out.push_str("</rdf:RDF>\n");
    out
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::panic)]
    use super::*;

    const FILE: &str = r#"<?xml version="1.0" ?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:sbol="http://sbols.org/v2#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:prov="http://www.w3.org/ns/prov#">
  <sbol:ComponentDefinition rdf:about="http://x.org/part">
    <sbol:displayId>part</sbol:displayId>
    <prov:wasDerivedFrom rdf:resource="http://somewhere.org/original"/>
    <dcterms:title>part</dcterms:title>
    <sbol:type rdf:resource="http://www.biopax.org/release/biopax-level3.owl#DnaRegion"/>
    <sbol:type rdf:resource="http://identifiers.org/so/SO:0000988"/>
    <sbol:sequenceAnnotation>
      <sbol:SequenceAnnotation rdf:about="http://x.org/part/prom">
        <sbol:displayId>prom</sbol:displayId>
        <dcterms:title>my promoter</dcterms:title>
        <sbol:location>
          <sbol:Range rdf:about="http://x.org/part/prom/r">
            <sbol:start>2</sbol:start>
            <sbol:end>7</sbol:end>
            <sbol:orientation rdf:resource="http://sbols.org/v2#inline"/>
          </sbol:Range>
        </sbol:location>
        <sbol:role rdf:resource="http://identifiers.org/so/SO:0000167"/>
      </sbol:SequenceAnnotation>
    </sbol:sequenceAnnotation>
    <sbol:sequenceAnnotation>
      <sbol:SequenceAnnotation rdf:about="http://x.org/part/wrap">
        <sbol:displayId>wrap</sbol:displayId>
        <sbol:location>
          <sbol:Range rdf:about="http://x.org/part/wrap/a">
            <sbol:start>18</sbol:start>
            <sbol:end>20</sbol:end>
            <sbol:orientation rdf:resource="http://sbols.org/v2#reverseComplement"/>
          </sbol:Range>
        </sbol:location>
        <sbol:location>
          <sbol:Range rdf:about="http://x.org/part/wrap/b">
            <sbol:start>1</sbol:start>
            <sbol:end>3</sbol:end>
            <sbol:orientation rdf:resource="http://sbols.org/v2#reverseComplement"/>
          </sbol:Range>
        </sbol:location>
        <sbol:role rdf:resource="http://identifiers.org/so/SO:0000316"/>
      </sbol:SequenceAnnotation>
    </sbol:sequenceAnnotation>
    <sbol:sequenceAnnotation>
      <sbol:SequenceAnnotation rdf:about="http://x.org/part/cut">
        <sbol:location>
          <sbol:Cut rdf:about="http://x.org/part/cut/c">
            <sbol:at>10</sbol:at>
          </sbol:Cut>
        </sbol:location>
      </sbol:SequenceAnnotation>
    </sbol:sequenceAnnotation>
    <sbol:sequence rdf:resource="http://x.org/seq"/>
  </sbol:ComponentDefinition>
  <sbol:Sequence rdf:about="http://x.org/seq">
    <sbol:elements>acgtacgtacgtacgtacgt</sbol:elements>
    <sbol:encoding rdf:resource="http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"/>
  </sbol:Sequence>
</rdf:RDF>
"#;

    fn source(i: &Imported) -> SbolSource {
        match &i.source {
            Source::Sbol(s) => s.as_ref().clone(),
            _ => panic!("not SBOL"),
        }
    }

    #[test]
    fn opens_ranges_strands_and_the_origin() {
        let i = import("part", FILE.as_bytes()).unwrap();
        let m = &i.molecule;
        assert_eq!((m.len(), m.topology()), (20, Topology::Circular));
        let p = m
            .features()
            .iter()
            .find(|f| f.name == "my promoter")
            .unwrap();
        assert_eq!((p.kind.as_str(), p.strand), ("promoter", Strand::Forward));
        assert_eq!(p.segments[0].range, Range::new(1, 7));
        let w = m.features().iter().find(|f| f.name == "wrap").unwrap();
        assert_eq!((w.kind.as_str(), w.strand), ("CDS", Strand::Reverse));
        assert_eq!(
            w.segments.iter().map(|s| s.range).collect::<Vec<_>>(),
            vec![Range::new(17, 3)]
        );
        // The cut site isn't a stretch: kept, reported.
        assert_eq!(m.features().len(), 2);
        assert!(i.warnings.iter().any(|w| w.contains("cut site")));
    }

    #[test]
    fn saves_unchanged_bytes_and_splices_edits() {
        let i = import("part", FILE.as_bytes()).unwrap();
        let src = source(&i);
        assert_eq!(
            export(&i.molecule, Some(&src)).unwrap().bytes,
            FILE.as_bytes()
        );
        let mut m = i.molecule.clone();
        let mut p = m
            .features()
            .iter()
            .find(|f| f.name == "my promoter")
            .unwrap()
            .clone();
        p.segments[0].range = Range::new(2, 8);
        m.update_feature(p).unwrap();
        let w = m.features().iter().find(|f| f.name == "wrap").unwrap().id;
        m.remove_feature(w).unwrap();
        m.add_feature(NewFeature {
            color: None,
            name: "new tag".into(),
            kind: "CDS".into(),
            strand: Strand::Forward,
            segments: vec![Segment {
                range: Range::new(10, 15),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        m.set_topology(Topology::Linear).unwrap();
        let out = export(&m, Some(&src)).unwrap();
        let text = String::from_utf8(out.bytes.clone()).unwrap();
        // Untouched content stays, byte for byte.
        assert!(text.contains(
            "    <prov:wasDerivedFrom rdf:resource=\"http://somewhere.org/original\"/>\n"
        ));
        assert!(text.contains("<sbol:Cut rdf:about=\"http://x.org/part/cut/c\">"));
        assert!(text.contains("SO:0000987") && !text.contains("SO:0000988"));
        assert!(!text.contains("http://x.org/part/wrap"));
        assert!(text.contains("<sbol:start>3</sbol:start>"));
        assert!(text.contains("rdf:about=\"http://x.org/part/new_tag\""));
        let again = import("part", text.as_bytes()).unwrap();
        assert_eq!(again.molecule.features().len(), 2);
        assert_eq!(again.molecule.topology(), Topology::Linear);
        // Saving again from the written source changes nothing.
        assert_eq!(
            export(
                &m,
                Some(match &out.source {
                    Source::Sbol(s) => s,
                    _ => panic!(),
                })
            )
            .unwrap()
            .bytes,
            out.bytes
        );
    }

    #[test]
    fn new_documents_read_back() {
        let i = import("part", FILE.as_bytes()).unwrap();
        let out = export(&i.molecule, None).unwrap();
        let again = import("part", &out.bytes).unwrap();
        assert_eq!(again.molecule.seq(), i.molecule.seq());
        assert_eq!(again.molecule.topology(), Topology::Circular);
        let names = |m: &Molecule| {
            let mut v: Vec<(String, String, Vec<Range>, Strand)> = m
                .features()
                .iter()
                .map(|f| {
                    (
                        f.name.clone(),
                        f.kind.clone(),
                        f.segments.iter().map(|s| s.range).collect(),
                        f.strand,
                    )
                })
                .collect();
            v.sort_by(|a, b| a.0.cmp(&b.0));
            v
        };
        assert_eq!(names(&again.molecule), names(&i.molecule));
    }
}
