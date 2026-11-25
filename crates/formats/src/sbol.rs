//! SBOL 2 (Synthetic Biology Open Language, version 2.3: Madsen et al., J. Integr. Bioinform.
//! 16(2), 2019) in its RDF/XML serialisation.
//!
//! The file is read into an element tree that keeps every element's byte span, with namespace
//! prefixes resolved, so writers can splice changed pieces into the original text and keep
//! everything else (provenance, sub-components, models, other designs) byte for byte. On top
//! of the tree, [`SbolDocument`] gives the parts a sequence editor uses: `ComponentDefinition`s
//! with their types, roles, sequence annotations (`Range`, `Cut` and `GenericLocation`
//! locations) and `Sequence`s.
//!
//! Only the RDF/XML forms that SBOL writers produce are understood: node elements with
//! `rdf:about`, property elements holding text, an `rdf:resource` reference or one nested node
//! element. Anything else is kept in the tree untouched.

use std::ops::Range;

use crate::xml::{Token, Tokenizer, decode};

pub const RDF: &str = "http://www.w3.org/1999/02/22-rdf-syntax-ns#";
pub const SBOL2: &str = "http://sbols.org/v2#";
pub const SBOL3: &str = "http://sbols.org/v3#";
pub const SBOL1: &str = "http://sbols.org/v1#";
pub const DCTERMS: &str = "http://purl.org/dc/terms/";
pub const PROV: &str = "http://www.w3.org/ns/prov#";
/// Topology types (Sequence Ontology).
pub const SO_LINEAR: &str = "http://identifiers.org/so/SO:0000987";
pub const SO_CIRCULAR: &str = "http://identifiers.org/so/SO:0000988";
/// Orientations.
pub const INLINE: &str = "http://sbols.org/v2#inline";
pub const REVERSE_COMPLEMENT: &str = "http://sbols.org/v2#reverseComplement";
/// Sequence encodings.
pub const IUPAC_DNA: &str = "http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html";
pub const IUPAC_PROTEIN: &str = "http://www.chem.qmul.ac.uk/iupac/AminoAcid/";
/// Component types (BioPAX).
pub const DNA_REGION: &str = "http://www.biopax.org/release/biopax-level3.owl#DnaRegion";
pub const RNA_REGION: &str = "http://www.biopax.org/release/biopax-level3.owl#RnaRegion";
pub const PROTEIN: &str = "http://www.biopax.org/release/biopax-level3.owl#Protein";

/// Deepest element nesting accepted.
const MAX_DEPTH: usize = 256;
/// Most elements accepted (the tree's size grows with them).
const MAX_ELEMENTS: usize = 2_000_000;

/// Why an SBOL file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SbolError {
    #[error("the file isn't valid UTF-8 text")]
    NotText,

    #[error("not an SBOL 2 file: {0}")]
    NotSbol2(String),

    #[error("malformed XML at byte {offset}: {reason}")]
    Xml { offset: usize, reason: String },
}

/// An attribute, namespace resolved and value decoded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attr {
    pub ns: String,
    pub local: String,
    pub value: String,
}

/// An XML element with its namespace resolved and its spans in the source text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Element {
    pub ns: String,
    pub local: String,
    pub attrs: Vec<Attr>,
    /// Namespace prefixes declared on this element: (prefix, URI).
    pub namespaces: Vec<(String, String)>,
    pub children: Vec<Element>,
    /// The element's own text (outside child elements), decoded.
    pub text: String,
    /// From `<` of the start tag to `>` of the end tag.
    pub span: Range<usize>,
    /// Between the start and end tags (empty, at the end of the start tag, when self-closing).
    pub inner: Range<usize>,
}

impl Element {
    pub fn is(&self, ns: &str, local: &str) -> bool {
        self.ns == ns && self.local == local
    }

    pub fn attr(&self, ns: &str, local: &str) -> Option<&str> {
        self.attrs
            .iter()
            .find(|a| a.ns == ns && a.local == local)
            .map(|a| a.value.as_str())
    }

    /// `rdf:about`.
    pub fn about(&self) -> Option<&str> {
        self.attr(RDF, "about")
    }

    /// `rdf:resource`.
    pub fn resource(&self) -> Option<&str> {
        self.attr(RDF, "resource")
    }

    pub fn children_named<'a>(
        &'a self,
        ns: &'a str,
        local: &'a str,
    ) -> impl Iterator<Item = &'a Element> + 'a {
        self.children.iter().filter(move |c| c.is(ns, local))
    }

    pub fn child(&self, ns: &str, local: &str) -> Option<&Element> {
        self.children.iter().find(|c| c.is(ns, local))
    }

    /// The trimmed text of a child property.
    pub fn text_of(&self, ns: &str, local: &str) -> Option<String> {
        self.child(ns, local).map(|c| c.text.trim().to_string())
    }

    /// The references of child properties (`rdf:resource`, or a nested node's `rdf:about`).
    pub fn references(&self, ns: &str, local: &str) -> Vec<String> {
        self.children_named(ns, local)
            .filter_map(|c| {
                c.resource()
                    .or_else(|| c.children.first().and_then(Element::about))
                    .map(str::to_string)
            })
            .collect()
    }

    /// The node elements nested in child properties (`<sbol:location><sbol:Range ...>`).
    pub fn nested<'a>(
        &'a self,
        ns: &'a str,
        local: &'a str,
    ) -> impl Iterator<Item = &'a Element> + 'a {
        self.children_named(ns, local)
            .filter_map(|c| c.children.first())
    }
}

/// Split `prefix:local` and resolve the prefix in the scopes (innermost last).
fn resolve(name: &str, scopes: &[Vec<(String, String)>], is_attr: bool) -> (String, String) {
    let (prefix, local) = match name.split_once(':') {
        Some((p, l)) => (p, l),
        None => ("", name),
    };
    if prefix.is_empty() && is_attr {
        return (String::new(), local.to_string());
    }
    let ns = scopes
        .iter()
        .rev()
        .flat_map(|s| s.iter())
        .find(|(p, _)| p == prefix)
        .map(|(_, uri)| uri.clone())
        .unwrap_or_default();
    (ns, local.to_string())
}

/// Read XML into an element tree (the document element), without recursion.
pub fn parse_tree(text: &str) -> Result<Element, SbolError> {
    let mut tokens = Tokenizer::new(text);
    let mut stack: Vec<Element> = Vec::new();
    let mut scopes: Vec<Vec<(String, String)>> = Vec::new();
    let mut root: Option<Element> = None;
    let mut count = 0usize;
    let xml_err = |offset: usize, reason: &str| SbolError::Xml {
        offset,
        reason: reason.to_string(),
    };
    while let Some(next) = tokens.next_token() {
        let (token, span) = next.map_err(|e| SbolError::Xml {
            offset: e.offset,
            reason: e.reason.to_string(),
        })?;
        match token {
            Token::Start {
                name,
                attrs,
                self_closing,
            } => {
                if root.is_some() {
                    return Err(xml_err(span.start, "content after the document element"));
                }
                count += 1;
                if count > MAX_ELEMENTS || stack.len() >= MAX_DEPTH {
                    return Err(xml_err(
                        span.start,
                        "too many or too deeply nested elements",
                    ));
                }
                let mut scope = Vec::new();
                for a in &attrs {
                    if a.name == "xmlns" {
                        scope.push((String::new(), decode(a.raw)));
                    } else if let Some(p) = a.name.strip_prefix("xmlns:") {
                        scope.push((p.to_string(), decode(a.raw)));
                    }
                }
                let namespaces = scope.clone();
                scopes.push(scope);
                let (ns, local) = resolve(name, &scopes, false);
                let attrs = attrs
                    .iter()
                    .filter(|a| a.name != "xmlns" && !a.name.starts_with("xmlns:"))
                    .map(|a| {
                        let (ns, local) = resolve(a.name, &scopes, true);
                        Attr {
                            ns,
                            local,
                            value: decode(a.raw),
                        }
                    })
                    .collect();
                let element = Element {
                    ns,
                    local,
                    attrs,
                    namespaces,
                    children: Vec::new(),
                    text: String::new(),
                    span: span.clone(),
                    inner: span.end..span.end,
                };
                if self_closing {
                    scopes.pop();
                    match stack.last_mut() {
                        Some(parent) => parent.children.push(element),
                        None => root = Some(element),
                    }
                } else {
                    stack.push(element);
                }
            }
            Token::End { name } => {
                let mut element = stack
                    .pop()
                    .ok_or_else(|| xml_err(span.start, "an end tag without a start tag"))?;
                let (ns, local) = resolve(name, &scopes, false);
                scopes.pop();
                if ns != element.ns || local != element.local {
                    return Err(xml_err(
                        span.start,
                        "the end tag doesn't match its start tag",
                    ));
                }
                element.inner = element.span.end..span.start;
                element.span = element.span.start..span.end;
                match stack.last_mut() {
                    Some(parent) => parent.children.push(element),
                    None => root = Some(element),
                }
            }
            Token::Text(t) => {
                if let Some(current) = stack.last_mut() {
                    current.text.push_str(&decode(t));
                } else if !t.trim().is_empty() {
                    return Err(xml_err(span.start, "text outside the document element"));
                }
            }
            Token::Other(raw) => {
                // CDATA sections carry text; comments, declarations and processing
                // instructions don't.
                if let (Some(current), Some(cdata)) = (
                    stack.last_mut(),
                    raw.strip_prefix("<![CDATA[")
                        .and_then(|r| r.strip_suffix("]]>")),
                ) {
                    current.text.push_str(cdata);
                }
            }
        }
    }
    if !stack.is_empty() {
        return Err(xml_err(text.len(), "the document ends inside an element"));
    }
    root.ok_or_else(|| xml_err(0, "no document element"))
}

/// Where a feature lies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Location {
    /// 1-based, inclusive.
    Range {
        start: u64,
        end: u64,
        orientation: Option<String>,
    },
    /// Between `at` and `at + 1` (0 is before the first residue).
    Cut {
        at: u64,
        orientation: Option<String>,
    },
    /// Somewhere in the sequence, position unknown.
    Generic { orientation: Option<String> },
}

impl Location {
    pub fn orientation(&self) -> Option<&str> {
        match self {
            Location::Range { orientation, .. }
            | Location::Cut { orientation, .. }
            | Location::Generic { orientation } => orientation.as_deref(),
        }
    }
}

/// A `SequenceAnnotation`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Annotation {
    pub uri: String,
    pub display_id: Option<String>,
    pub title: Option<String>,
    pub description: Option<String>,
    pub roles: Vec<String>,
    /// The `Component` it annotates, if any.
    pub component: Option<String>,
    pub locations: Vec<Location>,
    /// The `sbol:sequenceAnnotation` property element in the source.
    pub span: Range<usize>,
}

/// A `Component` (an instance of another definition inside this one).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Component {
    pub uri: String,
    pub definition: Option<String>,
}

/// A `ComponentDefinition`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ComponentDefinition {
    pub uri: String,
    /// The URI without version, which references may use for "the latest version".
    pub persistent_identity: Option<String>,
    pub display_id: Option<String>,
    pub title: Option<String>,
    pub description: Option<String>,
    pub types: Vec<String>,
    pub roles: Vec<String>,
    pub sequences: Vec<String>,
    pub components: Vec<Component>,
    pub annotations: Vec<Annotation>,
    /// The element in the source; `inner.end` is where new properties go.
    pub span: Range<usize>,
    pub inner: Range<usize>,
    /// Spans of the `dcterms:title`, `dcterms:description` and topology `sbol:type` properties.
    pub title_span: Option<Range<usize>>,
    pub description_span: Option<Range<usize>>,
    pub topology_span: Option<Range<usize>>,
}

impl ComponentDefinition {
    pub fn circular(&self) -> bool {
        self.types.iter().any(|t| t == SO_CIRCULAR)
    }
}

/// A `Sequence`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    pub uri: String,
    pub persistent_identity: Option<String>,
    pub display_id: Option<String>,
    pub elements: String,
    pub encoding: Option<String>,
    /// The text inside `sbol:elements` in the source.
    pub elements_span: Option<Range<usize>>,
}

/// The parts of an SBOL 2 document a sequence editor uses, over the element tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SbolDocument {
    pub root: Element,
    pub component_definitions: Vec<ComponentDefinition>,
    pub sequences: Vec<Sequence>,
}

fn location(node: &Element) -> Option<Location> {
    let orientation = node
        .child(SBOL2, "orientation")
        .and_then(Element::resource)
        .map(str::to_string);
    let number = |name: &str| {
        node.text_of(SBOL2, name)
            .and_then(|t| t.parse::<u64>().ok())
    };
    match node.local.as_str() {
        "Range" if node.ns == SBOL2 => Some(Location::Range {
            start: number("start")?,
            end: number("end")?,
            orientation,
        }),
        "Cut" if node.ns == SBOL2 => Some(Location::Cut {
            at: number("at")?,
            orientation,
        }),
        "GenericLocation" if node.ns == SBOL2 => Some(Location::Generic { orientation }),
        _ => None,
    }
}

fn annotation(property: &Element) -> Option<Annotation> {
    let node = property.children.first()?;
    if !node.is(SBOL2, "SequenceAnnotation") {
        return None;
    }
    Some(Annotation {
        uri: node.about().unwrap_or("").to_string(),
        display_id: node.text_of(SBOL2, "displayId"),
        title: node.text_of(DCTERMS, "title"),
        description: node.text_of(DCTERMS, "description"),
        roles: node.references(SBOL2, "role"),
        component: node.references(SBOL2, "component").into_iter().next(),
        locations: node
            .nested(SBOL2, "location")
            .filter_map(location)
            .collect(),
        span: property.span.clone(),
    })
}

fn component_definition(e: &Element) -> ComponentDefinition {
    let types = e.references(SBOL2, "type");
    ComponentDefinition {
        uri: e.about().unwrap_or("").to_string(),
        persistent_identity: e.references(SBOL2, "persistentIdentity").into_iter().next(),
        display_id: e.text_of(SBOL2, "displayId"),
        title: e.text_of(DCTERMS, "title"),
        description: e.text_of(DCTERMS, "description"),
        roles: e.references(SBOL2, "role"),
        sequences: e.references(SBOL2, "sequence"),
        components: e
            .nested(SBOL2, "component")
            .filter(|c| c.is(SBOL2, "Component"))
            .map(|c| Component {
                uri: c.about().unwrap_or("").to_string(),
                definition: c.references(SBOL2, "definition").into_iter().next(),
            })
            .collect(),
        annotations: e
            .children_named(SBOL2, "sequenceAnnotation")
            .filter_map(annotation)
            .collect(),
        span: e.span.clone(),
        inner: e.inner.clone(),
        title_span: e.child(DCTERMS, "title").map(|c| c.span.clone()),
        description_span: e.child(DCTERMS, "description").map(|c| c.span.clone()),
        topology_span: e
            .children_named(SBOL2, "type")
            .find(|c| matches!(c.resource(), Some(SO_LINEAR | SO_CIRCULAR)))
            .map(|c| c.span.clone()),
        types,
    }
}

fn sequence(e: &Element) -> Sequence {
    let elements = e.child(SBOL2, "elements");
    Sequence {
        uri: e.about().unwrap_or("").to_string(),
        persistent_identity: e.references(SBOL2, "persistentIdentity").into_iter().next(),
        display_id: e.text_of(SBOL2, "displayId"),
        elements: elements
            .map(|x| x.text.trim().to_string())
            .unwrap_or_default(),
        encoding: e
            .child(SBOL2, "encoding")
            .and_then(Element::resource)
            .map(str::to_string),
        elements_span: elements.map(|x| x.inner.clone()),
    }
}

/// Whether `bytes` look like SBOL (any version) in RDF/XML.
pub fn looks_like_sbol(bytes: &[u8]) -> bool {
    let head = bytes.get(..bytes.len().min(4096)).unwrap_or(bytes);
    let head = String::from_utf8_lossy(head);
    head.contains("rdf:RDF") && head.contains("http://sbols.org/")
}

impl SbolDocument {
    pub fn parse(bytes: &[u8]) -> Result<Self, SbolError> {
        let text = std::str::from_utf8(bytes).map_err(|_| SbolError::NotText)?;
        let root = parse_tree(text)?;
        if !root.is(RDF, "RDF") {
            return Err(SbolError::NotSbol2(
                "the document element isn't rdf:RDF".into(),
            ));
        }
        let uses = |ns: &str| root.children.iter().any(|c| c.ns == ns);
        if !uses(SBOL2) {
            let which = if uses(SBOL3) {
                "it is SBOL 3, which isn't supported yet"
            } else if uses(SBOL1) {
                "it is SBOL 1, which isn't supported"
            } else {
                "it has no SBOL 2 objects"
            };
            return Err(SbolError::NotSbol2(which.into()));
        }
        let component_definitions = root
            .children_named(SBOL2, "ComponentDefinition")
            .map(component_definition)
            .collect();
        let sequences = root
            .children_named(SBOL2, "Sequence")
            .map(sequence)
            .collect();
        Ok(SbolDocument {
            root,
            component_definitions,
            sequences,
        })
    }

    /// The sequence a reference names: by URI, or by persistent identity (the latest
    /// version, as SBOL 2 reads such references).
    /// The prefix the document element declares for a namespace.
    pub fn prefix(&self, ns: &str) -> Option<&str> {
        self.root
            .namespaces
            .iter()
            .find(|(_, uri)| uri == ns)
            .map(|(p, _)| p.as_str())
    }

    pub fn sequence(&self, uri: &str) -> Option<&Sequence> {
        self.sequences.iter().find(|s| s.uri == uri).or_else(|| {
            self.sequences
                .iter()
                .find(|s| s.persistent_identity.as_deref() == Some(uri))
        })
    }

    pub fn definition(&self, uri: &str) -> Option<&ComponentDefinition> {
        self.component_definitions
            .iter()
            .find(|c| c.uri == uri)
            .or_else(|| {
                self.component_definitions
                    .iter()
                    .find(|c| c.persistent_identity.as_deref() == Some(uri))
            })
    }

    /// The definition a sequence editor shows: one with a sequence that no other definition
    /// uses as a component (the top of the design), else the first with a sequence.
    pub fn main_definition(&self) -> Option<&ComponentDefinition> {
        let used: Vec<&str> = self
            .component_definitions
            .iter()
            .flat_map(|c| c.components.iter())
            .filter_map(|c| c.definition.as_deref())
            .collect();
        let with_sequence =
            |c: &&ComponentDefinition| c.sequences.iter().any(|s| self.sequence(s).is_some());
        let is_used = |c: &ComponentDefinition| {
            used.contains(&c.uri.as_str())
                || c.persistent_identity
                    .as_deref()
                    .is_some_and(|p| used.contains(&p))
        };
        self.component_definitions
            .iter()
            .filter(with_sequence)
            .find(|c| !is_used(c))
            .or_else(|| self.component_definitions.iter().find(with_sequence))
    }
}

/// Escape text for element content.
pub fn escape_text(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Escape text for an attribute value in double quotes.
pub fn escape_attr(s: &str) -> String {
    escape_text(s).replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::indexing_slicing)]
    use super::*;

    const SAMPLE: &str = r#"<?xml version="1.0" ?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:sbol="http://sbols.org/v2#" xmlns:dcterms="http://purl.org/dc/terms/">
  <sbol:ComponentDefinition rdf:about="http://x.org/part">
    <sbol:displayId>part</sbol:displayId>
    <dcterms:title>A &amp; B</dcterms:title>
    <sbol:type rdf:resource="http://www.biopax.org/release/biopax-level3.owl#DnaRegion"/>
    <sbol:type rdf:resource="http://identifiers.org/so/SO:0000988"/>
    <sbol:sequenceAnnotation>
      <sbol:SequenceAnnotation rdf:about="http://x.org/part/a1">
        <sbol:displayId>a1</sbol:displayId>
        <sbol:location>
          <sbol:Range rdf:about="http://x.org/part/a1/r">
            <sbol:start>3</sbol:start>
            <sbol:end>8</sbol:end>
            <sbol:orientation rdf:resource="http://sbols.org/v2#reverseComplement"/>
          </sbol:Range>
        </sbol:location>
        <sbol:role rdf:resource="http://identifiers.org/so/SO:0000167"/>
      </sbol:SequenceAnnotation>
    </sbol:sequenceAnnotation>
    <sbol:sequence rdf:resource="http://x.org/seq"/>
  </sbol:ComponentDefinition>
  <sbol:Sequence rdf:about="http://x.org/seq">
    <sbol:elements>acgtacgtac</sbol:elements>
    <sbol:encoding rdf:resource="http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"/>
  </sbol:Sequence>
</rdf:RDF>
"#;

    #[test]
    fn reads_definitions_annotations_and_sequences() {
        let d = SbolDocument::parse(SAMPLE.as_bytes()).unwrap();
        let cd = d.main_definition().unwrap();
        assert_eq!(cd.title.as_deref(), Some("A & B"));
        assert!(cd.circular());
        assert_eq!(d.sequence(&cd.sequences[0]).unwrap().elements, "acgtacgtac");
        let a = &cd.annotations[0];
        assert_eq!(
            a.roles,
            vec!["http://identifiers.org/so/SO:0000167".to_string()]
        );
        assert_eq!(
            a.locations,
            vec![Location::Range {
                start: 3,
                end: 8,
                orientation: Some(REVERSE_COMPLEMENT.into())
            }]
        );
        // Spans point at the source text.
        assert!(SAMPLE[a.span.clone()].starts_with("<sbol:sequenceAnnotation>"));
        assert!(SAMPLE[cd.topology_span.clone().unwrap()].contains("SO:0000988"));
        let s = d.sequence("http://x.org/seq").unwrap();
        assert_eq!(&SAMPLE[s.elements_span.clone().unwrap()], "acgtacgtac");
        assert!(SAMPLE[cd.inner.end..].starts_with("</sbol:ComponentDefinition>"));
    }

    #[test]
    fn other_prefixes_and_versions() {
        // Any prefix works: it's the namespace that counts.
        let renamed = SAMPLE
            .replace("sbol:", "s2:")
            .replace("xmlns:sbol=", "xmlns:s2=");
        let d = SbolDocument::parse(renamed.as_bytes()).unwrap();
        assert_eq!(d.component_definitions.len(), 1);
        let v3 = r#"<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:sbol="http://sbols.org/v3#"><sbol:Component rdf:about="x"/></rdf:RDF>"#;
        let e = SbolDocument::parse(v3.as_bytes()).unwrap_err();
        assert!(e.to_string().contains("SBOL 3"));
    }

    #[test]
    fn hostile_input_is_an_error() {
        for bad in [
            "",
            "<a>",
            "<a></b>",
            "</a>",
            "<a/><b/>",
            "text",
            "<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"/>",
        ] {
            assert!(SbolDocument::parse(bad.as_bytes()).is_err(), "{bad:?}");
        }
        let deep = "<a>".repeat(MAX_DEPTH + 5);
        assert!(parse_tree(&deep).is_err());
    }
}
