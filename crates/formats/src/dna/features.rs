//! The features packet: XML with one `<Feature>` per annotation.
//!
//! Every attribute and child element is kept in file order, including ones we don't model, and
//! each feature keeps its original text until it is modified. Unmodified features are therefore
//! written back exactly as read, and modified or new ones follow the layout observed in files
//! written by the reference application.

use std::ops::Range;

use super::{DnaError, DnaFile, tag};
use crate::xml::{self, RawAttr, Token, Tokenizer, XmlError};

/// Attributes with decoded values, in file order.
pub type Attrs = Vec<(String, String)>;

/// A child element of `<Feature>`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FeatureChild {
    /// `<Segment range="a-b" .../>`.
    Segment(Attrs),
    /// `<Q name="..."><V .../>...</Q>`.
    Qualifier { attrs: Attrs, values: Vec<Attrs> },
    /// Anything else, kept verbatim.
    Raw(String),
}

/// One `<Feature>` element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureXml {
    pub attrs: Attrs,
    pub children: Vec<FeatureChild>,
    pub self_closing: bool,
}

/// A segment's position: 1-based and inclusive, as stored. `start > end` wraps past the
/// origin of a circular molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SegmentInfo {
    pub start: u64,
    pub end: u64,
    pub name: Option<String>,
    pub color: Option<String>,
}

/// Default colour of a feature created in the reference application (observed).
pub const DEFAULT_COLOR: &str = "#a6acb3";

impl FeatureXml {
    /// A new feature laid out like one created in the reference application: `recentID`,
    /// `name`, optional `directionality` (1 forward, 2 reverse), `type`, then the two default
    /// flags; segments get `color` (or the default colour); qualifiers hold text values.
    pub fn new(
        recent_id: u64,
        color: Option<&str>,
        name: &str,
        kind: &str,
        directionality: Option<u8>,
        segments: &[(u64, u64)],
        qualifiers: &[(String, Vec<String>)],
    ) -> Self {
        let mut attrs: Attrs = vec![
            ("recentID".into(), recent_id.to_string()),
            ("name".into(), name.into()),
        ];
        if let Some(d) = directionality {
            attrs.push(("directionality".into(), d.to_string()));
        }
        attrs.push(("type".into(), kind.into()));
        attrs.push(("allowSegmentOverlaps".into(), "0".into()));
        attrs.push(("consecutiveTranslationNumbering".into(), "1".into()));
        let color = color.unwrap_or(DEFAULT_COLOR);
        let mut children: Vec<FeatureChild> = segments
            .iter()
            .map(|(a, b)| {
                FeatureChild::Segment(vec![
                    ("range".into(), format!("{a}-{b}")),
                    ("color".into(), color.into()),
                    ("type".into(), "standard".into()),
                ])
            })
            .collect();
        children.extend(qualifiers.iter().map(|(name, values)| {
            FeatureChild::Qualifier {
                attrs: vec![("name".into(), name.clone())],
                values: values
                    .iter()
                    .map(|v| vec![("text".into(), v.clone())])
                    .collect(),
            }
        }));
        FeatureXml {
            attrs,
            children,
            self_closing: false,
        }
    }

    pub fn attr(&self, name: &str) -> Option<&str> {
        get(&self.attrs, name)
    }

    /// Set an attribute, keeping its position if it exists, appending it otherwise.
    pub fn set_attr(&mut self, name: &str, value: &str) {
        set(&mut self.attrs, name, value);
    }

    pub fn remove_attr(&mut self, name: &str) {
        self.attrs.retain(|(n, _)| n != name);
    }

    pub fn name(&self) -> &str {
        self.attr("name").unwrap_or("")
    }

    /// The feature type (`CDS`, `promoter`, ...).
    pub fn kind(&self) -> Option<&str> {
        self.attr("type")
    }

    /// 1 = forward, 2 = reverse (observed); absent = not directional.
    pub fn directionality(&self) -> Option<&str> {
        self.attr("directionality")
    }

    /// The segments' positions, in file order (ascending position in every file so far,
    /// also for reverse-strand features).
    pub fn segments(&self) -> Result<Vec<SegmentInfo>, String> {
        self.children
            .iter()
            .filter_map(|c| match c {
                FeatureChild::Segment(attrs) => Some(attrs),
                _ => None,
            })
            .map(|attrs| {
                let range = get(attrs, "range").ok_or("a segment has no range")?;
                let (a, b) = range
                    .split_once('-')
                    .ok_or_else(|| format!("segment range {range:?} isn't \"start-end\""))?;
                let parse = |s: &str| {
                    s.trim()
                        .parse::<u64>()
                        .map_err(|_| format!("segment range {range:?} has a bad number"))
                };
                Ok(SegmentInfo {
                    start: parse(a)?,
                    end: parse(b)?,
                    name: get(attrs, "name").map(str::to_string),
                    color: get(attrs, "color").map(str::to_string),
                })
            })
            .collect()
    }

    /// Replace the segments' ranges (1-based, inclusive), keeping the other attributes of the
    /// existing segments in order. Extra ranges get new default segments; surplus segments are
    /// dropped.
    pub fn set_segment_ranges(&mut self, ranges: &[(u64, u64)]) {
        let mut ranges = ranges.iter();
        let mut kept = Vec::with_capacity(self.children.len());
        for child in self.children.drain(..) {
            match child {
                FeatureChild::Segment(mut attrs) => {
                    if let Some((a, b)) = ranges.next() {
                        set(&mut attrs, "range", &format!("{a}-{b}"));
                        kept.push(FeatureChild::Segment(attrs));
                    }
                }
                other => kept.push(other),
            }
        }
        let insert_at = kept
            .iter()
            .rposition(|c| matches!(c, FeatureChild::Segment(_)))
            .map_or(0, |i| i.saturating_add(1));
        let extra: Vec<FeatureChild> = ranges
            .map(|(a, b)| {
                FeatureChild::Segment(vec![
                    ("range".into(), format!("{a}-{b}")),
                    ("color".into(), DEFAULT_COLOR.into()),
                    ("type".into(), "standard".into()),
                ])
            })
            .collect();
        kept.splice(insert_at..insert_at, extra);
        self.children = kept;
        self.self_closing = false;
    }

    /// Qualifiers as (name, values). Each value is the first attribute of its `<V>` (`text`,
    /// `int`, ...), decoded.
    pub fn qualifiers(&self) -> Vec<(String, Vec<String>)> {
        self.children
            .iter()
            .filter_map(|c| match c {
                FeatureChild::Qualifier { attrs, values } => Some((
                    get(attrs, "name").unwrap_or("").to_string(),
                    values
                        .iter()
                        .filter_map(|v| v.first().map(|(_, value)| value.clone()))
                        .collect(),
                )),
                _ => None,
            })
            .collect()
    }

    /// Replace all qualifiers with text-valued ones, keeping the other children in place.
    pub fn set_qualifiers(&mut self, qualifiers: &[(String, Vec<String>)]) {
        let first_q = self
            .children
            .iter()
            .position(|c| matches!(c, FeatureChild::Qualifier { .. }))
            .unwrap_or(self.children.len());
        self.children
            .retain(|c| !matches!(c, FeatureChild::Qualifier { .. }));
        let at = first_q.min(self.children.len());
        let new = qualifiers
            .iter()
            .map(|(name, values)| FeatureChild::Qualifier {
                attrs: vec![("name".into(), name.clone())],
                values: values
                    .iter()
                    .map(|v| vec![("text".into(), v.clone())])
                    .collect(),
            });
        self.children.splice(at..at, new);
        if !self.children.is_empty() {
            self.self_closing = false;
        }
    }

    /// The element as the reference application writes it.
    pub fn to_xml(&self) -> String {
        let mut out = String::from("<Feature");
        push_attrs(&mut out, &self.attrs);
        if self.self_closing && self.children.is_empty() {
            out.push_str("/>");
            return out;
        }
        out.push('>');
        for child in &self.children {
            match child {
                FeatureChild::Segment(attrs) => {
                    out.push_str("<Segment");
                    push_attrs(&mut out, attrs);
                    out.push_str("/>");
                }
                FeatureChild::Qualifier { attrs, values } => {
                    out.push_str("<Q");
                    push_attrs(&mut out, attrs);
                    out.push('>');
                    for v in values {
                        out.push_str("<V");
                        push_attrs(&mut out, v);
                        out.push_str("/>");
                    }
                    out.push_str("</Q>");
                }
                FeatureChild::Raw(raw) => out.push_str(raw),
            }
        }
        out.push_str("</Feature>");
        out
    }
}

fn get<'a>(attrs: &'a Attrs, name: &str) -> Option<&'a str> {
    attrs
        .iter()
        .find(|(n, _)| n == name)
        .map(|(_, v)| v.as_str())
}

fn set(attrs: &mut Attrs, name: &str, value: &str) {
    match attrs.iter_mut().find(|(n, _)| n == name) {
        Some(slot) => slot.1 = value.to_string(),
        None => attrs.push((name.to_string(), value.to_string())),
    }
}

fn push_attrs(out: &mut String, attrs: &Attrs) {
    for (name, value) in attrs {
        out.push(' ');
        out.push_str(name);
        out.push_str("=\"");
        out.push_str(&xml::encode_attr(value));
        out.push('"');
    }
}

fn decode_attrs(attrs: &[RawAttr<'_>]) -> Attrs {
    attrs
        .iter()
        .map(|a| (a.name.to_string(), xml::decode(a.raw)))
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Item {
    /// `raw` is the original text while the feature is unmodified.
    Feature {
        xml: FeatureXml,
        raw: Option<String>,
    },
    /// Anything between features (whitespace, unknown elements), verbatim.
    Raw(String),
}

/// The whole features packet.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeaturesXml {
    /// Everything before the `<Features>` tag (the XML declaration).
    prolog: String,
    root_attrs: Attrs,
    /// The original `<Features ...>` tag while its attributes are unmodified.
    root_raw: Option<String>,
    root_self_closing: bool,
    items: Vec<Item>,
    /// `</Features>` and whatever follows it, verbatim (or what followed a self-closing root).
    tail: String,
}

impl FeaturesXml {
    /// An empty features packet in the layout the reference application uses.
    pub fn empty() -> Self {
        FeaturesXml {
            prolog: "<?xml version=\"1.0\"?>".into(),
            root_attrs: vec![("nextValidID".into(), "0".into())],
            root_raw: None,
            root_self_closing: false,
            items: Vec::new(),
            tail: "</Features>\n".into(),
        }
    }

    pub fn parse(text: &str) -> Result<Self, XmlError> {
        let mut tokens = Tokenizer::new(text);
        let mut next = || tokens.next_token().transpose();

        // Up to and including the root tag.
        let (root_attrs, root_span, root_self_closing) = loop {
            match next()? {
                None => {
                    return Err(XmlError {
                        offset: text.len(),
                        reason: "no <Features> element",
                    });
                }
                Some((
                    Token::Start {
                        name: "Features",
                        attrs,
                        self_closing,
                    },
                    span,
                )) => break (decode_attrs(&attrs), span, self_closing),
                Some((Token::Other(_) | Token::Text(_), _)) => {}
                Some((_, span)) => {
                    return Err(XmlError {
                        offset: span.start,
                        reason: "unexpected element before <Features>",
                    });
                }
            }
        };
        let slice = |r: Range<usize>| text.get(r).unwrap_or("").to_string();
        let mut parsed = FeaturesXml {
            prolog: slice(0..root_span.start),
            root_attrs,
            root_raw: Some(slice(root_span.clone())),
            root_self_closing,
            items: Vec::new(),
            tail: String::new(),
        };
        if root_self_closing {
            parsed.tail = slice(root_span.end..text.len());
            return Ok(parsed);
        }

        // Children of <Features>, each collected as a complete element.
        loop {
            let Some((token, span)) = next()? else {
                return Err(XmlError {
                    offset: text.len(),
                    reason: "missing </Features>",
                });
            };
            match token {
                Token::End { name: "Features" } => {
                    parsed.tail = slice(span.start..text.len());
                    return Ok(parsed);
                }
                Token::Start {
                    name,
                    attrs,
                    self_closing,
                } => {
                    let mut body = Vec::new();
                    let end = if self_closing {
                        span.end
                    } else {
                        collect_element(&mut next, &mut body, span.start)?
                    };
                    let raw = slice(span.start..end);
                    if name == "Feature" {
                        let xml = FeatureXml {
                            attrs: decode_attrs(&attrs),
                            children: children(&body, text),
                            self_closing,
                        };
                        parsed.items.push(Item::Feature {
                            xml,
                            raw: Some(raw),
                        });
                    } else {
                        parsed.items.push(Item::Raw(raw));
                    }
                }
                Token::Text(_) | Token::Other(_) => parsed.items.push(Item::Raw(slice(span))),
                Token::End { .. } => {
                    return Err(XmlError {
                        offset: span.start,
                        reason: "unexpected end tag inside <Features>",
                    });
                }
            }
        }
    }

    /// The packet text: unmodified parts exactly as read.
    pub fn to_text(&self) -> String {
        let mut out = self.prolog.clone();
        match &self.root_raw {
            Some(raw) if !(self.root_self_closing && !self.items.is_empty()) => out.push_str(raw),
            _ => {
                out.push_str("<Features");
                push_attrs(&mut out, &self.root_attrs);
                out.push('>');
            }
        }
        for item in &self.items {
            match item {
                Item::Feature { raw: Some(raw), .. } => out.push_str(raw),
                Item::Feature { xml, raw: None } => out.push_str(&xml.to_xml()),
                Item::Raw(raw) => out.push_str(raw),
            }
        }
        if self.root_self_closing && !self.items.is_empty() {
            out.push_str("</Features>");
        }
        out.push_str(&self.tail);
        out
    }

    pub fn features(&self) -> impl Iterator<Item = &FeatureXml> {
        self.items.iter().filter_map(|i| match i {
            Item::Feature { xml, .. } => Some(xml),
            Item::Raw(_) => None,
        })
    }

    pub fn len(&self) -> usize {
        self.features().count()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Whether feature `index` is still exactly as read.
    pub fn is_unmodified(&self, index: usize) -> bool {
        self.feature_item(index)
            .is_some_and(|i| matches!(i, Item::Feature { raw: Some(_), .. }))
    }

    /// Mutable access to feature `index`; it will be re-serialised on write.
    pub fn feature_mut(&mut self, index: usize) -> Option<&mut FeatureXml> {
        match self.feature_item_mut(index)? {
            Item::Feature { xml, raw } => {
                *raw = None;
                Some(xml)
            }
            Item::Raw(_) => None,
        }
    }

    pub fn push(&mut self, feature: FeatureXml) {
        self.items.push(Item::Feature {
            xml: feature,
            raw: None,
        });
    }

    /// Keep only the features for which `keep(index, feature)` is true; other content stays.
    pub fn retain(&mut self, mut keep: impl FnMut(usize, &FeatureXml) -> bool) {
        let mut index = 0usize;
        self.items.retain(|item| match item {
            Item::Feature { xml, .. } => {
                let k = keep(index, xml);
                index = index.saturating_add(1);
                k
            }
            Item::Raw(_) => true,
        });
    }

    /// Reorder features by their current indices; `order` must be a permutation of the
    /// feature indices (anything else is ignored). Non-feature items stay where they are.
    pub fn reorder(&mut self, order: &[usize]) {
        let mut taken: Vec<Option<Item>> = Vec::new();
        let mut slots = Vec::new();
        for (pos, item) in self.items.iter().enumerate() {
            if matches!(item, Item::Feature { .. }) {
                slots.push(pos);
            }
        }
        let mut sorted = order.to_vec();
        sorted.sort_unstable();
        if sorted != (0..slots.len()).collect::<Vec<_>>() {
            return;
        }
        for &pos in &slots {
            taken.push(self.items.get(pos).cloned());
        }
        for (slot, &from) in slots.iter().zip(order) {
            if let (Some(dst), Some(Some(item))) = (self.items.get_mut(*slot), taken.get(from)) {
                *dst = item.clone();
            }
        }
    }

    /// `nextValidID` on the root, if present and numeric.
    pub fn next_valid_id(&self) -> Option<u64> {
        get(&self.root_attrs, "nextValidID").and_then(|v| v.parse().ok())
    }

    pub fn set_next_valid_id(&mut self, id: u64) {
        if self.next_valid_id() != Some(id) {
            set(&mut self.root_attrs, "nextValidID", &id.to_string());
            self.root_raw = None;
        }
    }

    fn feature_item(&self, index: usize) -> Option<&Item> {
        self.items
            .iter()
            .filter(|i| matches!(i, Item::Feature { .. }))
            .nth(index)
    }

    fn feature_item_mut(&mut self, index: usize) -> Option<&mut Item> {
        self.items
            .iter_mut()
            .filter(|i| matches!(i, Item::Feature { .. }))
            .nth(index)
    }
}

/// Collect the tokens inside an element whose start tag was just read, up to its matching end
/// tag. Returns the byte offset just after the end tag. Iterative, so nesting depth can't
/// overflow the stack.
fn collect_element<'a>(
    next: &mut impl FnMut() -> Result<Option<(Token<'a>, Range<usize>)>, XmlError>,
    body: &mut Vec<(Token<'a>, Range<usize>)>,
    start: usize,
) -> Result<usize, XmlError> {
    let mut depth = 1usize;
    loop {
        let Some((token, span)) = next()? else {
            return Err(XmlError {
                offset: start,
                reason: "element is never closed",
            });
        };
        match &token {
            Token::Start {
                self_closing: false,
                ..
            } => depth = depth.saturating_add(1),
            Token::End { .. } => {
                depth = depth.saturating_sub(1);
                if depth == 0 {
                    return Ok(span.end);
                }
            }
            _ => {}
        }
        body.push((token, span));
    }
}

/// Interpret the tokens inside a `<Feature>` as children; anything unexpected is kept raw.
fn children(body: &[(Token<'_>, Range<usize>)], text: &str) -> Vec<FeatureChild> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while let Some((token, span)) = body.get(i) {
        let raw = |r: Range<usize>| text.get(r).unwrap_or("").to_string();
        match token {
            Token::Start {
                name: "Segment",
                attrs,
                self_closing: true,
            } => {
                out.push(FeatureChild::Segment(decode_attrs(attrs)));
                i = i.saturating_add(1);
            }
            Token::Start {
                name,
                self_closing: false,
                attrs,
            } => {
                // Find the matching end within the body.
                let mut depth = 1usize;
                let mut j = i.saturating_add(1);
                while let Some((t, _)) = body.get(j) {
                    match t {
                        Token::Start {
                            self_closing: false,
                            ..
                        } => depth = depth.saturating_add(1),
                        Token::End { .. } => {
                            depth = depth.saturating_sub(1);
                            if depth == 0 {
                                break;
                            }
                        }
                        _ => {}
                    }
                    j = j.saturating_add(1);
                }
                let end = body.get(j).map_or(span.end, |(_, s)| s.end);
                let inner = body.get(i.saturating_add(1)..j).unwrap_or(&[]);
                let values: Option<Vec<Attrs>> = inner
                    .iter()
                    .map(|(t, _)| match t {
                        Token::Start {
                            name: "V",
                            attrs,
                            self_closing: true,
                        } => Some(decode_attrs(attrs)),
                        _ => None,
                    })
                    .collect();
                match (*name, values) {
                    ("Q", Some(values)) => out.push(FeatureChild::Qualifier {
                        attrs: decode_attrs(attrs),
                        values,
                    }),
                    _ => out.push(FeatureChild::Raw(raw(span.start..end))),
                }
                i = j.saturating_add(1);
            }
            _ => {
                out.push(FeatureChild::Raw(raw(span.clone())));
                i = i.saturating_add(1);
            }
        }
    }
    out
}

impl DnaFile {
    pub fn features(&self) -> Result<FeaturesXml, DnaError> {
        let text = super::views::utf8(self.payload(tag::FEATURES)?, tag::FEATURES)?;
        FeaturesXml::parse(text).map_err(|e| DnaError::InvalidPacket {
            tag: tag::FEATURES,
            reason: e.to_string(),
        })
    }

    pub fn set_features(&mut self, features: &FeaturesXml) -> Result<(), DnaError> {
        self.set_payload(tag::FEATURES, features.to_text().into_bytes())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Observed in a file where the reference application added a feature at bp 1-3.
    const START_CODON: &str = r##"<Feature recentID="1" name="start-codon" type="misc_feature" allowSegmentOverlaps="0" consecutiveTranslationNumbering="1"><Segment range="1-3" color="#a6acb3" type="standard"/></Feature>"##;

    const PACKET: &str = r##"<?xml version="1.0"?><Features nextValidID="2"><Feature recentID="0" name="lac &quot;up&quot;" directionality="2" type="promoter" swappedSegmentNumbering="1"><Segment name="-10" range="10-16" color="#ffffff" type="standard"/><Segment range="17-30" color="#ffffff" type="standard"/><Q name="note"><V text="&lt;html&gt;&lt;body&gt;aph(3')-Ia ®&lt;/body&gt;&lt;/html&gt;"/></Q><Q name="codon_start"><V int="1"/></Q><Odd x="1"><Inner/></Odd></Feature><Feature recentID="1" name="x" type="misc_feature"/></Features>
"##;

    #[test]
    fn new_features_match_the_reference_layout() {
        let f = FeatureXml::new(1, None, "start-codon", "misc_feature", None, &[(1, 3)], &[]);
        assert_eq!(f.to_xml(), START_CODON);
    }

    #[test]
    fn parses_and_writes_back_exactly() {
        let p = FeaturesXml::parse(PACKET).unwrap();
        assert_eq!(p.to_text(), PACKET);
        assert_eq!(p.len(), 2);
        assert_eq!(p.next_valid_id(), Some(2));
        let f = p.features().next().unwrap();
        assert_eq!(f.name(), "lac \"up\"");
        assert_eq!(f.directionality(), Some("2"));
        let segs = f.segments().unwrap();
        assert_eq!((segs[0].start, segs[0].end), (10, 16));
        assert_eq!(segs[0].name.as_deref(), Some("-10"));
        let q = f.qualifiers();
        assert_eq!(q[0].1[0], "<html><body>aph(3')-Ia ®</body></html>");
        assert_eq!(q[1], ("codon_start".to_string(), vec!["1".to_string()]));
        assert!(
            matches!(f.children.last(), Some(FeatureChild::Raw(r)) if r == "<Odd x=\"1\"><Inner/></Odd>")
        );
        // Re-serialising every feature (not using the raw text) also reproduces the input.
        for (i, f) in p.features().enumerate() {
            assert!(p.is_unmodified(i));
            assert!(PACKET.contains(&f.to_xml()), "{}", f.to_xml());
        }
    }

    #[test]
    fn edits_touch_only_what_changed() {
        let mut p = FeaturesXml::parse(PACKET).unwrap();
        p.feature_mut(0)
            .unwrap()
            .set_segment_ranges(&[(11, 16), (17, 31)]);
        p.push(FeatureXml::new(
            2,
            None,
            "start-codon",
            "misc_feature",
            None,
            &[(1, 3)],
            &[],
        ));
        p.set_next_valid_id(3);
        let text = p.to_text();
        assert!(
            text.contains(
                r##"<Segment name="-10" range="11-16" color="#ffffff" type="standard"/>"##
            )
        );
        assert!(text.contains(r#"<Features nextValidID="3">"#));
        assert!(text.ends_with("</Feature></Features>\n"));
        let again = FeaturesXml::parse(&text).unwrap();
        assert_eq!(again.len(), 3);
        p.retain(|i, _| i != 1);
        assert_eq!(p.len(), 2);
        p.reorder(&[1, 0]);
        assert_eq!(p.features().next().unwrap().name(), "start-codon");
    }

    #[test]
    fn empty_and_self_closing_roots() {
        let mut e = FeaturesXml::empty();
        assert_eq!(
            e.to_text(),
            "<?xml version=\"1.0\"?><Features nextValidID=\"0\"></Features>\n"
        );
        e.push(FeatureXml::new(
            0,
            None,
            "a",
            "misc_feature",
            Some(1),
            &[(1, 2)],
            &[],
        ));
        assert_eq!(FeaturesXml::parse(&e.to_text()).unwrap().len(), 1);
        let mut s = FeaturesXml::parse("<Features nextValidID=\"0\"/>\n").unwrap();
        assert_eq!(s.to_text(), "<Features nextValidID=\"0\"/>\n");
        s.push(FeatureXml::new(
            0,
            None,
            "a",
            "misc_feature",
            None,
            &[(1, 2)],
            &[],
        ));
        let text = s.to_text();
        assert!(text.starts_with("<Features nextValidID=\"0\"><Feature"));
        assert!(text.ends_with("</Features>\n"));
    }

    #[test]
    fn malformed_packets_are_errors() {
        for bad in [
            "",
            "<Feature/>",
            "<Features>",
            "<Features><Feature>",
            "<Features></Feature>",
            "<Features><Feature><Segment range=\"1-2\"></Feature></Features>",
        ] {
            assert!(FeaturesXml::parse(bad).is_err(), "{bad:?}");
        }
        let f = FeatureXml {
            attrs: vec![],
            children: vec![FeatureChild::Segment(vec![("range".into(), "x-2".into())])],
            self_closing: false,
        };
        assert!(f.segments().is_err());
    }
}
