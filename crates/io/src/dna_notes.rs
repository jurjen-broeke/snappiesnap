//! Writing changed primers (packet 0x05) and notes (0x06) into `.dna` files, in the layout of
//! the reference application's files, changing only what changed.
//!
//! - Primers: each `<Primer>` the document didn't change is written back as read; a changed or
//!   new one gets its attributes (other attributes kept in place) and binding sites written
//!   anew. After a sequence edit every primer's binding sites are rewritten, since their
//!   positions moved. Binding sites: 0-based inclusive `location`, `boundStrand` 0 (forward) or
//!   1 (reverse), the annealed bases, a Tm, one `<Component>`, and the same site again with
//!   `simplified="1"`.
//! - Notes: one element per line; only the lines of changed fields are rewritten (attributes
//!   such as `UTC` on dates are kept), removed fields are dropped, new ones go before
//!   `</Notes>`. Notes text escapes `&` and `<` but not `>` (as observed).

use snappiesnap_bio::primers::{BindingRules, PrimerStrand, binding_sites};
use snappiesnap_bio::tm::{TmConditions, tm};
use snappiesnap_doc::{Molecule, Notes, Primer, PrimerId};
use snappiesnap_formats::xml::{Token, Tokenizer, encode_attr};

use crate::IoError;

/// The primers packet of a file with no primers (as the reference application writes it).
pub(crate) const EMPTY_PRIMERS: &str = "<?xml version=\"1.0\"?><Primers nextValidID=\"0\"><HybridizationParams minContinuousMatchLen=\"10\" allowMismatch=\"1\" minMeltingTemperature=\"40\" showAdditionalFivePrimeMatches=\"1\" minimumFivePrimeAnnealing=\"15\"/></Primers>\n";

/// Shortest exact annealed stretch written as a binding site.
const MIN_ANNEALED: usize = 10;

/// A primers packet cut into the `<Primer>` elements and the text around them.
enum Piece {
    Other(String),
    Primer {
        raw: String,
        /// Attributes of the start tag: name and raw (escaped) value, in order.
        attrs: Vec<(String, String)>,
    },
}

fn split(text: &str) -> Result<Vec<Piece>, IoError> {
    let bad = |e: String| IoError::Format(format!("primers packet: {e}"));
    let mut pieces = Vec::new();
    let mut tokens = Tokenizer::new(text);
    let mut depth = 0usize;
    let mut last = 0usize;
    let mut open: Option<(usize, Vec<(String, String)>)> = None;
    while let Some(next) = tokens.next_token() {
        let (token, range) =
            next.map_err(|e| bad(format!("XML at byte {}: {}", e.offset, e.reason)))?;
        match token {
            Token::Start {
                name,
                attrs,
                self_closing,
            } => {
                if depth == 1 && name == "Primer" {
                    let attrs: Vec<(String, String)> = attrs
                        .iter()
                        .map(|a| (a.name.to_string(), a.raw.to_string()))
                        .collect();
                    if self_closing {
                        pieces.push(Piece::Other(
                            text.get(last..range.start).unwrap_or_default().to_string(),
                        ));
                        pieces.push(Piece::Primer {
                            raw: text.get(range.clone()).unwrap_or_default().to_string(),
                            attrs,
                        });
                        last = range.end;
                    } else {
                        open = Some((range.start, attrs));
                    }
                }
                if !self_closing {
                    depth += 1;
                }
            }
            Token::End { name } => {
                depth = depth.saturating_sub(1);
                if depth == 1
                    && name == "Primer"
                    && let Some((start, attrs)) = open.take()
                {
                    pieces.push(Piece::Other(
                        text.get(last..start).unwrap_or_default().to_string(),
                    ));
                    pieces.push(Piece::Primer {
                        raw: text.get(start..range.end).unwrap_or_default().to_string(),
                        attrs,
                    });
                    last = range.end;
                }
            }
            _ => {}
        }
    }
    pieces.push(Piece::Other(
        text.get(last..).unwrap_or_default().to_string(),
    ));
    Ok(pieces)
}

/// A description as the file keeps it: HTML.
fn html(description: &str) -> String {
    if description.trim_start().starts_with("<html") {
        return description.to_string();
    }
    let escaped = description
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;");
    format!("<html><body>{escaped}</body></html>")
}

/// The binding sites of `p` on `m`, as XML: exact sites of at least 10 annealed bases (or the
/// whole primer), the `minContinuousMatchLen` the reference application's packets declare.
fn sites_xml(m: &Molecule, p: &Primer) -> String {
    let seq = p.seq.as_bytes();
    let len = m.len();
    let mut out = String::new();
    if len == 0 {
        return out;
    }
    for site in binding_sites(m, &[seq], BindingRules::default()) {
        if site.mismatches > 0 || site.annealed < MIN_ANNEALED.min(seq.len()) {
            continue;
        }
        let annealed = seq
            .get(seq.len().saturating_sub(site.annealed)..)
            .unwrap_or_default();
        let bases = String::from_utf8_lossy(annealed).into_owned();
        let last = (site.range.end + len - 1) % len;
        let location = format!("{}-{last}", site.range.start);
        let strand = if site.strand == PrimerStrand::Forward {
            0
        } else {
            1
        };
        let melting = tm(annealed, &TmConditions::default())
            .filter(|t| t.is_finite())
            .map_or(0, |t| t.round() as i64);
        let body = format!(
            "location=\"{location}\" boundStrand=\"{strand}\" annealedBases=\"{b}\" meltingTemperature=\"{melting}\"><Component hybridizedRange=\"{location}\" bases=\"{b}\"/></BindingSite>",
            b = encode_attr(&bases)
        );
        out.push_str(&format!(
            "<BindingSite {body}<BindingSite simplified=\"1\" {body}"
        ));
    }
    out
}

/// A `<Primer>` element: `attrs` (raw) with name, sequence and description set from `p`.
fn primer_xml(m: &Molecule, p: &Primer, attrs: &[(String, String)]) -> String {
    let values = [
        ("name", encode_attr(&p.name)),
        (
            "sequence",
            encode_attr(&String::from_utf8_lossy(p.seq.as_bytes())),
        ),
        ("description", encode_attr(&html(&p.description))),
    ];
    let mut out = String::from("<Primer");
    let mut done = [false; 3];
    for (name, raw) in attrs {
        match values.iter().position(|(n, _)| n == name) {
            Some(i) => {
                let value = values.get(i).map_or("", |(_, v)| v.as_str());
                out.push_str(&format!(" {name}=\"{value}\""));
                if let Some(d) = done.get_mut(i) {
                    *d = true;
                }
            }
            None => out.push_str(&format!(" {name}=\"{raw}\"")),
        }
    }
    for ((name, value), d) in values.iter().zip(done) {
        if !d {
            out.push_str(&format!(" {name}=\"{value}\""));
        }
    }
    out.push('>');
    out.push_str(&sites_xml(m, p));
    out.push_str("</Primer>");
    out
}

/// The `nextValidID` of the `<Primers>` start tag in `text`, and the text with it replaced.
fn set_next_id(text: &str, next: u64) -> String {
    let key = "nextValidID=\"";
    let Some(at) = text.find(key) else {
        return text.to_string();
    };
    let from = at + key.len();
    let to = text
        .get(from..)
        .and_then(|r| r.find('"'))
        .map_or(from, |n| from + n);
    format!(
        "{}{next}{}",
        text.get(..from).unwrap_or_default(),
        text.get(to..).unwrap_or_default()
    )
}

/// The primers packet for `now`, from the file's packet (`raw`; `None`: a file without one)
/// read into `then`, where `map` pairs document primers with the file's `<Primer>` elements
/// (by index). Returns the packet and, per `<Primer>` element written, the document primer it
/// holds (`None` for elements the document doesn't show).
pub(crate) fn primers_packet(
    now: &Molecule,
    then: &Molecule,
    raw: Option<&str>,
    map: &[(PrimerId, usize)],
) -> Result<(String, Vec<Option<PrimerId>>), IoError> {
    let mut pieces = split(raw.unwrap_or(EMPTY_PRIMERS))?;
    let sequence_moved = now.seq() != then.seq() || now.topology() != then.topology();
    let mut owners = Vec::new();
    let mut next_id = 0u64;
    let mut index = 0usize;
    let mut out_pieces: Vec<String> = Vec::new();
    for piece in pieces.drain(..) {
        match piece {
            Piece::Other(text) => out_pieces.push(text),
            Piece::Primer { raw, attrs } => {
                if let Some(id) = attrs
                    .iter()
                    .find(|(n, _)| n == "recentID")
                    .and_then(|(_, v)| v.parse::<u64>().ok())
                {
                    next_id = next_id.max(id.saturating_add(1));
                }
                let owner = map.iter().find(|(_, i)| *i == index).map(|(id, _)| *id);
                index += 1;
                match owner {
                    // Not shown in the document: kept as read.
                    None => {
                        out_pieces.push(raw);
                        owners.push(None);
                    }
                    Some(id) => match now.primers().iter().find(|p| p.id == id) {
                        // Deleted.
                        None => {}
                        Some(p) => {
                            let same = then.primers().iter().find(|q| q.id == id) == Some(p);
                            out_pieces.push(if same && !sequence_moved {
                                raw
                            } else {
                                primer_xml(now, p, &attrs)
                            });
                            owners.push(Some(id));
                        }
                    },
                }
            }
        }
    }
    // New primers go before the closing tag, in the last piece of text.
    let last = out_pieces.pop().unwrap_or_default();
    let cut = last.rfind("</Primers>").unwrap_or(last.len());
    out_pieces.push(last.get(..cut).unwrap_or_default().to_string());
    let closing = last.get(cut..).unwrap_or_default().to_string();
    let mapped: Vec<PrimerId> = map.iter().map(|(id, _)| *id).collect();
    if let Some(head) = out_pieces.first() {
        let declared = head
            .find("nextValidID=\"")
            .and_then(|at| head.get(at + 13..))
            .and_then(|r| r.split('"').next())
            .and_then(|v| v.parse::<u64>().ok())
            .unwrap_or(0);
        next_id = next_id.max(declared);
    }
    let mut added = false;
    for p in now.primers() {
        if mapped.contains(&p.id) {
            continue;
        }
        let attrs = vec![("recentID".to_string(), next_id.to_string())];
        out_pieces.push(primer_xml(now, p, &attrs));
        owners.push(Some(p.id));
        next_id = next_id.saturating_add(1);
        added = true;
    }
    out_pieces.push(closing);
    if added && let Some(head) = out_pieces.first_mut() {
        *head = set_next_id(head, next_id);
    }
    Ok((out_pieces.concat(), owners))
}

/// Notes text escaped as the reference application does: `&` and `<`, not `>`.
fn escape_notes(s: &str) -> String {
    s.replace('&', "&amp;").replace('<', "&lt;")
}

/// The element name of a notes line (`<Name ...>...`), if it is one.
fn line_name(line: &str) -> Option<&str> {
    let rest = line.strip_prefix('<')?;
    if rest.starts_with('/') {
        return None;
    }
    let end = rest.find(['>', ' ', '/']).unwrap_or(rest.len());
    rest.get(..end)
}

/// The notes packet for `now`, from the file's packet (`raw`, or `None`) read into `then`.
/// `None` if the notes didn't change.
pub(crate) fn notes_packet(
    now: &Notes,
    then: &Notes,
    raw: Option<&str>,
) -> Result<Option<String>, IoError> {
    if now == then {
        return Ok(None);
    }
    // Field name -> new text (None: removed). Description is a field of its own in the file.
    let mut changes: Vec<(String, Option<String>)> = Vec::new();
    if now.description != then.description {
        changes.push((
            "Description".to_string(),
            (!now.description.is_empty()).then(|| html(&now.description)),
        ));
    }
    for (k, v) in &then.fields {
        match now.fields.get(k) {
            None => changes.push((k.clone(), None)),
            Some(n) if n != v => changes.push((k.clone(), Some(n.clone()))),
            _ => {}
        }
    }
    for (k, v) in &now.fields {
        if !then.fields.contains_key(k) {
            changes.push((k.clone(), Some(v.clone())));
        }
    }
    let raw = raw.unwrap_or("<Notes>\n</Notes>\n");
    let mut lines: Vec<String> = raw.split('\n').map(str::to_string).collect();
    let mut pending: Vec<(String, String)> = Vec::new();
    for (name, value) in changes {
        let found = lines
            .iter()
            .position(|l| line_name(l) == Some(name.as_str()));
        match (found, value) {
            (Some(i), None) => {
                lines.remove(i);
            }
            (Some(i), Some(v)) => {
                let line = lines.get(i).cloned().unwrap_or_default();
                let close = format!("</{name}>");
                let Some(tag_end) = line.find('>') else {
                    return Err(IoError::Unsupported(format!(
                        "saving the notes field {name} (it isn't on one line in this file)"
                    )));
                };
                if !line.ends_with(&close) && !line.ends_with("/>") {
                    return Err(IoError::Unsupported(format!(
                        "saving the notes field {name} (it isn't on one line in this file)"
                    )));
                }
                // Keep the start tag (and its attributes); a self-closing one opens now.
                let start = line
                    .get(..tag_end)
                    .unwrap_or_default()
                    .trim_end_matches('/');
                if let Some(l) = lines.get_mut(i) {
                    *l = format!("{start}>{}{close}", escape_notes(&v));
                }
            }
            (None, Some(v)) => pending.push((name, v)),
            (None, None) => {}
        }
    }
    let closing = lines
        .iter()
        .rposition(|l| l.trim() == "</Notes>")
        .ok_or_else(|| IoError::Format("notes packet without </Notes>".into()))?;
    for (offset, (name, v)) in pending.into_iter().enumerate() {
        lines.insert(
            closing + offset,
            format!("<{name}>{}</{name}>", escape_notes(&v)),
        );
    }
    Ok(Some(lines.join("\n")))
}

#[cfg(test)]
mod tests {
    use super::*;
    use snappiesnap_doc::{Alphabet, Seq, Topology};

    const PACKET: &str = "<?xml version=\"1.0\"?><Primers nextValidID=\"3\"><HybridizationParams minContinuousMatchLen=\"10\" allowMismatch=\"1\" minMeltingTemperature=\"40\" showAdditionalFivePrimeMatches=\"1\" minimumFivePrimeAnnealing=\"15\"/><Primer recentID=\"1\" name=\"f\" sequence=\"ACGTTGCAAGGCTTAC\" description=\"&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;\" dateAdded=\"2024-04-11T11:06:53Z\"><BindingSite location=\"0-15\" boundStrand=\"0\" annealedBases=\"ACGTTGCAAGGCTTAC\" meltingTemperature=\"50\"><Component hybridizedRange=\"0-15\" bases=\"ACGTTGCAAGGCTTAC\"/></BindingSite></Primer><Primer recentID=\"2\" name=\"hidden\" sequence=\"GGGG\"/></Primers>\n";

    fn molecule() -> Molecule {
        let mut m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, "ACGTTGCAAGGCTTACGATCGGATCCTTAGGCATGCAAGCTT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        m.add_primer(
            "f",
            Seq::new(Alphabet::Dna, "ACGTTGCAAGGCTTAC").unwrap(),
            "<html><body></body></html>",
        )
        .unwrap();
        m
    }

    #[test]
    fn unchanged_primers_are_written_back_as_read() {
        let m = molecule();
        let map = [(m.primers()[0].id, 0)];
        let (text, owners) = primers_packet(&m, &m, Some(PACKET), &map).unwrap();
        assert_eq!(text, PACKET);
        assert_eq!(owners, vec![Some(m.primers()[0].id), None]);
    }

    #[test]
    fn edited_and_new_primers() {
        let then = molecule();
        let mut now = then.clone();
        let mut p = now.primers()[0].clone();
        p.name = "f & g".into();
        now.update_primer(p).unwrap();
        // A reverse primer on AAGCTTGCATGCC (the reverse complement of the end).
        now.add_primer(
            "r",
            Seq::new(Alphabet::Dna, "AAGCTTGCATGCCTAAGG").unwrap(),
            "made here",
        )
        .unwrap();
        let map = [(then.primers()[0].id, 0)];
        let (text, owners) = primers_packet(&now, &then, Some(PACKET), &map).unwrap();
        // The edited one keeps its other attributes, in place, and gets fresh sites.
        assert!(text.contains("<Primer recentID=\"1\" name=\"f &amp; g\" sequence=\"ACGTTGCAAGGCTTAC\" description=\"&lt;html&gt;&lt;body&gt;&lt;/body&gt;&lt;/html&gt;\" dateAdded=\"2024-04-11T11:06:53Z\"><BindingSite location=\"0-15\" boundStrand=\"0\""), "{text}");
        assert!(text.contains("<BindingSite simplified=\"1\" location=\"0-15\""));
        // The hidden one is kept; the new one comes last, with the next id.
        assert!(text.contains("<Primer recentID=\"2\" name=\"hidden\" sequence=\"GGGG\"/>"));
        assert!(text.contains("<Primer recentID=\"3\" name=\"r\" sequence=\"AAGCTTGCATGCCTAAGG\" description=\"&lt;html&gt;&lt;body&gt;made here&lt;/body&gt;&lt;/html&gt;\"><BindingSite location=\"24-41\" boundStrand=\"1\""), "{text}");
        assert!(text.contains("nextValidID=\"4\""));
        assert!(text.ends_with("</Primer></Primers>\n"));
        assert_eq!(owners.len(), 3);
        // Deleting.
        let mut gone = then.clone();
        gone.remove_primer(then.primers()[0].id).unwrap();
        let (text, _) = primers_packet(&gone, &then, Some(PACKET), &map).unwrap();
        assert!(!text.contains("name=\"f\""));
        // A file without a primers packet.
        let (text, _) = primers_packet(
            &now,
            &Molecule::new("e", now.seq().clone(), now.topology()).unwrap(),
            None,
            &[],
        )
        .unwrap();
        assert!(
            text.starts_with("<?xml version=\"1.0\"?><Primers nextValidID=\"2\">"),
            "{text}"
        );
    }

    #[test]
    fn notes_change_only_their_lines() {
        let raw = "<Notes>\n<UUID>abc</UUID>\n<Created UTC=\"14:57:19\">2026.3.6</Created>\n<Organism>E. coli</Organism>\n<Description>&lt;html>&lt;body>old&lt;/body>&lt;/html></Description>\n</Notes>\n";
        let mut then = Notes {
            description: "<html><body>old</body></html>".into(),
            ..Notes::default()
        };
        then.fields.insert("UUID".into(), "abc".into());
        then.fields.insert("Created".into(), "2026.3.6".into());
        then.fields.insert("Organism".into(), "E. coli".into());
        assert_eq!(notes_packet(&then, &then, Some(raw)).unwrap(), None);
        let mut now = then.clone();
        now.description = "new <b> & more".into();
        now.fields.remove("Organism");
        now.fields.insert("Created".into(), "2026.3.7".into());
        now.fields.insert("CreatedBy".into(), "Lab".into());
        let out = notes_packet(&now, &then, Some(raw)).unwrap().unwrap();
        assert_eq!(
            out,
            "<Notes>\n<UUID>abc</UUID>\n<Created UTC=\"14:57:19\">2026.3.7</Created>\n<Description>&lt;html>&lt;body>new &amp;lt;b&amp;gt; &amp;amp; more&lt;/body>&lt;/html></Description>\n<CreatedBy>Lab</CreatedBy>\n</Notes>\n"
        );
        // A file without notes.
        let out = notes_packet(&now, &Notes::default(), None)
            .unwrap()
            .unwrap();
        assert!(out.starts_with("<Notes>\n") && out.ends_with("</Notes>\n"));
        assert!(out.contains("<CreatedBy>Lab</CreatedBy>"));
    }
}
