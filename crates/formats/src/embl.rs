//! EMBL flat files (ENA/INSDC), as described in the ENA "EMBL Flat File Format" user manual.
//!
//! Each line starts with a two-letter code (`ID`, `AC`, `DE`, `FT`, `SQ`, ...); content starts
//! in column 6. The feature table (`FT` lines) has the same layout and syntax as GenBank's
//! (INSDC feature table definition), so its entries are parsed into [`FeatureEntry`] values
//! and written with the same rules, with `FT` in the first two columns. Sequence lines hold 60
//! residues in groups of ten, with the running count right-aligned to column 80.
//!
//! Like the GenBank reader, records keep their raw text (ID line, header, each feature, the
//! sequence), so unchanged parts are written back byte for byte.

use std::fmt::Write as _;

use crate::genbank::{Feature, FeatureEntry, QUALIFIER_INDENT, parse_entry};

/// Why an EMBL file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum EmblError {
    #[error("not an EMBL file: it must start with an ID line")]
    NotEmbl,

    #[error("the file isn't valid UTF-8 text")]
    NotText,

    #[error("line {line}: {reason}")]
    Malformed { line: usize, reason: String },
}

/// The fields of an ID line that the document uses.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdLine {
    /// The primary accession or entry name (the first word).
    pub name: String,
    pub circular: bool,
    /// The molecule type, such as `genomic DNA` or `other DNA`.
    pub molecule: String,
    /// The length in base pairs.
    pub length: u64,
}

impl IdLine {
    /// Read an ID line, in the current form (`ID   X56562; SV 1; circular; genomic DNA; STD;
    /// SYN; 2686 BP.`) or the pre-2006 one (`ID   X56562  standard; circular DNA; SYN; 2686
    /// BP.`).
    pub fn parse(line: &str) -> Option<IdLine> {
        let body = line.strip_prefix("ID")?.trim();
        let tokens: Vec<&str> = body
            .split(';')
            .map(str::trim)
            .filter(|t| !t.is_empty())
            .collect();
        let name = tokens.first()?.split_whitespace().next()?.to_string();
        let circular = tokens.iter().any(|t| {
            t.split_whitespace()
                .any(|w| w.eq_ignore_ascii_case("circular"))
        });
        let molecule = tokens
            .iter()
            .find(|t| t.contains("DNA") || t.contains("RNA"))
            .map(|t| {
                t.split_whitespace()
                    .filter(|w| {
                        !w.eq_ignore_ascii_case("circular") && !w.eq_ignore_ascii_case("linear")
                    })
                    .collect::<Vec<_>>()
                    .join(" ")
            })
            .unwrap_or_default();
        let length = tokens
            .iter()
            .rev()
            .find_map(|t| {
                let t = t.trim_end_matches('.');
                t.strip_suffix("BP").map(str::trim)
            })
            .and_then(|n| n.parse().ok())
            .unwrap_or(0);
        Some(IdLine {
            name,
            circular,
            molecule,
            length,
        })
    }

    /// The ID line in the current layout, ending with a newline.
    pub fn to_line(&self) -> String {
        let molecule = if self.molecule.is_empty() {
            "other DNA"
        } else {
            &self.molecule
        };
        format!(
            "ID   {}; SV 1; {}; {molecule}; STD; SYN; {} BP.\n",
            self.name,
            if self.circular { "circular" } else { "linear" },
            self.length
        )
    }
}

/// Patch the topology and length of an ID line in place, keeping everything else; `None` if
/// the line doesn't have the expected tokens.
fn patch_id_line(raw: &str, old: &IdLine, new: &IdLine) -> Option<String> {
    if old.name != new.name || old.molecule != new.molecule {
        return None;
    }
    let mut out = raw.to_string();
    if old.circular != new.circular {
        let (from, to) = if new.circular {
            ("linear", "circular")
        } else {
            ("circular", "linear")
        };
        let at = out.find(from)?;
        out.replace_range(at..at + from.len(), to);
    }
    if old.length != new.length {
        let from = format!(" {} BP.", old.length);
        let at = out.rfind(&from)?;
        out.replace_range(at..at + from.len(), &format!(" {} BP.", new.length));
    }
    Some(out)
}

/// One record (ID ... //).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EmblRecord {
    pub id: IdLine,
    id_raw: Option<String>,
    /// Lines after ID up to the feature table (AC, DE, KW, OS, references, XX, ...).
    pub header: String,
    /// The FH lines heading the feature table.
    fh: String,
    pub features: Vec<Feature>,
    /// Lines between the feature table and SQ (XX, CO, ...).
    pub middle: String,
    /// The SQ line (empty when the record had none), while the residues are unchanged.
    sq_raw: Option<String>,
    pub residues: Vec<u8>,
    residues_raw: Option<String>,
    /// The `//` line and anything up to the next record.
    pub tail: String,
}

impl EmblRecord {
    /// A new record with a minimal header.
    pub fn new(id: IdLine, header: String) -> Self {
        EmblRecord {
            id,
            id_raw: None,
            header,
            fh: "FH   Key             Location/Qualifiers\nFH\n".into(),
            features: Vec::new(),
            middle: "XX\n".into(),
            sq_raw: None,
            residues: Vec::new(),
            residues_raw: None,
            tail: "//\n".into(),
        }
    }

    /// Change the ID fields; the original line is patched in place when only the length or
    /// topology changed, so its layout is kept.
    pub fn set_id(&mut self, id: IdLine) {
        if id == self.id {
            return;
        }
        let patched = self
            .id_raw
            .as_deref()
            .and_then(|raw| patch_id_line(raw, &self.id, &id));
        self.id = id;
        self.id_raw = patched;
    }

    pub fn set_residues(&mut self, residues: Vec<u8>) {
        if residues != self.residues {
            self.residues = residues;
            self.residues_raw = None;
            self.sq_raw = None;
        }
    }

    pub fn to_text(&self) -> String {
        let mut out = String::new();
        match &self.id_raw {
            Some(raw) => out.push_str(raw),
            None => out.push_str(&self.id.to_line()),
        }
        out.push_str(&self.header);
        if !self.features.is_empty() || !self.fh.is_empty() {
            out.push_str(&self.fh);
        }
        for f in &self.features {
            match &f.raw {
                Some(raw) => out.push_str(raw),
                None => out.push_str(&ft_text(&f.entry)),
            }
        }
        out.push_str(&self.middle);
        match &self.sq_raw {
            Some(raw) => out.push_str(raw),
            None => out.push_str(&sq_line(&self.residues)),
        }
        match &self.residues_raw {
            Some(raw) => out.push_str(raw),
            None => out.push_str(&sequence_text(&self.residues)),
        }
        out.push_str(&self.tail);
        out
    }
}

/// A feature entry as FT lines.
pub fn ft_text(entry: &FeatureEntry) -> String {
    entry
        .to_text()
        .lines()
        .map(|l| format!("FT{}\n", l.get(2..).unwrap_or("")))
        .collect()
}

/// The SQ line: length and base counts.
pub fn sq_line(residues: &[u8]) -> String {
    let count = |b: u8| {
        residues
            .iter()
            .filter(|r| r.eq_ignore_ascii_case(&b))
            .count()
    };
    let (a, c, g, t) = (
        count(b'a'),
        count(b'c'),
        count(b'g'),
        count(b't') + count(b'u'),
    );
    let other = residues.len().saturating_sub(a + c + g + t);
    format!(
        "SQ   Sequence {} BP; {a} A; {c} C; {g} G; {t} T; {other} other;\n",
        residues.len()
    )
}

/// Sequence lines: five spaces, 60 lower-case residues in groups of ten, and the count of
/// residues so far right-aligned to column 80.
pub fn sequence_text(residues: &[u8]) -> String {
    let mut out = String::new();
    let mut done = 0usize;
    for line in residues.chunks(60) {
        done += line.len();
        let groups: Vec<String> = line
            .chunks(10)
            .map(|g| String::from_utf8_lossy(g).to_ascii_lowercase())
            .collect();
        let _ = writeln!(out, "     {:<65}{done:>10}", groups.join(" "));
    }
    out
}

/// An EMBL file: one or more records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EmblFile {
    pub records: Vec<EmblRecord>,
}

impl EmblFile {
    pub fn to_text(&self) -> String {
        self.records.iter().map(EmblRecord::to_text).collect()
    }

    pub fn parse(bytes: &[u8]) -> Result<Self, EmblError> {
        let text = std::str::from_utf8(bytes).map_err(|_| EmblError::NotText)?;
        let lines: Vec<&str> = text.split_inclusive('\n').collect();
        if !lines.first().is_some_and(|l| l.starts_with("ID   ")) {
            return Err(EmblError::NotEmbl);
        }
        let mut records = Vec::new();
        let mut i = 0usize;
        while i < lines.len() {
            let (record, next) = parse_record(&lines, i)?;
            records.push(record);
            i = next;
        }
        Ok(EmblFile { records })
    }
}

fn malformed(line: usize, reason: impl Into<String>) -> EmblError {
    EmblError::Malformed {
        line: line + 1,
        reason: reason.into(),
    }
}

/// The two-letter code of a line.
fn code(line: &str) -> &str {
    line.get(..2).unwrap_or("")
}

/// Parse the record starting at line `start` (an ID line). Returns it and the index of the next
/// record's ID line (or the end).
fn parse_record(lines: &[&str], start: usize) -> Result<(EmblRecord, usize), EmblError> {
    let id_raw = lines.get(start).copied().unwrap_or("");
    let id = IdLine::parse(id_raw).ok_or_else(|| malformed(start, "unreadable ID line"))?;
    let mut i = start + 1;
    let at = |i: usize| lines.get(i).copied();
    let mut header = String::new();
    while let Some(l) = at(i) {
        if matches!(code(l), "FH" | "FT" | "SQ" | "//") || l.starts_with("ID   ") {
            break;
        }
        header.push_str(l);
        i += 1;
    }
    let mut fh = String::new();
    while let Some(l) = at(i).filter(|l| code(l) == "FH") {
        fh.push_str(l);
        i += 1;
    }
    // Feature entries: a line with a key in columns 6-20 starts one.
    let mut features = Vec::new();
    let mut current: Option<(usize, String)> = None;
    let finish = |entry: Option<(usize, String)>, features: &mut Vec<Feature>| {
        if let Some((line, raw)) = entry {
            let gb: String = raw
                .lines()
                .map(|l| format!("  {}\n", l.get(2..).unwrap_or("")))
                .collect();
            let entry = parse_entry(&gb).map_err(|reason| malformed(line, reason))?;
            features.push(Feature {
                entry,
                raw: Some(raw),
            });
        }
        Ok::<(), EmblError>(())
    };
    while let Some(l) = at(i).filter(|l| code(l) == "FT") {
        let key = l.get(5..QUALIFIER_INDENT).unwrap_or("").trim();
        if !key.is_empty() {
            finish(current.take(), &mut features)?;
            current = Some((i, l.to_string()));
        } else if let Some((_, raw)) = current.as_mut() {
            raw.push_str(l);
        } else {
            return Err(malformed(i, "feature table line before any feature key"));
        }
        i += 1;
    }
    finish(current.take(), &mut features)?;
    let mut middle = String::new();
    while let Some(l) = at(i) {
        if matches!(code(l), "SQ" | "//") || l.starts_with("ID   ") {
            break;
        }
        middle.push_str(l);
        i += 1;
    }
    // No SQ block: nothing is written for it until residues are set.
    let mut sq_raw = Some(String::new());
    let mut residues = Vec::new();
    let mut residues_raw = String::new();
    if let Some(l) = at(i).filter(|l| code(l) == "SQ") {
        sq_raw = Some(l.to_string());
        i += 1;
        while let Some(l) = at(i).filter(|l| l.starts_with(' ')) {
            residues.extend(l.bytes().filter(u8::is_ascii_alphabetic));
            residues_raw.push_str(l);
            i += 1;
        }
    }
    let mut tail = String::new();
    while let Some(l) = at(i) {
        if l.starts_with("ID   ") && tail.contains("//") {
            break;
        }
        if l.starts_with("ID   ") {
            return Err(malformed(
                i,
                "a record must end with // before the next ID line",
            ));
        }
        tail.push_str(l);
        i += 1;
    }
    if !tail.starts_with("//") {
        return Err(malformed(
            i.saturating_sub(1),
            "the record doesn't end with //",
        ));
    }
    Ok((
        EmblRecord {
            id,
            id_raw: Some(id_raw.to_string()),
            header,
            fh,
            features,
            middle,
            sq_raw,
            residues,
            residues_raw: Some(residues_raw),
            tail,
        },
        i,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "ID   TEST01; SV 1; circular; other DNA; STD; SYN; 70 BP.
XX
AC   TEST01;
XX
DE   A made-up construct.
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..70
FT                   /organism=\"synthetic construct\"
FT                   /mol_type=\"other DNA\"
FT   CDS             complement(join(60..70,1..9))
FT                   /gene=\"x\"
FT                   /note=\"a long note that runs on and on so that it has to wrap
FT                   onto a second line\"
XX
SQ   Sequence 70 BP; 18 A; 17 C; 18 G; 17 T; 0 other;
     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60
     acgtacgtac                                                               70
//
";

    #[test]
    fn reads_and_writes_back_unchanged() {
        let f = EmblFile::parse(SAMPLE.as_bytes()).unwrap();
        let r = &f.records[0];
        assert_eq!(
            r.id,
            IdLine {
                name: "TEST01".into(),
                circular: true,
                molecule: "other DNA".into(),
                length: 70
            }
        );
        assert_eq!(r.residues.len(), 70);
        assert_eq!(r.features.len(), 2);
        assert_eq!(r.features[1].entry.key, "CDS");
        assert_eq!(
            r.features[1].entry.location,
            "complement(join(60..70,1..9))"
        );
        assert_eq!(
            r.features[1].entry.qualifier("note"),
            // Wrapped as the GenBank reader reads it (the same feature table rules).
            Some("a long note that runs on and on so that it has to wrap\nonto a second line")
        );
        assert_eq!(f.to_text(), SAMPLE);
    }

    #[test]
    fn edits_rewrite_only_what_changed() {
        let mut f = EmblFile::parse(SAMPLE.as_bytes()).unwrap();
        let r = &mut f.records[0];
        r.features[1].entry.location = "complement(join(61..70,1..9))".into();
        r.features[1].raw = None;
        let mut id = r.id.clone();
        id.circular = false;
        r.set_id(id);
        let text = f.to_text();
        assert!(text.starts_with("ID   TEST01; SV 1; linear; other DNA; STD; SYN; 70 BP.\n"));
        assert!(text.contains("FT   CDS             complement(join(61..70,1..9))\n"));
        // The untouched source feature and the sequence are as they were.
        assert!(text.contains("FT   source          1..70\nFT                   /organism"));
        assert!(text.contains(
            "acgtacgtac                                                               70\n"
        ));
        let again = EmblFile::parse(text.as_bytes()).unwrap();
        assert_eq!(
            again.records[0].features[1].entry.qualifiers,
            f.records[0].features[1].entry.qualifiers
        );
    }

    #[test]
    fn new_sequence_lines_follow_the_layout() {
        let residues = b"ACGT".repeat(16);
        assert_eq!(
            sequence_text(&residues),
            "     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt        60\n     acgt                                                                     64\n"
        );
        assert_eq!(
            sq_line(&residues),
            "SQ   Sequence 64 BP; 16 A; 16 C; 16 G; 16 T; 0 other;\n"
        );
    }

    #[test]
    fn a_record_without_sequence_stays_without_one() {
        // Found by the property test: writing back used to add an SQ line.
        let text = "ID   P; SV 1; circular; other DNA; STD; SYN; 10 BP.\nXX\n//\n";
        assert_eq!(EmblFile::parse(text.as_bytes()).unwrap().to_text(), text);
    }

    #[test]
    fn old_style_id_lines() {
        let id = IdLine::parse("ID   X56562     standard; circular DNA; SYN; 2686 BP.\n").unwrap();
        assert_eq!(
            (id.name.as_str(), id.circular, id.length),
            ("X56562", true, 2686)
        );
        assert_eq!(id.molecule, "DNA");
    }

    #[test]
    fn hostile_input_is_an_error() {
        for bad in [
            "",
            "LOCUS x",
            "ID   \n",
            "ID   A; linear; DNA; 5 BP.\nFT                   /note=\"x\"\n//\n",
            "ID   A; linear; DNA; 5 BP.\nSQ   Sequence 5 BP;\n     acgta 5\n",
            "ID   A; linear; DNA; 5 BP.\nID   B; linear; DNA; 5 BP.\n//\n",
        ] {
            assert!(EmblFile::parse(bad.as_bytes()).is_err(), "{bad:?}");
        }
    }
}
