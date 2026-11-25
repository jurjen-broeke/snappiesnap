//! GFF3 (Generic Feature Format version 3), as specified by the Sequence Ontology project
//! ("GFF3 specification", version 1.26).
//!
//! A file is a list of lines: directives (`##...`), comments (`#...`), blank lines, and
//! feature rows of nine tab-separated columns (seqid, source, type, start, end, score, strand,
//! phase, attributes). An optional `##FASTA` section at the end holds sequences. Coordinates
//! are 1-based and inclusive; a feature crossing the origin of a circular sequence has an end
//! beyond the sequence length. Attribute values are percent-encoded where they contain
//! `;`, `=`, `&`, `,`, tabs, line breaks or `%`.
//!
//! Every line keeps its raw text, so an unchanged file is written back byte for byte, and an
//! edit rewrites only the rows it touches.

use std::fmt::Write as _;

/// Why a GFF3 file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Gff3Error {
    #[error("not a GFF3 file: it must start with a \"##gff-version 3\" line")]
    NotGff3,

    #[error("the file isn't valid UTF-8 text")]
    NotText,

    #[error("line {line}: {reason}")]
    Malformed { line: usize, reason: String },
}

/// One feature row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Row {
    pub seqid: String,
    pub source: String,
    /// The feature type: a Sequence Ontology term or accession, such as `gene` or `CDS`.
    pub kind: String,
    /// 1-based, inclusive.
    pub start: u64,
    pub end: u64,
    /// `.` when not given.
    pub score: String,
    /// `+`, `-`, `.` (not stranded) or `?` (unknown).
    pub strand: char,
    /// `0`, `1`, `2` for CDS rows, else `.`.
    pub phase: String,
    /// Tag and values, in file order, decoded.
    pub attributes: Vec<(String, Vec<String>)>,
}

impl Row {
    pub fn attribute(&self, tag: &str) -> Option<&str> {
        self.attributes
            .iter()
            .find(|(t, _)| t == tag)
            .and_then(|(_, v)| v.first())
            .map(String::as_str)
    }

    /// The row as a line (with the newline), encoded.
    pub fn to_line(&self) -> String {
        let attributes = if self.attributes.is_empty() {
            ".".to_string()
        } else {
            self.attributes
                .iter()
                .map(|(t, v)| {
                    let values: Vec<String> = v.iter().map(|x| encode(x, true)).collect();
                    format!("{}={}", encode(t, true), values.join(","))
                })
                .collect::<Vec<_>>()
                .join(";")
        };
        let dot = |s: &str| {
            if s.is_empty() {
                ".".to_string()
            } else {
                encode(s, false)
            }
        };
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            dot(&self.seqid),
            dot(&self.source),
            dot(&self.kind),
            self.start,
            self.end,
            dot(&self.score),
            self.strand,
            dot(&self.phase),
            attributes
        )
    }
}

/// One line of the annotation part.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Line {
    /// A `##` directive, comment or blank line, kept as written.
    Other(String),
    /// A feature row; `raw` is its original text while unchanged.
    Feature { row: Row, raw: Option<String> },
}

/// A sequence of the `##FASTA` section.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sequence {
    /// The first word of the `>` header (the seqid it belongs to).
    pub id: String,
    pub header: String,
    pub residues: Vec<u8>,
    /// The original text while unchanged.
    raw: Option<String>,
}

/// A GFF3 file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gff3File {
    pub lines: Vec<Line>,
    /// The `##FASTA` line, if the file has a sequence section.
    fasta_line: Option<String>,
    /// Text in the sequence section before the first `>` (normally nothing).
    fasta_preamble: String,
    pub sequences: Vec<Sequence>,
}

/// Characters that must be escaped: in attributes `;=&,`, everywhere tab, line breaks, `%`
/// and other control characters.
fn encode(s: &str, attribute: bool) -> String {
    let mut out = String::new();
    for c in s.chars() {
        let escape =
            c == '%' || c.is_control() || (attribute && matches!(c, ';' | '=' | '&' | ','));
        if escape {
            let mut buf = [0u8; 4];
            for b in c.encode_utf8(&mut buf).bytes() {
                let _ = write!(out, "%{b:02X}");
            }
        } else {
            out.push(c);
        }
    }
    out
}

/// Undo `%XX` escapes (invalid ones are kept as written).
fn decode(s: &str) -> String {
    let bytes = s.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0usize;
    while let Some(&b) = bytes.get(i) {
        let hex = bytes
            .get(i + 1..i + 3)
            .and_then(|h| std::str::from_utf8(h).ok())
            .and_then(|h| u8::from_str_radix(h, 16).ok());
        match (b, hex) {
            (b'%', Some(v)) => {
                out.push(v);
                i += 3;
            }
            _ => {
                out.push(b);
                i += 1;
            }
        }
    }
    String::from_utf8_lossy(&out).into_owned()
}

fn malformed(line: usize, reason: impl Into<String>) -> Gff3Error {
    Gff3Error::Malformed {
        line: line + 1,
        reason: reason.into(),
    }
}

fn parse_row(text: &str, line: usize) -> Result<Row, Gff3Error> {
    let cols: Vec<&str> = text.trim_end_matches(['\n', '\r']).split('\t').collect();
    let [
        seqid,
        source,
        kind,
        start,
        end,
        score,
        strand,
        phase,
        attributes,
    ] = cols.as_slice()
    else {
        return Err(malformed(
            line,
            format!("expected 9 tab-separated columns, found {}", cols.len()),
        ));
    };
    let number = |s: &str, what: &str| {
        s.trim()
            .parse::<u64>()
            .map_err(|_| malformed(line, format!("{what} {s:?} isn't a number")))
    };
    let (start, end) = (number(start, "start")?, number(end, "end")?);
    if start == 0 || end < start {
        return Err(malformed(
            line,
            format!("start {start} and end {end}: positions are 1-based, with start <= end"),
        ));
    }
    let strand = match *strand {
        "+" => '+',
        "-" => '-',
        "." => '.',
        "?" => '?',
        other => {
            return Err(malformed(
                line,
                format!("strand {other:?} isn't +, -, . or ?"),
            ));
        }
    };
    let undot = |s: &str| if s == "." { String::new() } else { decode(s) };
    let mut attrs = Vec::new();
    if *attributes != "." && !attributes.trim().is_empty() {
        for part in attributes.split(';').filter(|p| !p.trim().is_empty()) {
            let (tag, values) = part
                .split_once('=')
                .ok_or_else(|| malformed(line, format!("attribute {part:?} has no '='")))?;
            attrs.push((decode(tag.trim()), values.split(',').map(decode).collect()));
        }
    }
    Ok(Row {
        seqid: undot(seqid),
        source: undot(source),
        kind: undot(kind),
        start,
        end,
        score: undot(score),
        strand,
        phase: undot(phase),
        attributes: attrs,
    })
}

/// Sequence lines for new or changed sequences: 60 residues per line.
fn fasta_text(header: &str, residues: &[u8]) -> String {
    let mut out = format!(">{header}\n");
    for line in residues.chunks(60) {
        out.push_str(&String::from_utf8_lossy(line));
        out.push('\n');
    }
    out
}

impl Gff3File {
    /// An empty file (just the version line).
    pub fn new() -> Self {
        Gff3File {
            lines: vec![Line::Other("##gff-version 3\n".into())],
            fasta_line: None,
            fasta_preamble: String::new(),
            sequences: Vec::new(),
        }
    }

    pub fn parse(bytes: &[u8]) -> Result<Self, Gff3Error> {
        let text = std::str::from_utf8(bytes).map_err(|_| Gff3Error::NotText)?;
        if !text.starts_with("##gff-version 3") {
            return Err(Gff3Error::NotGff3);
        }
        let mut lines = Vec::new();
        let mut fasta_line = None;
        let mut rest = text.split_inclusive('\n').enumerate();
        for (n, l) in rest.by_ref() {
            let trimmed = l.trim_end_matches(['\n', '\r']);
            if trimmed == "##FASTA" {
                fasta_line = Some(l.to_string());
                break;
            }
            if trimmed.starts_with('#') || trimmed.trim().is_empty() {
                lines.push(Line::Other(l.to_string()));
            } else {
                lines.push(Line::Feature {
                    row: parse_row(l, n)?,
                    raw: Some(l.to_string()),
                });
            }
        }
        let mut fasta_preamble = String::new();
        let mut sequences: Vec<Sequence> = Vec::new();
        for (n, l) in rest {
            if let Some(header) = l.strip_prefix('>') {
                let header = header.trim_end_matches(['\n', '\r']).to_string();
                let id = header.split_whitespace().next().unwrap_or("").to_string();
                if id.is_empty() {
                    return Err(malformed(n, "sequence header without a name"));
                }
                sequences.push(Sequence {
                    id,
                    header,
                    residues: Vec::new(),
                    raw: Some(l.to_string()),
                });
            } else if let Some(s) = sequences.last_mut() {
                s.residues
                    .extend(l.bytes().filter(|b| !b.is_ascii_whitespace()));
                if let Some(raw) = s.raw.as_mut() {
                    raw.push_str(l);
                }
            } else {
                fasta_preamble.push_str(l);
            }
        }
        Ok(Gff3File {
            lines,
            fasta_line,
            fasta_preamble,
            sequences,
        })
    }

    pub fn to_text(&self) -> String {
        let mut out = String::new();
        for l in &self.lines {
            match l {
                Line::Other(raw) => out.push_str(raw),
                Line::Feature { raw: Some(raw), .. } => out.push_str(raw),
                Line::Feature { row, raw: None } => out.push_str(&row.to_line()),
            }
        }
        if self.fasta_line.is_some() || !self.sequences.is_empty() {
            out.push_str(self.fasta_line.as_deref().unwrap_or("##FASTA\n"));
            out.push_str(&self.fasta_preamble);
            for s in &self.sequences {
                match &s.raw {
                    Some(raw) => out.push_str(raw),
                    None => out.push_str(&fasta_text(&s.header, &s.residues)),
                }
            }
        }
        out
    }

    /// The sequence for `seqid`, if the file has one.
    pub fn sequence(&self, seqid: &str) -> Option<&Sequence> {
        self.sequences.iter().find(|s| s.id == seqid)
    }

    /// Set (or add) the sequence for `seqid`.
    pub fn set_sequence(&mut self, seqid: &str, residues: Vec<u8>) {
        match self.sequences.iter_mut().find(|s| s.id == seqid) {
            Some(s) if s.residues == residues => {}
            Some(s) => {
                s.residues = residues;
                s.raw = None;
            }
            None => self.sequences.push(Sequence {
                id: seqid.to_string(),
                header: seqid.to_string(),
                residues,
                raw: None,
            }),
        }
    }

    /// The seqids of the feature rows, in order of first appearance.
    pub fn seqids(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for l in &self.lines {
            if let Line::Feature { row, .. } = l
                && !out.contains(&row.seqid)
            {
                out.push(row.seqid.clone());
            }
        }
        out
    }
}

impl Default for Gff3File {
    fn default() -> Self {
        Gff3File::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &str = "##gff-version 3
##sequence-region ctg1 1 30
# a comment
ctg1\t.\tregion\t1\t30\t.\t+\t.\tID=ctg1;Is_circular=true
ctg1\tlab\tCDS\t25\t36\t.\t+\t0\tID=cds1;Name=wrap;Note=crosses the origin%2C on purpose
ctg1\tlab\tpromoter\t3\t9\t.\t-\t.\tID=p1;Name=P%3B1;Dbxref=a:1,b:2

##FASTA
>ctg1 a made-up contig
ACGTACGTACGTACGTACGTACGTACGTAC
";

    #[test]
    fn reads_and_writes_back_unchanged() {
        let f = Gff3File::parse(SAMPLE.as_bytes()).unwrap();
        let rows: Vec<&Row> = f
            .lines
            .iter()
            .filter_map(|l| match l {
                Line::Feature { row, .. } => Some(row),
                Line::Other(_) => None,
            })
            .collect();
        assert_eq!(rows.len(), 3);
        assert_eq!((rows[1].start, rows[1].end, rows[1].strand), (25, 36, '+'));
        assert_eq!(
            rows[1].attribute("Note"),
            Some("crosses the origin, on purpose")
        );
        assert_eq!(rows[2].attribute("Name"), Some("P;1"));
        assert_eq!(
            rows[2].attributes[2].1,
            vec!["a:1".to_string(), "b:2".to_string()]
        );
        assert_eq!(f.sequence("ctg1").unwrap().residues.len(), 30);
        assert_eq!(f.seqids(), vec!["ctg1".to_string()]);
        assert_eq!(f.to_text(), SAMPLE);
    }

    #[test]
    fn edited_rows_are_encoded_again() {
        let mut f = Gff3File::parse(SAMPLE.as_bytes()).unwrap();
        if let Some(Line::Feature { row, raw }) = f.lines.get_mut(5) {
            row.attributes[1].1 = vec!["P;2".into()];
            *raw = None;
        }
        f.set_sequence("ctg1", b"ACGT".to_vec());
        let text = f.to_text();
        assert!(
            text.contains("ctg1\tlab\tpromoter\t3\t9\t.\t-\t.\tID=p1;Name=P%3B2;Dbxref=a:1,b:2\n")
        );
        assert!(text.ends_with("##FASTA\n>ctg1 a made-up contig\nACGT\n"));
        // Untouched rows keep their text.
        assert!(text.contains("Note=crosses the origin%2C on purpose\n"));
        assert_eq!(Gff3File::parse(text.as_bytes()).unwrap().to_text(), text);
    }

    #[test]
    fn hostile_input_is_an_error() {
        for bad in [
            "",
            "##gff-version 2\n",
            "##gff-version 3\nctg1\t.\tgene\t1\n",
            "##gff-version 3\nctg1\t.\tgene\t0\t5\t.\t+\t.\t.\n",
            "##gff-version 3\nctg1\t.\tgene\t9\t5\t.\t+\t.\t.\n",
            "##gff-version 3\nctg1\t.\tgene\tx\t5\t.\t+\t.\t.\n",
            "##gff-version 3\nctg1\t.\tgene\t1\t5\t.\t*\t.\t.\n",
            "##gff-version 3\nctg1\t.\tgene\t1\t5\t.\t+\t.\tnovalue\n",
            "##gff-version 3\n##FASTA\n> \nACGT\n",
        ] {
            assert!(Gff3File::parse(bad.as_bytes()).is_err(), "{bad:?}");
        }
        // Escapes that aren't hex stay as written; huge numbers are errors, not panics.
        assert_eq!(decode("100%zz"), "100%zz");
        assert!(
            Gff3File::parse(b"##gff-version 3\nc\t.\tg\t1\t99999999999999999999\t.\t+\t.\t.\n")
                .is_err()
        );
    }
}
