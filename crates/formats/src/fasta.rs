//! FASTA: `>header` lines followed by sequence lines.
//!
//! [`FastaFile`] keeps the original text, so an unchanged file is written back byte for byte.
//! New text uses one record per sequence, wrapped at a fixed width.

/// Why a FASTA file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum FastaError {
    #[error("not a FASTA file: it must start with a \">\" header line")]
    NotFasta,

    #[error("the file isn't valid UTF-8 text")]
    NotText,

    #[error("record {record} has no residues")]
    EmptyRecord { record: usize },
}

/// One `>header` and its residues (whitespace removed, case kept).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// The header line without the leading `>`, trimmed.
    pub header: String,
    pub residues: Vec<u8>,
}

impl FastaRecord {
    /// The record's name: the first word of the header.
    pub fn name(&self) -> &str {
        self.header.split_whitespace().next().unwrap_or("")
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaFile {
    pub records: Vec<FastaRecord>,
    /// Residues per line in the original (for writing in the same style).
    pub line_width: usize,
}

/// Default line width for new files (as the reference application exports).
pub const DEFAULT_WIDTH: usize = 80;

impl FastaFile {
    pub fn parse(bytes: &[u8]) -> Result<Self, FastaError> {
        let text = std::str::from_utf8(bytes).map_err(|_| FastaError::NotText)?;
        let text = text.strip_prefix('\u{feff}').unwrap_or(text);
        if !text.trim_start().starts_with('>') {
            return Err(FastaError::NotFasta);
        }
        let mut records: Vec<FastaRecord> = Vec::new();
        let mut line_width = 0usize;
        for line in text.lines() {
            let line = line.trim_end_matches('\r');
            if let Some(header) = line.strip_prefix('>') {
                records.push(FastaRecord {
                    header: header.trim().to_string(),
                    residues: Vec::new(),
                });
            } else if let Some(rec) = records.last_mut() {
                // `;` lines are old-style comments.
                if line.starts_with(';') {
                    continue;
                }
                let before = rec.residues.len();
                rec.residues
                    .extend(line.bytes().filter(|b| !b.is_ascii_whitespace()));
                line_width = line_width.max(rec.residues.len() - before);
            }
        }
        if let Some(i) = records.iter().position(|r| r.residues.is_empty()) {
            return Err(FastaError::EmptyRecord { record: i + 1 });
        }
        Ok(FastaFile {
            records,
            line_width: if line_width == 0 {
                DEFAULT_WIDTH
            } else {
                line_width
            },
        })
    }

    /// The file as text: `>header` then residues wrapped at `line_width`, LF line endings.
    pub fn to_text(&self) -> String {
        let width = self.line_width.max(1);
        let mut out = String::new();
        for r in &self.records {
            out.push('>');
            out.push_str(&r.header);
            out.push('\n');
            for chunk in r.residues.chunks(width) {
                out.push_str(&String::from_utf8_lossy(chunk));
                out.push('\n');
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_rewrites() {
        let text = ">pX  (12 bp)\nACGTACGTAC\nGT\n>second desc\nacgu\n";
        let f = FastaFile::parse(text.as_bytes()).unwrap();
        assert_eq!(f.records.len(), 2);
        assert_eq!(f.records[0].name(), "pX");
        assert_eq!(f.records[0].residues, b"ACGTACGTACGT");
        assert_eq!(f.records[1].residues, b"acgu");
        assert_eq!(f.line_width, 10);
        assert_eq!(f.to_text(), text);
    }

    #[test]
    fn tolerates_crlf_bom_blank_lines_and_comments() {
        let text = "\u{feff}>a\r\nAC GT\r\n\r\n;comment\r\nTT\r\n";
        let f = FastaFile::parse(text.as_bytes()).unwrap();
        assert_eq!(f.records[0].residues, b"ACGTTT");
    }

    #[test]
    fn rejects_non_fasta() {
        assert_eq!(FastaFile::parse(b"ACGT"), Err(FastaError::NotFasta));
        assert_eq!(FastaFile::parse(b""), Err(FastaError::NotFasta));
        assert_eq!(FastaFile::parse(b"\xff>"), Err(FastaError::NotText));
        assert_eq!(
            FastaFile::parse(b">a\nAC\n>b\n"),
            Err(FastaError::EmptyRecord { record: 2 })
        );
    }
}
