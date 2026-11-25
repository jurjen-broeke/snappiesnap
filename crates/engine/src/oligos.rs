//! The lab's oligo library: oligos (sequencing primers and the like) imported from a list once
//! and kept in the data folder ([`LIBRARY_FILE`]), so any construct can be searched for the ones
//! that bind it (`oligo_library_find`) and those added as its primers.
//!
//! Lists are read as plain text ([`parse_oligo_text`]): one oligo per line as name, sequence
//! and an optional description separated by tabs (the layout of a primer list exported from
//! SnapGene or a spreadsheet), or FASTA. A file holding nothing but one sequence has its name
//! nowhere but in the file name, so it is named after the file. Sequences keep their case,
//! since labs often write a 5' tail in lower case.

use serde::{Deserialize, Serialize};

/// The library, inside [`crate::Platform::data_dir`].
pub const LIBRARY_FILE: &str = "oligo_library.json";
/// Most oligos in the library.
pub const MAX_OLIGOS: usize = 5_000;
/// Longest oligo accepted, in bases.
pub const MAX_LENGTH: usize = 200;
/// Shortest oligo accepted, in bases (shorter ones can't be placed on a template).
pub const MIN_LENGTH: usize = 10;
/// Longest name or description, in characters.
const MAX_TEXT: usize = 200;
/// Largest list file read.
pub const MAX_FILE: usize = 8 << 20;
/// Most problems listed in a report (the count is always given).
const MAX_REPORTED: usize = 50;

/// One oligo of the library.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LabOligo {
    pub name: String,
    /// 5'->3', IUPAC nucleotide codes, case as given.
    pub residues: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub description: String,
}

impl LabOligo {
    /// Whether two oligos have the same sequence (ignoring case).
    pub fn same_sequence(&self, residues: &str) -> bool {
        self.residues.eq_ignore_ascii_case(residues)
    }
}

/// The library as stored.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OligoLibrary {
    pub oligos: Vec<LabOligo>,
}

/// A line that wasn't read, and why.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Skipped {
    /// 1-based line number in the file.
    pub line: usize,
    pub reason: String,
    /// The file the line is in, when a whole folder was imported.
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub file: String,
}

impl Skipped {
    fn new(line: usize, reason: impl Into<String>) -> Self {
        Skipped {
            line,
            reason: reason.into(),
            file: String::new(),
        }
    }

    /// Note which file of a folder import the line came from.
    pub fn in_file(mut self, file: &str) -> Self {
        self.file = file.to_string();
        self
    }
}

/// What reading a list (and adding it to the library) found.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct ImportReport {
    /// Oligos added to the library.
    pub added: usize,
    /// Oligos already there under the same name with the same sequence (or repeated in the
    /// file), left as they were.
    pub duplicates: usize,
    /// Names already used for a different sequence: the oligo already there is kept.
    pub conflicts: Vec<String>,
    /// Lines that aren't oligos (at most 50 listed; `skipped_count` has them all).
    pub skipped: Vec<Skipped>,
    pub skipped_count: usize,
}

/// Read an oligo list: tab-separated lines (name, sequence, optional description; a header
/// line such as "Name<TAB>Sequence" is passed over) or FASTA (`>name description` then the
/// sequence). UTF-8 or UTF-16. `named` is the name to give an oligo the file doesn't name,
/// normally the file's own name without its extension: labs keep one primer per file that way,
/// with the name only in the file name. Returns the oligos in file order and the lines that
/// weren't oligos. Errors only when nothing at all could be read.
pub fn parse_oligo_text(
    bytes: &[u8],
    named: Option<&str>,
) -> Result<(Vec<LabOligo>, Vec<Skipped>), String> {
    if bytes.len() > MAX_FILE {
        return Err(format!(
            "the file is larger than {} MB; an oligo list is a few lines per oligo",
            MAX_FILE >> 20
        ));
    }
    let text = crate::text::decode_text(bytes);
    let (oligos, skipped) = if text.trim_start().starts_with('>') {
        parse_fasta(&text)
    } else if let Some(bare) = parse_bare(&text, named) {
        bare
    } else {
        parse_table(&text)
    };
    if oligos.is_empty() {
        let why = skipped.first().map_or(String::new(), |s| {
            format!(" (line {}: {})", s.line, s.reason)
        });
        return Err(format!("the file holds no oligos{why}"));
    }
    Ok((oligos, skipped))
}

/// A file holding one bare sequence and nothing else, named after the file (`named`). `None`
/// when the file isn't of that shape, so it is read as a table instead.
///
/// Only a single line counts. Joining several bare lines into one oligo would invent a
/// sequence that is in no file, which matters more than reading a wrapped one: a list of
/// sequences without names has to say which name belongs to which, and FASTA is how a wrapped
/// sequence carries its name.
fn parse_bare(text: &str, named: Option<&str>) -> Option<(Vec<LabOligo>, Vec<Skipped>)> {
    let name = named?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let line = lines.next()?.trim();
    if lines.next().is_some() || line.contains('\t') {
        return None;
    }
    if !line
        .chars()
        .all(|c| "ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(c))
    {
        return None;
    }
    Some(match oligo(name, line, "") {
        Ok(o) => (vec![o], Vec::new()),
        Err(reason) => (Vec::new(), vec![Skipped::new(1, reason)]),
    })
}

fn parse_table(text: &str) -> (Vec<LabOligo>, Vec<Skipped>) {
    let mut oligos = Vec::new();
    let mut skipped = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line_no = i.saturating_add(1);
        if line.trim().is_empty() {
            continue;
        }
        let mut cells = line.split('\t').map(str::trim);
        let name = cells.next().unwrap_or_default();
        let residues = cells.next().unwrap_or_default();
        // Further cells: a description (empty cells from spreadsheet exports are dropped).
        let description = cells
            .filter(|c| !c.is_empty())
            .collect::<Vec<_>>()
            .join(" ");
        match oligo(name, residues, &description) {
            Ok(o) => oligos.push(o),
            // A first line that names columns is a header, not a bad oligo.
            Err(_) if oligos.is_empty() && skipped.is_empty() && is_header(name, residues) => {}
            Err(reason) => skipped.push(Skipped::new(line_no, reason)),
        }
    }
    (oligos, skipped)
}

fn parse_fasta(text: &str) -> (Vec<LabOligo>, Vec<Skipped>) {
    let mut oligos = Vec::new();
    let mut skipped = Vec::new();
    // The record being read: its header line number, name, description and sequence.
    let mut current: Option<(usize, String, String, String)> = None;
    let finish = |record: Option<(usize, String, String, String)>,
                  oligos: &mut Vec<LabOligo>,
                  skipped: &mut Vec<Skipped>| {
        if let Some((line, name, description, residues)) = record {
            match oligo(&name, &residues, &description) {
                Ok(o) => oligos.push(o),
                Err(reason) => skipped.push(Skipped::new(line, reason)),
            }
        }
    };
    for (i, line) in text.lines().enumerate() {
        let line_no = i.saturating_add(1);
        let line = line.trim();
        if let Some(header) = line.strip_prefix('>') {
            finish(current.take(), &mut oligos, &mut skipped);
            let header = header.trim();
            let (name, description) = header
                .split_once(char::is_whitespace)
                .map_or((header, ""), |(n, d)| (n, d.trim()));
            current = Some((
                line_no,
                name.to_string(),
                description.to_string(),
                String::new(),
            ));
        } else if let Some((_, _, _, residues)) = &mut current {
            if residues.len() <= MAX_LENGTH {
                residues.push_str(line);
            }
        } else if !line.is_empty() {
            skipped.push(Skipped::new(
                line_no,
                "sequence before the first >name line",
            ));
        }
    }
    finish(current.take(), &mut oligos, &mut skipped);
    (oligos, skipped)
}

/// "Name / Sequence"-style column titles.
fn is_header(name: &str, residues: &str) -> bool {
    let words = [
        "name", "oligo", "primer", "id", "sequence", "seq", "5'-3'", "5' to 3'",
    ];
    let is_title = |s: &str| {
        let s = s.trim().to_ascii_lowercase();
        words
            .iter()
            .any(|w| s == *w || s.starts_with(&format!("{w} ")))
    };
    is_title(name) || is_title(residues)
}

/// A checked oligo; errors say what is wrong, in words.
pub fn oligo(name: &str, residues: &str, description: &str) -> Result<LabOligo, String> {
    let name = name.trim();
    if name.is_empty() {
        return Err("no name".to_string());
    }
    if name.chars().count() > MAX_TEXT {
        return Err(format!("the name is longer than {MAX_TEXT} characters"));
    }
    let bases: String = residues.chars().filter(|c| !c.is_whitespace()).collect();
    if bases.is_empty() {
        return Err(format!("{name}: no sequence"));
    }
    if let Some(bad) = bases
        .chars()
        .find(|c| !"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(*c))
    {
        return Err(format!("{name}: {bad:?} is not a nucleotide code"));
    }
    let n = bases.len();
    if !(MIN_LENGTH..=MAX_LENGTH).contains(&n) {
        return Err(format!(
            "{name}: {n} bases (an oligo here is {MIN_LENGTH} to {MAX_LENGTH})"
        ));
    }
    Ok(LabOligo {
        name: name.to_string(),
        residues: bases,
        description: description.trim().chars().take(MAX_TEXT).collect(),
    })
}

impl OligoLibrary {
    /// Add `oligos`: one already there under the same name (ignoring case) with the same
    /// sequence is a duplicate; with another sequence a conflict, and the one already there is
    /// kept. If the library would grow past [`MAX_OLIGOS`], nothing is added and it errors.
    pub fn add(&mut self, oligos: Vec<LabOligo>) -> Result<ImportReport, String> {
        let mut report = ImportReport::default();
        let mut next = self.oligos.clone();
        for o in oligos {
            match next.iter().find(|x| x.name.eq_ignore_ascii_case(&o.name)) {
                Some(x) if x.same_sequence(&o.residues) => {
                    report.duplicates = report.duplicates.saturating_add(1);
                }
                Some(_) => {
                    if !report.conflicts.contains(&o.name) {
                        report.conflicts.push(o.name);
                    }
                }
                None => {
                    if next.len() >= MAX_OLIGOS {
                        return Err(format!(
                            "the library holds at most {MAX_OLIGOS} oligos; remove some first"
                        ));
                    }
                    next.push(o);
                    report.added = report.added.saturating_add(1);
                }
            }
        }
        self.oligos = next;
        Ok(report)
    }

    /// Check a library read from the data folder (it may have been edited by hand): invalid
    /// oligos and repeated names are dropped.
    pub fn checked(self) -> OligoLibrary {
        let mut out = OligoLibrary::default();
        let valid = self
            .oligos
            .into_iter()
            .filter_map(|o| oligo(&o.name, &o.residues, &o.description).ok())
            .take(MAX_OLIGOS)
            .collect();
        // Repeats are dropped quietly; nothing else can go wrong here.
        let _ = out.add(valid);
        out
    }
}

impl ImportReport {
    /// Note the lines a file couldn't use (all are counted, the first 50 listed).
    pub fn with_skipped(mut self, skipped: Vec<Skipped>) -> Self {
        self.skipped_count = skipped.len();
        self.skipped = skipped.into_iter().take(MAX_REPORTED).collect();
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tab_separated_lists() {
        // The layout of a list exported from SnapGene: name, tab, sequence, sometimes trailing
        // empty cells; Windows line ends; lower-case tails kept.
        let text = "FWD1\tacgttgcaggttacggaaga\r\nTAIL1\tgaattccATGACCGGTACCTAGCTAGC\t\t\r\n\r\nSEQ1\tACGTACGTACGT\tsequencing, M13 fwd\t\r\n";
        let (oligos, skipped) = parse_oligo_text(text.as_bytes(), None).unwrap();
        assert!(skipped.is_empty());
        assert_eq!(oligos.len(), 3);
        assert_eq!(oligos[1].residues, "gaattccATGACCGGTACCTAGCTAGC");
        assert_eq!(oligos[2].description, "sequencing, M13 fwd");
        // A header line is passed over; bad lines are reported with their line numbers.
        let text = "Name\tSequence\nok1\tACGTACGTACGT\nbad\tACGTXACGTACG\nshort\tACGT\n\tACGTACGTACGT\nonly name\n";
        let (oligos, skipped) = parse_oligo_text(text.as_bytes(), None).unwrap();
        assert_eq!(oligos.len(), 1);
        let lines: Vec<usize> = skipped.iter().map(|s| s.line).collect();
        assert_eq!(lines, vec![3, 4, 5, 6]);
        assert!(skipped[0].reason.contains("'X'"));
        assert!(skipped[1].reason.contains("4 bases"));
        // UTF-16 with a byte-order mark, as SnapGene exports on Windows.
        let mut utf16 = vec![0xFF, 0xFE];
        for u in "PROM1\tacgtggcaggtcgtatacg\r\n".encode_utf16() {
            utf16.extend(u.to_le_bytes());
        }
        assert_eq!(parse_oligo_text(&utf16, None).unwrap().0[0].name, "PROM1");
    }

    #[test]
    fn fasta_lists() {
        let text =
            ">M13F sequencing primer\nGTAAAACGAC\nGGCCAGT\n>bad\nACGT\n>T7\nTAATACGACTCACTATAGGG\n";
        let (oligos, skipped) = parse_oligo_text(text.as_bytes(), None).unwrap();
        assert_eq!(oligos.len(), 2);
        assert_eq!(
            (oligos[0].name.as_str(), oligos[0].residues.as_str()),
            ("M13F", "GTAAAACGACGGCCAGT")
        );
        assert_eq!(oligos[0].description, "sequencing primer");
        assert_eq!(skipped[0].line, 4);
    }

    #[test]
    fn one_bare_sequence_is_named_after_its_file() {
        // The lab's cloning primers: one file per primer, the name only in the file name.
        let (oligos, skipped) =
            parse_oligo_text(b"ggatccagatctacgtacgtacgttgcatgcatgcaag\n", Some("SP_201")).unwrap();
        assert!(skipped.is_empty());
        assert_eq!(oligos.len(), 1);
        assert_eq!(oligos[0].name, "SP_201");
        assert_eq!(oligos[0].residues, "ggatccagatctacgtacgtacgttgcatgcatgcaag");
        // No trailing newline, and surrounding blank lines, read the same.
        assert_eq!(
            parse_oligo_text(b"\r\n  ACGTACGTACGT  \r\n\r\n", Some("x"))
                .unwrap()
                .0[0]
                .residues,
            "ACGTACGTACGT"
        );
        // A named list is still read as a list, even when it holds one line.
        let (oligos, _) = parse_oligo_text(b"PROM1\tacgtggcaggtcgtatacg\n", Some("list")).unwrap();
        assert_eq!(oligos[0].name, "PROM1");
        // Without a name to give it there is nothing to call the oligo.
        assert!(parse_oligo_text(b"ACGTACGTACGT\n", None).is_err());
        // Several bare sequences are not joined into one invented oligo.
        assert!(parse_oligo_text(b"ACGTACGTACGT\nTGCATGCATGCA\n", Some("two")).is_err());
        // A bare sequence of the wrong length says so, with the file's name.
        let err = parse_oligo_text(b"ACGT\n", Some("short")).unwrap_err();
        assert!(err.contains("short: 4 bases"), "{err}");
    }

    #[test]
    fn nothing_readable_is_an_error_not_a_panic() {
        for bad in [
            &b""[..],
            b"\n\n",
            b">\n",
            b">x\n",
            b"\xFF\xFE\x00",
            b"x\ty",
            b"\xFE\xFF\xD8\x00",
        ] {
            assert!(parse_oligo_text(bad, None).is_err(), "{bad:?}");
        }
        assert!(parse_oligo_text(&vec![b'A'; MAX_FILE + 1], None).is_err());
        // A FASTA record with an endless sequence is cut short, not grown without limit.
        let long = format!(">x\n{}\n", "A".repeat(100_000));
        assert!(parse_oligo_text(long.as_bytes(), None).is_err());
    }

    #[test]
    fn adding_keeps_existing_oligos() {
        let o = |n: &str, s: &str| oligo(n, s, "").unwrap();
        let mut lib = OligoLibrary::default();
        let r = lib
            .add(vec![
                o("SP_301", "acgtacgtaccaggttacgt"),
                o("SP_302", "GACGTACCGGTTACCTAGCTA"),
                // The same oligo again, in capitals: a duplicate.
                o("SP_301", "ACGTACGTACCAGGTTACGT"),
            ])
            .unwrap();
        assert_eq!((r.added, r.duplicates), (2, 1));
        // The same name for another sequence: the library keeps its own.
        let r = lib.add(vec![o("sp_302", "ACGTACGTACGTACGT")]).unwrap();
        assert_eq!(r.conflicts, vec!["sp_302"]);
        assert_eq!(lib.oligos[1].residues, "GACGTACCGGTTACCTAGCTA");
        // Full: an error, and nothing changes.
        let mut full = OligoLibrary {
            oligos: (0..MAX_OLIGOS)
                .map(|i| o(&format!("o{i}"), "ACGTACGTACGT"))
                .collect(),
        };
        assert!(full.add(vec![o("new", "ACGTACGTACGT")]).is_err());
        assert_eq!(full.oligos.len(), MAX_OLIGOS);
        // A hand-edited file: bad and repeated entries are dropped.
        let edited = OligoLibrary {
            oligos: vec![
                o("a", "ACGTACGTACGT"),
                LabOligo {
                    name: "bad".into(),
                    residues: "ZZZ".into(),
                    description: String::new(),
                },
                o("A", "ACGTACGTACGT"),
            ],
        };
        assert_eq!(edited.checked().oligos.len(), 1);
    }
}
