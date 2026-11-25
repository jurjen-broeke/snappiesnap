//! FASTA ↔ document. The document is the first record; other records are kept and written
//! back unchanged. FASTA has no topology, features, primers or notes.

use snappiesnap_doc::{Alphabet, DocError, Molecule, Seq, Topology};
use snappiesnap_formats::fasta::{DEFAULT_WIDTH, FastaFile, FastaRecord};

use crate::{Exported, Format, Imported, IoError, Source};

#[derive(Debug, Clone, PartialEq)]
pub struct FastaSource {
    original: Vec<u8>,
    file: FastaFile,
    imported: Molecule,
}

pub(crate) fn looks_like_fasta(bytes: &[u8]) -> bool {
    bytes
        .iter()
        .find(|b| !b.is_ascii_whitespace() && **b != 0xEF && **b != 0xBB && **b != 0xBF)
        == Some(&b'>')
}

/// The alphabet a residue string most likely is: DNA, else RNA, else protein.
pub(crate) fn guess_alphabet(residues: &[u8]) -> Alphabet {
    [Alphabet::Dna, Alphabet::Rna]
        .into_iter()
        .find(|a| residues.iter().all(|&b| a.is_valid(b)))
        .unwrap_or(Alphabet::Protein)
}

pub(crate) fn import(name: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    let file = FastaFile::parse(bytes).map_err(|e| IoError::Format(e.to_string()))?;
    let first = file
        .records
        .first()
        .ok_or_else(|| IoError::Format("the file has no records".into()))?;
    let seq = Seq::new(guess_alphabet(&first.residues), first.residues.clone())
        .map_err(DocError::from)?;
    let molecule = Molecule::new(name, seq, Topology::Linear)?;
    let mut warnings = Vec::new();
    if file.records.len() > 1 {
        warnings.push(format!(
            "{} more record(s) in the file are kept but not shown",
            file.records.len() - 1
        ));
    }
    Ok(Imported {
        source: Source::Fasta(Box::new(FastaSource {
            original: bytes.to_vec(),
            file,
            imported: molecule.clone(),
        })),
        molecule,
        format: Format::Fasta,
        warnings,
    })
}

pub(crate) fn export(m: &Molecule, source: Option<&FastaSource>) -> Result<Exported, IoError> {
    let mut warnings = Vec::new();
    if !m.features().is_empty() || !m.primers().is_empty() {
        warnings.push(
            "FASTA stores only the sequence: features and primers were not saved in this file"
                .to_string(),
        );
    }
    if let Some(src) = source
        && *m == src.imported
    {
        return Ok(Exported {
            bytes: src.original.clone(),
            source: Source::Fasta(Box::new(src.clone())),
            warnings,
        });
    }
    let mut file = match source {
        Some(src) => src.file.clone(),
        None => FastaFile {
            records: Vec::new(),
            line_width: DEFAULT_WIDTH,
        },
    };
    let old_header = file.records.first().map(|r| r.header.clone());
    let header = header_for(m, source.map(|s| &s.imported), old_header.as_deref());
    let record = FastaRecord {
        header,
        residues: m.seq().as_bytes().to_vec(),
    };
    match file.records.first_mut() {
        Some(first) => *first = record,
        None => file.records.push(record),
    }
    let bytes = file.to_text().into_bytes();
    let source = FastaSource {
        original: bytes.clone(),
        file,
        imported: m.clone(),
    };
    Ok(Exported {
        bytes,
        source: Source::Fasta(Box::new(source)),
        warnings,
    })
}

/// The header for the first record. The reference application writes `name  (N bp)`; in that
/// form the length is kept up to date. If the document was renamed, the header's first word
/// becomes the new name; the rest of the header is kept.
fn header_for(m: &Molecule, imported: Option<&Molecule>, old: Option<&str>) -> String {
    let Some(old) = old else {
        return format!("{}  ({} bp)", m.name(), m.len());
    };
    let renamed = imported.is_some_and(|prev| prev.name() != m.name());
    if let Some(header_name) = generated_name(old) {
        let name = if renamed { m.name() } else { header_name };
        return format!("{name}  ({} bp)", m.len());
    }
    if renamed {
        let rest = old.split_once(char::is_whitespace).map_or("", |(_, r)| r);
        return if rest.is_empty() {
            m.name().to_string()
        } else {
            format!("{} {rest}", m.name())
        };
    }
    old.to_string()
}

/// The name in a header of the form `name  (N bp)`.
fn generated_name(header: &str) -> Option<&str> {
    let (name, rest) = header.split_once("  (")?;
    let digits = rest.strip_suffix(" bp)")?;
    (!name.is_empty()
        && !name.contains(char::is_whitespace)
        && !digits.is_empty()
        && digits.bytes().all(|b| b.is_ascii_digit()))
    .then_some(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open(text: &str) -> Imported {
        import("x", text.as_bytes()).unwrap()
    }

    #[test]
    fn unchanged_is_byte_exact_and_edits_keep_other_records() {
        let text = ">pX  (12 bp)\nACGTACGTAC\nGT\n>other\nMKV\n";
        let i = open(text);
        assert_eq!(i.molecule.seq().alphabet(), Alphabet::Dna);
        assert_eq!(i.warnings.len(), 1);
        let Source::Fasta(src) = &i.source else {
            panic!()
        };
        assert_eq!(
            export(&i.molecule, Some(src)).unwrap().bytes,
            text.as_bytes()
        );

        let mut m = i.molecule.clone();
        m.insert(0, b"TT").unwrap();
        let out = export(&m, Some(src)).unwrap();
        assert_eq!(
            String::from_utf8(out.bytes).unwrap(),
            ">pX  (14 bp)\nTTACGTACGT\nACGT\n>other\nMKV\n"
        );
    }

    #[test]
    fn alphabets_and_headers() {
        assert_eq!(guess_alphabet(b"ACGTN"), Alphabet::Dna);
        assert_eq!(guess_alphabet(b"ACGU"), Alphabet::Rna);
        assert_eq!(guess_alphabet(b"MKVLE*"), Alphabet::Protein);
        let m = Molecule::new(
            "new",
            Seq::new(Alphabet::Dna, "ACGT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert_eq!(header_for(&m, None, None), "new  (4 bp)");
        let prev = Molecule::new(
            "old",
            Seq::new(Alphabet::Dna, "ACGT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        assert_eq!(
            header_for(&m, Some(&prev), Some("old some description")),
            "new some description"
        );
        assert_eq!(
            header_for(&m, Some(&prev), Some("old  (4 bp)")),
            "new  (4 bp)"
        );
        assert_eq!(
            header_for(&prev, Some(&prev), Some("pX  (9 bp)")),
            "pX  (4 bp)"
        );
        assert_eq!(generated_name("a  (12 bp)"), Some("a"));
        assert_eq!(generated_name("a b  (12 bp)"), None);
        assert_eq!(generated_name("a  (x bp)"), None);
    }

    #[test]
    fn features_are_warned_about_not_dropped_silently() {
        let mut m = Molecule::new(
            "p",
            Seq::new(Alphabet::Dna, "ACGTACGT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        m.add_feature(snappiesnap_doc::NewFeature {
            color: None,
            name: "f".into(),
            kind: "misc_feature".into(),
            strand: Default::default(),
            segments: vec![snappiesnap_doc::Segment {
                range: snappiesnap_doc::Range::new(0, 2),
            }],
            qualifiers: vec![],
        })
        .unwrap();
        let out = export(&m, None).unwrap();
        assert_eq!(out.warnings.len(), 1);
        assert_eq!(out.bytes, b">p  (8 bp)\nACGTACGT\n");
    }
}
