//! EMBL ↔ document.
//!
//! EMBL's feature table follows the same INSDC rules as GenBank's, so features are read and
//! written exactly as for GenBank (see the `genbank` module): names from `/label`, `/gene`,
//! `/product` or `/locus_tag`, the same strands, and only edited entries re-encoded. The
//! document shows the first record; DE lines become the description.

use snappiesnap_doc::{Alphabet, DocError, FeatureId, Molecule, Notes, Seq, Topology};
use snappiesnap_formats::embl::{EmblFile, EmblRecord, IdLine};
use snappiesnap_formats::genbank::{Feature, FeatureEntry, Qualifier as GbQualifier};

use crate::genbank::{add_features, strip_tags, sync_features};
use crate::{Exported, Format, Imported, IoError, Source};

#[derive(Debug, Clone, PartialEq)]
pub struct EmblSource {
    original: Vec<u8>,
    file: EmblFile,
    imported: Molecule,
    /// Document feature id → index in the first record's features.
    feature_map: Vec<(FeatureId, usize)>,
}

pub(crate) fn looks_like_embl(bytes: &[u8]) -> bool {
    bytes.starts_with(b"ID   ")
}

pub(crate) fn import(name: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    let file = EmblFile::parse(bytes).map_err(|e| IoError::Format(e.to_string()))?;
    let record = file
        .records
        .first()
        .ok_or_else(|| IoError::Format("the file has no records".into()))?;
    let alphabet = if record.id.molecule.to_ascii_uppercase().contains("RNA") {
        Alphabet::Rna
    } else {
        Alphabet::Dna
    };
    let seq = Seq::new(alphabet, record.residues.to_ascii_uppercase()).map_err(DocError::from)?;
    let topology = if record.id.circular {
        Topology::Circular
    } else {
        Topology::Linear
    };
    let mut m = Molecule::new(name, seq, topology)?;
    let mut warnings = Vec::new();
    let description = description(&record.header);
    if !description.is_empty() {
        let notes = Notes {
            description,
            ..Notes::default()
        };
        if let Err(e) = m.set_notes(notes) {
            warnings.push(format!("description not shown ({e})"));
        }
    }
    let feature_map = add_features(&record.features, &mut m, topology, &mut warnings);
    if file.records.len() > 1 {
        warnings.push(format!(
            "{} more record(s) in the file are kept but not shown",
            file.records.len() - 1
        ));
    }
    Ok(Imported {
        source: Source::Embl(Box::new(EmblSource {
            original: bytes.to_vec(),
            file,
            imported: m.clone(),
            feature_map,
        })),
        molecule: m,
        format: Format::Embl,
        warnings,
    })
}

/// The DE lines' text, joined.
fn description(header: &str) -> String {
    header
        .lines()
        .filter_map(|l| l.strip_prefix("DE"))
        .map(str::trim)
        .filter(|t| !t.is_empty())
        .collect::<Vec<_>>()
        .join(" ")
}

/// An ID line name: the document name without spaces or semicolons.
fn id_name(name: &str) -> String {
    let n: String = name
        .chars()
        .map(|c| {
            if c.is_whitespace() || c == ';' {
                '_'
            } else {
                c
            }
        })
        .take(40)
        .collect();
    if n.is_empty() { "Untitled".into() } else { n }
}

pub(crate) fn export(m: &Molecule, source: Option<&EmblSource>) -> Result<Exported, IoError> {
    let mut warnings = Vec::new();
    if !m.primers().is_empty() {
        warnings.push("primers aren't saved in EMBL files".to_string());
    }
    if m.seq().alphabet() == Alphabet::Protein {
        return Err(IoError::Unsupported(
            "EMBL files hold nucleotide sequences only".into(),
        ));
    }
    if let Some(src) = source
        && *m == src.imported
    {
        return Ok(Exported {
            bytes: src.original.clone(),
            source: Source::Embl(Box::new(src.clone())),
            warnings,
        });
    }
    let (mut file, imported, feature_map) = match source {
        Some(src) => (
            src.file.clone(),
            Some(&src.imported),
            src.feature_map.clone(),
        ),
        None => (new_file(m), None, Vec::new()),
    };
    let record = file
        .records
        .first_mut()
        .ok_or_else(|| IoError::Format("the file has no records".into()))?;
    record.set_residues(m.seq().as_bytes().to_ascii_lowercase());
    let mut id = record.id.clone();
    id.length = m.len() as u64;
    id.circular = m.topology() == Topology::Circular;
    if imported.is_some_and(|i| i.name() != m.name()) {
        id.name = id_name(m.name());
    }
    record.set_id(id);
    let feature_map = sync_features(&mut record.features, m, imported, &feature_map);
    let bytes = file.to_text().into_bytes();
    let reread = EmblFile::parse(&bytes).map_err(|e| IoError::Format(e.to_string()))?;
    Ok(Exported {
        source: Source::Embl(Box::new(EmblSource {
            original: bytes.clone(),
            file: reread,
            imported: m.clone(),
            feature_map,
        })),
        bytes,
        warnings,
    })
}

/// A new file for a document that wasn't read from EMBL.
fn new_file(m: &Molecule) -> EmblFile {
    let description = strip_tags(&m.notes().description);
    let description = if description.trim().is_empty() {
        m.name().to_string()
    } else {
        description.split_whitespace().collect::<Vec<_>>().join(" ")
    };
    let header = format!(
        "XX\nAC   ;\nXX\nDE   {description}\nXX\nKW   .\nXX\nOS   synthetic construct\nOC   other sequences; artificial sequences.\nXX\n"
    );
    let id = IdLine {
        name: id_name(m.name()),
        circular: m.topology() == Topology::Circular,
        molecule: if m.seq().alphabet() == Alphabet::Rna {
            "other RNA".into()
        } else {
            "other DNA".into()
        },
        length: m.len() as u64,
    };
    let mut record = EmblRecord::new(id, header);
    let text = |name: &str, value: &str| GbQualifier {
        name: name.into(),
        value: Some(value.into()),
        quoted: true,
    };
    record.features.push(Feature {
        entry: FeatureEntry {
            key: "source".into(),
            location: format!("1..{}", m.len()),
            qualifiers: vec![
                text("organism", "synthetic construct"),
                text("mol_type", "other DNA"),
            ],
        },
        raw: None,
    });
    EmblFile {
        records: vec![record],
    }
}
