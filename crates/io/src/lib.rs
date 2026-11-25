//! SnappieSnap `io` (layer L3).
//!
//! Maps file formats onto the document model. Import keeps the parsed file as the document's
//! [`Source`], so a later save can write unchanged content back byte for byte and re-encode only
//! what was edited (AGENTS.md › File fidelity). Anything that can't be represented in the
//! document stays in the source and is reported as a warning; nothing is dropped silently.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

mod dna;
mod dna_notes;
mod embl;
mod fasta;
mod genbank;
mod gff3;
mod history;
mod sbol;
mod trace;

use snappiesnap_doc::{DocError, Molecule};
use snappiesnap_formats::dna::DnaError;

pub use dna::DnaSource;
pub use embl::EmblSource;
pub use fasta::FastaSource;
pub use genbank::GenBankSource;
pub use gff3::{Gff3Source, annotations as gff3_annotations};
pub use sbol::SbolSource;
pub use trace::{MatchedRead, TraceSource};

/// File formats we can read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// SnapGene `.dna`.
    Dna,
    /// FASTA (`.fa`, `.fasta`, ...).
    Fasta,
    /// GenBank flat file (`.gb`, `.gbk`, ...).
    GenBank,
    /// EMBL flat file (`.embl`).
    Embl,
    /// GFF3 annotations with their sequence (`.gff3`).
    Gff3,
    /// SBOL 2 in RDF/XML (`.sbol`, `.xml`).
    Sbol,
    /// Applied Biosystems sequencing trace (`.ab1`); read only.
    Abif,
    /// Standard Chromatogram Format sequencing trace (`.scf`); read only.
    Scf,
}

impl Format {
    /// Detect the format from the content, falling back to the file extension.
    pub fn detect(path: &str, bytes: &[u8]) -> Option<Format> {
        if trace::looks_like_abif(bytes) {
            return Some(Format::Abif);
        }
        if trace::looks_like_scf(bytes) {
            return Some(Format::Scf);
        }
        if dna::looks_like_dna(bytes) {
            return Some(Format::Dna);
        }
        if fasta::looks_like_fasta(bytes) {
            return Some(Format::Fasta);
        }
        if genbank::looks_like_genbank(bytes) {
            return Some(Format::GenBank);
        }
        if embl::looks_like_embl(bytes) {
            return Some(Format::Embl);
        }
        if gff3::looks_like_gff3(bytes) {
            return Some(Format::Gff3);
        }
        if sbol::looks_like_sbol(bytes) {
            return Some(Format::Sbol);
        }
        format_for_path(path)
    }

    pub fn name(self) -> &'static str {
        match self {
            Format::Dna => "SnapGene .dna",
            Format::Fasta => "FASTA",
            Format::GenBank => "GenBank",
            Format::Embl => "EMBL",
            Format::Gff3 => "GFF3",
            Format::Sbol => "SBOL 2",
            Format::Abif => "ABIF trace (.ab1)",
            Format::Scf => "SCF trace",
        }
    }
}

/// The original file behind a document, kept for faithful saving.
#[derive(Debug, Clone, PartialEq)]
pub enum Source {
    Dna(Box<DnaSource>),
    Fasta(Box<FastaSource>),
    GenBank(Box<GenBankSource>),
    Embl(Box<EmblSource>),
    Gff3(Box<Gff3Source>),
    Sbol(Box<SbolSource>),
    Trace(Box<TraceSource>),
}

impl Source {
    pub fn format(&self) -> Format {
        match self {
            Source::Dna(_) => Format::Dna,
            Source::Fasta(_) => Format::Fasta,
            Source::GenBank(_) => Format::GenBank,
            Source::Embl(_) => Format::Embl,
            Source::Gff3(_) => Format::Gff3,
            Source::Sbol(_) => Format::Sbol,
            Source::Trace(t) => t.format(),
        }
    }
}

/// The result of an export.
#[derive(Debug, Clone, PartialEq)]
pub struct Exported {
    pub bytes: Vec<u8>,
    /// What the written file now is, for the next save.
    pub source: Source,
    /// What the format couldn't hold (for example features in FASTA).
    pub warnings: Vec<String>,
}

/// The result of an import.
#[derive(Debug, Clone, PartialEq)]
pub struct Imported {
    pub molecule: Molecule,
    pub source: Source,
    pub format: Format,
    /// What the document doesn't show (it is kept in the file and saved back).
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum IoError {
    #[error(
        "unrecognised file format (supported: SnapGene .dna, GenBank, EMBL, GFF3, SBOL 2, FASTA, \
         .ab1 and .scf traces)"
    )]
    UnknownFormat,

    #[error("{0}")]
    Format(String),

    #[error(transparent)]
    Dna(#[from] DnaError),

    #[error(transparent)]
    Doc(#[from] DocError),

    #[error("not supported yet: {0}")]
    Unsupported(String),
}

/// Read a file. `path` is used for the document name (its file stem) and as a format hint.
pub fn import(path: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    let format = Format::detect(path, bytes).ok_or(IoError::UnknownFormat)?;
    match format {
        Format::Dna => dna::import(&file_stem(path), bytes),
        Format::Fasta => fasta::import(&file_stem(path), bytes),
        Format::GenBank => genbank::import(&file_stem(path), bytes),
        Format::Embl => embl::import(&file_stem(path), bytes),
        Format::Gff3 => gff3::import(&file_stem(path), bytes),
        Format::Sbol => sbol::import(&file_stem(path), bytes),
        Format::Abif | Format::Scf => trace::import(&file_stem(path), bytes, format),
    }
}

/// Settings for writing new content.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ExportOptions {
    /// Today's date for GenBank LOCUS lines, like `22-OCT-2026` (see [`genbank_date`]).
    pub date: Option<String>,
}

/// A date in the GenBank LOCUS style (`22-OCT-2026`). `None` for an invalid month.
pub fn genbank_date(year: i32, month: u32, day: u32) -> Option<String> {
    const MONTHS: [&str; 12] = [
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
    ];
    let m = MONTHS.get(usize::try_from(month).ok()?.checked_sub(1)?)?;
    Some(format!("{day:02}-{m}-{year:04}"))
}

/// Write `molecule` in `format`. `source` is what the document was read from (if anything):
/// when it is the same format, unchanged content is written back byte for byte.
pub fn export(
    molecule: &Molecule,
    source: Option<&Source>,
    format: Format,
) -> Result<Exported, IoError> {
    export_with(molecule, source, format, &ExportOptions::default())
}

/// [`export`] with options for new content.
pub fn export_with(
    molecule: &Molecule,
    source: Option<&Source>,
    format: Format,
    options: &ExportOptions,
) -> Result<Exported, IoError> {
    let mut out = export_any(molecule, source, format, options)?;
    if molecule.trace().is_some() && !matches!(format, Format::Abif | Format::Scf) {
        out.warnings.push(format!(
            "{} files don't hold the chromatogram: only the sequence of the read was saved",
            format.name()
        ));
    }
    Ok(out)
}

fn export_any(
    molecule: &Molecule,
    source: Option<&Source>,
    format: Format,
    options: &ExportOptions,
) -> Result<Exported, IoError> {
    match format {
        Format::Abif | Format::Scf => {
            let src = match source {
                Some(Source::Trace(s)) => Some(s.as_ref()),
                _ => None,
            };
            trace::export(molecule, src, format)
        }
        Format::Dna => {
            let src = match source {
                Some(Source::Dna(s)) => Some(s.as_ref()),
                _ => None,
            };
            let mut warnings = Vec::new();
            let (bytes, new, more) = dna::export(molecule, src)?;
            warnings.extend(more);
            Ok(Exported {
                bytes,
                source: Source::Dna(Box::new(new)),
                warnings,
            })
        }
        Format::Fasta | Format::GenBank | Format::Embl | Format::Gff3 | Format::Sbol
            if molecule.history().is_some() =>
        {
            let mut out = export_any(&without_history(molecule), source, format, options)?;
            out.warnings.push(format!(
                "{} files don't hold the construction history",
                format.name()
            ));
            Ok(out)
        }
        Format::Fasta => {
            let src = match source {
                Some(Source::Fasta(s)) => Some(s.as_ref()),
                _ => None,
            };
            fasta::export(molecule, src)
        }
        Format::GenBank => {
            let src = match source {
                Some(Source::GenBank(s)) => Some(s.as_ref()),
                _ => None,
            };
            genbank::export(molecule, src, options)
        }
        Format::Embl => {
            let src = match source {
                Some(Source::Embl(s)) => Some(s.as_ref()),
                _ => None,
            };
            embl::export(molecule, src)
        }
        Format::Gff3 => {
            let src = match source {
                Some(Source::Gff3(s)) => Some(s.as_ref()),
                _ => None,
            };
            gff3::export(molecule, src)
        }
        Format::Sbol => {
            let src = match source {
                Some(Source::Sbol(s)) => Some(s.as_ref()),
                _ => None,
            };
            sbol::export(molecule, src)
        }
    }
}

fn without_history(m: &Molecule) -> Molecule {
    let mut m = m.clone();
    m.take_history();
    m
}

/// The format to save `path` as, from its extension.
pub fn format_for_path(path: &str) -> Option<Format> {
    let ext = path.rsplit_once('.').map(|(_, e)| e.to_ascii_lowercase());
    match ext.as_deref() {
        Some("dna") => Some(Format::Dna),
        Some("fa" | "fasta" | "fna" | "fas" | "ffn" | "seq") => Some(Format::Fasta),
        Some("gb" | "gbk" | "genbank" | "gbff" | "ape") => Some(Format::GenBank),
        Some("embl" | "emb") => Some(Format::Embl),
        Some("gff3" | "gff") => Some(Format::Gff3),
        // SBOL is the only XML format written, so .xml saves as SBOL.
        Some("sbol" | "rdf" | "xml") => Some(Format::Sbol),
        Some("ab1" | "abi" | "abif" | "fsa") => Some(Format::Abif),
        Some("scf") => Some(Format::Scf),
        _ => None,
    }
}

/// The file name without directories and extension, or "Untitled".
pub fn file_stem(path: &str) -> String {
    let name = path.rsplit(['/', '\\']).next().unwrap_or(path);
    let stem = match name.rsplit_once('.') {
        Some((stem, _)) if !stem.is_empty() => stem,
        _ => name,
    };
    let stem: String = stem.chars().filter(|c| !c.is_control()).collect();
    if stem.trim().is_empty() {
        "Untitled".to_string()
    } else {
        stem
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stems() {
        assert_eq!(file_stem("/a/b/pUC19.dna"), "pUC19");
        assert_eq!(file_stem("C:\\x\\my.plasmid.dna"), "my.plasmid");
        assert_eq!(file_stem(".dna"), ".dna");
        assert_eq!(file_stem(""), "Untitled");
        assert_eq!(file_stem("/tmp/"), "Untitled");
    }
}
