//! Sequencing traces (ABIF `.ab1`, SCF) ↔ document. The read becomes a linear DNA molecule
//! with its chromatogram attached. We don't write trace files: an unchanged read is saved
//! from its original bytes, and an edited one goes to a sequence format.

use snappiesnap_doc::{Alphabet, Call, DocError, Molecule, Seq, Signal, Topology, Trace};
use snappiesnap_formats::trace::Chromatogram;
use snappiesnap_formats::{abif, scf};

use crate::{Exported, Format, Imported, IoError, Source};

#[derive(Debug, Clone, PartialEq)]
pub struct TraceSource {
    original: Vec<u8>,
    format: Format,
    imported: Molecule,
}

impl TraceSource {
    pub fn format(&self) -> Format {
        self.format
    }
}

pub(crate) fn looks_like_abif(bytes: &[u8]) -> bool {
    bytes.starts_with(b"ABIF")
}

pub(crate) fn looks_like_scf(bytes: &[u8]) -> bool {
    bytes.starts_with(b".scf")
}

pub(crate) fn import(name: &str, bytes: &[u8], format: Format) -> Result<Imported, IoError> {
    let c = match format {
        Format::Scf => scf::read_chromatogram(bytes).map_err(|e| IoError::Format(e.to_string()))?,
        _ => abif::read_chromatogram(bytes).map_err(|e| IoError::Format(e.to_string()))?,
    };
    let (molecule, warnings) = molecule(name, c)?;
    Ok(Imported {
        source: Source::Trace(Box::new(TraceSource {
            original: bytes.to_vec(),
            format,
            imported: molecule.clone(),
        })),
        molecule,
        format,
        warnings,
    })
}

/// The read as a molecule. Letters that aren't IUPAC bases (some basecallers write `-`) become
/// `N`, and the warnings say how many.
fn molecule(name: &str, c: Chromatogram) -> Result<(Molecule, Vec<String>), IoError> {
    let mut warnings = Vec::new();
    let mut odd = 0usize;
    let bases: Vec<u8> = c
        .bases
        .iter()
        .map(|&b| {
            if Alphabet::Dna.is_valid(b) {
                b
            } else {
                odd += 1;
                b'N'
            }
        })
        .collect();
    if odd > 0 {
        warnings.push(format!(
            "{odd} base call(s) weren't IUPAC letters and are shown as N"
        ));
    }
    let seq = Seq::new(Alphabet::Dna, bases).map_err(DocError::from)?;
    let [a, cc, g, t] = c.channels;
    let quality = c.quality.unwrap_or_default();
    if quality.is_empty() {
        warnings.push("the file has no quality values".to_string());
    }
    let calls = c
        .peaks
        .iter()
        .enumerate()
        .map(|(i, &p)| Call {
            peak: Some(p),
            quality: quality.get(i).copied(),
        })
        .collect();
    let (trace, dropped) = Trace::new(Signal { a, c: cc, g, t }, calls, c.sample_name)?;
    if dropped > 0 {
        warnings.push(format!(
            "{dropped} base call(s) have a peak outside the signal and are shown without one"
        ));
    }
    let mut m = Molecule::new(name, seq, Topology::Linear)?;
    m.set_trace(Some(trace))?;
    Ok((m, warnings))
}

pub(crate) fn export(
    m: &Molecule,
    source: Option<&TraceSource>,
    format: Format,
) -> Result<Exported, IoError> {
    match source {
        Some(src) if src.format == format && *m == src.imported => Ok(Exported {
            bytes: src.original.clone(),
            source: Source::Trace(Box::new(src.clone())),
            warnings: Vec::new(),
        }),
        Some(src) if src.format == format => Err(IoError::Unsupported(format!(
            "writing {} files; save the edited read as .dna, GenBank or FASTA",
            format.name()
        ))),
        _ => Err(IoError::Unsupported(format!(
            "writing {} files; only trace files opened unchanged are saved as such",
            format.name()
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn chromatogram() -> Chromatogram {
        Chromatogram {
            sample_name: Some("read 1".into()),
            bases: b"ACG-T".to_vec(),
            peaks: vec![2, 6, 10, 14, 99],
            quality: Some(vec![20, 30, 40, 5, 9]),
            channels: [vec![1; 20], vec![2; 20], vec![3; 20], vec![4; 20]],
        }
    }

    #[test]
    fn reads_become_documents_and_save_back_unchanged() {
        let bytes = abif::write_chromatogram(&chromatogram()).unwrap();
        let imported = crate::import("/x/read1.ab1", &bytes).unwrap();
        assert_eq!(imported.format, Format::Abif);
        let m = &imported.molecule;
        assert_eq!(m.seq().as_bytes(), b"ACGNT");
        assert_eq!(m.name(), "read1");
        let t = m.trace().unwrap();
        assert_eq!(t.sample_name(), Some("read 1"));
        assert_eq!(t.calls()[2].quality, Some(40));
        assert_eq!(t.calls()[4].peak, None);
        assert_eq!(imported.warnings.len(), 2, "{:?}", imported.warnings);
        let out = crate::export(m, Some(&imported.source), Format::Abif).unwrap();
        assert_eq!(out.bytes, bytes);
    }

    #[test]
    fn edited_reads_go_to_sequence_formats() {
        let bytes = scf::write_chromatogram(&chromatogram()).unwrap();
        let imported = crate::import("r.scf", &bytes).unwrap();
        assert_eq!(imported.format, Format::Scf);
        let mut m = imported.molecule.clone();
        m.delete(snappiesnap_doc::Range::new(0, 1)).unwrap();
        let err = crate::export(&m, Some(&imported.source), Format::Scf).unwrap_err();
        assert!(err.to_string().contains("save the edited read"), "{err}");
        assert!(crate::export(&m, None, Format::Abif).is_err());
        let fasta = crate::export(&m, Some(&imported.source), Format::Fasta).unwrap();
        assert!(fasta.bytes.ends_with(b"\nCGNT\n"));
        assert!(
            fasta.warnings.iter().any(|w| w.contains("chromatogram")),
            "{:?}",
            fasta.warnings
        );
    }

    #[test]
    fn damaged_traces_are_errors() {
        let mut bytes = abif::write_chromatogram(&chromatogram()).unwrap();
        bytes.truncate(140);
        assert!(crate::import("x.ab1", &bytes).is_err());
        assert!(crate::import("x.scf", b".scf").is_err());
    }
}

/// A sequencing read matched to a construct in a `.dna` file.
#[derive(Debug, Clone, PartialEq)]
pub struct MatchedRead {
    pub name: String,
    /// The stretch kept after trimming in the reference app: 0-based, end exclusive.
    pub trimmed: Option<(usize, usize)>,
    /// The read as a linear molecule (with its trace, if it is one), or why it can't be read.
    pub molecule: Result<Molecule, String>,
}

impl crate::DnaSource {
    /// The reads matched to this construct and stored in its file, in list order.
    pub fn matched_reads(&self) -> Vec<MatchedRead> {
        let reads = match self.file.matched_reads() {
            Ok(r) => r,
            Err(e) => {
                return vec![MatchedRead {
                    name: "matched reads".to_string(),
                    trimmed: None,
                    molecule: Err(e.to_string()),
                }];
            }
        };
        reads
            .into_iter()
            .map(|r| {
                let name = {
                    let n: String = r
                        .name
                        .chars()
                        .filter(|c| !c.is_control())
                        .take(200)
                        .collect();
                    if n.trim().is_empty() {
                        format!("read {}", r.id)
                    } else {
                        n
                    }
                };
                let molecule = r.read.and_then(|data| match data.trace {
                    Some(c) => molecule(&name, c)
                        .map(|(m, _)| m)
                        .map_err(|e| e.to_string()),
                    None => {
                        let bases: Vec<u8> = data
                            .bases
                            .iter()
                            .map(|&b| if Alphabet::Dna.is_valid(b) { b } else { b'N' })
                            .collect();
                        Seq::new(Alphabet::Dna, bases)
                            .map_err(|e| e.to_string())
                            .and_then(|s| {
                                Molecule::new(&name, s, Topology::Linear).map_err(|e| e.to_string())
                            })
                    }
                });
                MatchedRead {
                    name,
                    trimmed: r.trimmed,
                    molecule,
                }
            })
            .collect()
    }
}
