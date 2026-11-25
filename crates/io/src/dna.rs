//! SnapGene `.dna` ↔ document.
//!
//! Coordinates: the file stores segments 1-based and inclusive (`range="1-3"` is the first
//! three bases; `start > end` wraps past the origin); the document uses 0-based half-open
//! ranges. Segments are stored in ascending position; the document lists them 5'→3' along the
//! feature, so the order is reversed for reverse-strand features.

use snappiesnap_doc::{
    Alphabet, DocError, FeatureId, Molecule, NewFeature, Notes, PrimerId, Qualifier, Range,
    Segment, Seq, Strand, Topology,
};
use snappiesnap_formats::dna::{
    DnaFile, FeatureXml, FeaturesXml, MAGIC, Packet, SequencePacket, tag,
};
use snappiesnap_formats::xml::{Token, Tokenizer};

use snappiesnap_doc::Feature;

use crate::{Format, Imported, IoError, Source};

/// What a `.dna` document was read from, for faithful saving.
#[derive(Debug, Clone, PartialEq)]
pub struct DnaSource {
    /// The file exactly as read.
    pub(crate) file: DnaFile,
    /// The parsed features packet (empty if the file had none).
    pub(crate) features: FeaturesXml,
    pub(crate) had_features_packet: bool,
    /// The molecule as imported, to tell what changed since.
    pub(crate) imported: Molecule,
    /// Document feature id → index in `features`.
    pub(crate) feature_map: Vec<(FeatureId, usize)>,
    /// Document primer id → index of its `<Primer>` in the primers packet.
    pub(crate) primer_map: Vec<(PrimerId, usize)>,
}

impl DnaSource {
    pub fn file(&self) -> &DnaFile {
        &self.file
    }
}

pub(crate) fn looks_like_dna(bytes: &[u8]) -> bool {
    bytes.first() == Some(&tag::COOKIE) && bytes.get(5..13) == Some(MAGIC.as_slice())
}

pub(crate) fn import(name: &str, bytes: &[u8]) -> Result<Imported, IoError> {
    import_with(name, bytes, true)
}

/// A `.dna` file's molecule without its construction history (for the steps of a history).
pub(crate) fn import_plain(name: &str, bytes: &[u8]) -> Result<Molecule, IoError> {
    Ok(import_with(name, bytes, false)?.molecule)
}

fn import_with(name: &str, bytes: &[u8], with_history: bool) -> Result<Imported, IoError> {
    let file = DnaFile::parse(bytes)?;
    let mut warnings = Vec::new();

    let cookie = file.cookie()?;
    if cookie.sequence_type != 1 {
        return Err(IoError::Unsupported(format!(
            ".dna files of sequence type {} (only DNA, type 1, so far)",
            cookie.sequence_type
        )));
    }
    let sequence = file.sequence()?;
    let topology = if sequence.is_circular() {
        Topology::Circular
    } else {
        Topology::Linear
    };
    let seq = Seq::new(Alphabet::Dna, sequence.residues).map_err(DocError::from)?;
    let mut molecule = Molecule::new(name, seq, topology)?;

    if file.packet(tag::NOTES).is_some() {
        match file.notes() {
            Ok(fields) => {
                let mut notes = Notes::default();
                for f in fields {
                    if f.name == "Description" {
                        notes.description = f.text;
                    } else {
                        notes.fields.insert(f.name, f.text);
                    }
                }
                if let Err(e) = molecule.set_notes(notes) {
                    warnings.push(format!("notes not shown ({e}); they are kept in the file"));
                }
            }
            Err(e) => warnings.push(format!("notes not shown ({e}); they are kept in the file")),
        }
    }

    let had_features_packet = file.packet(tag::FEATURES).is_some();
    let features = if had_features_packet {
        file.features()?
    } else {
        FeaturesXml::empty()
    };
    let mut feature_map = Vec::new();
    for (index, xml) in features.features().enumerate() {
        match to_new_feature(xml, molecule.len())
            .and_then(|f| molecule.add_feature(f).map_err(|e| e.to_string()))
        {
            Ok(id) => feature_map.push((id, index)),
            Err(reason) => warnings.push(format!(
                "feature {:?} not shown ({reason}); it is kept in the file",
                xml.name()
            )),
        }
    }

    let mut primer_map = Vec::new();
    for (index, primer) in primers_in(&file, &mut warnings).into_iter().enumerate() {
        let added = Seq::new(Alphabet::Dna, primer.sequence.as_bytes())
            .map_err(|e| e.to_string())
            .and_then(|seq| {
                molecule
                    .add_primer(&primer.name, seq, &primer.description)
                    .map_err(|e| e.to_string())
            });
        match added {
            Ok(id) => primer_map.push((id, index)),
            Err(reason) => warnings.push(format!(
                "primer {:?} not shown ({reason}); it is kept in the file",
                primer.name
            )),
        }
    }

    if with_history {
        match crate::history::read(&file).and_then(|h| Ok(molecule.set_history(h)?)) {
            Ok(()) => {}
            Err(e) => warnings.push(format!(
                "construction history not shown ({e}); it is kept in the file"
            )),
        }
    }

    let source = DnaSource {
        file,
        features,
        had_features_packet,
        imported: molecule.clone(),
        feature_map,
        primer_map,
    };
    Ok(Imported {
        molecule,
        source: Source::Dna(Box::new(source)),
        format: Format::Dna,
        warnings,
    })
}

/// Write `molecule` as a `.dna` file. With the source it was read from, unchanged content is
/// written back byte for byte and only edited parts are re-encoded; without one, a new file is
/// built. Returns the bytes and the source describing them (for the next save).
pub(crate) fn export(
    molecule: &Molecule,
    source: Option<&DnaSource>,
) -> Result<(Vec<u8>, DnaSource, Vec<String>), IoError> {
    let Some(src) = source else {
        return export_new(molecule);
    };
    if *molecule == src.imported {
        let bytes = src.file.to_bytes()?;
        return Ok((bytes, src.clone(), Vec::new()));
    }
    let mut file = src.file.clone();
    // The packets as read (the copy above is being edited).
    let text = |t: u8| {
        src.file
            .packet(t)
            .and_then(|p| std::str::from_utf8(&p.payload).ok())
            .map(str::to_string)
    };
    // Primers: when they changed, or when the sequence moved under them (binding sites).
    let moved =
        molecule.seq() != src.imported.seq() || molecule.topology() != src.imported.topology();
    let mut primer_map = src.primer_map.clone();
    if molecule.primers() != src.imported.primers() || (moved && !src.primer_map.is_empty()) {
        let (packet, owners) = crate::dna_notes::primers_packet(
            molecule,
            &src.imported,
            text(tag::PRIMERS).as_deref(),
            &src.primer_map,
        )?;
        file.set_payload(tag::PRIMERS, packet.into_bytes())?;
        primer_map = owners
            .into_iter()
            .enumerate()
            .filter_map(|(i, o)| o.map(|id| (id, i)))
            .collect();
    }
    if let Some(packet) = crate::dna_notes::notes_packet(
        molecule.notes(),
        src.imported.notes(),
        text(tag::NOTES).as_deref(),
    )? {
        file.set_payload(tag::NOTES, packet.into_bytes())?;
    }

    if molecule.seq() != src.imported.seq() || molecule.topology() != src.imported.topology() {
        let mut packet = file.sequence()?;
        packet.residues = molecule.seq().as_bytes().to_vec();
        packet.set_circular(molecule.topology() == Topology::Circular);
        file.set_payload(tag::SEQUENCE, packet.to_payload())?;
    }

    let (features, owners) = updated_features(molecule, src)?;
    if src.had_features_packet || !features.is_empty() {
        file.set_features(&features)?;
    }
    let mut warnings = Vec::new();
    if molecule.history() != src.imported.history() {
        crate::history::write(&mut file, molecule, &mut warnings)?;
    }
    let bytes = file.to_bytes()?;
    let mut source = source_after_save(molecule, file, owners)?;
    source.primer_map = primer_map;
    Ok((bytes, source, warnings))
}

/// The features packet for `molecule`, plus which document feature owns each `<Feature>`
/// (`None` for features kept in the file but not shown in the document).
fn updated_features(
    molecule: &Molecule,
    src: &DnaSource,
) -> Result<(FeaturesXml, Vec<Option<FeatureId>>), IoError> {
    let mut xml = src.features.clone();
    let len = molecule.len();
    let mut owners: Vec<Option<FeatureId>> = vec![None; xml.len()];
    for (id, index) in &src.feature_map {
        if let Some(slot) = owners.get_mut(*index) {
            *slot = Some(*id);
        }
    }

    // Edited features: change only what the document models.
    for (id, index) in &src.feature_map {
        let (Some(now), Some(then)) = (molecule.feature(*id), src.imported.feature(*id)) else {
            continue;
        };
        if now == then {
            continue;
        }
        let f = xml
            .feature_mut(*index)
            .ok_or_else(|| IoError::Unsupported("feature table changed unexpectedly".into()))?;
        apply_edit(f, now, then, len);
    }

    // Removed features (shown in the document, then deleted); hidden ones stay.
    let mut index = 0usize;
    xml.retain(|_, _| {
        let keep = match owners.get(index).copied().flatten() {
            Some(id) => molecule.feature(id).is_some(),
            None => true,
        };
        index = index.saturating_add(1);
        keep
    });
    owners.retain(|o| o.is_none_or(|id| molecule.feature(id).is_some()));

    // New features, appended with fresh ids.
    let mapped: Vec<FeatureId> = src.feature_map.iter().map(|(id, _)| *id).collect();
    let mut next = xml.next_valid_id().unwrap_or(0).max(
        xml.features()
            .filter_map(|f| f.attr("recentID").and_then(|v| v.parse::<u64>().ok()))
            .max()
            .map_or(0, |m| m.saturating_add(1)),
    );
    let mut added = false;
    for f in molecule.features() {
        if mapped.contains(&f.id) {
            continue;
        }
        xml.push(FeatureXml::new(
            next,
            f.color.as_deref(),
            &f.name,
            &f.kind,
            directionality(f.strand),
            &file_segments(f, len),
            &f.qualifiers
                .iter()
                .map(|q| (q.name.clone(), q.values.clone()))
                .collect::<Vec<_>>(),
        ));
        owners.push(Some(f.id));
        next = next.saturating_add(1);
        added = true;
    }
    if added {
        xml.set_next_valid_id(next);
    }
    Ok((xml, owners))
}

fn apply_edit(xml: &mut FeatureXml, now: &Feature, then: &Feature, len: usize) {
    if now.name != then.name {
        xml.set_attr("name", &now.name);
    }
    if now.kind != then.kind {
        xml.set_attr("type", &now.kind);
    }
    if now.strand != then.strand {
        match directionality(now.strand) {
            Some(d) => xml.set_attr("directionality", &d.to_string()),
            None => xml.remove_attr("directionality"),
        }
    }
    if now.segments != then.segments || now.strand != then.strand {
        xml.set_segment_ranges(&file_segments(now, len));
    }
    if now.qualifiers != then.qualifiers {
        xml.set_qualifiers(
            &now.qualifiers
                .iter()
                .map(|q| (q.name.clone(), q.values.clone()))
                .collect::<Vec<_>>(),
        );
    }
}

fn directionality(strand: Strand) -> Option<u8> {
    match strand {
        Strand::Forward => Some(1),
        Strand::Reverse => Some(2),
        Strand::Both => Some(3),
        Strand::None => None,
    }
}

/// Segments as stored: 1-based inclusive, in ascending order along the sequence (the document
/// lists reverse-strand segments 5'→3', i.e. descending).
fn file_segments(f: &Feature, len: usize) -> Vec<(u64, u64)> {
    let mut segs: Vec<(u64, u64)> = f
        .segments
        .iter()
        .map(|s| range_to_segment(s.range, len))
        .collect();
    if f.strand == Strand::Reverse {
        segs.reverse();
    }
    segs
}

/// The source after a save: the written file, with features re-read so their original text is
/// the written text, mapped to the same document ids.
fn source_after_save(
    molecule: &Molecule,
    file: DnaFile,
    owners: Vec<Option<FeatureId>>,
) -> Result<DnaSource, IoError> {
    let had_features_packet = file.packet(tag::FEATURES).is_some();
    let features = if had_features_packet {
        file.features()?
    } else {
        FeaturesXml::empty()
    };
    let feature_map = owners
        .into_iter()
        .enumerate()
        .filter_map(|(i, o)| o.map(|id| (id, i)))
        .collect();
    Ok(DnaSource {
        file,
        features,
        had_features_packet,
        imported: molecule.clone(),
        feature_map,
        primer_map: Vec::new(),
    })
}

/// Cookie values seen in files written by the reference application (sequence type 1 = DNA).
const NEW_FILE_COOKIE: [u8; 6] = [0, 1, 0, 15, 0, 19];

/// A new `.dna` file for a document that wasn't read from one. Experimental: the minimal
/// packet set (cookie, sequence, features, primers) hasn't been confirmed to open in the
/// reference application yet.
fn export_new(molecule: &Molecule) -> Result<(Vec<u8>, DnaSource, Vec<String>), IoError> {
    if molecule.seq().alphabet() != Alphabet::Dna {
        return Err(IoError::Unsupported("saving RNA or protein as .dna".into()));
    }
    let mut cookie = MAGIC.to_vec();
    cookie.extend_from_slice(&NEW_FILE_COOKIE);
    // 0x02 as on the linear files we've seen, plus the circular bit.
    let mut sequence = SequencePacket {
        flags: 0x02,
        residues: molecule.seq().as_bytes().to_vec(),
    };
    sequence.set_circular(molecule.topology() == Topology::Circular);
    let file = DnaFile::from_packets(vec![
        Packet {
            tag: tag::COOKIE,
            payload: cookie,
        },
        Packet {
            tag: tag::SEQUENCE,
            payload: sequence.to_payload(),
        },
        Packet {
            tag: tag::FEATURES,
            payload: FeaturesXml::empty().to_text().into_bytes(),
        },
        Packet {
            tag: tag::PRIMERS,
            payload: crate::dna_notes::EMPTY_PRIMERS.as_bytes().to_vec(),
        },
    ])?;
    let empty = DnaSource {
        features: file.features()?,
        had_features_packet: true,
        // Nothing yet: the features, primers and notes are all written as new.
        imported: Molecule::new(molecule.name(), molecule.seq().clone(), molecule.topology())?,
        feature_map: Vec::new(),
        primer_map: Vec::new(),
        file,
    };
    export(molecule, Some(&empty))
}

fn to_new_feature(xml: &FeatureXml, len: usize) -> Result<NewFeature, String> {
    let strand = match xml.directionality() {
        Some("1") => Strand::Forward,
        Some("2") => Strand::Reverse,
        Some("3") => Strand::Both,
        _ => Strand::None,
    };
    let infos = xml.segments()?;
    // One colour per feature: the first segment's (they nearly always agree).
    let color = infos
        .iter()
        .find_map(|s| s.color.clone())
        .filter(|c| snappiesnap_doc::is_color(c))
        .map(|c| c.to_ascii_lowercase());
    let mut segments = infos
        .into_iter()
        .map(|s| segment_to_range(s.start, s.end, len).map(|range| Segment { range }))
        .collect::<Result<Vec<_>, _>>()?;
    if strand == Strand::Reverse {
        segments.reverse();
    }
    Ok(NewFeature {
        color,
        name: xml.name().to_string(),
        kind: xml.kind().unwrap_or("misc_feature").to_string(),
        strand,
        segments,
        qualifiers: xml
            .qualifiers()
            .into_iter()
            .map(|(name, values)| Qualifier { name, values })
            .collect(),
    })
}

/// 1-based inclusive (`start > end` wraps) → 0-based half-open.
pub(crate) fn segment_to_range(start: u64, end: u64, len: usize) -> Result<Range, String> {
    let to_usize = |v: u64| usize::try_from(v).map_err(|_| format!("position {v} is too large"));
    let (start, end) = (to_usize(start)?, to_usize(end)?);
    if start == 0 || end == 0 || start > len || end > len {
        return Err(format!(
            "segment {start}-{end} is outside the sequence (1-{len})"
        ));
    }
    Ok(Range::new(start.saturating_sub(1), end))
}

/// 0-based half-open → 1-based inclusive (`start > end` for a range across the origin).
pub(crate) fn range_to_segment(range: Range, len: usize) -> (u64, u64) {
    let start = range.start.saturating_add(1);
    let end = if range.end == 0 { len } else { range.end };
    (start as u64, end as u64)
}

/// A primer as stored in the primers packet (binding sites are recomputed, not read).
struct PrimerXml {
    name: String,
    sequence: String,
    description: String,
}

/// The primers in the primers packet, in file order. Problems become warnings; the packet
/// itself is always written back as read.
fn primers_in(file: &DnaFile, warnings: &mut Vec<String>) -> Vec<PrimerXml> {
    let mut out = Vec::new();
    let Some(packet) = file.packet(tag::PRIMERS) else {
        return out;
    };
    let Ok(text) = std::str::from_utf8(&packet.payload) else {
        warnings
            .push("primers not shown (the packet isn't UTF-8); they are kept in the file".into());
        return out;
    };
    let mut tokens = Tokenizer::new(text);
    while let Some(next) = tokens.next_token() {
        let token = match next {
            Ok((token, _)) => token,
            Err(e) => {
                warnings.push(format!(
                    "primers after byte {} not shown ({}); they are kept in the file",
                    e.offset, e.reason
                ));
                break;
            }
        };
        if let Token::Start {
            name: "Primer",
            attrs,
            ..
        } = token
        {
            let get = |n: &str| {
                attrs
                    .iter()
                    .find(|a| a.name == n)
                    .map(|a| snappiesnap_formats::xml::decode(a.raw))
            };
            out.push(PrimerXml {
                name: get("name").unwrap_or_default(),
                sequence: get("sequence").unwrap_or_default(),
                description: get("description").unwrap_or_default(),
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coordinates_convert_both_ways() {
        assert_eq!(segment_to_range(1, 3, 10), Ok(Range::new(0, 3)));
        assert_eq!(segment_to_range(8, 2, 10), Ok(Range::new(7, 2)));
        assert!(segment_to_range(0, 3, 10).is_err());
        assert!(segment_to_range(1, 11, 10).is_err());
        assert!(segment_to_range(u64::MAX, 3, 10).is_err());
        assert_eq!(range_to_segment(Range::new(0, 3), 10), (1, 3));
        assert_eq!(range_to_segment(Range::new(7, 2), 10), (8, 2));
        for (a, b) in [(1, 3), (8, 2), (1, 10), (10, 10), (10, 1)] {
            let r = segment_to_range(a, b, 10).unwrap();
            assert_eq!(range_to_segment(r, 10), (a, b));
        }
    }
}
