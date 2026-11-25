//! Gel lanes: what each lane holds, and its fragments. Shared by the `gel` command and the Gel
//! view, so agents and users see the same gel. A gel is either of the open molecule, or an
//! experiment ([`experiment`]): each lane a construct (an open document) and a method.

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use snappiesnap_bio::enzymes::digest;
use snappiesnap_bio::gel::{Band, Fragment, ladder, lane};
use snappiesnap_bio::pcr::{DEFAULT_MAX_PRODUCT, pcr};
use snappiesnap_doc::{Molecule, Topology};

use crate::Session;
use crate::enzymes::{EnzymeDb, SiteQuery, sites};

/// What a lane holds.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum LaneSpec {
    /// A size ladder by name ("1 kb", "100 bp", "High range").
    Ladder { name: String },
    /// The molecule itself (supercoiled if circular).
    Uncut,
    /// The molecule cut with these enzymes.
    Digest { enzymes: Vec<String> },
    /// The products of a PCR with two primers (a primer's name, or a sequence).
    Pcr { forward: String, reverse: String },
    /// Fragments of given sizes (bp).
    Sizes { sizes: Vec<usize> },
}

impl LaneSpec {
    /// A short title for the lane.
    pub fn title(&self) -> String {
        match self {
            LaneSpec::Ladder { name } => name.clone(),
            LaneSpec::Uncut => "Uncut".to_string(),
            LaneSpec::Digest { enzymes } => enzymes.join(" + "),
            LaneSpec::Pcr { forward, reverse } => format!("PCR {forward} / {reverse}"),
            LaneSpec::Sizes { .. } => "Fragments".to_string(),
        }
    }
}

/// Most fragments in a hand-made lane.
const MAX_SIZES: usize = 1000;

/// A primer by name in the molecule, or else the text as a sequence.
pub fn primer_residues(m: &Molecule, text: &str) -> Result<(String, Vec<u8>), String> {
    let t = text.trim();
    if let Some(p) = m.primers().iter().find(|p| p.name.eq_ignore_ascii_case(t)) {
        return Ok((p.name.clone(), p.seq.as_bytes().to_vec()));
    }
    let bases: Vec<u8> = t.bytes().filter(|b| !b.is_ascii_whitespace()).collect();
    if bases.len() >= 10
        && bases
            .iter()
            .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(b))
    {
        Ok((format!("{} nt", bases.len()), bases))
    } else {
        Err(format!(
            "{t:?} is neither a primer of this molecule nor a primer sequence"
        ))
    }
}

/// The fragments a lane holds.
pub fn fragments(m: &Molecule, db: &EnzymeDb, spec: &LaneSpec) -> Result<Vec<Fragment>, String> {
    let linear = |bp: usize, label: String| Fragment {
        bp,
        copies: 1.0,
        supercoiled: false,
        label,
    };
    Ok(match spec {
        LaneSpec::Ladder { .. } => Vec::new(),
        LaneSpec::Uncut => vec![Fragment {
            bp: m.len(),
            copies: 1.0,
            supercoiled: m.topology() == Topology::Circular,
            label: if m.topology() == Topology::Circular {
                "uncut, supercoiled (apparent size approximate)".to_string()
            } else {
                "uncut".to_string()
            },
        }],
        LaneSpec::Digest { enzymes } => {
            if enzymes.is_empty() {
                return Err("name at least one enzyme".to_string());
            }
            let q = SiteQuery {
                enzymes: Some(enzymes.clone()),
                ..SiteQuery::default()
            };
            let chosen = sites(db, m, &q)?;
            let list: Vec<_> = chosen
                .iter()
                .filter_map(|c| db.get(&c.name).map(|(_, e)| e.clone()))
                .collect();
            digest(m, &list)
                .into_iter()
                .map(|f| {
                    let ends = match (&f.left, &f.right) {
                        (Some(a), Some(b)) => format!("{a} \u{2013} {b}"),
                        (Some(a), None) | (None, Some(a)) => format!("{a} \u{2013} end"),
                        (None, None) => "uncut".to_string(),
                    };
                    // An uncut circle in a digest lane is still supercoiled.
                    let uncut_circle = f.left.is_none() && m.topology() == Topology::Circular;
                    Fragment {
                        bp: f.length,
                        copies: 1.0,
                        supercoiled: uncut_circle,
                        label: ends,
                    }
                })
                .collect()
        }
        LaneSpec::Pcr { forward, reverse } => {
            let (fname, f) = primer_residues(m, forward)?;
            let (rname, r) = primer_residues(m, reverse)?;
            pcr(m, [&f, &r], None, DEFAULT_MAX_PRODUCT)
                .into_iter()
                .map(|p| linear(p.length, format!("PCR product {fname} / {rname}")))
                .collect()
        }
        LaneSpec::Sizes { sizes } => {
            if sizes.len() > MAX_SIZES {
                return Err(format!("at most {MAX_SIZES} fragments per lane"));
            }
            sizes
                .iter()
                .map(|&bp| linear(bp, format!("{bp} bp")))
                .collect()
        }
    })
}

/// The bands of a lane at `agarose` %.
pub fn bands(
    m: &Molecule,
    db: &EnzymeDb,
    spec: &LaneSpec,
    agarose: f64,
) -> Result<Vec<Band>, String> {
    if let LaneSpec::Ladder { name } = spec {
        return ladder_bands(name, agarose);
    }
    Ok(lane(&fragments(m, db, spec)?, agarose))
}

/// The bands of a size ladder by name.
pub fn ladder_bands(name: &str, agarose: f64) -> Result<Vec<Band>, String> {
    ladder(name, agarose).ok_or_else(|| {
        let names: Vec<&str> = snappiesnap_bio::gel::LADDERS.iter().map(|l| l.0).collect();
        format!("unknown ladder {name:?} (known: {})", names.join(", "))
    })
}

/// One lane of a gel experiment: a construct (an open document, by its stable id; see the
/// `documents` command) and what is done with it. Without a method the lane stays empty.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ExperimentLane {
    pub document: u64,
    #[serde(default)]
    pub method: Option<LaneSpec>,
}

/// A lane as run: its title, and its bands or why it can't be run.
pub type RunLane = (String, Result<Vec<Band>, String>);

/// The lanes of an experiment gel: the ladder first, then one lane per construct in order.
/// Each lane runs exactly as on a gel of that construct alone ([`bands`]).
pub fn experiment(
    session: &Session,
    ladder: &str,
    agarose: f64,
    lanes: &[ExperimentLane],
) -> Vec<RunLane> {
    let mut out = vec![(ladder.to_string(), ladder_bands(ladder, agarose))];
    for l in lanes {
        let Some(doc) = session.document_by_id(l.document) else {
            out.push((
                "closed construct".to_string(),
                Err(format!(
                    "the construct of this lane (document {}) is no longer open",
                    l.document
                )),
            ));
            continue;
        };
        let m = &doc.molecule;
        out.push(match &l.method {
            Some(spec) => (
                format!("{}: {}", m.name(), spec.title()),
                bands(m, session.enzymes(), spec, agarose),
            ),
            None => (format!("{}: no method yet", m.name()), Ok(Vec::new())),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, Seq};

    use super::*;

    #[test]
    fn lanes() {
        let s = format!(
            "{}GAATTC{}GGATCC{}",
            "A".repeat(100),
            "C".repeat(1000),
            "T".repeat(2000)
        );
        let mut m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, s.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        let db = EnzymeDb::default();
        let digest = LaneSpec::Digest {
            enzymes: vec!["EcoRI".into(), "BamHI".into()],
        };
        let f = fragments(&m, &db, &digest).unwrap();
        let mut sizes: Vec<usize> = f.iter().map(|x| x.bp).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![1006, s.len() - 1006]);
        let uncut = fragments(&m, &db, &LaneSpec::Uncut).unwrap();
        assert!(uncut[0].supercoiled);
        m.add_primer(
            "F",
            Seq::new(Alphabet::Dna, "AAAAAAAAAAAAAAAAAAAAGAATTC").unwrap(),
            "",
        )
        .unwrap();
        assert!(primer_residues(&m, "f").is_ok());
        assert!(primer_residues(&m, "nonsense").is_err());
        assert!(
            bands(
                &m,
                &db,
                &LaneSpec::Ladder {
                    name: "1 kb".into()
                },
                1.0
            )
            .unwrap()
            .len()
                > 5
        );
        assert!(bands(&m, &db, &LaneSpec::Ladder { name: "x".into() }, 1.0).is_err());
        assert!(fragments(&m, &db, &LaneSpec::Digest { enzymes: vec![] }).is_err());
        let json = serde_json::to_value(&digest).unwrap();
        assert_eq!(json["kind"], "digest");
    }
}
