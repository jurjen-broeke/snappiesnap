//! SnappieSnap experiment files (`.sslab`): a gel experiment with its constructs, so it can be
//! closed, reopened and shared. The format is specified for other applications in
//! `docs/experiment-format.md`; keep the two in step.
//!
//! A file is one UTF-8 JSON object. Each construct is kept twice: as a reference to its file
//! (relative to the experiment file, and as originally opened) with a fingerprint of its sequence, and as an
//! embedded GenBank copy. On opening, the file is used when it is found and its sequence still
//! has the fingerprint; otherwise the copy, with a warning.
//!
//! Fields this version doesn't know are kept and written back, so a file edited by a newer
//! version or another application loses nothing when it is saved again here.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::gel::LaneSpec;

/// The value of `format`.
pub const FORMAT: &str = "snappiesnap-experiment";
/// The version written, and the newest one read.
pub const VERSION: u32 = 1;
/// The file extension, without the dot.
pub const EXTENSION: &str = "sslab";
/// Extensions of earlier versions, still recognised when opening (0.6.0 wrote `.ssexp`).
pub const OLD_EXTENSIONS: &[&str] = &["ssexp"];
/// Largest experiment file read.
pub const MAX_BYTES: usize = 256 << 20;
/// Most constructs and lanes.
pub const MAX_CONSTRUCTS: usize = 49;
pub const MAX_LANES: usize = 49;
/// Most reads in a reads experiment.
pub const MAX_READS: usize = 200;

/// A whole experiment file.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExperimentFile {
    pub format: String,
    pub version: u32,
    /// What wrote the file, such as "SnappieSnap 0.6.0" (informative).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub generator: Option<String>,
    /// The kind of experiment; version 1 knows "gel".
    pub kind: String,
    /// How the lanes were chosen, if by a design ("diagnostic_digest": digests chosen to tell
    /// the constructs apart). Informative; readers that don't know it show a gel experiment.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub design: Option<String>,
    pub name: String,
    #[serde(default)]
    pub gel: GelSettings,
    pub constructs: Vec<Construct>,
    pub lanes: Vec<Lane>,
    /// Kind "reads": the sequencing reads checked against the one construct.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reads: Vec<ReadFile>,
    /// Kind "reads": also check the reads stored with the construct in its .dna file.
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub matched_reads: bool,
    #[serde(flatten)]
    pub other: BTreeMap<String, Value>,
}

/// The gel the lanes run on.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GelSettings {
    /// Agarose, % (0.1 to 5).
    pub agarose: f64,
    /// The ladder in lane 1, by name ("1 kb", "100 bp", "High range").
    pub ladder: String,
    #[serde(flatten)]
    pub other: BTreeMap<String, Value>,
}

impl Default for GelSettings {
    fn default() -> Self {
        GelSettings {
            agarose: 1.0,
            ladder: "1 kb".to_string(),
            other: BTreeMap::new(),
        }
    }
}

/// One construct, used by one or more lanes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Construct {
    pub name: String,
    /// Its file relative to the experiment file's folder, with `/` between folders; absent
    /// when the construct had never been saved.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    /// The construct's path as the saving application knew it (absolute when the file was
    /// opened from a file dialog); tried when `file` isn't found.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub original_path: Option<String>,
    /// Residues.
    pub length: usize,
    /// "linear" or "circular".
    pub topology: String,
    /// `fnv1a64:` and 16 lowercase hex digits (see [`fingerprint`]).
    pub fingerprint: String,
    /// The construct as a GenBank flat file, as it was when the experiment was saved.
    pub genbank: String,
    #[serde(flatten)]
    pub other: BTreeMap<String, Value>,
}

/// A sequencing read, by its file (reads are raw data the user keeps; no copy is embedded).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReadFile {
    /// The read's name (its file name without the extension).
    pub name: String,
    /// The file relative to the experiment file's folder, with `/` between folders.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file: Option<String>,
    /// The path as the saving application knew it; tried when `file` isn't found.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub original_path: Option<String>,
    #[serde(flatten)]
    pub other: BTreeMap<String, Value>,
}

/// One lane after the ladder.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Lane {
    /// Index into `constructs`.
    pub construct: usize,
    /// What is done with the construct; null leaves the lane empty.
    #[serde(default)]
    pub method: Option<LaneSpec>,
    #[serde(flatten)]
    pub other: BTreeMap<String, Value>,
}

/// The fingerprint of a sequence: FNV-1a (64-bit; offset basis 0xcbf29ce484222325, prime
/// 0x100000001b3; Fowler, Noll and Vo, as published at www.isthe.com/chongo/tech/comp/fnv) over
/// the residues in upper case, written `fnv1a64:` plus 16 lowercase hex digits. It tells a
/// changed construct file from an unchanged one; it is not a security measure.
pub fn fingerprint(residues: &[u8]) -> String {
    let upper: Vec<u8> = residues.iter().map(u8::to_ascii_uppercase).collect();
    format!("fnv1a64:{:016x}", fnv1a64(&upper))
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

impl ExperimentFile {
    /// Read a file, checking what later steps rely on (indices, counts, numbers).
    pub fn parse(bytes: &[u8]) -> Result<ExperimentFile, String> {
        if bytes.len() > MAX_BYTES {
            return Err(format!("the file is over {} MB", MAX_BYTES >> 20));
        }
        let value: Value = serde_json::from_slice(bytes)
            .map_err(|e| format!("not a SnappieSnap experiment (invalid JSON: {e})"))?;
        if value.get("format").and_then(Value::as_str) != Some(FORMAT) {
            return Err(format!(
                "not a SnappieSnap experiment (\"format\" must be \"{FORMAT}\")"
            ));
        }
        let version = value.get("version").and_then(Value::as_u64).unwrap_or(0);
        if version == 0 || version > u64::from(VERSION) {
            return Err(format!(
                "experiment file version {version}; this SnappieSnap reads version {VERSION}"
            ));
        }
        let f: ExperimentFile =
            serde_json::from_value(value).map_err(|e| format!("damaged experiment file: {e}"))?;
        f.check()?;
        Ok(f)
    }

    /// The file's bytes: pretty-printed JSON with a final newline.
    pub fn to_bytes(&self) -> Result<Vec<u8>, String> {
        self.check()?;
        let mut out = serde_json::to_vec_pretty(self).map_err(|e| e.to_string())?;
        out.push(b'\n');
        Ok(out)
    }

    fn check(&self) -> Result<(), String> {
        match self.kind.as_str() {
            "gel" => {}
            "reads" => {
                if self.constructs.len() != 1 {
                    return Err("a reads experiment has exactly one construct".into());
                }
                if self.reads.len() > MAX_READS {
                    return Err(format!("at most {MAX_READS} reads"));
                }
            }
            other => {
                return Err(format!(
                    "experiment kind {other:?} isn't known to this SnappieSnap (it knows \"gel\" and \"reads\")"
                ));
            }
        }
        if !self.gel.agarose.is_finite() || !(0.1..=5.0).contains(&self.gel.agarose) {
            return Err("agarose must be 0.1-5 %".into());
        }
        if self.constructs.len() > MAX_CONSTRUCTS || self.lanes.len() > MAX_LANES {
            return Err(format!(
                "at most {MAX_CONSTRUCTS} constructs and {MAX_LANES} lanes"
            ));
        }
        if let Some(l) = self
            .lanes
            .iter()
            .find(|l| l.construct >= self.constructs.len())
        {
            return Err(format!(
                "a lane refers to construct {}, but there are {}",
                l.construct,
                self.constructs.len()
            ));
        }
        if let Some(c) = self
            .constructs
            .iter()
            .find(|c| !matches!(c.topology.as_str(), "linear" | "circular"))
        {
            return Err(format!(
                "construct {:?}: topology must be \"linear\" or \"circular\"",
                c.name
            ));
        }
        Ok(())
    }
}

/// The parts of a path, split at `/` or `\`, without empty parts and `.`.
fn parts(path: &str) -> Vec<&str> {
    path.split(['/', '\\'])
        .filter(|p| !p.is_empty() && *p != ".")
        .collect()
}

/// The folder a file is in (everything before its last separator; "" for a bare name).
pub fn folder_of(path: &str) -> &str {
    path.rfind(['/', '\\'])
        .map_or("", |i| path.get(..i).unwrap_or(""))
}

fn is_absolute(path: &str) -> bool {
    path.starts_with(['/', '\\']) || path.get(1..3) == Some(":\\") || path.get(1..3) == Some(":/")
}

/// `target` relative to the folder `base`, with `/` separators. `None` when there is no
/// relative path (one is absolute and the other isn't, or Windows drives differ).
pub fn relative_path(base: &str, target: &str) -> Option<String> {
    if is_absolute(base) != is_absolute(target) {
        return None;
    }
    let (b, t) = (parts(base), parts(target));
    // Windows paths on different drives have no relative path.
    let drive = |p: &[&str]| {
        p.first()
            .filter(|s| s.ends_with(':'))
            .map(|s| s.to_ascii_lowercase())
    };
    if drive(&b) != drive(&t) {
        return None;
    }
    let common = b.iter().zip(&t).take_while(|(x, y)| x == y).count();
    let mut out: Vec<&str> = vec![".."; b.len().saturating_sub(common)];
    out.extend(t.get(common..).unwrap_or(&[]));
    (!out.is_empty()).then(|| out.join("/"))
}

/// `relative` (with `/` or `\` separators) resolved against the folder `base`. `..` beyond the
/// top of `base` is an error rather than a path outside it.
pub fn resolve(base: &str, relative: &str) -> Option<String> {
    if is_absolute(relative) {
        return None;
    }
    let mut out: Vec<&str> = parts(base);
    for p in relative
        .split(['/', '\\'])
        .filter(|p| !p.is_empty() && *p != ".")
    {
        if p == ".." {
            out.pop()?;
        } else {
            out.push(p);
        }
    }
    let joined = out.join("/");
    Some(if base.starts_with(['/', '\\']) {
        format!("/{joined}")
    } else {
        joined
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fingerprint_by_hand() {
        // Published FNV-1a 64 test vectors: "" is the offset basis, "a" and "foobar".
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a64(b"foobar"), 0x8594_4171_f739_67e8);
        assert_eq!(
            fingerprint(b"acgt"),
            format!("fnv1a64:{:016x}", fnv1a64(b"ACGT"))
        );
        assert_ne!(fingerprint(b"ACGT"), fingerprint(b"ACGA"));
    }

    #[test]
    fn relative_paths() {
        assert_eq!(
            relative_path("/a/b", "/a/b/c.dna").as_deref(),
            Some("c.dna")
        );
        assert_eq!(
            relative_path("/a/b", "/a/x/c.dna").as_deref(),
            Some("../x/c.dna")
        );
        assert_eq!(
            relative_path("/a/b", "/c.dna").as_deref(),
            Some("../../c.dna")
        );
        assert_eq!(relative_path("/a", "c.dna"), None);
        assert_eq!(relative_path("C:\\x", "D:\\y\\c.dna"), None);
        assert_eq!(
            relative_path("C:\\x", "C:\\x\\y\\c.dna").as_deref(),
            Some("y/c.dna")
        );
        assert_eq!(relative_path("", "c.dna").as_deref(), Some("c.dna"));
        assert_eq!(resolve("/a/b", "../x/c.dna").as_deref(), Some("/a/x/c.dna"));
        assert_eq!(resolve("/a/b", "c.dna").as_deref(), Some("/a/b/c.dna"));
        assert_eq!(resolve("", "c.dna").as_deref(), Some("c.dna"));
        assert_eq!(resolve("/a", "../../c.dna"), None);
        assert_eq!(resolve("/a", "/etc/c.dna"), None);
        assert_eq!(folder_of("/a/b/e.sslab"), "/a/b");
        assert_eq!(folder_of("e.sslab"), "");
    }

    fn sample() -> ExperimentFile {
        ExperimentFile {
            format: FORMAT.into(),
            version: VERSION,
            generator: Some("test".into()),
            kind: "gel".into(),
            design: None,
            name: "Run".into(),
            gel: GelSettings::default(),
            constructs: vec![Construct {
                name: "p".into(),
                file: Some("p.dna".into()),
                original_path: Some("/x/p.dna".into()),
                length: 4,
                topology: "circular".into(),
                fingerprint: fingerprint(b"ACGT"),
                genbank: "LOCUS p\n//\n".into(),
                other: BTreeMap::new(),
            }],
            reads: Vec::new(),
            matched_reads: false,
            lanes: vec![Lane {
                construct: 0,
                method: Some(LaneSpec::Digest {
                    enzymes: vec!["EcoRI".into()],
                }),
                other: BTreeMap::new(),
            }],
            other: BTreeMap::new(),
        }
    }

    #[test]
    fn round_trip_keeps_unknown_fields() {
        let bytes = sample().to_bytes().unwrap();
        assert_eq!(ExperimentFile::parse(&bytes).unwrap(), sample());
        let mut v: Value = serde_json::from_slice(&bytes).unwrap();
        v["notes"] = Value::from("from a newer version");
        v["lanes"][0]["volume_ul"] = Value::from(10);
        let f = ExperimentFile::parse(&serde_json::to_vec(&v).unwrap()).unwrap();
        let again: Value = serde_json::from_slice(&f.to_bytes().unwrap()).unwrap();
        assert_eq!(again["notes"], "from a newer version");
        assert_eq!(again["lanes"][0]["volume_ul"], 10);
    }

    #[test]
    fn bad_files_are_errors() {
        let good: Value = serde_json::from_slice(&sample().to_bytes().unwrap()).unwrap();
        let bad = |f: &dyn Fn(&mut Value)| {
            let mut v = good.clone();
            f(&mut v);
            ExperimentFile::parse(&serde_json::to_vec(&v).unwrap()).unwrap_err()
        };
        assert!(bad(&|v| v["format"] = "x".into()).contains("not a SnappieSnap"));
        assert!(bad(&|v| v["version"] = 2.into()).contains("version 2"));
        assert!(bad(&|v| v["kind"] = "western".into()).contains("kind"));
        assert!(bad(&|v| v["lanes"][0]["construct"] = 5.into()).contains("construct 5"));
        assert!(bad(&|v| v["gel"]["agarose"] = 50.into()).contains("agarose"));
        assert!(bad(&|v| v["constructs"][0]["topology"] = "loop".into()).contains("topology"));
        assert!(
            bad(&|v| v["lanes"][0]["method"] = serde_json::json!({"kind": "boil"}))
                .contains("damaged")
        );
        assert!(ExperimentFile::parse(b"\xff{").is_err());
        assert!(ExperimentFile::parse(b"[]").is_err());
    }
}
