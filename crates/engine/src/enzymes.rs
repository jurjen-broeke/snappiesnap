//! The session's enzyme database: the built-in set, or a REBASE file the user downloaded or
//! loaded. REBASE is never bundled; a downloaded copy is cached in the platform's data folder
//! and loaded again at the next start ([`crate::Session::load_enzyme_cache`]).
//!
//! The database also carries the user's own enzyme sets ([`UserSet`]), such as the enzymes a
//! lab keeps in its freezer. They are kept in the data folder too
//! ([`crate::Session::load_enzyme_sets`]) and stay when the database is switched.

use snappiesnap_bio::enzymes::{Enzyme, builtin, parse_withrefm};

/// Where REBASE's `withrefm` file is published.
pub const REBASE_URL: &str = "http://rebase.neb.com/rebase/link_withrefm";
/// The cached download, inside [`crate::Platform::data_dir`].
pub const REBASE_CACHE: &str = "rebase_withrefm.txt";
/// The user's enzyme sets, inside [`crate::Platform::data_dir`].
pub const SETS_FILE: &str = "enzyme_sets.json";
/// Most enzyme sets kept.
pub const MAX_SETS: usize = 100;
/// Most enzymes in one set.
pub const MAX_SET_SIZE: usize = 2000;
/// Longest set name, in characters.
pub const MAX_SET_NAME: usize = 80;
/// Longest enzyme name accepted in a set (REBASE's longest is well under this).
const MAX_ENZYME_NAME: usize = 40;
/// Largest enzyme set file read.
pub const MAX_SET_FILE: usize = 1 << 20;
/// Credit shown wherever REBASE data is used.
pub const REBASE_CREDIT: &str = "REBASE, The Restriction Enzyme Database (rebase.neb.com), Roberts RJ, Vincze T, Posfai J, Macelis D.";

/// A set of enzymes and where it came from.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnzymeDb {
    /// `"built-in"` or `"REBASE"`.
    pub source: String,
    pub version: Option<String>,
    pub enzymes: Vec<Enzyme>,
    /// The user's own enzyme sets (names only; looked up in `enzymes`).
    pub sets: Vec<UserSet>,
    /// How many suppliers sell each enzyme, by name (from REBASE; empty for the built-in set).
    pub suppliers: std::collections::BTreeMap<String, usize>,
}

impl Default for EnzymeDb {
    fn default() -> Self {
        EnzymeDb {
            source: "built-in".to_string(),
            version: None,
            enzymes: builtin(),
            sets: Vec::new(),
            // The counts travel with the built-in table, so the name shown for a site that
            // several enzymes cut is the one most suppliers sell, before REBASE is downloaded.
            suppliers: snappiesnap_bio::enzymes::builtin_suppliers(),
        }
    }
}

impl EnzymeDb {
    /// Read a REBASE `withrefm` file. Returns the database and how many records were skipped.
    pub fn from_rebase(bytes: &[u8]) -> Result<(EnzymeDb, usize), String> {
        let r = parse_withrefm(bytes).map_err(|e| e.to_string())?;
        if r.enzymes.is_empty() {
            return Err("the file has no usable enzymes".to_string());
        }
        Ok((
            EnzymeDb {
                source: "REBASE".to_string(),
                version: r.version,
                enzymes: r.enzymes,
                sets: Vec::new(),
                suppliers: r.suppliers,
            },
            r.skipped,
        ))
    }

    /// Find an enzyme by name, ignoring case.
    pub fn get(&self, name: &str) -> Option<(usize, &Enzyme)> {
        let name = name.trim();
        self.enzymes
            .iter()
            .enumerate()
            .find(|(_, e)| e.name.eq_ignore_ascii_case(name))
    }

    pub fn is_builtin(&self) -> bool {
        self.source == "built-in"
    }

    /// The user's enzyme set called `name`, ignoring case.
    /// Order for naming enzymes that cut the same way: those in the user's own sets first (what
    /// the lab uses), then those sold by more suppliers, then the built-in list's.
    pub fn naming_rank(&self, name: &str) -> (bool, std::cmp::Reverse<usize>, bool) {
        let in_sets = self
            .sets
            .iter()
            .any(|s| s.enzymes.iter().any(|e| e.eq_ignore_ascii_case(name)));
        let suppliers = self.suppliers.get(name).copied().unwrap_or(0);
        let builtin = snappiesnap_bio::enzymes::builtin()
            .iter()
            .any(|e| e.name.eq_ignore_ascii_case(name));
        (!in_sets, std::cmp::Reverse(suppliers), !builtin)
    }

    pub fn set(&self, name: &str) -> Option<&UserSet> {
        let name = name.trim();
        self.sets.iter().find(|s| s.name.eq_ignore_ascii_case(name))
    }
}

/// An enzyme set of the user's own, such as the enzymes a lab has in stock. It holds names
/// only. They are looked up in whichever database is in use, so a name the built-in set lacks
/// is kept (REBASE may have it) and reported by [`UserSet::missing`].
#[derive(
    Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize, schemars::JsonSchema,
)]
#[serde(deny_unknown_fields)]
pub struct UserSet {
    pub name: String,
    pub enzymes: Vec<String>,
}

impl UserSet {
    /// A set from a name and enzyme names: names trimmed, duplicates (ignoring case) dropped.
    /// Names that can't be enzyme names are returned separately.
    pub fn new(name: &str, enzymes: &[String]) -> Result<(UserSet, Vec<String>), String> {
        let name = name.trim();
        if name.is_empty() {
            return Err("the enzyme set needs a name".to_string());
        }
        if name.chars().count() > MAX_SET_NAME {
            return Err(format!("set names are at most {MAX_SET_NAME} characters"));
        }
        let mut kept: Vec<String> = Vec::new();
        let mut rejected = Vec::new();
        for e in enzymes {
            let e = e.trim();
            if e.is_empty() {
                continue;
            }
            if !is_enzyme_name(e) {
                rejected.push(e.chars().take(MAX_ENZYME_NAME).collect());
            } else if !kept.iter().any(|k| k.eq_ignore_ascii_case(e)) {
                kept.push(e.to_string());
            }
        }
        if kept.is_empty() {
            return Err("the enzyme set names no enzymes".to_string());
        }
        if kept.len() > MAX_SET_SIZE {
            return Err(format!(
                "an enzyme set holds at most {MAX_SET_SIZE} enzymes"
            ));
        }
        Ok((
            UserSet {
                name: name.to_string(),
                enzymes: kept,
            },
            rejected,
        ))
    }

    /// The set's enzymes that `db` doesn't have, in set order.
    pub fn missing<'a>(&'a self, db: &EnzymeDb) -> Vec<&'a str> {
        self.enzymes
            .iter()
            .filter(|n| db.get(n).is_none())
            .map(String::as_str)
            .collect()
    }
}

/// Whether `s` can be an enzyme name: letters, digits and the few marks REBASE uses.
fn is_enzyme_name(s: &str) -> bool {
    s.chars().count() <= MAX_ENZYME_NAME
        && s.chars().any(|c| c.is_ascii_alphabetic())
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '.' | '_'))
}

/// An enzyme set exported as text: an optional title line followed by a blank line, then
/// enzyme names separated by commas, semicolons, tabs or line breaks (SnapGene's "Export
/// Enzyme Set" writes this layout). UTF-8 and UTF-16 (with a byte-order mark) are read.
/// Returns the title, if there is one, and the names in file order.
pub fn parse_set_text(bytes: &[u8]) -> Result<(Option<String>, Vec<String>), String> {
    if bytes.len() > MAX_SET_FILE {
        return Err(format!(
            "the file is larger than {} KB; an enzyme set file is a short list of names",
            MAX_SET_FILE >> 10
        ));
    }
    let text = crate::text::decode_text(bytes);
    let lines: Vec<&str> = text.lines().map(str::trim).collect();
    let first = lines.iter().position(|l| !l.is_empty());
    // A title is a first line without separators that a blank line sets apart.
    let title = first.and_then(|i| {
        let line = lines.get(i)?;
        let blank_next = lines.get(i.saturating_add(1)).is_some_and(|l| l.is_empty());
        let more = lines
            .iter()
            .skip(i.saturating_add(1))
            .any(|l| !l.is_empty());
        (blank_next && more && !line.contains([',', ';', '\t'])).then(|| (i, line.to_string()))
    });
    let body = lines
        .iter()
        .enumerate()
        .filter(|(i, _)| title.as_ref().is_none_or(|(t, _)| t != i))
        .map(|(_, l)| *l);
    let mut names = Vec::new();
    for line in body {
        for name in line.split([',', ';', '\t']).map(str::trim) {
            if !name.is_empty() {
                names.push(name.to_string());
            }
        }
    }
    if names.is_empty() {
        return Err("the file names no enzymes".to_string());
    }
    Ok((title.map(|(_, t)| t), names))
}

/// Which enzymes to report sites for.
#[derive(Debug, Clone, PartialEq, Eq, serde::Deserialize, schemars::JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct SiteQuery {
    /// Exactly these enzymes (case-insensitive), reported even if they don't cut. When set,
    /// the other filters are ignored.
    pub enzymes: Option<Vec<String>>,
    /// Only enzymes from this enzyme set of the user's own (see `enzyme_sets_list`), filtered
    /// by `min_length` and `max_cuts`. `commercial_only` doesn't apply: the set says what is
    /// at hand. Enzymes the database doesn't have are left out.
    pub set: Option<String>,
    /// Fewest specific (non-N) bases in the site.
    pub min_length: usize,
    /// Most sites an enzyme may have to be reported (1 = unique cutters).
    pub max_cuts: usize,
    /// Only enzymes that are sold commercially.
    pub commercial_only: bool,
}

impl Default for SiteQuery {
    /// Unique 6+ cutters, commercial enzymes only.
    fn default() -> Self {
        SiteQuery {
            enzymes: None,
            set: None,
            min_length: 6,
            max_cuts: 1,
            commercial_only: true,
        }
    }
}

/// One enzyme's sites.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize)]
pub struct EnzymeSites {
    pub name: String,
    pub notation: String,
    pub sites: Vec<snappiesnap_bio::enzymes::Site>,
}

/// Most enzymes one query may name.
pub const MAX_NAMED: usize = 1000;

/// The sites of the enzymes `q` selects, in database order. Unknown names are an error that
/// lists them.
pub fn sites(
    db: &EnzymeDb,
    m: &snappiesnap_doc::Molecule,
    q: &SiteQuery,
) -> Result<Vec<EnzymeSites>, String> {
    use snappiesnap_bio::enzymes::{find_sites, sites_by_enzyme};
    let chosen: Vec<Enzyme> = match &q.enzymes {
        Some(names) => {
            if names.len() > MAX_NAMED {
                return Err(format!("at most {MAX_NAMED} enzymes per query"));
            }
            let mut out = Vec::new();
            let mut unknown = Vec::new();
            for n in names {
                match db.get(n) {
                    Some((_, e)) if !out.contains(e) => out.push(e.clone()),
                    Some(_) => {}
                    None => unknown.push(n.as_str()),
                }
            }
            if !unknown.is_empty() {
                return Err(format!("unknown enzymes: {}", unknown.join(", ")));
            }
            out
        }
        None => {
            let members: Option<std::collections::HashSet<String>> = match &q.set {
                Some(name) => Some(
                    db.set(name)
                        .ok_or_else(|| format!("there is no enzyme set called {name:?}"))?
                        .enzymes
                        .iter()
                        .map(|n| n.to_ascii_lowercase())
                        .collect(),
                ),
                None => None,
            };
            db.enzymes
                .iter()
                .filter(|e| e.cuts() && e.specific_length() >= q.min_length)
                .filter(|e| match &members {
                    Some(m) => m.contains(&e.name.to_ascii_lowercase()),
                    None => e.commercial || !q.commercial_only,
                })
                .cloned()
                .collect()
        }
    };
    let mut grouped = sites_by_enzyme(find_sites(m, &chosen))
        .into_iter()
        .peekable();
    let mut out = Vec::new();
    for (i, e) in chosen.iter().enumerate() {
        let sites = match grouped.peek() {
            Some((g, _)) if *g == i => grouped.next().map(|(_, s)| s).unwrap_or_default(),
            _ => Vec::new(),
        };
        let keep = q.enzymes.is_some() || (!sites.is_empty() && sites.len() <= q.max_cuts);
        if keep {
            out.push(EnzymeSites {
                name: e.name.clone(),
                notation: e.notation(),
                sites,
            });
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use snappiesnap_doc::{Alphabet, Molecule, Seq, Topology};

    use super::*;

    #[test]
    fn queries() {
        let s = format!(
            "{}GAATTC{}GGATCC{}GGATCC",
            "A".repeat(10),
            "A".repeat(20),
            "A".repeat(5)
        );
        let m =
            Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Circular).unwrap();
        let db = EnzymeDb::default();
        let names = |q: &SiteQuery| -> Vec<String> {
            sites(&db, &m, q)
                .unwrap()
                .into_iter()
                .map(|e| e.name)
                .collect()
        };
        assert!(names(&SiteQuery::default()).contains(&"EcoRI".to_string()));
        assert!(!names(&SiteQuery::default()).contains(&"BamHI".to_string()));
        let two = SiteQuery {
            max_cuts: 2,
            ..SiteQuery::default()
        };
        assert!(names(&two).contains(&"BamHI".to_string()));
        let named = SiteQuery {
            enzymes: Some(vec!["ecori".into(), "NotI".into(), "EcoRI".into()]),
            ..SiteQuery::default()
        };
        let r = sites(&db, &m, &named).unwrap();
        assert_eq!(r.len(), 2);
        assert_eq!(
            (r[0].name.as_str(), r[0].sites.len(), r[1].sites.len()),
            ("EcoRI", 1, 0)
        );
        let bad = SiteQuery {
            enzymes: Some(vec!["Nope".into()]),
            ..SiteQuery::default()
        };
        assert!(sites(&db, &m, &bad).unwrap_err().contains("Nope"));

        // A user set narrows the pool; the filters still apply within it.
        let mut db = db;
        db.sets.push(
            UserSet::new("lab", &["bamhi".into(), "NotI".into(), "Missing1".into()])
                .unwrap()
                .0,
        );
        let from_set = |max_cuts| SiteQuery {
            set: Some("LAB".into()),
            max_cuts,
            ..SiteQuery::default()
        };
        let found = |q: &SiteQuery| -> Vec<String> {
            sites(&db, &m, q)
                .unwrap()
                .into_iter()
                .map(|e| e.name)
                .collect()
        };
        assert!(found(&from_set(1)).is_empty());
        assert_eq!(found(&from_set(2)), vec!["BamHI"]);
        let no_set = SiteQuery {
            set: Some("nope".into()),
            ..SiteQuery::default()
        };
        assert!(sites(&db, &m, &no_set).unwrap_err().contains("nope"));
        assert_eq!(db.set("lab").unwrap().missing(&db), vec!["Missing1"]);
    }

    #[test]
    fn user_sets_keep_valid_unique_names() {
        let names: Vec<String> = ["EcoRI", " ecori ", "", "Hpy188III", "Bad name", "x<y>"]
            .map(String::from)
            .to_vec();
        let (set, rejected) = UserSet::new(" Lab ", &names).unwrap();
        assert_eq!(set.name, "Lab");
        assert_eq!(set.enzymes, vec!["EcoRI", "Hpy188III"]);
        assert_eq!(rejected, vec!["Bad name", "x<y>"]);
        assert!(UserSet::new("", &names).is_err());
        assert!(UserSet::new("x", &[]).is_err());
        assert!(UserSet::new(&"x".repeat(MAX_SET_NAME + 1), &names).is_err());
        let many: Vec<String> = (0..=MAX_SET_SIZE).map(|i| format!("E{i}")).collect();
        assert!(UserSet::new("x", &many).is_err());
    }

    #[test]
    fn set_files_parse() {
        // The layout of an exported set: title, blank line, wrapped comma-separated names.
        let text = b"LAB ENZYMES\r\n\r\nPvuII, XmnI,\r\nEcoRI, BamHI\r\n";
        let (title, names) = parse_set_text(text).unwrap();
        assert_eq!(title.as_deref(), Some("LAB ENZYMES"));
        assert_eq!(names, vec!["PvuII", "XmnI", "EcoRI", "BamHI"]);
        // One name per line, no title.
        let (title, names) = parse_set_text(b"EcoRI\nBamHI\n").unwrap();
        assert_eq!((title, names.len()), (None, 2));
        // A lone first line followed by a blank line and nothing else is a name, not a title.
        let (title, names) = parse_set_text(b"EcoRI\n\n").unwrap();
        assert_eq!((title, names), (None, vec!["EcoRI".to_string()]));
        // UTF-16 LE and BE with a byte-order mark, and a UTF-8 BOM.
        let utf16 = |le: bool| {
            let mut out = if le {
                vec![0xFF, 0xFE]
            } else {
                vec![0xFE, 0xFF]
            };
            for u in "Set\r\n\r\nEcoRI;SmaI".encode_utf16() {
                out.extend(if le { u.to_le_bytes() } else { u.to_be_bytes() });
            }
            out
        };
        for le in [true, false] {
            let (title, names) = parse_set_text(&utf16(le)).unwrap();
            assert_eq!((title.as_deref(), names.len()), (Some("Set"), 2));
        }
        assert_eq!(
            parse_set_text(b"\xEF\xBB\xBFEcoRI").unwrap().1,
            vec!["EcoRI"]
        );
        // Hostile input: errors, never panics.
        for bad in [&b""[..], b"\n\n , ;", b"\xFF\xFE\x00", b"\xFE\xFF\xD8\x00"] {
            let _ = parse_set_text(bad);
        }
        assert!(parse_set_text(b"").is_err());
        assert!(parse_set_text(&vec![b'a'; MAX_SET_FILE + 1]).is_err());
    }
}
