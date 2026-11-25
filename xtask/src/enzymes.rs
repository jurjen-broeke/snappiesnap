//! `cargo xtask enzymes`: rebuild the built-in enzyme table from a REBASE `withrefm` file.
//!
//! The rule, applied mechanically and nothing else: **every restriction enzyme with a known cut
//! that at least one supplier sells**, according to the supplier codes in the REBASE file given.
//! Only the name
//! and the recognition sequence with its cut positions are written — published biochemical
//! facts in the standard notation. No REBASE text, references, methylation or other fields are
//! copied, and the REBASE file itself is never committed (ATTRIBUTION.md).
//!
//! The file is read from the path given, else from `$SNAPPIESNAP_REBASE`, else from the copy the
//! app caches when a user downloads REBASE. The table is parsed back with our own parser before
//! it is written, so a table that wouldn't load can't be generated.
//!
//! ```sh
//! cargo xtask enzymes                      # the app's downloaded copy
//! cargo xtask enzymes /path/withrefm.txt
//! ```

use std::fmt::Write as _;
use std::path::PathBuf;

use snappiesnap_bio::enzymes::{Enzyme, parse_withrefm};

use crate::util::{self, Result};

/// Where the generated table lives.
const TABLE_FILE: &str = "crates/bio/src/enzymes/builtin.rs";
/// The file the app caches a downloaded REBASE in (see `engine::enzymes::REBASE_CACHE`).
const CACHE_NAME: &str = "rebase_withrefm.txt";
/// Fewest suppliers an enzyme needs to be built in.
const MIN_SUPPLIERS: usize = 1;

pub fn run(arg: Option<&str>) -> Result {
    let root = util::root()?;
    let path = source(arg)?;
    let bytes = std::fs::read(&path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    let rebase = parse_withrefm(&bytes).map_err(|e| format!("{}: {e}", path.display()))?;
    let version = rebase.version.clone().unwrap_or_else(|| "?".to_string());

    let mut chosen: Vec<&Enzyme> = rebase
        .enzymes
        .iter()
        .filter(|e| e.cuts())
        .filter(|e| rebase.suppliers.get(&e.name).copied().unwrap_or(0) >= MIN_SUPPLIERS)
        .collect();
    chosen.sort_by(|a, b| {
        a.name
            .to_ascii_lowercase()
            .cmp(&b.name.to_ascii_lowercase())
    });
    chosen.dedup_by(|a, b| a.name.eq_ignore_ascii_case(&b.name));
    if chosen.len() < 200 {
        return Err(format!(
            "only {} enzymes passed the rule in {}: is this a REBASE withrefm file?",
            chosen.len(),
            path.display()
        ));
    }

    let text = render(&chosen, &version, &rebase.suppliers)?;
    // Our own parser has to accept every line before the table replaces the old one.
    for e in &chosen {
        let notation = e.notation();
        Enzyme::parse(&e.name, &notation, true)
            .map_err(|err| format!("{} ({notation}) would not parse back: {err}", e.name))?;
    }
    let out = root.join(TABLE_FILE);
    std::fs::write(&out, text).map_err(|e| format!("cannot write {}: {e}", out.display()))?;
    println!(
        "xtask enzymes: {} enzymes from REBASE {version} -> {TABLE_FILE} (run cargo fmt)",
        chosen.len()
    );
    Ok(())
}

/// The REBASE file to read: the argument, `$SNAPPIESNAP_REBASE`, or the app's cached download.
fn source(arg: Option<&str>) -> Result<PathBuf> {
    if let Some(p) = arg {
        return Ok(PathBuf::from(p));
    }
    if let Some(p) = std::env::var_os("SNAPPIESNAP_REBASE") {
        return Ok(PathBuf::from(p));
    }
    let cached = data_dir().map(|d| d.join(CACHE_NAME));
    match cached {
        Some(p) if p.is_file() => Ok(p),
        Some(p) => Err(format!(
            "no REBASE file: looked at {}. Download it in the app (Enzymes tab), set \
             SNAPPIESNAP_REBASE, or give the path: cargo xtask enzymes <withrefm.txt>",
            p.display()
        )),
        None => Err("no REBASE file: give the path, or set SNAPPIESNAP_REBASE".to_string()),
    }
}

/// The app's data folder, the same place the engine caches a downloaded REBASE in.
fn data_dir() -> Option<PathBuf> {
    let home = std::env::var_os("HOME").map(PathBuf::from);
    if cfg!(target_os = "macos") {
        return home.map(|h| h.join("Library/Application Support/SnappieSnap"));
    }
    if cfg!(target_os = "windows") {
        return std::env::var_os("APPDATA").map(|a| PathBuf::from(a).join("SnappieSnap"));
    }
    std::env::var_os("XDG_DATA_HOME")
        .map(PathBuf::from)
        .or_else(|| home.map(|h| h.join(".local/share")))
        .map(|d| d.join("snappiesnap"))
}

/// The generated source file.
fn render(
    chosen: &[&Enzyme],
    version: &str,
    suppliers: &std::collections::BTreeMap<String, usize>,
) -> Result<String> {
    let mut s = String::new();
    writeln!(
        s,
        "//! The built-in enzyme set, so sites and digests work offline and before REBASE has \
         been\n\
         //! downloaded. Recognition sequences and cut positions are the published specificities \
         (facts\n\
         //! of nature, written in the standard notation).\n\
         //!\n\
         //! Selection rule (ours, applied mechanically): every restriction enzyme with a known \
         cut that\n\
         //! is sold by at least one supplier, according to the supplier codes in REBASE version \
         {version}\n\
         //! (rebase.neb.com). No other curation. A lab's enzyme set names the enzyme it buys, \
         which is\n\
         //! often a single supplier's name for a site another enzyme also cuts, so leaving \
         those out\n\
         //! lost names the lab actually uses.\n\
         //!\n\
         //! Generated by `cargo xtask enzymes`; do not edit by hand.\n\
         \n\
         use super::Enzyme;\n\
         \n\
         /// (name, site, suppliers) for the built-in enzymes, sorted by name. `suppliers` is \
         how many\n\
         /// suppliers sell the enzyme under that name, so the name shown for a site that \
         several\n\
         /// enzymes cut is the one a user is most likely to have on the shelf.\n\
         const TABLE: &[(&str, &str, u16)] = &["
    )
    .map_err(|e| e.to_string())?;
    for e in chosen {
        let notation = e.notation();
        if e.name.contains('"') || notation.contains('"') {
            return Err(format!("{}: a quote in the name or site", e.name));
        }
        let suppliers = suppliers
            .get(&e.name)
            .copied()
            .unwrap_or(0)
            .min(u16::MAX as usize);
        writeln!(s, "    (\"{}\", \"{notation}\", {suppliers}),", e.name)
            .map_err(|e| e.to_string())?;
    }
    s.push_str(
        "];\n\
         \n\
         /// The built-in enzymes, parsed. Entries that don't parse are left out rather than \
         panicking;\n\
         /// the table's test fails if there is one.\n\
         pub fn builtin() -> Vec<Enzyme> {\n\
         \x20   TABLE\n\
         \x20       .iter()\n\
         \x20       .filter_map(|(name, site, _)| Enzyme::parse(name, site, true).ok())\n\
         \x20       .collect()\n\
         }\n\
         \n\
         /// How many suppliers sell each built-in enzyme. Used to pick the name a site is \
         labelled\n\
         /// with when several enzymes cut it the same way (see `EnzymeDb::naming_rank`), so \
         that\n\
         /// works before REBASE has been downloaded.\n\
         pub fn builtin_suppliers() -> std::collections::BTreeMap<String, usize> {\n\
         \x20   TABLE\n\
         \x20       .iter()\n\
         \x20       .map(|(name, _, suppliers)| ((*name).to_string(), *suppliers as usize))\n\
         \x20       .collect()\n\
         }\n",
    );
    s.push_str(TESTS);
    Ok(s)
}

/// The tests that travel with the generated table.
const TESTS: &str = r#"
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn table_is_valid_sorted_and_unique() {
        assert_eq!(builtin().len(), TABLE.len(), "every entry parses");
        let names: Vec<String> = TABLE
            .iter()
            .map(|(n, _, _)| n.to_ascii_lowercase())
            .collect();
        let mut sorted = names.clone();
        sorted.sort();
        sorted.dedup();
        assert_eq!(names, sorted, "sorted by name, no duplicates");
        for e in builtin() {
            assert!(e.cuts(), "{} has known cuts", e.name);
        }
    }

    /// The name shown for a site several enzymes cut is the one most suppliers sell, which
    /// needs the supplier counts to travel with the table.
    #[test]
    fn supplier_counts_travel_with_the_table() {
        let counts = builtin_suppliers();
        assert_eq!(counts.len(), TABLE.len());
        let of = |n: &str| counts.get(n).copied().unwrap_or_default();
        assert!(of("HindIII") >= 4, "HindIII is sold by several suppliers");
        assert!(
            of("HindIII") > of("BstFI"),
            "the widely sold name wins over its isoschizomer"
        );
        assert!(of("EcoRI") > of("Eco47III") || of("EcoRI") >= 4);
    }

    #[test]
    fn isoschizomers_agree() {
        let all = builtin();
        let get = |n: &str| all.iter().find(|e| e.name == n).expect("a built-in enzyme");
        for (a, b) in [
            ("BsmBI", "Esp3I"),
            ("XhoI", "PaeR7I"),
            ("MboI", "Sau3AI"),
            ("MboI", "DpnII"),
            // Names a lab buys, whose site another built-in enzyme cuts the same way. These
            // were missing while the set was limited to enzymes with two or more suppliers.
            ("MscI", "MlsI"),
            ("HindII", "HincII"),
            ("MunI", "MfeI"),
            ("KspI", "SacII"),
            ("BbrPI", "PmlI"),
            ("AviII", "FspI"),
        ] {
            let (a, b) = (get(a), get(b));
            assert_eq!(
                (&a.site, a.cut_top, a.cut_bottom),
                (&b.site, b.cut_top, b.cut_bottom)
            );
        }
    }

    /// Single-supplier names a lab's own enzyme set may use are built in, so a set imported
    /// before REBASE is downloaded doesn't lose them.
    #[test]
    fn single_supplier_catalogue_names_are_built_in() {
        let all = builtin();
        for name in ["AflIII", "BsiWI", "SexAI", "XcmI", "Hpy188III", "MscI"] {
            assert!(
                all.iter().any(|e| e.name == name),
                "{name} is missing from the built-in set"
            );
        }
        // A name nobody sells any more is not built in, by the same rule: REBASE 610 lists no
        // supplier for BlnI, and AvrII and XmaJI are the sold names for its site. Importing a
        // set that uses such a name still reports it, with the hint to download REBASE.
        assert!(!all.iter().any(|e| e.name == "BlnI"));
    }
}
"#;
