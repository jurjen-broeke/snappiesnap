//! Engine-level checks over the local corpus (`cargo xtask corpus-test`): oligo and enzyme
//! lists exported from other tools, read the way the import commands read them.
//!
//! The corpus is local-only and partial, and the files in it are a user's own working files, so
//! these tests name no file and pin no sequence. They discover whatever is present and check
//! the properties that must hold for any such list: nothing is dropped silently, a name is
//! never attached to the wrong sequence, and case is preserved.

#![cfg(feature = "corpus")]
// Test code: the never-crash lints apply to production code only (see clippy.toml).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::{Path, PathBuf};

use snappiesnap_engine::oligos::parse_oligo_text;

fn corpus_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../corpus")
}

/// A corpus is local and partial by nature, so a test whose material is absent says so and
/// skips instead of failing.
fn read(rel: &str) -> Option<Vec<u8>> {
    let path = corpus_dir().join(rel);
    match std::fs::read(&path) {
        Ok(bytes) => Some(bytes),
        Err(_) => {
            println!("no {rel} in the corpus: nothing to check");
            None
        }
    }
}

/// Every file directly under `corpus/<dir>` whose name ends in `ext`, sorted, so a run is
/// deterministic. Empty when the folder is absent.
fn files_in(dir: &str, ext: &str) -> Vec<PathBuf> {
    let Ok(entries) = std::fs::read_dir(corpus_dir().join(dir)) else {
        println!("no corpus/{dir}: nothing to check");
        return Vec::new();
    };
    let mut paths: Vec<PathBuf> = entries
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_file() && p.extension().is_some_and(|e| e.eq_ignore_ascii_case(ext)))
        .collect();
    paths.sort();
    paths
}

/// A file's name without its extension: the name the import gives a bare sequence.
fn stem(path: &Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default()
}

/// Cloning primers the way a lab keeps them: one file each, holding the bare sequence, named
/// only by the file name. Every one must come back under that name, with the file's own
/// sequence and its case intact — a name must never end up on another file's sequence.
#[test]
fn a_primer_file_is_read_under_its_own_name() {
    let mut checked = 0;
    for path in files_in("cloning", "txt") {
        let name = stem(&path);
        let bytes = std::fs::read(&path).unwrap();
        // Only single-sequence files belong to this check; a list in the same folder is read
        // as a list and is covered by the export checks below.
        let Ok((oligos, skipped)) = parse_oligo_text(&bytes, Some(&name)) else {
            continue;
        };
        if oligos.len() != 1 || !skipped.is_empty() {
            continue;
        }
        assert_eq!(oligos[0].name, name, "{}", path.display());
        // The sequence is exactly the file's letters, case and all.
        let expected: String = String::from_utf8_lossy(&bytes)
            .chars()
            .filter(|c| c.is_ascii_alphabetic())
            .collect();
        assert_eq!(oligos[0].residues, expected, "{}", path.display());
        checked += 1;
    }
    println!("{checked} primer files read under their own name");
}

/// An oligo list exported from another tool, in whatever encoding and line ending it uses.
/// Every oligo must carry a name and a sequence, the totals must add up between the file and
/// the library, and a name repeated in the file must keep the sequence it first had.
#[test]
fn an_exported_oligo_list_loses_no_oligo() {
    let Some(bytes) = read("import/exported-oligo-list.txt") else {
        return;
    };
    let (oligos, skipped) = parse_oligo_text(&bytes, Some("exported-oligo-list")).unwrap();
    assert!(!oligos.is_empty(), "the list holds oligos");

    // Every oligo has a name and a sequence, and the sequence is nothing but letters: a
    // mangled line would show up here rather than as a silently wrong oligo.
    for o in &oligos {
        assert!(!o.name.trim().is_empty(), "an oligo with no name");
        assert!(!o.residues.is_empty(), "{}: no sequence", o.name);
        assert!(
            o.residues.bytes().all(|b| b.is_ascii_alphabetic()),
            "{}: {} is not a sequence",
            o.name,
            o.residues
        );
    }

    // Names repeated in the file are kept once, with the sequence they first had, and the
    // totals add up: nothing vanishes between the file and the library.
    let mut library = snappiesnap_engine::oligos::OligoLibrary::default();
    let report = library.add(oligos.clone()).unwrap();
    assert!(report.conflicts.is_empty(), "{:?}", report.conflicts);
    assert_eq!(
        report.added + report.duplicates,
        oligos.len(),
        "every oligo is added or counted as a duplicate"
    );
    assert_eq!(library.oligos.len(), report.added);
    println!(
        "{} oligos read, {} added, {} duplicates, {} lines reported",
        oligos.len(),
        report.added,
        report.duplicates,
        skipped.len()
    );
}

/// An enzyme set exported from another tool: a title line, then the names wrapped over several
/// lines with trailing commas and whatever line ending the tool used. Every name in it must be
/// usable as an enzyme name — a name broken across a line end must not be split or joined,
/// because a mangled name would silently drop an enzyme from the user's own set.
#[test]
fn an_exported_enzyme_set_reads_every_name() {
    let Some(bytes) = read("import/exported-enzyme-set.txt") else {
        return;
    };
    let (name, enzymes) = snappiesnap_engine::enzymes::parse_set_text(&bytes).unwrap();
    assert!(!enzymes.is_empty(), "the set holds enzyme names");
    // A split or joined name shows up as a name that no enzyme could have.
    for e in &enzymes {
        assert!(
            !e.trim().is_empty()
                && e.chars().all(|c| c.is_ascii_alphanumeric() || c == '-')
                && e.starts_with(|c: char| c.is_ascii_alphabetic()),
            "{e} does not look like an enzyme name (a line end may have split it)"
        );
    }

    // None is rejected as malformed by the set the user would get.
    let set_name = name.as_deref().unwrap_or("imported set");
    let (set, rejected) = snappiesnap_engine::enzymes::UserSet::new(set_name, &enzymes).unwrap();
    assert!(rejected.is_empty(), "{rejected:?}");

    // The built-in database is a selection of REBASE by supplier count, so it knows one name of
    // a pair of isoschizomers and not the other (MscI/MluNI, HindII/HincII, SacII/KspI and the
    // like). A set naming the other one reports those as missing, with the hint to download
    // REBASE; that is expected, so only the count is printed.
    let db = snappiesnap_engine::EnzymeDb::default();
    assert!(db.is_builtin());
    let missing = set.missing(&db);
    println!(
        "{} enzymes in the set, {} of them in the built-in database",
        enzymes.len(),
        enzymes.len() - missing.len()
    );
}

/// The primer files and the constructs they belong to must agree: a file named after a primer
/// holds the sequence the construct stores under that name, case and all. This is what makes
/// the file name a usable name — if the two ever disagree, importing a folder of them would
/// put a wrong sequence in the library under a name the maps already use.
#[test]
fn primer_files_match_the_primers_of_the_constructs_that_use_them() {
    let mut compared = 0;
    for path in files_in("cloning", "dna") {
        let bytes = std::fs::read(&path).unwrap();
        let rel = format!(
            "cloning/{}",
            path.file_name().unwrap_or_default().to_string_lossy()
        );
        let doc = snappiesnap_io::import(&rel, &bytes).unwrap();
        for primer in doc.molecule.primers() {
            let Some(file) = read(&format!("cloning/{}.txt", primer.name)) else {
                continue;
            };
            let (oligos, _) = parse_oligo_text(&file, Some(primer.name.as_str())).unwrap();
            assert_eq!(oligos.len(), 1, "{}", primer.name);
            assert_eq!(
                oligos[0].residues,
                primer.seq.as_str(),
                "{} in {rel}",
                primer.name
            );
            compared += 1;
        }
    }
    println!("{compared} primer files checked against their constructs");
}
