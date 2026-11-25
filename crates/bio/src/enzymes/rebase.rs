//! Reader for REBASE's `withrefm` file (the "link_withrefm" download). REBASE is copyrighted
//! and is never bundled: users download it at runtime and we only parse their copy.
//!
//! Records are blocks of `<n>` fields: `<1>` name, `<2>` isoschizomers, `<3>` recognition
//! sequence, `<4>` methylation site, `<5>` organism, `<6>` source, `<7>` supplier codes,
//! `<8>` references. Records whose site we can't use (unknown, several specificities) and
//! methyltransferases (`M.` names, which don't cut) are skipped and counted.

use super::Enzyme;

/// Largest file we parse.
const MAX_BYTES: usize = 64 << 20;
/// Most records we keep.
const MAX_RECORDS: usize = 100_000;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RebaseError {
    #[error("not a REBASE withrefm file (no enzyme records found)")]
    NoRecords,
    #[error("file too large ({0} bytes)")]
    TooLarge(usize),
}

/// Enzymes read from a REBASE file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rebase {
    /// The `REBASE version NNN` number, if the header has one.
    pub version: Option<String>,
    pub enzymes: Vec<Enzyme>,
    /// How many suppliers sell each enzyme (distinct supplier codes in `<7>`), by name.
    pub suppliers: std::collections::BTreeMap<String, usize>,
    /// Records without a usable recognition sequence.
    pub skipped: usize,
}

#[derive(Default)]
struct Record {
    name: String,
    site: String,
    suppliers: String,
}

/// Parse a `withrefm` file. Never panics; unusable records are skipped.
pub fn parse_withrefm(bytes: &[u8]) -> Result<Rebase, RebaseError> {
    if bytes.len() > MAX_BYTES {
        return Err(RebaseError::TooLarge(bytes.len()));
    }
    let text = String::from_utf8_lossy(bytes);
    let mut version = None;
    let mut records = Vec::new();
    let mut current: Option<Record> = None;
    for line in text.lines() {
        let line = line.trim_end();
        if version.is_none()
            && let Some(rest) = line.trim_start().strip_prefix("REBASE version")
        {
            version = rest.split_whitespace().next().map(str::to_string);
        }
        let Some((field, value)) = field(line) else {
            continue;
        };
        match field {
            1 => {
                records.extend(current.take());
                if records.len() >= MAX_RECORDS {
                    break;
                }
                current = Some(Record {
                    name: value.trim().to_string(),
                    ..Record::default()
                });
            }
            3 => {
                if let Some(r) = current.as_mut() {
                    r.site = value.trim().to_string();
                }
            }
            7 => {
                if let Some(r) = current.as_mut() {
                    r.suppliers = value.trim().to_string();
                }
            }
            _ => {}
        }
    }
    records.extend(current);
    if records.is_empty() {
        return Err(RebaseError::NoRecords);
    }
    let mut enzymes = Vec::with_capacity(records.len());
    let mut skipped = 0;
    let mut suppliers = std::collections::BTreeMap::new();
    for r in records {
        let mut codes: Vec<char> = r
            .suppliers
            .chars()
            .filter(char::is_ascii_alphanumeric)
            .collect();
        codes.sort_unstable();
        codes.dedup();
        if !codes.is_empty() {
            suppliers.insert(r.name.clone(), codes.len());
        }
        let usable =
            !r.name.is_empty() && !r.name.starts_with("M.") && !r.site.contains(['?', ',']);
        match usable.then(|| Enzyme::parse(&r.name, &r.site, !r.suppliers.is_empty())) {
            Some(Ok(e)) => enzymes.push(e),
            _ => skipped += 1,
        }
    }
    Ok(Rebase {
        version,
        enzymes,
        suppliers,
        skipped,
    })
}

/// `<n>value` → (n, value).
fn field(line: &str) -> Option<(u8, &str)> {
    let rest = line.strip_prefix('<')?;
    let (n, value) = rest.split_once('>')?;
    Some((n.parse().ok()?, value))
}

#[cfg(test)]
mod tests {
    use super::*;

    // A made-up file in the withrefm layout (not REBASE data).
    const SAMPLE: &str = "
REBASE version 999                                              withrefm.999

<ENZYME NAME>   Restriction enzyme name.
<1>AaaI
<2>BbbI
<3>G^AATTC
<4>
<5>Some organism
<6>Someone
<7>XY
<8>A reference.

<1>CccI
<2>
<3>GGTCTC(1/5)
<4>
<5>
<6>
<7>
<8>

<1>DddI
<3>?
<7>

<1>EeeI
<3>GGCC,CCGG
<7>

<1>FffI
<3>(8/13)GACNNNNNNTGG(12/7)
<7>Z

<1>M.GggI
<3>GATC
<7>
";

    #[test]
    fn reads_records() {
        let r = parse_withrefm(SAMPLE.as_bytes()).unwrap();
        assert_eq!(r.version.as_deref(), Some("999"));
        assert_eq!(r.skipped, 3);
        let names: Vec<(&str, bool)> = r
            .enzymes
            .iter()
            .map(|e| (e.name.as_str(), e.commercial))
            .collect();
        assert_eq!(names, vec![("AaaI", true), ("CccI", false), ("FffI", true)]);
        assert_eq!(r.enzymes[1].cut_top, Some(7));
        assert_eq!(r.enzymes[2].cut_upstream, Some((-8, -13)));
    }

    #[test]
    fn hostile_input() {
        assert_eq!(parse_withrefm(b""), Err(RebaseError::NoRecords));
        assert_eq!(
            parse_withrefm(b"<html>not found</html>"),
            Err(RebaseError::NoRecords)
        );
        for cut in 0..SAMPLE.len() {
            let _ = parse_withrefm(&SAMPLE.as_bytes()[..cut]);
        }
        let _ = parse_withrefm(&[0xff, b'<', b'1', b'>', 0xfe, b'\n', b'<', b'3', b'>', 0x80]);
    }
}
