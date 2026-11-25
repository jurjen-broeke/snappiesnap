//! BLAST searches at NCBI: the web address of a BLAST form filled in with part of the
//! document. The command itself sends nothing anywhere; the app opens the address in the
//! user's browser, where the search starts when the user presses BLAST: a sequence leaves the
//! machine only on the user's own request.
//!
//! Programs (NCBI BLAST, Altschul et al. 1990 and 1997): DNA against a nucleotide database is
//! blastn, against a protein database blastx (the DNA translated), against a translated
//! nucleotide database tblastx; a protein against a protein database is blastp, against a
//! translated nucleotide database tblastn.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::translate::strand_residues;
use snappiesnap_doc::{Range, Strand};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

/// The BLAST form NCBI documents for links ("BLAST URL API").
const BLAST: &str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";
/// Longest query sent in a web address: servers commonly refuse request lines over 8 KB.
pub const MAX_QUERY: usize = 7_500;

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "blast_search",
            "BLAST Selection at NCBI",
            "The address of NCBI's BLAST form filled in with the selection (or `range`; the whole sequence if there is neither), for opening in a browser: nothing is sent by this command. `against`: \"nucleotide\" (DNA: blastn), \"protein\" (DNA: blastx; protein: blastp) or \"translated\" (a translated nucleotide database: tblastx for DNA, tblastn for protein); default nucleotide for DNA and protein for proteins. `database` names an NCBI database (for example nt, nr, refseq_rna); otherwise the form's default. Queries longer than 7,500 residues don't fit in a web address.",
            blast_search,
        )
        .menu("Sequence/BLAST Selection at NCBI")
        .example(json!({}))
        .example(json!({"against": "protein", "range": {"start": 0, "end": 9}}))
        .example(json!({"against": "translated", "database": "nt"})),
    ]
}

#[derive(Deserialize, JsonSchema, Clone, Copy, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
enum Against {
    Nucleotide,
    Protein,
    Translated,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Params {
    /// The database type to search (see the description).
    against: Option<Against>,
    /// The stretch to search; default the selection, else the whole sequence.
    range: Option<Range>,
    /// An NCBI database name, such as nt, nr or refseq_rna.
    database: Option<String>,
}

/// Percent-encode for a query string value (RFC 3986 unreserved characters kept).
fn encode(s: &str) -> String {
    let mut out = String::new();
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || matches!(b, b'-' | b'_' | b'.' | b'~') {
            out.push(char::from(b));
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}

fn blast_search(s: &Session, p: Params) -> Result<Value, CmdError> {
    let doc = s.document().ok_or_else(|| {
        CmdError::new(ErrorCode::NoDocument, "no document is open")
            .with_hint("open or create a sequence first")
    })?;
    let m = &doc.molecule;
    if m.is_empty() {
        return Err(CmdError::invalid_params("the sequence is empty"));
    }
    let nucleotide = m.seq().alphabet().is_nucleotide();
    let against = p.against.unwrap_or(if nucleotide {
        Against::Nucleotide
    } else {
        Against::Protein
    });
    let program = match (nucleotide, against) {
        (true, Against::Nucleotide) => "blastn",
        (true, Against::Protein) => "blastx",
        (true, Against::Translated) => "tblastx",
        (false, Against::Protein) => "blastp",
        (false, Against::Translated) => "tblastn",
        (false, Against::Nucleotide) => {
            return Err(CmdError::invalid_params(
                "a protein can't be searched against nucleotides directly",
            )
            .with_hint(
                "use `against: \"translated\"` (tblastn) to search translated nucleotide databases",
            ));
        }
    };
    let range = match p.range.or(doc.selection.filter(|r| !r.is_empty())) {
        Some(r) => r.validate(m.len(), m.topology())?,
        None => Range::new(0, m.len()),
    };
    if range.is_empty() {
        return Err(CmdError::invalid_params(
            "the range is empty: select some residues",
        ));
    }
    let query = strand_residues(m, range, Strand::Forward)
        .ok_or_else(|| CmdError::invalid_params("range outside the molecule"))?;
    let query: String = String::from_utf8_lossy(&query)
        .chars()
        .filter(|c| c.is_ascii_alphabetic() || *c == '*')
        .collect::<String>()
        .trim_end_matches('*')
        .to_ascii_uppercase();
    if query.len() > MAX_QUERY {
        return Err(CmdError::invalid_params(format!(
            "the query is {} residues; at most {MAX_QUERY} fit in a web address",
            query.len()
        ))
        .with_hint(
            "select a shorter stretch, or paste the sequence into the BLAST form yourself",
        ));
    }
    if let Some(db) = &p.database
        && (db.is_empty()
            || db.len() > 64
            || !db
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'/'))
    {
        return Err(CmdError::invalid_params(
            "`database` must be an NCBI database name (letters, digits, _ and /)",
        ));
    }
    let (a, b) = (
        range.start + 1,
        if range.end == 0 { m.len() } else { range.end },
    );
    let title = format!("{} {a}..{b}", m.name());
    let mut url = format!(
        "{BLAST}?PAGE_TYPE=BlastSearch&PROGRAM={program}&QUERY={}&JOB_TITLE={}",
        encode(&query),
        encode(&title)
    );
    if let Some(db) = &p.database {
        url.push_str(&format!("&DATABASE={}", encode(db)));
    }
    Ok(json!({
        "url": url,
        "program": program,
        "query_length": query.len(),
        "range": range,
        "title": title,
    }))
}

#[cfg(test)]
mod tests {
    use super::encode;

    #[test]
    fn query_values_are_encoded() {
        assert_eq!(encode("pUC19 1..10"), "pUC19%201..10");
        assert_eq!(encode("a&b=c/é"), "a%26b%3Dc%2F%C3%A9");
    }
}
