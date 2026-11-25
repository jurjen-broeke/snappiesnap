//! Sequencing reads matched to a construct and stored in its `.dna` file (see
//! `snappiesnap_io::MatchedRead`): listed, opened in tabs of their own, and checked against the
//! construct with `check_reads`.

use snappiesnap_io::{MatchedRead, Source};

use crate::{CmdError, ErrorCode, Session};

/// The reads matched to document `id` (default: the active one); empty if its file has none.
pub fn matched_reads(s: &Session, id: Option<u64>) -> Result<Vec<MatchedRead>, CmdError> {
    let doc = match id {
        Some(id) => s.document_by_id(id).ok_or_else(|| {
            CmdError::new(ErrorCode::NotFound, format!("document {id} is not open"))
                .with_hint("the `documents` command lists the open documents and their ids")
        })?,
        None => s.document().ok_or_else(CmdError::no_document)?,
    };
    Ok(match &doc.source {
        Some(Source::Dna(src)) => src.matched_reads(),
        _ => Vec::new(),
    })
}
