use std::fmt;

use serde::Serialize;
use snappiesnap_doc::DocError;

/// Machine-readable error categories, shared by the UI, CLI, control channel and MCP.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ErrorCode {
    /// No command with that id.
    UnknownCommand,
    /// No protocol method with that name.
    UnknownMethod,
    /// The params don't match the command's schema.
    InvalidParams,
    /// The command needs an open document.
    NoDocument,
    /// The command can't run in the current state.
    NotEnabled,
    /// The params are well-formed but the edit is impossible (out of range, invalid residue).
    InvalidInput,
    /// A feature, primer or other item with that id doesn't exist.
    NotFound,
    /// Not implemented yet.
    Unsupported,
    /// A file couldn't be read or written (missing, no permission, disk full, ...).
    Io,
    /// A script expectation (`expect`, `expect_error`) didn't hold.
    AssertionFailed,
    /// A bug: the command panicked. The document was restored.
    Internal,
}

/// An error a person or agent can act on: what went wrong, and where possible how to fix it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct CmdError {
    pub code: ErrorCode,
    pub message: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hint: Option<String>,
}

impl CmdError {
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        CmdError {
            code,
            message: message.into(),
            hint: None,
        }
    }

    pub fn with_hint(mut self, hint: impl Into<String>) -> Self {
        self.hint = Some(hint.into());
        self
    }

    pub fn invalid_params(message: impl Into<String>) -> Self {
        CmdError::new(ErrorCode::InvalidParams, message)
    }

    pub fn no_document() -> Self {
        CmdError::new(ErrorCode::NoDocument, "no document is open")
            .with_hint("create one with the `new` command")
    }

    pub fn unsupported(message: impl Into<String>) -> Self {
        CmdError::new(ErrorCode::Unsupported, message)
    }
}

impl fmt::Display for CmdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)?;
        if let Some(hint) = &self.hint {
            write!(f, " (hint: {hint})")?;
        }
        Ok(())
    }
}

impl std::error::Error for CmdError {}

impl From<DocError> for CmdError {
    fn from(e: DocError) -> Self {
        let code = match e {
            DocError::FeatureNotFound(_) | DocError::PrimerNotFound(_) => ErrorCode::NotFound,
            _ => ErrorCode::InvalidInput,
        };
        let err = CmdError::new(code, e.to_string());
        match e {
            DocError::FeatureSpansOrigin { .. } => {
                err.with_hint("use `rotate_origin` or edit the feature, then try again")
            }
            DocError::FeatureNotFound(_) => {
                err.with_hint("`doc.get` lists the features and their ids")
            }
            _ => err,
        }
    }
}

impl From<snappiesnap_doc::SeqError> for CmdError {
    fn from(e: snappiesnap_doc::SeqError) -> Self {
        DocError::from(e).into()
    }
}

impl From<snappiesnap_io::IoError> for CmdError {
    fn from(e: snappiesnap_io::IoError) -> Self {
        use snappiesnap_io::IoError;
        match e {
            IoError::Doc(d) => d.into(),
            IoError::Unsupported(_) => CmdError::new(ErrorCode::Unsupported, e.to_string()),
            IoError::UnknownFormat => CmdError::new(ErrorCode::InvalidInput, e.to_string())
                .with_hint("open a SnapGene .dna, GenBank or FASTA file"),
            IoError::Dna(_) | IoError::Format(_) => {
                CmdError::new(ErrorCode::InvalidInput, e.to_string())
            }
        }
    }
}
