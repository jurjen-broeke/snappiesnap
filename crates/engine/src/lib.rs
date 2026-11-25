//! SnappieSnap `engine` (layer L4).
//!
//! The [`Session`] holds the open document, undo/redo history and a log of recent problems.
//! Every user-visible action is a command ([`CommandSpec`]) in the [`Registry`]; the UI, CLI,
//! control channel and MCP all run commands through [`Session::execute`], which:
//!
//! 1. parses and checks the JSON params against the command's schema,
//! 2. refuses commands that aren't enabled in the current state,
//! 3. runs the command under `catch_unwind`, and
//! 4. on any error or panic, restores the document exactly as it was.
//!
//! Layering is enforced by `cargo xtask layers`; see `AGENTS.md`.

pub mod api;
mod cmd;
pub use cmd::copy_text;
pub mod enzymes;
mod error;
pub mod experiment;
pub mod gel;
#[cfg(not(target_arch = "wasm32"))]
mod http;
mod log;
pub mod matched;
pub mod oligos;
pub mod paths;
mod platform;
mod registry;
pub mod sample;
pub mod script;
mod session;
pub mod text;

pub use enzymes::EnzymeDb;
pub use error::{CmdError, ErrorCode};
pub use log::{Level, LogEntry};
#[cfg(not(target_arch = "wasm32"))]
pub use platform::NativeFs;
pub use platform::{DirEntry, MAX_ENTRIES, MAX_FILE, MemoryFs, Platform};
pub use registry::{CommandInfo, CommandSpec, Effect, NoParams, Registry};
pub use session::{Document, History, JournalEntry, MAX_UNDO, OpenDocument, Session};

pub use snappiesnap_bio as bio;
pub use snappiesnap_doc as doc;
pub use snappiesnap_io as io;
