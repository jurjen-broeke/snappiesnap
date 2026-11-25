//! Adding annotations from a file to the open document: GFF3 files that hold features
//! without a sequence (as NCBI and genome browsers export them).

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

use crate::{CmdError, CommandSpec, ErrorCode, Level, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::edit(
            "import_annotations",
            "Import Annotations",
            "Add the features of a GFF3 file (`path`) to the open document as one undoable edit. The rows of `seqid` are used; without it, those of the file's only sequence, or of the one whose region has the document's length. Features already in the document (same name, type and place) are skipped; rows outside the sequence are listed in `warnings`. GFF3 files with a ##FASTA section can also be opened as documents.",
            import_annotations,
        )
        .menu("Features/Import Annotations…")
        .example(json!({"path": "annotations.gff3"}))
        .example(json!({"path": "annotations.gff3", "seqid": "demo"})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct Params {
    /// The GFF3 file.
    path: String,
    /// Which sequence's rows to use, when the file describes several.
    seqid: Option<String>,
}

fn import_annotations(s: &mut Session, p: Params) -> Result<Value, CmdError> {
    let bytes = s
        .platform()
        .read(&p.path)
        .map_err(|e| CmdError::new(ErrorCode::Io, format!("cannot read {}: {e}", p.path)))?;
    if !bytes.starts_with(b"##gff-version 3") {
        return Err(CmdError::invalid_params(format!(
            "{} isn't a GFF3 file (it must start with ##gff-version 3)",
            p.path
        ))
        .with_hint("GenBank, EMBL and .dna files carry their own sequence: open them instead"));
    }
    let (len, topology) = {
        let m = s.molecule()?;
        (m.len(), m.topology())
    };
    let (seqid, features, warnings) =
        snappiesnap_io::gff3_annotations(&bytes, p.seqid.as_deref(), len, topology)?;
    let m = s.molecule_mut()?;
    let mut ids = Vec::new();
    let mut skipped = 0usize;
    for f in features {
        let present = m.features().iter().any(|e| {
            e.name == f.name
                && e.kind == f.kind
                && e.segments
                    .iter()
                    .map(|x| x.range)
                    .eq(f.segments.iter().map(|x| x.range))
        });
        if present {
            skipped += 1;
            continue;
        }
        ids.push(m.add_feature(f)?);
    }
    for w in &warnings {
        s.log(
            Level::Warn,
            "import_annotations",
            format!("{}: {w}", p.path),
        );
    }
    Ok(json!({
        "seqid": seqid,
        "added": ids.len(),
        "ids": ids,
        "skipped": skipped,
        "warnings": warnings,
    }))
}
