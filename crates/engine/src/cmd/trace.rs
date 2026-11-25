//! Sequencing reads opened from trace files (.ab1, .scf): their quality, calls and signal, and
//! trimming of low-quality ends.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::trace::{DEFAULT_CUTOFF, quality_trim, read_quality};
use snappiesnap_doc::{Channel, Molecule, Range, Trace};

use crate::{CmdError, CommandSpec, ErrorCode, Session};

/// Most calls listed by one `trace_info` request.
const MAX_CALLS: usize = 5_000;
/// Most sample points per channel returned by one `trace_info` request.
const MAX_SIGNAL: usize = 60_000;

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::query(
            "trace_info",
            "Read Quality",
            "The sequencing trace behind the open read (a document opened from .ab1 or .scf): sample name, length, sample points, whether it is reverse-complemented, quality (mean Phred, bases of Q20 and Q30) and the stretch that trimming would keep (`cutoff`: error probability, default 0.05). With `from`/`to` (0-based residues, end exclusive) it lists each residue's base, quality and peak (sample point); `signal` adds the four channels' values over those residues' peaks. Positions and channels are in the document's orientation.",
            info,
        )
        .example(json!({}))
        .example(json!({"from": 0, "to": 10, "signal": true})),
        CommandSpec::edit(
            "trace_trim",
            "Trim Low-Quality Ends",
            "Delete the low-quality ends of the open read: every base scores `cutoff` minus its error probability, and the best-scoring stretch is kept (phred's modified Mott trimming; `cutoff` default 0.05). Undoable.",
            trim,
        )
        .menu("Actions/Trim Low-Quality Ends")
        .example(json!({}))
        .example(json!({"cutoff": 0.01})),
    ]
}

fn trace_of(m: &Molecule) -> Result<&Trace, CmdError> {
    m.trace().ok_or_else(|| {
        CmdError::new(
            ErrorCode::NotEnabled,
            "the open document has no sequencing trace",
        )
        .with_hint("open a read from an .ab1 or .scf file")
    })
}

fn cutoff(c: Option<f64>) -> Result<f64, CmdError> {
    let c = c.unwrap_or(DEFAULT_CUTOFF);
    if !(c.is_finite() && c > 0.0 && c < 1.0) {
        return Err(CmdError::invalid_params(
            "`cutoff` is an error probability between 0 and 1, such as 0.05",
        ));
    }
    Ok(c)
}

fn qualities(t: &Trace) -> Vec<Option<u8>> {
    t.calls().iter().map(|c| c.quality).collect()
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct InfoParams {
    /// First residue whose call is listed (0-based).
    from: Option<usize>,
    /// End of the listed residues (exclusive); at most 5,000 residues at a time.
    to: Option<usize>,
    /// Also return the signal under the listed residues.
    #[serde(default)]
    signal: bool,
    /// Error-probability cutoff for the suggested trim (default 0.05).
    cutoff: Option<f64>,
}

fn info(s: &Session, p: InfoParams) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let t = trace_of(m)?;
    let cutoff = cutoff(p.cutoff)?;
    let q = qualities(t);
    let keep = quality_trim(&q, cutoff);
    let mut out = json!({
        "sample_name": t.sample_name(),
        "length": m.len(),
        "samples": t.samples(),
        "reversed": t.is_reversed(),
        "quality": read_quality(&q),
        "trim": keep.map(|r| json!({"start": r.start, "end": r.end})),
        "cutoff": cutoff,
    });
    if p.from.is_none() && p.to.is_none() {
        if p.signal {
            return Err(CmdError::invalid_params(
                "`signal` needs the residues to show: give `from` and `to`",
            ));
        }
        return Ok(out);
    }
    let from = p.from.unwrap_or(0);
    let to = p.to.unwrap_or(m.len());
    if from > to || to > m.len() {
        return Err(CmdError::new(
            ErrorCode::InvalidInput,
            format!(
                "residues {from}..{to} aren't within the read (0..{})",
                m.len()
            ),
        ));
    }
    if to - from > MAX_CALLS {
        return Err(CmdError::invalid_params(format!(
            "at most {MAX_CALLS} residues at a time"
        )));
    }
    let residues = m.residues(Range::new(from, to))?;
    let calls: Vec<Value> = (from..to)
        .zip(&residues)
        .map(|(i, b)| {
            json!({
                "position": i,
                "base": char::from(*b).to_string(),
                "quality": t.calls().get(i).and_then(|c| c.quality),
                "peak": t.peak(i),
            })
        })
        .collect();
    if let Some(o) = out.as_object_mut() {
        o.insert("calls".into(), json!(calls));
    }
    if p.signal {
        let peaks: Vec<usize> = (from..to).filter_map(|i| t.peak(i)).collect();
        let signal = match (peaks.iter().min(), peaks.iter().max()) {
            (Some(&lo), Some(&hi)) => {
                // A few samples either side, so the first and last peaks are whole.
                let start = lo.saturating_sub(6);
                let end = hi.saturating_add(7).min(t.samples());
                if end.saturating_sub(start) > MAX_SIGNAL {
                    return Err(CmdError::invalid_params(format!(
                        "the signal under these residues is over {MAX_SIGNAL} samples: ask for fewer residues"
                    )));
                }
                let channel = |c: Channel| -> Vec<u16> {
                    (start..end).map(|i| t.value(c, i).unwrap_or(0)).collect()
                };
                json!({
                    "start": start,
                    "a": channel(Channel::A),
                    "c": channel(Channel::C),
                    "g": channel(Channel::G),
                    "t": channel(Channel::T),
                })
            }
            _ => Value::Null,
        };
        if let Some(o) = out.as_object_mut() {
            o.insert("signal".into(), signal);
        }
    }
    Ok(out)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct TrimParams {
    /// Error-probability cutoff (default 0.05): bases more likely wrong than this count
    /// against the stretch kept.
    cutoff: Option<f64>,
}

fn trim(s: &mut Session, p: TrimParams) -> Result<Value, CmdError> {
    let cutoff = cutoff(p.cutoff)?;
    let m = s.molecule()?;
    let len = m.len();
    let keep = quality_trim(&qualities(trace_of(m)?), cutoff).ok_or_else(|| {
        CmdError::new(
            ErrorCode::InvalidInput,
            "no stretch of the read is good enough to keep at this cutoff",
        )
        .with_hint("try a higher cutoff, such as 0.1")
    })?;
    let m = s.molecule_mut()?;
    // The end first, so the start's positions stay put.
    if keep.end < len {
        m.delete(Range::new(keep.end, len))?;
    }
    if keep.start > 0 {
        m.delete(Range::new(0, keep.start))?;
    }
    Ok(json!({
        "removed_start": keep.start,
        "removed_end": len.saturating_sub(keep.end),
        "length": keep.len(),
    }))
}
