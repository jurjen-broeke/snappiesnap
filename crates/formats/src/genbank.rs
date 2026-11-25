//! GenBank flat files (INSDC feature table format).
//!
//! Records keep their raw text: the header (LOCUS to FEATURES), each feature entry, and the
//! sequence. Unchanged parts are written back verbatim; changed parts are written in the layout
//! of the INSDC specification as the reference application exports it (80-column lines,
//! qualifiers from column 22, text values wrapped at spaces, ORIGIN in groups of ten).

use std::fmt::Write as _;

/// Why a GenBank file couldn't be read.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum GenBankError {
    #[error("not a GenBank file: it must start with a LOCUS line")]
    NotGenBank,

    #[error("the file isn't valid UTF-8 text")]
    NotText,

    #[error("line {line}: {reason}")]
    Malformed { line: usize, reason: String },
}

/// Column where feature locations and qualifiers start (0-based).
pub const QUALIFIER_INDENT: usize = 21;
/// Longest line written.
pub const LINE_WIDTH: usize = 80;

/// The LOCUS line's fields.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub name: String,
    pub length: u64,
    /// For example `DNA`, `ds-DNA`, `RNA`, `mRNA`, `aa` (proteins use `aa` instead of `bp`).
    pub molecule: String,
    pub circular: bool,
    pub division: String,
    pub date: String,
    pub protein: bool,
}

impl Locus {
    pub fn parse(line: &str) -> Option<Locus> {
        let mut t = line.split_whitespace();
        if t.next()? != "LOCUS" {
            return None;
        }
        let name = t.next()?.to_string();
        let length = t.next()?.parse().ok()?;
        let unit = t.next()?;
        let protein = unit == "aa";
        let rest: Vec<&str> = t.collect();
        let circular = rest.iter().any(|w| w.eq_ignore_ascii_case("circular"));
        let molecule = rest
            .iter()
            .find(|w| !w.eq_ignore_ascii_case("circular") && !w.eq_ignore_ascii_case("linear"))
            .filter(|w| !is_date(w))
            .map_or(String::new(), |w| w.to_string());
        let date = rest
            .iter()
            .find(|w| is_date(w))
            .map_or(String::new(), |w| w.to_string());
        let division = rest
            .iter()
            .rev()
            .filter(|w| !is_date(w))
            .find(|w| w.len() == 3 && w.chars().all(|c| c.is_ascii_uppercase()) && **w != molecule)
            .map_or(String::new(), |w| w.to_string());
        Some(Locus {
            name,
            length,
            molecule,
            circular,
            division,
            date,
            protein,
        })
    }

    /// The LOCUS line in the column layout of the INSDC specification.
    pub fn to_line(&self) -> String {
        let (strand, molecule) = match self.molecule.split_once('-') {
            Some((s, m)) if matches!(s, "ss" | "ds" | "ms") => (format!("{s}-"), m.to_string()),
            _ => (String::new(), self.molecule.clone()),
        };
        let unit = if self.protein { "aa" } else { "bp" };
        let topology = if self.circular { "circular" } else { "linear" };
        format!(
            "LOCUS       {:<16} {:>11} {unit} {:<3}{:<6}  {:<8} {:<3} {}",
            self.name, self.length, strand, molecule, topology, self.division, self.date
        )
        .trim_end()
        .to_string()
    }
}

/// Patch the length and topology in a LOCUS line, keeping the line's own spacing: the length
/// stays right-aligned where it was (taking or giving spaces on its left), and the topology word
/// keeps an 8-column field. `None` if anything else changed or the line has an unexpected shape.
fn patch_locus_line(raw: &str, old: &Locus, new: &Locus) -> Option<String> {
    if old.name != new.name
        || old.molecule != new.molecule
        || old.division != new.division
        || old.date != new.date
        || old.protein != new.protein
    {
        return None;
    }
    let (line, newline) = match raw.strip_suffix('\n') {
        Some(l) => (l, "\n"),
        None => (raw, ""),
    };
    let unit = if old.protein { " aa" } else { " bp" };
    let old_len = format!("{}{unit}", old.length);
    let at = line.find(&format!(" {old_len}"))? + 1;
    let before = line.get(..at)?;
    let after = line.get(at + old_len.len()..)?;
    let new_len = format!("{}{unit}", new.length);
    let mut before = before.to_string();
    let grow = new_len.len() as isize - old_len.len() as isize;
    if grow > 0 {
        let spare = before.len() - before.trim_end().len();
        let take = (grow as usize).min(spare.saturating_sub(1));
        before.truncate(before.len() - take);
    } else {
        before.push_str(&" ".repeat((-grow) as usize));
    }
    let mut after = after.to_string();
    if old.circular != new.circular {
        let (from, to) = if new.circular {
            ("linear  ", "circular")
        } else {
            ("circular", "linear  ")
        };
        let from_short = from.trim_end();
        if let Some(i) = after.find(from) {
            after.replace_range(i..i + from.len(), to);
        } else {
            let i = after.find(from_short)?;
            after.replace_range(i..i + from_short.len(), to.trim_end());
        }
    }
    Some(format!("{before}{new_len}{after}{newline}"))
}

fn is_date(w: &str) -> bool {
    let b = w.as_bytes();
    b.len() == 11 && b.get(2) == Some(&b'-') && b.get(6) == Some(&b'-')
}

/// A qualifier: `/name=value`, `/name="value"` or `/name`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Qualifier {
    pub name: String,
    pub value: Option<String>,
    /// Whether the value was written in quotes.
    pub quoted: bool,
}

/// One feature table entry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureEntry {
    pub key: String,
    /// The location as written, without line breaks.
    pub location: String,
    pub qualifiers: Vec<Qualifier>,
}

impl FeatureEntry {
    pub fn qualifier(&self, name: &str) -> Option<&str> {
        self.qualifiers
            .iter()
            .find(|q| q.name == name)
            .and_then(|q| q.value.as_deref())
    }

    /// The entry as text, ending with a newline.
    pub fn to_text(&self) -> String {
        let mut out = String::new();
        let first = format!("     {:<15} ", self.key);
        let indent = " ".repeat(QUALIFIER_INDENT);
        for (i, line) in wrap_location(&self.location).iter().enumerate() {
            out.push_str(if i == 0 { &first } else { &indent });
            out.push_str(line);
            out.push('\n');
        }
        for q in &self.qualifiers {
            for line in qualifier_lines(q) {
                out.push_str(&indent);
                out.push_str(&line);
                out.push('\n');
            }
        }
        out
    }
}

/// Content width after the qualifier indent.
const CONTENT: usize = LINE_WIDTH - QUALIFIER_INDENT;

/// Break a location at commas so lines fit.
fn wrap_location(location: &str) -> Vec<String> {
    let mut lines = vec![String::new()];
    for (i, part) in location.split(',').enumerate() {
        let piece = if i == 0 {
            part.to_string()
        } else {
            format!(",{part}")
        };
        let fits = lines
            .last()
            .is_some_and(|l| l.len() + piece.len() <= CONTENT);
        if !fits && lines.last().is_some_and(|l| !l.is_empty()) {
            // The comma stays at the end of the line it follows.
            if let Some(l) = lines.last_mut()
                && i > 0
            {
                l.push(',');
            }
            lines.push(part.to_string());
        } else if let Some(l) = lines.last_mut() {
            l.push_str(&piece);
        }
    }
    lines
}

/// A qualifier as wrapped lines (without the indent). Text wraps after spaces (which stay at
/// the end of the line); text without spaces (such as translations) is cut at the line width;
/// embedded newlines start a new line.
fn qualifier_lines(q: &Qualifier) -> Vec<String> {
    let text = match (&q.value, q.quoted) {
        (None, _) => format!("/{}", q.name),
        (Some(v), true) => format!("/{}=\"{}\"", q.name, v.replace('"', "\"\"")),
        (Some(v), false) => format!("/{}={v}", q.name),
    };
    let mut lines = Vec::new();
    for paragraph in text.split('\n') {
        let mut line = String::new();
        for word in split_keep_spaces(paragraph) {
            if line.len() + word.len() <= CONTENT {
                line.push_str(word);
                continue;
            }
            if !line.is_empty() {
                lines.push(std::mem::take(&mut line));
            }
            let mut rest = word;
            while rest.len() > CONTENT {
                let cut = floor_char_boundary(rest, CONTENT);
                if cut == 0 {
                    break;
                }
                lines.push(rest.get(..cut).unwrap_or("").to_string());
                rest = rest.get(cut..).unwrap_or("");
            }
            line.push_str(rest);
        }
        lines.push(line);
    }
    lines
}

fn floor_char_boundary(s: &str, max: usize) -> usize {
    (0..=max.min(s.len()))
        .rev()
        .find(|&i| s.is_char_boundary(i))
        .unwrap_or(0)
}

/// Split into words, each keeping the spaces that follow it.
fn split_keep_spaces(s: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut start = 0usize;
    let bytes = s.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        if bytes.get(i) == Some(&b' ') {
            while bytes.get(i) == Some(&b' ') {
                i += 1;
            }
            out.push(s.get(start..i).unwrap_or(""));
            start = i;
        } else {
            i += 1;
        }
    }
    if start < s.len() {
        out.push(s.get(start..).unwrap_or(""));
    }
    out
}

/// The ORIGIN block lines (without the ORIGIN line): 60 residues per line in groups of ten,
/// numbered, lower case.
pub fn origin_text(residues: &[u8]) -> String {
    let mut out = String::new();
    for (i, line) in residues.chunks(60).enumerate() {
        let _ = write!(out, "{:>9}", i * 60 + 1);
        for group in line.chunks(10) {
            out.push(' ');
            out.push_str(&String::from_utf8_lossy(group).to_ascii_lowercase());
        }
        out.push('\n');
    }
    out
}

/// A part of a record's text.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub entry: FeatureEntry,
    /// The original text while unmodified.
    pub raw: Option<String>,
}

/// One record (LOCUS ... //).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub locus: Locus,
    locus_raw: Option<String>,
    /// Lines after LOCUS up to (not including) the FEATURES line.
    pub header: String,
    /// The FEATURES line, if any.
    features_line: Option<String>,
    pub features: Vec<Feature>,
    /// Lines between the features and ORIGIN (BASE COUNT, CONTIG, ...).
    pub middle: String,
    /// The ORIGIN line, if any.
    origin_line: Option<String>,
    pub residues: Vec<u8>,
    residues_raw: Option<String>,
    /// The `//` line and anything up to the next record.
    pub tail: String,
}

impl Record {
    /// A new record with a minimal header.
    pub fn new(locus: Locus, header: String) -> Self {
        Record {
            locus,
            locus_raw: None,
            header,
            features_line: Some("FEATURES             Location/Qualifiers\n".into()),
            features: Vec::new(),
            middle: String::new(),
            origin_line: Some("ORIGIN\n".into()),
            residues: Vec::new(),
            residues_raw: None,
            tail: "//\n".into(),
        }
    }

    /// Change the LOCUS fields. When only the length or topology changed, the original line is
    /// patched in place so its spacing (which varies between writers) is kept.
    pub fn set_locus(&mut self, locus: Locus) {
        if locus == self.locus {
            return;
        }
        let patched = self
            .locus_raw
            .as_deref()
            .and_then(|raw| patch_locus_line(raw, &self.locus, &locus));
        self.locus = locus;
        self.locus_raw = patched;
    }

    pub fn set_residues(&mut self, residues: Vec<u8>) {
        if residues != self.residues {
            self.residues = residues;
            self.residues_raw = None;
            // A base count would be stale.
            self.middle = self
                .middle
                .lines()
                .filter(|l| !l.starts_with("BASE COUNT"))
                .map(|l| format!("{l}\n"))
                .collect();
        }
    }

    pub fn to_text(&self) -> String {
        let mut out = String::new();
        match &self.locus_raw {
            Some(raw) => out.push_str(raw),
            None => {
                out.push_str(&self.locus.to_line());
                out.push('\n');
            }
        }
        out.push_str(&self.header);
        if self.features_line.is_some() || !self.features.is_empty() {
            out.push_str(
                self.features_line
                    .as_deref()
                    .unwrap_or("FEATURES             Location/Qualifiers\n"),
            );
        }
        for f in &self.features {
            match &f.raw {
                Some(raw) => out.push_str(raw),
                None => out.push_str(&f.entry.to_text()),
            }
        }
        out.push_str(&self.middle);
        if let Some(origin) = &self.origin_line {
            out.push_str(origin);
        }
        match &self.residues_raw {
            Some(raw) => out.push_str(raw),
            None => out.push_str(&origin_text(&self.residues)),
        }
        out.push_str(&self.tail);
        out
    }
}

/// A GenBank file: one or more records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenBankFile {
    pub records: Vec<Record>,
    /// Anything before the first LOCUS line.
    pub preamble: String,
}

impl GenBankFile {
    pub fn to_text(&self) -> String {
        let mut out = self.preamble.clone();
        for r in &self.records {
            out.push_str(&r.to_text());
        }
        out
    }

    pub fn parse(bytes: &[u8]) -> Result<Self, GenBankError> {
        let text = std::str::from_utf8(bytes).map_err(|_| GenBankError::NotText)?;
        let lines: Vec<&str> = text.split_inclusive('\n').collect();
        let first_locus = lines
            .iter()
            .position(|l| l.starts_with("LOCUS"))
            .ok_or(GenBankError::NotGenBank)?;
        let preamble: String = lines.get(..first_locus).unwrap_or(&[]).concat();
        if !preamble.trim().is_empty() {
            return Err(GenBankError::NotGenBank);
        }
        let mut records = Vec::new();
        let mut i = first_locus;
        while i < lines.len() {
            let (record, next) = parse_record(&lines, i)?;
            records.push(record);
            i = next;
        }
        Ok(GenBankFile { records, preamble })
    }
}

fn malformed(line: usize, reason: impl Into<String>) -> GenBankError {
    GenBankError::Malformed {
        line: line + 1,
        reason: reason.into(),
    }
}

/// Parse the record starting at line `start` (a LOCUS line). Returns it and the index of the
/// next record's LOCUS line (or the end).
fn parse_record(lines: &[&str], start: usize) -> Result<(Record, usize), GenBankError> {
    let get = |i: usize| lines.get(i).copied().unwrap_or("");
    let locus_line = get(start);
    let locus =
        Locus::parse(locus_line).ok_or_else(|| malformed(start, "unreadable LOCUS line"))?;
    let mut i = start + 1;
    let mut header = String::new();
    while i < lines.len()
        && !get(i).starts_with("FEATURES")
        && !get(i).starts_with("ORIGIN")
        && !get(i).starts_with("//")
        && !get(i).starts_with("LOCUS")
    {
        header.push_str(get(i));
        i += 1;
    }
    let mut features_line = None;
    let mut features = Vec::new();
    if get(i).starts_with("FEATURES") {
        features_line = Some(get(i).to_string());
        i += 1;
        while i < lines.len() && get(i).starts_with("     ") {
            let entry_start = i;
            i += 1;
            while i < lines.len() && get(i).starts_with(&" ".repeat(QUALIFIER_INDENT)) {
                i += 1;
            }
            let raw: String = lines.get(entry_start..i).unwrap_or(&[]).concat();
            let entry = parse_entry(&raw).map_err(|e| malformed(entry_start, e))?;
            features.push(Feature {
                entry,
                raw: Some(raw),
            });
        }
    }
    let mut middle = String::new();
    while i < lines.len()
        && !get(i).starts_with("ORIGIN")
        && !get(i).starts_with("//")
        && !get(i).starts_with("LOCUS")
    {
        middle.push_str(get(i));
        i += 1;
    }
    let mut origin_line = None;
    let mut residues = Vec::new();
    let mut residues_raw = String::new();
    if get(i).starts_with("ORIGIN") {
        origin_line = Some(get(i).to_string());
        i += 1;
        while i < lines.len() && !get(i).starts_with("//") && !get(i).starts_with("LOCUS") {
            residues_raw.push_str(get(i));
            residues.extend(get(i).bytes().filter(u8::is_ascii_alphabetic));
            i += 1;
        }
    }
    let mut tail = String::new();
    while i < lines.len() && !get(i).starts_with("LOCUS") {
        tail.push_str(get(i));
        i += 1;
    }
    Ok((
        Record {
            locus,
            locus_raw: Some(locus_line.to_string()),
            header,
            features_line,
            features,
            middle,
            origin_line,
            residues,
            residues_raw: Some(residues_raw),
            tail,
        },
        i,
    ))
}

/// Parse one feature entry's text.
pub(crate) fn parse_entry(raw: &str) -> Result<FeatureEntry, String> {
    let mut lines = raw.lines();
    let first = lines.next().ok_or("empty feature")?;
    let key = first
        .get(5..QUALIFIER_INDENT)
        .unwrap_or("")
        .trim()
        .to_string();
    if key.is_empty() {
        return Err("feature without a key".into());
    }
    let mut location = first
        .get(QUALIFIER_INDENT..)
        .unwrap_or("")
        .trim()
        .to_string();
    // Continuation lines: more location until the first qualifier.
    let content: Vec<&str> = lines
        .map(|l| l.get(QUALIFIER_INDENT..).unwrap_or(""))
        .collect();
    let mut k = 0usize;
    while let Some(line) = content.get(k) {
        if line.starts_with('/') {
            break;
        }
        location.push_str(line.trim());
        k += 1;
    }
    let mut qualifiers = Vec::new();
    let mut current: Option<(String, String)> = None;
    let push = |q: Option<(String, String)>, qualifiers: &mut Vec<Qualifier>| {
        if let Some((name, value)) = q {
            qualifiers.push(qualifier(name, value));
        }
    };
    for line in content.get(k..).unwrap_or(&[]) {
        // A new qualifier starts with '/' unless we're inside an open quoted value.
        let inside_quotes = current
            .as_ref()
            .is_some_and(|(_, v)| v.starts_with('"') && !closed(v));
        if line.starts_with('/') && !inside_quotes {
            push(current.take(), &mut qualifiers);
            let body = line.get(1..).unwrap_or("");
            let (name, value) = match body.split_once('=') {
                Some((n, v)) => (n.to_string(), v.to_string()),
                None => (body.to_string(), String::new()),
            };
            current = Some((
                name,
                if body.contains('=') {
                    value
                } else {
                    "\u{0}".into()
                },
            ));
        } else if let Some((_, value)) = current.as_mut() {
            // Join continuation lines: a line that ended in a space or was cut at the width
            // continues directly; otherwise the value had a line break.
            let prev_full = value.ends_with(' ')
                || value
                    .lines()
                    .last()
                    .is_some_and(|l| l.len() + QUALIFIER_INDENT >= LINE_WIDTH - 1)
                || !value.contains(' ');
            if !prev_full {
                value.push('\n');
            }
            value.push_str(line);
        }
    }
    push(current.take(), &mut qualifiers);
    Ok(FeatureEntry {
        key,
        location,
        qualifiers,
    })
}

/// Whether a value starting with `"` has its closing quote (doubled quotes are escapes).
fn closed(v: &str) -> bool {
    let inner = v.get(1..).unwrap_or("");
    let mut quotes = 0usize;
    for c in inner.chars() {
        if c == '"' {
            quotes += 1;
        }
    }
    inner.ends_with('"') && quotes % 2 == 1
}

fn qualifier(name: String, value: String) -> Qualifier {
    if value == "\u{0}" {
        return Qualifier {
            name,
            value: None,
            quoted: false,
        };
    }
    if let Some(inner) = value.strip_prefix('"') {
        let inner = inner.strip_suffix('"').unwrap_or(inner);
        return Qualifier {
            name,
            value: Some(inner.replace("\"\"", "\"")),
            quoted: true,
        };
    }
    Qualifier {
        name,
        value: Some(value),
        quoted: false,
    }
}

/// A location, as far as we understand it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Location {
    /// `a..b`, with optional `<` / `>` partial ends.
    Span {
        start: u64,
        end: u64,
    },
    /// `a` (one base).
    Point(u64),
    /// `a^b` (between two bases).
    Between(u64, u64),
    Complement(Box<Location>),
    Join(Vec<Location>),
    Order(Vec<Location>),
    /// A location on another record (`J00194.1:100..202`), kept as text.
    Remote(String),
}

/// A flattened piece of a location: 1-based inclusive, in biological (5'→3') order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocSegment {
    pub start: u64,
    pub end: u64,
    pub reverse: bool,
}

/// Deepest nesting we accept in a location.
const MAX_DEPTH: usize = 32;

impl Location {
    pub fn parse(text: &str) -> Result<Location, String> {
        let compact: String = text.chars().filter(|c| !c.is_whitespace()).collect();
        let (loc, rest) = parse_loc(&compact, 0)?;
        if !rest.is_empty() {
            return Err(format!("unexpected {rest:?} in location"));
        }
        Ok(loc)
    }

    /// The segments in biological order, or `None` for locations that point elsewhere or are
    /// between bases.
    pub fn segments(&self) -> Option<Vec<LocSegment>> {
        let mut out = Vec::new();
        flatten(self, false, &mut out, 0)?;
        Some(out)
    }
}

fn flatten(loc: &Location, reverse: bool, out: &mut Vec<LocSegment>, depth: usize) -> Option<()> {
    if depth > MAX_DEPTH {
        return None;
    }
    match loc {
        Location::Span { start, end } => out.push(LocSegment {
            start: *start,
            end: *end,
            reverse,
        }),
        Location::Point(p) => out.push(LocSegment {
            start: *p,
            end: *p,
            reverse,
        }),
        Location::Between(..) | Location::Remote(_) => return None,
        Location::Complement(inner) => {
            let mut sub = Vec::new();
            flatten(inner, !reverse, &mut sub, depth + 1)?;
            sub.reverse();
            out.extend(sub);
        }
        Location::Join(parts) | Location::Order(parts) => {
            for p in parts {
                flatten(p, reverse, out, depth + 1)?;
            }
        }
    }
    Some(())
}

fn parse_loc(s: &str, depth: usize) -> Result<(Location, &str), String> {
    if depth > MAX_DEPTH {
        return Err("location nested too deeply".into());
    }
    for (word, kind) in [("complement(", 0), ("join(", 1), ("order(", 2)] {
        if let Some(rest) = s.strip_prefix(word) {
            let mut parts = Vec::new();
            let mut rest = rest;
            loop {
                let (loc, r) = parse_loc(rest, depth + 1)?;
                parts.push(loc);
                if let Some(r) = r.strip_prefix(',') {
                    rest = r;
                } else if let Some(r) = r.strip_prefix(')') {
                    rest = r;
                    break;
                } else {
                    return Err(format!("expected ',' or ')' in {word}...)"));
                }
            }
            let loc = match kind {
                0 => {
                    let inner = if parts.len() == 1 {
                        parts.pop().ok_or("empty complement()")?
                    } else {
                        Location::Join(parts)
                    };
                    Location::Complement(Box::new(inner))
                }
                1 => Location::Join(parts),
                _ => Location::Order(parts),
            };
            return Ok((loc, rest));
        }
    }
    // A simple location runs to the next ',' or ')'.
    let end = s.find([',', ')']).unwrap_or(s.len());
    let (simple, rest) = s.split_at(end);
    if simple.contains(':') {
        return Ok((Location::Remote(simple.to_string()), rest));
    }
    let number = |t: &str| -> Result<u64, String> {
        t.trim_start_matches(['<', '>'])
            .parse::<u64>()
            .map_err(|_| format!("bad position {t:?}"))
    };
    let loc = if let Some((a, b)) = simple.split_once("..") {
        Location::Span {
            start: number(a)?,
            end: number(b)?,
        }
    } else if let Some((a, b)) = simple.split_once('^') {
        Location::Between(number(a)?, number(b)?)
    } else if let Some((a, b)) = simple.split_once('.') {
        // `a.b`: one base somewhere in the range; treat as the range.
        Location::Span {
            start: number(a)?,
            end: number(b)?,
        }
    } else {
        Location::Point(number(simple)?)
    };
    Ok((loc, rest))
}

/// Write segments (1-based inclusive, in biological order) as a location, the way the reference
/// application does: `a..b`, `complement(a..b)`, `join(a..b,c..d)`, and for reverse-strand
/// joins `complement(join(...))` with the pieces in ascending order.
pub fn format_location(segments: &[LocSegment]) -> String {
    let piece = |s: &LocSegment| {
        if s.start == s.end {
            s.start.to_string()
        } else {
            format!("{}..{}", s.start, s.end)
        }
    };
    let all_reverse = !segments.is_empty() && segments.iter().all(|s| s.reverse);
    let all_forward = segments.iter().all(|s| !s.reverse);
    if all_reverse {
        let pieces: Vec<String> = segments.iter().rev().map(piece).collect();
        return match pieces.as_slice() {
            [one] => format!("complement({one})"),
            _ => format!("complement(join({}))", pieces.join(",")),
        };
    }
    let pieces: Vec<String> = segments
        .iter()
        .map(|s| {
            if all_forward || !s.reverse {
                piece(s)
            } else {
                format!("complement({})", piece(s))
            }
        })
        .collect();
    match pieces.as_slice() {
        [one] => one.clone(),
        _ => format!("join({})", pieces.join(",")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const RECORD: &str =
        "LOCUS       pTest                     30 bp    DNA     circular SYN 06-MAR-2024
DEFINITION  A test record.
ACCESSION   .
FEATURES             Location/Qualifiers
     source          1..30
                     /mol_type=\"other DNA\"
     CDS             complement(join(1..5,10..20))
                     /codon_start=1
                     /label=my gene
                     /note=\"E. coli lac promoter with an \"\"up\"\" mutation\"
                     /translation=\"MSHIQRETSRPRLNSNMDADLYGYKWARDNVGQSGATIYRLYGKP
                     DAPELFLKHGKGSV\"
                     /pseudo
     misc_feature    25..5
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac
//
";

    #[test]
    fn parses_records_features_and_qualifiers() {
        let f = GenBankFile::parse(RECORD.as_bytes()).unwrap();
        let r = &f.records[0];
        assert_eq!(r.locus.name, "pTest");
        assert_eq!(r.locus.length, 30);
        assert!(r.locus.circular);
        assert_eq!(r.locus.molecule, "DNA");
        assert_eq!(r.locus.division, "SYN");
        assert_eq!(r.locus.date, "06-MAR-2024");
        assert_eq!(r.residues.len(), 30);
        assert_eq!(r.features.len(), 3);
        let cds = &r.features[1].entry;
        assert_eq!(cds.key, "CDS");
        assert_eq!(cds.location, "complement(join(1..5,10..20))");
        assert_eq!(cds.qualifier("label"), Some("my gene"));
        assert_eq!(cds.qualifier("codon_start"), Some("1"));
        assert_eq!(
            cds.qualifier("note"),
            Some("E. coli lac promoter with an \"up\" mutation")
        );
        assert_eq!(
            cds.qualifier("translation"),
            Some("MSHIQRETSRPRLNSNMDADLYGYKWARDNVGQSGATIYRLYGKPDAPELFLKHGKGSV")
        );
        assert!(!cds.qualifier("translation").unwrap().contains(' '));
        assert_eq!(
            cds.qualifiers
                .last()
                .map(|q| (q.name.as_str(), q.value.is_none())),
            Some(("pseudo", true))
        );
        assert_eq!(f.to_text(), RECORD);
    }

    #[test]
    fn writers_reproduce_the_layout() {
        let f = GenBankFile::parse(RECORD.as_bytes()).unwrap();
        let r = &f.records[0];
        for feat in &r.features {
            assert_eq!(Some(feat.entry.to_text()), feat.raw, "{}", feat.entry.key);
        }
        assert_eq!(
            r.locus.to_line() + "\n",
            RECORD.lines().next().unwrap().to_string() + "\n"
        );
        assert_eq!(
            origin_text(&r.residues),
            "        1 acgtacgtac gtacgtacgt acgtacgtac\n"
        );
    }

    #[test]
    fn edits_rewrite_only_what_changed() {
        let mut f = GenBankFile::parse(RECORD.as_bytes()).unwrap();
        let r = &mut f.records[0];
        r.set_residues(b"ACGT".to_vec());
        let mut locus = r.locus.clone();
        locus.length = 4;
        r.set_locus(locus);
        let text = r.to_text();
        assert!(text.starts_with(
            "LOCUS       pTest                      4 bp    DNA     circular SYN 06-MAR-2024\n"
        ));
        assert!(text.contains("ORIGIN\n        1 acgt\n//\n"));
        assert!(text.contains("DEFINITION  A test record.\n"));
    }

    #[test]
    fn locus_edits_keep_the_original_spacing() {
        let raw =
            "LOCUS       pTEST-synthetic-vector01    6495 bp    DNA     circular SYN 20-OCT-2026\n";
        let old = Locus::parse(raw).unwrap();
        let mut new = old.clone();
        new.length = 6498;
        assert_eq!(
            patch_locus_line(raw, &old, &new).unwrap(),
            "LOCUS       pTEST-synthetic-vector01    6498 bp    DNA     circular SYN 20-OCT-2026\n"
        );
        new.length = 16498;
        assert_eq!(
            patch_locus_line(raw, &old, &new).unwrap(),
            "LOCUS       pTEST-synthetic-vector01   16498 bp    DNA     circular SYN 20-OCT-2026\n"
        );
        new.length = 649;
        new.circular = false;
        assert_eq!(
            patch_locus_line(raw, &old, &new).unwrap(),
            "LOCUS       pTEST-synthetic-vector01     649 bp    DNA     linear   SYN 20-OCT-2026\n"
        );
        new.name = "other".into();
        assert!(patch_locus_line(raw, &old, &new).is_none());
    }

    #[test]
    fn locations() {
        let seg = |start, end, reverse| LocSegment {
            start,
            end,
            reverse,
        };
        let l = Location::parse("complement(join(1926..1932,1933..1950,1951..1956))").unwrap();
        assert_eq!(
            l.segments().unwrap(),
            vec![
                seg(1951, 1956, true),
                seg(1933, 1950, true),
                seg(1926, 1932, true)
            ]
        );
        assert_eq!(
            format_location(&l.segments().unwrap()),
            "complement(join(1926..1932,1933..1950,1951..1956))"
        );
        let l = Location::parse("join(complement(20..30),complement(1..5))").unwrap();
        assert_eq!(
            l.segments().unwrap(),
            vec![seg(20, 30, true), seg(1, 5, true)]
        );
        let l = Location::parse("join(<1..5, 10..>20)").unwrap();
        assert_eq!(
            l.segments().unwrap(),
            vec![seg(1, 5, false), seg(10, 20, false)]
        );
        assert_eq!(format_location(&l.segments().unwrap()), "join(1..5,10..20)");
        assert_eq!(Location::parse("42").unwrap(), Location::Point(42));
        assert!(Location::parse("12^13").unwrap().segments().is_none());
        assert!(
            Location::parse("J00194.1:100..202")
                .unwrap()
                .segments()
                .is_none()
        );
        for bad in [
            "",
            "join(",
            "join(1..2",
            "1..x",
            "complement()",
            "1..2)",
            "join(1..2;3..4)",
        ] {
            assert!(Location::parse(bad).is_err(), "{bad:?}");
        }
        let deep = format!("{}1..2{}", "join(".repeat(100), ")".repeat(100));
        assert!(Location::parse(&deep).is_err());
    }

    #[test]
    fn rejects_non_genbank() {
        assert_eq!(
            GenBankFile::parse(b">fasta\nACGT\n"),
            Err(GenBankError::NotGenBank)
        );
        assert_eq!(GenBankFile::parse(b""), Err(GenBankError::NotGenBank));
        assert!(GenBankFile::parse(b"LOCUS\n").is_err());
    }

    #[test]
    fn long_text_wraps_after_spaces() {
        let q = Qualifier {
            name: "note".into(),
            value: Some("transcription terminator T1 from the E. coli rrnB gene".into()),
            quoted: true,
        };
        assert_eq!(
            qualifier_lines(&q),
            vec![
                "/note=\"transcription terminator T1 from the E. coli rrnB ".to_string(),
                "gene\"".to_string()
            ]
        );
    }
}
