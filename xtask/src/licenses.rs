//! Dependency licences (AGENTS.md › Clean-room).
//!
//! Every non-workspace package in the dependency graph must declare an SPDX licence expression
//! that is satisfiable using only [`ALLOWED`] licences.

use crate::util::{self, Result};

/// Permissive licences we accept.
const ALLOWED: &[&str] = &[
    "MIT",
    "MIT-0",
    "Apache-2.0",
    "BSD-2-Clause",
    "BSD-3-Clause",
    "0BSD",
    "ISC",
    "Zlib",
    "BSL-1.0",
    "Unlicense",
    "CC0-1.0",
    "Unicode-3.0",
    "Unicode-DFS-2016",
    "MPL-2.0",
];

/// Extra licences accepted for one specific package only, each with its reason below. The gate
/// still rejects these licences everywhere else.
const PACKAGE_EXCEPTIONS: &[(&str, &[&str])] = &[
    // egui's bundled UI fonts: font data, freely redistributable with software.
    ("epaint_default_fonts", &["OFL-1.1", "Ubuntu-font-1.0"]),
];

/// Guard against pathological nesting in a licence string.
const MAX_DEPTH: usize = 32;

pub fn run() -> Result {
    let meta = util::metadata()?;
    let members = util::workspace_members(&meta);
    let mut problems = Vec::new();
    for pkg in util::packages(&meta)
        .iter()
        .filter(|p| !members.contains(&p.id))
    {
        match pkg.license {
            None => problems.push(format!(
                "`{}` declares no SPDX licence expression; review it before depending on it",
                pkg.name
            )),
            Some(expr) => match allowed(expr, exceptions_for(pkg.name)) {
                Ok(true) => {}
                Ok(false) => problems.push(format!(
                    "`{}` has licence `{expr}`, which this project does not allow",
                    pkg.name
                )),
                Err(e) => problems.push(format!(
                    "`{}` has an unreadable licence `{expr}`: {e}",
                    pkg.name
                )),
            },
        }
    }
    util::verdict("licenses", problems)
}

fn exceptions_for(package: &str) -> &'static [&'static str] {
    PACKAGE_EXCEPTIONS
        .iter()
        .find(|(name, _)| *name == package)
        .map_or(&[], |(_, extra)| *extra)
}

/// Whether an SPDX expression (with `AND`, `OR`, `WITH`, parentheses, and the legacy `/` for
/// `OR`) can be satisfied using only allowed licences (plus `extra` for this package).
fn allowed(expr: &str, extra: &'static [&'static str]) -> Result<bool> {
    let tokens = tokenize(expr);
    let mut parser = Parser {
        tokens,
        pos: 0,
        extra,
    };
    let value = parser.or_expr(0)?;
    if parser.pos != parser.tokens.len() {
        return Err("unexpected trailing tokens".to_string());
    }
    Ok(value)
}

fn tokenize(expr: &str) -> Vec<String> {
    let spaced = expr
        .replace('(', " ( ")
        .replace(')', " ) ")
        .replace('/', " OR ");
    spaced.split_whitespace().map(str::to_string).collect()
}

struct Parser {
    tokens: Vec<String>,
    pos: usize,
    extra: &'static [&'static str],
}

impl Parser {
    fn peek(&self) -> Option<&str> {
        self.tokens.get(self.pos).map(String::as_str)
    }

    fn next(&mut self) -> Option<String> {
        let tok = self.tokens.get(self.pos).cloned();
        self.pos = self.pos.saturating_add(1);
        tok
    }

    fn or_expr(&mut self, depth: usize) -> Result<bool> {
        let mut value = self.and_expr(depth)?;
        while self.peek() == Some("OR") {
            self.next();
            let rhs = self.and_expr(depth)?;
            value = value || rhs;
        }
        Ok(value)
    }

    fn and_expr(&mut self, depth: usize) -> Result<bool> {
        let mut value = self.atom(depth)?;
        while self.peek() == Some("AND") {
            self.next();
            let rhs = self.atom(depth)?;
            value = value && rhs;
        }
        Ok(value)
    }

    fn atom(&mut self, depth: usize) -> Result<bool> {
        if depth > MAX_DEPTH {
            return Err("licence expression nested too deeply".to_string());
        }
        match self.next().as_deref() {
            Some("(") => {
                let value = self.or_expr(depth.saturating_add(1))?;
                if self.next().as_deref() != Some(")") {
                    return Err("missing `)`".to_string());
                }
                Ok(value)
            }
            Some(id) if !matches!(id, ")" | "AND" | "OR" | "WITH") => {
                let id = id.trim_end_matches('+');
                let ok = ALLOWED.contains(&id) || self.extra.contains(&id);
                if self.peek() == Some("WITH") {
                    self.next();
                    // An exception only adds permissions; it needs a name.
                    if self.next().is_none() {
                        return Err("`WITH` needs an exception name".to_string());
                    }
                }
                Ok(ok)
            }
            Some(tok) => Err(format!("unexpected `{tok}`")),
            None => Err("unexpected end of expression".to_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{allowed, exceptions_for};

    #[test]
    fn accepts_common_permissive_expressions() {
        assert!(allowed("MIT", &[]).unwrap());
        assert!(allowed("MIT OR Apache-2.0", &[]).unwrap());
        assert!(allowed("MIT/Apache-2.0", &[]).unwrap());
        assert!(allowed("(MIT OR Apache-2.0) AND Unicode-3.0", &[]).unwrap());
        assert!(allowed("Apache-2.0 WITH LLVM-exception OR MIT", &[]).unwrap());
    }

    #[test]
    fn rejects_copyleft() {
        assert!(!allowed("GPL-3.0-only", &[]).unwrap());
        assert!(!allowed("MIT AND LGPL-2.1-or-later", &[]).unwrap());
        assert!(allowed("GPL-2.0 OR MIT", &[]).unwrap());
    }

    #[test]
    fn malformed_expressions_are_errors_not_panics() {
        for bad in ["", "(", ")", "MIT AND", "OR MIT", "(MIT", "MIT WITH", "AND"] {
            assert!(allowed(bad, &[]).is_err(), "{bad:?} should be an error");
        }
        let deep = format!("{}MIT{}", "(".repeat(1000), ")".repeat(1000));
        assert!(allowed(&deep, &[]).is_err());
    }

    #[test]
    fn exceptions_apply_to_their_package_only() {
        let fonts = "(MIT OR Apache-2.0) AND OFL-1.1 AND Ubuntu-font-1.0";
        assert!(allowed(fonts, exceptions_for("epaint_default_fonts")).unwrap());
        assert!(!allowed(fonts, exceptions_for("some_other_crate")).unwrap());
        assert!(!allowed("OFL-1.1", &[]).unwrap());
    }
}
