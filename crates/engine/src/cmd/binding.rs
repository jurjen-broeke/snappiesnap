//! The rules for where primers bind, as set in Preferences: for short and for long template
//! sequences, matched 3' bases, an isolated mismatch, and a minimum Tm (`bio::primers::BindingPrefs`). They
//! are the session's, used by `primer_sites`, `oligo_library_find`, PCR's weaker sites and the
//! map; not part of any document.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::primers::{BindingPrefs, ClassRules};

use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::session(
            "primer_binding",
            "Primer Binding Rules",
            "Show or change which primer binding sites count. On a template sequence of `long_from_bp` bases or more (default 1,000,000) primers follow the `long` rules, on a shorter one the `short` rules. Each class has its own: the primer's 3'-terminal `matched` bases must pair; with `isolated_mismatches`, one mismatch whose neighbours both pair is tolerated (never in the last 3 bases); with `min_tm_on`, sites whose paired stretch melts below `min_tm` °C don't count. Primary sites (exact over 15 bases or the whole primer) always count. For both classes, with `five_prime_on` sites where the primer's 5'-terminal `five_prime_matched` bases pair but its 3' end doesn't also count (the primer anneals there but can't be extended; listed with `five_prime_only`). Only the values given change; returns all of them. Used by primer_sites, oligo_library_find, pcr with `secondary`, and the map.",
            primer_binding,
        )
        .needs_document(false)
        .example(json!({}))
        .example(json!({"long_from_bp": 50000, "short": {"matched": 12, "isolated_mismatches": false, "min_tm_on": true, "min_tm": 38}, "long": {"matched": 18, "isolated_mismatches": true, "min_tm_on": false, "min_tm": 50}}))
        .example(json!({"five_prime_on": true, "five_prime_matched": 18})),
    ]
}

#[derive(Default, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct BindingParams {
    /// Template sequences of at least this many bases are long (1,000-1,000,000,000).
    long_from_bp: Option<usize>,
    /// Rules on templates shorter than `long_from_bp`.
    short: ClassParams,
    /// Rules on templates of `long_from_bp` bases or more.
    long: ClassParams,
    /// Also count sites where only the primer's 5' end pairs.
    five_prime_on: Option<bool>,
    /// 5'-terminal bases that must pair for such a site (10-60).
    five_prime_matched: Option<usize>,
}

#[derive(Default, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct ClassParams {
    /// 3'-terminal bases that must pair (10-60).
    matched: Option<usize>,
    /// Tolerate one mismatch whose neighbours both pair.
    isolated_mismatches: Option<bool>,
    /// Leave out sites whose paired stretch melts below `min_tm`.
    min_tm_on: Option<bool>,
    /// Lowest Tm of a site, °C (0-90).
    min_tm: Option<f64>,
}

fn out_of_range(what: &str, range: &str) -> CmdError {
    CmdError::invalid_params(format!("{what} must be {range}"))
        .with_hint("call primer_binding with {} to see the current values")
}

impl ClassParams {
    /// `rules` with the values given; `class` names the class in errors.
    fn apply(&self, class: &str, mut rules: ClassRules) -> Result<ClassRules, CmdError> {
        if let Some(n) = self.matched {
            if !ClassRules::MATCH.contains(&n) {
                return Err(out_of_range(&format!("{class}.matched"), "10-60"));
            }
            rules.matched = n;
        }
        if let Some(t) = self.min_tm {
            if !ClassRules::MIN_TM.contains(&t) {
                return Err(out_of_range(&format!("{class}.min_tm"), "0-90 °C"));
            }
            rules.min_tm = t;
        }
        rules.isolated_mismatches = self
            .isolated_mismatches
            .unwrap_or(rules.isolated_mismatches);
        rules.min_tm_on = self.min_tm_on.unwrap_or(rules.min_tm_on);
        Ok(rules)
    }
}

fn primer_binding(s: &mut Session, p: BindingParams) -> Result<Value, CmdError> {
    let mut prefs = *s.primer_binding();
    if let Some(n) = p.long_from_bp {
        if !BindingPrefs::LONG_FROM.contains(&n) {
            return Err(out_of_range("long_from_bp", "1,000-1,000,000,000"));
        }
        prefs.long_from_bp = n;
    }
    if let Some(n) = p.five_prime_matched {
        if !ClassRules::MATCH.contains(&n) {
            return Err(out_of_range("five_prime_matched", "10-60"));
        }
        prefs.five_prime_matched = n;
    }
    prefs.five_prime_on = p.five_prime_on.unwrap_or(prefs.five_prime_on);
    prefs.short = p.short.apply("short", prefs.short)?;
    prefs.long = p.long.apply("long", prefs.long)?;
    s.set_primer_binding(prefs);
    serde_json::to_value(s.primer_binding())
        .map_err(|e| CmdError::invalid_params(format!("primer_binding: {e}")))
}
