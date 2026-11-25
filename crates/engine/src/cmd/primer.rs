//! Primers: adding and removing them, and where they bind. Tm arrives in M5.

use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};
use snappiesnap_bio::design::{DesignOptions, design};
use snappiesnap_bio::primers::{BindingRules, BindingSite, binding_sites, sites_with_prefs};
use snappiesnap_bio::structure;
use snappiesnap_bio::tm::{TmConditions, gc_fraction, tm};
use snappiesnap_doc::{Alphabet, PrimerId, Range, Seq};

use super::residue_bytes;
use crate::{CmdError, CommandSpec, Session};

pub(super) fn specs() -> Vec<CommandSpec> {
    vec![
        CommandSpec::edit(
            "add_primer",
            "Add Primer",
            "Store a primer with the molecule. Returns the new primer's id.",
            add_primer,
        )
        .menu("Primers/Add Primer…")
        .example(json!({"name": "fwd", "residues": "ACGTACGTAC", "description": "example"})),
        CommandSpec::edit(
            "edit_primer",
            "Edit Primer",
            "Change a primer's name, sequence (5'->3') or description, by id; what isn't given stays.",
            edit_primer,
        )
        .menu("Primers/Edit Primer…")
        .example(json!({"id": 1, "name": "renamed"})),
        CommandSpec::edit(
            "remove_primer",
            "Delete Primer",
            "Remove a primer by id.",
            remove_primer,
        )
        .menu("Primers/Delete Primer")
        .example(json!({"id": 1})),
        CommandSpec::query(
            "primer_sites",
            "Primer Binding Sites",
            "Where each primer binds, under the session's binding rules (primer_binding: matched 3' bases for short and long primers, isolated mismatches, minimum Tm). Giving any of `window`, `mismatches`, `clamp` uses that rule instead: the 3'-terminal `window` bases pair with at most `mismatches` mismatches, none in the last `clamp` bases (others default to 10, 1, 3); the annealed stretch extends 5'-ward while bases pair. `primary` sites pair exactly over at least 15 bases (or the whole primer).",
            primer_sites,
        )
        .example(json!({}))
        .example(json!({"window": 15, "mismatches": 0, "clamp": 15})),
        CommandSpec::query(
            "tm",
            "Melting Temperature",
            "Tm of an oligo with its exact complement (nearest-neighbour, SantaLucia 1998; salt correction with Mg2+ as sodium equivalent, von Ahsen 2001), and its GC content. Conditions default to 100 nM oligo, 100 mM Na+, no Mg2+.",
            tm_command,
        )
        .needs_document(false)
        .example(json!({"residues": "CGTTGCATGCAAGCTT"}))
        .example(json!({"residues": "CGTTGCATGCAAGCTT", "conditions": {"oligo_nm": 250, "na_mm": 50, "mg_mm": 1.5, "dntp_mm": 0.2}})),
        CommandSpec::query(
            "design_primers",
            "Design Primers",
            "Design forward and reverse primers that amplify `range` (best pairs first). Primers start at the range ends, or up to `search` bases outside them. Scored on Tm (target 60 °C, nearest-neighbour), GC 40-60 %, 3' GC clamp, runs and repeats, 3' self- and cross-complementarity, Tm difference and other binding sites. Use add_primer to keep a pair.",
            design_primers,
        )
        .menu("Primers/Design Primers…")
        .example(json!({"range": {"start": 0, "end": 12}, "min_len": 5, "max_len": 6}))
        .example(json!({"range": {"start": 2, "end": 12}, "min_len": 5, "max_len": 6, "target_tm": 20, "search": 1, "count": 3})),
        CommandSpec::query(
            "primer_checks",
            "Check Primer Structure",
            "Hairpins and dimers of an oligo by free energy (ΔG at 37 °C, kcal/mol; nearest-neighbour stacks with salt correction, hairpin loops after SantaLucia & Hicks 2004; ungapped structures only). With `partner`, also the cross-dimer. Warnings follow common limits: hairpin -3 (-2 at the 3' end), dimer -9 (-5 at a 3' end).",
            primer_checks,
        )
        .needs_document(false)
        .example(json!({"residues": "GTAAAACGACGGCCAGTG"}))
        .example(json!({"residues": "TTTTTTGAATTC", "partner": "CCCCGAATTC"})),
    ]
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct AddParams {
    name: String,
    /// Primer sequence, 5'→3' (DNA IUPAC letters; whitespace is ignored).
    residues: String,
    #[serde(default)]
    description: String,
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct EditParams {
    id: PrimerId,
    name: Option<String>,
    /// New sequence, 5'->3' (whitespace is ignored).
    residues: Option<String>,
    description: Option<String>,
}

fn edit_primer(s: &mut Session, p: EditParams) -> Result<Value, CmdError> {
    let m = s.molecule_mut()?;
    let mut primer = m
        .primers()
        .iter()
        .find(|x| x.id == p.id)
        .cloned()
        .ok_or_else(|| {
            CmdError::new(crate::ErrorCode::NotFound, format!("no primer {}", p.id.0))
                .with_hint("primer_sites lists the primers and their ids")
        })?;
    if let Some(name) = p.name {
        primer.name = name;
    }
    if let Some(r) = p.residues {
        primer.seq = Seq::new(Alphabet::Dna, residue_bytes(&r))?;
    }
    if let Some(d) = p.description {
        primer.description = d;
    }
    m.update_primer(primer)?;
    Ok(json!({ "id": p.id }))
}

fn add_primer(s: &mut Session, p: AddParams) -> Result<Value, CmdError> {
    let seq = Seq::new(Alphabet::Dna, residue_bytes(&p.residues))?;
    let id = s.molecule_mut()?.add_primer(&p.name, seq, &p.description)?;
    Ok(json!({ "id": id }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct IdParams {
    id: PrimerId,
}

fn remove_primer(s: &mut Session, p: IdParams) -> Result<Value, CmdError> {
    s.molecule_mut()?.remove_primer(p.id)?;
    Ok(json!({}))
}

#[derive(Default, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
struct SitesParams {
    /// Conditions for the Tm values.
    conditions: Conditions,
    /// 3'-terminal bases that must pair (at least 10).
    window: Option<usize>,
    /// Mismatches allowed in the window.
    mismatches: Option<usize>,
    /// 3'-terminal bases that must pair exactly.
    clamp: Option<usize>,
}

/// Binding sites of `primers` on `m`: under explicit `window`/`mismatches`/`clamp` (missing
/// ones from [`BindingRules::default`]) when any is given, else under the session's rules.
pub(crate) fn find_sites(
    s: &Session,
    m: &snappiesnap_doc::Molecule,
    primers: &[&[u8]],
    explicit: [Option<usize>; 3],
    conditions: &TmConditions,
) -> Vec<BindingSite> {
    if explicit.iter().all(Option::is_none) {
        return sites_with_prefs(m, primers, s.primer_binding(), conditions);
    }
    let d = BindingRules::default();
    let [window, mismatches, clamp] = explicit;
    let rules = BindingRules {
        window: window.unwrap_or(d.window).min(10_000),
        mismatches: mismatches.unwrap_or(d.mismatches).min(10_000),
        clamp: clamp.unwrap_or(d.clamp).min(10_000),
        isolated: false,
    };
    binding_sites(m, primers, rules)
}

fn primer_sites(s: &Session, p: SitesParams) -> Result<Value, CmdError> {
    let m = s.molecule()?;
    let conditions = p.conditions.checked()?;
    let seqs: Vec<&[u8]> = m.primers().iter().map(|p| p.seq.as_bytes()).collect();
    let sites = find_sites(s, m, &seqs, [p.window, p.mismatches, p.clamp], &conditions);
    let primers: Vec<Value> = m
        .primers()
        .iter()
        .enumerate()
        .map(|(i, primer)| {
            let mine: Vec<Value> = sites
                .iter()
                .filter(|b| b.primer == i)
                .map(|b| {
                    json!({
                        "strand": b.strand,
                        "range": b.range,
                        "annealed": b.annealed,
                        "mismatches": b.mismatches,
                        "primary": b.is_primary(primer.seq.len()),
                        "five_prime_only": b.five_prime_only,
                        "tm": round1(b.tm(primer.seq.as_bytes(), &conditions)),
                    })
                })
                .collect();
            json!({
                "id": primer.id,
                "name": primer.name,
                "length": primer.seq.len(),
                "tm": round1(tm(primer.seq.as_bytes(), &conditions)),
                "warnings": structure::warnings(primer.seq.as_bytes(), &conditions),
                "sites": mine,
            })
        })
        .collect();
    Ok(json!({ "primers": primers }))
}

/// Reaction conditions as command params.
#[derive(Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub(crate) struct Conditions {
    /// Oligo concentration, nM (default 100).
    oligo_nm: f64,
    /// Monovalent cations, mM (default 100).
    na_mm: f64,
    /// Mg2+, mM (default 0).
    mg_mm: f64,
    /// dNTPs, mM (default 0).
    dntp_mm: f64,
}

impl Default for Conditions {
    fn default() -> Self {
        let c = TmConditions::default();
        Conditions {
            oligo_nm: c.oligo_nm,
            na_mm: c.na_mm,
            mg_mm: c.mg_mm,
            dntp_mm: c.dntp_mm,
        }
    }
}

impl Conditions {
    pub(crate) fn checked(&self) -> Result<TmConditions, CmdError> {
        let ok = |v: f64, max: f64| v.is_finite() && (0.0..=max).contains(&v);
        if !(ok(self.oligo_nm, 1e9) && self.oligo_nm > 0.0)
            || !ok(self.na_mm, 10_000.0)
            || !ok(self.mg_mm, 10_000.0)
            || !ok(self.dntp_mm, 10_000.0)
            || self.na_mm + self.mg_mm <= 0.0
        {
            return Err(CmdError::invalid_params(
                "conditions must be finite, non-negative, with some salt and oligo",
            ));
        }
        Ok(TmConditions {
            oligo_nm: self.oligo_nm,
            na_mm: self.na_mm,
            mg_mm: self.mg_mm,
            dntp_mm: self.dntp_mm,
        })
    }
}

fn round1(v: Option<f64>) -> Option<f64> {
    v.map(|t| (t * 10.0).round() / 10.0)
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct TmParams {
    /// Oligo sequence, 5'->3' (A, C, G, T; whitespace is ignored).
    residues: String,
    #[serde(default)]
    conditions: Conditions,
}

fn tm_command(_: &Session, p: TmParams) -> Result<Value, CmdError> {
    let conditions = p.conditions.checked()?;
    let bases = residue_bytes(&p.residues);
    let t = tm(&bases, &conditions).ok_or_else(|| {
        CmdError::invalid_params("Tm needs 2 to 1000 unambiguous bases (A, C, G, T)")
    })?;
    Ok(json!({
        "tm": round1(Some(t)),
        "length": bases.len(),
        "gc_percent": gc_fraction(&bases).map(|g| (g * 1000.0).round() / 10.0),
    }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct DesignParams {
    /// The region to amplify.
    range: Range,
    /// Target Tm, °C (default 60).
    target_tm: Option<f64>,
    /// Primer length limits (defaults 18 and 30).
    min_len: Option<usize>,
    max_len: Option<usize>,
    /// How far outside the region primers may start (default 0).
    search: Option<usize>,
    /// Pairs to return (default 5).
    count: Option<usize>,
    #[serde(default)]
    conditions: Conditions,
}

fn design_primers(s: &Session, p: DesignParams) -> Result<Value, CmdError> {
    let d = DesignOptions::default();
    let opts = DesignOptions {
        target_tm: p.target_tm.unwrap_or(d.target_tm),
        min_len: p.min_len.unwrap_or(d.min_len),
        max_len: p.max_len.unwrap_or(d.max_len),
        search: p.search.unwrap_or(d.search),
        count: p.count.unwrap_or(d.count),
        conditions: p.conditions.checked()?,
    };
    let m = s.molecule()?;
    let pairs = design(m, p.range, &opts).map_err(|e| {
        CmdError::new(crate::ErrorCode::InvalidInput, e.to_string())
            .with_hint("select a longer region, widen the length limits or allow a search margin")
    })?;
    let round = |v: f64| (v * 10.0).round() / 10.0;
    let primer = |d: &snappiesnap_bio::design::Designed| {
        json!({
            "residues": d.seq,
            "range": d.range,
            "length": d.seq.len(),
            "tm": round(d.tm),
            "gc_percent": round(d.gc_percent),
            "warnings": d.warnings,
        })
    };
    let pairs: Vec<Value> = pairs
        .iter()
        .map(|p| {
            json!({
                "forward": primer(&p.forward),
                "reverse": primer(&p.reverse),
                "product": p.product,
                "penalty": round(p.penalty),
                "warnings": p.warnings,
            })
        })
        .collect();
    Ok(json!({ "pairs": pairs }))
}

#[derive(Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
struct ChecksParams {
    /// The oligo, 5'->3' (A, C, G, T; whitespace is ignored).
    residues: String,
    /// Another oligo used in the same reaction, for the cross-dimer.
    #[serde(default)]
    partner: Option<String>,
    #[serde(default)]
    conditions: Conditions,
}

fn dimer_json(d: Option<structure::Dimer>) -> Value {
    d.map_or(Value::Null, |d| {
        json!({
            "dg": round1(Some(d.dg)),
            "pairs": d.pairs,
            "start": d.a_start,
            "partner_start": d.b_start,
            "three_prime": d.a_three_prime || d.b_three_prime,
        })
    })
}

fn primer_checks(_: &Session, p: ChecksParams) -> Result<Value, CmdError> {
    let conditions = p.conditions.checked()?;
    let seq = residue_bytes(&p.residues);
    let plain =
        |s: &[u8]| !s.is_empty() && s.len() <= 200 && s.iter().all(|b| b"ACGTacgt".contains(b));
    if !plain(&seq) {
        return Err(CmdError::invalid_params(
            "the oligo must be 1 to 200 bases of A, C, G, T",
        ));
    }
    let hairpin = structure::hairpin(&seq, &conditions).map_or(Value::Null, |h| {
        json!({
            "dg": round1(Some(h.dg)),
            "stem": h.stem,
            "loop": h.loop_len,
            "start": h.start,
            "end": h.end,
            "three_prime": h.three_prime,
        })
    });
    let (any, three) = structure::dimers(&seq, &seq, &conditions).unwrap_or((None, None));
    let mut warnings = structure::warnings(&seq, &conditions);
    let mut out = json!({
        "hairpin": hairpin,
        "self_dimer": dimer_json(any),
        "self_dimer_3prime": dimer_json(three),
    });
    if let Some(partner) = p.partner {
        let other = residue_bytes(&partner);
        if !plain(&other) {
            return Err(CmdError::invalid_params(
                "the partner must be 1 to 200 bases of A, C, G, T",
            ));
        }
        let (any, three) = structure::dimers(&seq, &other, &conditions).unwrap_or((None, None));
        warnings.extend(structure::pair_warnings(&seq, &other, &conditions));
        if let Some(obj) = out.as_object_mut() {
            obj.insert("cross_dimer".into(), dimer_json(any));
            obj.insert("cross_dimer_3prime".into(), dimer_json(three));
        }
    }
    if let Some(obj) = out.as_object_mut() {
        obj.insert("warnings".into(), json!(warnings));
    }
    Ok(out)
}
