//! Restriction sites shown on the map: the session's enzymes (unique 6+ cutters), with
//! enzymes that recognise the same site and cut at the same place merged into one label.
//! Computed once per document change, not per frame.

use std::sync::Arc;

use snappiesnap_engine::Session;
use snappiesnap_engine::doc::{Molecule, Range, Topology};
use snappiesnap_engine::enzymes::{EnzymeDb, SiteQuery, sites};
use snappiesnap_layout::SiteMark;

/// One labelled site on the map.
#[derive(Debug, Clone, PartialEq)]
pub struct MapSite {
    /// Named in the cut list of the step being set up: drawn to stand out from the enzymes
    /// that are merely offered. Always false for the document views.
    pub chosen: bool,
    pub mark: SiteMark,
    /// The recognition site, selected when the label is clicked.
    pub range: Range,
    /// Bottom-strand cut (top-strand coordinates); the top-strand cut is `mark.pos`.
    pub cut_bottom: usize,
    /// What the hover card shows.
    pub detail: SiteDetail,
}

/// The hover card for a site: enzymes, how often they cut, overhang and the cut drawn on both
/// strands of the actual sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteDetail {
    /// Enzymes with this site and cut, joined with ", ".
    pub names: String,
    /// The one name maps and the sequence view label the site with: the name the user's own
    /// enzyme set uses, else the one most suppliers sell (see `EnzymeDb::naming_rank`), with
    /// ties broken by name so the choice is the same on every platform and every run. The
    /// other names are in `names`, which the hover card shows.
    pub label: String,
    pub notation: String,
    /// The cut as shown to the user (residues before it).
    pub position: usize,
    /// How many sites the enzyme has in the molecule.
    pub count: usize,
    /// Specific (non-N) bases in the recognition site.
    pub specific: usize,
    /// Bottom-strand cut minus top-strand cut: positive is a 5' overhang, negative 3'.
    pub overhang: i64,
    /// The residues around the site and cuts, if the window is small enough to draw.
    pub window: Option<CutWindow>,
}

impl SiteDetail {
    pub fn overhang_text(&self) -> String {
        match self.overhang {
            0 => "Blunt".to_string(),
            n if n > 0 => format!("5' overhang, {n} nt"),
            n => format!("3' overhang, {} nt", n.unsigned_abs()),
        }
    }
}

/// A stretch of both strands around a site. Indices count characters from the window start.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CutWindow {
    /// Top strand 5'→3' (a space past the end of a linear molecule).
    pub top: String,
    /// Bottom strand, complementary and aligned (so it reads 3'→5').
    pub bottom: String,
    /// The recognition site, as `start..end`.
    pub site: (usize, usize),
    /// Characters before the top- and bottom-strand cuts.
    pub cut_top: usize,
    pub cut_bottom: usize,
}

/// Widest window drawn in the hover card.
const MAX_WINDOW: i64 = 48;
/// Residues shown either side of the site and cuts.
const PAD: i64 = 2;

fn complement(b: u8) -> char {
    match b.to_ascii_uppercase() {
        b'A' => 'T',
        b'T' | b'U' => 'A',
        b'C' => 'G',
        b'G' => 'C',
        b'R' => 'Y',
        b'Y' => 'R',
        b'K' => 'M',
        b'M' => 'K',
        b'B' => 'V',
        b'V' => 'B',
        b'D' => 'H',
        b'H' => 'D',
        b' ' => ' ',
        other => other as char,
    }
}

/// The residues around a site with its cuts (top-strand coordinates).
fn cut_window(m: &Molecule, range: Range, cut: usize, cut_bottom: usize) -> Option<CutWindow> {
    let len = i64::try_from(m.len()).ok()?;
    if len == 0 {
        return None;
    }
    let circular = m.topology() == Topology::Circular;
    let start = i64::try_from(range.start).ok()?;
    let k = i64::try_from(range.len_on(m.len())).ok()?;
    // Signed distance from the site start; on a circle the shorter way round.
    let offset = |x: usize| -> Option<i64> {
        let d = i64::try_from(x).ok()? - start;
        Some(if circular && d < -len / 2 {
            d + len
        } else if circular && d > len / 2 {
            d - len
        } else {
            d
        })
    };
    let (dt, db) = (offset(cut)?, offset(cut_bottom)?);
    let lo = 0.min(dt).min(db) - PAD;
    let hi = k.max(dt).max(db) + PAD;
    if hi - lo > MAX_WINDOW {
        return None;
    }
    let seq = m.seq().as_bytes();
    let mut top = String::new();
    let mut bottom = String::new();
    for i in lo..hi {
        let pos = start + i;
        let pos = if circular { pos.rem_euclid(len) } else { pos };
        let b = usize::try_from(pos)
            .ok()
            .and_then(|p| seq.get(p))
            .map_or(b' ', u8::to_ascii_uppercase);
        top.push(b as char);
        bottom.push(complement(b));
    }
    let idx = |d: i64| usize::try_from(d - lo).ok();
    Some(CutWindow {
        top,
        bottom,
        site: (idx(0)?, idx(k)?),
        cut_top: idx(dt)?,
        cut_bottom: idx(db)?,
    })
}

/// The sites to label on `m`: the enzymes `query` picks (names the database doesn't know are
/// left out).
pub fn map_sites(db: &EnzymeDb, m: &Molecule, query: &SiteQuery) -> Vec<MapSite> {
    let mut query = query.clone();
    if let Some(names) = &mut query.enzymes {
        names.retain(|n| db.get(n).is_some());
        if names.is_empty() {
            return Vec::new();
        }
    }
    let Ok(found) = sites(db, m, &query) else {
        return Vec::new();
    };
    struct Group {
        cut: usize,
        cut_bottom: usize,
        notation: String,
        range: Range,
        names: Vec<String>,
        count: usize,
        specific: usize,
    }
    let mut groups: Vec<Group> = Vec::new();
    for e in found {
        let count = e.sites.len();
        let specific = db.get(&e.name).map_or(0, |(_, x)| x.specific_length());
        for s in e.sites {
            let (Some(cut), Some(cut_bottom)) = (s.cut, s.cut_bottom) else {
                continue;
            };
            match groups
                .iter_mut()
                .find(|g| g.cut == cut && g.notation == e.notation && g.range == s.range)
            {
                Some(g) => g.names.push(e.name.clone()),
                None => groups.push(Group {
                    cut,
                    cut_bottom,
                    notation: e.notation.clone(),
                    range: s.range,
                    names: vec![e.name.clone()],
                    count,
                    specific,
                }),
            }
        }
    }
    groups.sort_by(|a, b| a.cut.cmp(&b.cut).then_with(|| a.names.cmp(&b.names)));
    let circular = m.topology() == Topology::Circular;
    let len = m.len();
    groups
        .into_iter()
        .map(|mut g| {
            // A cut at the origin of a circular molecule is after its last residue.
            let position = if circular && g.cut == 0 { len } else { g.cut };
            // One name on the map, all of them in the hover card.
            g.names
                .sort_by(|a, b| db.naming_rank(a).cmp(&db.naming_rank(b)).then(a.cmp(b)));
            let label = g.names.first().cloned().unwrap_or_default();
            let names = g.names.join(", ");
            // Overhang in top-strand coordinates, the short way round a circle.
            let mut overhang = g.cut_bottom as i64 - g.cut as i64;
            if circular && len > 0 {
                let l = len as i64;
                if overhang > l / 2 {
                    overhang -= l;
                } else if overhang < -l / 2 {
                    overhang += l;
                }
            }
            MapSite {
                chosen: false,
                mark: SiteMark {
                    pos: g.cut,
                    text: format!("{label} ({position})"),
                    keep: false,
                },
                range: g.range,
                cut_bottom: g.cut_bottom,
                detail: SiteDetail {
                    names,
                    label,
                    notation: g.notation,
                    position,
                    count: g.count,
                    specific: g.specific,
                    overhang,
                    window: cut_window(m, g.range, g.cut, g.cut_bottom),
                },
            }
        })
        .collect()
}

/// A primer binding site on the map.
#[derive(Debug, Clone, PartialEq)]
pub struct PrimerMark {
    /// Index of the primer in the molecule's primer list.
    pub primer: usize,
    pub name: String,
    /// The annealed stretch (selected when the primer is clicked).
    pub range: Range,
    pub forward: bool,
    /// Where the primer is meant to bind, as opposed to a weaker secondary site: it pairs
    /// exactly over at least 15 bases, or along its whole length with the changes a repair or
    /// mutagenesis oligo carries. Only the first of those primes a PCR, which is why this is
    /// wider than `BindingSite::is_primary`.
    pub primary: bool,
    pub annealed: usize,
    pub mismatches: usize,
    /// Only the primer's 5' end (`annealed` bases) pairs here: it anneals but can't be
    /// extended.
    pub five_prime_only: bool,
    /// The whole primer, 5'->3'.
    pub seq: String,
    /// Stacking lane, so overlapping primers don't overlap on the map.
    pub lane: usize,
    /// The primer drawn against both strands, for the hover card.
    pub window: Option<PrimerWindow>,
    /// Tm (°C) of the annealed part (its longest perfectly paired stretch) and of the whole
    /// primer with its exact complement, at the default conditions.
    pub tm_site: Option<f64>,
    pub tm_full: Option<f64>,
    /// The primer's own structure: most stable hairpin and self-dimer ΔG (kcal/mol), and
    /// warnings where they pass the usual limits.
    pub hairpin_dg: Option<f64>,
    pub dimer_dg: Option<f64>,
    pub checks: Vec<String>,
}

/// How a primer base sits on the template.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pairing {
    Paired,
    Mismatch,
    /// Part of the 5' tail, not annealed.
    Tail,
}

/// A primer aligned to a stretch of both strands. Columns count characters from the window
/// start; the primer row is written 5'->3' for forward primers and 3'->5' (left to right) for
/// reverse ones, so every primer base sits over or under the base it pairs with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrimerWindow {
    pub top: String,
    pub bottom: String,
    /// (column, base, pairing) in column order.
    pub primer: Vec<(usize, char, Pairing)>,
    /// Primer bases left out at its 5' end because the window is too narrow.
    pub hidden_5prime: usize,
}

/// The template around a site where only the primer's 5' end pairs, with the whole primer laid
/// on it: the paired 5' part, then the 3' part overhanging where it doesn't pair.
fn five_prime_window(m: &Molecule, p: &PrimerMark) -> Option<PrimerWindow> {
    let len = i64::try_from(m.len()).ok()?;
    if len == 0 {
        return None;
    }
    let circular = m.topology() == Topology::Circular;
    let primer: Vec<char> = p.seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    let n = primer.len();
    let shown = n.min(usize::try_from(MAX_WINDOW - 2 * PAD).ok()?);
    let n_i = i64::try_from(n).ok()?;
    // Top-strand position of the primer's leftmost base: its 5' end for forward primers,
    // its 3' end for reverse ones (written 3'->5' left to right).
    let left = if p.forward {
        i64::try_from(p.range.start).ok()?
    } else {
        let end = if p.range.end == 0 {
            len
        } else {
            i64::try_from(p.range.end).ok()?
        };
        end - n_i
    };
    let skip = i64::try_from(n - shown).ok()?;
    // Forward: show the 5' end and as much of the primer as fits; reverse: the same, from the right.
    let w0 = if p.forward {
        left - PAD
    } else {
        left + skip - PAD
    };
    let width = shown + 2 * usize::try_from(PAD).ok()?;
    let seq = m.seq().as_bytes();
    let mut top = String::new();
    let mut bottom = String::new();
    for i in 0..i64::try_from(width).ok()? {
        let pos = w0 + i;
        let pos = if circular { pos.rem_euclid(len) } else { pos };
        let b = usize::try_from(pos)
            .ok()
            .and_then(|q| seq.get(q))
            .map_or(b' ', u8::to_ascii_uppercase);
        top.push(b as char);
        bottom.push(complement(b));
    }
    let pad = usize::try_from(PAD).ok()?;
    let mut row = Vec::new();
    for k in 0..shown {
        let col = pad + k;
        let t = top.chars().nth(col).filter(|c| *c != ' ');
        let (base, pairs) = if p.forward {
            let base = *primer.get(k)?;
            (base, t.is_some_and(|t| base == t || base == 'N'))
        } else {
            let i = n - 1 - (k + (n - shown));
            let base = *primer.get(i)?;
            let c = complement(base as u8);
            (base, t.is_some_and(|t| c == t || base == 'N'))
        };
        row.push((
            col,
            base,
            if pairs {
                Pairing::Paired
            } else {
                Pairing::Mismatch
            },
        ));
    }
    Some(PrimerWindow {
        top,
        bottom,
        primer: row,
        hidden_5prime: 0,
    })
}

/// The template around a primer site with the primer aligned to it.
fn primer_window(m: &Molecule, p: &PrimerMark) -> Option<PrimerWindow> {
    if p.five_prime_only {
        return five_prime_window(m, p);
    }
    let len = i64::try_from(m.len()).ok()?;
    if len == 0 {
        return None;
    }
    let circular = m.topology() == Topology::Circular;
    let primer: Vec<char> = p.seq.chars().map(|c| c.to_ascii_uppercase()).collect();
    let n = primer.len();
    let a = p.annealed.min(n);
    let budget = usize::try_from(MAX_WINDOW - 2 * PAD).ok()?;
    let shown_ann = a.min(budget);
    let tail_shown = (n - a).min(budget - shown_ann);
    let hidden_5prime = n - shown_ann - tail_shown;
    let pad = usize::try_from(PAD).ok()?;
    let width = shown_ann + tail_shown + 2 * pad;
    let start = i64::try_from(p.range.start).ok()?;
    let (a_i, shown_i, tail_i) = (
        i64::try_from(a).ok()?,
        i64::try_from(shown_ann).ok()?,
        i64::try_from(tail_shown).ok()?,
    );
    // Top-strand position of the window's first column.
    let w0 = if p.forward {
        start + (a_i - shown_i) - tail_i - PAD
    } else {
        start - PAD
    };
    let seq = m.seq().as_bytes();
    let mut top = String::new();
    let mut bottom = String::new();
    for i in 0..i64::try_from(width).ok()? {
        let pos = w0 + i;
        let pos = if circular { pos.rem_euclid(len) } else { pos };
        let b = usize::try_from(pos)
            .ok()
            .and_then(|q| seq.get(q))
            .map_or(b' ', u8::to_ascii_uppercase);
        top.push(b as char);
        bottom.push(complement(b));
    }
    let pairs_with = |base: char, template: Option<char>| {
        template.is_some_and(|t| t != ' ' && (base == t || base == 'N'))
    };
    let mut row = Vec::new();
    if p.forward {
        // 5' tail, then the annealed 3' part, ending at the window's last padded column.
        let first = n - shown_ann - tail_shown;
        for (k, &base) in primer.iter().enumerate().skip(first) {
            let col = pad + (k - first);
            let pairing = if k < n - a {
                Pairing::Tail
            } else if pairs_with(base, top.chars().nth(col)) {
                Pairing::Paired
            } else {
                Pairing::Mismatch
            };
            row.push((col, base, pairing));
        }
    } else {
        // Written 3'->5': the 3' end pairs at the site start, the tail runs to the right.
        let complement_of = |c: char| complement(c as u8);
        for k in 0..shown_ann + tail_shown {
            let Some(&base) = primer.get(n - 1 - k) else {
                continue;
            };
            let col = pad + k;
            let pairing = if k >= a {
                Pairing::Tail
            } else if pairs_with(complement_of(base), top.chars().nth(col)) {
                Pairing::Paired
            } else {
                Pairing::Mismatch
            };
            row.push((col, base, pairing));
        }
    }
    Some(PrimerWindow {
        top,
        bottom,
        primer: row,
        hidden_5prime,
    })
}

impl PrimerMark {
    /// The middle of the annealed stretch (where the label points).
    pub fn mid(&self, len: usize) -> usize {
        let mid = self.range.start.saturating_add(self.annealed / 2);
        if len == 0 { 0 } else { mid % len }
    }

    pub fn detail(&self, len: usize) -> String {
        let end = if self.range.end == 0 {
            len
        } else {
            self.range.end
        };
        let mut text = format!(
            "{}  ({} nt)\n5'-{}-3'\n{} primer, anneals at {}..{} ({} nt",
            self.name,
            self.seq.len(),
            self.seq,
            if self.forward { "Forward" } else { "Reverse" },
            self.range.start + 1,
            end,
            self.annealed
        );
        if self.mismatches > 0 {
            text.push_str(&format!(", {} mismatch", self.mismatches));
            if self.mismatches > 1 {
                text.push_str("es");
            }
        }
        text.push(')');
        if let Some(t) = self.tm_site {
            text.push_str(&format!("\nTm {t:.1} \u{b0}C"));
        }
        if self.five_prime_only {
            text.push_str("\nOnly the 5' end pairs here: the primer anneals but can't be extended");
        } else if !self.primary {
            text.push_str("\nSecondary site");
        }
        text
    }
}

/// A free energy for display: one decimal, never "-0.0"; "none" without a structure.
pub fn format_dg(v: Option<f64>) -> String {
    match v {
        Some(g) => {
            let g = (g * 10.0).round() / 10.0;
            format!("{:.1}", if g == 0.0 { 0.0 } else { g })
        }
        None => "none".to_string(),
    }
}

/// Hairpin ΔG, self-dimer ΔG and warnings for a primer (5'->3').
pub fn primer_structure(
    seq: &[u8],
    c: &snappiesnap_engine::bio::tm::TmConditions,
) -> (Option<f64>, Option<f64>, Vec<String>) {
    use snappiesnap_engine::bio::structure::{DIMER_3PRIME_WARN, dimers, hairpin, warnings};
    let hp = hairpin(seq, c).map(|h| h.dg);
    let dimer = dimers(seq, seq, c).and_then(|(any, three)| match (any, three) {
        (_, Some(t)) if t.dg <= DIMER_3PRIME_WARN => Some(t.dg),
        (Some(a), _) => Some(a.dg),
        (None, t) => t.map(|t| t.dg),
    });
    (hp, dimer, warnings(seq, c))
}

/// The primers' binding sites on `m` under the binding rules of Preferences, in stacking
/// lanes.
pub fn primer_marks(
    m: &Molecule,
    rules: &snappiesnap_engine::bio::primers::BindingPrefs,
) -> Vec<PrimerMark> {
    use snappiesnap_engine::bio::primers::{PrimerStrand, sites_with_prefs};
    use snappiesnap_engine::bio::tm::{TmConditions, tm};
    let primers = m.primers();
    let seqs: Vec<&[u8]> = primers.iter().map(|p| p.seq.as_bytes()).collect();
    let c = TmConditions::default();
    // Per primer: hairpin ΔG, self-dimer ΔG (the 3'-end one when it is the concern), warnings.
    let structure: Vec<(Option<f64>, Option<f64>, Vec<String>)> =
        seqs.iter().map(|s| primer_structure(s, &c)).collect();
    let mut marks: Vec<PrimerMark> = sites_with_prefs(m, &seqs, rules, &c)
        .into_iter()
        .filter_map(|b| {
            let p = primers.get(b.primer)?;
            Some(PrimerMark {
                primer: b.primer,
                name: p.name.clone(),
                range: b.range,
                forward: b.strand == PrimerStrand::Forward,
                primary: b.is_primary(p.seq.len())
                    || b.is_primary_with_one_mismatch(p.seq.len())
                    || b.anneals_whole_primer(p.seq.len()),
                annealed: b.annealed,
                mismatches: b.mismatches,
                five_prime_only: b.five_prime_only,
                seq: String::from_utf8_lossy(p.seq.as_bytes()).into_owned(),
                lane: 0,
                window: None,
                tm_site: b.tm(p.seq.as_bytes(), &TmConditions::default()),
                tm_full: tm(p.seq.as_bytes(), &TmConditions::default()),
                hairpin_dg: structure.get(b.primer).and_then(|s| s.0),
                dimer_dg: structure.get(b.primer).and_then(|s| s.1),
                checks: structure
                    .get(b.primer)
                    .map(|s| s.2.clone())
                    .unwrap_or_default(),
            })
        })
        .collect();
    marks.sort_by_key(|p| (p.range.start, p.range.end));
    // Greedy lanes over the unwrapped extent (a wrapping site occupies start..start+annealed).
    let mut lane_ends: Vec<usize> = Vec::new();
    for p in &mut marks {
        let (start, end) = (p.range.start, p.range.start.saturating_add(p.annealed));
        let lane = match lane_ends.iter().position(|&e| e + 2 <= start) {
            Some(l) => l,
            None => {
                lane_ends.push(0);
                lane_ends.len() - 1
            }
        };
        if let Some(e) = lane_ends.get_mut(lane) {
            *e = end;
        }
        p.lane = lane;
    }
    for p in &mut marks {
        p.window = primer_window(m, p);
    }
    marks
}

/// An open reading frame on the map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OrfMark {
    pub orf: snappiesnap_engine::bio::translate::Orf,
    /// Stacking lane, so overlapping ORFs don't overlap on the map.
    pub lane: usize,
}

/// ORFs of at least 75 amino acids (ATG starts, standard code), in stacking lanes.
pub fn orf_marks(m: &Molecule) -> Vec<OrfMark> {
    use snappiesnap_engine::bio::translate::{OrfOptions, find_orfs};
    let mut orfs = find_orfs(m, &OrfOptions::default());
    // A map with thousands of ORFs is unreadable; keep the longest.
    orfs.truncate(500);
    orfs.sort_by_key(|o| o.range.start);
    let len = m.len();
    let mut lane_ends: Vec<usize> = Vec::new();
    orfs.into_iter()
        .map(|orf| {
            let start = orf.range.start;
            let end = start + orf.range.len_on(len);
            let lane = match lane_ends.iter().position(|&e| e + 2 <= start) {
                Some(l) => l,
                None => {
                    lane_ends.push(0);
                    lane_ends.len() - 1
                }
            };
            if let Some(e) = lane_ends.get_mut(lane) {
                *e = end;
            }
            OrfMark { orf, lane }
        })
        .collect()
}

/// What the map shows besides features.
pub struct Marks<'a> {
    pub sites: &'a [MapSite],
    pub primers: &'a [PrimerMark],
    pub orfs: &'a [OrfMark],
}

/// The map's sites, primers and ORFs, recomputed only when the molecule, the enzyme
/// database or the primer binding rules changed.
#[derive(Default)]
pub struct SiteCache {
    /// Document id and session revision the marks are for.
    revision: Option<(Option<u64>, u64)>,
    query: Option<SiteQuery>,
    db: Option<Arc<EnzymeDb>>,
    topology: Option<Topology>,
    seq: Vec<u8>,
    primer_list: Vec<snappiesnap_engine::doc::Primer>,
    /// The binding rules the primer marks were found with.
    primer_rules: Option<snappiesnap_engine::bio::primers::BindingPrefs>,
    sites: Vec<MapSite>,
    primers: Vec<PrimerMark>,
    orfs: Vec<OrfMark>,
}

impl SiteCache {
    /// The marks of the open document, with the sites `query` picks.
    pub fn get(&mut self, session: &Session, query: &SiteQuery) -> Marks<'_> {
        let id = session.active_document_id();
        self.get_for(session, id, session.document(), query)
    }

    /// The marks of `doc` (the open document with tab id `id`, active or not).
    pub fn get_for(
        &mut self,
        session: &Session,
        id: Option<u64>,
        doc: Option<&snappiesnap_engine::Document>,
        query: &SiteQuery,
    ) -> Marks<'_> {
        if let Some(doc) = doc {
            let new_query = self.query.as_ref() != Some(query);
            let key = (id, session.revision());
            if self.revision != Some(key) || new_query {
                self.revision = Some(key);
                let db = session.enzymes_arc();
                let m = &doc.molecule;
                let same_molecule =
                    self.topology == Some(m.topology()) && self.seq == m.seq().as_bytes();
                let same_db = self.db.as_ref().is_some_and(|d| Arc::ptr_eq(d, &db));
                if !same_molecule || !same_db || new_query {
                    self.sites = map_sites(&db, m, query);
                    self.query = Some(query.clone());
                }
                let rules = *session.primer_binding();
                if !same_molecule
                    || self.primer_list != m.primers()
                    || self.primer_rules != Some(rules)
                {
                    self.primers = primer_marks(m, &rules);
                    self.primer_list = m.primers().to_vec();
                    self.primer_rules = Some(rules);
                }
                if !same_molecule {
                    self.orfs = orf_marks(m);
                    self.seq = m.seq().as_bytes().to_vec();
                }
                self.db = Some(db);
                self.topology = Some(m.topology());
            }
        } else {
            self.revision = None;
            self.sites.clear();
            self.primers.clear();
            self.orfs.clear();
        }
        Marks {
            sites: &self.sites,
            primers: &self.primers,
            orfs: &self.orfs,
        }
    }
}

#[cfg(test)]
mod tests {
    use snappiesnap_engine::doc::{Alphabet, Seq};

    use super::*;

    #[test]
    fn isoschizomers_share_a_label() {
        // AvaI and its isoschizomers (C^YCGRG) share one site, so they share one label: the
        // name most suppliers sell, with the others in the hover card. XmaI (C^CCGGG) cuts in
        // the same place but recognises a different site; SmaI (CCC^GGG) cuts elsewhere.
        let s = format!("{}CCCGGG{}", "A".repeat(20), "T".repeat(30));
        let m = Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Linear).unwrap();
        let sites = map_sites(&EnzymeDb::default(), &m, &SiteQuery::default());
        let texts: Vec<&str> = sites.iter().map(|s| s.mark.text.as_str()).collect();
        assert_eq!(texts, vec!["AvaI (21)", "XmaI (21)", "SmaI (23)"]);
        // The label is the name a user is most likely to have, not the first alphabetically.
        assert_eq!(sites[0].detail.label, "AvaI");
        assert!(
            sites[0].detail.names.contains("Ama87I") && sites[0].detail.names.contains("AvaI"),
            "the hover card lists them all: {}",
            sites[0].detail.names
        );
        assert_eq!(sites[1].range, Range::new(20, 26));
    }

    /// With an enzyme set of the user's own, a shared site is labelled with the name that set
    /// uses, so the map speaks the lab's vocabulary.
    #[test]
    fn a_users_own_set_names_the_site() {
        let s = format!("{}CCCGGG{}", "A".repeat(20), "T".repeat(30));
        let m = Molecule::new("m", Seq::new(Alphabet::Dna, s).unwrap(), Topology::Linear).unwrap();
        let mut db = EnzymeDb::default();
        db.sets.push(snappiesnap_engine::enzymes::UserSet {
            name: "Fridge".into(),
            enzymes: vec!["Ama87I".into()],
        });
        let sites = map_sites(&db, &m, &SiteQuery::default());
        assert_eq!(sites[0].detail.label, "Ama87I");
    }

    #[test]
    fn five_prime_only_sites_on_the_map() {
        use snappiesnap_engine::bio::primers::BindingPrefs;
        let primer = "ACGTTGCATGCAAGCTTGGCACTG";
        let comp = |c: char| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        };
        // Forward: the primer's first 18 bases at 10, then bases its last 6 don't pair with.
        let t = format!("TTTTTTTTTT{}AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", &primer[..18]);
        for reverse in [false, true] {
            let template: String = if reverse {
                t.chars().rev().map(comp).collect()
            } else {
                t.clone()
            };
            let mut m = Molecule::new(
                "m",
                Seq::new(Alphabet::Dna, template).unwrap(),
                Topology::Linear,
            )
            .unwrap();
            m.add_primer("p", Seq::new(Alphabet::Dna, primer).unwrap(), "")
                .unwrap();
            let marks = primer_marks(&m, &BindingPrefs::default());
            let mark = marks.iter().find(|p| p.five_prime_only).expect("a 5' site");
            assert_eq!(mark.forward, !reverse);
            assert!(mark.detail(m.len()).contains("can't be extended"));
            let w = mark.window.as_ref().unwrap();
            // Written 5'->3' (forward) or 3'->5' (reverse): the paired 18 are at the 5' end.
            let pairing: Vec<Pairing> = w.primer.iter().map(|(_, _, p)| *p).collect();
            let paired = pairing.iter().filter(|p| **p == Pairing::Paired).count();
            assert!(paired >= 18, "{pairing:?}");
            let five_prime_side = if reverse {
                &pairing[6..]
            } else {
                &pairing[..18]
            };
            assert!(
                five_prime_side.iter().all(|p| *p == Pairing::Paired),
                "{pairing:?}"
            );
            let off = BindingPrefs {
                five_prime_on: false,
                ..BindingPrefs::default()
            };
            assert!(!primer_marks(&m, &off).iter().any(|p| p.five_prime_only));
        }
    }

    #[test]
    fn primer_lanes_and_details() {
        let mut m = Molecule::new(
            "m",
            Seq::new(
                Alphabet::Dna,
                "AAAAAAAAAAACGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAAAAAAAAAA",
            )
            .unwrap(),
            Topology::Circular,
        )
        .unwrap();
        for (n, p) in [
            ("f", "CGTTGCATGCAAGCTT"),
            ("f2", "TTGCATGCAAGCTTGG"),
            ("r", "GTAAAACGACGGCCAGTGCCA"),
        ] {
            m.add_primer(n, Seq::new(Alphabet::Dna, p).unwrap(), "")
                .unwrap();
        }
        let marks = primer_marks(&m, &Default::default());
        let lanes: Vec<(&str, usize)> = marks.iter().map(|p| (p.name.as_str(), p.lane)).collect();
        assert_eq!(lanes, vec![("f", 0), ("f2", 1), ("r", 2)]);
        assert!(marks[0].primary && marks[0].forward && !marks[2].forward);
        assert!(
            marks[2]
                .detail(m.len())
                .contains("Reverse primer, anneals at 27..47 (21 nt)")
        );
    }

    #[test]
    fn primer_windows() {
        let template = "AAAAAAAAAAACGTTGCATGCAAGCTTGGCACTGGCCGTCGTTTTACAAAAAAAAAA";
        let mut m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, template).unwrap(),
            Topology::Linear,
        )
        .unwrap();
        // Forward with a 5' tail and a mismatch 3 bases from the 3' end; reverse plain.
        for (n, p) in [
            ("f", "gaattcCGTTGCATGCAAACTT"),
            ("r", "GTAAAACGACGGCCAGTGCCA"),
        ] {
            m.add_primer(n, Seq::new(Alphabet::Dna, p).unwrap(), "")
                .unwrap();
        }
        let marks = primer_marks(&m, &Default::default());
        let f = marks
            .iter()
            .find(|p| p.name == "f")
            .unwrap()
            .window
            .clone()
            .unwrap();
        let primer: String = f.primer.iter().map(|x| x.1).collect();
        assert_eq!(primer, "GAATTCCGTTGCATGCAAACTT");
        assert_eq!(f.primer.first().unwrap().0, 2);
        assert_eq!(&f.top[8..24], "CGTTGCATGCAAGCTT");
        let kinds: Vec<Pairing> = f.primer.iter().map(|x| x.2).collect();
        assert_eq!(kinds.iter().filter(|k| **k == Pairing::Tail).count(), 6);
        assert_eq!(kinds.iter().filter(|k| **k == Pairing::Mismatch).count(), 1);
        assert_eq!(f.hidden_5prime, 0);
        let r = marks
            .iter()
            .find(|p| p.name == "r")
            .unwrap()
            .window
            .clone()
            .unwrap();
        // Written 3'->5' under the bottom strand: it reads like the top strand's complement.
        let primer: String = r.primer.iter().map(|x| x.1).collect();
        assert_eq!(primer, "ACCGTGACCGGCAGCAAAATG");
        assert_eq!(&r.bottom[2..23], "ACCGTGACCGGCAGCAAAATG");
        assert!(r.primer.iter().all(|x| x.2 == Pairing::Paired));
    }

    #[test]
    fn free_energy_format() {
        assert_eq!(format_dg(Some(-0.04)), "0.0");
        assert_eq!(format_dg(Some(-5.26)), "-5.3");
        assert_eq!(format_dg(None), "none");
    }

    #[test]
    fn cut_windows() {
        // EcoRI G^AATTC: 5' overhang of 4, drawn with two residues either side.
        let m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, "TTTTGAATTCTTTT").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let sites = map_sites(&EnzymeDb::default(), &m, &SiteQuery::default());
        let eco = sites
            .iter()
            .find(|s| s.detail.names.contains("EcoRI"))
            .unwrap();
        assert_eq!(eco.detail.overhang_text(), "5' overhang, 4 nt");
        let w = eco.detail.window.as_ref().unwrap();
        assert_eq!(
            (w.top.as_str(), w.bottom.as_str()),
            ("TTGAATTCTT", "AACTTAAGAA")
        );
        assert_eq!((w.site, w.cut_top, w.cut_bottom), ((2, 8), 3, 7));

        // BsaI on the reverse strand (GAGACC) cuts upstream of the site, across the origin.
        let s = format!("GAGACC{}ACCCCCC", "A".repeat(20));
        let m = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, s.as_str()).unwrap(),
            Topology::Circular,
        )
        .unwrap();
        let w = cut_window(&m, Range::new(0, 6), s.len() - 5, s.len() - 1).unwrap();
        assert_eq!(w.top, "ACCCCCCGAGACCAA");
        assert_eq!((w.site, w.cut_top, w.cut_bottom), ((7, 13), 2, 6));
        // Windows at the end of a linear molecule are padded with spaces.
        let lin = Molecule::new(
            "m",
            Seq::new(Alphabet::Dna, "GAATTC").unwrap(),
            Topology::Linear,
        )
        .unwrap();
        let w = cut_window(&lin, Range::new(0, 6), 1, 5).unwrap();
        assert_eq!(w.top, "  GAATTC  ");
    }
}
