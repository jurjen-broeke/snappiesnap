# Gaps

**Updated:** 2026-11-13

Known shortfalls, ranked. Pick work from the top. Remove or update a row when the work lands, and record the date.

| Rank | Gap | Dimension | Milestone | Notes |
|---|---|---|---|---|
| 1 | No methylation sensitivity (dam/dcm/CpG) | enzymes | M4 | Sensitivity needs a per-enzyme dataset (the source is still to be chosen) |
| 2 | Built-in DNA elements lack promoters, origins and polyA signals (their source records annotate them only as points or not at all) | features | — | 42 elements plus 976 FPbase proteins. Literature boundaries allowed; T7 promoter and the eight Gateway att sites done, more to add. Vector library download landed (system curl for HTTPS). Next: more literature-boundary elements |
| 6 | Derived caches in `.dna` (e.g. packet 0x03, `translationMW`, notes `LastModified`) aren't updated after edits | file formats | M2/M4 | Primer binding sites are rewritten since 2026-10-24 |
| 10 | No SBOL 3; a read opened from .ab1/.scf loses its chromatogram when saved as .dna; overlapping reads aren't assembled into a consensus | file formats / alignment | M2/M7 | EMBL, GFF3 and SBOL 2 done 2026-10-24. ABIF/SCF reads open, show and trim since 2026-11-02; reads matched to constructs in .dna (ZTR) are read since 0.11.0. Writing traces into .dna needs a ZTR writer and a check in the reference app |

## Revision history

- 2026-11-13: the Cloning experiment builds from the constructs open in tabs, so the web build can clone and an edited tab needs no saving first.
- 2026-11-13: the cloning experiment's enzyme fields offer what the lab holds, and say when a typed enzyme isn't in the set or in the database.
- 2026-11-02: near-match feature detection allows inserted and missing bases.
- 2026-11-02: maps exported as SVG (vector) and PNG.
- 2026-11-02: sequencing reads matched to a construct in its .dna file are read and checked.
- 2026-11-02: primer binding sites that match only at the 5' end.
- 2026-11-02: Check Reads: several reads against a construct, with the traces at each difference.
- 2026-11-02: diagnostic digest designer in the Experiments menu.
- 2026-11-02: go to a position or range by its numbers; the status bar describes the selection (GC %, Tm, molecular weight).
- 2026-11-02: gel experiments are saved to and reopened from `.ssexp` files (format documented for other applications).
- 2026-11-02: sequencing reads (.ab1, .scf) open with their chromatogram in a Trace tab; low-quality ends are trimmed (also before aligning to the construct).
- 2026-11-05: the cloning experiment: technique, vector, insert and product in four steps, from files, with the product opened as a new construct (the Cloning tab is gone).
- 2026-11-06: a folder of one-primer files imports into the oligo library, each named after its file (how labs keep cloning primers).
- 2026-11-06: Diagnostic Digest searches the chosen enzyme set (what the lab holds), so suggestions are runnable; widening is explicit and marks what the set hasn't got.
- 2026-11-05: the eight Gateway att sites detected out of the box (attB1/2, attL1/2, attR1/2, attP1/2); promoters, origins and polyA signals still open.
- 2026-10-29: lab oligo library: import a list once, find the oligos that bind a construct, add them as primers.
- 2026-10-28: feature colours set with a colour picker, and shown (stored colours, from the picker or `.dna` files, now win over the automatic colour by name and type); map/sequence background in Preferences.
- 2026-10-28: gel experiments: one gel across several open constructs, each lane with its own method.
- 2026-10-28: the map and sequence view can draw their enzymes from an enzyme set of the user's own (imported from a text file, kept between sessions).
- 2026-10-24: SBOL 2 (open, save, convert).
- 2026-10-24: EMBL and GFF3 (open, save, convert; GFF3 annotations onto an open sequence).
- 2026-10-24: user-started vector library download (HTTPS via the system curl).
- 2026-10-24: users add their own and designed features to their library.
- 2026-10-24: 14 more built-in elements (markers, regulators, CMV enhancer, ccdA/B, gIII).
- 2026-10-24: FPbase fluorescent proteins detected by translation.
- 2026-10-24: near-match feature detection (substitutions, default 2%).
- 2026-10-24: built-in features from public records (exact detection); near matches and a larger set still open.
- 2026-10-24: MCP attaches to the running app (`mcp --attach PORT`).
- 2026-10-24: undo keeps only the changed residues per step (10 Mb, 100 edits: 1.05 GB to 70 MB).
- 2026-10-24: map zoom (pinch, Cmd/Ctrl+scroll, buttons; linear maps along the sequence) and two documents side by side on the map.
- 2026-10-24: saving without construction history (a copy, or removing it from the document).
- 2026-10-24: primers and notes are edited and saved into `.dna` (binding sites follow sequence edits); confirmed in the reference application.
- 2026-10-24: the enzymes shown on the map and in the sequence view can be chosen; methylation remains.
- 2026-10-24: find in the sequence view; the copy/paste/find row is closed.
- 2026-10-24: copy, cut and paste in the sequence view; find remains.
- 2026-10-24: closed feature colours (colours must stay stable through processing rather than match the reference app; file colours are carried and our maps colour by name and type) and new-file confirmation (files built from scratch, with and without history, open in the reference app).
- 2026-10-24: save-without-history option added to the list.
- 2026-10-23: several documents open in tabs.
- 2026-10-23: sites and translations in the sequence view; methylation noted.
- 2026-10-23: `.dna` primers are imported and their binding sites computed.
- 2026-10-22: GenBank and FASTA landed; added feature colours and remaining formats.
- 2026-10-22: M3 landed (UI); reranked.
- 2026-10-22: `.dna` read/write landed; reranked file-format gaps.
- 2026-10-21: M1 done. Infrastructure and MCP gaps closed; added snapshot-undo and single-document gaps.
- 2026-10-21: page created (M0).
