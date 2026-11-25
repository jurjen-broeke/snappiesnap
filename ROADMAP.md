# SnappieSnap roadmap

**Stage:** M7 done (DNA, protein and multiple alignment; chromatograms: .ab1/.scf reads, Trace tab, quality trimming); M6 cloning done (restriction, Gateway, Gibson/In-Fusion, Golden Gate, mutagenesis, TA/TOPO; history read, shown, recorded and written to .dna, checked in the reference app); M4 done except methylation (deferred by owner); M5 done · **Updated:** 2026-11-13

| Parity | Value | Source |
|---|---|---|
| Version | 0.28.0 (0.x.y: features raise x, fixes raise y; see AGENTS.md › Versioning) | `Cargo.toml`, [CHANGELOG.md](CHANGELOG.md) |
| Overall parity (weighted) | 57.5% | [docs/target-app-parity.md](docs/target-app-parity.md) |
| Byte-exact round trip (open → save) | 8 / 8 local corpus files (.dna, GenBank, FASTA) | `cargo xtask corpus-test` (local corpus) |

## Milestones

| | Milestone | State |
|---|---|---|
| M0 | Workspace, never-crash lints, layering, licence/asset/wasm gates, CI | done (licence files wait on the licence decision) |
| M1 | Sequence model, command registry, undo, panic guard, **MCP**, CLI, hostile-input and scenario tests | done (control channel arrives with the app in M3) |
| M2 | File formats: **`.dna` byte-exact**, GenBank, EMBL, FASTA, GFF3 | done for `.dna`, GenBank, EMBL, GFF3, SBOL 2 and FASTA (open, save, convert); GFF3 annotations onto an open sequence |
| M3 | Desktop UI: sequence view, circular/linear maps, panels, session recording | done: shell, maps, sequence view with editing, features/history panels, control channel with `ui.render`, CLI image export, `--record` |
| M4 | Enzymes and digests, ORFs, translation | nearly done: enzymes, digests, REBASE, map labels, Enzymes view, translation (NCBI codes), ORFs, protein properties; sites in the sequence view; left: methylation sensitivity (needs a dataset decision) |
| M5 | Primers, Tm, PCR | done: binding sites, Tm, hairpins/dimers, design, PCR simulation; saving primers into `.dna` waits on a check in the reference app |
| M6 | Cloning simulation, history tree, gel simulation | done: restriction cloning, Gateway (checked against the reference app), Gibson/In-Fusion, Golden Gate and mutagenesis done; gel done; construction history read, recorded and written to `.dna` (opens in the reference app); TA/TOPO. Open: primers of earlier steps in written histories |
| M7 | Alignment, chromatograms | done: DNA and protein alignment, alignment to the construct, multiple alignment; chromatograms from .ab1 and .scf with a Trace tab, quality trimming and reads aligned to the construct. Later: traces under alignments, assembling several reads |
| M8 | Web build, performance budgets, installers | in progress: web version runs in the browser with file open and download (docs/web.md, 0.15.0); packaging on GitHub Actions (docs/packaging.md): Linux .deb/.rpm, Windows .exe/.msi, the web archive and a universal macOS .dmg, published from a version tag (unsigned until certificates are set up) |
| M9 | Parity long tail ([docs/gaps.md](docs/gaps.md)) | in progress: saving without construction history, map zoom and side-by-side maps, features table and editor, feature detection (own files and a built-in set from public records, near matches, FPbase proteins by translation) done |

## Planned

| Item | Notes |
|---|---|
| More experiment designs | The Experiments menu holds Simulate Gel, Diagnostic Digest and Check Reads. Any workflow for one construct run across several files or constructs belongs here. |
