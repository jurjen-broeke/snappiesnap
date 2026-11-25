# SnappieSnap — instructions for contributors and agents

SnappieSnap is a clean-room, open-source, pure-Rust molecular biology application for planning,
visualising and documenting DNA, RNA and protein work: sequence and plasmid maps, features,
primers, restriction and cloning simulation, alignments and chromatograms. It aims at feature
parity with established tools, with speed, openness and agent control as its priorities. It runs
natively on macOS, Windows and Linux, and on the web via WASM.

The product name is always written **SnappieSnap** in user-facing text (UI, window titles, About,
installers, docs prose). Machine names stay lowercase: crates (`snappiesnap-*`), binaries, file
names, app ids. SnappieSnap is an independent project. It is not affiliated with or endorsed by
SnapGene or Dotmatics, and nothing in it may suggest otherwise.

This file applies to human and AI contributors alike. Read it first, then `ROADMAP.md` and
`docs/`.

## Where to start

1. `ROADMAP.md` is the one-page summary: stage, parity numbers, what's next.
2. `docs/target-app-parity.md` assesses parity dimension by dimension; `docs/gaps.md` ranks the
   shortfalls, so pick work from its top rows.
3. `docs/control-protocol.md` documents the command protocol that MCP and the control channel
   share.

## Never crash (outranks feature work)

Scientists trust SnappieSnap with constructs that represent weeks of bench work. A crash loses
that work, and a silently corrupted sequence is worse. Any input must produce an error that a
user or agent can act on, never a panic: a malformed GenBank file, a truncated chromatogram, an
absurd primer length, a feature spanning the origin of a circular molecule, a bad MCP parameter,
a corrupt settings file or a full disk. Fix a crash before building on top of it.

- **No panics in non-test code.** No `unwrap()`, `expect()`, `panic!`, `unreachable!`, `todo!` or
  `unimplemented!`. Return the crate's error type and propagate with `?`. Use `ok_or(..)?`,
  `let .. else`, `if let`, or `unwrap_or*` only where the fallback is genuinely correct, never
  where it would silently change a sequence. An unfinished feature returns an "unsupported"
  error. The single exception is a provably infallible literal: `#[allow(clippy::expect_used)]`
  plus `.expect("why this cannot fail")`.
- **No `unsafe`.** The workspace sets `unsafe_code = "forbid"`. If a platform integration truly
  needs it, isolate the code in its own crate, with a `SAFETY:` comment on every block.
- **Input-derived numbers are hostile.** This covers file contents, coordinates, lengths, enzyme
  sites, trace offsets and command params. Use `get()` instead of `[i]` or `[a..b]`, slice
  strings only at char boundaries, use `checked_*` or `saturating_*` for offsets and lengths,
  guard against division by zero and NaN/inf casts, and cap allocations and loop counts that
  input sizes control.
- **Circular coordinates are a classic trap.** Ranges on circular molecules wrap past the origin.
  Every range operation (insert, delete, rotate, reverse-complement, feature mapping) has tests
  for wrap-around, zero-length and full-length ranges, and for sequences of length 0 and 1.
- **Bound recursion and iteration.** Use depth limits or seen-sets for nested history, cloning
  trees and recursive file structures.
- **Don't cascade.** Recover from lock poisoning
  (`lock().unwrap_or_else(PoisonError::into_inner)`) and treat thread joins as `Result`s.
- **Last-resort guard.** The engine runs command dispatch and file import/export under
  `catch_unwind`. An escaped panic restores the document and is reported as an error. This is a
  safety net, not permission to panic. Keep `panic = "unwind"`.
- **Prove it.** Every crash fix lands with a small synthetic regression test that panicked before
  the fix. The `hostile_input` suite (`crates/engine/tests/hostile_input.rs`) runs every
  registered command in every kind of session state, with its examples mutated value by value
  into hostile values. It must stay green, and parsers get fuzz tests of their own.
- **Enforced by clippy.** The workspace `Cargo.toml` denies `unwrap_used`, `expect_used`,
  `panic`, `unimplemented`, `todo`, `unreachable` and `indexing_slicing` and forbids
  `unsafe_code`. Every crate opts in with `[lints] workspace = true`, and `cargo xtask lints`
  fails if one doesn't or if the levels are weakened. `clippy.toml` relaxes the panic lints in
  tests only.

## Scientific correctness

A wrong answer that looks right is a lab failure: a primer that won't bind, an enzyme that cuts
where the map says it doesn't.

- **Sequences are never altered silently.** Every edit goes through a command and lands in the
  undo history. Importers keep data they don't understand (unknown qualifiers, extra records) and
  report what they couldn't represent.
- **Algorithms come from cited public sources.** Each implementation has a doc comment naming its
  reference: the paper, the standard, or the public-domain dataset (e.g. SantaLucia 1998
  nearest-neighbour Tm, IUPAC ambiguity codes, the NCBI genetic code tables, Needleman–Wunsch and
  Smith–Waterman alignment). Tests check results against values worked out by hand or published
  in those sources.
- **Results are deterministic.** The same input gives the same output on every platform, natively
  and in WASM. Any randomness takes a fixed, recorded seed.

## File fidelity: `.dna` is byte-exact

`.dna` is our users' working format, so SnappieSnap must read and write it **byte for byte**.

- **Reading then writing an unmodified `.dna` file gives back identical bytes.** No reordering,
  re-encoding, whitespace normalisation or XML reformatting. `formats::dna` keeps every packet's
  raw bytes, including types we don't understand, alongside the parsed view.
- **Edits are minimal.** A change re-encodes only the packets it affects. Every other packet is
  written back verbatim and in its original order. A re-encoded packet follows the conventions
  observed in reference files (attribute order, number formatting, encodings), so files we create
  look like files the reference app creates.
- **Files we produce open in the reference app** without warnings and show the same sequence,
  features, primers and notes.
- **The fidelity tests are the gate.** For every file in the corpus: (1) read then write gives
  identical bytes; (2) read, apply one edit, write, read again gives the expected document, and
  every packet the edit didn't touch is unchanged; (3) arbitrary bytes never panic the reader
  (fuzzing). The corpus lives in `corpus/`, which is gitignored and never committed: the files in
  it are working files belonging to whoever runs the tests. Every file has its source and licence
  recorded in `corpus/MANIFEST.md`, and corpus tests must name no file and pin no sequence —
  they discover what is present and check properties that hold for any such file.
- The same approach applies to every other format as far as it can: keep what you don't
  understand, and write back what you didn't change.

## Clean-room (non-negotiable)

- **Observe behaviour, never internals.** Another application may be studied black-box only: run
  it and use its UI on synthetic sequences we made ourselves.
- **Never decompile, disassemble, debug, patch or read the contents of another application's
  bundle, binaries, scripts or resources.** File and folder names are the one exception. Don't
  use any tool that does this on our behalf. Don't read leaked or decompiled code from anywhere.
- **Never copy another product's data or content.** That includes feature and primer databases,
  enzyme sets and their curation, plasmid and vector collections, templates, help text,
  tutorials, icons, colours chosen as a brand, or any file it produced. We build our own
  databases from primary public sources and record their provenance.
- **File formats come from public documentation and our own files.** That means published
  standards (INSDC/GenBank, EMBL, FASTA, GFF3, SBOL, Applied Biosystems ABIF, SCF) and public
  third-party descriptions of the `.dna` format. Byte-level analysis is limited to `.dna` files
  we authored ourselves from synthetic sequences: comparing files before and after one known
  change. That covers data files only, never the program. Facts about a format may be learned
  from these sources, but their text and code are never copied.
- **Licences.** Dependencies must be permissively licensed (MIT, Apache-2.0, BSD, ISC, Zlib, CC0,
  or equivalent) and checked by `cargo xtask licenses`. Never copy, link or read the source of
  GPL, LGPL or AGPL bioinformatics code (for example EMBOSS or GPL'd plasmid editors). When an
  algorithm is only available in such code, implement it from the paper.

## Assets and brand (absolutely essential)

- SnappieSnap contains **no SnapGene, Dotmatics or other third-party iconography, logos, images,
  fonts, colour palettes or artwork, ever.** Screenshots of another application are never
  committed. Breaking this rule is the most serious mistake you can make in this repo.
- **Our logo and brand are original work** and live in `assets/brand/` with their own licence
  file. Don't derive them from, or make them resemble, any existing product's mark.
- **Icons are drawn in code** (`crates/ui-egui/src/icons.rs`) where practical. Any file asset
  (icon, image, font, sample sequence, enzyme or feature dataset) must be original, public
  domain/CC0, OSI-licensed or redistributable Creative Commons. It must have a row in
  `ATTRIBUTION.md` (path, title, author, source, licence) in the same change, with its licence
  file next to it. `cargo xtask assets` enforces this.
- **Sample sequences** come only from public-domain or openly licensed sources (e.g. NCBI/INSDC
  records, Addgene's openly shared sequences where the licence allows it), with the accession
  number recorded. Never use personal or unpublished lab data, in code, tests or fixtures.
- **Requests for external data** (NCBI, REBASE and so on) use a generic `SnappieSnap-dev`
  User-Agent and never contain anyone's name, email or other personal details. Check a dataset's
  licence before vendoring it.

## Architecture rules

- **Rust only.** No handwritten JavaScript or TypeScript, and no Tauri, Electron or webview
  shells. The desktop app is native (egui/eframe on wgpu). The web build is the same Rust
  compiled to WASM (trunk + wasm-bindgen). Avoid C/C++ dependencies.
- **Everything is a command.** User-visible behaviour is a command (`CommandSpec`) in
  `crates/engine/src/cmd/*`. Each command has a stable `snake_case` id, label, description, menu
  path, shortcut, typed params whose JSON Schema is generated from the Rust type, and tests. It
  also needs **at least one `.example(...)`** of valid params: examples document the command for
  agents (`commands.describe`) and seed the `hostile_input` suite, which mutates them value by
  value. A test fails if a command has no working example. Menus, toolbars, the CLI, the control
  channel and MCP all dispatch the same commands by id. Pure view state (zoom, panels, selection
  highlight) belongs to the UI crate.
- **Programmatic calls never open dialogs.** `engine.execute` runs a command's JSON form with
  defaults. Only menu or toolbar invocation opens interactive dialogs.
- **Layering** is enforced by `cargo xtask layers`. A crate depends only on lower layers, and a
  new crate must be registered there. The planned map (update it as crates land):

  ```text
  L0  seq                   alphabets, IUPAC, sequences, linear/circular ranges (no workspace deps)
  L0  formats               standalone parsers/writers: GenBank, EMBL, FASTA, GFF3, SBOL, ABIF/SCF, .dna
  L1  doc                   molecule document: sequence, features, primers, notes, history (pure data)
  L2  bio                   enzymes & digests, Tm, primers, ORFs/translation, alignment, cloning simulation
  L2  layout                map geometry: circular/linear maps, sequence view, feature lanes (no drawing)
  L3  render  io            CPU renderer (the oracle) + GPU path; document <-> file formats, native format
  L4  engine                session, command registry, undo/redo, panic guard
  L5  ui-egui  mcp          thin egui shell; MCP server and JSON control channel
  apps/ snappiesnap  snappiesnap-cli  snappiesnap-web
  ```

  Nothing below L5 depends on egui, eframe, winit, wgpu or rfd. The UI crate is swappable.
- **The UI is thin and data-driven.** UI state is serde-serialisable so the control channel can
  read and drive it. Colours, spacing and radii come from `theme::Tokens`, never hard-coded.
- **Never break WASM.** L0–L5 must pass `cargo check --target wasm32-unknown-unknown`
  (`cargo xtask wasm`). File-system and network code sits behind `cfg(not(target_arch = "wasm32"))`
  or platform services.
- **Performance is a feature.** Benchmark on large inputs in release mode: a 10 Mb chromosome, a
  1,000-enzyme digest, a 200-sequence alignment. Use rayon where it pays off, cache per document
  revision, and never rescan a whole sequence every frame.
- **Keep the native format complete.** A new `doc` field must also be added to the native format
  with `#[serde(default)]`, so older files still load. Where a field has an equivalent in an
  interchange format, map it in `io` too.

## Agent control: MCP and the control channel

Agent control is part of the product, and it is also the main way SnappieSnap gets tested. Agents
drive the real app to develop features, check their own work, reproduce bugs and write regression
tests. Keep it complete.

- **Two entry points, one protocol.** (1) `snappiesnap-cli mcp` is a headless MCP server over
  stdio that works on documents without a window. (2) `snappiesnap --control PORT` is the running
  GUI app with a JSON-lines control channel on `127.0.0.1:PORT`.
  `snappiesnap-cli mcp --attach PORT` exposes that running app over MCP. Both entry points share
  their method names and schemas, documented in `docs/control-protocol.md`.
- **Every command is reachable automatically.** MCP tools and control methods are generated from
  the command registry, with names, descriptions and JSON schemas taken from `CommandSpec`. A new
  command needs no extra wiring. A test fails if any registered command is missing from the MCP
  tool list or lacks a params schema.
- **Agents can see everything a user can.** Required introspection methods:
  - `commands.list` and `commands.describe`
  - `doc.get` (the sequence, features, primers and history as JSON, with ranges for large molecules)
  - `doc.open`, `doc.save` and `doc.export`
  - `ui.state` (open panels, selection, zoom)
  - `ui.screenshot` (the live window) and `ui.render` (offscreen, no window or focus needed)
  - `log.tail` (recent errors and warnings)
  - `app.version`
- **Never crash, never hang, never prompt.** Every method returns a structured error (`code`,
  `message`, and `hint` where there is one) for any input. Programmatic calls never open dialogs.
  Long jobs (large alignments, big digests) report progress and can be cancelled.
- **Record and replay.** A UI session can be recorded as a command script (`--record
  session.jsonl`) and replayed headless (`snappiesnap-cli run --script session.jsonl`). A
  reproduced bug becomes a scenario test under `tests/scenarios/` that CI replays, with assertions
  on `doc.get` output and, for UI work, on offscreen renders.
- **Safe by default.** The control channel is off unless `--control` is passed and binds only to
  `127.0.0.1`. File access over MCP stays inside the paths the user or agent opened. The web build
  has no control channel.

Day-to-day use:

- `cargo run --release -p snappiesnap -- --sample --control PORT` opens a sample plasmid with
  control enabled. Pick a free port.
- Send
  `{"id":1,"method":"engine.execute","params":{"command":"digest","params":{"enzymes":["EcoRI"]}}}`,
  then `{"id":2,"method":"ui.screenshot","params":{"path":"/tmp/shot.png"}}`.
- **For UI work, look at the result.** Use `ui.render` (offscreen, identical drawing code, no
  display needed) or `ui.screenshot` (the live window), and read the PNG.
- **macOS pitfall:** when the display sleeps or the screen is locked, the window isn't drawn.
  Control requests are still answered after `caffeinate -u -t 1`, but `ui.screenshot` times out.
  Use `ui.render` instead.
- Headless:
  - `snappiesnap-cli info file.dna`
  - `snappiesnap-cli render file.dna out.png [--view sequence] [--theme dark]`: the app's exact
    drawing, no window
  - `snappiesnap-cli run --script s.jsonl --export out.png`
  - `snappiesnap-cli mcp`
- **Turn manual sessions into tests:** `snappiesnap --record session.jsonl` writes every
  successful edit as a script step. Replay it with
  `snappiesnap-cli run --script session.jsonl`, add `expect` checks, and save it under
  `tests/scenarios/`.

## Before you finish a task

```sh
cargo xtask ci          # fmt, clippy -D warnings, tests (incl. hostile_input, scenarios, MCP coverage),
                        # layers, lints, licenses, assets, hygiene, wasm
cargo xtask corpus-test # if you touched formats or io: byte-exact round trips over the local corpus
cargo xtask web         # if the feature is visible in the app: dist/web + the .tar.gz for a server
```

`cargo xtask ci` is the same gate CI runs (`.github/workflows/ci.yml`). While iterating, run the
narrower `cargo test -p <crate>`, and keep `cargo xtask ci` as the final gate.

If you added or changed a command, also drive it once through MCP or the control channel and check
the result with `doc.get` or a render. Then add or update a scenario under `tests/scenarios/`.

### Committing

- **Commit every feature once it is developed and `cargo xtask ci` passes.** One feature or fix
  per commit, with a message saying what changed and why. Keep the tree green at every commit.
- **Never commit files that exist only for development and testing:**
  - `corpus/` (test files, including anything produced by another application) and `dist/`
  - `target/` and scratch or temporary files
  - local tool or editor settings

  Stage files by name (`git add <paths>`), never `git add -A` or `git add .` blindly, and never
  `git add -f`. Check `git status` before committing. `cargo xtask hygiene` fails if a local-only
  path is tracked.
- Code that tests need, such as tests, test helpers and small synthetic fixtures written by us,
  *is* committed. Test **data files** that come from outside (including anything in `corpus/`)
  are not.

### Versioning

The version is `0.x.y`, set once for the whole workspace (`[workspace.package] version` in
`Cargo.toml`). It appears in the window title, `--version` and `app.version`, so you can tell
which build is running.

- **A commit that adds a feature** (or several) raises `x` by one and sets `y` to 0: 0.3.2 becomes
  0.4.0.
- **A commit that fixes something broken or not working** raises `y` by one: 0.4.0 becomes 0.4.1.
- Commits that change only docs, tests or internal structure, with nothing a user can see, keep
  the version.
- The bump goes in the same commit as the change, with `Cargo.lock` updated (any `cargo` build
  does it), and the commit message names the new version.
- Every version gets an entry at the top of `CHANGELOG.md`, in the same commit: the version and
  date, then **Added** (features) and **Fixed** (fixes) as separate lists, written for the people
  who use the app. Leave out a list that would be empty.
- Pushing a `v<version>` tag publishes a release (`.github/workflows/release.yml`); the tag must
  match the workspace version. See `docs/packaging.md`.

## Tracking progress

When work moves a number or closes a gap, update the row, the doc's status line and its dated
revision history in the same change (`ROADMAP.md`, `docs/target-app-parity.md`, `docs/gaps.md`).

## Parallel agents

- Use your own target dir (`CARGO_TARGET_DIR=target/agent-<name>`) and edit only the crates you
  own. Make small, surgical edits to shared files (command registry list, menus), and re-read them
  just before editing.
- Put new commands in a new module (`crates/engine/src/cmd/<area>.rs`) rather than growing a
  shared file.
- If someone else's in-progress edit breaks the build, wait and retry; don't fix their files.
- Keep every `Cargo.toml` valid at all times: write manifests to a temp file outside `crates/`,
  then move them into place.
- Delete your `target/agent-*` directory when you're done.
