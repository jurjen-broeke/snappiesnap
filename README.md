# SnappieSnap

SnappieSnap is an open-source molecular biology application written from scratch in pure Rust. It's for planning, visualising and documenting DNA, RNA and protein work: plasmid maps, features, primers, restriction and cloning simulation, alignments and chromatograms.

- **Never crashes.** Bad input produces an error message, not a crash.
- **Byte-exact `.dna` files.** Opening and saving an unchanged file gives back identical bytes.
- **Built for agents.** Every action is a command that agents can drive through MCP or a local control channel.
- **Runs everywhere.** Natively on macOS, Windows and Linux, and in the browser via WebAssembly.

> **Status:** early development. You can open, edit, save and convert SnapGene `.dna`, GenBank and FASTA files in a desktop app with map and sequence views. See [ROADMAP.md](ROADMAP.md), and [CHANGELOG.md](CHANGELOG.md) for what each version added and fixed.

SnappieSnap is an independent clean-room project. It is not affiliated with, endorsed by or derived from SnapGene or Dotmatics. "SnapGene" is used only to describe compatibility and the parity target.

## Running

```sh
cargo run --release -p snappiesnap -- path/to/plasmid.dna   # the desktop app
cargo run --release -p snappiesnap -- --sample              # a synthetic sample plasmid
cargo run --release -p snappiesnap-cli -- render plasmid.dna map.png   # an image, no window
cargo run --release -p snappiesnap-cli -- mcp               # MCP server for AI agents
(cd apps/snappiesnap-web && trunk serve --release)           # in the browser, on 127.0.0.1:8080
```

The web version is a folder of static files that any web server can host; see [docs/web.md](docs/web.md).

## Where SnappieSnap keeps its files

SnappieSnap writes its own files (settings, your feature library, downloaded libraries and
enzyme data) only to a folder in your user account. It never writes next to the application,
so it needs no administrator rights and doesn't trigger security warnings.

| | Data | Settings |
|---|---|---|
| macOS | `~/Library/Application Support/SnappieSnap` | the same folder |
| Windows | `%APPDATA%\SnappieSnap` | the same folder |
| Linux | `~/.local/share/snappiesnap` (or `$XDG_DATA_HOME`) | `~/.config/snappiesnap` (or `$XDG_CONFIG_HOME`) |

Set `SNAPPIESNAP_DATA` or `SNAPPIESNAP_SETTINGS` to an absolute path to use another location.
Relative paths are ignored.

The web version keeps nothing: open files stay in the page's memory and are saved as downloads.

## Building

Requires Rust (the version is pinned in `rust-toolchain.toml`; `rustup` installs it automatically).

```sh
cargo build --workspace
cargo xtask ci      # every quality gate: fmt, clippy, tests, layers, lints, licences, assets, wasm
```

## Contributing

Read [AGENTS.md](AGENTS.md) first. It applies to human and AI contributors alike: never crash, clean-room, pure Rust.

## Licence

SnappieSnap is licensed under either of

- the Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE)), or
- the MIT licence ([LICENSE-MIT](LICENSE-MIT)),

at your option. Unless you explicitly state otherwise, any contribution you intentionally
submit for inclusion in SnappieSnap, as defined in the Apache-2.0 licence, is dual licensed
as above, without any additional terms or conditions.

Bundled data keeps its own licence, listed in [ATTRIBUTION.md](ATTRIBUTION.md). For example,
the FPbase proteins are CC BY-SA 4.0.
