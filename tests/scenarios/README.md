# Scenario tests

Each `*.jsonl` file here is a command script with expectations. `cargo test` replays every file in a fresh session (`crates/engine/tests/scenarios.rs`), and every step must pass. The format is described in `crates/engine/src/script.rs` and `docs/control-protocol.md`.

Turn a bug into a scenario: reproduce it through MCP, the control channel or `snappiesnap-cli run`, save the steps here, and add `expect` / `expect_error` checks that fail until the bug is fixed.

Use only synthetic sequences (or short public-domain ones) here, never corpus files.
