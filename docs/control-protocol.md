# Control protocol

Agents, scripts and tests drive SnappieSnap through one protocol. It is available over two transports:

| Transport | Start it with | Status |
|---|---|---|
| **MCP** (Model Context Protocol, stdio) | `snappiesnap-cli mcp` | available (M1) |
| **Control channel** (JSON lines on `127.0.0.1:PORT`, inside the running app) | `snappiesnap --control PORT` | available (M3) |

Both expose the same commands and methods with the same parameter schemas.

## Conventions

- **Positions are 0-based.** A range is half-open: `{"start": 2, "end": 5}` covers positions 2, 3 and 4.
- **Ranges can wrap.** On a circular molecule, `start > end` wraps past the origin: `{"start": 8, "end": 2}` on a 10 bp plasmid covers 8, 9, 0 and 1.
- **Residues** are IUPAC letters. Case is preserved, and whitespace in parameters is ignored.
- **Failures change nothing.** A command that fails, for any reason, leaves the document exactly as it was.
- **Errors are structured:**

  ```json
  {"code": "invalid_input", "message": "position 99 is outside the sequence (length 12)", "hint": "..."}
  ```

  | Code | Meaning |
  |---|---|
  | `unknown_command` | No command with that id. |
  | `unknown_method` | No protocol method with that name. |
  | `invalid_params` | The params don't match the command's schema. |
  | `no_document` | The command needs an open document. |
  | `not_enabled` | The command can't run now, e.g. nothing to undo or unsaved changes. |
  | `invalid_input` | The edit is impossible: out of range, invalid residue, and so on. |
  | `not_found` | No item with that id. |
  | `unsupported` | Not implemented yet. |
  | `io` | A file couldn't be read or written (missing, no permission, disk full). |
  | `assertion_failed` | A script expectation didn't hold (scripts only). |
  | `internal` | A bug. The document was restored; please report it. |

## Methods

| Method | MCP tool | Params | Result |
|---|---|---|---|
| `engine.execute` | one tool per command, named after the command id | `{"command": id, "params": {...}}` | the command's result |
| `commands.list` | `commands_list` | none | `{"commands": [{id, label, description, menu, shortcut, effect, needs_document}]}` |
| `commands.describe` | `commands_describe` | `{"id": "insert"}` | `{"command": {..., "params": <JSON Schema>}}` |
| `doc.get` | `doc_get` | `{"residues": true}` | `{"open": bool, "path", "dirty", "selection", "can_undo", "can_redo", "history", "construction", "trace", "molecule"}` of the active document (`history`: the session's edits; `construction`: how the molecule was made, as in `history_tree`, or null; `trace`: for a sequencing read, its sample name, sample points and orientation, or null; its calls and signal come from `trace_info`). Several documents can be open (one per tab): the `documents` command lists them and `switch_document` changes the active one; `open`, `new` and commands that open a result add a tab. |
| `log.tail` | `log_tail` | `{"n": 50}` | `{"entries": [{seq, level, source, message}]}` |
| `app.version` | `app_version` | none | `{"name", "version", "protocol"}` |

The control channel also has methods for the running app:

| Method | Params | Result |
|---|---|---|
| `ui.state` | none | `{"state": {view, theme, map_style, folder, sidebar, features, enzymes, primers, positions, ruler, orfs, gel, experiments, active_experiment, history, translation, map_background, preferences, preferences_page, ui_scale, sequence_font, restore_window, check_updates, genetic_code, enzyme_set, map_zoom, status}, "compare", "dark", "title", "document_open"}` |
| `ui.set` | any of `{"view": "map" \| "sequence" \| "features" \| "primers" \| "enzymes" \| "gel" \| "align" \| "history", "theme": "system" \| "light" \| "dark", "cloning_sequence": true \| false (the Cloning experiment shows its construct as a sequence rather than a map; the layer flags below apply to it too), "map_style": "auto" \| "circular" \| "linear", "folder": "/path", "sidebar": true \| false, "sequence_single_line": true \| false (the sequence view runs the whole molecule along one line scrolled sideways, instead of wrapping it into rows), "align_panel": true \| false (the sequence view's Align with panel and, under the sequence, the pane aligning the ticked sequences to the document with their traces; opening it puts the sequence on one line, and turning `sequence_single_line` off closes it), "align_with": [{"document": id, "entries": [{"source": {"file": "/reads/r1.ab1"} \| {"document": tab id} \| {"matched": index}, "on": true, "collapsed": false}]}] (what each document is aligned with, in the order shown; `on` ticks it, `collapsed` folds its trace away; at most 200 entries; reads stored in a document's .dna file are listed by themselves; the alignment is `check_reads` with `rows`), "features" \| "enzymes" \| "primers" \| "positions" \| "ruler" \| "orfs": true \| false, "history": "list" \| "diagram", "translation": "orfs" \| "plus1" \| "plus2" \| "plus3" \| "minus1" \| "minus2" \| "minus3" \| "top" \| "bottom" \| "six" (what the ORF layer shows in the sequence view), "genetic_code": 1-33 (NCBI table of the frame translations), "enzyme_set": {"max_cuts": 1-3, "min_length": 0-12, "commercial_only": bool, "pool": "set name" | null, "only": ["EcoRI", ...]} (enzymes on the map and in the sequence view; `pool` draws them from one of the user's enzyme sets (`enzyme_sets_list`) instead of the whole database, and a non-empty `only` list overrides the rest), "experiments": [{"name": "Gel experiment", "kind": "gel", "step": "simulation" \| "not_started" \| "constructs" \| "enzymes" \| "lanes" \| "results", "agarose": 0.3-2, "ladder": "1 kb", "lanes": [{"document": id, "method": {"kind": "digest", "enzymes": ["EcoRI"]} \| null}], "selected": lane index \| null}] (experiment tabs, opened from the toolbar's Experiments menu beside the document tabs; kind `gel`: one gel across open documents, a lane per construct in order after the ladder; `gel_experiment` computes the same gel headless. Simulate Gel is one menu item covering kinds `gel` and `diagnostic`: it opens on step `simulation`, which chooses between them, and the steps after it follow that choice. Kind `reads` (Check Reads) has its own menu item and does not use the `simulation` step. Kind `cloning` (Experiments > Cloning) is bound to no document and has four steps of its own, `technique` > `vector` > `insert` > `product`, with `technique` (`restriction`, `gateway_lr`, `gateway_bp`, `gibson`, `infusion`, `golden_gate`, `topo_ta`, `topo_blunt`, `topo_directional`, `mutagenesis`), `vector`, `insert` and `extra` (assembly parts after the second). Each side is `{"file": "/lab/vector.dna", "document": tab id | null, "using": "whole" | "cut" | "pcr", "enzymes": "EcoRI", "forward": "", "reverse": "", "position": ""}` (a side is a file or one of the open tabs; `document` wins over `file`, which is kept beside it so a saved experiment reopens, and is the only way to choose a construct in the web build, where the browser gives the page no file paths; `document` is never written to a `.sslab` file); the step runs the same `restriction_clone`, `gateway_clone`, `assemble`, `topo_clone` or `mutagenesis` call an agent would make), "active_experiment": index (show that experiment tab; negative: the active document), "map_background": {"light": "#rrggbb" \| null, "dark": "#rrggbb" \| null} (the Map and Sequence views' background per theme, as set in Preferences; null: the theme's own), "preferences": true \| false (the Preferences window is open), "preferences_page": "general" \| "files_and_tabs" \| "export" \| "dna" \| "rna" \| "gels" \| "enzymes" \| "primers" (the category it shows), "ui_scale": 0.75-2 (interface size, 1 = normal), "sequence_font": 9-24 (size of the residues in the sequence view and the Align with pane, in points; 13 by default), "primer_tails": "simplified" \| "full" (how the sequence view shows a primer's 5' tail: a loop with its length, or written out), "restore_window": true \| false (reopen the window at its last size), "check_updates": true \| false (installed builds look for a new version once at startup), "map_zoom": 1-50, "compare": tab index (show that document beside the active one on the map; negative: off), "gel": {"mode": "digest" \| "pcr", "digest": {"agarose": 1.0, "ladder": "1 kb", "lanes": [{"kind": "uncut"}, {"kind": "digest", "enzymes": ["EcoRI"]}]}, "pcr": {"agarose": 1.5, "ladder": "100 bp", "lanes": [{"kind": "pcr", "forward": "F", "reverse": "R"}]}}}` (map layers: features, restriction sites, primer binding sites, cut positions in site labels, ruler) | as `ui.state` |
| `ui.screenshot` | `{"path": "/tmp/shot.png"}` | `{"path", "width", "height"}`: a PNG of the window as shown |
| `ui.render` | `{"path": "/tmp/render.png", "width": 1280, "height": 820, "scale": 1.0}` (size and scale optional) | `{"path", "width", "height"}`: the same UI drawn offscreen in software. It needs no display, so it works with a locked screen or a sleeping display, where `ui.screenshot` can't. |
| `ui.export_gel` | `{"path": "/tmp/gel.png", "invert": false, "experiment": 0}` (`.png` or `.jpg`; `invert`: dark bands on white; `experiment`: an experiment tab's index, else the open document's gel in the Gel tab) | `{"path", "width", "height", "inverted"}`: the gel as the app draws it, with a title and a legend of the lanes, at twice screen resolution (as **Save image…** in the Gel tab or an experiment tab, which asks whether to invert) |
| `ui.export_map` | `{"path": "/tmp/map.svg", "width": 1000, "height": 1000}` (`.svg` or `.png`; size in points, 200-4000, default square for circular maps and 1400x500 for linear ones) | `{"path", "format", "width", "height"}`: the whole map of the active document as the Map view draws it (layers, style, enzymes, background; unzoomed, no selection), as SVG with text kept as text, or as PNG at twice the size (as **Export…** in the Map view). Without a window: `snappiesnap-cli render file.dna out.svg`. |
| `ui.open_experiment` | `{"path": "/lab/run.sslab"}` | `{"experiment", "warnings"}`: opens a saved experiment ([format](experiment-format.md)) in a tab of its own, with its constructs, as **Experiments › Open Experiment…** does; `experiment` is the tab index. An experiment already open from that file is shown instead. |
| `ui.save_experiment` | `{"experiment": 0, "path": "/lab/run.sslab"}` (both optional: the shown experiment, its own file) | `{"experiment", "path", "warnings"}`: saves an experiment tab as the **Save** button does; `warnings` lists constructs kept only as copies (never saved, or with unsaved changes). Agents without a window use the `experiment_save` and `experiment_open` commands. |

Control-channel requests and responses are one JSON object per line: `{"id": 1, "method": "doc.get", "params": {}}` gets `{"id": 1, "result": ...}` or `{"id": 1, "error": {code, message, hint}}`. With `--control 0` the app picks a free port and prints it on stderr (`control channel listening on 127.0.0.1:PORT`).

`commands.list` and `commands.describe` are the source of truth for commands. Every command is listed with a JSON Schema generated from the parameter types the command actually parses.

## MCP specifics

- **Transport:** newline-delimited JSON-RPC 2.0 on stdin/stdout. Nothing else is written to stdout; diagnostics go to stderr.
- **Protocol revisions:** `2025-06-18` (preferred), `2025-03-26` and `2024-11-05`.
- **Tool results** carry the result object as `structuredContent` and as JSON text in `content`. A failed command is a tool result with `isError: true` and `structuredContent: {"error": {code, message, hint}}`, not a JSON-RPC error.
- **JSON-RPC errors** are reserved for protocol problems:
  - `-32700`: unparseable, non-UTF-8, or over 16 MiB
  - `-32600`: not a request object (batches aren't supported)
  - `-32601`: unknown method
  - `-32602`: unknown tool
- **Attached to the running app:** `snappiesnap-cli mcp --attach PORT` serves the same tools, but every call goes to the app started with `snappiesnap --control PORT` (on 127.0.0.1): commands act on its open document, and four more tools drive the window: `ui_state`, `ui_set`, `ui_render` and `ui_screenshot`. If the app goes away, calls fail with error code `io` (with a hint) instead of hanging; requests time out after 60 s.

Example session (one JSON object per line):

```json
{"jsonrpc":"2.0","id":0,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"me","version":"1"}}}
{"jsonrpc":"2.0","method":"notifications/initialized"}
{"jsonrpc":"2.0","id":1,"method":"tools/call","params":{"name":"new","arguments":{"name":"pDemo","residues":"ACGTACGT","topology":"circular"}}}
{"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":"add_feature","arguments":{"feature":{"name":"site","segments":[{"range":{"start":6,"end":2}}]}}}}
{"jsonrpc":"2.0","id":3,"method":"tools/call","params":{"name":"doc_get","arguments":{}}}
```

## Scripts and scenario tests

`snappiesnap-cli run --script <file|->` replays a JSON-lines script. The same format is used for scenario tests in `tests/scenarios/`, which `cargo test` replays.

- **One step per line:** `{"command": id, "params": {...}}` or `{"method": name, "params": {...}}`.
- **Comments:** blank lines and lines starting with `#` or `//` are skipped.
- **`"expect"`** checks the step's result by JSON pointer. Give one check or a list. Each check has a `path` and one of:
  - `eq`: the value is equal
  - `len`: the length of an array, string or object
  - `contains`: a substring, or an array element
- **`"expect_error"`** requires the step to fail with that code.
- **`"note"`** is free text.
- Any other key is an error, so typos don't pass silently.

```json
{"command": "new", "params": {"residues": "ACGTACGT", "topology": "circular"}, "expect": {"path": "/length", "eq": 8}}
{"command": "insert", "params": {"at": 99, "residues": "A"}, "expect_error": "invalid_input"}
{"method": "doc.get", "expect": [{"path": "/molecule/length", "eq": 8}, {"path": "/molecule/features", "len": 0}]}
```

Recording: `snappiesnap --record session.jsonl` writes every successful edit made in the app (by a user, the control channel or a script) as a script step, ready for `run --script`.

## Security

- The control channel is off by default and binds to `127.0.0.1` only.
- MCP runs only when a client starts `snappiesnap-cli mcp`.
- The web build has no control channel.
