# SnappieSnap experiment file format (`.sslab`), version 1

**Status:** stable · **Updated:** 2026-11-02

A `.sslab` file holds one SnappieSnap experiment: a plan for a lab experiment across one or more DNA constructs. Version 1 has two kinds:

- a **gel experiment** (`"gel"`): an agarose gel with a size ladder in lane 1, followed by one lane per entry of `lanes`. Each lane runs a construct uncut, digested with restriction enzymes, or as the products of a PCR.
- a **reads check** (`"reads"`): one construct checked with sequencing reads (see `reads[]`).

This document is enough to read and write the format without SnappieSnap's source. Its reference implementation is `crates/engine/src/experiment.rs`, and `tests/scenarios/data/digest_check.sslab` is a complete example.

## Encoding

- One JSON object (RFC 8259), encoded as UTF-8 without a byte order mark. SnappieSnap writes it indented, with a final newline. Readers must accept any valid JSON layout.
- The file name ends in `.sslab`. SnappieSnap 0.6.0 wrote `.ssexp` for the same format; readers should accept both. There is no registered media type; use `application/json` where one is needed.
- Readers should refuse files over 256 MiB.

## Top level

| Field | Type | Required | Meaning |
|---|---|---|---|
| `format` | string | yes | Always `"snappiesnap-experiment"`. Identifies the file. |
| `version` | integer | yes | Format version, `1` for this document. |
| `generator` | string | no | What wrote the file, for example `"SnappieSnap 0.6.0"`. Informative only. |
| `kind` | string | yes | `"gel"` or `"reads"`. |
| `design` | string | no | How the lanes were chosen, if by a design. `"diagnostic_digest"`: the digests were chosen to tell the constructs apart (SnappieSnap's Diagnostic Digest). Informative: a reader that doesn't know the value shows an ordinary gel experiment. |
| `name` | string | yes | The experiment's name, as shown on its tab. |
| `gel` | object | no | The gel (see below). If absent: 1 % agarose, ladder `"1 kb"`. |
| `constructs` | array | yes | The constructs, each listed once (see below). At most 49. |
| `lanes` | array | yes | The lanes after the ladder, in gel order (see below). At most 49. Empty for `"reads"`. |
| `reads` | array | no | Kind `"reads"`: the sequencing reads (see below). At most 200. |
| `matched_reads` | boolean | no | Kind `"reads"`: also check the reads stored with the construct in its `.dna` file (SnapGene keeps reads matched to a construct inside the file). Default false. |

## `gel`

| Field | Type | Required | Meaning |
|---|---|---|---|
| `agarose` | number | yes | Agarose concentration in percent, from 0.1 to 5. |
| `ladder` | string | yes | The size ladder in lane 1, by name. SnappieSnap knows `"1 kb"`, `"100 bp"` and `"High range"`. |

## `constructs[]`

Each construct is stored twice: as a **reference** to its file, and as an **embedded copy**. The copy makes the experiment complete on its own, so it can be moved or sent to someone else. The reference lets later changes to the construct file show up in the experiment.

| Field | Type | Required | Meaning |
|---|---|---|---|
| `name` | string | yes | The construct's name. |
| `file` | string | no | The construct's file, relative to the folder holding the `.sslab` file, with `/` between folder names (for example `../plasmids/pUC19.dna`). Absent if the construct had never been saved to a file. |
| `original_path` | string | no | The construct's path as the saving application knew it, usually absolute (for example `/Users/lab/plasmids/pUC19.dna` or `C:\Lab\pUC19.dna`). Used when `file` isn't found. |
| `length` | integer | yes | Length of the sequence in bases. |
| `topology` | string | yes | `"linear"` or `"circular"`. |
| `fingerprint` | string | yes | Fingerprint of the sequence when the experiment was saved (see below). |
| `genbank` | string | yes | The construct as a complete GenBank flat file (INSDC format, `LOCUS` line to `//`), exactly as it was when the experiment was saved, including unsaved edits. |

The `file` it points to may be in any sequence format SnappieSnap reads: SnapGene `.dna`, GenBank, EMBL, FASTA, GFF3 with sequence, SBOL 2, or a sequencing trace (`.ab1`, `.scf`).

### Fingerprint

The fingerprint tells an unchanged construct file from a changed one. It is not a security measure. To compute it:

1. Take the sequence's residues as ASCII letters and convert them to upper case. Nothing else: no spaces, no numbers, no line breaks.
2. Hash the bytes with 64-bit FNV-1a: start from `h = 0xcbf29ce484222325`; for each byte `b`, `h = (h XOR b) × 0x100000001b3`, keeping the low 64 bits.
3. Write `fnv1a64:` followed by `h` as 16 lowercase hexadecimal digits, with leading zeros.

For example, the empty sequence gives `fnv1a64:cbf29ce484222325`. In Python:

```python
def fingerprint(residues: str) -> str:
    h = 0xCBF29CE484222325
    for b in residues.upper().encode("ascii"):
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return f"fnv1a64:{h:016x}"
```

## `lanes[]`

| Field | Type | Required | Meaning |
|---|---|---|---|
| `construct` | integer | yes | Index into `constructs`, from 0. Several lanes may use the same construct. |
| `method` | object or null | no | What is done with the construct. `null` or absent means the lane is empty for now. |

A `method` is an object whose `kind` says which of these it is:

| `kind` | Other fields | The lane shows |
|---|---|---|
| `"uncut"` | none | The construct itself; a circular construct runs as supercoiled DNA. |
| `"digest"` | `enzymes`: array of enzyme names, such as `["EcoRI", "BamHI"]` | The fragments after cutting with all the enzymes together. Names follow REBASE spelling. |
| `"pcr"` | `forward`, `reverse`: strings | The products of a PCR on the construct. Each primer is given by the name of one of the construct's primers, or as its sequence (5'→3'). |

## `reads[]`

Kind `"reads"` has exactly one construct (the one the reads are checked against) and a list of reads. Reads are raw data files that the user keeps, such as `.ab1` or `.scf` traces or FASTA files, so they are stored as references only, without a copy.

| Field | Type | Required | Meaning |
|---|---|---|---|
| `name` | string | yes | The read's name (its file name without the extension). |
| `file` | string | no | The read file, relative to the folder holding the `.sslab` file, with `/` between folder names. |
| `original_path` | string | no | The read's path as the saving application knew it. Used when `file` isn't found. |

When opening, a read is looked for at `file`, then at `original_path`. A read that isn't found is left out, and the user is told.

## Opening a file

SnappieSnap opens each construct as follows. Other applications should do the same, so that users see the same thing everywhere:

1. Try `file`, resolved against the folder of the `.sslab` file. Then try `original_path`.
2. If a file is found and its sequence has the stored `fingerprint`, use the file.
3. Otherwise use the embedded `genbank` copy, and tell the user why: the file wasn't found, it has changed since the experiment was saved, or it can't be read. A construct that never had a file opens from its copy without a warning.

A relative `file` must not climb above the top of the folder tree (for example `../../..` from `/a`). Treat that as a file that wasn't found.

## Compatibility

- Readers must refuse a file whose `format` isn't `"snappiesnap-experiment"`, or whose `version` is greater than the newest version they know. Each version will raise `version` for changes that older readers would misread.
- Readers must ignore fields they don't know, at every level except inside a lane's `method`. When writing a file back, SnappieSnap keeps unknown fields as they were, so data added by newer versions or other applications isn't lost.
- A `method` has exactly the fields its `kind` lists. A method with another `kind` or extra fields is an error. New kinds of method will come with a new `version`.
- A `kind` a reader doesn't know is an error for that reader. SnappieSnap 0.6 to 0.8 know only `"gel"`.
- Writers should write every required field. They must not change the meaning of the fields above.

## Example

A gel at 2 % agarose with the 100 bp ladder. pExp is found as `pExp.gb` next to the experiment; pLinear has no file and lives only in its copy. The GenBank text is shortened here; files hold the whole record.

```json
{
  "format": "snappiesnap-experiment",
  "version": 1,
  "generator": "SnappieSnap 0.6.0",
  "kind": "gel",
  "name": "Digest check",
  "gel": { "agarose": 2.0, "ladder": "100 bp" },
  "constructs": [
    {
      "name": "pExp",
      "file": "pExp.gb",
      "original_path": "/Users/lab/experiments/pExp.gb",
      "length": 60,
      "topology": "circular",
      "fingerprint": "fnv1a64:7466b1f004b64a12",
      "genbank": "LOCUS       pExp                      60 bp    DNA     circular SYN\n…\n//\n"
    },
    {
      "name": "pLinear",
      "length": 41,
      "topology": "linear",
      "fingerprint": "fnv1a64:482ac0cae691a3ed",
      "genbank": "LOCUS       pLinear                   41 bp    DNA     linear   SYN\n…\n//\n"
    }
  ],
  "lanes": [
    { "construct": 0, "method": { "kind": "uncut" } },
    { "construct": 0, "method": { "kind": "digest", "enzymes": ["EcoRI"] } },
    { "construct": 1, "method": { "kind": "digest", "enzymes": ["EcoRI"] } },
    { "construct": 1, "method": null }
  ]
}
```

## Revision history

- 2026-11-02: optional `matched_reads` for kind `"reads"` (SnappieSnap 0.11.0).
- 2026-11-02: kind `"reads"` with a `reads` list (SnappieSnap 0.9.0); still version 1.
- 2026-11-02: optional `design` field (SnappieSnap 0.8.0); still version 1, since older readers ignore it.
- 2026-11-02: the extension is `.sslab` (SnappieSnap 0.6.1); the format is unchanged.
- 2026-11-02: version 1 (SnappieSnap 0.6.0): gel experiments.
