# histo_contacts — Design & Implementation Plan

## 1. Purpose

`histo_contacts` computes a **bond-typed contact map** between one part of a
3D biological structure (the **query**, the "from" side) and one or more
other chains (the **target**, the "to" side) — the same query→target
contract as the sibling tool `histo_neighbours`, but using
[PDBe Arpeggio](https://github.com/PDBeurope/arpeggio) as the contact-finding
engine instead of a raw distance cutoff. Every contacting pair is classified
by Arpeggio's CREDO-derived physicochemical interaction rules (`hbond`,
`ionic`, `vdw`, `hydrophobic`, `aromatic`, `covalent`, etc.) — the **bond
types** — in addition to distance. It ships as:

1. A Python library (`import histo_contacts`)
2. A CLI tool (`histo-contacts`), built with Click, with Rich-formatted
   console output
3. A Claude Code / Claude Desktop skill that wraps the CLI

## 2. Tooling

- Python **3.14**, managed with **uv** (`uv venv`, `uv sync`, `uv run`, `uv build`)
- **Biopython** for structure parsing and query/target selector resolution —
  the same role it plays in `histo_neighbours`
- **`pdbe-arpeggio`** (imported as `arpeggio.core.InteractionComplex`) as the
  contact-classification backend. Pulls in `gemmi`, `numpy`, `biopython`.
- **`openbabel`** (Python bindings) — required by `pdbe-arpeggio` for
  hydrogen addition and bond-order perception. Not declared in
  `pdbe-arpeggio`'s own `install_requires`, so it must be declared explicitly
  as a direct dependency here. Verified installing cleanly from PyPI wheels
  (`openbabel==3.2.0`) on Python 3.12 and 3.14, macOS arm64 — no system
  package (Conda/Homebrew) required, despite Arpeggio's own README
  recommending a Conda install.
- **Click** for the CLI, **Rich** for console summary output
- `pytest` for tests

## 3. Structure loading — two parallel representations

**Biopython** (`histo_contacts.core.load_structure`) — identical to
`histo_neighbours`: picks `MMCIFParser`/`PDBParser` from the file extension
(`.cif`/`.mmcif` vs `.pdb`/`.ent`), `QUIET=True`, first model only. Used
exclusively to resolve query/target selectors into concrete chains/residues
(reuses `histo_neighbours`' selector grammar and resolution logic, reimplemented
locally per house convention).

**Arpeggio's `InteractionComplex`** needs its own, separate mmCIF file, and
has a hard requirement Biopython does not: the file must contain a populated
`_chem_comp.` category, or `InteractionComplex.__init__` raises
`ValueError: Missing _chem_comp. category in mmcif`. Verified empirically:

- A raw legacy PDB file fails even earlier — Arpeggio reads the file with
  `gemmi.cif.read(path)`, which cannot parse PDB-format text at all
  (`expected block header (data_)`).
- A coordinate-only mmCIF (e.g. `coordinates.histo.fyi`'s `_aligned.cif`
  files, which contain only an `_atom_site` loop) fails identically to a raw
  PDB — no `_chem_comp.` category present.
- A full mmCIF downloaded directly from RCSB/PDBe for the same entry works
  immediately.

Since this tool's canonical, provenance-preferred test fixture is the
coordinate-only **PDB** file from `coordinates.histo.fyi` (see §11), every
input must be normalized before reaching Arpeggio. `histo_contacts.arpeggio_backend.normalize_to_mmcif(path)`
unconditionally round-trips the input through **gemmi**:

```python
structure = gemmi.read_structure(path)
structure.setup_entities()
doc = structure.make_mmcif_document()
doc.write_file(tmp_path)
```

This always produces a `_chem_comp.` category, regardless of whether the
source file already had one, and is applied to every input uniformly (not
conditionally) so behaviour doesn't branch on source format.

Verified empirically against the 8GVI fixture:

- `auth_asym_id` and `auth_seq_id` survive the round-trip unchanged, so
  chain/residue numbers stay consistent between Biopython-resolved selectors
  and Arpeggio's output rows.
- Contact-finding results (row count, contact-type distribution) are
  statistically identical whether Arpeggio runs against this
  gemmi-normalized coordinate-only file or against a fully-populated
  deposited mmCIF for the same entry (1465 vs. 1464 total atom-atom rows for
  an identical query/target/distance, same type-count breakdown).
- The one visible side effect: when the source lacks real `_chem_comp` type
  data, gemmi's generated category has an empty `type` field per residue, so
  every residue's internal `label_comp_type` comes back `"M"` (missing)
  rather than `"P"`/`"W"` etc. This field is **not** surfaced in our output
  rows, so it's invisible to callers — documented here so it isn't
  mistaken for a bug if anyone inspects Arpeggio's raw JSON directly.

The normalized temp mmCIF is created once per `ContactMapper` (cached
alongside the Biopython parse) and reused across multiple `contact_map()`
calls — see §5 for why a *fresh* `InteractionComplex` is still built per call
even though the underlying file is reused.

## 4. Selector grammar

Identical grammar and semantics to `histo_neighbours` — reused verbatim (see
`histo_neighbours/docs/PLAN.md` §4 for the full grammar table):

- Query (`histo_contacts.selectors.parse_query`): comma-separated chains
  and/or chain:residue / chain:range tokens, unioned into one "from" atom set.
- Target (`histo_contacts.selectors.parse_target`): comma-separated list of
  chain ids only, no residue-level restriction.

### Translating query selectors into Arpeggio selection strings

Arpeggio's own selection syntax (`-s`) is `/<chain_id>/<res_num>/<atom_name>`
with fields omittable, but **does not support residue ranges** — confirmed
empirically: `/P/1-3/` raises `SelectionError`. `histo_contacts.arpeggio_backend`
translates each *resolved* `QueryPart` (already expanded against the actual
structure via the same residue-lookup logic `histo_neighbours` uses) into one
or more Arpeggio selection strings:

- A whole-chain part (`QueryPart(chain="H")`) → one string: `/H//`.
- A single-residue or range part → one `/H/<n>/` string **per residue
  actually present** in that chain within `[start, end]` (mirrors
  `histo_neighbours`' existing range-to-residue-list resolution, so gaps in
  numbering resolve correctly rather than guessing at a blind numeric range).

Selections are additive (Arpeggio's own documented behaviour), so unioning
many per-residue strings for a wide range is safe, just verbose.

## 5. Contact-finding algorithm

`ContactMapper` (in `core.py`):

1. Parses the structure once via Biopython in `__init__` (`self.structure`/
   `self.model`), exactly like `histo_neighbours`.
2. Normalizes the input to a temp mmCIF via gemmi once in `__init__`
   (`self._arpeggio_cif_path`), per §3.
3. On each `contact_map(query, target, distance=5.0)` call:
   a. Resolves `query` via Biopython into a set of **query residues**
      (`(chain_id, resnum)` keys) — residue-, not atom-, granularity, because
      Arpeggio's plane/group-level rows (see §6) represent a whole ring or
      amide group per side, not a single atom. A query token always expands
      to whole residues already (same as `histo_neighbours`), so
      residue-level membership is sufficient and exact.
   b. Resolves `target` via Biopython into a set of target chain ids
      (reuses `histo_neighbours`' chain validation — unknown chain raises
      `StructureError` listing available chains).
   c. Builds Arpeggio selection strings from the query (§4).
   d. Builds a **fresh** `InteractionComplex(self._arpeggio_cif_path, vdw_comp=0.1, interacting=distance, ph=7.4)`
      each call — `.structure_checks()`, `.initialize()`,
      `.run_arpeggio(selection_strings, distance, 0.1, False)`,
      `.get_contacts()`. A fresh instance is required per call: `run_arpeggio`
      appends to internal result lists that are only initialized once in
      `__init__`, so reusing one `InteractionComplex` across two different
      queries would silently accumulate/duplicate rows from the first query
      into the second. This is a deliberate deviation from `histo_neighbours`'
      "parse once, query many times" invariant for the Arpeggio side
      specifically (the Biopython structure and the normalized mmCIF *are*
      still reused, cheaply) — see `CLAUDE.md`.
   e. Filters Arpeggio's raw contact list: keep only rows where
      `interacting_entities == "INTER"` (Arpeggio's own flag for "one side in
      the selection, one side not" — the selection *is* our query), then
      determines which side (`bgn`/`end`) is the query side by checking
      `(auth_asym_id, auth_seq_id)` against the resolved query-residue set
      (order is not guaranteed — confirmed empirically that `bgn` is not
      always the selection side), and keeps the row only if the *other*
      side's `auth_asym_id` is in the resolved target chain set (Arpeggio's
      "not selection" side is "everything else in the structure", broader
      than our declared target chains, so this second filter is required to
      match `histo_neighbours`' target-scoping contract).
4. Arpeggio's constants are fixed, not exposed as options (matching "keep the
   CLI surface deliberately minimal" — see `CLAUDE.md` Scope): `vdw_comp=0.1`,
   `ph=7.4`, `include_sequence_adjacent=False` — Arpeggio's own CLI defaults.

**Overlap handling — a forced behavioural difference from `histo_neighbours`:**
because Arpeggio classifies every atom pair as `INTER`/`INTRA_SELECTION`/
`INTRA_NON_SELECTION` at the whole-selection level, when query and target
overlap, contacts *within* the query selection are always `INTRA_SELECTION`
and therefore always excluded — including distinct, non-identical atom pairs.
`histo_neighbours` (which does its own pairwise `NeighborSearch`) reports
those distinct pairs; this tool cannot, because Arpeggio itself never labels
them `INTER`. This is an inherent consequence of delegating contact
classification to Arpeggio, not a free design choice — documented here and in
`CLAUDE.md` rather than worked around, since working around it would mean
reimplementing Arpeggio's own pairwise geometry, defeating the point of using
it.

## 6. Interaction types covered

Arpeggio produces contacts at five granularities, all sharing the same
`bgn`/`end`/`type`/`distance`/`contact`/`interacting_entities` JSON shape
(verified directly against `arpeggio/core/interactions.py`'s `get_contacts()`
and its `_prepare_*_contact_for_export` helpers):

| `type` | Meaning | `auth_atom_id` shape |
|---|---|---|
| `atom-atom` | Single atom vs. single atom (CREDO contact types) | one atom name |
| `atom-plane` | Single atom vs. an aromatic/amide ring or group | comma-joined atom names on one side |
| `plane-plane` | Ring vs. ring (stacking geometries) | comma-joined atom names on both sides |
| `group-group` | Amide/functional group vs. group | comma-joined atom names on both sides |
| `group-plane` | Group vs. ring | comma-joined atom names on one/both sides |

Per the explicit choice made when this plan was written: **all five are
reported**, not just `atom-atom` — broader scope than `histo_neighbours`'
pure atom-pair contract, and the reason `from_atom`/`to_atom` in this tool's
output can hold a comma-joined multi-atom string rather than always a single
atom name (see `CLAUDE.md` Key Invariants). Each output row carries a
`type` field so callers can tell which shape they're looking at.

Note: the committed fixture's default-cutoff query/target combinations
exercised during development produced only `atom-atom` rows in practice
(ring-stacking and amide-group proximity are genuinely rare, precise
geometric events — not every interface has one) — the other four types are
implemented and exercised against Arpeggio's documented JSON shape, but
aren't guaranteed to appear for every query. This is normal Arpeggio
behaviour, not a gap in this tool.

## 7. Output format & sorting

Each row is a dict:

```json
{
  "type": "atom-atom",
  "from_chain": "P", "from_residue": 5, "from_atom": "CA", "from_aa": "PRO",
  "to_chain": "H", "to_residue": 63, "to_atom": "OE1", "to_aa": "GLU",
  "distance": 3.84,
  "bond_types": ["hbond", "polar"]
}
```

- `bond_types` is Arpeggio's own `contact` list (its CREDO interaction-type
  vocabulary — see `histo_neighbours`-style README table for the full set:
  `clash`, `covalent`, `vdw_clash`, `vdw`, `proximal`, `hbond`, `weak_hbond`,
  `xbond`, `ionic`, `metal_complex`, `aromatic`, `hydrophobic`, `carbonyl`,
  `polar`, `weak_polar`, plus ring/group-level types like `FT`/`ET`/
  `AMIDEAMIDE`/`AMIDERING` for the non-atom-atom rows), sorted alphabetically
  for deterministic output.
- `distance` is Arpeggio's own computed value, used as-is (confirmed
  Arpeggio rounds to 2dp internally) rather than recomputed independently —
  unlike `histo_neighbours`, which computes its own Euclidean distance,
  because plane/group-level distances aren't a simple atom-atom Euclidean
  distance we could recompute ourselves.
- The full contact map is a **JSON list** of these row dicts, no wrapping
  metadata object — same contract as `histo_neighbours`.

Sort order: same convention as `histo_neighbours` — `(to_chain, to_residue,
to_atom, from_chain, from_residue, from_atom)` ascending. Only
`(to_chain, to_residue)` is part of the requested contract; the rest exists
purely for deterministic output.

## 8. Library API (`histo_contacts/core.py`)

```python
from histo_contacts import ContactMapper

cm = ContactMapper("structure.pdb")           # parses once, normalizes once
rows = cm.contact_map("P", ["H", "L"], distance=5.0)   # -> list[dict]
cm.write_json(rows, "contacts.json")
```

Module-level convenience function `contact_map(path, query, target,
distance=5.0)` wraps `ContactMapper` for one-shot scripting use — same shape
as `histo_neighbours`.

## 9. CLI (`histo_contacts/cli.py`)

```
histo-contacts FILENAME --query SPEC --target SPEC --output PATH [--distance FLOAT]
```

Same four-plus-one option surface as `histo_neighbours` — `FILENAME`,
`--query`/`-q`, `--target`/`-t`, `--output`/`-o`, `--distance`/`-d` (default
`5.0`, forwarded to Arpeggio's `interacting` cutoff). No Arpeggio-specific
knobs (`vdw-comp`, `ph`, `minimise-hydrogens`, `include-sequence-adjacent`)
are exposed — fixed per §5.

After computing, the CLI prints a Rich table (truncated with a "N more rows"
note beyond a display cap) with a `bond_types` column (joined with `, `)
alongside the same columns `histo_neighbours` shows, plus a one-line summary,
then confirms the file written.

## 10. Claude skill

`skills/histo-contacts/SKILL.md` — describes when/how to invoke the
`histo-contacts` CLI (computing a bond-typed contact map between a chain/
residue selection and one or more other chains), cross-referencing
`histo-neighbours`' skill as the plain-distance sibling.

## 11. Test fixture

Per the explicit choice made when this plan was written: the same PDB entry
`histo_neighbours` uses (8GVI, a pMHC-TCR complex), but the **coordinate-only
PDB** download (not mmCIF) from its original documented source
(`coordinates.histo.fyi`), matching `histo_neighbours`' data-provenance
convention:

```
https://coordinates.histo.fyi/structures/downloads/class_i/without_solvent/8gvi_1_aligned.pdb
```

Chains: `H` (1-274, MHC heavy chain), `L` (1-99, β2m), `P` (1-8, peptide),
`A` (3-205), `B` (3-245) (TCR α/β) — identical numbering to
`histo_neighbours`' fixture (it's the same structure, just legacy-PDB format
instead of coordinate-only mmCIF). Verified this file requires the gemmi
normalization step of §3 to be usable with Arpeggio at all (a raw feed to
`pdbe-arpeggio` fails); confirmed the normalized version produces
contact-finding results equivalent to running Arpeggio against a fully
deposited mmCIF for the same entry.

## 12. Package layout

```
histo_contacts/
  pyproject.toml
  README.md
  CLAUDE.md
  CHANGELOG.md
  docs/PLAN.md
  .gitignore
  src/histo_contacts/
    __init__.py
    core.py              # load_structure(), ContactMapper, contact_map()
    selectors.py         # query/target selector parsing (same grammar as histo_neighbours)
    arpeggio_backend.py  # gemmi normalization, selector translation, Arpeggio invocation, row mapping
    cli.py                # Click CLI (entry point: histo-contacts)
    py.typed
  skills/histo-contacts/SKILL.md
  tests/
    fixtures/
      8gvi_1_aligned.pdb   # pMHC-TCR complex, coordinate-only PDB from coordinates.histo.fyi
    test_core.py
    test_selectors.py
    test_arpeggio_backend.py
    test_cli.py
  tmp/
    .gitkeep
```

## 13. Testing plan

Unit tests with pytest against the committed `8gvi_1_aligned.pdb` fixture —
real Arpeggio runs throughout, no mocking of Arpeggio's output:

- `test_selectors.py`: same coverage as `histo_neighbours` (same grammar).
- `test_arpeggio_backend.py`:
  - `normalize_to_mmcif()` produces a file with a `_chem_comp.` category, and
    preserves `auth_asym_id`/`auth_seq_id`/coordinates for a known atom.
  - Selector-string translation: whole chain → `["/H//"]`; single residue →
    `["/P/5/"]`; range → one string per residue actually present.
- `test_core.py` (`ContactMapper`):
  - Query a whole chain (`"P"`) against a single target chain (`["H"]`) —
    non-empty result, correct row shape/keys, `from_chain == "P"`,
    `to_chain == "H"`, `0 < distance <= cutoff`, `bond_types` non-empty.
  - Multi-target sort order (`to_chain`, then `to_residue`), matching
    `histo_neighbours`.
  - Single-residue and residue-range queries — only the requested residues
    appear as `from_residue`.
  - Every row's `type` is one of the five known Arpeggio interaction types.
  - Distance cutoff: a larger cutoff strictly increases (or maintains) the
    contact count relative to a smaller one.
  - Error paths: unknown target chain, unknown query chain, bad selector
    syntax, residue range selector (rejected by Arpeggio's own selection
    parser, surfaced as a clear error rather than a raw traceback).
  - Two sequential `contact_map()` calls on the same `ContactMapper` with
    different queries return independent, non-accumulating results (guards
    the fresh-`InteractionComplex`-per-call design in §5).
- `test_cli.py`: writes valid JSON to `--output`, matches the library result
  for the same arguments, non-zero exit + usage error when a required option
  is missing.

## 14. Workflow

1. Write this plan, commit it.
2. Scaffold project with `uv init`/`uv add`.
3. Implement `selectors.py`, `arpeggio_backend.py`, `core.py`, `cli.py`.
4. Write and run tests against the committed fixture.
5. Write `README.md`, `CLAUDE.md`, `CHANGELOG.md`, the skill file.
6. Manually exercise the documented use cases end-to-end (`uv run
   histo-contacts ...`), inspect the JSON output in `tmp/`.
7. Present for approval, then commit.
