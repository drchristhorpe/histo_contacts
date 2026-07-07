# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`histo_contacts` computes a **bond-typed contact map** between a query
selection (the "from" side — chains and/or specific residues) and one or
more target chains (the "to" side) of a 3D biological structure (PDB/mmCIF),
using [PDBe Arpeggio](https://github.com/PDBeurope/arpeggio) to classify
every contact by CREDO-derived physicochemical interaction type (`hbond`,
`ionic`, `vdw`, `hydrophobic`, `aromatic`, `covalent`, etc.) — the "bond
types" — in addition to distance. It's the sibling of `histo_neighbours`
(plain distance-cutoff contacts, no bond typing), sharing the same
query/target selector contract. It's a Python library, a Click-based CLI
(`histo-contacts`) with Rich console output, and a Claude skill wrapping
that CLI (`skills/histo-contacts/`). See [README.md](README.md) for
user-facing usage and [docs/PLAN.md](docs/PLAN.md) for the full design
rationale, including the empirical findings behind several non-obvious
decisions below.

## Environment

- Python 3.14, managed with `uv`. Use `uv sync`, `uv run <cmd>`, `uv run pytest`.
- Don't invoke a bare `python`/`pip` — always go through `uv run` /
  `uv add` so the lockfile stays authoritative.
- Dependencies include `pdbe-arpeggio` and `openbabel` (Python bindings) —
  both install cleanly from PyPI wheels on Python 3.14/macOS arm64, no
  Conda/Homebrew system package required despite Arpeggio's own README
  recommending a Conda install.

## Layout

```
src/histo_contacts/
  core.py              # load_structure(), ContactMapper, contact_map() convenience wrapper
  selectors.py         # query ("from") and target ("to") selector parsing — same grammar as histo_neighbours
  arpeggio_backend.py  # gemmi mmCIF normalization, selector translation, Arpeggio invocation, row mapping
  cli.py                # Click CLI + Rich table output (entry point: histo-contacts)
tests/
  fixtures/            # real structure file used by the tests, keep committed
  test_core.py
  test_selectors.py
  test_arpeggio_backend.py
  test_cli.py
skills/histo-contacts/SKILL.md
```

## Key invariants — don't break these

- **Every input is normalized through gemmi before reaching Arpeggio**, via
  `arpeggio_backend.normalize_to_mmcif()` — unconditionally, not just when
  the source lacks `_chem_comp.`. Arpeggio raises `ValueError: Missing
  _chem_comp. category in mmcif` on any file without that category
  (confirmed: this includes `histo_neighbours`-style coordinate-only mmCIF
  downloads), and cannot parse legacy PDB syntax at all (`gemmi.cif.read`
  fails immediately). This is not defensive/optional — don't skip it or
  make it conditional.
- The gemmi round-trip leaves `auth_asym_id`/`auth_seq_id`/coordinates
  unchanged (verified), so it's safe to build Arpeggio selection strings and
  filter results using the *original* chain ids and residue numbers.
- **A fresh `InteractionComplex` is built per `contact_map()` call.**
  `ContactMapper` still parses the structure and normalizes the mmCIF once
  in `__init__` (like `histo_neighbours`), but `arpeggio_backend.run_arpeggio()`
  constructs a new `InteractionComplex` every call. This is because
  `run_arpeggio()`'s internal result lists (`atom_contacts`,
  `plane_plane_contacts`, etc.) are only initialized in `__init__` and are
  never reset — reusing one instance across two different queries silently
  accumulates/duplicates results from the first query into the second
  (confirmed empirically). Don't "optimize" this into a single shared
  `InteractionComplex` without re-verifying that risk.
- Arpeggio's own tunables (`vdw_comp=0.1`, `ph=7.4`,
  `include_sequence_adjacent=False`) are **hardcoded to Arpeggio's own CLI
  defaults**, not exposed as options — see Scope below.
- Arpeggio's selection syntax has **no residue-range support**
  (`/P/1-3/` raises `SelectionError`, confirmed) — a query range is always
  expanded into one `/chain/resnum/` string per residue *actually present*
  in that chain within the requested range before being handed to Arpeggio.
  A whole-chain query token becomes a single `/chain//` string.
- **Query/target membership is residue-level, not atom-level**, because
  Arpeggio's plane/group-level rows represent a whole ring or amide group
  per side (comma-joined `auth_atom_id`), not a single atom. A query token
  always expands to whole residues already, so this loses nothing.
- **Output rows cover all five of Arpeggio's interaction-type granularities**
  (`atom-atom`, `atom-plane`, `plane-plane`, `group-group`, `group-plane` —
  see the `type` field), not just atom-atom. This is a deliberate,
  confirmed choice — broader scope than `histo_neighbours`' pure atom-pair
  contract. Consequence: `from_atom`/`to_atom` can be a comma-joined
  multi-atom string (a ring or group) rather than always a single atom name
  — don't assume it's always one atom.
- **Only `interacting_entities == "INTER"` rows are kept**, mapped into
  `from_*`/`to_*` by checking which side's `(chain, residue)` is in the
  resolved query set (Arpeggio doesn't guarantee `bgn` is always the
  selection side — confirmed empirically), then requiring the other side's
  chain to be in the resolved target chain set.
- **Overlap handling differs from `histo_neighbours`, and cannot be made to
  match it**: when query and target overlap, Arpeggio classifies every pair
  within the query selection as `INTRA_SELECTION`, which this tool always
  excludes — including distinct, non-identical atom pairs that
  `histo_neighbours` *would* report (it does its own pairwise
  `NeighborSearch` and only excludes literal self-pairs). This is a forced
  consequence of delegating contact classification to Arpeggio, not a
  configurable behaviour — don't try to "fix" it without reimplementing
  Arpeggio's own pairwise geometry.
- `distance` in output rows is **Arpeggio's own computed value**, used
  as-is (already rounded to 2dp internally), not independently recomputed —
  unlike `histo_neighbours`, which computes its own Euclidean distance.
  Ring/group-level distances aren't a simple atom-atom Euclidean distance,
  so there's nothing to recompute against.
- Only the first model (`structure[0]`) is used for Biopython-side selector
  resolution, same as `histo_neighbours`.
- Sort order is `(to_chain, to_residue, to_atom, from_chain, from_residue,
  from_atom)` ascending — same convention as `histo_neighbours`. Only the
  first two fields are part of the requested contract.
- The JSON output (`write_json()` / `--output`) is a **plain list** of row
  dicts — no wrapping metadata object, same as `histo_neighbours`.

## Testing

- `uv run pytest` — the fixture is already committed in `tests/fixtures/`;
  tests don't hit the network, but they do run real Arpeggio computations
  (no mocking of Arpeggio's output), so the suite takes longer than
  `histo_neighbours`' (tens of seconds, not milliseconds).
- `8gvi_1_aligned.pdb` is the same pMHC-TCR complex as `histo_neighbours`'
  fixture (chains `H`, `L`, `P`, `A`, `B`), but the coordinate-only **PDB**
  download from `coordinates.histo.fyi` rather than coordinate-only mmCIF —
  Arpeggio needs the gemmi-normalization path exercised regardless of
  format, and this is the tool's own documented data source.
- The committed fixture's default queries only ever produced `atom-atom`
  rows during development — `atom-plane`/`plane-plane`/`group-group`/
  `group-plane` are implemented and exercised against Arpeggio's documented
  JSON shape (see `test_arpeggio_backend.py`'s `map_contacts` unit tests,
  which use literal Arpeggio-shaped dicts rather than requiring a real ring
  contact to exist in the fixture), but don't assume the fixture will ever
  produce one — ring-stacking/amide-proximity are genuinely rare geometric
  events, not a coverage gap.
- When changing selector parsing or the contact-mapping logic, add a case to
  `test_selectors.py`/`test_arpeggio_backend.py`/`test_core.py` rather than
  only eyeballing CLI output.

## Scope

The CLI intentionally exposes exactly four options: `filename`, `--query`,
`--target`, `--output`, `--distance` — same surface as `histo_neighbours`.
Don't add further options (e.g. Arpeggio's own `vdw-comp`, `ph`,
`minimise-hydrogens`, `include-sequence-adjacent`, `use-ambiguities` flags,
residue-level target selectors, or a toggle to restrict output to
`atom-atom` only) without checking with the user first — these were
deliberate constraints from the initial design conversation, not
oversights.
