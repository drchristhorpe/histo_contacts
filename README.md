# histo-contacts

Compute a **bond-typed contact map** between one part of a 3D biological
structure (PDB/mmCIF) and one or more other chains, using
[PDBe Arpeggio](https://github.com/PDBeurope/arpeggio) to classify every
contact by CREDO-derived physicochemical interaction type (`hbond`, `ionic`,
`vdw`, `hydrophobic`, `aromatic`, `covalent`, and more) in addition to
distance. It's the sibling of
[`histo-neighbours`](../histo_neighbours) (plain distance-cutoff contacts,
no bond typing) — same query/target selector contract, richer rows.

It ships as:

- a Python library — `import histo_contacts`
- a CLI tool — `histo-contacts`
- a [Claude Code / Claude Desktop skill](skills/histo-contacts/SKILL.md)

Requires Python 3.14+.

## Install

```bash
uv sync                 # dev environment, from a checkout
uv tool install .       # install the `histo-contacts` CLI globally
# or
pip install .
```

## CLI usage

```
histo-contacts FILENAME --query SPEC --target SPEC --output PATH [--distance FLOAT]
```

- `FILENAME` — a `.cif`/`.mmcif` or `.pdb`/`.ent` structure file.
- `--query`, `-q` — the **"from"** side: comma-separated chains and/or
  chain:residue / chain:range tokens.
- `--target`, `-t` — the **"to"** side: comma-separated list of chain ids.
- `--output`, `-o` — path to write the contact map as a JSON file.
- `--distance`, `-d` — contact distance cutoff in Ångströms (default `5.0`),
  forwarded to Arpeggio's own interacting-distance cutoff.

### Example: a peptide residue against an MHC chain, with bond types

```bash
$ histo-contacts 8gvi_1_aligned.pdb --query P:1 --target H,L --output contacts.json --distance 4.0
                       Contacts: 'P:1' -> 'H,L' (≤ 4.0Å)
┏━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ from_chain ┃ from_residue ┃ from_atom ┃ to_chain ┃ to_residue ┃ to_atom ┃ distance ┃ bond_types               ┃
┡━━━━━━━━━━━╇━━━━━━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━╇━━━━━━━━━━━╇━━━━━━━━╇━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ P          │ 1            │ NH2       │ H        │ 62         │ OE1     │ 3.03     │ ionic, polar, vdw_clash  │
│ P          │ 1            │ N         │ H        │ 7          │ OH      │ 2.99     │ hbond, polar, vdw_clash  │
│ ...        │ ...          │ ...       │ ...      │ ...        │ ...     │ ...      │ ...                      │
└───────────┴─────────────┴──────────┴─────────┴───────────┴────────┴──────────┴──────────────────────────┘
31 contact(s) within 4.0Å written to contacts.json
```

### Example: whole chain against two other chains

```bash
$ histo-contacts 8gvi_1_aligned.pdb --query P --target H,L --output contacts.json
```

Output rows are sorted by target chain ascending (when the target spans
multiple chains), then target residue number ascending.

## Selector grammar

Identical to `histo-neighbours`. Query (`--query`) tokens:

| Token | Meaning |
|---|---|
| `A` | every atom in chain `A` |
| `A:12` | residue 12 on chain `A` |
| `A:1-180` | residues 1–180 on chain `A` |
| `12` / `1-180` | residue/range on the *default chain* — only valid for single-chain structures |

Multiple comma-separated tokens are unioned into one query set, e.g.
`--query A:1-9,B:12` selects chain A residues 1–9 **and** chain B residue
12 as the "from" side.

Target (`--target`) is simpler by design: a comma-separated list of chain
ids only, e.g. `--target H,L`. There is no residue-level restriction on
the target side.

## JSON output format

A plain JSON list of contact rows:

```json
[
  {
    "type": "atom-atom",
    "from_chain": "P", "from_residue": 1, "from_atom": "NH2", "from_aa": "ARG",
    "to_chain": "H", "to_residue": 62, "to_atom": "OE1", "to_aa": "GLU",
    "distance": 3.03,
    "bond_types": ["ionic", "polar", "vdw_clash"]
  }
]
```

- `type` is one of Arpeggio's five interaction-type granularities:
  `atom-atom`, `atom-plane`, `plane-plane`, `group-group`, `group-plane`.
  For anything other than `atom-atom`, `from_atom`/`to_atom` holds a
  comma-joined list of atom names (a ring or functional group), not a
  single atom.
- `bond_types` is Arpeggio's own CREDO interaction-type classification for
  that pair — the "bond types". Common values: `covalent`, `clash`,
  `vdw_clash`, `vdw`, `proximal`, `hbond`, `weak_hbond`, `xbond`, `ionic`,
  `metal_complex`, `aromatic`, `hydrophobic`, `carbonyl`, `polar`,
  `weak_polar` (atom-level), plus ring/group-level types like `FT`, `ET`,
  `AMIDEAMIDE`, `AMIDERING` for the non-atom-atom rows. See the
  [Arpeggio README](https://github.com/PDBeurope/arpeggio#output) for the
  full reference.
- `distance` is Arpeggio's own computed value (2dp), not independently
  recomputed.

## Library usage

`ContactMapper` parses a structure file once and normalizes it for Arpeggio
once; the same instance can be reused for multiple `contact_map()` calls
(each call still runs a fresh Arpeggio computation — see
[CLAUDE.md](CLAUDE.md)).

```python
from histo_contacts import ContactMapper

cm = ContactMapper("8gvi_1_aligned.pdb")
rows = cm.contact_map("P", ["H", "L"], distance=5.0)   # -> list[dict]
cm.write_json(rows, "contacts.json")
```

A one-shot convenience function is also available:

```python
from histo_contacts import contact_map

contact_map("8gvi_1_aligned.pdb", "P", ["H", "L"], distance=5.0)
```

`target` accepts either a comma-separated chain-id string (`"H,L"`) or an
iterable of chain ids (`["H", "L"]`).

## Notes and limitations

- Every input (`.cif`/`.mmcif`/`.pdb`/`.ent`) is round-tripped through
  [gemmi](https://gemmi.readthedocs.io/) before reaching Arpeggio, since
  Arpeggio requires a populated `_chem_comp.` mmCIF category that
  coordinate-only downloads (and all legacy PDB files) lack. This doesn't
  change contact-finding results but does mean each `ContactMapper(...)`
  call does a small amount of extra I/O up front.
- Only the **first model** in a file is used — sufficient for X-ray/cryo-EM
  structures; NMR ensembles are not averaged/iterated across models.
- If the query and target selections overlap, contacts *within* the query
  selection are never reported (Arpeggio classifies them as
  `INTRA_SELECTION`, not `INTER`) — even distinct, non-identical atom pairs.
  This differs from `histo-neighbours`, which reports those pairs; see
  [CLAUDE.md](CLAUDE.md) for why this can't be made to match.
- Whole-chain or wide-range queries are noticeably slower than
  `histo-neighbours` (seconds to tens of seconds, not milliseconds) — the
  cost is Arpeggio's own per-atom geometric classification (hydrogenation,
  valence/ring perception, hydrogen-bond donor/acceptor/angle checks), not
  this tool's own logic.
- The target chain must exist in the structure; the query's chain-less
  tokens (bare residue numbers/ranges) are only valid for single-chain
  structures.

## Residue-level aggregation for TCR:pMHC structures

This library produces atom-level contacts with bond-type classification. For
pipelines that need residue-level contact maps grouped by CDR loops and MHC
interface regions, the `residue_aggregator` utility module provides the
canonical pattern:

```python
from histo_contacts import ContactMapper, residue_aggregator

cm = ContactMapper("structure.pdb")
atom_rows = cm.contact_map(
    "D:27-38,D:56-65,D:105-117,E:27-38,E:56-65,E:105-117",  # All TCR CDRs
    ["A", "C"],  # MHC + peptide
    distance=5.0,
)

residue_rows = residue_aggregator.aggregate_contacts(atom_rows)
residue_aggregator.write_residue_contacts_json(
    {"file": "structure.pdb", "ct_total_atom_pairs": len(atom_rows)},
    residue_rows,
    "structure_residue_contacts.json",
)
```

Key features:

- **CDR loops at fixed IMGT positions** (chains D/E, not string-matched) —
  requires IMGT-renumbered structures.
- **Helix-only MHC regions** (α1: A 50–86, α2: A 137–180, peptide: C all) —
  the TCR reads the helical surfaces, not the β-sheet floor.
- **Bond-type counts, not unions** — one atom pair carries multiple types
  simultaneously (e.g., both `proximal` and `vdw`), so the output is
  `{"proximal": 3, "vdw": 1}`, not `["proximal", "vdw"]`. This is required
  for acceptance testing: re-aggregating residue-level rows must exactly
  match the atom-pair counts they came from.
- **Validation functions** — `reconcile_against_aggregates()` checks that
  residue rows re-aggregate correctly, and `validate_no_regression()`
  ensures helix-only counts ≤ whole-domain baselines.

See [examples/aggregate_residue_contacts.py](examples/aggregate_residue_contacts.py)
for a complete working example.

## Development

```bash
uv sync
uv run pytest
```

Test fixtures under `tests/fixtures/` are real structure files from
[coordinates.histo.fyi](https://coordinates.histo.fyi/).

See [docs/PLAN.md](docs/PLAN.md) for the design rationale (including the
empirical findings behind the gemmi-normalization step and the
fresh-`InteractionComplex`-per-call design, and a detailed design for
residue-level aggregation) and [CHANGELOG.md](CHANGELOG.md)
for release history.
