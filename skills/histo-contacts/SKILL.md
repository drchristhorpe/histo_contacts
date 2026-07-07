---
name: histo-contacts
description: Compute a bond-typed contact map (hydrogen bonds, ionic, van der Waals, hydrophobic, aromatic, covalent, etc., via PDBe Arpeggio) between a query selection (chains and/or specific residues) and one or more target chains of a 3D biological structure (PDB/mmCIF file). Use when asked what kind of interaction/bond exists between residues or chains, to classify contacts by physicochemical type, to find hydrogen bonds/salt bridges/hydrophobic contacts at an interface, or when a plain distance-based contact map (no bond typing) isn't enough.
---

# histo-contacts

`histo-contacts` is a CLI tool (installed from the `histo_contacts` package)
that computes a bond-typed contact map between a query selection and one or
more target chains in a PDB or mmCIF structure file, using
[PDBe Arpeggio](https://github.com/PDBeurope/arpeggio) to classify each
contact by physicochemical interaction type. Invoke it with the Bash tool.

It's the sibling of `histo-neighbours` (same query/target selector grammar,
plain distance-cutoff contacts, no bond typing). Use `histo-contacts` when
the *type* of interaction matters (hydrogen bond vs. ionic vs. hydrophobic
vs. just "nearby"); use `histo-neighbours` when only distance/proximity
matters and speed matters more (it's much faster — no Arpeggio overhead).

## When to use this skill

The user provides (or references) a `.cif`/`.mmcif` or `.pdb`/`.ent`
structure file and asks about the *nature* of contacts between one part of
the structure (a chain, a residue range, individual residues) and one or
more other chains — e.g. "what hydrogen bonds does this peptide make with
the MHC groove", "are there any salt bridges at this interface", "classify
the contacts between chain A and chain B".

## Checking availability

```bash
histo-contacts --help
```

If this fails with "command not found", install it first:

```bash
uv tool install histo_contacts   # or: pip install histo_contacts
```

(If working from a checkout of the `histo_contacts` source repo instead of
an installed package, use `uv run histo-contacts ...` there instead.)

## Usage

```
histo-contacts FILENAME --query SPEC --target SPEC --output PATH [--distance FLOAT]
```

- `--query`/`-q` (required): the **"from"** side — comma-separated chains
  and/or chain:residue / chain:range tokens, e.g. `P`, `A:1-9`,
  `A:1-9,B:12`. A bare chain letter means "every atom in that chain".
- `--target`/`-t` (required): the **"to"** side — a comma-separated list
  of chain ids only, e.g. `H,L`. No residue-level restriction on this side.
- `--output`/`-o` (required): path to write the contact map as a JSON file.
- `--distance`/`-d` (optional, default `5.0`): contact distance cutoff in
  Ångströms, forwarded to Arpeggio's own interacting-distance cutoff.

Examples:

- **Whole chain vs. one other chain**: `histo-contacts structure.cif --query P --target H --output contacts.json`
- **A single residue, tighter cutoff**: `histo-contacts structure.cif --query P:5 --target H,L --output contacts.json --distance 4.0`
- **A residue range vs. a chain**: `histo-contacts structure.cif --query A:1-9 --target B --output contacts.json`

Note: whole-chain or wide-range queries are noticeably slower than
`histo-neighbours` (Arpeggio does per-atom hydrogenation and geometric
classification) — seconds to tens of seconds rather than milliseconds. A
single residue or small range typically finishes in a couple of seconds.

## Output

The JSON file at `--output` is a plain list of contact rows:

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

- `bond_types` is the key addition over `histo-neighbours`: Arpeggio's
  CREDO-derived interaction classification for that pair. Common values —
  `hbond` (hydrogen bond), `ionic` (salt bridge), `hydrophobic`, `aromatic`
  (ring-ring/ring-atom), `vdw`/`vdw_clash` (van der Waals), `covalent`,
  `polar`/`weak_polar`, `proximal` (within cutoff but no stronger
  classification applies). A row can have several bond types at once.
- `type` distinguishes single-atom pairs (`atom-atom`, the common case)
  from ring/group-level interactions (`atom-plane`, `plane-plane`,
  `group-group`, `group-plane`) — for the latter, `from_atom`/`to_atom` may
  be a comma-joined list of atom names rather than a single atom.
- Rows are sorted by target chain (`to_chain`) ascending, then target
  residue number (`to_residue`) ascending. The CLI also prints a Rich table
  preview (capped at 25 rows) including a `bond_types` column, and a total
  contact count, to the console.

## Example

```bash
$ histo-contacts 8gvi_1_aligned.pdb --query P:1 --target H,L --output contacts.json --distance 4.0
                       Contacts: 'P:1' -> 'H,L' (≤ 4.0Å)
┏━━━━━━━━━━━┳━━━━━━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━┳━━━━━━━━━━━┳━━━━━━━━┳━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ from_chain ┃ from_residue ┃ from_atom ┃ to_chain ┃ to_residue ┃ to_atom ┃ distance ┃ bond_types               ┃
┗━━━━━━━━━━━┻━━━━━━━━━━━━━┻━━━━━━━━━━┻━━━━━━━━━┻━━━━━━━━━━━┻━━━━━━━━┻━━━━━━━━━━┻━━━━━━━━━━━━━━━━━━━━━━━━━━┛
31 contact(s) within 4.0Å written to contacts.json
```

Report the result back to the user in whatever form they asked for (a
summary count, specific hydrogen bonds/salt bridges filtered from the JSON
by `bond_types`, the full table, etc.) — this skill only tells you how to
obtain the bond-typed contact map. If the user only cares about distance,
not bond type, and speed matters, use `histo-neighbours` instead.
