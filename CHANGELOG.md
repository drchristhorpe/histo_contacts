# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.1.0] - 2026-07-07

### Added

- `histo_contacts` Python library: `ContactMapper` class that parses a
  PDB/mmCIF structure once, normalizes it once for
  [PDBe Arpeggio](https://github.com/PDBeurope/arpeggio) via a gemmi
  round-trip, and computes a bond-typed contact map from a query selection
  (chains and/or specific residues, the "from" side) to one or more target
  chains (the "to" side).
- Every contact row carries Arpeggio's own CREDO-derived `bond_types` list
  (`hbond`, `ionic`, `vdw`, `hydrophobic`, `aromatic`, `covalent`, etc.) in
  addition to distance, at atom-atom, atom-plane, plane-plane, group-group,
  and group-plane granularity (`type` field per row).
- Selector grammar (`histo_contacts.selectors`): identical to
  `histo_neighbours` — query tokens support bare chains, chain-scoped single
  residues, and chain-scoped ranges (`A`, `A:12`, `A:1-180`), comma-joined
  into one atom set; target is a simple comma-separated list of chain ids.
- `histo_contacts.arpeggio_backend`: gemmi-based mmCIF normalization
  (`normalize_to_mmcif`), query-selector-to-Arpeggio-selection-string
  translation, Arpeggio invocation (`InteractionComplex`), and mapping of
  Arpeggio's raw JSON into this tool's row shape (`map_contacts`), filtering
  to `INTER` contacts between the resolved query residues and target chains.
- A fresh `InteractionComplex` is built per `contact_map()` call (unlike
  `histo_neighbours`' single-parse-many-queries `ContactMapper`) because
  Arpeggio's `run_arpeggio()` accumulates into per-instance result lists
  that are never reset — reusing one instance across two different queries
  would silently merge their results.
- `histo-contacts` CLI (Click-based, Rich console table output with a
  `bond_types` column) with `--query`, `--target`, `--output`, and
  `--distance` options — same four-plus-one surface as `histo_neighbours`.
- Claude Code / Claude Desktop skill (`skills/histo-contacts/`) wrapping the
  CLI.
- Test suite (pytest) against a committed pMHC-TCR structure
  (`8gvi_1_aligned.pdb`, the coordinate-only PDB download from
  `coordinates.histo.fyi` — the same entry as `histo_neighbours`' fixture,
  in legacy PDB rather than coordinate-only mmCIF format).
- `README.md`, `CLAUDE.md`, and design plan (`docs/PLAN.md`).
