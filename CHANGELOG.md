# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- `histo_contacts.residue_aggregator` — utility module for aggregating
  atom-level contacts to residue-level pairs, with TCR:pMHC grouping by CDR
  loops (IMGT-fixed positions on chains D/E) and MHC interface regions
  (helix-only α1/α2 on chain A, peptide on chain C). Particularly useful for
  pipeline consumers that need residue-level contact maps.
- `aggregate_contacts()` — main entry point: takes atom-level rows from
  `contact_map()` and produces residue-level rows with bond-type counts
  (not unions), required for acceptance testing (self-consistency across
  re-aggregated totals).
- `categorize_contact()` — classifies each atom-level contact into a
  (CDR loop, MHC region) pair, or None if out of scope.
- `reconcile_against_aggregates()` — validates residue-level rows re-aggregate
  to expected atom-pair totals per (CDR, region) cell.
- `validate_no_regression()` — ensures helix-only counts ≤ whole-domain counts
  (strict subset, no expansion allowed).
- `write_residue_contacts_json()` — outputs residue-level contacts in the
  standard JSON structure (keyed by structure ID, with per-structure metadata).
- Example script `examples/aggregate_residue_contacts.py` — demonstrates the
  canonical pattern: compute atom-level contacts, aggregate to residue level,
  write JSON.
- Comprehensive test suite for aggregation (`tests/test_residue_aggregator.py`):
  26 unit tests covering definitions, categorization, aggregation, validation,
  and an acceptance test against sample 1AO7 data.

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
