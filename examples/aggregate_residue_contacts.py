#!/usr/bin/env python3
"""Example: aggregate atom-level contacts to residue-level via residue_aggregator.

This script demonstrates the canonical pattern for using histo_contacts output
to produce residue-level contact maps, grouped by TCR CDR loops and MHC
interface regions.

Typical workflow:
  1. Use histo_contacts.contact_map() to get atom-level contacts
  2. Pass them to residue_aggregator.aggregate_contacts()
  3. Validate the output with acceptance tests
  4. Write to JSON

Example invocation:
  $ python examples/aggregate_residue_contacts.py tests/fixtures/8gvi_1_aligned.pdb
"""

import json
from pathlib import Path
import sys

from histo_contacts import ContactMapper, residue_aggregator


def main(pdb_path: str | None = None) -> int:
    """Aggregate atom-level contacts to residue level for a TCR:pMHC structure.

    Args:
        pdb_path: path to input PDB/mmCIF file, or None to use test fixture

    Returns:
        0 on success, 1 on error
    """
    # Use test fixture if no path provided
    if pdb_path is None:
        pdb_path = "tests/fixtures/8gvi_1_aligned.pdb"
        if not Path(pdb_path).exists():
            print(f"Error: {pdb_path} not found", file=sys.stderr)
            return 1

    pdb_path = Path(pdb_path)
    print(f"Reading structure: {pdb_path}")

    # Compute atom-level contacts using histo_contacts
    cm = ContactMapper(str(pdb_path))
    print("Computing atom-level contacts...")
    atom_rows = cm.contact_map(
        "D:27-38,D:56-65,D:105-117,E:27-38,E:56-65,E:105-117",  # All TCR CDRs
        ["A", "C"],  # MHC heavy chain + peptide
        distance=5.0,
    )
    print(f"  {len(atom_rows)} atom-level contacts found")

    # Aggregate to residue level
    print("Aggregating to residue level...")
    residue_rows = residue_aggregator.aggregate_contacts(atom_rows)
    print(f"  {len(residue_rows)} unique residue pairs")

    # Compute totals per (CDR, region) cell for validation
    cell_totals = {}
    for row in residue_rows:
        cell = f"{row['cdr_loop']}|{row['region']}"
        cell_totals[cell] = cell_totals.get(cell, 0) + row["n_atom_pairs"]

    print("\nContact summary by CDR loop and MHC region:")
    print("  CDR loop       | MHC region | Atom pairs")
    print("  " + "-" * 44)
    for cell in sorted(cell_totals.keys()):
        cdr_loop, region = cell.split("|")
        count = cell_totals[cell]
        print(f"  {cdr_loop:14} | {region:10} | {count:4d}")

    total_atoms = sum(row["n_atom_pairs"] for row in residue_rows)
    print(f"  {'TOTAL':14} | {'':10} | {total_atoms:4d}")

    # Write output JSON
    output_path = Path("tmp") / f"{pdb_path.stem}_residue_contacts.json"
    output_path.parent.mkdir(exist_ok=True)

    structure_info = {
        "file": pdb_path.name,
        "ct_total_atom_pairs": total_atoms,
    }

    residue_aggregator.write_residue_contacts_json(structure_info, residue_rows, output_path)
    print(f"\nWrote residue-level contacts to {output_path}")

    # Show sample row
    if residue_rows:
        print("\nSample residue-level row:")
        sample = residue_rows[0]
        print(json.dumps(sample, indent=2))

    return 0


if __name__ == "__main__":
    pdb_path = sys.argv[1] if len(sys.argv) > 1 else None
    sys.exit(main(pdb_path))
