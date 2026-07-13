"""Residue-level aggregation utility for TCR:pMHC contact maps.

This module aggregates atom-level contacts from histo_contacts output to
residue-level pairs, grouped by CDR loops (IMGT-fixed positions on TCR chains D/E)
and MHC interface regions (helix-only α1/α2 on chain A, or peptide on chain C).

Design decisions:
- Bond types are COUNTED per residue pair (not unioned), because one atom pair
  carries multiple types simultaneously. This is required for acceptance testing:
  SUM(bond_type_counts[type]) over residue pairs in a cell must equal the
  atom-pair count for that cell.
- CDR loops are at fixed IMGT positions (not string-matched), requiring
  IMGT-renumbered structures.
- MHC regions are helix-only (scientific choice: TCR reads the helices flanking
  the groove, not the β-sheet floor).
"""

from pathlib import Path
from collections import defaultdict
import json


def define_cdr_loops() -> dict[str, dict[str, tuple[int, int]]]:
    """Return CDR loop definitions per chain, IMGT-renumbered positions.

    Returns:
        dict: {'alpha': {'cdr1': (27, 38), ...}, 'beta': {...}}
    """
    return {
        "alpha": {
            "cdr1": (27, 38),
            "cdr2": (56, 65),
            "cdr3": (105, 117),
        },
        "beta": {
            "cdr1": (27, 38),
            "cdr2": (56, 65),
            "cdr3": (105, 117),
        },
    }


def define_mhc_regions() -> dict[str, tuple[str, int, int]]:
    """Return MHC region definitions: helix-only, Class I numbering.

    Returns:
        dict: region_name -> (chain, start_resnum, end_resnum)
              e.g., 'alpha1' -> ('A', 50, 86)
    """
    return {
        "alpha1": ("A", 50, 86),
        "alpha2": ("A", 137, 180),
        "peptide": ("C", 1, 1000),  # All of chain C; max 1000 is arbitrary upper bound
    }


def categorize_contact(
    atom_row: dict,
    cdr_loops: dict,
    mhc_regions: dict,
) -> tuple[str, str, str, str] | None:
    """Assign an atom-level contact to a (CDR loop, MHC region) pair.

    Args:
        atom_row: dict from histo_contacts, with from_chain, from_residue, to_chain, to_residue
        cdr_loops: from define_cdr_loops()
        mhc_regions: from define_mhc_regions()

    Returns:
        (cdr_chain, cdr_loop_name, mhc_region_name, from_or_to)
            where cdr_loop_name includes the alpha/beta suffix (e.g., "cdr1_alpha")
            and from_or_to is either "from" or "to" (which side is the TCR)
        or None if the row doesn't fall into CDR-region scope.
    """
    from_chain = atom_row["from_chain"]
    from_resnum = atom_row["from_residue"]
    to_chain = atom_row["to_chain"]
    to_resnum = atom_row["to_residue"]

    # Check if from_side is a TCR CDR
    if from_chain in ("D", "E"):
        tcr_label = "alpha" if from_chain == "D" else "beta"
        cdr_loop = _find_cdr_loop(from_resnum, cdr_loops[tcr_label])
        if cdr_loop:
            mhc_region = _find_mhc_region(to_resnum, to_chain, mhc_regions)
            if mhc_region:
                cdr_loop_with_label = f"{cdr_loop}_{tcr_label}"
                return (from_chain, cdr_loop_with_label, mhc_region, "from")

    # Check if to_side is a TCR CDR
    if to_chain in ("D", "E"):
        tcr_label = "alpha" if to_chain == "D" else "beta"
        cdr_loop = _find_cdr_loop(to_resnum, cdr_loops[tcr_label])
        if cdr_loop:
            mhc_region = _find_mhc_region(from_resnum, from_chain, mhc_regions)
            if mhc_region:
                cdr_loop_with_label = f"{cdr_loop}_{tcr_label}"
                return (to_chain, cdr_loop_with_label, mhc_region, "to")

    return None


def _find_cdr_loop(resnum: int, cdr_ranges: dict[str, tuple[int, int]]) -> str | None:
    """Find which CDR loop contains resnum, or None."""
    for loop_name, (start, end) in cdr_ranges.items():
        if start <= resnum <= end:
            return loop_name
    return None


def _find_mhc_region(resnum: int, chain: str, mhc_regions: dict) -> str | None:
    """Find which MHC region contains (chain, resnum), or None."""
    for region_name, (region_chain, start, end) in mhc_regions.items():
        if chain == region_chain and start <= resnum <= end:
            return region_name
    return None


def aggregate_contacts(atom_rows: list[dict]) -> dict[str, list[dict]]:
    """Aggregate atom-level contacts to residue-level pairs.

    Groups atom pairs by (cdr_loop, mhc_region), then by (tcr_residue, mhc_residue),
    and for each residue pair:
    - counts bond types (one atom pair may carry multiple types)
    - counts atom pairs
    - records min distance
    - records residue names

    Args:
        atom_rows: list of dicts from histo_contacts.contact_map()

    Returns:
        dict[residue_pair_key] -> list of residue-level contact dicts, where
        residue_pair_key = (cdr_loop, mhc_region, tcr_chain, tcr_resnum, mhc_chain, mhc_resnum)
    """
    cdr_loops = define_cdr_loops()
    mhc_regions = define_mhc_regions()

    # Group by (cdr_loop, region, tcr_chain, tcr_resnum, mhc_chain, mhc_resnum)
    aggregated: dict = defaultdict(lambda: {
        "atom_pairs": [],
        "bond_type_counts": defaultdict(int),
        "distances": [],
    })

    for atom_row in atom_rows:
        result = categorize_contact(atom_row, cdr_loops, mhc_regions)
        if not result:
            continue

        tcr_chain, cdr_loop, mhc_region, from_or_to = result

        # Extract TCR and MHC sides
        if from_or_to == "from":
            tcr_side = atom_row["from_chain"], atom_row["from_residue"], atom_row["from_aa"]
            mhc_side = atom_row["to_chain"], atom_row["to_residue"], atom_row["to_aa"]
        else:
            tcr_side = atom_row["to_chain"], atom_row["to_residue"], atom_row["to_aa"]
            mhc_side = atom_row["from_chain"], atom_row["from_residue"], atom_row["from_aa"]

        tcr_chain_id, tcr_resnum, tcr_resname = tcr_side
        mhc_chain_id, mhc_resnum, mhc_resname = mhc_side

        key = (cdr_loop, mhc_region, tcr_chain_id, tcr_resnum, mhc_chain_id, mhc_resnum)

        # Track the atom pair and its bond types
        aggregated[key]["atom_pairs"].append(atom_row)
        aggregated[key]["distances"].append(atom_row["distance"])
        aggregated[key]["tcr_resname"] = tcr_resname
        aggregated[key]["mhc_resname"] = mhc_resname

        for bond_type in atom_row["bond_types"]:
            aggregated[key]["bond_type_counts"][bond_type] += 1

    # Convert to final residue-level rows
    residue_rows = []
    for (cdr_loop, mhc_region, tcr_chain, tcr_resnum, mhc_chain, mhc_resnum), data in aggregated.items():
        row = {
            "cdr_loop": cdr_loop,
            "region": mhc_region,
            "tcr_chain": tcr_chain,
            "tcr_resnum": tcr_resnum,
            "tcr_resname": data["tcr_resname"],
            "mhc_chain": mhc_chain,
            "mhc_resnum": mhc_resnum,
            "mhc_resname": data["mhc_resname"],
            "n_atom_pairs": len(data["atom_pairs"]),
            "bond_type_counts": dict(data["bond_type_counts"]),
            "min_distance": round(min(data["distances"]), 2),
        }
        residue_rows.append(row)

    return sorted(residue_rows, key=lambda r: (
        r["cdr_loop"],
        r["region"],
        r["tcr_resnum"],
        r["mhc_resnum"],
    ))


def write_residue_contacts_json(
    structure_info: dict,
    residue_rows: list[dict],
    path: Path,
) -> None:
    """Write residue-level contacts to JSON.

    Args:
        structure_info: dict with pdb_id, complex, file, ct_total_atom_pairs
        residue_rows: list of residue-level contact dicts
        path: output file path
    """
    output = {**structure_info, "contacts": residue_rows}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(output, f, indent=2)


def reconcile_against_aggregates(
    residue_rows: list[dict],
    expected_atom_pair_counts: dict[str, int],
) -> bool:
    """Validate residue-level counts against expected atom-pair totals.

    For each (cdr_loop, region) cell, sums bond_type_counts over all residue
    pairs in that cell and checks it equals the expected count.

    Args:
        residue_rows: list of residue-level dicts
        expected_atom_pair_counts: dict[f"{cdr_loop}|{region}"] -> count

    Returns:
        True if all cells reconcile, False otherwise.
    """
    cells: dict[str, int] = defaultdict(int)

    # Aggregate residue-level counts by bond type, per cell
    for row in residue_rows:
        cell = f"{row['cdr_loop']}|{row['region']}"
        for bond_type, count in row["bond_type_counts"].items():
            cells[cell] += count

    # Check against expected
    for cell, count in expected_atom_pair_counts.items():
        if cells[cell] != count:
            print(f"Reconciliation FAILED for {cell}: expected {count}, got {cells[cell]}")
            return False

    return True


def validate_no_regression(
    helix_only_counts: dict[str, int],
    whole_domain_counts: dict[str, int],
) -> bool:
    """Ensure helix-only counts <= whole-domain counts (strict subset).

    Args:
        helix_only_counts: dict[cell_key] -> count (new, helix-only)
        whole_domain_counts: dict[cell_key] -> count (old, whole-domain baseline)

    Returns:
        True if all helix-only <= whole-domain, False otherwise.
    """
    for cell, new_count in helix_only_counts.items():
        old_count = whole_domain_counts.get(cell, 0)
        if new_count > old_count:
            print(f"Regression for {cell}: old={old_count}, new={new_count} (INCREASED)")
            return False
    return True
