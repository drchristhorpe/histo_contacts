"""Arpeggio invocation: mmCIF normalization, selector translation, and
mapping Arpeggio's raw contact JSON into this tool's row shape.

See docs/PLAN.md (sections 3, 4, 5, 6) for the rationale behind each step
here — in particular why every input is round-tripped through gemmi before
reaching Arpeggio, and why a fresh InteractionComplex is built per call.
"""

from __future__ import annotations

from pathlib import Path

import gemmi
from arpeggio.core import InteractionComplex

# Arpeggio's own fixed parameters (its own CLI defaults) — not exposed as
# options here, see docs/PLAN.md section 5 and CLAUDE.md Scope.
_VDW_COMP = 0.1
_PH = 7.4
_INCLUDE_SEQUENCE_ADJACENT = False

INTERACTION_TYPES = ("atom-atom", "atom-plane", "plane-plane", "group-group", "group-plane")


def normalize_to_mmcif(path: str | Path, dest_dir: str | Path) -> Path:
    """Round-trip ``path`` (PDB or mmCIF, with or without a ``_chem_comp``
    category) through gemmi into a fresh mmCIF file Arpeggio can load.

    Arpeggio hard-requires a populated ``_chem_comp.`` mmCIF category and
    cannot read legacy PDB syntax at all; this always produces one,
    regardless of whether the source already had it, so behaviour doesn't
    branch on input format.
    """
    structure = gemmi.read_structure(str(path))
    structure.setup_entities()
    dest = Path(dest_dir) / "structure.cif"
    structure.make_mmcif_document().write_file(str(dest))
    return dest


def selection_string_for_chain(chain_id: str) -> str:
    """Arpeggio selection string for an entire chain."""
    return f"/{chain_id}//"


def selection_strings_for_residues(chain_id: str, resnums) -> list[str]:
    """Arpeggio selection strings, one per residue.

    Arpeggio's own selection syntax has no residue-range support (confirmed:
    ``/P/1-3/`` raises ``SelectionError``), so a range must be expanded to
    one string per residue; selections are additive, so this is safe.
    """
    return [f"/{chain_id}/{resnum}/" for resnum in resnums]


def run_arpeggio(cif_path: str | Path, selections: list[str], distance: float) -> list[dict]:
    """Run Arpeggio against ``cif_path`` restricted to ``selections`` and
    return its raw contact list (all five interaction-type shapes, every
    ``interacting_entities`` category, unfiltered)."""
    complex_ = InteractionComplex(str(cif_path), vdw_comp=_VDW_COMP, interacting=distance, ph=_PH)
    complex_.structure_checks()
    complex_.initialize()
    complex_.run_arpeggio(selections, distance, _VDW_COMP, _INCLUDE_SEQUENCE_ADJACENT)
    return complex_.get_contacts()


def map_contacts(
    raw_contacts: list[dict], query_keys: set[tuple[str, int]], target_chain_ids: set[str]
) -> list[dict]:
    """Filter Arpeggio's raw contacts to query -> target ``INTER`` rows and
    map them into this tool's output row shape (see docs/PLAN.md section 7).

    ``query_keys`` is a set of ``(auth_asym_id, auth_seq_id)`` for every
    residue in the resolved query selection. Arpeggio's own ``bgn``/``end``
    order is not consistent about which side is the selection side, so
    membership is checked explicitly rather than assumed.
    """
    rows = []
    for contact in raw_contacts:
        if contact["interacting_entities"] != "INTER":
            continue

        bgn, end = contact["bgn"], contact["end"]
        bgn_key = (bgn["auth_asym_id"], bgn["auth_seq_id"])
        end_key = (end["auth_asym_id"], end["auth_seq_id"])

        if bgn_key in query_keys:
            frm, to = bgn, end
        elif end_key in query_keys:
            frm, to = end, bgn
        else:
            continue

        if to["auth_asym_id"] not in target_chain_ids:
            continue

        rows.append(
            {
                "type": contact["type"],
                "from_chain": frm["auth_asym_id"],
                "from_residue": frm["auth_seq_id"],
                "from_atom": frm["auth_atom_id"],
                "from_aa": frm["label_comp_id"],
                "to_chain": to["auth_asym_id"],
                "to_residue": to["auth_seq_id"],
                "to_atom": to["auth_atom_id"],
                "to_aa": to["label_comp_id"],
                "distance": float(contact["distance"]),
                "bond_types": sorted(contact["contact"]),
            }
        )
    return rows
