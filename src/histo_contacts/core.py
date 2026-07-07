"""Structure loading and Arpeggio-backed, bond-typed contact-map calculations."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

from Bio.PDB.Chain import Chain
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

from histo_contacts import arpeggio_backend
from histo_contacts.selectors import QueryPart, SelectorError, parse_query, parse_target

_CIF_SUFFIXES = {".cif", ".mmcif"}
_PDB_SUFFIXES = {".pdb", ".ent"}


class StructureError(ValueError):
    """Raised for problems loading or querying a structure."""


def load_structure(path: str | Path, structure_id: str | None = None) -> Structure:
    """Parse a PDB or mmCIF file into a Bio.PDB Structure.

    Format is chosen from the file extension (case-insensitive):
    ``.cif``/``.mmcif`` -> mmCIF, ``.pdb``/``.ent`` -> legacy PDB.
    """
    path = Path(path)
    if not path.is_file():
        raise StructureError(f"No such file: {path}")

    suffix = path.suffix.lower()
    sid = structure_id or path.stem

    if suffix in _CIF_SUFFIXES:
        parser = MMCIFParser(QUIET=True)
    elif suffix in _PDB_SUFFIXES:
        parser = PDBParser(QUIET=True)
    else:
        raise StructureError(
            f"Unrecognised structure file extension {suffix!r} for {path}; "
            "expected one of .cif, .mmcif, .pdb, .ent"
        )

    structure = parser.get_structure(sid, str(path))
    if len(structure) == 0:
        raise StructureError(f"No models found in {path}")
    return structure


def _get_chain(model, chain_id: str) -> Chain:
    try:
        return model[chain_id]
    except KeyError:
        available = ", ".join(c.id for c in model)
        raise StructureError(
            f"No such chain {chain_id!r}; available chains: {available}"
        ) from None


def _default_chain(model) -> Chain:
    chains = list(model)
    if len(chains) != 1:
        ids = ", ".join(c.id for c in chains)
        raise StructureError(
            "A chain-less selector is only valid for single-chain structures; "
            f"this structure has {len(chains)} chains ({ids}). "
            "Prefix the selector with a chain id, e.g. 'A:1-180'."
        )
    return chains[0]


def _residues_for_part(model, part: QueryPart) -> tuple[Chain, list]:
    chain = _get_chain(model, part.chain) if part.chain else _default_chain(model)

    if part.start is None:
        residues = list(chain)
    else:
        residues = [r for r in chain if part.start <= r.id[1] <= part.end]
        if not residues:
            range_text = str(part.start) if part.start == part.end else f"{part.start}-{part.end}"
            raise StructureError(
                f"No residues in range {range_text} found in chain {chain.id!r}"
            )
    return chain, residues


class ContactMapper:
    """Loads a structure once, normalizes it for Arpeggio once, and computes
    bond-typed contact maps from it.

    >>> cm = ContactMapper("structure.pdb")
    >>> cm.contact_map("P", ["H", "L"], distance=5.0)
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.structure = load_structure(self.path)
        self.model = self.structure[0]
        self._tmpdir = tempfile.TemporaryDirectory(prefix="histo_contacts_")
        self._arpeggio_cif_path = arpeggio_backend.normalize_to_mmcif(self.path, self._tmpdir.name)

    def _resolve_query(self, query: str) -> tuple[set[tuple[str, int]], list[str]]:
        parts = parse_query(query)
        query_keys: set[tuple[str, int]] = set()
        selections: list[str] = []
        for part in parts:
            chain, residues = _residues_for_part(self.model, part)
            query_keys.update((chain.id, r.id[1]) for r in residues)
            if part.start is None:
                selections.append(arpeggio_backend.selection_string_for_chain(chain.id))
            else:
                selections.extend(
                    arpeggio_backend.selection_strings_for_residues(chain.id, [r.id[1] for r in residues])
                )
        if not query_keys:
            raise StructureError(f"Query selector {query!r} matched no atoms")
        return query_keys, selections

    def _resolve_target_chains(self, target) -> list[str]:
        chain_ids = parse_target(target) if isinstance(target, str) else list(target)
        if not chain_ids:
            raise SelectorError("Empty target selector")
        for chain_id in chain_ids:
            _get_chain(self.model, chain_id)
        return chain_ids

    def contact_map(self, query: str, target, distance: float = 5.0) -> list[dict]:
        """Compute the bond-typed contact map from ``query`` (the "from"
        selection: chains and/or specific residues) to ``target`` (the "to"
        selection: a list of chain ids), for every contact Arpeggio reports
        within ``distance`` Angstroms.

        Returns a list of row dicts (``type``, ``from_chain``,
        ``from_residue``, ``from_atom``, ``from_aa``, ``to_chain``,
        ``to_residue``, ``to_atom``, ``to_aa``, ``distance``,
        ``bond_types``), sorted by target chain ascending, then target
        residue number ascending.
        """
        if distance <= 0:
            raise ValueError(f"distance must be positive, got {distance}")

        query_keys, selections = self._resolve_query(query)
        target_chain_ids = self._resolve_target_chains(target)

        raw_contacts = arpeggio_backend.run_arpeggio(self._arpeggio_cif_path, selections, distance)
        rows = arpeggio_backend.map_contacts(raw_contacts, query_keys, set(target_chain_ids))

        rows.sort(
            key=lambda r: (
                r["to_chain"],
                r["to_residue"],
                r["to_atom"],
                r["from_chain"],
                r["from_residue"],
                r["from_atom"],
            )
        )
        return rows

    def write_json(self, rows: list[dict], output_path: str | Path) -> Path:
        """Write contact-map rows as a JSON list to ``output_path``."""
        output_path = Path(output_path)
        output_path.write_text(json.dumps(rows, indent=2))
        return output_path


def contact_map(path: str | Path, query: str, target, distance: float = 5.0) -> list[dict]:
    """Convenience wrapper: compute a contact map for a structure file."""
    return ContactMapper(path).contact_map(query, target, distance=distance)
