import tempfile
from pathlib import Path

import gemmi

from histo_contacts import arpeggio_backend

FIXTURE = Path(__file__).parent / "fixtures" / "8gvi_1_aligned.pdb"


def test_normalize_to_mmcif_produces_chem_comp_category(tmp_path):
    dest = arpeggio_backend.normalize_to_mmcif(FIXTURE, tmp_path)
    assert dest.exists()

    block = gemmi.cif.read(str(dest)).sole_block()
    assert "_chem_comp." in block.get_mmcif_category_names()


def test_normalize_to_mmcif_preserves_chain_and_residue_numbering(tmp_path):
    dest = arpeggio_backend.normalize_to_mmcif(FIXTURE, tmp_path)

    doc = gemmi.cif.read(str(dest))
    block = doc.sole_block()
    atom_site = block.get_mmcif_category("_atom_site.")

    # First atom of the fixture: ATOM 1 N GLN A 3 (auth chain A, auth resnum 3).
    assert atom_site["auth_asym_id"][0] == "A"
    assert atom_site["auth_seq_id"][0] == "3"
    assert atom_site["label_atom_id"][0] == "N"


def test_normalize_to_mmcif_is_arpeggio_loadable(tmp_path):
    from arpeggio.core import InteractionComplex

    dest = arpeggio_backend.normalize_to_mmcif(FIXTURE, tmp_path)
    # Would raise ValueError("Missing _chem_comp. category in mmcif") otherwise.
    complex_ = InteractionComplex(str(dest), vdw_comp=0.1, interacting=5.0, ph=7.4)
    assert complex_ is not None


def test_selection_string_for_chain():
    assert arpeggio_backend.selection_string_for_chain("H") == "/H//"


def test_selection_strings_for_residues():
    assert arpeggio_backend.selection_strings_for_residues("P", [1, 2, 3]) == [
        "/P/1/",
        "/P/2/",
        "/P/3/",
    ]


def test_map_contacts_picks_query_side_regardless_of_bgn_end_order():
    query_keys = {("P", 5)}
    target_chain_ids = {"H"}

    # bgn is the non-query side here, matching Arpeggio's observed behaviour
    # of not always putting the selection atom first.
    raw = [
        {
            "type": "atom-atom",
            "bgn": {
                "auth_asym_id": "H",
                "auth_seq_id": 63,
                "auth_atom_id": "OE1",
                "label_comp_id": "GLU",
            },
            "end": {
                "auth_asym_id": "P",
                "auth_seq_id": 5,
                "auth_atom_id": "CA",
                "label_comp_id": "THR",
            },
            "distance": 3.84,
            "contact": ["hbond", "polar"],
            "interacting_entities": "INTER",
        }
    ]

    rows = arpeggio_backend.map_contacts(raw, query_keys, target_chain_ids)
    assert len(rows) == 1
    row = rows[0]
    assert row["from_chain"] == "P"
    assert row["from_residue"] == 5
    assert row["to_chain"] == "H"
    assert row["to_residue"] == 63
    assert row["bond_types"] == ["hbond", "polar"]


def test_map_contacts_excludes_non_inter():
    query_keys = {("P", 5)}
    target_chain_ids = {"H"}
    raw = [
        {
            "type": "atom-atom",
            "bgn": {"auth_asym_id": "P", "auth_seq_id": 5, "auth_atom_id": "CA", "label_comp_id": "THR"},
            "end": {"auth_asym_id": "P", "auth_seq_id": 6, "auth_atom_id": "CA", "label_comp_id": "TRP"},
            "distance": 3.8,
            "contact": ["proximal"],
            "interacting_entities": "INTRA_SELECTION",
        }
    ]
    assert arpeggio_backend.map_contacts(raw, query_keys, target_chain_ids) == []


def test_map_contacts_excludes_contacts_outside_target_chains():
    query_keys = {("P", 5)}
    target_chain_ids = {"H"}
    raw = [
        {
            "type": "atom-atom",
            "bgn": {"auth_asym_id": "P", "auth_seq_id": 5, "auth_atom_id": "CA", "label_comp_id": "THR"},
            "end": {"auth_asym_id": "L", "auth_seq_id": 10, "auth_atom_id": "CB", "label_comp_id": "SER"},
            "distance": 4.5,
            "contact": ["vdw"],
            "interacting_entities": "INTER",
        }
    ]
    assert arpeggio_backend.map_contacts(raw, query_keys, target_chain_ids) == []
