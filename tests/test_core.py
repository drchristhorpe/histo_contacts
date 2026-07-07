import json
from pathlib import Path

import pytest

from histo_contacts import ContactMapper, StructureError, contact_map
from histo_contacts.arpeggio_backend import INTERACTION_TYPES
from histo_contacts.selectors import SelectorError

FIXTURE = Path(__file__).parent / "fixtures" / "8gvi_1_aligned.pdb"

ROW_KEYS = {
    "type",
    "from_chain",
    "from_residue",
    "from_atom",
    "from_aa",
    "to_chain",
    "to_residue",
    "to_atom",
    "to_aa",
    "distance",
    "bond_types",
}


def test_contact_map_row_shape_and_distance_bound():
    cm = ContactMapper(FIXTURE)
    rows = cm.contact_map("P", ["H"], distance=5.0)
    assert rows, "expected at least one contact between peptide P and MHC chain H"
    for row in rows:
        assert set(row) == ROW_KEYS
        assert row["from_chain"] == "P"
        assert row["to_chain"] == "H"
        assert 0 < row["distance"] <= 5.0
        assert row["type"] in INTERACTION_TYPES
        assert row["bond_types"], "every row should carry at least one Arpeggio bond type"


def test_contact_map_sorted_by_target_chain_then_residue():
    cm = ContactMapper(FIXTURE)
    rows = cm.contact_map("P", ["H", "L"], distance=5.0)
    assert rows

    chains_seen = [r["to_chain"] for r in rows]
    assert chains_seen == sorted(chains_seen)

    for chain_id in set(chains_seen):
        residues = [r["to_residue"] for r in rows if r["to_chain"] == chain_id]
        assert residues == sorted(residues)


def test_contact_map_single_residue_query():
    cm = ContactMapper(FIXTURE)
    rows = cm.contact_map("P:5", ["H", "L"], distance=5.0)
    assert rows
    assert {r["from_residue"] for r in rows} == {5}
    assert {r["from_chain"] for r in rows} == {"P"}


def test_contact_map_residue_range_query():
    cm = ContactMapper(FIXTURE)
    rows = cm.contact_map("P:1-3", ["H", "L"], distance=5.0)
    assert rows
    assert {r["from_residue"] for r in rows} <= {1, 2, 3}


def test_contact_map_distance_cutoff_monotonic():
    cm = ContactMapper(FIXTURE)
    small = cm.contact_map("P:1-3", ["H", "L"], distance=3.0)
    large = cm.contact_map("P:1-3", ["H", "L"], distance=6.0)
    assert len(large) >= len(small)


def test_contact_map_rejects_non_positive_distance():
    cm = ContactMapper(FIXTURE)
    with pytest.raises(ValueError):
        cm.contact_map("P", ["H"], distance=0)


def test_contact_map_unknown_target_chain():
    cm = ContactMapper(FIXTURE)
    with pytest.raises(StructureError):
        cm.contact_map("P", ["Z"], distance=5.0)


def test_contact_map_unknown_query_chain():
    cm = ContactMapper(FIXTURE)
    with pytest.raises(StructureError):
        cm.contact_map("Z", ["H"], distance=5.0)


def test_contact_map_bad_selector_syntax():
    cm = ContactMapper(FIXTURE)
    with pytest.raises(SelectorError):
        cm.contact_map("P:9-1", ["H"], distance=5.0)


def test_write_json_roundtrip(tmp_path):
    cm = ContactMapper(FIXTURE)
    rows = cm.contact_map("P:1-3", ["H"], distance=5.0)
    out = tmp_path / "contacts.json"
    cm.write_json(rows, out)

    loaded = json.loads(out.read_text())
    assert loaded == rows


def test_module_level_convenience_function():
    rows = contact_map(FIXTURE, "P:1-3", ["H"], distance=5.0)
    assert rows
    assert set(rows[0]) == ROW_KEYS


def test_target_accepts_string_or_iterable():
    cm = ContactMapper(FIXTURE)
    rows_str = cm.contact_map("P:1-3", "H,L", distance=5.0)
    rows_list = cm.contact_map("P:1-3", ["H", "L"], distance=5.0)
    assert rows_str == rows_list


def test_repeated_calls_on_same_mapper_do_not_accumulate():
    cm = ContactMapper(FIXTURE)
    first = cm.contact_map("P", ["H"], distance=5.0)
    cm.contact_map("P:5", ["H", "L"], distance=5.0)  # a different query in between
    third = cm.contact_map("P", ["H"], distance=5.0)
    assert first == third
