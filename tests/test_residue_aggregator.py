"""Tests for residue_aggregator: residue-level contact aggregation."""

import json
from pathlib import Path
import pytest

from histo_contacts import residue_aggregator


class TestCDRAndMHCDefinitions:
    """Test CDR loop and MHC region definitions."""

    def test_cdr_loops_structure(self):
        loops = residue_aggregator.define_cdr_loops()
        assert "alpha" in loops
        assert "beta" in loops
        assert loops["alpha"]["cdr1"] == (27, 38)
        assert loops["alpha"]["cdr2"] == (56, 65)
        assert loops["alpha"]["cdr3"] == (105, 117)
        assert loops["beta"]["cdr1"] == (27, 38)
        assert loops["beta"]["cdr3"] == (105, 117)

    def test_mhc_regions_structure(self):
        regions = residue_aggregator.define_mhc_regions()
        assert regions["alpha1"] == ("A", 50, 86)
        assert regions["alpha2"] == ("A", 137, 180)
        assert regions["peptide"] == ("C", 1, 1000)


class TestFindCDRLoop:
    """Test _find_cdr_loop internal function."""

    def test_resnum_in_cdr1(self):
        loops = residue_aggregator.define_cdr_loops()
        result = residue_aggregator._find_cdr_loop(30, loops["alpha"])
        assert result == "cdr1"

    def test_resnum_in_cdr3(self):
        loops = residue_aggregator.define_cdr_loops()
        result = residue_aggregator._find_cdr_loop(110, loops["alpha"])
        assert result == "cdr3"

    def test_resnum_outside_any_loop(self):
        loops = residue_aggregator.define_cdr_loops()
        result = residue_aggregator._find_cdr_loop(50, loops["alpha"])
        assert result is None


class TestFindMHCRegion:
    """Test _find_mhc_region internal function."""

    def test_resnum_in_alpha1(self):
        regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator._find_mhc_region(60, "A", regions)
        assert result == "alpha1"

    def test_resnum_in_alpha2(self):
        regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator._find_mhc_region(160, "A", regions)
        assert result == "alpha2"

    def test_chain_c_is_peptide(self):
        regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator._find_mhc_region(5, "C", regions)
        assert result == "peptide"

    def test_wrong_chain(self):
        regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator._find_mhc_region(50, "B", regions)
        assert result is None

    def test_alpha1_boundary_inclusive(self):
        regions = residue_aggregator.define_mhc_regions()
        assert residue_aggregator._find_mhc_region(50, "A", regions) == "alpha1"
        assert residue_aggregator._find_mhc_region(86, "A", regions) == "alpha1"

    def test_gap_between_alpha1_alpha2(self):
        regions = residue_aggregator.define_mhc_regions()
        # 87-136 is the gap (β-sheet floor), excluded
        assert residue_aggregator._find_mhc_region(87, "A", regions) is None
        assert residue_aggregator._find_mhc_region(136, "A", regions) is None


class TestCategorizeContact:
    """Test contact categorization into (CDR, region) pairs."""

    def test_tcr_alpha_to_mhc_alpha1(self):
        """TCR CDR1 on chain D to MHC α1 on chain A."""
        atom_row = {
            "from_chain": "D",
            "from_residue": 30,
            "to_chain": "A",
            "to_residue": 60,
        }
        cdr_loops = residue_aggregator.define_cdr_loops()
        mhc_regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator.categorize_contact(atom_row, cdr_loops, mhc_regions)
        assert result == ("D", "cdr1_alpha", "alpha1", "from")

    def test_tcr_beta_to_peptide(self):
        """TCR CDR3 on chain E to peptide on chain C."""
        atom_row = {
            "from_chain": "E",
            "from_residue": 110,
            "to_chain": "C",
            "to_residue": 5,
        }
        cdr_loops = residue_aggregator.define_cdr_loops()
        mhc_regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator.categorize_contact(atom_row, cdr_loops, mhc_regions)
        assert result == ("E", "cdr3_beta", "peptide", "from")

    def test_mhc_to_tcr_reversed(self):
        """MHC α1 to TCR CDR1 (reversed sides)."""
        atom_row = {
            "from_chain": "A",
            "from_residue": 70,
            "to_chain": "D",
            "to_residue": 32,
        }
        cdr_loops = residue_aggregator.define_cdr_loops()
        mhc_regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator.categorize_contact(atom_row, cdr_loops, mhc_regions)
        assert result == ("D", "cdr1_alpha", "alpha1", "to")

    def test_out_of_scope(self):
        """Contact outside CDR and MHC regions."""
        atom_row = {
            "from_chain": "B",
            "from_residue": 100,
            "to_chain": "L",
            "to_residue": 50,
        }
        cdr_loops = residue_aggregator.define_cdr_loops()
        mhc_regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator.categorize_contact(atom_row, cdr_loops, mhc_regions)
        assert result is None

    def test_tcr_outside_cdr_loops(self):
        """TCR on D but outside any CDR loop."""
        atom_row = {
            "from_chain": "D",
            "from_residue": 50,  # Between CDR1 and CDR2
            "to_chain": "A",
            "to_residue": 60,
        }
        cdr_loops = residue_aggregator.define_cdr_loops()
        mhc_regions = residue_aggregator.define_mhc_regions()
        result = residue_aggregator.categorize_contact(atom_row, cdr_loops, mhc_regions)
        assert result is None


class TestAggregateContacts:
    """Test residue-level aggregation."""

    def test_single_atom_pair(self):
        """Aggregate a single atom pair."""
        atom_rows = [
            {
                "from_chain": "D",
                "from_residue": 110,
                "from_atom": "OG",
                "from_aa": "THR",
                "to_chain": "C",
                "to_residue": 5,
                "to_atom": "CB",
                "to_aa": "TYR",
                "bond_types": ["vdw"],
                "distance": 3.5,
            }
        ]
        residue_rows = residue_aggregator.aggregate_contacts(atom_rows)
        assert len(residue_rows) == 1
        row = residue_rows[0]
        assert row["cdr_loop"] == "cdr3_alpha"
        assert row["region"] == "peptide"
        assert row["tcr_chain"] == "D"
        assert row["tcr_resnum"] == 110
        assert row["mhc_chain"] == "C"
        assert row["mhc_resnum"] == 5
        assert row["n_atom_pairs"] == 1
        assert row["bond_type_counts"] == {"vdw": 1}
        assert row["min_distance"] == 3.5

    def test_multiple_atom_pairs_same_residue_pair(self):
        """Multiple atom pairs connecting the same two residues."""
        atom_rows = [
            {
                "from_chain": "D",
                "from_residue": 30,
                "from_atom": "CB",
                "from_aa": "ALA",
                "to_chain": "A",
                "to_residue": 60,
                "to_atom": "OE",
                "to_aa": "GLU",
                "bond_types": ["hbond", "polar"],
                "distance": 3.0,
            },
            {
                "from_chain": "D",
                "from_residue": 30,
                "from_atom": "N",
                "from_aa": "ALA",
                "to_chain": "A",
                "to_residue": 60,
                "to_atom": "O",
                "to_aa": "GLU",
                "bond_types": ["vdw"],
                "distance": 3.8,
            },
        ]
        residue_rows = residue_aggregator.aggregate_contacts(atom_rows)
        assert len(residue_rows) == 1
        row = residue_rows[0]
        assert row["n_atom_pairs"] == 2
        # hbond: 1, polar: 1, vdw: 1 (one pair has multiple types)
        assert row["bond_type_counts"] == {"hbond": 1, "polar": 1, "vdw": 1}
        assert row["min_distance"] == 3.0

    def test_overlapping_bond_types(self):
        """Verify bond types overlap correctly (counts don't sum to n_atom_pairs)."""
        atom_rows = [
            {
                "from_chain": "D",
                "from_residue": 30,
                "from_atom": "CB",
                "from_aa": "LEU",
                "to_chain": "A",
                "to_residue": 60,
                "to_atom": "CD",
                "to_aa": "LEU",
                "bond_types": ["proximal", "vdw"],
                "distance": 4.0,
            },
            {
                "from_chain": "D",
                "from_residue": 30,
                "from_atom": "CD",
                "from_aa": "LEU",
                "to_chain": "A",
                "to_residue": 60,
                "to_atom": "CE",
                "to_aa": "LEU",
                "bond_types": ["proximal"],
                "distance": 4.2,
            },
        ]
        residue_rows = residue_aggregator.aggregate_contacts(atom_rows)
        row = residue_rows[0]
        assert row["n_atom_pairs"] == 2
        # proximal: 2 (both pairs), vdw: 1 (first pair only)
        assert row["bond_type_counts"] == {"proximal": 2, "vdw": 1}
        assert sum(row["bond_type_counts"].values()) == 3  # NOT equal to n_atom_pairs

    def test_sorting_by_cdr_loop_region_residue(self):
        """Verify output is sorted by CDR, region, then residue numbers."""
        atom_rows = [
            # CDR3_alpha, peptide
            {
                "from_chain": "D",
                "from_residue": 110,
                "from_atom": "CA",
                "from_aa": "GLY",
                "to_chain": "C",
                "to_residue": 2,
                "to_atom": "CA",
                "to_aa": "ALA",
                "bond_types": ["vdw"],
                "distance": 4.0,
            },
            # CDR1_alpha, alpha1
            {
                "from_chain": "D",
                "from_residue": 30,
                "from_atom": "CA",
                "from_aa": "GLY",
                "to_chain": "A",
                "to_residue": 60,
                "to_atom": "CA",
                "to_aa": "ALA",
                "bond_types": ["vdw"],
                "distance": 4.0,
            },
        ]
        residue_rows = residue_aggregator.aggregate_contacts(atom_rows)
        assert len(residue_rows) == 2
        assert residue_rows[0]["cdr_loop"] == "cdr1_alpha"
        assert residue_rows[1]["cdr_loop"] == "cdr3_alpha"


class TestReconciliation:
    """Test acceptance test: reconciliation against expected totals."""

    def test_reconcile_pass(self):
        """Reconciliation passes when counts match."""
        residue_rows = [
            {
                "cdr_loop": "cdr1_alpha",
                "region": "alpha1",
                "tcr_resnum": 30,
                "mhc_resnum": 60,
                "n_atom_pairs": 5,
                "bond_type_counts": {"proximal": 5},
                "min_distance": 3.5,
                "tcr_chain": "D",
                "tcr_resname": "ALA",
                "mhc_chain": "A",
                "mhc_resname": "GLU",
            },
            {
                "cdr_loop": "cdr1_alpha",
                "region": "alpha1",
                "tcr_resnum": 31,
                "mhc_resnum": 60,
                "n_atom_pairs": 3,
                "bond_type_counts": {"vdw": 3},
                "min_distance": 4.0,
                "tcr_chain": "D",
                "tcr_resname": "GLY",
                "mhc_chain": "A",
                "mhc_resname": "GLU",
            },
        ]
        expected = {"cdr1_alpha|alpha1": 8}  # 5 proximal + 3 vdw
        assert residue_aggregator.reconcile_against_aggregates(residue_rows, expected)

    def test_reconcile_fail_mismatch(self):
        """Reconciliation fails when counts don't match."""
        residue_rows = [
            {
                "cdr_loop": "cdr1_alpha",
                "region": "alpha1",
                "tcr_resnum": 30,
                "mhc_resnum": 60,
                "bond_type_counts": {"proximal": 5},
                "tcr_chain": "D",
                "tcr_resname": "ALA",
                "mhc_chain": "A",
                "mhc_resname": "GLU",
            },
        ]
        expected = {"cdr1_alpha|alpha1": 10}  # Expected 10, got 5
        assert not residue_aggregator.reconcile_against_aggregates(residue_rows, expected)


class TestNoRegression:
    """Test no-regression check: helix-only <= whole-domain."""

    def test_no_regression_pass(self):
        """No regression when helix-only <= whole-domain."""
        new_counts = {"CDR1_alpha|alpha1": 20, "CDR1_alpha|alpha2": 30}
        old_counts = {"CDR1_alpha|alpha1": 25, "CDR1_alpha|alpha2": 40}
        assert residue_aggregator.validate_no_regression(new_counts, old_counts)

    def test_no_regression_equal(self):
        """No regression when helix-only == whole-domain."""
        new_counts = {"CDR1_alpha|alpha1": 25}
        old_counts = {"CDR1_alpha|alpha1": 25}
        assert residue_aggregator.validate_no_regression(new_counts, old_counts)

    def test_regression_detected(self):
        """Regression detected when helix-only > whole-domain."""
        new_counts = {"CDR1_alpha|alpha1": 30}
        old_counts = {"CDR1_alpha|alpha1": 25}
        assert not residue_aggregator.validate_no_regression(new_counts, old_counts)


class TestAcceptanceTest1AO7:
    """Acceptance test with sample 1AO7 data."""

    @pytest.fixture
    def sample_atoms(self):
        """Load 1AO7 sample atom-level contacts from the brief."""
        sample_dir = Path(__file__).parent / "fixtures" / ".." / ".." / "briefs" / "residue_contacts" / "sample"
        totals_path = sample_dir / "expected_totals__1AO7_aligned_1.json"
        if not totals_path.exists():
            pytest.skip(f"Sample data not found at {totals_path}")
        with open(totals_path) as f:
            expected = json.load(f)
        return expected

    @pytest.fixture
    def sample_structure(self):
        """Load sample 1AO7 PDB for contact extraction."""
        sample_dir = Path(__file__).parent / "fixtures" / ".." / ".." / "briefs" / "residue_contacts" / "sample"
        pdb_path = sample_dir / "1ao7_aligned_1.pdb"
        if not pdb_path.exists():
            pytest.skip(f"Sample PDB not found at {pdb_path}")
        return pdb_path

    def test_aggregate_1ao7_sample(self, sample_atoms, sample_structure):
        """Test aggregation on real 1AO7 sample data."""
        from histo_contacts import ContactMapper

        # Generate contacts using histo_contacts
        cm = ContactMapper(str(sample_structure))
        atom_rows = cm.contact_map(
            "D:27-38,D:56-65,D:105-117,E:27-38,E:56-65,E:105-117",
            ["A", "C"],
            distance=5.0
        )

        # Aggregate to residue level
        residue_rows = residue_aggregator.aggregate_contacts(atom_rows)

        # Compute aggregates
        cell_totals = {}
        for row in residue_rows:
            cell = f"{row['cdr_loop']}|{row['region']}"
            cell_totals[cell] = cell_totals.get(cell, 0) + row["n_atom_pairs"]

        # Reconcile against expected
        expected_counts = sample_atoms["ct_atom_pairs"]
        # Normalize expected keys (brief uses mixed cases)
        expected_norm = {k.upper(): v for k, v in expected_counts.items()}
        computed_norm = {k.upper(): v for k, v in cell_totals.items()}

        assert computed_norm == expected_norm, f"Mismatch: expected {expected_norm}, got {computed_norm}"

        # Check total atom pairs
        total_computed = sum(row["n_atom_pairs"] for row in residue_rows)
        total_expected = sample_atoms["ct_total_atom_pairs"]
        assert total_computed == total_expected, f"Total mismatch: expected {total_expected}, got {total_computed}"
