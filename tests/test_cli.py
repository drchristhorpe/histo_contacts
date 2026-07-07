import json
from pathlib import Path

from click.testing import CliRunner

from histo_contacts.cli import main
from histo_contacts.core import contact_map

FIXTURE = Path(__file__).parent / "fixtures" / "8gvi_1_aligned.pdb"


def test_cli_writes_json_matching_library(tmp_path):
    output = tmp_path / "contacts.json"
    runner = CliRunner()
    result = runner.invoke(
        main,
        [str(FIXTURE), "--query", "P:1-3", "--target", "H,L", "--output", str(output)],
    )
    assert result.exit_code == 0, result.output
    assert output.exists()

    rows = json.loads(output.read_text())
    expected = contact_map(FIXTURE, "P:1-3", ["H", "L"], distance=5.0)
    assert rows == expected
    assert str(len(rows)) in result.output


def test_cli_respects_distance_option(tmp_path):
    output = tmp_path / "contacts.json"
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            str(FIXTURE),
            "--query",
            "P:1-3",
            "--target",
            "H",
            "--output",
            str(output),
            "--distance",
            "3.0",
        ],
    )
    assert result.exit_code == 0, result.output
    rows = json.loads(output.read_text())
    assert all(row["distance"] <= 3.0 for row in rows)


def test_cli_missing_required_option_errors():
    runner = CliRunner()
    result = runner.invoke(main, [str(FIXTURE), "--target", "H", "--output", "out.json"])
    assert result.exit_code != 0
    assert "query" in result.output.lower()


def test_cli_unknown_chain_errors_gracefully(tmp_path):
    output = tmp_path / "contacts.json"
    runner = CliRunner()
    result = runner.invoke(
        main,
        [str(FIXTURE), "--query", "P:1-3", "--target", "ZZ", "--output", str(output)],
    )
    assert result.exit_code != 0
    assert not output.exists()
