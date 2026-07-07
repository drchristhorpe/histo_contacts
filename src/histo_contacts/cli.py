"""Command line interface for histo_contacts."""

from __future__ import annotations

import click
from rich.console import Console
from rich.table import Table

from histo_contacts.core import ContactMapper, StructureError
from histo_contacts.selectors import SelectorError

_DISPLAY_CAP = 25


@click.command()
@click.argument("filename", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--query",
    "-q",
    required=True,
    help="Query selector (the 'from' side): comma-separated chains and/or "
    "chain:residue / chain:range tokens, e.g. 'P', 'A:1-9', 'A:1-9,B:12'.",
)
@click.option(
    "--target",
    "-t",
    required=True,
    help="Target selector (the 'to' side): comma-separated list of chain ids, e.g. 'H,L'.",
)
@click.option(
    "--output",
    "-o",
    required=True,
    type=click.Path(dir_okay=False, writable=True),
    help="Path to write the contact map as a JSON file.",
)
@click.option(
    "--distance",
    "-d",
    type=float,
    default=5.0,
    show_default=True,
    help="Contact distance cutoff in Angstroms (forwarded to Arpeggio's interacting-distance cutoff).",
)
def main(filename: str, query: str, target: str, output: str, distance: float) -> None:
    """Compute a bond-typed contact map, using PDBe Arpeggio, between a
    query selection and one or more target chains in a 3D biological
    structure (PDB/mmCIF).

    FILENAME is the path to a .cif/.mmcif or .pdb/.ent structure file.
    """
    console = Console()

    try:
        mapper = ContactMapper(filename)
        rows = mapper.contact_map(query, target, distance=distance)
        mapper.write_json(rows, output)
    except (StructureError, SelectorError) as exc:
        raise click.ClickException(str(exc)) from exc

    table = Table(title=f"Contacts: {query!r} -> {target!r} (≤ {distance}Å)")
    for column in ("from_chain", "from_residue", "from_atom", "to_chain", "to_residue", "to_atom", "distance", "bond_types"):
        table.add_column(column)

    for row in rows[:_DISPLAY_CAP]:
        table.add_row(
            row["from_chain"],
            str(row["from_residue"]),
            row["from_atom"],
            row["to_chain"],
            str(row["to_residue"]),
            row["to_atom"],
            f"{row['distance']:.2f}",
            ", ".join(row["bond_types"]),
        )

    console.print(table)
    if len(rows) > _DISPLAY_CAP:
        console.print(f"... and {len(rows) - _DISPLAY_CAP} more rows (see {output})")

    console.print(f"[bold]{len(rows)}[/bold] contact(s) within {distance}Å written to {output}")


if __name__ == "__main__":
    main()
