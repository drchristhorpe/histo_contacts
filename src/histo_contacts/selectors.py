"""Parsing of query/target selector strings.

Same grammar as histo_neighbours.selectors, reimplemented locally per house
convention (sibling tools are independent packages, not a shared dependency).
See docs/PLAN.md (section 4) for the selector grammar.
"""

from __future__ import annotations

from dataclasses import dataclass


class SelectorError(ValueError):
    """Raised for malformed or ambiguous query/target selectors."""


@dataclass(frozen=True)
class QueryPart:
    """A chain, or a residue range within a chain — one token of a query."""

    chain: str | None
    start: int | None = None
    end: int | None = None


def _split_chain(token: str) -> tuple[str | None, str]:
    """Split ``"A:1-180"`` into ``("A", "1-180")``; ``"1-180"`` into ``(None, "1-180")``."""
    if ":" in token:
        chain, _, rest = token.partition(":")
        chain = chain.strip()
        if not chain:
            raise SelectorError(f"Empty chain id in selector token {token!r}")
        return chain, rest.strip()
    return None, token.strip()


def _is_range(text: str) -> bool:
    if "-" not in text:
        return False
    left, _, right = text.partition("-")
    return left.strip().lstrip("+").isdigit() and right.strip().isdigit()


def _tokenize(spec: str) -> list[str]:
    tokens = [t.strip() for t in spec.split(",")]
    tokens = [t for t in tokens if t]
    if not tokens:
        raise SelectorError("Empty selector")
    return tokens


def _parse_query_token(token: str) -> QueryPart:
    """Parse one query token: a bare chain (``"A"``), a bare residue range
    (``"1-180"``), a chain-scoped range (``"A:1-180"``), or a chain-scoped
    single residue (``"A:12"``)."""
    text = token.strip()
    if not text:
        raise SelectorError("Empty query token")

    chain, rest = _split_chain(text)

    if not rest:
        if chain is None:
            raise SelectorError(f"Empty query token {token!r}")
        return QueryPart(chain=chain)

    if _is_range(rest):
        start_s, _, end_s = rest.partition("-")
        start, end = int(start_s), int(end_s)
        if start > end:
            raise SelectorError(f"Invalid residue range {rest!r}: start > end")
        return QueryPart(chain=chain, start=start, end=end)

    if rest.lstrip("+").isdigit():
        resnum = int(rest)
        return QueryPart(chain=chain, start=resnum, end=resnum)

    if chain is None:
        # A bare alphabetic token with no ':' is a chain letter, e.g. "A".
        return QueryPart(chain=text)

    raise SelectorError(f"Could not parse query token {token!r}")


def parse_query(spec: str) -> list[QueryPart]:
    """Parse a query ("from") selector into its constituent parts.

    A comma-separated string of tokens, each a chain, chain+range, or
    chain+single-residue (see docs/PLAN.md section 4). Every part's atoms
    are unioned into a single "from" atom set by the caller.
    """
    if not isinstance(spec, str):
        raise SelectorError(f"Query selector must be a string, got {type(spec).__name__}")
    return [_parse_query_token(t) for t in _tokenize(spec)]


def parse_target(spec: str) -> list[str]:
    """Parse a target ("to") selector: a comma-separated list of chain ids."""
    if not isinstance(spec, str):
        raise SelectorError(f"Target selector must be a string, got {type(spec).__name__}")
    tokens = _tokenize(spec)
    for t in tokens:
        if ":" in t or _is_range(t):
            raise SelectorError(
                f"Target selector takes chain ids only, not residues/ranges: {t!r}"
            )
    return tokens
