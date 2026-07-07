import pytest

from histo_contacts.selectors import QueryPart, SelectorError, parse_query, parse_target


def test_parse_query_bare_chain():
    assert parse_query("A") == [QueryPart(chain="A")]


def test_parse_query_chain_residue():
    assert parse_query("A:12") == [QueryPart(chain="A", start=12, end=12)]


def test_parse_query_chain_range():
    assert parse_query("A:1-180") == [QueryPart(chain="A", start=1, end=180)]


def test_parse_query_multiple_tokens_union():
    parts = parse_query("A:1-9,B:12")
    assert parts == [
        QueryPart(chain="A", start=1, end=9),
        QueryPart(chain="B", start=12, end=12),
    ]


def test_parse_query_bare_residue_no_chain():
    assert parse_query("1-9") == [QueryPart(chain=None, start=1, end=9)]
    assert parse_query("12") == [QueryPart(chain=None, start=12, end=12)]


def test_parse_query_invalid_range():
    with pytest.raises(SelectorError):
        parse_query("A:9-1")


def test_parse_query_empty():
    with pytest.raises(SelectorError):
        parse_query("")


def test_parse_query_non_string():
    with pytest.raises(SelectorError):
        parse_query(["A"])


def test_parse_target_list_of_chains():
    assert parse_target("H,L") == ["H", "L"]


def test_parse_target_single_chain():
    assert parse_target("H") == ["H"]


def test_parse_target_rejects_residue_tokens():
    with pytest.raises(SelectorError):
        parse_target("H:1-10")


def test_parse_target_rejects_range_token():
    with pytest.raises(SelectorError):
        parse_target("1-10")


def test_parse_target_empty():
    with pytest.raises(SelectorError):
        parse_target("")
