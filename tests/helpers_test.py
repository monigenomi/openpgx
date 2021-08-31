from src.openpgx.dpwg import get_dpwg_recommendations
from src.openpgx.helpers import *


def test_index_table_by_something():
    result = index_items_by_key(
        [
            {"Drug": "Tamsulosin", "Gene": "CYP2D6"},
            {"Drug": "Trimipramine", "Gene": "CYP2D6"},
        ],
        "Gene",
    )

    expected = {
        "CYP2D6": [
            {"Drug": "Tamsulosin", "Gene": "CYP2D6"},
            {"Drug": "Trimipramine", "Gene": "CYP2D6"},
        ]
    }

    assert result == expected


def test_get_phenoconversion_data_from_recommendations():
    assert (
        len(get_phenoconversion_data_from_recommendations(get_dpwg_recommendations()))
        > 0
    )


def test_normalize_gene_and_factor():
    assert normalize_hla_gene_and_factor("HLA-B", "*57:01 positive") == (
        "HLA-B*57:01",
        "positive",
    )
    assert normalize_hla_gene_and_factor("HLA-A*31:01", "*31:01 positive") == (
        "HLA-A*31:01",
        "positive",
    )


def test_words_to_sentence():
    assert words_to_sentence(["foo"]) == "foo"
    assert words_to_sentence(["foo", "bar"]) == "foo and bar"
    assert words_to_sentence(["foo", "bar", "baz"]) == "foo, bar and baz"


def test_format_with_populations():
    assert format_with_populations(
        {
            "adults": {"recommendation": "hello", "foo": "bar"},
            "CBZ naive": {"recommendation": "world", "foo": "bar"},
        }
    ) == {
        "recommendation": "Adults: hello\n\nIf patient has not previously used carbamazepine: world",
        "foo": "bar",
    }
