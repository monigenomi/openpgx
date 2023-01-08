from openpgx.dpwg import create_dpwg_database
from openpgx.helpers import *

cwd = os.path.dirname(os.path.realpath(__file__))


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

# TODO - check this
# def test_format_with_populations():
#     assert format_with_populations(
#         {
#             "adults": {"recommendation": "hello", "foo": "bar"},
#             "CBZ naive": {"recommendation": "world", "foo": "bar"},
#         }
#     ) == {
#         "recommendation": "Adults: hello\n\nIf patient has not previously used carbamazepine: world",
#         "foo": "bar",
#     }


def test_url_to_dir():
    assert (
        url_to_cache_dir(
            "https://github.com/cpicpgx/cpic-data/releases/download/v1.10/inserts.sql.gz"
        )
        == "github.com/cpicpgx/cpic-data/releases/download/v1.10"
    )
    assert (
        url_to_cache_dir("https://github.com/cpicpgx/example.zip")
        == "github.com/cpicpgx/example"
    )


def test_extract_usage():
    usage = extract_usage(repository_path("README.md")).split("\n")
    assert usage[0][0:9] == "$ openpgx"
    assert usage[-1] == "https://github.com/monigenomi/openpgx"
