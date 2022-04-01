import os

from openpgx.dpwg import get_dpwg_recommendations
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


def test_get_normalizations():
    print(get_dpwg_recommendations())
    assert len(get_normalizations(get_dpwg_recommendations())) != 0


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


def test_import_parse_copy():
    # TODO:
    # 1. parse first line with regexp to extract table name and column names
    # 2. parse next lines up until \. to extract data. use tsv reader
    #    https://www.pythonpool.com/read-tsv-file-python/
    
    input = "COPY cpic.allele (id, version, genesymbol, name, functionalstatus, clinicalfunctionalstatus, clinicalfunctionalsubstrate, activityvalue, definitionid, citations, strength, functioncomments, findings, frequency) FROM stdin;"
    
    table, columns = parse_copy(input)
    
    assert table == "allele"
    assert columns == ['id', 'version', 'genesymbol', 'name', 'functionalstatus', 'clinicalfunctionalstatus',
                       'clinicalfunctionalsubstrate', 'activityvalue', 'definitionid', 'citations', 'strength',
                       'functioncomments', 'findings', 'frequency']


def test_yield_rows_from_sql_file():
    database = defaultdict(list)
    for table, row in yield_rows_from_sql_file(os.path.join(cwd, "fixtures/cpic.sql")):
        database[table].append(row)
    assert dict(database) == {
        'allele': [{
            'activityvalue': None,
            'citations': '{}',
            'clinicalfunctionalstatus': 'Normal Function',
            'clinicalfunctionalsubstrate': None,
            'definitionid': '777262',
            'findings': None,
            'frequency': None,
            'functionalstatus': None,
            'functioncomments': None,
            'genesymbol': 'CACNA1S',
            'id': '777263',
            'name': 'Reference',
            'strength': None,
            'version': '25'
        },
            {
                'activityvalue': None,
                'citations': '{22232210}',
                'clinicalfunctionalstatus': 'No function',
                'clinicalfunctionalsubstrate': None,
                'definitionid': '1357093',
                'findings': 'SLCO1B1*48 is assigned no function due to evidence '
                            'supporting a partial gene deletion (22232210). '
                            'Therefore, consensus among experts was no function '
                            'due to limited evidence.',
                'frequency': None,
                'functionalstatus': None,
                'functioncomments': None,
                'genesymbol': 'SLCO1B1',
                'id': '1357094',
                'name': '*49',
                'strength': 'Limited',
                'version': '9'
            }],
        'allele_definition': [{
            'genesymbol': 'HLA-A',
            'id': '9000',
            'name': '*31:01',
            'pharmvarid': None,
            'reference': False,
            'structuralvariation': False,
            'version': '1'
        },
            {
                'genesymbol': 'HLA-B',
                'id': '9001',
                'name': '*15:02',
                'pharmvarid': None,
                'reference': False,
                'structuralvariation': False,
                'version': '1'
            }]
    }


def test_url_to_dir():
    assert url_to_cache_dir("https://github.com/cpicpgx/cpic-data/releases/download/v1.10/inserts.sql.gz") == \
           "github.com/cpicpgx/cpic-data/releases/download/v1.10"
    assert url_to_cache_dir("https://github.com/cpicpgx/example.zip") == \
           "github.com/cpicpgx/example"


def test_extract_usage():
    usage = extract_usage(repository_path('README.md')).split("\n")
    assert usage[0][0:9] == "$ openpgx"
    assert usage[-1] == "https://github.com/monigenomi/openpgx"
