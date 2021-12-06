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

def test_simple_sql_parse():
    file = open(os.path.join(cwd, "fixtures/cpic_inserts.sql"), 'r')
    database = defaultdict(list)
    for query in yield_inserts_from_file(file):
        table_name, row = sql_to_data(query)
        database[table_name].append(row)
    assert dict(database) == {
        'allele_definition': [
            {'genesymbol': 'HLA-A',
             'id': 9000,
             'name': '*31:01',
             'pharmvarid': None,
             'reference': False,
             'structuralvariation': False,
             'version': 1},
            {'genesymbol': 'HLA-B',
             'id': 9001,
             'name': '*31:02',
             'pharmvarid': None,
             'reference': True,
             'structuralvariation': True,
             'version': 2}
        ],
        'gene': [
            {'chr': 'chrM',
             'chromosequenceid': 'NC_012920.1',
             'ensemblid': 'ENSG00000211459',
             'frequencymethods': 'Methods\nMultiline',
             'functionmethods': None,
             'genesequenceid': 'NG_042193.1',
             'hgncid': 'HGNC:7470',
             'lookupmethod': 'PHENOTYPE',
             'mrnasequenceid': 'NM_001042351.2',
             'ncbiid': '4549',
             'notesondiplotype': None,
             'pharmgkbid': 'PA31274',
             'proteinsequenceid': 'NP_000101.2',
             'symbol': 'MT-RNR1',
             'url': 'https://cpicpgx.org/gene/mt-rnr1/',
             'version': 43
             }
        ]
    }

