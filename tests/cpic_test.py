#!/usr/bin/env python3
import os

from openpgx.cpic import *
from openpgx.helpers import save_json

CPIC_DATABASE = create_cpic_database()
CPIC_RECOMMENDATIONS = CPIC_DATABASE["recommendations"]

cwd = os.path.dirname(os.path.realpath(__file__))

cached_sql_gz = download_to_cache_dir(CPIC_DEFAULT_URL)

DATA = load_cpic_dump(cached_sql_gz)


def test_get_alleles():
    assert get_alleles([
        {
            'id': 828130, 'version': 11, 'genesymbol': 'MT-RNR1', 'name': '1095T>C', 'functionalstatus': None,
            'clinicalfunctionalstatus': 'increased risk of aminoglycoside-induced hearing loss',
            'clinicalfunctionalsubstrate': None, 'activityvalue': 'n/a', 'definitionid': 828129,
            'citations': '{11079536,11313749,15555598,15841390,16875663,21205314}', 'strength': 'Moderate',
            'functioncomments': None,
            'findings': 'MT-RNR1 1095T>C is assigned an increased risk of aminoglycoside-induced hearing loss based on moderate evidence. Screening of individuals or probands of individuals with aminoglycoside-hearing looss identified 1095T>C in a number of cases (11079536, 11313749, 15555598, 15841390, 16875663, 21205314). Therefore, consensus among experts was this variant is associated with increased risk of aminoglycoside-induced hearing loss based on moderate evidence.'
        },
        {
            'id': 777350, 'version': 11, 'genesymbol': 'CFTR', 'name': 'G551D', 'functionalstatus': None,
            'clinicalfunctionalstatus': 'ivacaftor responsive', 'clinicalfunctionalsubstrate': None,
            'activityvalue': None, 'definitionid': 777349, 'citations': '{}', 'strength': None,
            'functioncomments': 'ACMG screening panel', 'findings': None
        },
        {
            'id': 777316, 'version': 11, 'genesymbol': 'CFTR', 'name': '711+3A->G', 'functionalstatus': None,
            'clinicalfunctionalstatus': 'ivacaftor responsive', 'clinicalfunctionalsubstrate': None,
            'activityvalue': None, 'definitionid': 777315, 'citations': '{}', 'strength': None,
            'functioncomments': 'ACMG screening panel', 'findings': None
        }
    
    ]) == {
               "CFTR": [{"allele": "G551D"}, {"allele": "711+3A->G"}],
               "MT-RNR1": [{"allele": "1095T>C"}]
           }


def test_normalize_cpic_factor():
    assert normalize_cpic_factor("HLA-A", "No *12 Result") == ("HLA-A*12", None)
    assert normalize_cpic_factor("AAA", "No Result") == ("AAA", None)
    assert normalize_cpic_factor("FOO", "n/a") == ("FOO", None)


def test_normalize_genename():
    assert normalize_genename("CYP2D6", "*1/*1") == "CYP2D6"
    assert normalize_genename("HLA-A", "*15:02 positive") == "HLA-A*15:02"


def test_normalize_activityscore():
    # For factors
    assert normalize_activityscore("1", True) == "== 1.00"
    assert normalize_activityscore("n/a", True) == None
    assert normalize_activityscore("≥4", True) == ">= 4.00"
    # For phenotype table
    assert normalize_activityscore("1", False) == 1.00
    assert normalize_activityscore("1.00", False) == 1.00
    assert normalize_activityscore("> 1.5", False) == 1.50
    assert normalize_activityscore("≥4", False) == 4.00


def test_create_phenotype_and_activityscore_table():
    activityscore, phenotype = create_phenotype_and_activityscore_table(DATA["gene_result_diplotype"],
                                                                        DATA["gene_result_lookup"],
                                                                        DATA["gene_result"])
    
    activity_two = [i for i in activityscore["CYP2C9"] if i["activityscore"] == 2.00][0]
    assert activity_two == {'activityscore': 2.0, 'genotypes': [['*1', '*1'], ['*9', '*9'], ['*1', '*9']]}
    assert phenotype['CACNA1S'][0] == {
        'genotypes': [['c.3257G>A', 'c.3257G>A'],
                      ['c.520C>T', 'c.520C>T'],
                      ['Reference', 'c.520C>T'],
                      ['Reference', 'c.3257G>A'],
                      ['c.3257G>A', 'c.520C>T']],
        'phenotype': 'Malignant Hyperthermia Susceptibility'
    }
    
    assert [i for i in phenotype["HLA-B*15:02"] if i["phenotype"] == "negative"] == {
        # TODO - maybe it should have different format?
        "phenotype": "negative",
        "genotypes": [
            ["negative"]
        ]
    }


def test_normalize_cpic_factors():
    assert normalize_cpic_factors({'HLA-B': '*57:01 positive'}) == {'HLA-B*57:01': 'positive'}
    assert normalize_cpic_factors({"SOME_GENE": "2"}) == {"SOME_GENE": "== 2.00"}


def test_create_cpic_recommendations():
    assert CPIC_RECOMMENDATIONS["abacavir"] == [{
        'factors': {
            'HLA-B*57:01': 'negative', 'population': 'general'
        },
        'guideline': 'https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/',
        'recommendation': 'Use abacavir per standard dosing guidelines',
        'strength': 'strong'
    },
        {
            'factors': {
                'HLA-B*57:01': 'positive', 'population': 'general'
            },
            'guideline': 'https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/',
            'recommendation': 'Abacavir is not recommended',
            'strength': 'strong'
        }]


def test_check_guideline():
    def help():
        rec = CPIC_RECOMMENDATIONS
        result = []
        for drugname, values in rec.items():
            guidelines = set()
            for factor in values:
                guidelines.add(factor["guideline"])
            is_drug_in_guideline_url = False
            for i in guidelines:
                if drugname in i:
                    is_drug_in_guideline_url = True
            result.append({
                "drug": drugname,
                "guidelines": guidelines,
                "count_guidelines": len(guidelines),
                "is_drug_in_guideline_url": is_drug_in_guideline_url
            })
        return result
    
    for i in help():
        if i["is_drug_in_guideline_url"] == False:
            print(i["drug"])
            print(i["guidelines"])


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
    with open(os.path.join(cwd, "fixtures/cpic.sql"), 'r') as sql_file:
        for table, row in yield_rows_from_sql_file(sql_file):
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
