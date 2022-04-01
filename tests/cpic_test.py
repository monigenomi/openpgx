#!/usr/bin/env python3

from openpgx.cpic import *
from openpgx.helpers import save_json

load_cpic_database_from_url(CPIC_DEFAULT_URL)


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
    assert normalize_cpic_factor("HLA-A", "No *12 Result") == ("HLA-A", None)
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
        'phenotype': 'malignant hyperthermia susceptibility'
    }
    
    assert [i for i in phenotype["HLA-B*15:02"] if i["phenotype"] == "negative"][0] == {
        # TODO - maybe it should have different format?
        "phenotype": "negative",
        "genotypes": [
            ["negative"]
        ]
    }


def test_normalize_factors():
    assert normalize_factors_for_recommendation({'HLA-B': '*57:01 positive'}) == {'HLA-B*57:01': 'positive'}
    assert normalize_factors_for_recommendation({
        "CYP2D6": "Ultrarapid Metabolizer", "CYP2C19": "Ultrarapid Metabolizer"
    }) == {
               "CYP2D6": "ultrarapid metabolizer", "CYP2C19": "ultrarapid metabolizer"
           }
    assert normalize_factors_for_recommendation({
        "CYP2D6": "0.5", "CYP2C19": "Likely Poor Metabolizer"
    }) == {"CYP2D6": "== 0.50", "CYP2C19": "likely poor metabolizer"}
    assert normalize_factors_for_recommendation({"SOME_GENE": "2"}) == {"SOME_GENE": "== 2.00"}


def test_create_cpic_recommendations():
    assert create_cpic_recommendations()["abacavir"] == [
        {
            "factors": {
                "HLA-B*57:01": "negative"
            },
            "recommendation": "Use abacavir per standard dosing guidelines",
            "strength": "strong",
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/",
            "population": "general"
        },
        {
            "factors": {
                "HLA-B*57:01": "positive"
            },
            "recommendation": "Abacavir is not recommended",
            "strength": "strong",
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/",
            "population": "general"
        }
    ]


def test_write_to_file():
    save_json("/tmp/raw_cpic.json", create_cpic_recommendations())


def test_check_guideline():
    def help():
        rec = create_cpic_recommendations()
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
