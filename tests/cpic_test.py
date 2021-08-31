#!/usr/bin/env python3

from numpy.testing import assert_equal

from src.openpgx.cpic import *
from src.openpgx.helpers import *

CPIC_RECOMMENDATIONS = get_cpic_recommendations()


def test_normalize_activityscore():
    assert_equal(normalize_activityscore("No result"), None)
    assert_equal(normalize_activityscore("n/a"), None)
    assert_equal(normalize_activityscore("1"), "== 1.00")
    assert_equal(normalize_activityscore("≥4"), ">= 4.00")
    assert_equal(normalize_activityscore("4.25"), "== 4.25")


def test_get():
    assert select("gene", "symbol", "CYP2D6")[0]["chr"] == "chr22"


def test_get_factors_for_recommendation():
    result = get_factors_for_recommendation(
        {
            "id": 872151,
            "guidelineid": 100421,
            "drugid": "RxNorm:190521",
            "implications": {
                "HLA-B": "Low or reduced risk of abacavir hypersensitivity"
            },
            "drugrecommendation": "Use abacavir per standard dosing guidelines",
            "classification": "Strong",
            "phenotypes": {},
            "activityscore": {},
            "allelestatus": {"HLA-B": "HLA-B*57:01 negative"},
            "lookupkey": {"HLA-B": "*57:01 negative"},
            "population": "general",
            "comments": "n/a",
            "version": 1,
        }
    )
    assert result == {"HLA-B*57:01": "negative"}


def test_get_cpic_recommendations():
    assert CPIC_RECOMMENDATIONS["rasburicase"] == [
        {
            "factors": {"G6PD": "Normal"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-rasburicase-and-g6pd/",
            "recommendation": "No reason to withhold rasburicase based on G6PD status.",
            "strength": "strong",
        },
        {
            "factors": {"G6PD": "Variable"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-rasburicase-and-g6pd/",
            "recommendation": "To ascertain that G6PD status is normal, enzyme activity "
            "must be measured; alternatives include allopurinol.",
            "strength": "moderate",
        },
        {
            "factors": {"G6PD": "Deficient"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-rasburicase-and-g6pd/",
            "recommendation": "Rasburicase is contraindicated; alternatives include "
            "allopurinol.",
            "strength": "strong",
        },
    ]
    assert CPIC_RECOMMENDATIONS["abacavir"] == [
        {
            "factors": {"HLA-B*57:01": "negative"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/",
            "recommendation": "Use abacavir per standard dosing guidelines",
            "strength": "strong",
        },
        {
            "factors": {"HLA-B*57:01": "positive"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/",
            "recommendation": "Abacavir is not recommended",
            "strength": "strong",
        },
    ]
    assert len(CPIC_RECOMMENDATIONS) == 63


def test_check_lookup_is_valid():
    def create_all_possible_ascores():
        result = []
        for i in range(0, 1000, 5):
            result.append("== " + str("{:.2f}".format(i / 100)))
            result.append(">= " + str("{:.2f}".format(i / 100)))
        return result + [None]

    activity_scores = create_all_possible_ascores()
    cpic_factors = list(PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC.keys())
    for drug, recommendations in CPIC_RECOMMENDATIONS.items():
        for recommendation in recommendations:
            for gene, value in recommendation["factors"].items():
                assert (
                    value in cpic_factors or value in activity_scores
                ), f"{drug} {gene} {value}"


# TODO: what happens if factors: {gene1: None, gene2: Ultrarapid metabolizer} - None is when gene was not provided
def test_phenoconversion():
    result = get_cpic_phenoconversion_data([])

    assert result != {}
    assert result["CYP2D6:*1/*1"] == ["Normal Metabolizer", 2.0]
    assert result["HLA-B*57:01:negative"] == ["negative", None]


def test_normalize_cpic_factor():
    assert normalize_cpic_factor("HLA-A", "No *12 Result") == ("HLA-A", None)
    assert normalize_cpic_factor("AAA", "No Result") == ("AAA", None)
    assert normalize_cpic_factor("FOO", "n/a") == ("FOO", None)
