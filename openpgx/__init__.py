import os
from pkgutil import get_data
import re
from collections import defaultdict
from typing import Optional

from loguru import logger

from openpgx.cpic import create_cpic_database
from openpgx.dpwg import create_dpwg_database
from openpgx.fda import create_fda_database

from openpgx.helpers import words_to_sentence, get_database


def index_recommendations(all_recommendations: list) -> dict:
    result = defaultdict(lambda: {"cpic": [], "dpwg": [], "fda": []})

    for recommendation in all_recommendations:
        result[recommendation["drug"]][recommendation["source"]].append(recommendation)

    return result


result = {}


def get_best_recommendation(recommendations: list) -> dict:
    # TODO: probably more factors that number of factors needs to be considered
    def score(recommendation):
        return len(recommendation["factors"].keys())

    return max(recommendations, key=score)


def does_encoding_match_factor(encoding: str, factor: str) -> bool:
    # Case with activity score: "== 2.00" and ">= 2.00"
    factor_operator, factor_value = factor[0:2], factor[2:]
    if factor_operator == "==":
        return float(encoding) == float(factor_value)
    elif factor_operator == ">=":
        return float(encoding) >= float(factor_value)

    return encoding == factor


def recommendation_matches_genotype(recommendation: dict, genotype: dict) -> bool:
    if len(genotype) == 0:
        return len(recommendation["factors"]) == 0

    for gene, factor in recommendation["factors"].items():
        if gene not in genotype:
            return False

        if not does_encoding_match_factor(genotype[gene], factor):
            return False

    return True


def get_recommendations_for_drug(drug: str, genotype: str) -> dict:
    result = {}
    database = get_database()

    drug_recommendations = database[drug]

    for source, recommendations in drug_recommendations.items():
        matched_recommendations = []

        for recommendation in recommendations:
            if recommendation_matches_genotype(recommendation, genotype):
                matched_recommendations.append(recommendation)

        if len(matched_recommendations) > 0:
            result[source] = get_best_recommendation(matched_recommendations)

        elif len(recommendations) > 0:
            factors_of_recommendations = set(
                [f for r in recommendations for f in r["factors"]]
            )
            genes_missing = sorted(
                list(factors_of_recommendations - set(genotype.keys()))
            )
            if len(genes_missing) > 0:
                result[source] = {
                    "factors": {},
                    "recommendation": f"Recommendations are available, but they require genotypes of following genes: {words_to_sentence(genes_missing)}",
                    "guideline": recommendations[0]["guideline"],
                }

    return result


def create_database(*, cpic_url=None, dpwg_url=None, fda_url=None):
    result = {}

    result["cpic"] = create_cpic_database(cpic_url)
    result["dpwg"] = create_dpwg_database(dpwg_url)
    result["fda"] = create_fda_database(fda_url)

    return result


def get_recommendations(genotype: dict) -> dict:
    result = {}

    database = get_database()

    for drug in database:
        cpic_recommendations = get_recommendations_for_drug(drug, genotype)
        dpwg_recommendations = get_recommendations_for_drug(drug, genotype)
        fda_recommendations = get_recommendations_for_drug(drug, genotype)

        result[drug] = cpic_recommendations + dpwg_recommendations + fda_recommendations

    return result
