import os
import re
from collections import defaultdict
from typing import Optional

from loguru import logger

from .cpic import create_cpic_database
from .dpwg import create_dpwg_database
from .fda import create_fda_database

from .helpers import (
    words_to_sentence,
    load_json, cache_path, save_json, repository_path
)


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


def compare_factor(
        factor: Optional[str], activityscore: float, factor_range: Optional[str]
) -> bool:
    if factor_range is None:
        return factor is None
    
    # Case with activity score: "== 2.00" and ">= 2.00"
    operator, value = factor_range[0:2], factor_range[2:]
    
    if operator == "==":
        return activityscore == float(value)
    elif operator == ">=":
        return activityscore >= float(value)
    
    return factor == factor_range


def recommendation_matches_factors(
        source: str, recommendation: dict, factors: dict
) -> bool:
    if len(factors) == 0:
        return False
    
    factor_key = "cpic_factor" if source == "cpic" else "factor"
    
    for gene, factor_range in recommendation["factors"].items():
        if gene not in factors:
            return False
        
        if not compare_factor(
                factors[gene][factor_key], factors[gene]["activityscore"], factor_range
        ):
            return False
    
    return True


def get_recommendations_for_drug(source: str, recommendations: list, factors: dict) -> dict:
    result = {}
    
    matched_recommendations = []
    for recommendation in recommendations:
        if recommendation_matches_factors(source, recommendation, factors):
            matched_recommendations.append(recommendation)
    if len(matched_recommendations) > 0:
        result[source] = get_best_recommendation(matched_recommendations)
    elif len(recommendations) > 0:
        factors_of_recommendations = set(
            [f for r in recommendations for f in r["factors"]]
        )
        genes_missing = sorted(
            list(factors_of_recommendations - set(factors.keys()))
        )
        if len(genes_missing) > 0:
            result[source] = {
                "factors": {},
                "recommendation": f"Recommendations are available, but they require genotypes of following genes: {words_to_sentence(genes_missing)}",
                "guideline": recommendations[0]["guideline"],
            }
    
    return result


def create_database(*, cpic_url=None, dpwg_url=None, fda_url=None):
    return {
        "cpic": create_cpic_database(cpic_url),
        "dpwg": create_dpwg_database(dpwg_url),
        "fda": create_fda_database(fda_url)
    }


DATABASE_PATH = repository_path('database.json')


def load_database() -> dict:
    if not os.path.exists(DATABASE_PATH):
        logger.error('No database present. Please use "openpgx update".')
    
    return load_json(DATABASE_PATH)


def save_database(data: dict) -> dict:
    save_json(DATABASE_PATH, data)


# @with_logs
def get_recommendations(genotype: dict) -> dict:
    database = load_database()
    
    factors = phenotyping(genotype)
    
    recommendations = {}
    
    for drug in database:
        recommendations[drug] = get_recommendations_for_drug(drug, factors)
    
    return recommendations
