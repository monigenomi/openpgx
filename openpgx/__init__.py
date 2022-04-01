import os
import re
from collections import defaultdict
from typing import Optional

from loguru import logger

from .cpic import CPIC_DEFAULT_URL  # , get_cpic_normalizations, get_cpic_recommendations,
from .dpwg import get_dpwg_recommendations, DPWG_DEFAULT_URL, get_dpwg_normalizations
from .fda import get_fda_recommendations, FDA_DEFAULT_URL, get_fda_normalizations
from .helpers import (
    PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC,
    words_to_sentence,
    get_normalizations, load_json, cache_path, save_json, repository_path
)


def index_recommendations(all_recommendations: list) -> dict:
    result = defaultdict(lambda: {"cpic": [], "dpwg": [], "fda": []})
    
    for recommendation in all_recommendations:
        result[recommendation["drug"]][recommendation["source"]].append(recommendation)
    
    return result


result = {}


def prepare_range(allele: str):
    result = [allele]
    # Converting ranges to sets
    match = re.match(r"^(\*\d+[A-Z]?)x(\d{1,2})$", allele)
    if match:
        allele_name = match.groups()[0]
        how_many = int(match.groups()[1])
        result.extend([f"{allele_name}≥{str(i)}" for i in range(how_many, 0, -1)])
    
    return result


# Generates all possible genotype indexes for gene and its genotype
def get_genotype_indexes(gene: str, genotype: str) -> str:
    result = []
    
    if "/" in genotype:
        first_allele, second_allele = genotype.split("/")
        for first_index in prepare_range(first_allele):
            for second_index in prepare_range(second_allele):
                result.append(
                    gene + ":" + "/".join(sorted([first_index, second_index]))
                )
    else:
        for index in prepare_range(genotype):
            result.append(gene + ":" + index)
    
    return result


def genotype_to_phenotype(genotype: dict) -> dict:
    factors = {}
    
    phenotypings = get_all_phenotypes()
    
    for gene, genotype in input_factors.items():
        factor, cpic_factor, activityscore = None, None, None
        
        for index in get_genotype_indexes(gene, genotype):
            # If we have normalized factor for both CPIC and other databases, go to next gene
            if factor and cpic_factor:
                break
            
            if cpic_factor is None:
                if index in phenotypings["cpic"]:
                    cpic_factor, activityscore = phenotypings["cpic"][index]
                    factor = PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC.get(
                        cpic_factor, None
                    )
            
            if factor is None:
                if index in phenotypings["dpwg"]:
                    factor = phenotypings["dpwg"][index]
                
                elif index in phenotypings["fda"]:
                    factor = phenotypings["fda"][index]
        
        factors[gene] = {
            "factor": factor,
            "cpic_factor": cpic_factor,
            "activityscore": activityscore,
        }
    
    return factors


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


def create_database(*, cpic=None, dpwg=None, fda=None):
    cpic_url = cpic or CPIC_DEFAULT_URL
    dpwg_url = dpwg or DPWG_DEFAULT_URL
    fda_url = fda or FDA_DEFAULT_URL
    
    cpic_recommendations = get_cpic_recommendations(cpic_url)
    cpic_normalizations = get_cpic_normalizations(cpic_recommendations)
    
    dpwg_recommendations = get_dpwg_recommendations(dpwg_url)
    dpwg_normalizations = get_dpwg_normalizations(dpwg_recommendations)
    
    fda_recommendations = get_fda_recommendations(fda_url)
    fda_normalizations = get_fda_normalizations(fda_recommendations)
    
    return {
        "sources": {
            "cpic": cpic_url,
            "dpwg": dpwg_url,
            "fda": fda_url
        },
        "recommendations": {
            "cpic": cpic_recommendations,
            "dpwg": dpwg_recommendations,
            "fda": fda_recommendations
        },
        "normalizations": {
            "cpic": cpic_normalizations,
            "dpwg": dpwg_normalizations,
            "fda": fda_normalizations,
        }
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
