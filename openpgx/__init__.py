import os
from pkgutil import get_data
import re
from collections import defaultdict
from typing import Optional

from loguru import logger
from numpy import source

from openpgx.cpic import create_cpic_database
from openpgx.dpwg import create_dpwg_database
from openpgx.fda import create_fda_database

from openpgx.helpers import words_to_sentence, get_database

DATABASES = {
    "cpic": create_cpic_database,
    "dpwg": create_dpwg_database,
    "fda": create_fda_database,
}


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
    """
    Checks if encoding matches factor

    encoding is the value of factor, e.g. "5.25", "poor metabilizer", "positive"
    factor is:
        - activity score: a range for which encoding matches factor, e.g. ">= 2.0"
        or
        - phenotype
        or
        - genotype (in HLA cases)
    
    """
    # activity score: "== 2.00" and ">= 2.00"
    if '= ' in encoding:
        factor_operator, factor_value = factor[0:2], factor[2:]
        if factor_operator == "==":
            return float(encoding) == float(factor_value)
        elif factor_operator == ">=":
            return float(encoding) >= float(factor_value)

    # In factors other than activity score (phenotype, genotype)
    return encoding == factor



def recommendation_matches_genotype(recommendation: dict, genotype: dict) -> bool:
    """
    Checks if all factors for specific drug z
    Usage:

    - recommendation["factors"]
      dictionary for specific drug, where keys are factor names and values are allowed factor range values:
      for example:
        {
            'HLA-B*57:01': 'negative',
            # 'population': 'general'
            }
    - genotype
        dictionary where key is genename and value is genotype (haplotype or diplotype)
        for example: 
            {'CYP2D6': '*1/*1'}

    
    """
    if len(genotype) == 0:
        return len(recommendation["factors"]) == 0

    for gene, factor in recommendation["factors"].items():
        if gene not in genotype.keys():
            return False # filter out other factors than genes (example: population)

        if factor is None:
            return False

        if not does_encoding_match_factor(genotype[gene], factor):
            return False

    return True


def get_recommendation_for_drug(database: dict, drug: str, genotype: str):
    """
    Gets best matched recommention for specific drug in specific database (cpic, dpwd, fda)
    database:
        for example:
    {
        'recommendations': {
            'abacavir': {
                'factors': {'HLA-B*57:01': 'negative', 'population': 'general'},
                'recommendation': 'Use abacavir per standard dosing guidelines',
                'strength': 'strong',
                'guideline': 'https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/'
                    },
                'allopurinol': [...]
                },
        "encodings" : {...}
    
    """
    if drug not in database["recommendations"]:
        return None

    drug_recommendations = database["recommendations"][drug]

    matched_recommendations = []

    for recommendation in drug_recommendations:
        if recommendation_matches_genotype(recommendation, genotype):
            matched_recommendations.append(recommendation)

    if len(matched_recommendations) > 0:
        recommendation = get_best_recommendation(matched_recommendations)

        return recommendation

    return None


def verify_vendor_database(data):
    recommendation_gene_names = data.keys()
    recommendation_factor_names = [d["factors"] for d in data.values()]


def create_database(sources):
    # TODO make default sources
    result = {}

    for source, source_url in sources.items():
        result[source] = DATABASES[source](source_url)

    return result


def get_drugs(database) -> list:
    drugs = []
    for source_database in database.values():
        drugs.extend(source_database["recommendations"].keys())
    return drugs


def get_recommendations_for_person(genotype: dict) -> dict:
    """
    1. Creates database with all databases data (cpic, fda, dpwg)
    2. Creates recommendation dictionary for each drug in database. Recommedation matches genotype

    """
    recommendations = defaultdict(dict)

    database = get_database()

    drugs = get_drugs(database)

    for source, source_database in database.items():
        for drug in drugs:
            recommendation = get_recommendation_for_drug(
                source_database, drug, genotype
            )

            if recommendation != None:
                recommendations[drug].append(recommendation)

    return dict(recommendations)
