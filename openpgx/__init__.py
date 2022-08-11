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
    "cpic": create_cpic_database(),
    "dpwg": create_dpwg_database(),
    "fda": create_fda_database(),
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

    encoding is the value of factor, e.g. "== 5.25", "poor metabilizer", "positive"
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


def get_recommendation_for_drug(database: dict, drug: str, encodings: str):
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
        "encodings" : {"CYP2D6": ["poor metabolizer", 0.0, ...}
    
    """
    if drug not in database["recommendations"]:
        return None

    drug_recommendations = database["recommendations"][drug]

    matched_recommendations = []

    for recommendation in drug_recommendations:
        if recommendation_matches_genotype(recommendation, encodings):
            matched_recommendations.append(recommendation)

    if len(matched_recommendations) > 0:
        recommendation = get_best_recommendation(matched_recommendations)

        return recommendation

    return None


def verify_vendor_database(data):
    recommendation_gene_names = data.keys()
    recommendation_factor_names = [d["factors"] for d in data.values()]


def create_database():
    # TODO make default sources
    result = {}

    for data in ["cpic", "dpwg", "fda"]:
        result[data] = DATABASES[data]

    return result


def get_drugs(database) -> list:
    drugs = []
    for source_database in database.values():
        drugs.extend(source_database["recommendations"].keys())
    return drugs

def phenotyping(genotypes: dict, database: dict ) -> dict:
    """
    Performs translation, changing genotype to encoding according to encodings taken from databases.
    genotype: according to main input example.json
    database: dictionary with databases names as keys (cpic, fda, dpwg) and "recommendations" and "encodings"
    """
    cpic_encodings = database["cpic"]["encodings"] #TODO implement encodings from DPWG and FDA also
    phenotyping_result = {}
    for gene, genotype in genotypes.items():
        sorted_genotype = "/".join(sorted(genotype.split("/")))
        phenotyping_result[gene] = []
        if gene in cpic_encodings and sorted_genotype in cpic_encodings[gene]:
            phenotyping_result[gene] = cpic_encodings[gene][sorted_genotype]
    return phenotyping_result
    
    
def get_recommendations_for_patient(genotypes: dict) -> dict:
    """
    1. Creates database with all databases data (cpic, fda, dpwg). Including recommendations + encodings for each
    2. Performs phenotyping (using encodings). Example:
        "encodings": {"NUDT15": {
            "*1/*1": [
              "normal metabolizer"
            ],
            "*1/*3": [
              "intermediate metabolizer"
            ]}
    2. Creates recommendation dictionary for each drug in database
        if genotype translated to encoding matches factors in recomendation

    genotype: dictionary with all genes that were genotyped for specific patient, according to example.json
    """
    recommendations = defaultdict(dict)

    database = get_database()

    drugs = get_drugs(database)
    
    genotypes_translated_to_encodings = phenotyping(genotypes, database)

    
    for source, source_database in database.items():
        for drug in drugs:
            recommendations[drug] == defaultdict(dict)
            recommendation = get_recommendation_for_drug(
                source_database, drug, genotypes_translated_to_encodings
            )

            if recommendation != None:
                recommendations[drug][source] = recommendation

    return dict(recommendations)
