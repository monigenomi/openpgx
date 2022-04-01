import gzip
import json
import os
import pickle
import re
from collections import defaultdict
from typing import Optional, Any

from openpgx.cpic import *
from openpgx.helpers import index_items_by_key
from .helpers import (
    download_to_cache_dir,
    normalize_hla_gene_and_factor,
    yield_rows_from_sql_file,
    logger
)

CPIC_DEFAULT_URL = (
    "https://github.com/cpicpgx/cpic-data/releases/download/v1.15.1/cpic_db_dump-v1.15.1_inserts.sql.gz")

DATA_TABLES = [
    'allele',
    'allele_definition',
    'gene',
    'pair',
    'guideline',
    'drug',
    'publication',
    'recommendation',
    'gene_result_diplotype',
    'gene_result_lookup',
    'test_alert',
    'gene_result'
]

DATA = {}
INDEXES = {}


# normalizes row from cpic database, given table name and row
def normalize(table: str, record: dict) -> Optional[dict]:
    if table == "recommendation":
        if record["drugrecommendation"] == "No recommendation":
            return None
        
        if record["lookupkey"]:
            parsed_lookupkey = json.loads(record["lookupkey"])
            for key in list(parsed_lookupkey.keys()):
                match = re.match(r"^No (.+) result$", parsed_lookupkey[key])
                if match:
                    parsed_lookupkey[key] = match[1] + " n/a"
            
            record["lookupkey"] = parsed_lookupkey
        else:
            record["lookupkey"] = {}
    
    return record


def load_cpic_database_from_descriptor(file):
    global DATA, INDEXES
    
    # keys are table names, values are lists of records in them
    DATA = defaultdict(list)
    
    for table_name, row in yield_rows_from_sql_file(file):
        if table_name in DATA_TABLES:
            normalized_row = normalize(table_name, row)
            if normalized_row is not None:
                if table_name == 'recommendation' and normalized_row['id'] == 1767274:
                    print(row)
                    print(normalized_row)
                DATA[table_name].append(normalized_row)
    
    INDEXES = {}
    DATA = dict(DATA)


def load_cpic_database(sql_gz_path):
    with gzip.open(sql_gz_path, 'rt') as file:
        load_cpic_database_from_descriptor(file)


def load_cpic_database_cached(cached_sql_gz: str, force: bool = False):
    global DATA, INDEXES
    cached_file_path = cached_sql_gz + ".pkl"
    
    if force or not os.path.exists(cached_file_path):
        logger.info("Creating cached cpic database")
        load_cpic_database(cached_sql_gz)
        with open(cached_file_path, "wb") as cpicdb:
            pickle.dump(DATA, cpicdb)
            return
    
    logger.info("Loading cached cpic database")
    with open(cached_file_path, "rb") as cpicdb:
        print("loading")
        for key, value in pickle.load(cpicdb).items():
            DATA[key] = value
        
        INDEXES.clear()


def load_cpic_database_from_url(url: str, force: bool = False):
    cached_sql_gz = download_to_cache_dir(url)
    load_cpic_database_cached(cached_sql_gz, force)


def normalize_activityscore(activityscore):
    if activityscore == "n/a" or activityscore == "No result" or activityscore is None:
        return None
    if "≥" not in activityscore:
        return "== " + "{0:.2f}".format(round(float(activityscore) * 4) / 4)
    return ">= " + "{0:.2f}".format(
        round(float(activityscore.replace("≥", "")) * 4) / 4
    )


def normalize_cpic_factor(genename: str, factor: str) -> tuple:
    if re.match(r"≥?\d+(\.\d+)?", factor):
        return genename, normalize_activityscore(factor)
    
    if re.match(r"No (.*?) Result|No Result|n/a", factor):
        return genename, None
    
    if "HLA-" in genename:
        genename, factor = normalize_hla_gene_and_factor(genename, factor)
    
    return genename, factor


def get_alleles(allele_table: list) -> dict:
    alleles = defaultdict(list)
    for raw in allele_table:
        alleles[raw["genesymbol"]].append({"allele": raw["name"]})
    return dict(alleles)


def normalize_genotype(genotype: str) -> str:
    return validate_hla_genotype(genotype)


def validate_hla_genotype(genotype: str) -> str:
    if "positive" in genotype:
        return "positive"
    if "negative" in genotype:
        return "negative"
    return genotype


def normalize_cpic_phenotype(phenotype: str) -> Any:
    phenotype_lower = phenotype.lower()
    # for item in ["normal", "intermediate", "poor"]:
    #     if item in phenotype_lower:
    #         return item + " metabolizer"
    # if "rapid" in phenotype_lower:
    #     return "ultrarapid metabolizer"
    if "decreased" in phenotype_lower:
        return "poor metabolizer"
    if "negative" in phenotype_lower:
        return "negative"
    if "positive" in phenotype_lower:
        return "positive"
    if re.match(r"No (.*?) Result|No Result|n/a", phenotype):
        return None
    return phenotype_lower


def normalize_genename(genename: str, genotype: str) -> str:
    if "HLA-" in genename:
        if "positive" in genotype:
            if "*" not in genename:
                genename = (genename + genotype.replace("positive", "")).strip()
        if "negative" in genotype:
            if "*" not in genename:
                genename = (genename + genotype.replace("negative", "")).strip()
    return genename


def diplotype_to_phenotype_allele_tables(diplotype):
    return sorted(diplotype.split("/"))


def normalize_activityscore(activityscore: str, is_for_factors: bool):
    if activityscore == "n/a" or activityscore == "No result" or activityscore is None:
        return None
    if is_for_factors:
        if "≥" not in activityscore:
            return "== " + "{0:.2f}".format(round(float(activityscore) * 4) / 4)
        return ">= " + "{0:.2f}".format(
            round(float(activityscore.replace("≥", "")) * 4) / 4
        )
    
    return round(float(re.sub(r"[^\d+.]", "", activityscore)) * 4) / 4


def collect_genotypes(genesymbol: str, factor: str, factors: dict, genotype_table: dict) -> dict:
    """factor can be activityscore or phenotype and factors is a dictionary with them """
    diplotype = genotype_table["diplotype"]
    normalized_genesymbol = normalize_genename(genesymbol, diplotype)
    normalized_genotype = sorted(diplotype_to_phenotype_allele_tables(normalize_genotype(diplotype)))
    if factor not in factors[normalized_genesymbol]:
        factors[normalized_genesymbol][factor] = []
    factors[normalized_genesymbol][factor].append(normalized_genotype)
    return factors


def construct_factors_dictionary(factors, key):
    result = {}
    for gene in factors:
        result[gene] = []
        for factor in factors[gene]:
            result[gene].append({
                key: factor,
                "genotypes": factors[gene][factor]
            })
    return result


def create_phenotype_and_activityscore_table(gene_result_diplotype_table: list, gene_result_lookup_table: list,
                                             gene_result_table: list) -> Any:
    indexed_gene_result_lookup = index_items_by_key(gene_result_lookup_table, "id")
    indexed_gene_result = index_items_by_key(gene_result_table, "id")
    
    activityscores = defaultdict(dict)
    phenotypes = defaultdict(dict)
    for genotype_row in gene_result_diplotype_table:
        
        lookup_id_value = indexed_gene_result_lookup[genotype_row["functionphenotypeid"]][0]["phenotypeid"]
        lookup_row = indexed_gene_result[lookup_id_value][0]
        genesymbol = lookup_row["genesymbol"]
        
        activityscore = normalize_activityscore(lookup_row["activityscore"], False)
        if activityscore:
            activityscores = collect_genotypes(genesymbol, activityscore, activityscores, genotype_row)
        
        phenotype = normalize_cpic_phenotype(lookup_row["result"])
        phenotypes = collect_genotypes(genesymbol, phenotype, phenotypes, genotype_row)
    
    return construct_factors_dictionary(activityscores, "activityscore"), construct_factors_dictionary(phenotypes,
                                                                                                       "phenotype")


def normalize_factors_for_recommendation(factors) -> dict:
    """factor_value can be: phenotype, allele (in CPIC only HLA) or activityscore"""
    
    result = {}
    for genesymbol, factor_value in factors.items():
        normalized_genesymbol = normalize_genename(genesymbol, factor_value)
        if re.match(r"\d+(.[0-9])?|No Result|n/a", factor_value):
            result[normalized_genesymbol] = normalize_activityscore(factor_value, True)
        elif "HLA-" in genesymbol:
            result[normalized_genesymbol] = validate_hla_genotype(factor_value)
        else:
            result[normalized_genesymbol] = normalize_cpic_phenotype(factor_value)
    return result


def create_cpic_recommendations() -> dict:
    drug_table = DATA["drug"]
    recommendation_table = DATA["recommendation"]
    guideline_table = DATA["guideline"]
    
    guideline_indexed_by_id = index_items_by_key(guideline_table, "id")
    drug_indexed_by_id = index_items_by_key(drug_table, "drugid")
    
    dupicates = defaultdict(lambda: defaultdict(lambda: 0))
    result = defaultdict(list)
    for raw in recommendation_table:
        drug_name = drug_indexed_by_id[raw["drugid"]][0]["name"]
        
        result[drug_name].append({
            "factors": normalize_factors_for_recommendation(raw["lookupkey"]),
            "recommendation": raw["drugrecommendation"],
            "strength": raw["classification"].lower(),
            "guideline": guideline_indexed_by_id[raw["guidelineid"]][0]["url"],
            "population": raw["population"]
        })
        key = str(result[drug_name][-1])
        dupicates[drug_name][key] += 1
    
    return dict(result)
