import csv
import gzip
import json
import re
from collections import defaultdict
from typing import Optional, Any

from openpgx.cpic import *
from openpgx.helpers import index_items_by_key, normalize_hla_gene_and_factor
from .helpers import download_to_cache_dir

CPIC_DEFAULT_URL = (
    "https://github.com/cpicpgx/cpic-data/releases/download/v1.15.1/cpic_db_dump-v1.15.1.sql.gz")


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


def normalize_cpic_factor(genename: str, factor: str) -> tuple:
    if re.match(r"≥?\d+(\.\d+)?", factor):
        return genename, normalize_activityscore(factor, True)
    
    if "HLA-" in genename:
        genename, factor = normalize_hla_gene_and_factor(genename, factor)
        return genename, factor
    
    if re.match(r"No (.*?) Result|No Result|n/a", factor):
        return genename, None
    
    return genename, factor


# normalizes record from cpic database dump
def normalize(table: str, record: dict) -> Optional[dict]:
    for key, value in record.items():
        if value is not None:
            if value == '\\N':
                record[key] = None
            elif value == 'f':
                record[key] = False
            elif value == 't':
                record[key] = True
            elif value[0:1] == "{":
                try:
                    record[key] = json.loads(value)
                except json.JSONDecodeError:
                    # We don't need sets like {"asdfa sdfgs","asdfas asdfa"}
                    record[key] = None
    
    if table == "recommendation":
        if record["drugrecommendation"] == "No recommendation":
            return None
        
        record["lookupkey"] = normalize_cpic_factors(record["lookupkey"] or {})
    
    return record


def normalize_cpic_factors(factors: dict) -> dict:
    result = {}
    for gene, factor in factors.items():
        gene, factor = normalize_cpic_factor(gene, factor)
        result[gene] = factor
    return result


def parse_copy(sql: str) -> dict:
    table_name, columns = re.search(r"COPY\s+(?:(?:[^(\.]+)\.)?([^(\.]+)\s+\(([^)]+)\)\s+FROM\s+stdin", sql).groups()
    column_names = [c.strip() for c in columns.split(r",")]
    
    return table_name, column_names


def yield_rows_from_sql_file(sql_file: any):
    table, columns, lines, reading = None, None, [], False
    
    for line in sql_file:
        if re.match("\s*COPY", line):
            table, columns = parse_copy(line)
            reading = True
            continue
        
        if line[0:2] == "\\.":
            reading = False
            reader = csv.DictReader(lines, fieldnames=columns, dialect='excel-tab')
            for record in reader:
                record = normalize(table, record)
                if record is not None:
                    yield table, record
            
            lines = []
            continue
        
        if reading:
            lines.append(line)


def load_cpic_database_from_descriptor(file) -> dict:
    # keys are table names, values are lists of records in them
    data = defaultdict(list)
    
    for table, record in yield_rows_from_sql_file(file):
        data[table].append(record)
    
    return dict(data)


def load_cpic_dump(sql_gz_path) -> dict:
    with gzip.open(sql_gz_path, 'rt') as file:
        return load_cpic_database_from_descriptor(file)


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
        
        phenotype = lookup_row["result"]
        phenotypes = collect_genotypes(genesymbol, phenotype, phenotypes, genotype_row)
    
    return construct_factors_dictionary(activityscores, "activityscore"), construct_factors_dictionary(phenotypes,
                                                                                                       "phenotype")


def create_cpic_database(url: Optional[str] = None) -> dict:
    if url is None:
        url = CPIC_DEFAULT_URL
    
    cached_sql_gz = download_to_cache_dir(url)
    
    data = load_cpic_dump(cached_sql_gz)
    
    guideline_indexed_by_id = index_items_by_key(data["guideline"], "id")
    drug_indexed_by_id = index_items_by_key(data["drug"], "drugid")
    
    dupicates = defaultdict(lambda: defaultdict(lambda: 0))
    recommendations = defaultdict(list)
    for raw in data["recommendation"]:
        drug_name = drug_indexed_by_id[raw["drugid"]][0]["name"]
        
        factors = raw["lookupkey"]
        factors["population"] = raw["population"]
        
        recommendations[drug_name].append({
            "factors": factors,
            "recommendation": raw["drugrecommendation"],
            "strength": raw["classification"].lower(),
            "guideline": guideline_indexed_by_id[raw["guidelineid"]][0]["url"]
        })
        key = str(recommendations[drug_name][-1])
        dupicates[drug_name][key] += 1
    
    recommendations = dict(recommendations)
    
    return {
        "recommendations": recommendations
    }
