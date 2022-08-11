import csv
import gzip
import json
import re
from collections import defaultdict
from typing import Optional, Any

from openpgx.cpic import *
from openpgx.helpers import (
    index_items_by_key,
    normalize_hla_gene_and_factor,
    download_to_cache_dir,
)

CPIC_DEFAULT_URL = "https://github.com/cpicpgx/cpic-data/releases/download/v1.15.1/cpic_db_dump-v1.15.1.sql.gz"


def normalize_activityscore(activityscore: str, is_for_factors: bool):
    if activityscore == "n/a" or activityscore == "No result" or activityscore is None:
        return None
    # Check if activityscore is as. If is genotype or phenotype - return this
    if not any(char.isdigit() for char in activityscore):
        return activityscore
    
    if is_for_factors:
        if "≥" not in activityscore:
            return "== " + "{0:.2f}".format(round(float(activityscore) * 4) / 4)
        return ">= " + "{0:.2f}".format(
            round(float(activityscore.replace("≥", "")) * 4) / 4
        )

    return round(float(re.sub(r"[^\d+.]", "", activityscore)) * 4) / 4


def normalize_cpic_factor(genename: str, factor: str) -> tuple:
    if "HLA-" in genename:
        genename, factor = normalize_hla_gene_and_factor(genename, factor)

    if factor:
        if re.match(r"No (.*?) Result|No Result|n/a", factor):
            factor = None

        if factor:
            factor = normalize_activityscore(factor.lower(), True)

    return genename, factor


def normalize_cpic_factors(factors: dict) -> dict:
    result = {}
    for gene, factor in factors.items():
        gene, factor = normalize_cpic_factor(gene, factor)
        result[gene] = factor
    return result


# normalizes record from cpic database dump
def normalize(table: str, record: dict) -> Optional[dict]:
    for key, value in record.items():
        if value is not None:
            if value == "\\N":
                record[key] = None
            elif value == "f":
                record[key] = False
            elif value == "t":
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


def parse_copy(sql: str) -> dict:
    table_name, columns = re.search(
        r"COPY\s+(?:(?:[^(\.]+)\.)?([^(\.]+)\s+\(([^)]+)\)\s+FROM\s+stdin", sql
    ).groups()
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
            reader = csv.DictReader(lines, fieldnames=columns, dialect="excel-tab")
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
    with gzip.open(sql_gz_path, "rt") as file:
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


def create_cpic_encodings(data) -> dict:
    result = defaultdict(lambda: defaultdict(list))

    indexed_gene_result_lookup = index_items_by_key(data["gene_result_lookup"], "id")
    indexed_gene_result = index_items_by_key(data["gene_result"], "id")

    for diplotype_row in data["gene_result_diplotype"]:
        diplotype = "/".join(sorted(diplotype_row["diplotype"].split("/")))
        phenotypes = indexed_gene_result_lookup[diplotype_row["functionphenotypeid"]]
        gene_result = indexed_gene_result[phenotypes[0]["phenotypeid"]][0]
        genename = normalize_genename(gene_result["genesymbol"], diplotype)

        # First, encoding can be a phenotype name or genotype (mainly in case of HLA)
        normalized_genename, factor = normalize_cpic_factor(
            genename, gene_result["result"]
        )
        result[normalized_genename][diplotype].append(factor)

        # Then optionally gene can be represented by an activity score
        activityscore = normalize_activityscore(gene_result["activityscore"], False)
        if activityscore is not None:
            result[normalized_genename][diplotype].append(activityscore)

    return {k: dict(v) for k, v in result.items()}


def create_cpic_recommendations(data: dict) -> dict:
    """
    Creates dictionary with all cpic data needed to match recommendation for every drug existing in database.
    data:
        authomaticly created dictionary from sql file. Each key corresponds to sql table name
    """
    guideline_indexed_by_id = index_items_by_key(data["guideline"], "id")
    drug_indexed_by_id = index_items_by_key(data["drug"], "drugid")

    recommendations = defaultdict(list)
    for raw in data["recommendation"]:
        drug_name = drug_indexed_by_id[raw["drugid"]][0]["name"]
        drug_recommendations = recommendations[drug_name]

        factors = raw["lookupkey"]
        # factors["population"] = raw["population"]

        drug_recommendations.append(
            {
                "factors": factors,
                "recommendation": raw["drugrecommendation"],
                "strength": raw["classification"].lower(),
                "guideline": guideline_indexed_by_id[raw["guidelineid"]][0]["url"],
            }
        )

    return dict(recommendations)


def create_cpic_database(url: Optional[str] = None) -> dict:
    if url is None:
        url = CPIC_DEFAULT_URL

    cached_sql_gz = download_to_cache_dir(url)
    data = load_cpic_dump(cached_sql_gz)
    recommendations = create_cpic_recommendations(data)
    encodings = create_cpic_encodings(data)

    return {"recommendations": recommendations, "encodings": encodings}
