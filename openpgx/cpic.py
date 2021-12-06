import os
from collections import defaultdict
import gzip
import re
import pickle
import json
from typing import Any, Optional

from loguru import logger

from .helpers import (
    download_to_cache_dir,
    normalize_hla_gene_and_factor,
    sql_to_data,
    yield_inserts_from_file,
)

CPIC_DEFAULT_URL = (
    "https://github.com/cpicpgx/cpic-data/releases/download/v1.10/cpic_db_dump-v1.10_inserts.sql.gz"
)

data_tables = [
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

data_keys = {
    "drug": "drugid",
    "pair": "pairid",
    "gene": "symbol"
}


def get_index(columns) -> Any:
    if callable(columns):
        return columns

    if type(columns) == tuple:
        return lambda item: tuple([item[column] for column in columns])

    return lambda item: item[columns]


def create_index(rows: list, columns: Any) -> dict:
    index_fn = get_index(columns)
    index = defaultdict(list)
    for item in rows:
        key_or_keys = index_fn(item)
        if type(key_or_keys) == list:
            for key in key_or_keys:
                index[key].append(item)
        else:
            index[key_or_keys].append(item)

    return dict(index)


class CpicDB:
    def __init__(self, sql_gz_path):
        # keys are table names, values are lists of records in them
        self.data = defaultdict(list)

        with gzip.open(sql_gz_path, 'rt') as file:
            for query in yield_inserts_from_file(file):
                table_name, row = sql_to_data(query)
                if table_name in data_tables:
                    normalized_row = normalize(table_name, row)
                    if normalized_row is not None:
                        self.data[table_name].append(normalized_row)

        self.indexes = {}

    def all(self, table: str):
        return self.data[table]

    def select(self, table: str, columns: Any, values: Any) -> list:
        if type(values) != list:
            values = [values]

        rows = self.data[table]

        index_name = (table, columns)

        if hasattr(columns, "__name__"):
            index_name = (table, columns.__name__)

        if index_name not in self.indexes:
            self.indexes[index_name] = create_index(rows, columns)

        index = self.indexes[index_name]

        results = {}
        for value in values:
            if value in index:
                for record in index[value]:
                    key_fn = get_index(data_keys.get(table, "id"))
                    results[key_fn(record)] = record

        return list(results.values())

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




def find(table: str, columns: Any, values: Any) -> dict:
    results = select(table, columns, values)

    return results[0]



def by_gene_and_diplotype(item):
    return list(item["diplotypekey"].keys())[0], item["diplotype"]


def get_allele_by_haplotype(gene, query):
    haplotypes = [(gene, r) for r in prepare_range(query)]
    result = select("allele", by_allele_range, haplotypes)
    return result[0] if len(result) > 0 else None


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


def get_alleles(gene: str, genotype: list) -> list:
    result = []
    for allele in genotype:
        allele_info = get_allele_info(gene, allele)
        if not allele_info:
            continue
        result.append(allele_info)
    return result


def get_cpic_info(name: str, genotype: list) -> Optional[dict]:
    cpic_genes = select("gene", "symbol", name)

    if len(cpic_genes) == 0:
        return None

    gene = cpic_genes[0]

    alleles = []

    if pass_diplotype(name, genotype):
        alleles = get_alleles(name, genotype)

    if len(alleles) == 0:
        logger.warning(
            "Gene found in CPIC but none of the alleles were",
            gene=name,
            genotype=genotype,
        )
        return None

    activity_score = normalize_activityscore(
        get_activity_score(name, "/".join(genotype))
    )

    phenotype = get_normalized_phenotype(name, genotype)

    if not phenotype:
        logger.error(
            "Phenotype for genotype does not exist in CPIC database",
            genotype=genotype,
            gene=name,
        )

    gene_result = select(
        "gene_result",
        ("genesymbol", "result"),
        (gene, phenotype),
    )[0]

    return {
        "alleles": alleles,
        "activityscore": activity_score,
        "chromosome": gene["chr"],
        "sequence_id": gene["genesequenceid"],
        "chromosome_sequence_id": gene["chromosequenceid"],
        "mRNA_sequence_id": gene["mrnasequenceid"],
        "hgnc_id": gene["hgncid"],
        "ncbi_id": gene["ncbiid"],
        "ensembl_id": gene["ensemblid"],
        "pharmgkb_id": gene["pharmgkbid"],
        "phenotype": phenotype,
        "phenotype_description": gene_result["consultationtext"],
    }


def get_normalized_phenotype(name, genotype):
    phenotype = get_phenotype(name, "/".join(genotype))
    if phenotype is None:
        return None
    return normalize_cpic_phenotype(phenotype)


def normalize_genotype(gene: str, genotype: list) -> list:
    result = []
    for allele_name in genotype:
        normalized_allele_name = get_allele_by_haplotype(gene, allele_name)
        if normalized_allele_name:
            result.append(normalized_allele_name["name"])
        else:
            result.append(allele_name)
    if len(result) == 2:
        diplotype = select("gene_result_diplotype", "diplotype", "/".join(result[::-1]))
        if len(diplotype) > 0:
            return result[::-1]
    return result


def get_allele_info(gene: str, allele_name: str) -> Optional[dict]:
    queried_allele_name = allele_name
    # TODO: Do not treat them specially
    if "positive" in allele_name or "negative" in allele_name:
        queried_allele_name = allele_name.split(" ")[0]

    allele = get_allele_by_haplotype(gene, queried_allele_name)

    if allele:
        allele_definition = select("allele_definition", "id", allele["definitionid"])[0]

        result = {
            "gene": allele["genesymbol"],
            "name": allele["name"],
            "functionalstatus": allele["functionalstatus"],
            "clinicalfunctionalstatus": allele["clinicalfunctionalstatus"],
            "activity_value": allele["activityvalue"],
            "pmid_s": allele["citations"],
            "pharmvar_id": allele_definition["pharmvarid"],
            "structural_variation": allele_definition["structuralvariation"],
        }
    else:
        logger.error(
            "Allele for gene does not exist in CPIC database",
            gene=gene,
            allele=allele_name,
        )
        return None

    if "HLA-" in gene:
        result["variant"] = allele_name.split(" ")[1]

    return result


def is_haplo_or_diplo():
    return [g["symbol"] for g in select("gene", "chr", ["chrX", "chrM"])] + ["VKORC1"]


def get_url_for_publication(publication):
    if publication["url"]:
        return publication["url"]
    if publication["pmid"]:
        return "https://pubmed.ncbi.nlm.nih.gov/" + publication["pmid"] + "/"


def get_factors_for_recommendation(recommendation) -> dict:
    factors = {}
    for gene, factor in recommendation["lookupkey"].items():
        gene, factor = normalize_cpic_factor(gene, factor)
        if factor:
            factors[gene] = factor
    return factors


CLASSIFICATION_TO_STRENGTH = {
    "Optional": "optional",
    "Moderate": "moderate",
    "Strong": "strong",
    "No Recommendation": None,
    "n/a": None,
}


def get_genotype_index(genesymbol: str, diplotype: str) -> str:
    if "HLA-" in genesymbol:
        if " positive" in diplotype:
            return genesymbol + ":positive"
        if " negative" in diplotype:
            return genesymbol + ":negative"

    return genesymbol + ":" + "/".join(sorted(diplotype.split("/")))

def get_cpic_phenoconversion_data():
    pass

def get_records(cached_sql_gz):
    if os.path.exists("./cpicdb.pkl"):
        with open("./cpicdb.pkl", "rb") as cpicdb:
            return pickle.load(cpicdb)
    
    db = CpicDB(cached_sql_gz)
    with open("./cpicdb.pkl", "wb") as cpicdb:
        pickle.dump(db, cpicdb)
    return db

def get_cpic_recommendations(url: str = CPIC_DEFAULT_URL) -> dict:
    result = defaultdict(list)

    cached_sql_gz = download_to_cache_dir(url, "cpic")

    db = get_records(cached_sql_gz)

    for recommendation in db.all("recommendation"):
        # TODO: check all populations values
        # if recommendation["population"] not in ["general", "adults"]:
        #     continue

        drug = db.select("drug", "drugid", recommendation["drugid"])[0]
        guideline = db.select("guideline", "id", recommendation["guidelineid"])[0]

        result[drug["name"]].append(
            {
                "factors": get_factors_for_recommendation(recommendation),
                "recommendation": recommendation["drugrecommendation"],
                "strength": CLASSIFICATION_TO_STRENGTH[
                    recommendation["classification"]
                ],
                "guideline": guideline["url"],
                # "population": recommendation["population"]
            }
        )

    return result


def cpic_main():
    return get_cpic_recommendations()
