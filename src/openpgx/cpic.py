import re
from collections import defaultdict
from typing import Optional

import psycopg2.extras
from loguru import logger

from .helpers import normalize_hla_gene_and_factor

conn = psycopg2.connect(
    "postgresql://monika@localhost/cpicv1.8?options=-c%20search_path%3Dcpic",
    cursor_factory=psycopg2.extras.DictCursor,
)

data_queries = {
    "allele": "SELECT id, genesymbol, name, functionalstatus, definitionid, clinicalfunctionalstatus, activityvalue, citations FROM "
    "allele",
    "allele_definition": "SELECT * FROM allele_definition",
    "gene": "SELECT * FROM gene",
    "pair": "SELECT genesymbol, pairid, drugid, guidelineid, cpiclevel, pgkbcalevel, pgxtesting, citations FROM pair",
    "guideline": "SELECT * FROM guideline",
    "drug": "SELECT drugid, name, pharmgkbid, rxnormid, drugbankid, atcid, guidelineid FROM drug",
    "publication": "SELECT * FROM publication",
    "recommendation": "SELECT * FROM recommendation",
    "gene_result_diplotype": "SELECT * FROM gene_result_diplotype",
    "gene_result_lookup": "SElECT * FROM gene_result_lookup",
    "test_alert": "SELECT * FROM test_alert",
    "gene_result": "SELECT * FROM gene_result",
}

data_keys = {"drug": "drugid", "pair": "pairid", "gene": "symbol"}

DATA = {}
INDEXES = {}


def fetch(query: str) -> list:
    with conn.cursor() as curs:
        curs.execute(query)
        return [dict(result) for result in curs.fetchall()]


def get_index(columns) -> any:
    if callable(columns):
        return columns

    if type(columns) == tuple:
        return lambda item: tuple([item[column] for column in columns])

    return lambda item: item[columns]


def create_index(data: list, columns: any) -> dict:
    index_fn = get_index(columns)
    index = defaultdict(list)
    for item in data:
        key_or_keys = index_fn(item)
        if type(key_or_keys) == list:
            for key in key_or_keys:
                index[key].append(item)
        else:
            index[key_or_keys].append(item)

    return dict(index)


def normalize(table: str, record: any) -> any:
    if table == "recommendation":
        if record["drugrecommendation"] == "No recommendation":
            return None

        for key in list(record["lookupkey"].keys()):
            match = re.match(r"^No (.+) result$", record["lookupkey"][key])
            if match:
                record["lookupkey"][key] = match[1] + " n/a"

    return record


def fetch_all(table):
    if table not in DATA:
        data = []
        for row in fetch(data_queries[table]):
            normalized_row = normalize(table, row)
            if normalized_row:
                data.append(normalized_row)
        DATA[table] = data

    return DATA[table]


def select(table: str, columns: any, values: any) -> list:
    if type(values) != list:
        values = [values]

    rows = fetch_all(table)

    index_name = (table, columns)

    if hasattr(columns, "__name__"):
        index_name = (table, columns.__name__)

    if index_name not in INDEXES:
        INDEXES[index_name] = create_index(rows, columns)

    index = INDEXES[index_name]

    results = {}
    for value in values:
        if value in index:
            for record in index[value]:
                key_fn = get_index(data_keys.get(table, "id"))
                results[key_fn(record)] = record

    return list(results.values())


def find(table: str, columns: any, values: any) -> Optional[dict]:
    results = select(table, columns, values)

    if len(results) > 0:
        return results[0]

    return results


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


def get_factors_for_recommendation(recommendation) -> list:
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


def get_cpic_recommendations() -> dict:
    result = defaultdict(list)

    for recommendation in fetch_all("recommendation"):
        # TODO: check all populations values
        # if recommendation["population"] not in ["general", "adults"]:
        #     continue

        drug = select("drug", "drugid", recommendation["drugid"])[0]
        guideline = select("guideline", "id", recommendation["guidelineid"])[0]

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

    return dict(result)


def get_genotype_index(genesymbol: str, diplotype: str) -> str:
    sorted_diplotype = "/".join(sorted(diplotype.split("/")))
    if "HLA-" in genesymbol:
        for i in [" positive", " negative"]:
            if i in diplotype:
                return genesymbol + ":" + i[1:]
    return genesymbol + ":" + sorted_diplotype


def get_cpic_phenoconversion_data(recommendations):
    result = {}
    for gene_result_diplotype in fetch_all("gene_result_diplotype"):
        gene_result_lookup = find(
            "gene_result_lookup", "id", gene_result_diplotype["functionphenotypeid"]
        )
        gene_result = find("gene_result", "id", gene_result_lookup["phenotypeid"])

        activity_score = normalize_activityscore(gene_result["activityscore"])

        gene_name, phenotype = normalize_hla_gene_and_factor(
            gene_result["genesymbol"], gene_result["result"]
        )

        index = get_genotype_index(gene_name, gene_result_diplotype["diplotype"])

        result[index] = [
            phenotype,
            float(activity_score[3:]) if activity_score else None,
        ]

    return result


# weź funkcje≤ która będzie przyjmowała input od użytkownika i zwracała fenotypy
