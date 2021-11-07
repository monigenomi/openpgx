import re
from collections import defaultdict
from typing import Optional

from .cpic import get_cpic_phenoconversion_data, get_cpic_recommendations, is_haplo_or_diplo
from .dpwg import get_dpwg_phenoconversion_data, get_dpwg_recommendations
from .fda import get_fda_phenoconversion_data, get_fda_recommendations
from .helpers import PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC, words_to_sentence


def index_recommendations(all_recommendations: list) -> dict:
    result = defaultdict(lambda: {"cpic": [], "dpwg": [], "fda": []})

    for recommendation in all_recommendations:
        result[recommendation["drug"]][recommendation["source"]].append(recommendation)

    return result


RECOMMENDATIONS = {}


def get_all_recommendations():
    if len(RECOMMENDATIONS) == 0:
        RECOMMENDATIONS["cpic"] = get_cpic_recommendations()
        RECOMMENDATIONS["fda"] = get_fda_recommendations()
        RECOMMENDATIONS["dpwg"] = get_dpwg_recommendations()
    return RECOMMENDATIONS


PHENOCONVERSIONS = {}


def get_all_phenoconversions():
    recommendations = get_all_recommendations()
    if len(PHENOCONVERSIONS) == 0:
        PHENOCONVERSIONS["cpic"] = get_cpic_phenoconversion_data(
            recommendations["cpic"]
        )
        PHENOCONVERSIONS["fda"] = get_dpwg_phenoconversion_data(recommendations["dpwg"])
        PHENOCONVERSIONS["dpwg"] = get_fda_phenoconversion_data(recommendations["fda"])
    return PHENOCONVERSIONS


# def get_gene_info(name: str, genotype: list) -> Optional[dict]:
#     genotype = normalize_genotype(name, genotype)
#
#     result = {"name": name, "genotype": genotype}
#
#     cpic = get_cpic_info(name, genotype)
#     if cpic:
#         result["cpic"] = cpic
#
#     if "cpic" not in result and "dpwg" not in result and "fda" not in result:
#         logger.error("Gene does not exist in any database", gene=name)
#         return None
#
#     for allele in genotype:
#         if (
#                 "cpic" not in result
#                 and ("fda" not in result or normalize_hla_allele(allele) not in fda["allele"])
#                 and ("dpwg" not in result or normalize_hla_allele(allele) not in dpwg["allele"])
#         ):
#             logger.error(
#                 "Allele does not exist in any database",
#                 genotype=genotype,
#                 allele=allele,
#                 gene=name,
#             )
#             return result
#
#     return result
#

# def get_genes(inputs: dict) -> list:
#     result = []
#
#     for gene, diplotype in inputs.items():
#         genotype = diplotype.split("/")
#         gene_info = get_gene_info(gene, genotype)
#         if gene_info is None:
#             continue
#
#         result.append(gene_info)
#     return result



def single_database_records(gene_allele_table):
    result = []
    for gene in gene_allele_table:
        min_value = 2
        max_value = 2
        input_type = "select"
        if gene["gene"] in is_haplo_or_diplo():
            min_value = 1
        if "HLA-" in gene["gene"]:
            min_value = 1
            max_value = 1
            input_type = "radio"

        alleles = [
            {"name": record["allele"]}
            for record in gene_allele_table
            if record["gene"] == gene["gene"]
        ]

        single = {
            "name": gene["gene"],
            "values": alleles,
            "min_values": min_value,
            "max_values": max_value,
            "input": input_type,
        }
        if single not in result:
            result.append(single)
    return result


def add_database(table: list, name: str) -> list:
    raw_data = single_database_records(table)
    new = []
    for record in raw_data:
        record["database"] = name
        new.append(record)
    return new


def prepare_range(allele: str):
    result = [allele]
    # Converting ranges to sets
    match = re.match(r"^(\*\d+[A-Z]?)x(\d{1,2})$", allele)
    if match:
        allele_name = match.groups()[0]
        how_many = int(match.groups()[1])
        result.extend([f"{allele_name}≥{str(i)}" for i in range(how_many, 0, -1)])

    return result


def get_genotype_indexes(genesymbol: str, genotype: str) -> str:
    result = []

    if "/" in genotype:
        first_allele, second_allele = genotype.split("/")
        for first_index in prepare_range(first_allele):
            for second_index in prepare_range(second_allele):
                result.append(
                    genesymbol + ":" + "/".join(sorted([first_index, second_index]))
                )
    else:
        for index in prepare_range(genotype):
            result.append(genesymbol + ":" + index)

    return result


def phenoconversion(input_factors: dict) -> dict:
    factors = {}

    for gene, genotype in input_factors.items():
        factor, cpic_factor, activityscore = None, None, None

        for index in get_genotype_indexes(gene, genotype):
            if factor and cpic_factor:
                break

            if cpic_factor is None:
                if index in PHENOCONVERSIONS["cpic"]:
                    cpic_factor, activityscore = PHENOCONVERSIONS["cpic"][index]
                    factor = PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC.get(
                        cpic_factor, None
                    )

            if factor is None:
                if index in PHENOCONVERSIONS["dpwg"]:
                    factor = PHENOCONVERSIONS["dpwg"][index]

                elif index in PHENOCONVERSIONS["fda"]:
                    factor = PHENOCONVERSIONS["fda"][index]

        factors[gene] = {
            "factor": factor,
            "cpic_factor": cpic_factor,
            "activityscore": activityscore,
        }

    return factors


def get_all_drugs() -> list:
    drugs = []
    for recommendations_by_drug in RECOMMENDATIONS.values():
        drugs.extend(list(recommendations_by_drug.keys()))
    return set(drugs)


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


def get_recommendations_for_drug(drug: str, factors: dict) -> dict:
    result = {}

    for source, recommendations_by_drug in get_all_recommendations().items():
        recommendations = recommendations_by_drug.get(drug, [])
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

    # TODO: tell that no recommendation was found for given factors
    # TODO: tell user that some factor needs to be provided to select single recommendation:
    # for example: recommendation is for 2 genes and second needs to be provided

    # wypisać dla każdego leku jakie faktory są potrzebne żeby zrobić rekomendację
    # "nie znaleziono rekomendacji ale są możliwe do znalezienia jeżeli przetestuesz ten gen, ten gen i podasz populację

    return result


# @with_logs
def get_recommendations(main_input: list) -> dict:
    drugs = get_all_drugs()
    factors = phenoconversion(main_input)

    return {drug: get_recommendations_for_drug(drug, factors) for drug in drugs}


def test_no_numbers_in_factors():
    for source, recommendations_by_drug in get_all_recommendations().items():
        for drug, recommendations in recommendations_by_drug.items():
            for recommendation in recommendations:
                for gene, factor in recommendation["factors"].items():
                    if factor:
                        match = bool(re.match(r"^\d+(\.\d+)?$", factor))
                        assert not match, f"{source} {drug} {gene} {factor}"
