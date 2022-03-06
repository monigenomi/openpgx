from openpgx.cpic import *
from openpgx.helpers import index_items_by_key


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


def normalize_cpic_phenotype(phenotype: str) -> str:
    phenotype_lower = phenotype.lower()
    for item in ["normal", "intermediate", "poor"]:
        if item in phenotype_lower:
            return item + " metabolizer"
    if "rapid" in phenotype_lower:
        return "ultrarapid metabolizer"
    if "decreased" in phenotype_lower:
        return "poor metabolizer"
    if "negative" in phenotype_lower:
        return "negative"
    if "positive" in phenotype_lower:
        return "positive"
    if re.match(r"No (.*?) Result|No Result|n/a", phenotype):
        return None
    return phenotype


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
    
    return round(float(re.sub(r"[^\d+\.]", "", activityscore)) * 4) / 4


def collect_genotypes(genesymbol: str, factor: str, factors: dict, genotype_table: dict) -> dict:
    """factor can be activityscore or phenotype and factors is a dictionary with them """
    diplotype = genotype_table["diplotype"]
    normalized_genesymbol = normalize_genename(genesymbol, diplotype)
    normalized_genotype = diplotype_to_phenotype_allele_tables(normalize_genotype(diplotype))
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
    
    return construct_factors_dictionary(activityscores, "activityscore"), \
           construct_factors_dictionary(phenotypes, "phenotype")


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
    
    result = defaultdict(list)
    for raw in recommendation_table:
        drug_name = drug_indexed_by_id[raw["drugid"]][0]["name"]
        
        result[drug_name].append({
            "factors": normalize_factors_for_recommendation(raw["lookupkey"]),
            "recommendation": raw["drugrecommendation"],
            "strength": raw["classification"].lower(),
            "guideline": guideline_indexed_by_id[raw["guidelineid"]][0]["url"]
        })
    return result
