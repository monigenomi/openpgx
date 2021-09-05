from .helpers import *

STRENGTH_MAPPING = {"recommendation": "strong", "impact": "moderate", "pk": "optional"}

VARIANT_MAPPING_FDA = {
    "-1639G>A": "rs9923231 reference (C)",
    "V433M": "*3 (rs2108622 T, V433M)",
    "521 TC": "521 TC",
    "521 CC": "521 CC",
    "*02:01 positive": "*02:01 positive",
    "*07:01 positive": "*07:01 positive",
}


def normalize_strength(strength: str) -> str:
    if strength in STRENGTH_MAPPING:
        return STRENGTH_MAPPING[strength]

    raise Exception(f"Unknown strength: {strength}")


def subgroup_to_factor(subgroups: str, subgroup: str) -> str:
    # *57:01 allele positive
    if " positive" in subgroup:
        return subgroup.split(" ")[0] + " positive"

    # *57:01 allele negative
    if " negative" in subgroup:
        return subgroup.split(" ")[0] + " negative"

    # -1639G>A variant carriers
    if "variant carrier" in subgroup:
        return VARIANT_MAPPING_FDA[subgroup.split(" ")[0]]

    # *28/*28 (poor metabolizers)
    for possible in ["ultrarapid", "intermediate", "poor", "normal"]:
        if possible in subgroup.lower():
            if " metabolizer" in subgroups:
                return possible + " metabolizer"
            if " function" in subgroups:
                return possible + " function"

    if is_star(subgroup):
        return None

    return VARIANT_MAPPING_FDA[subgroup]


def subgroups_to_factors(subgroups: str) -> list:
    # This matches sentences that have additional parenthesis
    match = re.match(r"(.+) \((.+)\)", subgroups)
    if match:
        subgroups = " or ".join(match.groups())

    factors = []
    for subgroup in re.split(", or | or |,", subgroups):
        factor = subgroup_to_factor(subgroups, subgroup)
        if factor:
            factors.append(factor)

    return factors


def read_fda_entries(fda_json_path):
    return load_json(fda_json_path)["table"]


def get_fda_phenoconversion_data(recommendations):
    return get_phenoconversion_data_from_recommendations(recommendations)


FDA_DEFAULT_URL = "https://raw.githubusercontent.com/PharmGKB/fda-biomarker/master/fda_pgx_associations_table.json"


def get_fda_recommendations(url=FDA_DEFAULT_URL) -> list:
    fda_json_path = download_to_cache_dir(url, "fda")

    entries = read_fda_entries(fda_json_path)

    result = defaultdict(list)

    for entry in entries:
        for factor in subgroups_to_factors(entry["Affected Subgroups+"]):
            gene, factor = normalize_hla_gene_and_factor(entry["Gene"], factor)
            drug = entry["Drug"].lower()

            result[drug].append(
                {
                    "factors": {gene: factor},
                    "recommendation": entry["Description of Gene-Drug Interaction"],
                    "strength": normalize_strength(entry["table"]),
                    "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
                }
            )

    return dict(result)
