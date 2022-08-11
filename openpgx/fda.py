import re
from collections import defaultdict
from typing import Optional

from openpgx.helpers import (
    is_star,
    load_json,
    download_to_cache_dir,
    normalize_hla_gene_and_factor,
)

STRENGTH_MAPPING = {"recommendation": "strong", "impact": "moderate", "pk": "optional"}

VARIANT_MAPPING_FDA = {
    "-1639G>A": "rs9923231 reference (C)",
    "V433M": "*3 (rs2108622 T, V433M)",
    "521 TC": "521 TC",
    "521 CC": "521 CC",
    "*02:01 positive": "positive",
    "*07:01 positive": "positive",
}


def read_fda_entries(fda_json_path):
    return load_json(fda_json_path)["table"]


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


# TODO: change source to use directly https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations
# Scrap with requests and parse with BeautifulSoup
FDA_DEFAULT_URL = "https://raw.githubusercontent.com/PharmGKB/fda-biomarker/master/fda_pgx_associations_table.json"


def create_fda_database(url: Optional[str] = None) -> dict:
    if url is None:
        url = FDA_DEFAULT_URL

    fda_json_path = download_to_cache_dir(url)

    entries = read_fda_entries(fda_json_path)

    recommendations = defaultdict(list)
    encodings = defaultdict(set)

    for entry in entries:
        gene_name = entry["Gene"]
        for factor in subgroups_to_factors(entry["Affected Subgroups+"]):
            gene, factor = normalize_hla_gene_and_factor(gene_name, factor)
            drug = entry["Drug"].lower()

            recommendations[drug].append(
                {
                    "factors": {
                        gene: factor,
                        # "population": "general",
                        },
                    "recommendation": entry["Description of Gene-Drug Interaction"],
                    "strength": normalize_strength(entry["table"]),
                    "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
                }
            )

            encodings[drug].add(factor)

    encodings = {k: sorted(list(v)) for k, v in encodings.items()}

    return {"recommendations": dict(recommendations), "encodings": dict(encodings)}
