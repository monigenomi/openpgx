import glob
import ntpath
from collections import defaultdict
from os import path

import bs4

from .helpers import (
    load_json,
    is_star,
    get_phenoconversion_data_from_recommendations,
    format_with_populations,
    download_to_cache_dir,
)


def extract_gene_name(dpwg_filename: str) -> str:
    return dpwg_filename.split("_and_")[1].split(".json")[0].replace("_", "-")


def extract_drug_name(dpwg_filename: str) -> str:
    return dpwg_filename.split("_and_")[0].split("_")[-1]


def detect_engine():
    try:
        import lxml
    except ImportError:
        engine = "html.parser"
    else:
        engine = "lxml"
    return engine


class Converter:
    def __init__(self, **kwargs):
        engine = kwargs.get("engine")
        if engine is None:
            self.engine = detect_engine()
        else:
            self.engine = engine
        self.params = kwargs

    def convert(self, html_doc):
        soup = bs4.BeautifulSoup(html_doc, self.engine)
        output = []
        # TODO do not miss new line \n
        for table_num, table in enumerate(soup.find_all("table"), start=1):
            for tr in table.find_all("tr"):
                row = ["".join(cell.strings) for cell in tr.find_all(["td", "th"])]
                output.append([table_num, *row])

        return output


def table_from_html(html_text: str) -> list:
    c = Converter()
    return c.convert(html_text)


def tables_to_dicts(tables: list):
    if len(tables) < 1:
        return tables

    temporary = {}
    result = {}
    temporary[1] = []
    result[1] = []
    for line in tables:
        if line[0] == 1:
            temporary[1].append(line)
        continue
    header = temporary[1][0]
    for line in temporary[1][1:]:
        line_data = {}
        for cell in range(1, len(line)):
            line_data[header[cell]] = line[cell]
        result[1].append(line_data)

    return result


def html_to_table_of_recommendations(html_text: str) -> list:
    dicts = tables_to_dicts(table_from_html(html_text))

    if not dicts:
        raise Exception(f"No table in {html_text}")

    return list(dicts.values())[0]


FACTOR_NORMALIZATION = {
    "PM": "poor metabolizer",
    "UM": "ultrarapid metabolizer",
    "NM": "normal metabolizer",
    "IM": "intermidiate metabolizer",
    "VKORC1 -1639 AA": "rs9923231 variant (T)",
    "VKORC1 rs9923231 AA": "rs9923231 variant (T)",
    "VKORC1 rs9923231 AG": "rs9923231 reference (C)",
    "VKORC1 -1639 AG": "rs9923231 reference (C)",
    # TODO: check if is possible to use in this recommendations poor function and indeterminate function for SLCO1B1
    #  (phenoconversion from CPIC, the same as is named in fda: "521 TC or 521 CC (intermediate or poor function transporters)"
    "SLCO1B1 521 CC": "521 CC",
    "SLCO1B1 521 TC": "521 TC",
    # TYPO in database
    "SLCO1B1 512 CC": "521 CC",
    "SLCO1B1 512 TC": "521 TC",
    # DPYD
    "Activity Score 0": "== 0.00",
    "Activity Score 1": "== 1.00",
    "Activity Score 1.5": "== 1.50",
    # problem with DPYD and only Activity Score
    "DPD AS 0": "== 0.00",
    "DPD AS 1": "== 1.00",
    "DPD AS 1.5": "== 1.50",
    "DPD FENO": None,
    "FENO": None,
    "HLA-B*44": "*44 positive",
    ##
    "Factor V Leiden heterozygous": "Factor V Leiden heterozygous",
    "Factor V Leiden homozygous": "Factor V Leiden homozygous",
    "CYP3A5 heterozygote expressor": "CYP3A5 heterozygote expressor",
    "CYP3A5 homozygous expressor": "CYP3A5 homozygous expressor",
}


def normalize_dpwg_factor(factor: str) -> str:
    for key, value in FACTOR_NORMALIZATION.items():
        if key in factor:
            return value

    if "HLA-" in factor:
        if factor[-2::1].isdigit():
            factor = factor.replace(":", "")
            return factor[-5:-2:1] + ":" + factor[-2::1] + " positive"

    raise Exception("Unknown factor: " + factor)


def get_recommendations_by_factors(html_text: str) -> dict:
    result = {}

    # TODO: make table for rasburicase: https://www.pharmgkb.org/guidelineAnnotation/PA166119846
    for row in html_to_table_of_recommendations(html_text):
        factor = row.get("ALLELE/GENOTYPE/PHENOTYPE") or row.get(
            "Allele/Genotype/Phenotype"
        )

        if factor is None or is_star(factor):
            continue

        factor = normalize_dpwg_factor(factor)

        recommendation = row.get("RECOMMENDATION") or row.get("Recommendation")

        if "Recommendation LIVER transplantation" in row:
            recommendation = " ".join(
                [
                    "Recommendation Indications OTHER than liver transplantation: ",
                    row["Recommendation Indications OTHER than liver transplantation"],
                    "Recommendation LIVER transplantation: ",
                    row["Recommendation LIVER transplantation"],
                ]
            )

        result[factor] = recommendation

    return result


def load_dpwg_entry(gene_drug_filename: str) -> dict:
    filename = ntpath.basename(gene_drug_filename)
    gene_drug_dict = load_json(gene_drug_filename)

    # There are invalid links in dwpg json data. All of them "No page was found. Exception: link to guideline where everything is described.
    # Additionally in "citations" and "literature" commented below there are no all resources:
    # Even though recommendations are based on 2018, 2019 or 2020 DPWG update in citations is old publication from 2011 (... bench... etc.)
    # Unfortunately in secion guideline.literature there is no year and journal, but link to PHAMRGKB website is included, where all data is available.

    publications = []
    # if "citations" in gene_drug_dict:
    #     for citation in gene_drug_dict["citations"]:
    #         publications.append(citation["title"] + " (" + str(citation["year"]) + ")")

    guideline = None
    if "guideline" in gene_drug_dict and "@id" in gene_drug_dict["guideline"]:
        guideline = gene_drug_dict["guideline"]["@id"]
        if "literature" in gene_drug_dict["guideline"]:
            for item in gene_drug_dict["guideline"]["literature"]:
                title = (
                    item["title"] if item["title"][-1] != "." else item["title"][:-2:1]
                )
                publications.append(title)

    drugname = extract_drug_name(filename)

    recommendations_by_factor = {}

    if gene_drug_dict["guideline"]["recommendation"]:
        html_text = gene_drug_dict["guideline"]["textMarkdown"]["html"]
        recommendations_by_factor = get_recommendations_by_factors(html_text)

    return {
        "gene": extract_gene_name(filename),
        "drug": drugname,
        "recommendations_by_factor": recommendations_by_factor,
        "guideline": guideline,
        "summary": gene_drug_dict["guideline"]["summaryMarkdown"]["html"],
        "population": "general"
        if gene_drug_dict["guideline"]["pediatric"]
        else "adults",
    }


# TODO - saved json contains only default recommendations with empty factors, but shoudln't

DPWG_DEFAULT_URL = (
    "https://api.pharmgkb.org/v1/download/file/data/dosingGuidelines.json.zip"
)


def get_dpwg_recommendations(url=DPWG_DEFAULT_URL) -> dict:
    cache_dir = download_to_cache_dir(url, "dpwg")

    wildcard = path.join(
        path.dirname(path.realpath(__file__)),
        f"{cache_dir}/Annotation_of_DPWG_*.json",
    )

    entries = [load_dpwg_entry(filename) for filename in glob.glob(wildcard)]

    recommendations_by_drug_and_population = defaultdict(lambda: defaultdict(list))

    for entry in entries:
        if len(entry["recommendations_by_factor"]):
            default_recommendation = {
                "factors": {},
                # Taken from https://api.pharmgkb.org/v1/download/file/attachment/DPWG_November_2018.pdf
                "recommendation": "NO action is needed for this gene-drug interaction.",
                "guideline": entry["guideline"],
            }
            recommendations_by_drug_and_population[entry["drug"]]["general"].append(
                default_recommendation
            )
            continue

        for factor, raw_recommendation in entry["recommendations_by_factor"].items():
            gene = entry["gene"]
            if "HLA-" in gene:
                if " positive" in factor:
                    gene = gene + factor.replace(" positive", "")
                    factor = "positive"
                elif " negative" in factor:
                    gene = gene + factor.replace(" negative", "")
                    factor = "negative"

            recommendation = {
                "factors": {gene: factor},
                "recommendation": raw_recommendation,
                "guideline": entry["guideline"],
            }

            assert (
                entry["population"]
                not in recommendations_by_drug_and_population[entry["drug"]]
            )

            recommendations_by_drug_and_population[entry["drug"]][
                entry["population"]
            ] = recommendation

    result = defaultdict(list)
    for (
        drug,
        recommendations_by_population,
    ) in recommendations_by_drug_and_population.items():
        full_recommendation = format_with_populations(recommendations_by_population)
        result[drug].append(full_recommendation)

    return dict(result)


def get_dpwg_phenoconversion_data(recommendations):
    return get_phenoconversion_data_from_recommendations(recommendations)
