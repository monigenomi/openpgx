import json
import os
import re
import tempfile
import traceback
import zipfile
from typing import Any
from collections import defaultdict
from os import path
from urllib.request import Request, urlopen

import appdirs
from loguru import logger
from termcolor import colored


def load_json(json_path: str) -> dict:
    with open(json_path) as f:
        return json.load(f)


def save_json(json_path: str, data: Any):
    with open(json_path, "w") as f:
        return json.dump(data, f, indent=2)


def index_items_by_key(items: list, key: str) -> dict:
    result = defaultdict(list)
    for item in items:
        result[item[key]].append(item)
    return dict(result)


def add_traceback(record):
    if record["level"].name == "ERROR":
        record["msg"] = colored(record["message"], "red")
    elif record["level"].name == "WARNING":
        record["msg"] = colored(record["message"], "yellow")
    else:
        record["msg"] = record["message"]

    tb = traceback.extract_stack()
    tb = [f"{t[1]}: {t[3]}" for t in tb[::-1] if t.filename == tb[-1].filename]
    record["stacktrace"] = "\n".join(list(dict.fromkeys(tb[2:])))


logger.configure(
    handlers=[{"sink": lambda x: x, "format": "{line}: {msg} {extra}\n{stacktrace}\n"}],
    patcher=add_traceback,
)


def with_logs(fn):
    def fn_with_logs(*args, **kwargs):
        warnings = []
        errors = []

        def log(record):
            entry = {"message": record["message"], **record["extra"]}

            if record["level"].name == "WARNING":
                warnings.append(entry)

            if record["level"].name == "ERROR":
                errors.append(entry)

            return False

        handler_id = logger.add(lambda x: x, filter=log)
        result = fn(*args, **kwargs)
        result["warnings"] = warnings
        result["errors"] = errors
        logger.remove(handler_id)
        return result

    return fn_with_logs


def normalize_hla_allele(allele: str) -> str:
    return allele.split(" ")[0]


def is_star(allele):
    return bool(re.search(r"\*\d+[A-Z]?/\*\d+[A-Z]?", allele))


PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC = {
    "Ultrarapid Metabolizer": "ultrarapid metabolizer",
    "Rapid Metabolizer": "ultrarapid metabolizer",
    "Likely Intermediate Metabolizer": "intermediate metabolizer",
    "Possible Intermediate Metabolizer": "intermediate metabolizer",
    "Intermediate Metabolizer": "intermediate metabolizer",
    "Likely Poor Metabolizer": "poor metabolizer",
    "Poor Metabolizer": "poor metabolizer",
    "Normal Metabolizer": "normal metabolizer",
    "uncertain risk of aminoglycoside-induced hearing loss": None,
    "normal risk of aminoglycoside-induced hearing loss": None,
    "ivacaftor responsive in CF patients": None,
    "ivacaftor non-responsive in CF patients": None,
    "increased risk of aminoglycoside-induced hearing loss": None,
    "Uncertain Susceptibility": None,
    "Malignant Hyperthermia Susceptibility": None,
    "Variable": "variable",  # G6PD
    "Deficient": "deficient",  # G6PD
    "Normal": "normal",  # G6PD
    "Decreased Function": "intermediate function",  # SLCO1B1
    "Possible Increased Function": "intermediate function",  # SLCO1B1
    "Possible Decreased Function": "intermediate function",  # SLCO1B1
    "Possible Poor Function": "poor function",  # SLCO1B1
    "Poor Function": "poor function",  # SLCO1B1
    "Normal Function": "normal function",  # SLCO1B1
    "Indeterminate": None,
    # all HLA-A and HLA-B in CPIC are only genes with allele lookup method
    "positive": "positive",
    "negative": "negative",
    None: None,  # "n/a" as value for HLA genes in recommendations
}

PHENOTYPE_AND_ALLELE_NORMALIZATIONS_DPWG_FDA = {}


def get_phenoconversion_data_from_recommendations(recommendations: list) -> dict:
    result = {}
    for recommendations in recommendations.values():
        for recommendation in recommendations:
            for genename, factor in recommendation["factors"].items():
                result[f"{genename}:{factor}"] = factor
    return result


def normalize_hla_gene_and_factor(genename: str, factor: str) -> str:
    if "HLA-" in genename:
        for i in [" positive", " negative"]:
            if i in factor:
                if "*" not in genename:
                    genename = genename + factor.replace(i, "")
                factor = i[1:]

    return genename, factor


def words_to_sentence(words):
    if len(words) == 1:
        return words[0]

    return ", ".join(words[0:-1]) + " and " + words[-1]


POPULATIONS = {
    "PHT use >3mos": "If patient has not previously used phenytoin",
    "PHT naive": "If patient has previously used phenytoin for longer than 3 months without incidence of cutaneous adverse reactions",
    "CBZ use >3mos": "If patient has previously used carbamazepine for longer than 3 months without incidence of cutaneous adverse reactions",
    "CBZ naive": "If patient has not previously used carbamazepine",
    "CBZ-no alternatives": "If there are no alternatives for carbamazepine for patient",
    "OXC naive": "If patient has not previously used oxcarbazepine",
    "OXC use >3 mos": "If patient has previously used oxcarbazepine for longer than 3 months without incidence of cutaneous adverse reactions",
    "child >40kg_adult": "Adults and children > 40 kg",
    "adults": "Adults",
    "general": "General population",
}


def _key_without_description(recommendation: dict) -> str:
    r = recommendation.copy()
    del r["recommendation"]
    del r["population"]
    return ":::".join([f"{a}::{str(b)}" for a, b in r.items()])


def format_with_populations(recommendations_by_population: dict) -> str:
    if len(recommendations_by_population) == 1:
        return list(recommendations_by_population.values())[0]

    result = []

    for key, recommendation in recommendations_by_population.items():
        result.append(POPULATIONS[key] + ": " + recommendation["recommendation"])

    recommendation = list(recommendations_by_population.values())[0]

    return {**recommendation, "recommendation": "\n\n".join(result)}


# assert normalize_gene_and_factor("HLA-A*31:01", "*31:01 positive") == ("HLA-A*31:01", "positive")
def download_url(url: str, save_path: str):
    request = Request(
        url=url,
        headers={
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:91.0) Gecko/20100101 Firefox/91.0"
        },
    )
    response = urlopen(request)
    with open(save_path, "wb") as file:
        file.write(response.read())


def get_cache_dir(subdir=None):
    cache_dir = path.join(appdirs.user_cache_dir("openpgx"), subdir)

    if not path.exists(cache_dir):
        os.makedirs(cache_dir)

    return cache_dir


def download_to_cache_dir(url, subdir_name):
    cache_dir = get_cache_dir(subdir_name)

    if url.endswith(".zip"):
        with tempfile.TemporaryDirectory(prefix="openpgx") as tmpdirname:
            filename_path = path.join(tmpdirname, path.basename(url))
            download_url(url, filename_path)

            with zipfile.ZipFile(filename_path, "r") as zip_ref:
                zip_ref.extractall(cache_dir)

            return cache_dir
    else:
        download_path = path.join(cache_dir, path.basename(url))
        download_url(url, download_path)
        return download_path
