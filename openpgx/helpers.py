import json
import os
import re
import tempfile
import traceback
import zipfile
from collections import defaultdict
from os import path
from pathlib import Path
from typing import Any, Tuple
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from io import StringIO
from html.parser import HTMLParser

from loguru import logger
from termcolor import colored


def load_json(json_path: str) -> dict:
    with open(json_path) as f:
        return json.load(f)


def save_json(json_path: str, data: Any):
    with open(json_path, "w") as f:
        return json.dump(data, f, indent=2)


def extract_usage(readme_path):
    with open(readme_path, "r") as f:
        contents = f.read()
    
    usage_section = re.search(r"(?<=## Usage)[^#]+", contents).group(0)
    return re.search(r"(?<=```sh\n).+(?=\n```)", usage_section, re.DOTALL).group(0)


def index_items_by_key(items: list, key: str) -> dict:
    result = defaultdict(list)
    for item in items:
        result[item[key]].append(item)
    return dict(result)


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


def format_with_populations(recommendations_by_population: dict) -> dict:
    if len(recommendations_by_population) == 1:
        return list(recommendations_by_population.values())[0]
    
    result = []
    
    for key, recommendation in recommendations_by_population.items():
        result.append(POPULATIONS[key] + ": " + recommendation["recommendation"])
    
    recommendation = list(recommendations_by_population.values())[0]
    
    return {**recommendation, "recommendation": "\n\n".join(result)}


# assert normalize_gene_and_factor("HLA-A*31:01", "*31:01 positive") == ("HLA-A*31:01", "positive")
def download_url(url: str, save_path: str):
    logger.info("Downloading file", url=url, path=save_path)
    request = Request(
        url=url,
        headers={
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:91.0) Gecko/20100101 Firefox/91.0"
        },
    )
    response = urlopen(request)
    with open(save_path, "wb") as file:
        file.write(response.read())


def url_to_cache_dir(url: str) -> str:
    parsed = urlparse(url)
    if url.endswith(".zip"):
        return parsed.hostname + parsed.path[0:-4]
    return parsed.hostname + os.path.dirname(parsed.path)


def repository_path(path: str) -> str:
    return str(Path(__file__).joinpath('../../' + path).resolve())


def cache_path(path: str) -> str:
    return repository_path('.cache/' + path)


def get_cache_dir_for_url(url: str) -> str:
    cache_dir = repository_path('.cache/' + url_to_cache_dir(url))
    
    if not path.exists(cache_dir):
        os.makedirs(cache_dir)
    
    return cache_dir


def download_to_cache_dir(url, force=False):
    if url.endswith(".zip"):
        cache_dir = get_cache_dir_for_url(url)
        
        if not force and len(os.listdir(cache_dir)) > 0:
            return cache_dir
        
        with tempfile.TemporaryDirectory(prefix="openpgx") as tmpdirname:
            filename_path = path.join(tmpdirname, path.basename(url))
            
            download_url(url, filename_path)
            
            with zipfile.ZipFile(filename_path, "r") as zip_ref:
                zip_ref.extractall(cache_dir)
            
            return cache_dir
    else:
        cache_dir = get_cache_dir_for_url(url)
        download_path = path.join(cache_dir, path.basename(url))
        if not force and os.path.exists(download_path):
            return download_path
        download_url(url, download_path)
        return download_path


def add_traceback(record):
    if record["level"].name == "ERROR":
        record["message"] = colored(record["message"], "red")
    elif record["level"].name == "WARNING":
        record["message"] = colored(record["message"], "yellow")
    else:
        record["message"] = record["message"]
    
    tb = traceback.extract_stack()
    tb = [f"{t[1]}: {t[3]}" for t in tb[::-1] if t.filename == tb[-1].filename]
    record["stacktrace"] = "\n".join(list(dict.fromkeys(tb[2:])))


logger.configure(
    handlers=[{"sink": lambda x: x, "format": "{line}: {message} {extra}\n{stacktrace}\n"}],
    patcher=add_traceback,
)


def normalize_hla_gene_and_factor(genename: str, factor: str) -> Tuple[str, str]:
    match = re.match(r"No (.*?) Result|No Result|n/a", factor)
    if match:
        if match.groups()[0] is not None:
            return genename + match.groups()[0], None
        else:
            return genename, None
    
    if " positive" in factor:
        if "*" not in genename:
            genename = genename + factor.replace(" positive", "")
        factor = "positive"
    if " negative" in factor:
        if "*" not in genename:
            genename = genename + factor.replace(" negative", "")
        factor = "negative"
    
    return genename, factor


class MLStripper(HTMLParser):
    def __init__(self):
        super().__init__()
        self.reset()
        self.strict = False
        self.convert_charrefs = True
        self.text = StringIO()
    
    def handle_data(self, d):
        self.text.write(d)
    
    def get_data(self):
        return self.text.getvalue()


def strip_tags(html):
    s = MLStripper()
    s.feed(html)
    return s.get_data()
