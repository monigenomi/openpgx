if __name__ == "__main__":
    from os.path import dirname

    __path__ = [dirname(__file__)]
    del dirname

import argparse
import json
from typing import Optional, IO

from .openpgx import get_genes, get_all_recommendations
from .helpers import save_json


def convert_pgx_input(inputs: list) -> dict:
    result = {}
    for item in inputs:
        result[item["gene"]] = item["genotype"]
    return result


class ArgumentParser(argparse.ArgumentParser):
    def print_help(self, file: Optional[IO[str]] = ...) -> None:
        print(
            """
        Gets data about genotype from CPIC database.
        """
        )


# recommendations = get_all_recommendations()
# save_json('recommendations.json', recommendations)

if __name__ == "__main__":
    parser = ArgumentParser(prog="openpgx")

    parser.add_argument(
        "-i",
        "--input",
        help="Path to json with genotypes example: 'gene': '*1/*1'",
        required=True,
    )
    args = vars(parser.parse_args())

    def load_json(path):
        with open(path) as f:
            return json.load(f)

    input = convert_pgx_input(load_json(args["input"]))
    genes = get_genes(input)
    result = get_recommendations(genes)

    print(json.dumps(result, indent=2))
