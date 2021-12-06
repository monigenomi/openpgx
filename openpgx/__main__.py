import argparse
import json
from typing import IO, Optional

from . import get_recommendations
from .cpic import cpic_main


class ArgumentParser(argparse.ArgumentParser):
    def print_help(self, _: Optional[IO[str]] = ...) -> None:
        print(
            """
        Calculates pharmacogenomic recommendations for each supported drug
        """
        )


# recommendations = get_all_recommendations()
# save_json('recommendations.json', recommendations)


def main():
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

    input = load_json(args["input"])
    result = get_recommendations(input)

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    cpic_main()
