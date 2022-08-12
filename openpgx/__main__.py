import sys
from argparse import ArgumentParser

from openpgx import create_database, get_recommendations_for_patient
from openpgx.helpers import (
    load_json,
    save_json,
    save_database,
    repository_path,
    logger,
    extract_usage,
    with_logs,
    load_database,
)


def main(args: dict):
    if "positional" not in args or len(args["positional"]) == 0:
        help = extract_usage(repository_path("README.md"))
        print(help)
        return

    command = args["positional"][0]

    if command == "update":
        db = create_database(sources=args)
        save_database(db)

    else:
        genotype = load_json(args["positional"][0])
        recommendations = get_recommendations_for_patient(genotype)
        save_json(args["output"], recommendations)


if __name__ == "__main__":
    logger.add(
        sys.stderr, level="INFO", format="<level>{level: <8}</level> {message} {extra}"
    )

    parser = ArgumentParser(prog="openpgx")
    parser.add_argument("positional", nargs="*")
    parser.add_argument("-o", "--output", default="recommendations.json")
    parser.add_argument("--cpic")
    parser.add_argument("--dpwg")
    parser.add_argument("--fda")
    args = vars(parser.parse_args())

    main(args)
