import sys
from argparse import ArgumentParser

from openpgx import create_database, save_database, get_recommendations
from openpgx.helpers import (
    load_json,
    save_json,
    repository_path,
    logger,
    extract_usage
)

if __name__ == "__main__":
    logger.add(sys.stderr,
               level="INFO",
               format="<level>{level: <8}</level> {message} {extra}")
    
    parser = ArgumentParser(prog="openpgx")
    parser.add_argument('command', nargs='*')
    parser.add_argument("-g", "--genotype")
    parser.add_argument("-o", "--output", default="recommendations.json")
    parser.add_argument("--cpic")
    parser.add_argument("--dpwg")
    parser.add_argument("--fda")
    args = vars(parser.parse_args())
    command = args["command"]
    
    if len(command) == 0:
        help = extract_usage(repository_path('README.md'))
    
    if command[0] == "update":
        db = create_database(
            cpic_url=args["cpic"],
            dpwg_url=args["dpwg"],
            fda_url=args["fda"]
        )
        save_database(db)
    
    else:
        genotype = {}
        if "input" in args:
            genotype = load_json(args["genotype"])
        recommendations = get_recommendations(genotype)
        save_json(args["output"], recommendations)
    
    # input = load_json(args["input"])
    # recommendations = get_recommendations(input)
    # save_json(args["output"], recommendations)
