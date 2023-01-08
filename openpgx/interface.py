from collections import defaultdict
import json
from helpers import load_json, save_json


def main():

    databases = load_json("./database.json")
    genes_alleles_dict = defaultdict(set)
    for database_name, values in databases.items():
        print(database_name)
        encodings = values["encodings"]
        for gene, alleles in encodings.items():
            for haplotype, enc in alleles.items():
                for allele in haplotype.split("/"):
                    genes_alleles_dict[gene].add(allele)

    result = []
    for gene, alleles in genes_alleles_dict.items():
        type = "diplotype"
        if "HLA" in gene:
            type = "select"
        result.append({
            "gene": gene,
            "alleles": sorted(list(alleles)),
            "type": type
        })
    print(json.dumps(result, indent=4))

    save_json("./interfejs.json", result)

main()