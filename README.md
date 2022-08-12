# OpenPGx [![build](https://github.com/monigenomi/openpgx/workflows/CI/badge.svg)](https://github.com/monigenomi/openpgx/actions) 

OpenPGx is software with useful pharmacogenomics utilities.

It now implements to convert human's genotype to phenotype and annotate with gene-drug recommendations (for all drugs having recommendations in CPIC, DPWG, or FDA pharmacogenomics databases)

## Terminal Usage

```

$ openpgx <input> [-o <output>]
  
  <input> is a path to JSON file with genotypes to filter recommendations
  <output> is a path to where results will be put, in JSON format
  
  Here is an example <input> file that describes person's genotype:
  
    {
      "SLCO1B1", "*1A/*1B",
      "CYP2D6": "*2≥3/*1≥3",
      "HLA-A*31:01": "positive",
      "HLA-B*15:02": "negative",
      "DPYD": "c.601A>C/c.2194G>A (*6)",
      "G6PD": "B (wildtype)"
    }

   
Thank you for using OpenPGx! We really appreciate contrubitions and discussions:
https://github.com/monigenomi/openpgx

```

## Python API Usage

```python
import openpgx

openpgx.recommendations({ "SLCO1B1", "*1A/*1B" })
```

## Development

Create virtual environment and install dependencies with:

```
bash setup.sh
```

Then it's recommended that you use VSCode for development. Tests in CLI can be run with `pytest -vv`.

If you wish to compute recommendations from raw databases, you can use `openpgx update` command:

```sh
$ openpgx update
  
  Re-computes recommendations database for CPIC, DPWG, and FDA.
  
  You can use this command for example to test and implement new kinds of recommendation.
  
  Options:
    --cpic   Link or path from which to fetch CPIC recommendations, default:
       https://github.com/cpicpgx/cpic-data/releases/download/v1.10/cpic_db_dump-v1.10_inserts.sql.gz
    --dpwg   Link or path from which to fetch DPWG recommendations, default:
       https://api.pharmgkb.org/v1/download/file/data/dosingGuidelines.json.zip
    --fda    Link or path from which to fetch FDA recommendations, default:
       https://raw.githubusercontent.com/PharmGKB/fda-biomarker/master/fda_pgx_associations_table.json
```

Some tips:

- Please add tests for each change you make
- Recommendations database generation can be created with `openpgx update` command
- It's best to use `python -m openpgx` during development to handle module loading properly

## License

EUPL v1.2
