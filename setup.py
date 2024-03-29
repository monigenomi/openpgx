# -*- coding: utf-8 -*-
import sys
from pathlib import Path

from setuptools import setup

assert sys.version_info >= (3, 8, 0), "openpgx requires Python 3.8.0+"

CURRENT_DIR = Path(__file__).parent
sys.path.insert(0, str(CURRENT_DIR))

def package(name, authors, **args):
    long_description = (CURRENT_DIR / "README.md").read_text(encoding="utf8")

    setup(
        name=name,
        entry_points={"console_scripts": ["{} = {}.__main__:main".format(name, name)]},
        long_description=long_description,
        long_description_content_type="text/markdown",
        author=", ".join(authors),
        setup_requires=["setuptools>=45", "wheel"],
        package_dir={name: name},
        py_modules=[name],
        packages=[name],
        include_package_data=True,
        zip_safe=False,
        **args
    )


package(
    name="openpgx",
    version="0.0.1",
    description="PGx tool for converting star alleles of pharmacogenes to gene-drug recommendations",
    keywords=["pgx", "pharmacogenomics", "cpic", "dpwg", "fda", "bioinformatics"],
    authors=[
        "Monika Krzyżanowska <monigenomi@gmail.com>",
        "Adam Stankiewicz <sheerun@sher.pl>",
    ],
    classifiers=[
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)",
        "Programming Language :: Python :: 3.8",
    ],
    url="https://openpgx.com/",
    project_urls={
        "Bug Tracker": "https://github.com/monigenomi/openpgx/issues",
        "Source Code": "https://github.com/monigenomi/openpgx/",
    },
    install_requires=[
        "loguru",
        "termcolor",
        "bs4",
        "appdirs",
        "numpy"
    ],
    extras_require={
        "dev": [
            "black==22.6.0",
            "pre-commit",
            "pylint",
            "pytest",
            "isort==5.9.3",
            "jupyter",
            "pandas"
        ]
    },
    python_requires="~=3.8",
)
