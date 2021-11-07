# -*- coding: utf-8 -*-
import sys
from pathlib import Path  # noqa E402

from setuptools import setup

assert sys.version_info >= (3, 8, 0), "openpgx requires Python 3.8.0+"

try:
    # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:
    # for pip <= 9.0.3
    from pip.req import parse_requirements

def load_requirements(fname):
    reqs = parse_requirements(fname, session="test")
    return [str(ir.req) for ir in reqs]

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
    install_requires=load_requirements("requirements.txt"),
    extras_require={
        "dev": ["black==21.8b0", "pre-commit", "pylint", "pytest", "isort==5.9.3"]
    },
    python_requires="~=3.8",
)
