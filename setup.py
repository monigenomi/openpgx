# -*- coding: utf-8 -*-
import os

from setuptools import setup


def package(name, authors, **args):
    def get_long_description():
        """Retrieves package description from README.md"""
        try:
            base_path = os.path.dirname(os.path.realpath(__file__))

            with open(os.path.join(base_path, "README.md")) as file:
                description = file.read()
        except FileNotFoundError:
            description = ""
        return description

    setup(
        name=name,
        entry_points={"console_scripts": [f"{name} = {name}.__main__:main"]},
        long_description=get_long_description(),
        long_description_content_type="text/markdown",
        author=", ".join(authors),
        setup_requires=["setuptools>=45", "wheel", "setuptools_scm>=6.2"],
        package_dir={name: f"src/{name}"},
        packages=[name],
        include_package_data=True,
        zip_safe=False,
        **args,
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
    install_requires=["loguru"],
    extras_require={
        "dev": ["black==21.8b0", "pre-commit", "pylint", "pytest", "isort==5.9.3"]
    },
    python_requires="~=3.8",
)
