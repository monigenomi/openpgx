# def test_gene_not_in_database():
#     result = get_all({"some_gene": "x"})
#
#     assert result["errors"] == [{
#         'gene': 'some_gene',
#         'message': 'Gene does not exist in any database'
#     }]
#     assert result["warnings"] == []
#
#
# def test_allele_does_not_exists_in_any_database():
#     result = get_all({"CYP2C9": "*1/*200"})
#     assert result["errors"] == [{
#         'gene': 'CYP2C9',
#         'allele': '*200',
#         'message': 'Allele does not exist in any database'
#     }]
#     assert result["warnings"] == []
#
# def test_allele_exists_in_cpic():
#     result = get_all({"HLA-A": "*31:01 allele positive"})
#     assert result["errors"] == []
#     assert result["warnings"] == []
#
#
# def test_allele_does_not_exists_in_CPIC():
#     result = get_all({"HLA-DRB1": "*07:01 allele negative"})
#     assert result["errors"] == []
#     assert result["warnings"] == [{
#         'gene': 'HLA-DRB1',
#         'allele': "*07:01 allele negative",
#         'message': 'Allele does not exist in CPIC'
#     }]
#
#
# def test_phenotype_does_not_exist():
#     result = get_all({"CYP2C9": "*1/*62"})
#     assert result["errors"] == []
#     assert result["warnings"] == [{
#         'gene': 'CYP2C9',
#         'genotype': ['*1', '*62'],
#         'message': 'Phenotype for genotype does not exist in CPIC database'
#     }]
#

# def test_recommendation_present_no_phenotype():
#     # There is no recommendation for this phenotype
#     result = get_all({"HLA-B": "*15:02 allele negative"})
