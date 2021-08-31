#
# def test_fetch_cpic_genes():
#     assert "CYP2D6" in fetch_cpic_genes()
#
#
# def test_fetch_cpic_alleles_for_gene():
#     assert fetch_cpic_alleles_for_gene("VKORC1") == [{'name': 'rs9923231 reference (C)'},
#                                                      {'name': 'rs9923231 variant (T)'}]
#
# def test_single_record_mit():
#     result = single_record("G6PD")
#     assert result['max_values'] == 2
#     assert result['min_values'] == 1


# def test_fda_records_hla():
#     assert single_database_records(
#         [{"allele": "*57:01", "gene": "HLA-B", "variant": "positive"}]
#     ) == [
#         {
#             "name": "HLA-B",
#             "values": [{"name": "*57:01"}],
#             "min_values": 1,
#             "max_values": 1,
#             "input": "radio",
#         },
#     ]
#
#
# def test_G6PD():
#     assert single_database_records(
#         [{"allele": "something", "gene": "G6PD", "variant": None}]
#     ) == [
#         {
#             "name": "G6PD",
#             "values": [{"name": "something"}],
#             "min_values": 1,
#             "max_values": 2,
#             "input": "select",
#         },
#     ]


# def test_gene_with_2_alleles():
#     assert single_database_records(
#         [
#             {"allele": "*1", "gene": "X", "variant": "None"},
#             {"allele": "*2", "gene": "X", "variant": "None"},
#         ]
#     ) == [
#         {
#             "name": "X",
#             "values": [{"name": "*1"}, {"name": "*2"}],
#             "min_values": 2,
#             "max_values": 2,
#             "input": "select",
#         }
#     ]
#
