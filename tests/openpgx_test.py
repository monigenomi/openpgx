from openpgx import *

database = get_database()

# def test_create_database():
#     create_database(cpic_url=args["cpic"], dpwg_url=args["dpwg"], fda_url=args["fda"])
#     #TODO create this test. This is only for "update" option in main function


def test_get_all_drugs():
    drugs = get_drugs(database)
    assert len(drugs) == 253


def test_recommendation_matches_genotype():
    
    assert recommendation_matches_genotype({
        'factors': {
            'HLA-B*57:01': 'negative', 'population': 'general'
            },
        'recommendation': 'Use abacavir per standard dosing guidelines',
        'strength': 'strong',
        'guideline': 'https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/'
        },{
        "HLA-B*57:01": "negative",
        "population": "general"
        }) == True


def test_get_recommendations_for_drug_HLA_ABACAVIR():
    # Abacavir exists in every base with different recommendation in each
    abacavir = get_recommendation_for_drug(
                database["cpic"],
                "abacavir",
                {
                    "HLA-B*57:01": "positive"
                },
            )
    assert abacavir == {'factors': {'HLA-B*57:01': 'positive'},
     'guideline': 'https://cpicpgx.org/guidelines/guideline-for-abacavir-and-hla-b/',
     'recommendation': 'Abacavir is not recommended',
     'strength': 'strong'}


# def test_get_recommendation_for_drug_empty():
#     #TODO implement message for needed genotyping
#     assert get_recommendation_for_drug(database["cpic"], "allopurinol", {}) == {
#             "factors": {},
#             "guideline": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
#             "recommendation": "Recommendations are available, but they require "
#             "genotypes of following genes: HLA-B*58:01",
#         }
#
#
#     assert get_recommendation_for_drug(database["fda"], "allopurinol", {}) == {
#         "factors": {},
#         "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
#         "recommendation": "Recommendations are available, but they require "
#         "genotypes of following genes: HLA-B*58:01",
#     }


def test_get_recommendation_for_drug_not_empty_cpic():
    assert get_recommendation_for_drug(database["cpic"],"allopurinol", {"HLA-B*58:01": "positive"}) == {
            "factors": {"HLA-B*58:01": "positive"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
            "recommendation": "Allopurinol is contraindicated",
            "strength": "strong",
        }


def test_get_recommendation_for_drug_not_empty_fda():
    assert get_recommendation_for_drug(database["fda"],"allopurinol", {"HLA-B*58:01": "positive"}) == {
            "factors": {"HLA-B*58:01": "positive"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher adverse reaction risk (severe "
            "skin reactions).",
            "strength": "moderate",
        }


def test_get_recommendation_for_drug_not_empty_cpic():
   assert get_recommendation_for_drug(database["cpic"], "allopurinol", {"HLA-B*58:01": "negative"}) == {
            "factors": {"HLA-B*58:01": "negative"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
            "recommendation": "Use allopurinol per standard dosing guidelines",
            "strength": "strong",
        }

# # ma działać na genotype, jednak w tstach dodaję fenotyp, bo na nim jeszcze nie działa! fenotyping
# def test_get_recommendation_for_drug_double():
#     assert get_recommendation_for_drug(database["cpic"], "escitalopram", {"CYP2D6": "== 2.00", "CYP2C19": "poor metabolizer"}) == []


def test_phenotyping():
    # example of input and output to this function
    assert phenotyping({"CYP2D6": "*2≥3/*1≥3"}, database) == {
        "CYP2D6": [ "ultrarapid metabolizer", 6.0 ] }

    def short_(gene, genotype):
        
        result = phenotyping({gene: genotype}, database)[gene]
        return result
    
    # "normal genes"
    assert short_("G6PD", "B (wildtype)") == ["normal", "normal metabolizer"]
    assert short_("CYP2D6", "*7/*7") == ["poor metabolizer", 0.00]
    assert short_("G6PD", "B (wildtype)") == ["normal", 'normal metabolizer']
    assert short_("TPMT", "*4/*10") == [
        "possible intermediate metabolizer", "intermediate metabolizer"
    ]
    assert short_("F5", "Factor V Leiden heterozygous") == ["Factor V Leiden heterozygous"]

    # Activity score obligatory for DPYD
    # TODO: index is not sorted properly: 'DPYD:Reference/c.1898delC (*3)' why result is empty
    # assert short_("DPYD", "c.1898delC (*3)/Reference") == [
    #     "intermediate metabolizer",
    #     1.0,
    # ]
    # assert short_("DPYD", "Reference/c.1905+1G>A") == []

    # Phenotype in CPIC
    #TODO translate *1A/*1B to other genotype - there are no this one in CPIC
    # assert short_("SLCO1B1", "*1A/*1B") == ["normal function", "normal metabolizer"]

    # CPIC phenotyping impossible (allele does not exists) but 521 CC exists in dpwg and fda both
    # TODO implement dwpg and fda
    # assert short_("SLCO1B1", "521 CC") == ["521 CC"]

    # Allele based recommendation - "rs9923231 reference (C)" exists in CPIC but no phenotyping possible
    # assert short_("VKORC1", "rs9923231 reference (C)") == (
    #     #TODO check this case
    #     "rs9923231 reference (C)",
    # )

    # HLA-B*44: Gene and allele exists only in dpwg but recommendation is default whether genotype is given or not
    # assert short_("HLA-B*44", "positive") == ["positive"]

    # "negative and positive is always converted to dpwg and fda for the same name
    # assert short_("HLA-B*57:01", "negative") == ["negative"]
    # assert short_("HLA-B*57:01", "positive") == ["positive"]

    # If allele or gene does not exists return None everywhere
    assert short_("CYP2D6", "*150/*190") == []

    # TODO: What to do with not existing gene?
    assert short_("FOO", "bar") == []

    # assert short_("CYP2D6", "*104/*1x5") ==["indeterminate"] # TODO fix utf in encodings
    assert short_("DPYD", "Reference/c.1905+1G>A (*2A)") == [
        "intermediate metabolizer", 1.0
    ]


def test_get_recommendations_CYPS():
    # Test if "multiple gene" factors works
    recommendations = get_recommendations_for_patient({"CYP2D6": "*7/*7", "CYP2C19": "*1/*2"})

    assert recommendations["trimipramine"]["cpic"] == {
        "factors": {"CYP2D6": "== 0.00", "CYP2C19": "intermediate metabolizer"},
        "recommendation": "Avoid trimipramine use. If a trimipramine is warranted, consider a 50% reduction of recommended starting dose. Utilizing therapeutic drug monitoring to guide dose adjustments is strongly recommended.",
        "strength": "optional",
        "guideline": "https://cpicpgx.org/guidelines/guideline-for-tricyclic-antidepressants-and-cyp2d6-and-cyp2c19/",
    }


def test_get_recommendations_with_multiple_factors():
    recommendations = get_recommendations_for_patient(
        {"HLA-A*31:01": "positive", "HLA-B*15:02": "negative"}
    )
    assert recommendations["carbamazepine"]["cpic"] == {
        "factors": {"HLA-A*31:01": "positive", "HLA-B*15:02": "negative"},
        "guideline": "https://cpicpgx.org/guidelines/guideline-for-carbamazepine-and-hla-b/",
        "recommendation": "If patient is carbamazepine-naïve and alternative agents "
        "are available, do not use carbamazepine.",
        "strength": "strong",
    }


def test_get_recommendations_dpwg_by_activity_score():
    recommendations = get_recommendations_for_patient({"DPYD": "c.601A>C/c.2194G>A (*6)"})
    assert recommendations["capecitabine"]["dpwg"]["factors"] == {"DPYD": "== 1.00"}


def test_compare_activity_score():
    for n in [2.0, 2]:
        assert does_encoding_match_factor("Normal Metabolizer", n, ">= 1.5") == True
        assert does_encoding_match_factor("Normal Metabolizer", n, ">= 2.0") == True
        assert does_encoding_match_factor("Normal Metabolizer", n, ">= 2") == True
        assert does_encoding_match_factor("Normal Metabolizer", n, ">= 2.5") == False
        assert does_encoding_match_factor("Normal Metabolizer", n, "== 2") == True
        assert does_encoding_match_factor("Normal Metabolizer", n, "== 2.0") == True

    does_encoding_match_factor("Normal Metabolizer", 1.0, "Normal Metabolizer") == True
    does_encoding_match_factor("*57:01 negative", 1.0, "*57:01 negative") == True
    does_encoding_match_factor("Ultra Metabolizer", 1.0, "Normal Metabolizer") == False
    does_encoding_match_factor(None, 1.0, None) == True
    does_encoding_match_factor("Ultra Metabolizer", None, None) == False
    does_encoding_match_factor(None, 1.0, "Ultra Metabolizer") == False


def test_prepare_range():
    assert prepare_range("*1x5") == ["*1x5", "*1≥5", "*1≥4", "*1≥3", "*1≥2", "*1≥1"]


def test_no_no_results_in_recommendations():
    for source, recommendations_by_drug in get_all_recommendations().items():
        for drug, recommendations in recommendations_by_drug.items():
            for r in recommendations:
                assert r["recommendation"] != "No recommendation", r


def test_no_duplicate_factors_in_recommendations():
    for source, recommendations_by_drug in get_all_recommendations().items():
        for drug, recommendations in recommendations_by_drug.items():
            existing = {}
            for recommendation in recommendations:
                key = str(recommendation["factors"])
                if key in existing:
                    assert "Duplicate factors in recommendation", (
                        str(recommendation) + "\n" + str(existing[key])
                    )
                existing[key] = recommendation


def test_get_genotype_indexes():
    assert get_genotype_indexes("CYP2D6", "*2≥3/*1≥3") == []
    

def test_check_if_database_contains_proper_vkorc():
    ace = create_database()["dpwg"]["recommendations"]["acenocoumarol"]
    assert len(ace) == 3
    assert type(ace) == list


def test_issue3():
    inputs = {
      "F5": "Factor V Leiden heterozygous",
      "VKORC1": "rs9923231 reference (C)",
    }
    issue3_reco = get_recommendations_for_patient(inputs)
    assert issue3_reco["acenocoumarol"]["dpwg"] == [
         {'factors': {'VKORC1': 'rs9923231 reference (C)'},
          'guideline': 'https://www.pharmgkb.org/guidelineAnnotation/PA166104938',
          'recommendation': 'NO action is needed for this gene-drug interaction'}]

    assert issue3_reco["hormonal contraceptives for systemic use"]["dpwg"] == [{
            'factors': {'F5': 'Factor V Leiden heterozygous'},
            'guideline': 'https://www.pharmgkb.org/guidelineAnnotation/PA166104955',
            'recommendation': '- If the patient has a FAMILY HISTORY WITH A LOT OF '
                              'THROMBOSIS, or has had a PREVIOUS THROMBOSIS:1. Advise the '
                              'prescriber to avoid the use of contraceptives that contain '
                              'oestrogens and prescribe an on-hormone contraceptive-such '
                              'as a copper IUD - as an alternative. One could also opt '
                              'for a progestogen-only contraceptive method, such as the '
                              'depot injection, an IUD with levonorgestrel or an implant '
                              'with etonogestrel.- OTHER CASES:1. Advise the patient to '
                              'avoid additional risk factors for thrombosis (obesity, '
                              'smoking, etc.).'
        }]
    
    
