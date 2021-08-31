from src.openpgx.dpwg import *
from src.openpgx.helpers import *

DPWG_RECOMMENDATIONS = get_dpwg_recommendations()


def test_extract_gene_name():
    name = "Annotation_of_DPWG_Guideline_for_acenocoumarol_and_VKORC1.json"
    assert extract_gene_name(name) == "VKORC1"


def test_extract_drug_name():
    name = "Annotation_of_DPWG_Guideline_for_acenocoumarol_and_VKORC1.json"
    assert extract_drug_name(name) == "acenocoumarol"


def test_normalize_dpwg_factor():
    assert normalize_dpwg_factor("CYP2D6 UM") == "ultrarapid metabolizer"
    assert normalize_dpwg_factor("HLA-B*57:01") == "*57:01 positive"
    assert normalize_dpwg_factor("HLA-B*44") == "*44 positive"


def test_converter():
    base_dir = path.join(path.dirname(path.realpath(__file__)), "../data/dpwg/")
    raw_data = load_json(
        path.join(
            base_dir, "Annotation_of_DPWG_Guideline_for_aripiprazole_and_CYP2D6.json"
        )
    )
    html_text = raw_data["guideline"]["textMarkdown"]["html"]
    output = Converter().convert(html_text)
    assert type(output) == list
    assert output[0] == [
        1,
        "Allele/Genotype/Phenotype",
        "Drug",
        "Description",
        "Recommendation",
    ]
    assert output != []


def test_therapy_table():
    base_dir = path.join(path.dirname(path.realpath(__file__)), "../data/dpwg/")
    irinotecan = load_json(
        path.join(
            base_dir, "Annotation_of_DPWG_Guideline_for_irinotecan_and_UGT1A1.json"
        )
    )

    assert get_recommendations_by_factors(
        irinotecan["guideline"]["textMarkdown"]["html"]
    ) == {
        "intermidiate metabolizer": "NO action is needed for this gene-drug "
        "interaction.",
        "poor metabolizer": "Start with 70% of the standard dose If the patient "
        "tolerates this initial dose, the dose can be increased, "
        "guided by the neutrophil count.",
    }
    # It should be ampty because "recommendation": false
    gliclazide = load_json(
        path.join(
            base_dir, "Annotation_of_DPWG_Guideline_for_gliclazide_and_CYP2C9.json"
        )
    )
    assert (
        get_recommendations_by_factors(gliclazide["guideline"]["textMarkdown"]["html"])
        == {}
    )


# TODO: Jeżeli nie ma rekomendacji to powinno się dać jakiś tekst z jsona na ten temat.
def test_message_for_no_recommendations():
    pass


def test_load_dpwg_entry():
    base_dir = path.join(path.dirname(path.realpath(__file__)), "../data/dpwg/")
    assert (
        load_dpwg_entry(
            path.join(
                base_dir, "Annotation_of_DPWG_Guideline_for_ribavirin_and_HLA_B.json"
            )
        )["recommendations_by_factor"]
        == {}
    )


def test_get_dpwg_recommendations():
    assert DPWG_RECOMMENDATIONS["irinotecan"] == [
        {
            "factors": {"UGT1A1": "intermidiate metabolizer"},
            "guideline": "https://pharmgkb.org/guidelineAnnotation/PA166104951",
            "recommendation": "NO action is needed for this gene-drug interaction.",
        },
        {
            "factors": {"UGT1A1": "poor metabolizer"},
            "guideline": "https://pharmgkb.org/guidelineAnnotation/PA166104951",
            "recommendation": "Start with 70% of the standard dose If the patient "
            "tolerates this initial dose, the dose can be increased, "
            "guided by the neutrophil count.",
        },
    ]
    assert DPWG_RECOMMENDATIONS["ribavirin"] == [
        {
            "factors": {},
            "guideline": "https://pharmgkb.org/guidelineAnnotation/PA166104947",
            "recommendation": "NO action is needed for this gene-drug interaction.",
        }
    ]
    assert DPWG_RECOMMENDATIONS["abacavir"] == [
        {
            "factors": {"HLA-B*57:01": "positive"},
            "guideline": "https://pharmgkb.org/guidelineAnnotation/PA166104991",
            "recommendation": "Abacavir is contra-indicated for HLA-B*5701-positive "
            "patients.1. Avoid abacavir.",
        }
    ]
    assert DPWG_RECOMMENDATIONS["carbamazepine"] == []
    # TODO implement https://www.pharmgkb.org/guidelineAnnotation/PA166119846
    # assert DPWG_RECOMMENDATIONS["rasburicase"] == ???


def test_check_cpic_factor_conversion():
    cpic_factors = set(PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC.values())
    dpwg_factors = set(FACTOR_NORMALIZATION.values())

    for drug, recommendations in DPWG_RECOMMENDATIONS.items():
        for recommendation in recommendations:
            for gene, factor in recommendation["factors"].items():
                assert (
                    factor in cpic_factors or factor in dpwg_factors
                ), f"{drug} {gene} {factor}"


# todo: handle pediatric field

# todo: check CYP4F2
