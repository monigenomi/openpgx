from openpgx.dpwg import *
from openpgx.helpers import *

DATA_DIRECTORY = (
    "../.cache/api.pharmgkb.org/v1/download/file/data/dosingGuidelines.json/"
)
DATA = create_dpwg_database()
DPWG_ENCODINGS = DATA["encodings"]
DPWG_RECOMMENDATIONS = DATA["recommendations"]


def test_extract_gene_name():
    name = "Annotation_of_DPWG_Guideline_for_acenocoumarol_and_VKORC1.json"
    assert extract_gene_name(name) == "VKORC1"


def test_extract_drug_name():
    name = "Annotation_of_DPWG_Guideline_for_acenocoumarol_and_VKORC1.json"
    assert extract_drug_name(name) == "acenocoumarol"


def test_extract_drug_name_contraceptives():
    assert extract_drug_name(
        "Annotation_of_DPWG_Guideline_for_hormonal_contraceptives_for_systemic_use_and_F5")\
           == "hormonal contraceptives for systemic use"


def test_normalize_dpwg_factor():
    assert normalize_dpwg_factor("CYP2D6 UM") == "ultrarapid metabolizer"
    assert normalize_dpwg_factor("HLA-B*57:01") == "*57:01 positive"
    assert normalize_dpwg_factor("HLA-B*44") == "*44 positive"


def test_therapy_table():
    base_dir = path.join(path.dirname(path.realpath(__file__)), DATA_DIRECTORY)
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
    base_dir = path.join(path.dirname(path.realpath(__file__)), DATA_DIRECTORY)
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
            "factors": {
                "UGT1A1": "intermidiate metabolizer",
                # "population": "adults"
            },
            "guideline": "https://www.pharmgkb.org/guidelineAnnotation/PA166104951",
            "recommendation": "NO action is needed for this gene-drug interaction.",
        },
        {
            "factors": {
                "UGT1A1": "poor metabolizer",
                # "population": "adults"
            },
            "guideline": "https://www.pharmgkb.org/guidelineAnnotation/PA166104951",
            "recommendation": "Start with 70% of the standard dose If the patient "
            "tolerates this initial dose, the dose can be increased, "
            "guided by the neutrophil count.",
        },
    ]
    assert DPWG_RECOMMENDATIONS["ribavirin"] == [
        {
            "factors": {
                # "population": "adults"
            },
            "guideline": "https://pharmgkb.org/guidelineAnnotation/PA166104947",
            "recommendation": "Although there is some evidence for lower treatment "
            "response in HLA-B*44 negative patients,  there are no "
            "dosing recommendations for ribavirin at this time.\n",
        }
    ]

    assert DPWG_RECOMMENDATIONS["abacavir"] == [
        {
            "factors": {
                "HLA-B*57:01": "positive",
                # "population": "adults"
            },
            "guideline": "https://www.pharmgkb.org/guidelineAnnotation/PA166104991",
            "recommendation": "Abacavir is contra-indicated for HLA-B*5701-positive "
            "patients.1. Avoid abacavir.",
        }
    ]

def test_issue3_dpwg_f5():
    assert DPWG_RECOMMENDATIONS["hormonal contraceptives for systemic use"] == [{'factors': {'F5': 'Factor V Leiden heterozygous'},
      'guideline': 'https://www.pharmgkb.org/guidelineAnnotation/PA166104955',
      'recommendation': '- If the patient has a FAMILY HISTORY WITH A LOT OF '
                        'THROMBOSIS, or has had a PREVIOUS THROMBOSIS:1. Advise '
                        'the prescriber to avoid the use of contraceptives that '
                        'contain oestrogens and prescribe an on-hormone '
                        'contraceptive-such as a copper IUD - as an alternative. '
                        'One could also opt for a progestogen-only contraceptive '
                        'method, such as the depot injection, an IUD with '
                        'levonorgestrel or an implant with etonogestrel.- OTHER '
                        'CASES:1. Advise the patient to avoid additional risk '
                        'factors for thrombosis (obesity, smoking, etc.).'},
     {'factors': {'F5': 'Factor V Leiden homozygous'},
      'guideline': 'https://www.pharmgkb.org/guidelineAnnotation/PA166104955',
      'recommendation': '- If the patient has a FAMILY HISTORY WITH A LOT OF '
                        'THROMBOSIS, or has had a PREVIOUS THROMBOSIS:1. Advise '
                        'the prescriber to avoid the use of contraceptives that '
                        'contain oestrogens and prescribe a non-hormone '
                        'contraceptive-such as a copper IUD - as an alternative. '
                        'One could also opt for a progestogen-only contraceptive '
                        'method, such as the depot injection, an IUD with '
                        'levonorgestrel or an implant with etonogestrel.- OTHER '
                        'CASES:1. Advise the patient to avoid additional risk '
                        'factors for thrombosis (obesity, smoking, etc.).'}]

def test_issue3_vcorc():
    assert DPWG_RECOMMENDATIONS["acenocoumarol"] == [{'factors': {'VKORC1': 'rs9923231 variant (T)'},
      'guideline': 'https://www.pharmgkb.org/guidelineAnnotation/PA166104938',
      'recommendation': 'Monitoring by the ANTICOAGULATION CLINIC (National INR '
                        'Monitoring Service): recommend to use 50% of the standard '
                        'initial dose. OTHERWISE: recommend to use 50% of the '
                        'standard initial dose and recommend more frequent '
                        'monitoring of the INR .The initial dose and the '
                        'maintenance dose can be calculated using an algorithm. '
                        'However, for patients with two or more VKORC1 and/or '
                        'CYP2C9 variations, the algorithm used in EU-PACT (see '
                        'footnote for a link to a calculation tool in the form of '
                        'an Excel file 1) did not result in a significant '
                        'reduction in the incidence of INRs above the target range '
                        'when compared to an algorithm without genetic '
                        'information. We are therefore unable to recommend the use '
                        'of this algorithm at this time. A (non-validated) '
                        'algorithm has been prescribed for children that should '
                        'result in a better prediction of the maintenance dose for '
                        'AA than the current guideline used by the Anticoagulation '
                        'Clinic [Article:29935043].'},
     {'factors': {'VKORC1': 'rs9923231 reference (C)'},
      'guideline': 'https://www.pharmgkb.org/guidelineAnnotation/PA166104938',
      'recommendation': 'NO action is needed for this gene-drug interaction'},
     {'factors': {},
      'guideline': 'https://pharmgkb.org/guidelineAnnotation/PA166104979',
      'recommendation': 'There are currently no recommendations for acenocoumarol '
                        'dosing based on CYP2C9 genotypes.\n'}]
    
    
    # TODO: assert DPWG_RECOMMENDATIONS["carbamazepine"] == []
    # TODO: assert DPWG_RECOMMENDATIONS["rasburicase"] == [] .... only cpic?
    # TODO: warfarin... it should include CYP4F2 as factor

def test_check_encodings_for_f5():
    assert DPWG_ENCODINGS["F5"] == {
        "Factor V Leiden heterozygous": ["Factor V Leiden heterozygous"],
        "Factor V Leiden homozygous": ["Factor V Leiden homozygous"]
    }

def test_check_encodings_for_vkorc():
    assert DPWG_ENCODINGS["VKORC1"] == {
        'rs9923231 reference (C)': ['rs9923231 reference (C)'],
     'rs9923231 variant (T)': ['rs9923231 variant (T)']}