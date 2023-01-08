from openpgx.fda import *

FDA_DATABASE = create_fda_database()
FDA_RECOMMENDATIONS = FDA_DATABASE["recommendations"]
FDA_ENCODINGS = FDA_DATABASE["encodings"]


def test_subgroups_to_factors():
    assert subgroups_to_factors(
        "ultrarapid, normal, intermediate, or poor metabolizers"
    ) == [
        "ultrarapid metabolizer",
        "normal metabolizer",
        "intermediate metabolizer",
        "poor metabolizer",
    ]
    assert subgroups_to_factors("intermediate or poor metabolizers") == [
        "intermediate metabolizer",
        "poor metabolizer",
    ]
    assert subgroups_to_factors("*57:01 allele positive") == ["*57:01 positive"]
    assert subgroups_to_factors("*28/*28 (poor metabolizers)") == ["poor metabolizer"]
    assert subgroups_to_factors(
        "521 TC or 521 CC (intermediate or poor function transporters)"
    ) == ["521 TC", "521 CC", "intermediate function", "poor function"]
    assert subgroups_to_factors("521 CC (poor function transporters)") == [
        "521 CC",
        "poor function",
    ]
    assert subgroups_to_factors("-1639G>A variant carriers") == [
        "rs9923231 reference (C)"
    ]
    assert subgroups_to_factors("Ultrarapid metabolizers") == ["ultrarapid metabolizer"]
    assert subgroups_to_factors("V433M variant carriers") == ["*3 (rs2108622 T, V433M)"]


def test_get_fda_recommendations():
    assert FDA_RECOMMENDATIONS["abacavir"] == [
        {
            "factors": {
                "HLA-B*57:01": "positive",
#                "population": "general"
            },
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher adverse reaction risk (hypersensitivity "
            "reactions). Do not use abacavir in patients positive for "
            "HLA-B*57:01.",
            "strength": "strong",
        }
    ]

    assert FDA_RECOMMENDATIONS["doxepin"] == [
        {
            "factors": {
                "CYP2C19": "intermediate metabolizer",
                # "population": "general"
            },
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {
                "CYP2C19": "poor metabolizer",
                # "population": "general"
            },
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {
                "CYP2D6": "ultrarapid metabolizer",
                # "population": "general"
            },
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "May alter systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {
                "CYP2D6": "intermediate metabolizer",
                # "population": "general"
            },
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "May alter systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {
                "CYP2D6": "poor metabolizer",
                # "population": "general"
            },
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "May alter systemic concentrations.",
            "strength": "optional",
        },
    ]

    assert len(FDA_RECOMMENDATIONS) == 104


def test_fda_encodings_not():
    assert "BCHE" not in FDA_ENCODINGS
    
    assert "NAT2" not in FDA_ENCODINGS
    # https: // www.pharmgkb.org / vip / PA166170337
    

def test_fda_encodings_exists():
    assert FDA_ENCODINGS["HLA-DRB1*07:01"]["positive"] == ["positive"]
    assert FDA_ENCODINGS["SLCO1B1"]["521 TC"] == ["521 TC"]
    assert FDA_ENCODINGS["VKORC1"]["rs9923231 reference (C)"] == ["rs9923231 reference (C)"]

# TODO Check which genes are in FDA but no phenotyping exists in CPIC and DPWG