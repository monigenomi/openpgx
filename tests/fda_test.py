from src.openpgx.fda import *
from src.openpgx.helpers import *

FDA_RECOMMENDATIONS = get_fda_recommendations()


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
            "factors": {"HLA-B*57:01": "positive"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher adverse reaction risk (hypersensitivity "
            "reactions). Do not use abacavir in patients positive for "
            "HLA-B*57:01.",
            "strength": "strong",
        }
    ]

    assert FDA_RECOMMENDATIONS["doxepin"] == [
        {
            "factors": {"CYP2C19": "intermediate metabolizer"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {"CYP2C19": "poor metabolizer"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {"CYP2D6": "ultrarapid metabolizer"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "May alter systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {"CYP2D6": "intermediate metabolizer"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "May alter systemic concentrations.",
            "strength": "optional",
        },
        {
            "factors": {"CYP2D6": "poor metabolizer"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "May alter systemic concentrations.",
            "strength": "optional",
        },
    ]

    assert len(FDA_RECOMMENDATIONS) == 104


def test_factors_conversion():
    cpic = set(PHENOTYPE_AND_ALLELE_NORMALIZATIONS_CPIC.values())
    fda = set(VARIANT_MAPPING_FDA.values())

    for drug, recommendations in FDA_RECOMMENDATIONS.items():
        for recommendation in recommendations:
            for gene, factor in recommendation["factors"].items():
                assert factor in cpic or factor in fda, f"{drug} {gene} {factor}"
