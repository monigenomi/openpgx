from src.openpgx.openpgx import *


def test_get_all_drugs():
    drugs = get_all_drugs()
    assert len(set(["doxepin", "abacavir", "rasburicase"]) & drugs) == 3
    assert len(drugs) == 139


def test_get_recommendations_for_drug_HLA_ABACAVIR():
    # Abacavir exists in every base with different recommendation in each
    assert (
        sorted(
            list(
                get_recommendations_for_drug(
                    "abacavir",
                    {
                        "HLA-B*57:01": {
                            "factor": "positive",
                            "cpic_factor": "positive",
                            "activityscore": None,
                        }
                    },
                ).keys()
            )
        )
        == ["cpic", "dpwg", "fda"]
    )


def test_get_recommendation_for_drug():
    assert get_recommendations_for_drug("allopurinol", {}) == {
        "cpic": {
            "factors": {},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
            "recommendation": "Recommendations are available, but they require "
            "genotypes of following genes: HLA-B*58:01",
        },
        "fda": {
            "factors": {},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Recommendations are available, but they require "
            "genotypes of following genes: HLA-B*58:01",
        },
    }
    assert get_recommendations_for_drug(
        "allopurinol", phenoconversion({"HLA-B*58:01": "positive"})
    ) == {
        "cpic": {
            "factors": {"HLA-B*58:01": "positive"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
            "recommendation": "Allopurinol is contraindicated",
            "strength": "strong",
        },
        "fda": {
            "factors": {"HLA-B*58:01": "positive"},
            "guideline": "https://www.fda.gov/medical-devices/precision-medicine/table-pharmacogenetic-associations",
            "recommendation": "Results in higher adverse reaction risk (severe "
            "skin reactions).",
            "strength": "moderate",
        },
    }
    assert get_recommendations_for_drug(
        "allopurinol", phenoconversion({"HLA-B*58:01": "negative"})
    ) == {
        "cpic": {
            "factors": {"HLA-B*58:01": "negative"},
            "guideline": "https://cpicpgx.org/guidelines/guideline-for-allopurinol-and-hla-b/",
            "recommendation": "Use allopurinol per standard dosing guidelines",
            "strength": "strong",
        }
    }


def test_phenoconversion():
    # example of input and output to this function
    assert phenoconversion({"CYP2D6": "*2≥3/*1≥3"}) == {
        "CYP2D6": {
            "activityscore": 6.0,
            "cpic_factor": "Ultrarapid Metabolizer",
            "factor": "ultrarapid metabolizer",
        }
    }

    def short_(gene, genotype):
        result = phenoconversion({gene: genotype})[gene]
        return result["activityscore"], result["cpic_factor"], result["factor"]

    # "normal genes"
    assert short_("G6PD", "B (wildtype)") == (None, "Normal", "normal")
    assert short_("CYP2D6", "*7/*7") == (0.00, "Poor Metabolizer", "poor metabolizer")
    assert short_("G6PD", "B (wildtype)") == (None, "Normal", "normal")
    assert short_("TPMT", "*4/*10") == (
        None,
        "Possible Intermediate Metabolizer",
        "intermediate metabolizer",
    )

    # Activity score obligatory for DPYD
    # TODO: index is not sorted properly: 'DPYD:Reference/c.1898delC (*3)' why result is empty
    assert short_("DPYD", "c.1898delC (*3)/Reference") == (
        1.0,
        "Intermediate Metabolizer",
        "intermediate metabolizer",
    )
    assert short_("DPYD", "Reference/c.1905+1G>A") == (None, None, None)

    # Phenotype in CPIC
    assert short_("SLCO1B1", "*1A/*1B") == (None, "Normal Function", "normal function")

    # CPIC Phenoconversion impossible (allele does not exists) but 521 CC exists in dpwg and fda both
    assert short_("SLCO1B1", "521 CC") == (None, None, "521 CC")

    # Allele based recommendation - "rs9923231 reference (C)" exists in CPIC but no phenoconversion possible
    assert short_("VKORC1", "rs9923231 reference (C)") == (
        None,
        None,
        "rs9923231 reference (C)",
    )

    # HLA-B*44: Gene and allele exists only in dpwg but recommendation is default whether genotype is given or not
    assert short_("HLA-B*44", "positive") == (None, None, None)

    # "negative and positive is always converted to dpwg and fda for the same name
    assert short_("HLA-B*57:01", "negative") == (None, "negative", "negative")
    assert short_("HLA-B*57:01", "positive") == (None, "positive", "positive")

    # If allele or gene does not exists return None everywhere
    assert short_("CYP2D6", "*150/*190") == (None, None, None)

    # TODO: What to do with not existing gene?
    assert short_("FOO", "bar") == (None, None, None)

    assert short_("CYP2D6", "*104/*1x5") == (None, "Indeterminate", None)
    assert short_("DPYD", "Reference/c.1905+1G>A (*2A)") == (
        1.0,
        "Intermediate Metabolizer",
        "intermediate metabolizer",
    )


def test_get_recommendations_CYPS():
    # Test if "multiple gene" factors works
    recommendations = get_recommendations({"CYP2D6": "*7/*7", "CYP2C19": "*1/*2"})

    assert recommendations["trimipramine"]["cpic"] == {
        "factors": {"CYP2D6": "== 0.00", "CYP2C19": "Intermediate Metabolizer"},
        "recommendation": "Avoid trimipramine use. If a trimipramine is warranted, consider a 50% reduction of recommended starting dose. Utilizing therapeutic drug monitoring to guide dose adjustments is strongly recommended.",
        "strength": "optional",
        "guideline": "https://cpicpgx.org/guidelines/guideline-for-tricyclic-antidepressants-and-cyp2d6-and-cyp2c19/",
    }


def test_get_recommendations_with_multiple_factors():
    recommendations = get_recommendations(
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
    recommendations = get_recommendations({"DPYD": "c.601A>C/c.2194G>A (*6)"})
    assert recommendations["capecitabine"]["dpwg"]["factors"] == {"DPYD": "== 1.00"}


def test_compare_activity_score():
    for n in [2.0, 2]:
        assert compare_factor("Normal Metabolizer", n, ">= 1.5") == True
        assert compare_factor("Normal Metabolizer", n, ">= 2.0") == True
        assert compare_factor("Normal Metabolizer", n, ">= 2") == True
        assert compare_factor("Normal Metabolizer", n, ">= 2.5") == False
        assert compare_factor("Normal Metabolizer", n, "== 2") == True
        assert compare_factor("Normal Metabolizer", n, "== 2.0") == True

    compare_factor("Normal Metabolizer", 1.0, "Normal Metabolizer") == True
    compare_factor("*57:01 negative", 1.0, "*57:01 negative") == True
    compare_factor("Ultra Metabolizer", 1.0, "Normal Metabolizer") == False
    compare_factor(None, 1.0, None) == True
    compare_factor("Ultra Metabolizer", None, None) == False
    compare_factor(None, 1.0, "Ultra Metabolizer") == False


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
