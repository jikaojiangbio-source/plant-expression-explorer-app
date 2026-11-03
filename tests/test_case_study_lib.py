"""Unit tests for the pure helpers behind the tomato jasmonate case study."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "analysis" / "tomato_ja_hairy_roots"))

import case_study_lib as lib  # noqa: E402

CONFIG = """<?xml version="1.0"?>
<configuration experimentType="rnaseq_mrna_differential">
  <analytics>
    <assay_groups>
      <assay_group id="g1" label="jasmonate; 50 micromolar">
        <assay>ERR3</assay><assay>ERR1</assay><assay>ERR2</assay>
      </assay_group>
      <assay_group id="g2" label="none">
        <assay>ERR6</assay><assay>ERR4</assay><assay>ERR5</assay>
      </assay_group>
    </assay_groups>
    <contrasts>
      <contrast id="g2_g1">
        <reference_assay_group>g2</reference_assay_group>
        <test_assay_group>g1</test_assay_group>
      </contrast>
    </contrasts>
  </analytics>
</configuration>"""


# --- design -----------------------------------------------------------------


def test_parse_atlas_configuration_reads_reference_and_test_groups() -> None:
    contrast = lib.parse_atlas_configuration(CONFIG)
    assert contrast.reference_label == "none"
    assert contrast.reference_assays == ("ERR4", "ERR5", "ERR6")
    assert contrast.test_assays == ("ERR1", "ERR2", "ERR3")


def test_parse_atlas_configuration_rejects_assay_in_both_groups() -> None:
    broken = CONFIG.replace("<assay>ERR6</assay>", "<assay>ERR1</assay>")
    with pytest.raises(ValueError, match="both groups"):
        lib.parse_atlas_configuration(broken)


def test_parse_atlas_configuration_requires_exactly_one_contrast() -> None:
    doubled = CONFIG.replace("</contrasts>", '<contrast id="x"><reference_assay_group>g1</reference_assay_group><test_assay_group>g2</test_assay_group></contrast></contrasts>')
    with pytest.raises(ValueError, match="exactly one contrast"):
        lib.parse_atlas_configuration(doubled)


def test_build_sample_metadata_puts_reference_first_and_keeps_accessions() -> None:
    metadata = lib.build_sample_metadata(lib.parse_atlas_configuration(CONFIG))
    assert list(metadata["sample_id"]) == [
        "mock_1", "mock_2", "mock_3", "jasmonate_1", "jasmonate_2", "jasmonate_3",
    ]
    assert list(metadata["run_accession"]) == ["ERR4", "ERR5", "ERR6", "ERR1", "ERR2", "ERR3"]


# --- counts -----------------------------------------------------------------


def _raw_counts() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Gene ID": ["ENSRNA1", "Solyc01g000010.1", "Solyc01g000020.2"],
            "Gene Name": ["tRNA", None, "JAZ2"],
            "ERR1": [5, 10, 100], "ERR2": [5, 12, 90], "ERR3": [5, 8, 110],
            "ERR4": [5, 20, 10], "ERR5": [5, 22, 12], "ERR6": [5, 18, 8],
        }
    )


def test_prepare_counts_keeps_nuclear_genes_and_renames_samples() -> None:
    metadata = lib.build_sample_metadata(lib.parse_atlas_configuration(CONFIG))
    counts = lib.prepare_counts(_raw_counts(), metadata)
    assert list(counts.index) == ["Solyc01g000010.1", "Solyc01g000020.2"]
    assert list(counts.columns) == list(metadata["sample_id"])
    assert counts.loc["Solyc01g000020.2", "jasmonate_1"] == 100


def test_prepare_counts_rejects_missing_runs_and_negative_values() -> None:
    metadata = lib.build_sample_metadata(lib.parse_atlas_configuration(CONFIG))
    with pytest.raises(ValueError, match="missing"):
        lib.prepare_counts(_raw_counts().drop(columns="ERR5"), metadata)
    negative = _raw_counts()
    negative.loc[1, "ERR1"] = -1
    with pytest.raises(ValueError, match="non-negative"):
        lib.prepare_counts(negative, metadata)


def test_prepare_counts_rejects_fractional_counts() -> None:
    metadata = lib.build_sample_metadata(lib.parse_atlas_configuration(CONFIG))
    fractional = _raw_counts().astype({"ERR1": float})
    fractional.loc[1, "ERR1"] = 10.5
    with pytest.raises(ValueError, match="integers"):
        lib.prepare_counts(fractional, metadata)


def test_filter_low_counts_keeps_genes_expressed_in_one_condition_only() -> None:
    counts = pd.DataFrame(
        {"a": [0, 50, 9], "b": [0, 50, 9], "c": [0, 50, 9], "d": [12, 0, 9], "e": [15, 0, 9], "f": [11, 0, 9]},
        index=["only_test", "only_ref", "everywhere_low"],
    )
    kept = lib.filter_low_counts(counts, min_count=10, min_samples=3)
    assert list(kept.index) == ["only_test", "only_ref"]


def test_size_factors_match_hand_calculation_and_scale_with_depth() -> None:
    counts = pd.DataFrame({"s1": [10, 20, 30], "s2": [20, 40, 60]}, index=list("abc"))
    factors = lib.size_factors_median_of_ratios(counts)
    # geometric mean per gene = x * sqrt(2); ratios are 1/sqrt(2) and sqrt(2).
    assert factors["s1"] == pytest.approx(1 / math.sqrt(2))
    assert factors["s2"] == pytest.approx(math.sqrt(2))
    assert factors["s2"] / factors["s1"] == pytest.approx(2.0)


def test_size_factors_ignore_genes_with_a_zero() -> None:
    counts = pd.DataFrame({"s1": [10, 0], "s2": [20, 999]}, index=["a", "zero"])
    factors = lib.size_factors_median_of_ratios(counts)
    assert factors["s2"] / factors["s1"] == pytest.approx(2.0)


def test_log2_normalised_expression_formula() -> None:
    counts = pd.DataFrame({"s1": [0, 3]}, index=pd.Index(["g0", "g3"], name="gene_id"))
    table = lib.log2_normalised_expression(counts, pd.Series({"s1": 0.5}))
    assert list(table.columns) == ["gene_id", "s1"]
    assert table["s1"].tolist() == [0.0, pytest.approx(math.log2(7))]


# --- gene families ------------------------------------------------------------


def _hits(*rows: tuple[str, str, str]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["gene_id", "protein_id", "pfam"])


def test_myb_repeat_class_names() -> None:
    assert lib.myb_repeat_class(1) == "1R-MYB (MYB-related)"
    assert lib.myb_repeat_class(2) == "R2R3-MYB"
    assert lib.myb_repeat_class(3) == "3R-MYB"
    with pytest.raises(ValueError):
        lib.myb_repeat_class(0)


def test_classify_families_calls_each_family_by_its_rule() -> None:
    hits = _hits(
        ("tps_full", "tps_full.1", lib.PFAM_TERPENE_SYNTH_N),
        ("tps_full", "tps_full.1", lib.PFAM_TERPENE_SYNTH_C),
        ("tps_frag", "tps_frag.1", lib.PFAM_TERPENE_SYNTH_C),
        ("myb2", "myb2.1", lib.PFAM_MYB),
        ("myb2", "myb2.1", lib.PFAM_MYB),
        ("myb_iso", "myb_iso.1", lib.PFAM_MYB),
        ("myb_iso", "myb_iso.2", lib.PFAM_MYB),
        ("myb_iso", "myb_iso.2", lib.PFAM_MYB),
        ("myb_iso", "myb_iso.2", lib.PFAM_MYB),
        ("jaz", "jaz.1", lib.PFAM_TIFY),
        ("jaz", "jaz.1", lib.PFAM_JAS),
        ("tify_only", "tify_only.1", lib.PFAM_TIFY),
        ("split", "split.1", lib.PFAM_TIFY),
        ("split", "split.2", lib.PFAM_JAS),
    )
    calls = lib.classify_families(hits).set_index("gene_id")

    assert calls.loc["tps_full", "family"] == "TPS"
    assert "N-terminal + metal-binding" in calls.loc["tps_full", "evidence"]
    assert "possibly partial" in calls.loc["tps_frag", "evidence"]
    assert calls.loc["myb2", "family"] == "R2R3-MYB"
    # The gene takes the isoform with the most repeats (3), not the first one.
    assert calls.loc["myb_iso", "family"] == "3R-MYB"
    assert calls.loc["jaz", "family"] == "JAZ"
    # tify alone is not JAZ, and the two domains must be on the same isoform.
    assert "tify_only" not in calls.index
    assert "split" not in calls.index


def test_classify_families_requires_columns() -> None:
    with pytest.raises(ValueError, match="columns"):
        lib.classify_families(pd.DataFrame({"gene_id": ["a"]}))


# --- differential-expression tables -------------------------------------------


def test_app_de_table_keeps_missing_values_and_four_columns() -> None:
    results = pd.DataFrame(
        {
            "log2FoldChange": [1.123456789123, -2.0, 0.1],
            "pvalue": [1e-200, np.nan, 0.5],
            "padj": [3.14159265e-180, np.nan, 1.0],
            "baseMean": [1, 2, 3],
        },
        index=["a", "b", "c"],
    )
    table = lib.app_de_table(results)
    assert list(table.columns) == ["gene_id", "log2FoldChange", "pvalue", "padj"]
    assert table.loc[0, "padj"] == pytest.approx(3.14159e-180)
    assert math.isnan(table.loc[1, "pvalue"]) and math.isnan(table.loc[1, "padj"])
    assert table["padj"].dropna().between(0, 1).all()


def test_classify_regulation_uses_padj_and_sign_only() -> None:
    padj = pd.Series([0.01, 0.01, 0.2, np.nan], index=list("abcd"))
    lfc = pd.Series([0.4, -3.0, 5.0, 1.0], index=list("abcd"))
    labels = lib.classify_regulation(padj, lfc)
    assert labels.tolist() == ["up", "down", "not significant", "not tested"]


def test_hypergeometric_enrichment_matches_exact_small_case() -> None:
    # Urn of 10 genes, 4 induced; draw a family of 3; P(at least 2 induced).
    expected = (math.comb(4, 2) * math.comb(6, 1) + math.comb(4, 3)) / math.comb(10, 3)
    assert lib.hypergeometric_enrichment(2, 3, 4, 10) == pytest.approx(expected)
    assert lib.hypergeometric_enrichment(0, 3, 4, 10) == pytest.approx(1.0)


def test_hypergeometric_enrichment_rejects_impossible_counts() -> None:
    with pytest.raises(ValueError):
        lib.hypergeometric_enrichment(5, 3, 4, 10)


def test_family_response_summary_counts_by_class() -> None:
    regulation = pd.Series({"t1": "up", "t2": "not significant", "m1": "down"})
    families = pd.DataFrame({"gene_id": ["t1", "t2", "t3", "m1"], "family": ["TPS", "TPS", "TPS", "R2R3-MYB"], "evidence": ""})
    summary = lib.family_response_summary(regulation, families, {"TPS": 3, "R2R3-MYB": 1}).set_index("family")
    assert summary.loc["TPS", "in_proteome"] == 3
    assert summary.loc["TPS", "expressed_and_tested"] == 2
    assert summary.loc["TPS", "up"] == 1
    assert summary.loc["R2R3-MYB", "down"] == 1


# --- identifiers and reference sets -------------------------------------------


def test_strip_version_and_unversioned_matching() -> None:
    assert lib.strip_version("Solyc12g009220.2") == "Solyc12g009220"
    assert lib.strip_version("ENSRNA049444660") == "ENSRNA049444660"
    mapping = lib.match_unversioned_ids(["Solyc12g009220", "Solyc99g000000"], ["Solyc12g009220.2", "Solyc01g1.1"])
    assert mapping == {"Solyc12g009220": "Solyc12g009220.2", "Solyc99g000000": None}
    with pytest.raises(ValueError, match="several"):
        lib.match_unversioned_ids(["Solyc1"], ["Solyc1.1", "Solyc1.2"])


def test_reference_response_table_reports_filtered_genes_instead_of_dropping_them() -> None:
    metadata = pd.DataFrame({"sample_id": ["m1", "j1"], "condition": ["mock", "jasmonate"]})
    counts_all = pd.DataFrame({"m1": [10, 0], "j1": [80, 3]}, index=["G1.1", "G2.1"])
    results = pd.DataFrame({"log2FoldChange": [3.0], "padj_lfc1": [1e-5], "regulation_lfc1": ["up"]}, index=["G1.1"])
    reference = pd.DataFrame({"symbol": ["A", "B", "C"], "gene_id": ["G1", "G2", "G3"]})
    table = lib.reference_response_table(reference, counts_all, metadata, results, ["G1.1"]).set_index("symbol")
    assert table.loc["A", "status"] == "up" and bool(table.loc["A", "in_domain_call"])
    assert table.loc["B", "status"] == "below expression filter"
    assert table.loc["B", "mean_count_jasmonate"] == 3
    assert table.loc["C", "status"] == "absent from annotation"


def test_compare_with_atlas_like_for_like() -> None:
    ours = pd.DataFrame({"log2FoldChange": [2.0, -1.5, 0.1, 3.0], "padj": [0.001, 0.01, 0.9, 0.2]}, index=list("abcd"))
    atlas = pd.DataFrame({"Gene ID": list("abce"), "g2_g1.p-value": [0.001, 0.5, 0.9, 0.0], "g2_g1.log2foldchange": [2.0, -1.4, 0.1, 9.0]})
    summary = lib.compare_with_atlas(ours, atlas)
    assert summary["shared_genes"] == 3
    assert summary["deg_ours"] == 2 and summary["deg_atlas"] == 1 and summary["deg_both"] == 1
    assert summary["deg_jaccard"] == pytest.approx(0.5)


def test_first_present_skips_missing_and_blank() -> None:
    assert lib.first_present([None, float("nan"), "  ", "JAZ2"]) == "JAZ2"
    assert lib.first_present([None]) is None


def test_atlas_pvalue_scale_check_recognises_an_adjusted_column() -> None:
    ours = pd.DataFrame({"pvalue": [1e-6, 1e-3, 0.01], "padj": [1e-5, 1e-2, 0.05]}, index=list("abc"))
    atlas = pd.DataFrame({"Gene ID": list("abc"), "g2_g1.p-value": [1e-5, 1e-2, 0.05]})
    check = lib.atlas_pvalue_scale_check(ours, atlas)
    assert check["median_abs_log10_diff_vs_padj"] == 0
    assert check["median_abs_log10_diff_vs_pvalue"] > 0.5
    assert check["fraction_atlas_ge_pvalue"] == 1.0
