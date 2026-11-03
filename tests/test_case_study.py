"""The bundled real-data case study loads, validates and stays scientifically sane.

These tests read only the committed CSVs, so they run in the application's
pinned CI environment without pyDESeq2 or network access.
"""

from __future__ import annotations

import shutil
from dataclasses import fields
from pathlib import Path

import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from plant_expression_explorer.annotations import parse_custom_annotation_table
from plant_expression_explorer.dataset import (
    CASE_STUDY_DIRECTORY,
    CASE_STUDY_SOURCE_LABEL,
    CURRENT_DATASET_KEY,
    DatasetBundle,
    build_dataset_bundle,
    load_case_study_candidate,
)
from plant_expression_explorer.provenance import CASE_STUDY_PROVENANCE
from plant_expression_explorer.validation import IssueCode, Severity

REPOSITORY = Path(__file__).resolve().parents[1]
ANALYSIS = REPOSITORY / "analysis" / "tomato_ja_hairy_roots"


@pytest.fixture(scope="module")
def candidate():
    return load_case_study_candidate()


def test_case_study_tables_pass_application_validation(candidate) -> None:
    assert candidate.status == "valid"
    assert not candidate.report.has_errors
    expression, metadata, de_results = candidate.tables
    assert expression.shape == (19_522, 7)
    assert list(metadata["sample_id"]) == [
        "mock_1", "mock_2", "mock_3", "jasmonate_1", "jasmonate_2", "jasmonate_3",
    ]
    assert list(metadata["condition"]) == ["mock"] * 3 + ["jasmonate"] * 3
    assert metadata["run_accession"].str.fullmatch(r"ERR\d+").all()
    assert de_results is not None
    assert list(de_results["gene_id"]) == list(expression["gene_id"])


def test_case_study_reports_only_expected_non_blocking_issues(candidate) -> None:
    codes = {issue.code for issue in candidate.report.issues}
    # DESeq2 sets p-values to NA for Cook's-distance outliers; padj = 0 is
    # floating-point underflow for extremely strong responses. Both are
    # documented in the case-study README and must stay non-blocking.
    assert codes <= {IssueCode.MISSING_DE_STATISTIC, IssueCode.ZERO_ADJUSTED_P_VALUE}
    assert all(issue.severity is not Severity.ERROR for issue in candidate.report.issues)


def test_case_study_expression_is_log2_normalised_scale(candidate) -> None:
    expression = candidate.tables[0].set_index("gene_id")
    assert (expression >= 0).all().all()
    # log2(count/size_factor + 1) for ~20M-read libraries stays well below 25.
    assert expression.max().max() < 25


def test_jasmonate_positive_controls_are_induced(candidate) -> None:
    """Literature JAZ genes (Chini et al. 2017) must be up in jasmonate."""

    de_results = candidate.tables[2].set_index("gene_id")
    reference = pd.read_csv(ANALYSIS / "reference" / "chini2017_tomato_jaz.csv")
    by_base = {gene_id.rsplit(".", 1)[0]: gene_id for gene_id in de_results.index}
    tested = [by_base[g] for g in reference["gene_id"] if g in by_base]
    assert len(tested) == 11
    lfc = de_results.loc[tested, "log2FoldChange"]
    padj = de_results.loc[tested, "padj"]
    assert (lfc > 0).all()
    assert ((padj < 0.05) & (lfc > 1)).sum() >= 10


def test_case_study_gene_annotation_is_uploadable_and_matches_expression(candidate) -> None:
    table = pd.read_csv(CASE_STUDY_DIRECTORY / "gene_annotation.csv")
    annotations = parse_custom_annotation_table(table)
    expressed = set(candidate.tables[0]["gene_id"])
    assert set(annotations) <= expressed
    families = table["description"].str.split(":").str[0]
    assert {"TPS", "R2R3-MYB", "JAZ"} <= set(families)


def test_case_study_provenance_is_fully_described() -> None:
    for field in fields(CASE_STUDY_PROVENANCE):
        value = getattr(CASE_STUDY_PROVENANCE, field.name)
        assert isinstance(value, str) and value.strip(), field.name
    assert "E-MTAB-11876" in CASE_STUDY_PROVENANCE.dataset_title


def test_case_study_bundle_uses_its_own_source(candidate) -> None:
    bundle = build_dataset_bundle(
        candidate.tables,
        source="case_study",
        source_label=CASE_STUDY_SOURCE_LABEL,
        report=candidate.report,
        provenance=CASE_STUDY_PROVENANCE,
        checksums=candidate.checksums,
    )
    assert bundle.source == "case_study"
    assert bundle.gene_count == 19_522


def test_missing_case_study_file_is_a_clear_blocking_error(tmp_path) -> None:
    copied = tmp_path / "case"
    shutil.copytree(CASE_STUDY_DIRECTORY, copied)
    (copied / "deg_results.csv").unlink()

    result = load_case_study_candidate(copied)

    assert result.status == "invalid"
    (issue,) = result.report.issues
    assert issue.code is IssueCode.DEMO_FILE_MISSING
    assert "run_case_study.py" in issue.message


def test_upload_page_loads_the_real_case_study() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py", default_timeout=60).run()

    app.button(key="pee_load_case_study").click().run()

    assert not app.exception
    bundle = app.session_state[CURRENT_DATASET_KEY]
    assert isinstance(bundle, DatasetBundle)
    assert bundle.source == "case_study"
    assert bundle.provenance == CASE_STUDY_PROVENANCE
    assert bundle.sample_count == 6


def test_differential_expression_page_explains_case_study_statistics(candidate) -> None:
    bundle = build_dataset_bundle(
        candidate.tables,
        source="case_study",
        source_label=CASE_STUDY_SOURCE_LABEL,
        report=candidate.report,
        provenance=CASE_STUDY_PROVENANCE,
        checksums=candidate.checksums,
    )
    app = AppTest.from_file("pages/5_Differential_Expression.py", default_timeout=120)
    app.session_state[CURRENT_DATASET_KEY] = bundle
    app.run()

    assert not app.exception
    info_text = " ".join(str(element.value) for element in app.info)
    assert "real pyDESeq2 results computed offline" in info_text
    assert not any("synthetic" in str(element.value) for element in app.warning)


@pytest.mark.parametrize(
    "page",
    [
        "pages/2_Sample_Quality_Control.py",
        "pages/3_PCA.py",
        "pages/4_Sample_Correlation.py",
        "pages/6_Gene_Expression.py",
        "pages/7_Limitations.py",
        "pages/8_Report_Export.py",
    ],
)
def test_every_analysis_page_runs_on_the_real_case_study(candidate, page) -> None:
    bundle = build_dataset_bundle(
        candidate.tables,
        source="case_study",
        source_label=CASE_STUDY_SOURCE_LABEL,
        report=candidate.report,
        provenance=CASE_STUDY_PROVENANCE,
        checksums=candidate.checksums,
    )
    app = AppTest.from_file(page, default_timeout=120)
    app.session_state[CURRENT_DATASET_KEY] = bundle
    app.run()

    assert not app.exception
