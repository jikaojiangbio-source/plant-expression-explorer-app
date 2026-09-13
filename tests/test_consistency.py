"""Regression tests for Phase 2 cross-file validation."""

import pandas as pd
from pandas.testing import assert_frame_equal
import pytest

from plant_expression_explorer.consistency import (
    validate_expression_de_consistency,
    validate_expression_metadata_consistency,
    validate_input_tables,
)
from plant_expression_explorer.validation import IssueCode


def _codes(report) -> list[IssueCode]:
    return [issue.code for issue in report.issues]


def test_sample_mismatch_set_differences_have_correct_direction() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0], "s2": [2.0]})
    metadata = pd.DataFrame(
        {"sample_id": ["s2", "s3"], "condition": ["a", "b"]}
    )

    report = validate_expression_metadata_consistency(expression, metadata)

    missing = next(
        issue for issue in report.issues if issue.code is IssueCode.MISSING_METADATA_SAMPLE
    )
    unexpected = next(
        issue
        for issue in report.issues
        if issue.code is IssueCode.UNEXPECTED_METADATA_SAMPLE
    )
    assert missing.example_values == ("s1",)
    assert unexpected.example_values == ("s3",)
    assert report.has_errors


def test_no_sample_overlap_emits_only_specific_error() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    metadata = pd.DataFrame({"sample_id": ["s2"], "condition": ["a"]})

    report = validate_expression_metadata_consistency(expression, metadata)

    assert _codes(report) == [IssueCode.NO_SAMPLE_OVERLAP]
    assert report.has_errors


def test_sample_order_difference_is_information_and_does_not_reorder() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s2": [1.0], "s1": [2.0]})
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["a", "b"]}
    )
    original = metadata.copy(deep=True)

    report = validate_expression_metadata_consistency(expression, metadata)

    assert not report.has_errors
    assert _codes(report) == [IssueCode.SAMPLE_ORDER_DIFFERS]
    assert_frame_equal(metadata, original)


def test_partial_gene_mismatches_warn_in_both_directions() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "s1": [1.0, 2.0]}
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g2", "g3"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.1, 0.2],
            "padj": [0.2, 0.3],
        }
    )

    report = validate_expression_de_consistency(expression, de_results)

    assert not report.has_errors
    assert _codes(report) == [
        IssueCode.DE_GENE_NOT_IN_EXPRESSION,
        IssueCode.EXPRESSION_GENE_NOT_IN_DE,
    ]
    assert report.warnings[0].example_values == ("g3",)
    assert report.warnings[1].example_values == ("g1",)
    for issue in report.warnings:
        assert "annotation releases" in issue.message
        assert "gene-versus-transcript" in issue.message
        assert "did not rewrite any identifier" in issue.message


def test_no_gene_overlap_emits_only_specific_blocking_error() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    de_results = pd.DataFrame(
        {
            "gene_id": ["g2"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )

    report = validate_expression_de_consistency(expression, de_results)

    assert _codes(report) == [IssueCode.NO_GENE_OVERLAP]
    assert report.has_errors
    assert "isoform or version suffixes" in report.errors[0].message
    assert "did not rewrite any identifier" in report.errors[0].message


def test_input_validation_without_de_is_informational_and_non_mutating() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "s1": [1.0, 2.0], "s2": [2.0, 3.0]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["control", "treated"]}
    )
    expression_original = expression.copy(deep=True)
    metadata_original = metadata.copy(deep=True)

    report = validate_input_tables(expression, metadata, None)

    assert not report.has_errors
    assert _codes(report) == [
        IssueCode.CONDITION_WITH_SINGLE_SAMPLE,
        IssueCode.DE_RESULTS_NOT_SUPPLIED,
    ]
    assert report.information[0].table == "Differential-expression results"
    assert "remain available" in report.information[0].message
    assert_frame_equal(expression, expression_original)
    assert_frame_equal(metadata, metadata_original)


def test_input_validation_aggregates_table_and_cross_file_reports() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "s1": [1.0, 2.0], "s2": [2.0, 3.0]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["control", "control"]}
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g3"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.1, 0.2],
            "padj": [0.2, None],
        }
    )

    report = validate_input_tables(expression, metadata, de_results)

    assert not report.has_errors
    assert IssueCode.SINGLE_CONDITION in _codes(report)
    assert IssueCode.MISSING_DE_STATISTIC in _codes(report)
    assert IssueCode.DE_GENE_NOT_IN_EXPRESSION in _codes(report)
    assert IssueCode.EXPRESSION_GENE_NOT_IN_DE in _codes(report)


def test_input_validation_blocks_single_sample_expression() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    metadata = pd.DataFrame({"sample_id": ["s1"], "condition": ["control"]})
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )

    report = validate_input_tables(expression, metadata, de_results)

    assert IssueCode.LOW_SAMPLE_COUNT in [issue.code for issue in report.errors]
    assert IssueCode.LOW_SAMPLE_COUNT not in [
        issue.code for issue in report.warnings
    ]
    assert report.has_errors is True


@pytest.mark.parametrize(
    "duplicate_table", ["expression", "metadata", "de_results"]
)
def test_input_validation_blocks_duplicate_identifier_in_each_table(
    duplicate_table: str,
) -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "s1": [1.0, 2.0],
            "s2": [2.0, 3.0],
            "s3": [3.0, 4.0],
            "s4": [4.0, 5.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3", "s4"],
            "condition": ["control", "control", "treated", "treated"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.1, 0.2],
            "padj": [0.2, 0.3],
        }
    )

    if duplicate_table == "expression":
        expression = pd.concat([expression, expression.iloc[[0]]], ignore_index=True)
    elif duplicate_table == "metadata":
        metadata = pd.concat([metadata, metadata.iloc[[0]]], ignore_index=True)
    else:
        de_results = pd.concat([de_results, de_results.iloc[[0]]], ignore_index=True)

    report = validate_input_tables(expression, metadata, de_results)

    assert IssueCode.DUPLICATE_IDENTIFIER in [
        issue.code for issue in report.errors
    ]
    assert report.has_errors is True


def test_cross_file_validators_preserve_all_dataframes() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s2": [1.0], "s1": [2.0]})
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["a", "b"]}
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )
    originals = [table.copy(deep=True) for table in (expression, metadata, de_results)]

    validate_input_tables(expression, metadata, de_results)

    for actual, expected in zip((expression, metadata, de_results), originals):
        assert_frame_equal(actual, expected)


def test_input_validation_reports_duplicate_required_columns_without_crashing() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    metadata = pd.DataFrame(
        [["s1", "s1", "control"]],
        columns=["sample_id", "sample_id", "condition"],
    )
    de_results = pd.DataFrame(
        [["g1", "g1", 1.0, 0.1, 0.2]],
        columns=["gene_id", "gene_id", "log2FoldChange", "pvalue", "padj"],
    )

    report = validate_input_tables(expression, metadata, de_results)

    assert report.has_errors
    assert _codes(report).count(IssueCode.DUPLICATE_COLUMN_NAME) == 2
