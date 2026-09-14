"""Regression tests for Phase 2 single-table validation."""

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest

from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationReport,
    validate_de_results,
    validate_expression_matrix,
    validate_sample_metadata,
)


def _codes(report: ValidationReport) -> list[IssueCode]:
    return [issue.code for issue in report.issues]


def test_report_separates_errors_warnings_and_information() -> None:
    expression = pd.DataFrame(
        {"gene_id": [" g1 "], "sample1": ["-1.5"], "sample2": [0.0]}
    )

    report = validate_expression_matrix(expression)

    assert report.has_errors
    assert {issue.code for issue in report.errors} == {
        IssueCode.WHITESPACE_IN_IDENTIFIER
    }
    assert {issue.code for issue in report.warnings} == {
        IssueCode.NEGATIVE_EXPRESSION_VALUE,
    }
    assert {issue.code for issue in report.information} == {
        IssueCode.COERCIBLE_NUMERIC_STRING,
    }
    assert all(isinstance(issue.code, IssueCode) for issue in report.issues)
    assert all(isinstance(issue.severity, Severity) for issue in report.issues)


@pytest.mark.parametrize(
    ("table", "expected_code"),
    [
        (pd.DataFrame(), IssueCode.EMPTY_TABLE),
        (pd.DataFrame({"sample1": [1.0]}), IssueCode.MISSING_REQUIRED_COLUMN),
    ],
)
def test_expression_structural_errors_block(
    table: pd.DataFrame, expected_code: IssueCode
) -> None:
    report = validate_expression_matrix(table)

    assert report.has_errors
    assert expected_code in [issue.code for issue in report.errors]


def test_expression_rejects_zero_sample_columns() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"]})

    report = validate_expression_matrix(expression)

    assert IssueCode.NO_SAMPLE_COLUMNS in [issue.code for issue in report.errors]
    assert report.has_errors is True


def test_expression_rejects_single_sample_matrix_without_mutating_it() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "sample1": [1.0]})
    original = expression.copy(deep=True)

    assert len([column for column in expression.columns if column != "gene_id"]) == 1

    report = validate_expression_matrix(expression)

    assert IssueCode.LOW_SAMPLE_COUNT in [issue.code for issue in report.errors]
    assert IssueCode.LOW_SAMPLE_COUNT not in [
        issue.code for issue in report.warnings
    ]
    assert report.has_errors is True
    assert_frame_equal(expression, original)


def test_expression_accepts_exactly_two_sample_columns() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "sample1": [1.0], "sample2": [2.0]}
    )

    report = validate_expression_matrix(expression)

    assert IssueCode.NO_SAMPLE_COLUMNS not in _codes(report)
    assert IssueCode.LOW_SAMPLE_COUNT not in _codes(report)
    assert not report.has_errors


def test_expression_rejects_duplicate_sample_columns() -> None:
    expression = pd.DataFrame(
        [["g1", 1.0, 2.0]], columns=["gene_id", "sample1", "sample1"]
    )

    report = validate_expression_matrix(expression)

    assert report.has_errors
    assert IssueCode.DUPLICATE_COLUMN_NAME in _codes(report)


def test_expression_distinguishes_numeric_text_invalid_text_and_missing_values() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample1": ["1.5", "high", None],
        }
    )

    report = validate_expression_matrix(expression)

    assert IssueCode.COERCIBLE_NUMERIC_STRING in _codes(report)
    assert IssueCode.NON_NUMERIC_VALUE in _codes(report)
    assert IssueCode.MISSING_EXPRESSION_VALUE in _codes(report)


def test_expression_flags_pandas_auto_named_blank_header_as_empty_column_name() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "sample1": [1.0],
            "Unnamed: 2": [np.nan],
        }
    )

    report = validate_expression_matrix(expression)

    empty_name_issues = [
        issue for issue in report.issues if issue.code is IssueCode.EMPTY_COLUMN_NAME
    ]
    assert len(empty_name_issues) == 1
    assert "Unnamed: 2" in empty_name_issues[0].message
    assert "trailing empty column" in empty_name_issues[0].message


def test_expression_hints_at_missing_value_placeholders() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample1": [1.0, 2.0], "sample2": ["-", "3.0"]}
    )

    report = validate_expression_matrix(expression)

    non_numeric = next(
        issue for issue in report.issues if issue.code is IssueCode.NON_NUMERIC_VALUE
    )
    assert "leave the cell blank" in non_numeric.message


def test_expression_does_not_hint_at_missing_value_placeholders_for_ordinary_text() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "sample1": [1.0], "sample2": ["not_a_number"]}
    )

    report = validate_expression_matrix(expression)

    non_numeric = next(
        issue for issue in report.issues if issue.code is IssueCode.NON_NUMERIC_VALUE
    )
    assert "leave the cell blank" not in non_numeric.message


def test_expression_accepts_negative_values_as_warning() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "sample1": [-1.5], "sample2": [0.0]}
    )

    report = validate_expression_matrix(expression)

    assert not report.has_errors
    assert [issue.code for issue in report.warnings] == [
        IssueCode.NEGATIVE_EXPRESSION_VALUE
    ]
    assert "transformed or centred expression data" in report.warnings[0].message
    assert "confirm that the matrix data type is appropriate" in report.warnings[0].message


def test_expression_rejects_boolean_values_as_non_numeric() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "sample1": [True]})

    report = validate_expression_matrix(expression)

    assert IssueCode.NON_NUMERIC_VALUE in _codes(report)
    assert report.has_errors


def test_expression_rejects_identifier_and_sample_name_whitespace() -> None:
    expression = pd.DataFrame({"gene_id": [" g1 "], " sample1 ": [1.0]})

    report = validate_expression_matrix(expression)

    assert IssueCode.WHITESPACE_IN_IDENTIFIER in _codes(report)
    assert IssueCode.WHITESPACE_IN_SAMPLE_NAME in _codes(report)
    assert report.has_errors


def test_expression_rejects_missing_and_duplicate_gene_ids() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g1", None],
            "sample1": [1.0, 2.0, 3.0],
            "sample2": [2.0, 3.0, 4.0],
        }
    )

    report = validate_expression_matrix(expression)

    assert IssueCode.MISSING_IDENTIFIER in _codes(report)
    assert IssueCode.DUPLICATE_IDENTIFIER in _codes(report)
    assert report.has_errors


def test_expression_duplicate_gene_id_is_error_without_mutating_input() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g1"],
            "sample1": [1.0, 2.0],
            "sample2": [2.0, 3.0],
        }
    )
    original = expression.copy(deep=True)

    report = validate_expression_matrix(expression)

    assert IssueCode.DUPLICATE_IDENTIFIER in [
        issue.code for issue in report.errors
    ]
    assert IssueCode.DUPLICATE_IDENTIFIER not in [
        issue.code for issue in report.warnings
    ]
    assert IssueCode.DUPLICATE_IDENTIFIER not in [
        issue.code for issue in report.information
    ]
    assert report.has_errors is True
    assert_frame_equal(expression, original)


def test_expression_rejects_infinite_values() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "sample1": [float("inf")], "sample2": [1.0]}
    )

    report = validate_expression_matrix(expression)

    assert IssueCode.NON_FINITE_VALUE in _codes(report)
    assert report.has_errors


def test_expression_phase2_does_not_flag_constant_or_zero_variance_values() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample1": [1.0, 1.0],
            "sample2": [1.0, 1.0],
        }
    )

    report = validate_expression_matrix(expression)

    assert report.issues == ()


def test_metadata_reports_condition_limitations_without_blocking() -> None:
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["control", "control"]}
    )

    report = validate_sample_metadata(metadata)

    assert not report.has_errors
    assert [issue.code for issue in report.warnings] == [IssueCode.SINGLE_CONDITION]


def test_metadata_warns_for_single_sample_conditions() -> None:
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["control", "treated"]}
    )

    report = validate_sample_metadata(metadata)

    assert not report.has_errors
    assert IssueCode.CONDITION_WITH_SINGLE_SAMPLE in _codes(report)


def test_metadata_rejects_identifier_and_condition_whitespace() -> None:
    metadata = pd.DataFrame(
        {"sample_id": [" s1 "], "condition": [" control "]}
    )

    report = validate_sample_metadata(metadata)

    assert IssueCode.WHITESPACE_IN_IDENTIFIER in _codes(report)
    assert IssueCode.WHITESPACE_IN_CATEGORY_VALUE in _codes(report)
    assert report.has_errors


@pytest.mark.parametrize(
    "metadata",
    [
        pd.DataFrame(),
        pd.DataFrame({"condition": ["control"]}),
        pd.DataFrame({"sample_id": ["s1"]}),
        pd.DataFrame({"sample_id": [None], "condition": ["control"]}),
        pd.DataFrame(
            {"sample_id": ["s1", "s1"], "condition": ["a", "b"]}
        ),
        pd.DataFrame({"sample_id": ["s1"], "condition": [None]}),
    ],
)
def test_metadata_contract_errors_block(metadata: pd.DataFrame) -> None:
    report = validate_sample_metadata(metadata)

    assert report.has_errors


def test_metadata_duplicate_sample_id_is_error_without_mutating_input() -> None:
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s1", "s2", "s3"],
            "condition": ["control", "control", "treated", "treated"],
        }
    )
    original = metadata.copy(deep=True)

    report = validate_sample_metadata(metadata)

    assert IssueCode.DUPLICATE_IDENTIFIER in [
        issue.code for issue in report.errors
    ]
    assert IssueCode.DUPLICATE_IDENTIFIER not in [
        issue.code for issue in report.warnings
    ]
    assert IssueCode.DUPLICATE_IDENTIFIER not in [
        issue.code for issue in report.information
    ]
    assert report.has_errors is True
    assert_frame_equal(metadata, original)


def test_de_partial_missing_padj_warns_and_remains_missing() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.01, 0.2],
            "padj": [0.02, None],
        }
    )
    original = de_results.copy(deep=True)

    report = validate_de_results(de_results)

    assert not report.has_errors
    assert [issue.code for issue in report.warnings] == [
        IssueCode.MISSING_DE_STATISTIC
    ]
    assert pd.isna(de_results.loc[1, "padj"])
    assert_frame_equal(de_results, original)


def test_de_hints_at_missing_value_placeholders() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, "nd"],
            "pvalue": [0.01, 0.2],
            "padj": [0.02, 0.3],
        }
    )

    report = validate_de_results(de_results)

    non_numeric = next(
        issue for issue in report.issues if issue.code is IssueCode.NON_NUMERIC_VALUE
    )
    assert "leave the cell blank" in non_numeric.message


def test_de_rejects_completely_empty_table_without_mutating_it() -> None:
    de_results = pd.DataFrame()
    original = de_results.copy(deep=True)

    report = validate_de_results(de_results)

    empty_issue = next(
        issue for issue in report.issues if issue.code is IssueCode.EMPTY_TABLE
    )
    assert empty_issue.severity is Severity.ERROR
    assert report.has_errors
    assert_frame_equal(de_results, original)


@pytest.mark.parametrize("column", ["log2FoldChange", "pvalue", "padj"])
def test_de_rejects_required_statistic_with_all_values_missing(column: str) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )
    de_results[column] = None

    report = validate_de_results(de_results)

    issue = next(issue for issue in report.errors if issue.column == column)
    assert issue.code is IssueCode.NO_USABLE_VALUES
    assert report.has_errors


def test_de_distinguishes_existing_missing_from_invalid_coercion() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "log2FoldChange": [1.0, None, "not-a-number"],
            "pvalue": [0.1, 0.2, 0.3],
            "padj": [0.2, 0.3, 0.4],
        }
    )

    report = validate_de_results(de_results)

    logfc_issues = [
        issue for issue in report.issues if issue.column == "log2FoldChange"
    ]
    assert {issue.code for issue in logfc_issues} == {
        IssueCode.NON_NUMERIC_VALUE,
        IssueCode.MISSING_DE_STATISTIC,
    }
    assert next(
        issue for issue in logfc_issues if issue.code is IssueCode.NON_NUMERIC_VALUE
    ).count == 1
    assert next(
        issue for issue in logfc_issues if issue.code is IssueCode.MISSING_DE_STATISTIC
    ).count == 1


@pytest.mark.parametrize("column", ["pvalue", "padj"])
@pytest.mark.parametrize("value", [-0.1, 1.1])
def test_de_rejects_out_of_range_probabilities(column: str, value: float) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )
    de_results.loc[0, column] = value

    report = validate_de_results(de_results)

    assert any(
        issue.code is IssueCode.VALUE_OUT_OF_RANGE and issue.column == column
        for issue in report.errors
    )


def test_de_accepts_zero_padj_without_classifying_significance() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.0],
            "padj": [0.0],
        }
    )

    report = validate_de_results(de_results)

    assert not report.has_errors
    assert IssueCode.ZERO_ADJUSTED_P_VALUE in _codes(report)
    assert "significant" not in de_results.columns


@pytest.mark.parametrize(
    "missing_column", ["gene_id", "log2FoldChange", "pvalue", "padj"]
)
def test_de_rejects_missing_required_columns(missing_column: str) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    ).drop(columns=missing_column)

    report = validate_de_results(de_results)

    assert any(
        issue.code is IssueCode.MISSING_REQUIRED_COLUMN
        and issue.column == missing_column
        for issue in report.errors
    )


def test_de_rejects_missing_and_duplicate_gene_ids() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g1", None],
            "log2FoldChange": [1.0, 2.0, 3.0],
            "pvalue": [0.1, 0.2, 0.3],
            "padj": [0.2, 0.3, 0.4],
        }
    )

    report = validate_de_results(de_results)

    assert IssueCode.MISSING_IDENTIFIER in _codes(report)
    assert IssueCode.DUPLICATE_IDENTIFIER in _codes(report)


def test_de_duplicate_gene_id_is_error_without_mutating_input() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g1"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.1, 0.2],
            "padj": [0.2, 0.3],
        }
    )
    original = de_results.copy(deep=True)

    report = validate_de_results(de_results)

    assert IssueCode.DUPLICATE_IDENTIFIER in [
        issue.code for issue in report.errors
    ]
    assert IssueCode.DUPLICATE_IDENTIFIER not in [
        issue.code for issue in report.warnings
    ]
    assert IssueCode.DUPLICATE_IDENTIFIER not in [
        issue.code for issue in report.information
    ]
    assert report.has_errors is True
    assert_frame_equal(de_results, original)


def test_de_reports_numeric_strings_and_rejects_infinite_values() -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": ["1.0", float("inf")],
            "pvalue": ["0.1", "0.2"],
            "padj": ["0.2", "0.3"],
        }
    )

    report = validate_de_results(de_results)

    assert IssueCode.COERCIBLE_NUMERIC_STRING in _codes(report)
    assert IssueCode.NON_FINITE_VALUE in _codes(report)
    assert report.has_errors


@pytest.mark.parametrize("column", ["log2FoldChange", "pvalue", "padj"])
def test_de_rejects_complex_statistics_without_raising_or_mutating(column: str) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )
    de_results[column] = [1 + 2j]
    original = de_results.copy(deep=True)

    report = validate_de_results(de_results)

    matching = [issue for issue in report.errors if issue.column == column]
    assert [issue.code for issue in matching] == [IssueCode.NON_NUMERIC_VALUE]
    assert report.has_errors is True
    assert_frame_equal(de_results, original)


@pytest.mark.parametrize("column", ["log2FoldChange", "pvalue", "padj"])
@pytest.mark.parametrize("value", [np.bool_(False), np.bool_(True)])
def test_de_rejects_numpy_boolean_statistics_without_coercing_or_mutating(
    column: str, value: np.bool_
) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, 2.0],
            "pvalue": [0.1, 0.2],
            "padj": [0.2, 0.3],
        },
        dtype="object",
    )
    de_results.loc[1, column] = value
    original = de_results.copy(deep=True)

    report = validate_de_results(de_results)

    matching = [issue for issue in report.errors if issue.column == column]
    assert [issue.code for issue in matching] == [IssueCode.NON_NUMERIC_VALUE]
    assert matching[0].row_positions == (2,)
    assert report.has_errors is True
    assert_frame_equal(de_results, original)


@pytest.mark.parametrize("value", [0.0, 1.0])
def test_de_accepts_probability_boundaries(value: float) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [1.0],
            "pvalue": [value],
            "padj": [value],
        }
    )

    report = validate_de_results(de_results)

    assert not report.has_errors


def test_all_single_table_validators_preserve_original_dataframes() -> None:
    expression = pd.DataFrame(
        {"gene_id": [" g1 "], "sample1": ["1.0"], "sample2": [-1.0]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s2", "s1"], "condition": ["a", "b"]}
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, None],
            "pvalue": [0.1, 0.2],
            "padj": [0.2, None],
        }
    )
    originals = [table.copy(deep=True) for table in (expression, metadata, de_results)]

    validate_expression_matrix(expression)
    validate_sample_metadata(metadata)
    validate_de_results(de_results)

    for actual, expected in zip((expression, metadata, de_results), originals):
        assert_frame_equal(actual, expected)
