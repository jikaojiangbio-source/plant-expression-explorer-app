"""Regression tests for Phase 8 differential-expression exploration."""

from __future__ import annotations

import ast
import math
from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest

import plant_expression_explorer as package
from plant_expression_explorer.differential_expression import (
    MA_PLOT_COLUMNS,
    STATUS_COLUMN,
    VOLCANO_PLOT_COLUMNS,
    DifferentialExpressionComputationError,
    DifferentialExpressionErrorReason,
    DifferentialExpressionResult,
    DifferentialExpressionStatus,
    build_category_summary,
    build_ma_plot_data,
    build_volcano_plot_data,
    classify_differential_expression_results,
    select_rows_by_status,
)
from plant_expression_explorer.validation import IssueCode


def _valid_results() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "gene_id": ["positive", "negative", "other", "missing"],
            "log2FoldChange": [1.0, -1.0, 0.5, None],
            "pvalue": [0.01, 0.02, 0.5, 0.8],
            "padj": [0.05, 0.05, 0.5, 0.8],
        },
        index=pd.Index([8, 3, 8, 1], name="source_index"),
    )


def _classify(table: pd.DataFrame, **kwargs: object) -> DifferentialExpressionResult:
    return classify_differential_expression_results(table, **kwargs)


def _blocking_error(table: object) -> DifferentialExpressionComputationError:
    original = table.copy(deep=True) if isinstance(table, pd.DataFrame) else None
    with pytest.raises(DifferentialExpressionComputationError) as captured:
        _classify(table)  # type: ignore[arg-type]
    assert captured.value.reason is DifferentialExpressionErrorReason.BLOCKING_DE_VALIDATION
    assert captured.value.validation_report is not None
    if original is not None:
        assert_frame_equal(table, original, check_exact=True)
    return captured.value


def test_module_has_no_streamlit_or_session_state_dependency() -> None:
    module_path = Path(__file__).parents[1] / "plant_expression_explorer" / "differential_expression.py"
    tree = ast.parse(module_path.read_text(encoding="utf-8"))
    imported_roots = {
        alias.name.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imported_roots.update(
        (node.module or "").split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    )

    assert "streamlit" not in imported_roots
    assert "session_state" not in module_path.read_text(encoding="utf-8")


def test_package_exports_the_approved_public_api() -> None:
    assert package.STATUS_COLUMN == STATUS_COLUMN
    assert package.DifferentialExpressionStatus is DifferentialExpressionStatus
    assert package.DifferentialExpressionErrorReason is DifferentialExpressionErrorReason
    assert package.DifferentialExpressionComputationError is DifferentialExpressionComputationError
    assert package.DifferentialExpressionResult is DifferentialExpressionResult
    assert package.classify_differential_expression_results is classify_differential_expression_results
    assert package.build_category_summary is build_category_summary
    assert package.select_rows_by_status is select_rows_by_status
    assert package.VOLCANO_PLOT_COLUMNS == VOLCANO_PLOT_COLUMNS
    assert package.build_volcano_plot_data is build_volcano_plot_data
    assert package.MA_PLOT_COLUMNS == MA_PLOT_COLUMNS
    assert package.build_ma_plot_data is build_ma_plot_data


def test_all_four_statuses_and_counts_are_mutually_exclusive_and_exhaustive() -> None:
    result = _classify(_valid_results())

    assert result.annotated_results[STATUS_COLUMN].tolist() == [
        DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH.value,
        DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH.value,
        DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS.value,
        DifferentialExpressionStatus.NOT_EVALUABLE.value,
    ]
    assert result.total_row_count == 4
    assert result.evaluable_row_count == 3
    assert result.positive_threshold_match_count == 1
    assert result.negative_threshold_match_count == 1
    assert result.does_not_meet_combined_thresholds_count == 1
    assert result.not_evaluable_count == 1
    assert sum(build_category_summary(result)["row_count"]) == result.total_row_count


def test_build_volcano_plot_data_excludes_not_evaluable_rows() -> None:
    result = _classify(_valid_results())

    volcano = build_volcano_plot_data(result)

    assert list(volcano.plot_rows.columns) == list(VOLCANO_PLOT_COLUMNS)
    assert volcano.plot_rows["gene_id"].tolist() == ["positive", "negative", "other"]
    assert volcano.excluded_zero_padj_count == 0
    positive_row = volcano.plot_rows.set_index("gene_id").loc["positive"]
    assert positive_row["log2FoldChange"] == 1.0
    assert positive_row["padj"] == 0.05
    assert positive_row["neg_log10_padj"] == pytest.approx(-math.log10(0.05))
    assert (
        positive_row[STATUS_COLUMN]
        == DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH.value
    )


def test_build_volcano_plot_data_excludes_and_counts_zero_padj_rows() -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [2.0, -2.0],
            "pvalue": [0.0, 0.01],
            "padj": [0.0, 0.01],
        }
    )
    result = _classify(table)

    volcano = build_volcano_plot_data(result)

    assert volcano.plot_rows["gene_id"].tolist() == ["g2"]
    assert volcano.excluded_zero_padj_count == 1


def test_build_volcano_plot_data_does_not_mutate_the_annotated_results() -> None:
    result = _classify(_valid_results())
    original = result.annotated_results.copy(deep=True)

    build_volcano_plot_data(result)

    assert_frame_equal(result.annotated_results, original, check_exact=True)


def _matching_expression() -> pd.DataFrame:
    # Matches _valid_results()'s evaluable gene_ids: "positive", "negative",
    # "other". "missing" is intentionally absent (its row is NOT_EVALUABLE
    # anyway, per _valid_results()).
    return pd.DataFrame(
        {
            "gene_id": ["positive", "negative", "other"],
            "s1": [10.0, 20.0, 30.0],
            "s2": [12.0, 22.0, 34.0],
        }
    )


def test_build_ma_plot_data_excludes_not_evaluable_rows() -> None:
    result = _classify(_valid_results())

    ma_plot = build_ma_plot_data(result, _matching_expression())

    assert list(ma_plot.plot_rows.columns) == list(MA_PLOT_COLUMNS)
    assert ma_plot.plot_rows["gene_id"].tolist() == ["positive", "negative", "other"]
    assert ma_plot.excluded_no_expression_match_count == 0
    assert ma_plot.excluded_all_missing_expression_count == 0
    positive_row = ma_plot.plot_rows.set_index("gene_id").loc["positive"]
    assert positive_row["mean_expression"] == pytest.approx(11.0)
    assert positive_row["log2FoldChange"] == 1.0
    assert (
        positive_row[STATUS_COLUMN]
        == DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH.value
    )


def test_build_ma_plot_data_excludes_and_counts_genes_absent_from_expression() -> None:
    result = _classify(_valid_results())
    expression = pd.DataFrame(
        {"gene_id": ["positive", "negative"], "s1": [10.0, 20.0], "s2": [12.0, 22.0]}
    )

    ma_plot = build_ma_plot_data(result, expression)

    assert ma_plot.plot_rows["gene_id"].tolist() == ["positive", "negative"]
    assert ma_plot.excluded_no_expression_match_count == 1
    assert ma_plot.excluded_all_missing_expression_count == 0


def test_build_ma_plot_data_excludes_a_duplicated_ambiguous_gene_id() -> None:
    result = _classify(_valid_results())
    expression = pd.DataFrame(
        {
            "gene_id": ["positive", "positive", "negative", "other"],
            "s1": [10.0, 999.0, 20.0, 30.0],
            "s2": [12.0, 999.0, 22.0, 34.0],
        }
    )

    ma_plot = build_ma_plot_data(result, expression)

    assert ma_plot.plot_rows["gene_id"].tolist() == ["negative", "other"]
    assert ma_plot.excluded_no_expression_match_count == 1


def test_build_ma_plot_data_excludes_and_counts_all_missing_expression() -> None:
    result = _classify(_valid_results())
    expression = pd.DataFrame(
        {
            "gene_id": ["positive", "negative", "other"],
            "s1": [None, 20.0, 30.0],
            "s2": [None, 22.0, 34.0],
        }
    )

    ma_plot = build_ma_plot_data(result, expression)

    assert ma_plot.plot_rows["gene_id"].tolist() == ["negative", "other"]
    assert ma_plot.excluded_no_expression_match_count == 0
    assert ma_plot.excluded_all_missing_expression_count == 1


def test_build_ma_plot_data_skips_a_partially_missing_gene_value_only() -> None:
    result = _classify(_valid_results())
    expression = pd.DataFrame(
        {
            "gene_id": ["positive", "negative", "other"],
            "s1": [None, 20.0, 30.0],
            "s2": [12.0, 22.0, 34.0],
        }
    )

    ma_plot = build_ma_plot_data(result, expression)

    positive_row = ma_plot.plot_rows.set_index("gene_id").loc["positive"]
    assert positive_row["mean_expression"] == pytest.approx(12.0)


def test_build_ma_plot_data_does_not_mutate_inputs() -> None:
    result = _classify(_valid_results())
    original_annotated = result.annotated_results.copy(deep=True)
    expression = _matching_expression()
    original_expression = expression.copy(deep=True)

    build_ma_plot_data(result, expression)

    assert_frame_equal(result.annotated_results, original_annotated, check_exact=True)
    assert_frame_equal(expression, original_expression, check_exact=True)


def test_result_is_frozen_and_category_summary_has_stable_order() -> None:
    result = _classify(_valid_results())

    with pytest.raises(FrozenInstanceError):
        result.total_row_count = 99  # type: ignore[misc]
    assert build_category_summary(result).to_dict("records") == [
        {"status": "NOT_EVALUABLE", "row_count": 1},
        {"status": "POSITIVE_THRESHOLD_MATCH", "row_count": 1},
        {"status": "NEGATIVE_THRESHOLD_MATCH", "row_count": 1},
        {"status": "DOES_NOT_MEET_COMBINED_THRESHOLDS", "row_count": 1},
    ]


def test_threshold_comparisons_are_inclusive_without_tolerance() -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["padj_edge", "positive_edge", "negative_edge", "padj_over", "fc_under"],
            "log2FoldChange": [2.0, 1.0, -1.0, 2.0, 0.999999999],
            "pvalue": [0.01] * 5,
            "padj": [0.05, 0.01, 0.01, 0.050000001, 0.01],
        }
    )

    result = _classify(table)

    assert result.annotated_results[STATUS_COLUMN].tolist() == [
        "POSITIVE_THRESHOLD_MATCH",
        "POSITIVE_THRESHOLD_MATCH",
        "NEGATIVE_THRESHOLD_MATCH",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
    ]


def test_zero_adjusted_p_value_and_zero_adjusted_threshold_are_valid() -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.0, 0.1],
            "padj": [0.0, 0.1],
        }
    )

    result = _classify(table, adjusted_p_value_threshold=0)

    assert result.adjusted_p_value_threshold == 0.0
    assert result.annotated_results[STATUS_COLUMN].tolist() == [
        "POSITIVE_THRESHOLD_MATCH",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
    ]
    assert table.loc[0, "padj"] == 0.0


@pytest.mark.parametrize(
    ("log2_fold_change", "padj"),
    [(1.0, None), (None, 0.01), (None, None), ("", 0.01), (1.0, "   ")],
)
def test_missing_or_blank_classification_values_are_not_evaluable(
    log2_fold_change: object, padj: object
) -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [2.0, log2_fold_change],
            "pvalue": [0.01, 0.02],
            "padj": [0.01, padj],
        }
    )

    result = _classify(table)

    assert result.annotated_results.iloc[1][STATUS_COLUMN] == "NOT_EVALUABLE"
    assert result.not_evaluable_count == 1


def test_missing_pvalue_alone_does_not_affect_classification() -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [2.0, 0.0],
            "pvalue": [0.01, None],
            "padj": [0.01, 0.01],
        }
    )

    result = _classify(table)

    assert result.annotated_results[STATUS_COLUMN].tolist() == [
        "POSITIVE_THRESHOLD_MATCH",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
    ]
    assert result.not_evaluable_count == 0


def test_numeric_strings_are_temporary_and_original_values_and_dtypes_are_preserved() -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": ["1.0", "-1.0"],
            "pvalue": ["0.01", "0.02"],
            "padj": ["0.05", "0.05"],
        }
    )
    original = table.copy(deep=True)

    result = _classify(table)

    assert result.annotated_results[STATUS_COLUMN].tolist() == [
        "POSITIVE_THRESHOLD_MATCH",
        "NEGATIVE_THRESHOLD_MATCH",
    ]
    assert_frame_equal(table, original, check_exact=True)
    assert_frame_equal(
        result.annotated_results.iloc[:, :-1], original, check_exact=True
    )


@pytest.mark.parametrize("column", ["log2FoldChange", "pvalue", "padj"])
@pytest.mark.parametrize(
    "value", ["not-numeric", float("inf"), float("-inf"), True, 1 + 2j]
)
def test_invalid_statistic_values_are_controlled_and_non_mutating(
    column: str, value: object
) -> None:
    table = pd.DataFrame(
        {"gene_id": ["g1"], "log2FoldChange": [1.0], "pvalue": [0.01], "padj": [0.02]}
    )
    table[column] = [value]
    original = table.copy(deep=True)

    error = _blocking_error(table)

    expected_code = (
        IssueCode.NON_FINITE_VALUE
        if isinstance(value, float) and not pd.notna(value - value)
        else IssueCode.NON_NUMERIC_VALUE
    )
    assert expected_code in {issue.code for issue in error.validation_report.errors}
    assert_frame_equal(table, original, check_exact=True)


@pytest.mark.parametrize("column", ["log2FoldChange", "pvalue", "padj"])
@pytest.mark.parametrize("value", [np.bool_(False), np.bool_(True)])
def test_numpy_boolean_statistics_block_classification_without_mutation(
    column: str, value: np.bool_
) -> None:
    table = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [1.0, 2.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.02, 0.03],
        },
        dtype="object",
    )
    table.loc[1, column] = value

    error = _blocking_error(table)

    matching = [
        issue
        for issue in error.validation_report.errors
        if issue.column == column
    ]
    assert [issue.code for issue in matching] == [IssueCode.NON_NUMERIC_VALUE]


@pytest.mark.parametrize(
    "value", [None, "0.05", True, 0.05 + 0j, [0.05], float("nan"), float("inf"), -0.01, 1.01]
)
def test_invalid_adjusted_p_value_thresholds_are_controlled(value: object) -> None:
    original = _valid_results()
    table = original.copy(deep=True)

    with pytest.raises(DifferentialExpressionComputationError) as captured:
        _classify(table, adjusted_p_value_threshold=value)

    assert captured.value.reason is DifferentialExpressionErrorReason.INVALID_ADJUSTED_P_VALUE_THRESHOLD
    assert_frame_equal(table, original, check_exact=True)


@pytest.mark.parametrize(
    "value", [None, "1", True, 1 + 0j, [1], float("nan"), float("inf"), float("-inf"), 0, -0.1]
)
def test_invalid_fold_change_thresholds_are_controlled(value: object) -> None:
    original = _valid_results()
    table = original.copy(deep=True)

    with pytest.raises(DifferentialExpressionComputationError) as captured:
        _classify(table, absolute_log2_fold_change_threshold=value)

    assert captured.value.reason is DifferentialExpressionErrorReason.INVALID_ABSOLUTE_LOG2_FOLD_CHANGE_THRESHOLD
    assert_frame_equal(table, original, check_exact=True)


@pytest.mark.parametrize("missing", ["gene_id", "log2FoldChange", "pvalue", "padj"])
def test_every_missing_required_column_is_blocking(missing: str) -> None:
    _blocking_error(_valid_results().drop(columns=missing))


def test_empty_results_are_blocking() -> None:
    _blocking_error(pd.DataFrame(columns=["gene_id", "log2FoldChange", "pvalue", "padj"]))


def test_duplicate_required_columns_are_blocking() -> None:
    table = pd.DataFrame(
        [["g1", 1.0, 1.1, 0.01, 0.02]],
        columns=["gene_id", "log2FoldChange", "log2FoldChange", "pvalue", "padj"],
    )

    error = _blocking_error(table)

    assert IssueCode.DUPLICATE_COLUMN_NAME in {issue.code for issue in error.validation_report.errors}


@pytest.mark.parametrize("gene_ids", [[None, "g2"], ["", "g2"], [" g1 ", "g2"], ["g1", "g1"]])
def test_invalid_gene_identifiers_are_blocking_without_normalization(gene_ids: list[object]) -> None:
    table = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "log2FoldChange": [1.0, -1.0],
            "pvalue": [0.01, 0.02],
            "padj": [0.02, 0.03],
        }
    )
    original = table.copy(deep=True)

    _blocking_error(table)

    assert_frame_equal(table, original, check_exact=True)


@pytest.mark.parametrize("column", ["log2FoldChange", "pvalue", "padj"])
def test_entirely_missing_statistic_columns_are_blocking(column: str) -> None:
    table = _valid_results()
    table[column] = None

    error = _blocking_error(table)

    assert IssueCode.NO_USABLE_VALUES in {issue.code for issue in error.validation_report.errors}


@pytest.mark.parametrize("column", ["pvalue", "padj"])
@pytest.mark.parametrize("value", [-0.1, 1.1])
def test_out_of_range_probabilities_are_blocking(column: str, value: float) -> None:
    table = _valid_results()
    table.loc[table.index[0], column] = value

    error = _blocking_error(table)

    assert IssueCode.VALUE_OUT_OF_RANGE in {issue.code for issue in error.validation_report.errors}


def test_extra_columns_index_order_values_and_dtypes_are_preserved() -> None:
    table = _valid_results()
    table.insert(1, "annotation", pd.Series(["a", None, "c", "d"], index=table.index, dtype="object"))
    original = table.copy(deep=True)

    result = _classify(table)

    assert result.annotated_results.columns.tolist() == [*original.columns, STATUS_COLUMN]
    assert result.annotated_results.index.equals(original.index)
    assert_frame_equal(result.annotated_results.iloc[:, :-1], original, check_exact=True)
    assert_frame_equal(table, original, check_exact=True)


def test_repeated_computation_is_deterministic_and_does_not_mutate_input() -> None:
    table = _valid_results()
    original = table.copy(deep=True)

    first = _classify(table)
    second = _classify(table)

    assert first.adjusted_p_value_threshold == second.adjusted_p_value_threshold
    assert first.absolute_log2_fold_change_threshold == second.absolute_log2_fold_change_threshold
    assert first.total_row_count == second.total_row_count
    assert first.evaluable_row_count == second.evaluable_row_count
    assert first.not_evaluable_count == second.not_evaluable_count
    assert first.positive_threshold_match_count == second.positive_threshold_match_count
    assert first.negative_threshold_match_count == second.negative_threshold_match_count
    assert (
        first.does_not_meet_combined_thresholds_count
        == second.does_not_meet_combined_thresholds_count
    )
    assert_frame_equal(first.annotated_results, second.annotated_results, check_exact=True)
    assert_frame_equal(table, original, check_exact=True)


def test_status_subsets_preserve_order_index_and_are_independent_copies() -> None:
    result = _classify(_valid_results())

    positive = select_rows_by_status(result, DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH)
    negative = select_rows_by_status(result, DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH)
    other = select_rows_by_status(result, DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS)
    missing = select_rows_by_status(result, DifferentialExpressionStatus.NOT_EVALUABLE)

    assert positive["gene_id"].tolist() == ["positive"]
    assert negative["gene_id"].tolist() == ["negative"]
    assert other["gene_id"].tolist() == ["other"]
    assert missing["gene_id"].tolist() == ["missing"]
    assert [positive.index[0], negative.index[0], other.index[0], missing.index[0]] == [8, 3, 8, 1]
    positive.iloc[0, 0] = "changed"
    assert result.annotated_results.iloc[0]["gene_id"] == "positive"


def test_reserved_status_column_collision_is_controlled_and_non_mutating() -> None:
    table = _valid_results()
    table[STATUS_COLUMN] = "supplied"
    original = table.copy(deep=True)

    with pytest.raises(DifferentialExpressionComputationError) as captured:
        _classify(table)

    assert captured.value.reason is DifferentialExpressionErrorReason.STATUS_COLUMN_CONFLICT
    assert_frame_equal(table, original, check_exact=True)


def test_multiple_validation_errors_are_retained_in_controlled_report() -> None:
    table = pd.DataFrame(
        {
            "gene_id": [None, None],
            "log2FoldChange": ["bad", float("inf")],
            "pvalue": [-1.0, 2.0],
            "padj": [None, None],
        }
    )

    error = _blocking_error(table)
    codes = {issue.code for issue in error.validation_report.errors}

    assert {
        IssueCode.MISSING_IDENTIFIER,
        IssueCode.NON_NUMERIC_VALUE,
        IssueCode.NON_FINITE_VALUE,
        IssueCode.VALUE_OUT_OF_RANGE,
        IssueCode.NO_USABLE_VALUES,
    } <= codes
