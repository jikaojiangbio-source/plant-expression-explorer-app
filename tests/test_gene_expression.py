"""Regression tests for Phase 9 descriptive gene-expression lookup."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest

import plant_expression_explorer as package
from plant_expression_explorer.dataset import load_demo_candidate
from plant_expression_explorer.gene_expression import (
    CONDITION_EXPRESSION_SUMMARY_COLUMNS,
    GENE_EXPRESSION_CHART_COLUMNS,
    SAMPLE_EXPRESSION_COLUMNS,
    GeneExpressionComputationError,
    GeneExpressionErrorReason,
    GeneExpressionResult,
    build_gene_expression_chart_data,
    build_gene_expression_observations,
    list_gene_ids,
    lookup_gene_expression,
)
from plant_expression_explorer.validation import ValidationReport


_METADATA_NOT_SUPPLIED = object()


def _known_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": [101, 202],
            "sample_b": ["2.0", "8.0"],
            "sample_a": ["1.0", "4.0"],
            "sample_c": ["3.0", "6.0"],
        },
        index=pd.Index([7, 3], name="source_index"),
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_c", "sample_b"],
            "condition": ["Control", "Treated", "Control"],
            "unused": [1, 2, 3],
        },
        index=[20, 30, 10],
    )
    return expression, metadata


def _capture_error(
    expression: object,
    metadata: object = _METADATA_NOT_SUPPLIED,
    gene_id: object = "g1",
) -> GeneExpressionComputationError:
    with pytest.raises(GeneExpressionComputationError) as captured:
        if metadata is _METADATA_NOT_SUPPLIED:
            list_gene_ids(expression)  # type: ignore[arg-type]
        else:
            lookup_gene_expression(
                expression,  # type: ignore[arg-type]
                metadata,  # type: ignore[arg-type]
                gene_id,
            )
    return captured.value


def test_module_has_no_streamlit_or_session_state_dependency() -> None:
    module_path = (
        Path(__file__).parents[1]
        / "plant_expression_explorer"
        / "gene_expression.py"
    )
    source = module_path.read_text(encoding="utf-8")
    tree = ast.parse(source)
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
    assert "session_state" not in source


def test_package_exports_the_approved_public_api() -> None:
    assert package.SAMPLE_EXPRESSION_COLUMNS == SAMPLE_EXPRESSION_COLUMNS
    assert (
        package.CONDITION_EXPRESSION_SUMMARY_COLUMNS
        == CONDITION_EXPRESSION_SUMMARY_COLUMNS
    )
    assert package.GENE_EXPRESSION_CHART_COLUMNS == GENE_EXPRESSION_CHART_COLUMNS
    assert package.GeneExpressionErrorReason is GeneExpressionErrorReason
    assert package.GeneExpressionComputationError is GeneExpressionComputationError
    assert package.GeneExpressionResult is GeneExpressionResult
    assert package.list_gene_ids is list_gene_ids
    assert package.lookup_gene_expression is lookup_gene_expression
    assert (
        package.build_gene_expression_chart_data
        is build_gene_expression_chart_data
    )
    assert (
        package.build_gene_expression_observations
        is build_gene_expression_observations
    )


def test_error_reason_values_are_unique() -> None:
    values = [reason.value for reason in GeneExpressionErrorReason]
    assert len(values) == len(set(values))


def test_demo_options_and_lookup_preserve_exact_source_order_and_values() -> None:
    candidate = load_demo_candidate()
    assert candidate.tables is not None
    expression, metadata, _de_results = candidate.tables
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    gene_ids = list_gene_ids(expression)
    result = lookup_gene_expression(expression, metadata, gene_ids[0])

    assert len(gene_ids) == 120
    assert gene_ids[:2] == ("SYN_Solyc_0001", "SYN_Solyc_0002")
    assert gene_ids[-1] == "SYN_Solyc_0120"
    assert result.gene_id == "SYN_Solyc_0001"
    assert result.sample_expression["sample_id"].tolist() == list(
        expression.columns[1:]
    )
    assert result.sample_expression["condition"].tolist() == metadata[
        "condition"
    ].tolist()
    assert result.sample_expression["expression_value"].tolist() == expression.iloc[
        0, 1:
    ].tolist()
    assert result.condition_summary["condition"].tolist() == [
        "Control",
        "High_nitrate",
    ]
    assert_frame_equal(expression, original_expression, check_exact=True)
    assert_frame_equal(metadata, original_metadata, check_exact=True)


def test_identifier_matching_uses_exact_string_representations_without_source_changes() -> None:
    expression, metadata = _known_tables()
    original_expression = expression.copy(deep=True)

    assert list_gene_ids(expression) == ("101", "202")
    from_numeric_query = lookup_gene_expression(expression, metadata, 101)
    from_string_query = lookup_gene_expression(expression, metadata, "101")

    assert from_numeric_query.gene_id == from_string_query.gene_id == "101"
    assert from_numeric_query.sample_expression["expression_value"].tolist() == [
        "2.0",
        "1.0",
        "3.0",
    ]
    assert all(
        isinstance(value, str)
        for value in from_numeric_query.sample_expression["expression_value"]
    )
    assert expression["gene_id"].tolist() == [101, 202]
    assert_frame_equal(expression, original_expression, check_exact=True)


def test_sample_order_and_condition_mapping_follow_expression_columns() -> None:
    expression, metadata = _known_tables()

    result = lookup_gene_expression(expression, metadata, "202")

    assert result.metadata_order_matches_expression is False
    assert result.sample_expression.to_dict("records") == [
        {"sample_id": "sample_b", "condition": "Control", "expression_value": "8.0"},
        {"sample_id": "sample_a", "condition": "Control", "expression_value": "4.0"},
        {"sample_id": "sample_c", "condition": "Treated", "expression_value": "6.0"},
    ]
    assert result.condition_summary["condition"].tolist() == ["Control", "Treated"]
    assert result.condition_summary["sample_ids"].tolist() == [
        ("sample_b", "sample_a"),
        ("sample_c",),
    ]


def test_condition_summary_uses_disclosed_arithmetic_conventions() -> None:
    expression, metadata = _known_tables()

    result = lookup_gene_expression(expression, metadata, "202")
    control = result.condition_summary.iloc[0]
    treated = result.condition_summary.iloc[1]

    assert result.condition_summary.columns.tolist() == list(
        CONDITION_EXPRESSION_SUMMARY_COLUMNS
    )
    assert control["sample_count"] == 2
    assert control["minimum_expression"] == 4.0
    assert control["median_expression"] == 6.0
    assert control["mean_expression"] == 6.0
    assert control["maximum_expression"] == 8.0
    assert control["standard_deviation"] == pytest.approx(np.sqrt(8.0))
    assert treated["sample_count"] == 1
    assert treated["mean_expression"] == 6.0
    assert pd.isna(treated["standard_deviation"])


@pytest.mark.parametrize(
    ("values", "expected_median", "expected_mean", "expected_sd"),
    [
        ([1e308, 1e308], 1e308, 1e308, 0.0),
        ([1e308, -1e308], 0.0, 0.0, np.sqrt(2.0) * 1e308),
    ],
)
def test_condition_summaries_are_stable_for_large_finite_values(
    values: list[float],
    expected_median: float,
    expected_mean: float,
    expected_sd: float,
) -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "s1": [values[0]], "s2": [values[1]]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["A", "A"]}
    )
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    result = lookup_gene_expression(expression, metadata, "g1")
    summary = result.condition_summary.iloc[0]

    assert summary["median_expression"] == expected_median
    assert summary["mean_expression"] == expected_mean
    assert summary["standard_deviation"] == pytest.approx(expected_sd)
    assert np.isfinite(
        summary[
            [
                "minimum_expression",
                "median_expression",
                "mean_expression",
                "maximum_expression",
                "standard_deviation",
            ]
        ].astype(float)
    ).all()
    assert_frame_equal(expression, original_expression, check_exact=True)
    assert_frame_equal(metadata, original_metadata, check_exact=True)


def test_unrepresentable_derived_standard_deviation_is_controlled_and_non_mutating(
) -> None:
    maximum = np.finfo(float).max
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "s1": [maximum], "s2": [-maximum]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["A", "A"]}
    )
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    assert list_gene_ids(expression) == ("g1",)
    error = _capture_error(expression, metadata)

    assert error.reason is GeneExpressionErrorReason.NUMERICAL_RANGE_ERROR
    assert "no partial summary was returned" in str(error)
    assert_frame_equal(expression, original_expression, check_exact=True)
    assert_frame_equal(metadata, original_metadata, check_exact=True)


def test_chart_data_has_explicit_positions_and_is_an_independent_numeric_copy() -> None:
    expression, metadata = _known_tables()
    result = lookup_gene_expression(expression, metadata, "101")

    chart_data = build_gene_expression_chart_data(result)

    assert chart_data.columns.tolist() == list(GENE_EXPRESSION_CHART_COLUMNS)
    assert chart_data["sample_position"].tolist() == [0, 1, 2]
    assert chart_data["sample_id"].tolist() == ["sample_b", "sample_a", "sample_c"]
    assert chart_data["expression_value"].tolist() == [2.0, 1.0, 3.0]
    assert pd.api.types.is_numeric_dtype(chart_data["expression_value"])
    chart_data.loc[0, "expression_value"] = 999.0
    assert result.sample_expression.loc[0, "expression_value"] == "2.0"


def test_one_sample_one_condition_and_equal_value_boundaries_are_explicit() -> None:
    expression = pd.DataFrame({"gene_id": ["g1", "g2"], "only": [5.0, 5.0]})
    metadata = pd.DataFrame({"sample_id": ["only"], "condition": ["Control"]})

    result = lookup_gene_expression(expression, metadata, "g1")
    observations = build_gene_expression_observations(result, ValidationReport())

    assert list_gene_ids(expression) == ("g1", "g2")
    assert result.sample_count == 1
    assert result.condition_count == 1
    assert result.all_values_equal is True
    assert pd.isna(result.condition_summary.loc[0, "standard_deviation"])
    assert any("Only one sample" in message for message in observations)
    assert any("Only one supplied condition" in message for message in observations)
    assert any("exactly identical" in message for message in observations)
    assert any("replication is not inferred" in message for message in observations)


def test_multiple_equal_values_are_retained_and_reported() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "s1": [-2.0], "s2": [-2.0], "s3": [-2.0]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2", "s3"], "condition": ["A", "A", "A"]}
    )

    result = lookup_gene_expression(expression, metadata, "g1")

    assert result.all_values_equal is True
    assert result.sample_expression["expression_value"].tolist() == [-2.0] * 3
    assert result.condition_summary.loc[0, "standard_deviation"] == 0.0


@pytest.mark.parametrize("query", [None, pd.NA, "", "   ", " g1 "])
def test_invalid_selected_gene_id_is_controlled_without_trimming(query: object) -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    metadata = pd.DataFrame({"sample_id": ["s1"], "condition": ["A"]})

    error = _capture_error(expression, metadata, query)

    assert error.reason is GeneExpressionErrorReason.INVALID_GENE_ID


def test_exact_case_difference_is_unknown_without_fuzzy_fallback() -> None:
    expression = pd.DataFrame({"gene_id": ["GeneA"], "s1": [1.0]})
    metadata = pd.DataFrame({"sample_id": ["s1"], "condition": ["A"]})

    error = _capture_error(expression, metadata, "genea")

    assert error.reason is GeneExpressionErrorReason.UNKNOWN_GENE_ID
    assert "genea" in str(error)


@pytest.mark.parametrize(
    ("table", "expected_reason"),
    [
        (None, GeneExpressionErrorReason.INVALID_TABLE),
        (
            pd.DataFrame(columns=["gene_id", "s1"]),
            GeneExpressionErrorReason.EMPTY_EXPRESSION_TABLE,
        ),
        (
            pd.DataFrame({"gene_id": ["g1"]}),
            GeneExpressionErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
        ),
        (
            pd.DataFrame({"wrong": ["g1"], "s1": [1.0]}),
            GeneExpressionErrorReason.MISSING_REQUIRED_COLUMN,
        ),
    ],
)
def test_list_gene_ids_defensively_validates_structure(
    table: object,
    expected_reason: GeneExpressionErrorReason,
) -> None:
    error = _capture_error(table)
    assert error.reason is expected_reason


def test_duplicate_gene_id_column_is_controlled_by_both_public_lookups() -> None:
    expression = pd.DataFrame(
        [["g1", "g1", 1.0]],
        columns=["gene_id", "gene_id", "s1"],
    )
    metadata = pd.DataFrame({"sample_id": ["s1"], "condition": ["A"]})

    assert _capture_error(expression).reason is GeneExpressionErrorReason.DUPLICATE_REQUIRED_COLUMN
    assert (
        _capture_error(expression, metadata).reason
        is GeneExpressionErrorReason.DUPLICATE_REQUIRED_COLUMN
    )


@pytest.mark.parametrize(
    ("gene_ids", "reason"),
    [
        (["g1", None], GeneExpressionErrorReason.MISSING_REQUIRED_IDENTIFIER),
        (["g1", " g2 "], GeneExpressionErrorReason.WHITESPACE_IN_IDENTIFIER),
        (["g1", "g1"], GeneExpressionErrorReason.DUPLICATE_REQUIRED_IDENTIFIER),
        ([1, "1"], GeneExpressionErrorReason.DUPLICATE_REQUIRED_IDENTIFIER),
    ],
)
def test_invalid_gene_identifiers_block_option_building_and_lookup(
    gene_ids: list[object],
    reason: GeneExpressionErrorReason,
) -> None:
    expression = pd.DataFrame({"gene_id": gene_ids, "s1": [1.0, 2.0]})
    metadata = pd.DataFrame({"sample_id": ["s1"], "condition": ["A"]})
    original = expression.copy(deep=True)

    assert _capture_error(expression).reason is reason
    assert _capture_error(expression, metadata).reason is reason
    assert_frame_equal(expression, original, check_exact=True)


@pytest.mark.parametrize(
    ("value", "reason"),
    [
        (None, GeneExpressionErrorReason.MISSING_EXPRESSION_VALUE),
        ("", GeneExpressionErrorReason.MISSING_EXPRESSION_VALUE),
        ("not-numeric", GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE),
        (True, GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE),
        (False, GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE),
        (np.bool_(True), GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE),
        (np.bool_(False), GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE),
        (1 + 2j, GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE),
        (float("inf"), GeneExpressionErrorReason.INFINITE_EXPRESSION_VALUE),
        (float("-inf"), GeneExpressionErrorReason.INFINITE_EXPRESSION_VALUE),
    ],
)
def test_invalid_expression_values_block_both_public_lookups_without_mutation(
    value: object,
    reason: GeneExpressionErrorReason,
) -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "s1": [1.0, value]},
        dtype="object",
    )
    metadata = pd.DataFrame({"sample_id": ["s1"], "condition": ["A"]})
    original = expression.copy(deep=True)

    assert _capture_error(expression).reason is reason
    assert _capture_error(expression, metadata).reason is reason
    assert_frame_equal(expression, original, check_exact=True)


@pytest.mark.parametrize(
    ("metadata", "reason"),
    [
        (None, GeneExpressionErrorReason.INVALID_TABLE),
        (
            pd.DataFrame({"condition": ["A"]}),
            GeneExpressionErrorReason.MISSING_REQUIRED_COLUMN,
        ),
        (
            pd.DataFrame({"sample_id": ["s1"]}),
            GeneExpressionErrorReason.MISSING_REQUIRED_COLUMN,
        ),
        (
            pd.DataFrame({"sample_id": [None], "condition": ["A"]}),
            GeneExpressionErrorReason.MISSING_REQUIRED_IDENTIFIER,
        ),
        (
            pd.DataFrame({"sample_id": ["s1"], "condition": [None]}),
            GeneExpressionErrorReason.MISSING_REQUIRED_VALUE,
        ),
        (
            pd.DataFrame({"sample_id": ["s1"], "condition": [" A "]}),
            GeneExpressionErrorReason.WHITESPACE_IN_REQUIRED_VALUE,
        ),
    ],
)
def test_lookup_defensively_validates_metadata(
    metadata: object,
    reason: GeneExpressionErrorReason,
) -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    error = _capture_error(expression, metadata)
    assert error.reason is reason


def test_duplicate_metadata_and_sample_column_identifiers_are_controlled() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0], "s2": [2.0]})
    duplicate_metadata = pd.DataFrame(
        {"sample_id": ["s1", "s1"], "condition": ["A", "B"]}
    )
    duplicate_samples = expression.copy(deep=True)
    duplicate_samples.columns = ["gene_id", "s1", "s1"]
    valid_metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["A", "B"]}
    )

    assert (
        _capture_error(expression, duplicate_metadata).reason
        is GeneExpressionErrorReason.DUPLICATE_REQUIRED_IDENTIFIER
    )
    assert (
        _capture_error(duplicate_samples).reason
        is GeneExpressionErrorReason.DUPLICATE_REQUIRED_IDENTIFIER
    )
    assert (
        _capture_error(duplicate_samples, valid_metadata).reason
        is GeneExpressionErrorReason.DUPLICATE_REQUIRED_IDENTIFIER
    )


def test_sample_mismatch_reports_both_directions_and_takes_no_intersection() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "s1": [1.0], "s2": [2.0]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s2", "s3"], "condition": ["A", "B"]}
    )
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    error = _capture_error(expression, metadata)

    assert error.reason is GeneExpressionErrorReason.SAMPLE_MISMATCH
    assert "missing metadata IDs: s1" in str(error)
    assert "metadata-only IDs: s3" in str(error)
    assert "No intersection was taken" in str(error)
    assert_frame_equal(expression, original_expression, check_exact=True)
    assert_frame_equal(metadata, original_metadata, check_exact=True)


def test_result_is_frozen_and_repeated_calls_are_equivalent_and_non_mutating() -> None:
    expression, metadata = _known_tables()
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    first = lookup_gene_expression(expression, metadata, "101")
    second = lookup_gene_expression(expression, metadata, "101")

    with pytest.raises(FrozenInstanceError):
        first.sample_count = 99  # type: ignore[misc]
    assert first is not second
    assert first.sample_expression is not second.sample_expression
    assert first.condition_summary is not second.condition_summary
    assert_frame_equal(first.sample_expression, second.sample_expression, check_exact=True)
    assert_frame_equal(first.condition_summary, second.condition_summary, check_exact=True)
    assert_frame_equal(expression, original_expression, check_exact=True)
    assert_frame_equal(metadata, original_metadata, check_exact=True)
