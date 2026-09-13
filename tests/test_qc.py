"""Tests for Phase 5 descriptive sample quality-control contracts."""

from __future__ import annotations

import math

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from plant_expression_explorer.dataset import load_demo_candidate
from plant_expression_explorer.qc import (
    CONDITION_SUMMARY_COLUMNS,
    SAMPLE_SUMMARY_COLUMNS,
    QcComputationError,
    QcErrorReason,
    build_condition_chart_data,
    build_grouped_condition_chart_data,
    build_grouped_condition_summary,
    build_grouped_sample_summary,
    build_qc_observations,
    build_sample_chart_data,
    compute_sample_qc,
    count_zero_variance_genes,
    find_constant_samples,
    list_additional_metadata_columns,
    should_display_count_chart,
)
from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
)


def _known_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "s2": [0.0, 1.0, 2.0, 3.0],
            "s1": [-2.0, -1.0, 0.0, 1.0],
            "s3": [5.0, 5.0, 5.0, 5.0],
        },
        index=[40, 10, 30, 20],
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3"],
            "condition": ["control", "treated", "treated"],
        },
        index=[8, 5, 3],
    )
    return expression, metadata


def _simple_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "s1": [1.0, 2.0],
            "s2": [3.0, 4.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "condition": ["control", "treated"],
        }
    )
    return expression, metadata


def _assert_inputs_unchanged(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    original_expression: pd.DataFrame,
    original_metadata: pd.DataFrame,
) -> None:
    assert_frame_equal(expression, original_expression)
    assert_frame_equal(metadata, original_metadata)


def test_demo_qc_has_exact_gene_and_sample_counts() -> None:
    candidate = load_demo_candidate()
    assert candidate.tables is not None
    expression, metadata, _ = candidate.tables

    result = compute_sample_qc(expression, metadata)

    assert result.gene_count == 120
    assert result.sample_count == 6
    assert result.expression_cell_count == 720


def test_sample_order_and_metadata_condition_mapping_are_exact() -> None:
    expression, metadata = _known_tables()

    result = compute_sample_qc(expression, metadata)

    assert result.sample_summary["sample_id"].tolist() == ["s2", "s1", "s3"]
    assert result.sample_summary["condition"].tolist() == [
        "treated",
        "control",
        "treated",
    ]
    assert result.metadata_order_matches_expression is False


def test_sample_summary_has_exact_columns() -> None:
    expression, metadata = _known_tables()

    result = compute_sample_qc(expression, metadata)

    assert tuple(result.sample_summary.columns) == SAMPLE_SUMMARY_COLUMNS


def test_known_sample_statistics_use_documented_pandas_conventions() -> None:
    expression, metadata = _known_tables()

    result = compute_sample_qc(expression, metadata)
    sample = result.sample_summary.set_index("sample_id").loc["s2"]

    assert sample["minimum"] == pytest.approx(0.0)
    assert sample["first_quartile"] == pytest.approx(0.75)
    assert sample["median"] == pytest.approx(1.5)
    assert sample["mean"] == pytest.approx(1.5)
    assert sample["third_quartile"] == pytest.approx(2.25)
    assert sample["maximum"] == pytest.approx(3.0)
    assert sample["standard_deviation"] == pytest.approx(math.sqrt(5 / 3))
    assert sample["interquartile_range"] == pytest.approx(1.5)


def test_zero_negative_missing_and_non_finite_counts_are_exact() -> None:
    expression, metadata = _known_tables()

    result = compute_sample_qc(expression, metadata)
    samples = result.sample_summary.set_index("sample_id")

    assert result.zero_value_count == 2
    assert samples.loc["s2", "zero_value_count"] == 1
    assert samples.loc["s1", "zero_value_count"] == 1
    assert result.negative_value_count == 2
    assert samples.loc["s1", "negative_value_count"] == 2
    assert result.missing_value_count == 0
    assert result.non_finite_value_count == 0
    assert samples["missing_value_count"].tolist() == [0, 0, 0]
    assert samples["non_finite_value_count"].tolist() == [0, 0, 0]


def test_numeric_strings_are_calculated_on_a_temporary_copy() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "s1": ["1.0", "2.0", "3.0"],
            "s2": ["4", "5", "6"],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "condition": ["a", "b"],
        }
    )
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    result = compute_sample_qc(expression, metadata)

    assert result.sample_summary["mean"].tolist() == pytest.approx([2.0, 5.0])
    _assert_inputs_unchanged(
        expression,
        metadata,
        original_expression,
        original_metadata,
    )
    assert expression["s1"].dtype == object
    assert expression["s2"].dtype == object


def test_successful_calculation_preserves_values_dtypes_indices_and_order() -> None:
    expression, metadata = _known_tables()
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    compute_sample_qc(expression, metadata)

    _assert_inputs_unchanged(
        expression,
        metadata,
        original_expression,
        original_metadata,
    )


def test_constant_samples_and_zero_variance_genes_use_exact_equality() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "s1": [2.0, 2.0, 2.0],
            "s2": [2.0, 3.0, 4.0],
            "s3": [2.0, 5.0, 6.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s2", "s3"],
            "condition": ["a", "a", "b"],
        }
    )

    result = compute_sample_qc(expression, metadata)

    assert result.constant_samples == ("s1",)
    assert result.zero_variance_gene_count == 1


def test_constant_sample_and_zero_variance_gene_axes_are_distinguished() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "s_constant": [5.0, 5.0, 5.0],
            "s_variable_1": [5.0, 5.0, 6.0],
            "s_variable_2": [5.0, 5.0, 7.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s_constant", "s_variable_1", "s_variable_2"],
            "condition": ["control", "control", "treated"],
        }
    )
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)
    numeric_expression = expression.drop(columns="gene_id")

    assert find_constant_samples(numeric_expression) == ("s_constant",)
    assert count_zero_variance_genes(numeric_expression) == 2

    result = compute_sample_qc(expression, metadata)

    assert result.constant_samples == ("s_constant",)
    assert result.zero_variance_gene_count == 2
    _assert_inputs_unchanged(
        expression,
        metadata,
        original_expression,
        original_metadata,
    )


def test_condition_summary_uses_first_appearance_and_metadata_row_order() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "s3": [1.0, 2.0],
            "s1": [2.0, 3.0],
            "s2": [3.0, 4.0],
            "s4": [4.0, 5.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s3", "s2", "s4"],
            "condition": ["control", "treated", "control", "single"],
        }
    )

    result = compute_sample_qc(expression, metadata)

    assert tuple(result.condition_summary.columns) == CONDITION_SUMMARY_COLUMNS
    assert result.condition_summary["condition"].tolist() == [
        "control",
        "treated",
        "single",
    ]
    assert result.condition_summary["sample_count"].tolist() == [2, 1, 1]
    assert result.condition_summary["sample_ids"].tolist() == [
        ("s1", "s2"),
        ("s3",),
        ("s4",),
    ]


def test_list_additional_metadata_columns_excludes_required_columns() -> None:
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1"],
            "condition": ["control"],
            "genotype": ["WT"],
            "batch": ["1"],
        }
    )

    assert list_additional_metadata_columns(metadata) == ("genotype", "batch")


def test_list_additional_metadata_columns_returns_empty_for_non_dataframe() -> None:
    assert list_additional_metadata_columns(None) == ()


def test_build_grouped_sample_summary_matches_default_for_condition() -> None:
    expression, metadata = _known_tables()
    result = compute_sample_qc(expression, metadata)

    assert_frame_equal(
        build_grouped_sample_summary(result, metadata, "condition"),
        result.sample_summary,
    )


def test_build_grouped_sample_summary_relabels_without_recomputing_statistics() -> None:
    expression, metadata = _known_tables()
    metadata = metadata.assign(genotype=["mutant", "WT", "WT"])
    result = compute_sample_qc(expression, metadata)

    grouped = build_grouped_sample_summary(result, metadata, "genotype")

    assert "genotype" in grouped.columns
    assert "condition" not in grouped.columns
    expected_genotype = {"s1": "mutant", "s2": "WT", "s3": "WT"}
    for sample_id, genotype in zip(
        grouped["sample_id"], grouped["genotype"], strict=True
    ):
        assert genotype == expected_genotype[sample_id]
    # Numeric statistics are untouched by the relabelling.
    for column in ("minimum", "median", "mean", "maximum"):
        pd.testing.assert_series_equal(
            grouped[column], result.sample_summary[column], check_names=False
        )


def test_build_grouped_condition_summary_matches_default_for_condition() -> None:
    expression, metadata = _known_tables()
    result = compute_sample_qc(expression, metadata)

    assert_frame_equal(
        build_grouped_condition_summary(result, metadata, "condition"),
        result.condition_summary,
    )


def test_build_grouped_condition_summary_reuses_build_condition_summary() -> None:
    expression, metadata = _known_tables()
    metadata = metadata.assign(genotype=["mutant", "WT", "WT"])
    result = compute_sample_qc(expression, metadata)

    grouped = build_grouped_condition_summary(result, metadata, "genotype")

    assert list(grouped.columns) == ["genotype", "sample_count", "sample_ids"]
    by_genotype = grouped.set_index("genotype")
    assert by_genotype.loc["WT", "sample_count"] == 2
    assert sorted(by_genotype.loc["WT", "sample_ids"]) == ["s2", "s3"]
    assert by_genotype.loc["mutant", "sample_count"] == 1


def test_build_grouped_condition_chart_data_uses_the_group_column_name() -> None:
    expression, metadata = _known_tables()
    metadata = metadata.assign(genotype=["mutant", "WT", "WT"])
    result = compute_sample_qc(expression, metadata)

    chart_data = build_grouped_condition_chart_data(result, metadata, "genotype")

    assert list(chart_data.columns) == ["genotype", "sample_count"]


def test_build_grouped_summary_rejects_unknown_column() -> None:
    expression, metadata = _known_tables()
    result = compute_sample_qc(expression, metadata)

    with pytest.raises(ValueError, match="does not contain column"):
        build_grouped_sample_summary(result, metadata, "tissue")


def test_single_condition_and_one_replicate_observations_reuse_report() -> None:
    expression, metadata = _simple_tables()
    metadata["condition"] = ["one", "one"]
    report = ValidationReport(
        (
            ValidationIssue(
                code=IssueCode.SINGLE_CONDITION,
                severity=Severity.WARNING,
                table="Sample metadata",
                message="Existing single-condition observation.",
            ),
            ValidationIssue(
                code=IssueCode.CONDITION_WITH_SINGLE_SAMPLE,
                severity=Severity.WARNING,
                table="Sample metadata",
                message="Existing limited-replication observation.",
            ),
        )
    )

    result = compute_sample_qc(expression, metadata)
    observations = build_qc_observations(result, report)

    assert result.condition_count == 1
    assert result.condition_summary["sample_count"].tolist() == [2]
    assert observations.count(
        "Existing single-condition observation. This may warrant further review."
    ) == 1
    assert observations.count(
        "Existing limited-replication observation. "
        "This may warrant further review."
    ) == 1


def test_empty_expression_table_is_a_controlled_failure() -> None:
    expression = pd.DataFrame(columns=["gene_id", "s1", "s2"])
    _, metadata = _simple_tables()

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.EMPTY_EXPRESSION_TABLE


def test_zero_sample_columns_is_a_controlled_failure() -> None:
    expression = pd.DataFrame({"gene_id": ["g1"]})
    metadata = pd.DataFrame(
        {"sample_id": [], "condition": []},
        dtype="object",
    )

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS


def test_non_coercible_expression_value_is_a_controlled_failure() -> None:
    expression, metadata = _simple_tables()
    expression["s1"] = expression["s1"].astype("object")
    expression.loc[0, "s1"] = "not-a-number"

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.NON_COERCIBLE_EXPRESSION_VALUE
    assert "1 non-coercible" in str(raised.value)


@pytest.mark.parametrize("missing_value", [None, pd.NA, "", "   "])
def test_missing_expression_value_is_a_controlled_failure(
    missing_value: object,
) -> None:
    expression, metadata = _simple_tables()
    expression["s1"] = expression["s1"].astype("object")
    expression.loc[0, "s1"] = missing_value

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.MISSING_EXPRESSION_VALUE
    assert "1 missing or blank" in str(raised.value)


@pytest.mark.parametrize(
    ("infinite_value", "expected_counts"),
    [
        (float("inf"), "1 positive and 0 negative"),
        (float("-inf"), "0 positive and 1 negative"),
    ],
)
def test_infinite_expression_values_are_controlled_failures(
    infinite_value: float,
    expected_counts: str,
) -> None:
    expression, metadata = _simple_tables()
    expression.loc[0, "s1"] = infinite_value

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.INFINITE_EXPRESSION_VALUE
    assert expected_counts in str(raised.value)


@pytest.mark.parametrize(
    ("table_name", "column"),
    [
        ("expression", "gene_id"),
        ("metadata", "sample_id"),
        ("metadata", "condition"),
    ],
)
def test_missing_required_columns_are_controlled_failures(
    table_name: str,
    column: str,
) -> None:
    expression, metadata = _simple_tables()
    if table_name == "expression":
        expression = expression.drop(columns=column)
    else:
        metadata = metadata.drop(columns=column)

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.MISSING_REQUIRED_COLUMN


def test_duplicate_required_column_is_a_controlled_failure() -> None:
    expression = pd.DataFrame(
        [
            ["g1", "g1-copy", 1.0, 2.0],
            ["g2", "g2-copy", 3.0, 4.0],
        ],
        columns=["gene_id", "gene_id", "s1", "s2"],
    )
    _, metadata = _simple_tables()

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.DUPLICATE_REQUIRED_COLUMN


def test_duplicate_expression_sample_identifier_is_a_controlled_failure() -> None:
    expression = pd.DataFrame(
        [
            ["g1", 1.0, 2.0],
            ["g2", 3.0, 4.0],
        ],
        columns=["gene_id", "s1", "s1"],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1"], "condition": ["control"]}
    )

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.DUPLICATE_REQUIRED_IDENTIFIER


def test_missing_expression_sample_identifier_is_a_controlled_failure() -> None:
    expression = pd.DataFrame(
        [
            ["g1", 1.0, 2.0],
            ["g2", 3.0, 4.0],
        ],
        columns=["gene_id", None, "s2"],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["missing", "s2"], "condition": ["a", "b"]}
    )

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.MISSING_REQUIRED_IDENTIFIER


def test_missing_metadata_identifier_is_a_controlled_failure() -> None:
    expression, metadata = _simple_tables()
    metadata.loc[0, "sample_id"] = None

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.MISSING_REQUIRED_IDENTIFIER


def test_missing_metadata_condition_is_a_controlled_failure() -> None:
    expression, metadata = _simple_tables()
    metadata.loc[0, "condition"] = " "

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.MISSING_REQUIRED_VALUE


def test_duplicate_metadata_identifier_is_a_controlled_failure() -> None:
    expression, metadata = _simple_tables()
    metadata.loc[1, "sample_id"] = "s1"

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.DUPLICATE_REQUIRED_IDENTIFIER


@pytest.mark.parametrize(
    ("gene_ids", "reason"),
    [
        (["g1", None], QcErrorReason.MISSING_REQUIRED_IDENTIFIER),
        (["g1", "g1"], QcErrorReason.DUPLICATE_REQUIRED_IDENTIFIER),
    ],
)
def test_invalid_gene_identifiers_are_controlled_failures(
    gene_ids: list[object],
    reason: QcErrorReason,
) -> None:
    expression, metadata = _simple_tables()
    expression["gene_id"] = gene_ids

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is reason


def test_metadata_expression_mismatch_is_a_controlled_failure() -> None:
    expression, metadata = _simple_tables()
    metadata.loc[1, "sample_id"] = "other"

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is QcErrorReason.SAMPLE_MISMATCH
    assert "missing metadata IDs: s2" in str(raised.value)
    assert "metadata-only IDs: other" in str(raised.value)


def test_metadata_order_difference_is_accepted_and_observed_once() -> None:
    expression, metadata = _simple_tables()
    metadata = metadata.iloc[::-1].copy()

    result = compute_sample_qc(expression, metadata)
    report = ValidationReport(
        (
            ValidationIssue(
                code=IssueCode.SAMPLE_ORDER_DIFFERS,
                severity=Severity.INFORMATION,
                table="Expression matrix",
                message="Existing order observation.",
            ),
        )
    )
    observations = build_qc_observations(result, report)

    assert result.metadata_order_matches_expression is False
    assert result.sample_summary["sample_id"].tolist() == ["s1", "s2"]
    assert sum("order differs" in item for item in observations) == 1
    assert all("Existing order observation" not in item for item in observations)


def test_repeated_calculation_produces_equivalent_new_results() -> None:
    expression, metadata = _known_tables()

    first = compute_sample_qc(expression, metadata)
    second = compute_sample_qc(expression, metadata)

    assert first is not second
    assert first.sample_summary is not second.sample_summary
    assert first.condition_summary is not second.condition_summary
    assert_frame_equal(first.sample_summary, second.sample_summary)
    assert_frame_equal(first.condition_summary, second.condition_summary)
    assert first.constant_samples == second.constant_samples
    assert first.zero_variance_gene_count == second.zero_variance_gene_count


def test_one_gene_standard_deviation_and_observation_are_exact() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "s1": [2.0],
            "s2": [3.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["s1", "s2"],
            "condition": ["a", "b"],
        }
    )

    result = compute_sample_qc(expression, metadata)
    observations = build_qc_observations(result, ValidationReport())

    assert result.gene_count == 1
    assert result.sample_summary["standard_deviation"].isna().all()
    assert result.missing_value_count == 0
    assert result.zero_variance_gene_count == 0
    assert result.constant_samples == ("s1", "s2")
    assert observations[0] == (
        "The matrix contains only one gene, so sample-spread and "
        "constant-column diagnostics are not informative."
    )
    assert all("constant sample column" not in item for item in observations)


def test_one_sample_zero_variance_gene_boundary_is_explicit() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "only_sample": [1.0, 2.0, 3.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["only_sample"],
            "condition": ["control"],
        }
    )
    numeric_expression = expression.drop(columns="gene_id")

    assert find_constant_samples(numeric_expression) == ()
    assert count_zero_variance_genes(numeric_expression) == 3

    result = compute_sample_qc(expression, metadata)

    assert result.constant_samples == ()
    assert result.zero_variance_gene_count == 3
    assert result.gene_count == 3
    assert result.sample_count == 1


@pytest.mark.parametrize(
    ("mutator", "reason"),
    [
        (
            lambda expression, metadata: expression.__setitem__(
                "s1", ["bad", 2.0]
            ),
            QcErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
        ),
        (
            lambda expression, metadata: metadata.__setitem__(
                "sample_id", ["s1", "other"]
            ),
            QcErrorReason.SAMPLE_MISMATCH,
        ),
    ],
)
def test_inputs_are_not_mutated_when_calculation_fails(
    mutator,
    reason: QcErrorReason,
) -> None:
    expression, metadata = _simple_tables()
    mutator(expression, metadata)
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    with pytest.raises(QcComputationError) as raised:
        compute_sample_qc(expression, metadata)

    assert raised.value.reason is reason
    _assert_inputs_unchanged(
        expression,
        metadata,
        original_expression,
        original_metadata,
    )


def test_chart_data_preserves_sample_and_condition_order() -> None:
    expression, metadata = _known_tables()
    result = compute_sample_qc(expression, metadata)

    sample_chart = build_sample_chart_data(result, "median")
    condition_chart = build_condition_chart_data(result)

    assert sample_chart["sample_id"].tolist() == ["s2", "s1", "s3"]
    assert condition_chart["condition"].tolist() == ["control", "treated"]
    assert sample_chart is not result.sample_summary
    assert condition_chart is not result.condition_summary


def test_count_chart_decisions_reflect_non_zero_values() -> None:
    expression, metadata = _known_tables()
    result = compute_sample_qc(expression, metadata)
    positive_expression = expression.copy(deep=True)
    positive_expression.loc[:, ["s2", "s1", "s3"]] = [
        [1.0, 2.0, 3.0],
        [2.0, 3.0, 4.0],
        [3.0, 4.0, 5.0],
        [4.0, 5.0, 6.0],
    ]
    all_positive = compute_sample_qc(positive_expression, metadata)

    assert should_display_count_chart(result, "zero_value_count") is True
    assert should_display_count_chart(result, "negative_value_count") is True
    assert (
        should_display_count_chart(all_positive, "zero_value_count") is False
    )
    assert (
        should_display_count_chart(all_positive, "negative_value_count")
        is False
    )


def test_pure_helpers_do_not_create_session_state_or_heuristic_labels() -> None:
    expression, metadata = _known_tables()
    unrelated_state = {"keep": "unchanged"}

    result = compute_sample_qc(expression, metadata)
    build_sample_chart_data(result, "interquartile_range")
    build_condition_chart_data(result)
    observations = build_qc_observations(result, ValidationReport())

    assert unrelated_state == {"keep": "unchanged"}
    observation_text = " ".join(observations).lower()
    assert "outlier" not in observation_text
    assert "quality score" not in observation_text
    assert "pass" not in observation_text
    assert "fail" not in observation_text
