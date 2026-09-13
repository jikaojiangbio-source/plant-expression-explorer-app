"""Tests for descriptive sample-to-sample Pearson correlations."""

from dataclasses import fields

import numpy as np
import pandas as pd
import pytest

import plant_expression_explorer.correlation as correlation_module
from plant_expression_explorer.correlation import (
    CONDITION_CORRELATION_SUMMARY_COLUMNS,
    HEATMAP_DATA_COLUMNS,
    PAIR_SUMMARY_COLUMNS,
    SAMPLE_CORRELATION_SUMMARY_COLUMNS,
    CorrelationComputationError,
    CorrelationErrorReason,
    SampleCorrelationResult,
    build_correlation_observations,
    build_grouped_condition_correlation_summary,
    build_grouped_pair_summary,
    build_grouped_sample_correlation_summary,
    build_heatmap_data,
    compute_sample_correlation,
    order_samples_by_group,
)
from plant_expression_explorer.dataset import load_demo_candidate
from plant_expression_explorer.qc import list_additional_metadata_columns
from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
)


@pytest.fixture
def asymmetric_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return five gene rows and four deliberately reordered samples."""

    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4", "g5"],
            "sample_b": [1, 2, 3, 4, 5],
            "sample_a": [2, 4, 6, 8, 10],
            "sample_c": [5, 4, 3, 2, 1],
            "sample_d": [2, 1, 4, 3, 5],
        },
        index=[11, 13, 17, 19, 23],
    )
    metadata = pd.DataFrame(
        {
            "sample_id": [
                "sample_d",
                "sample_a",
                "sample_c",
                "sample_b",
            ],
            "condition": ["treated", "control", "treated", "control"],
        },
        index=[8, 6, 4, 2],
    )
    return expression, metadata


def _constant_tables(
    *,
    constant_count: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4", "g5"],
            "sample_1": [1, 2, 3, 4, 5],
            "sample_2": [2, 4, 6, 8, 10],
            "sample_3": [5, 4, 3, 2, 1],
        }
    )
    metadata_rows = [
        {"sample_id": "sample_3", "condition": "treated"},
        {"sample_id": "sample_1", "condition": "control"},
        {"sample_id": "sample_2", "condition": "control"},
    ]
    for position in range(constant_count):
        sample_id = f"constant_{position + 1}"
        expression[sample_id] = 7 + position
        metadata_rows.insert(
            1,
            {"sample_id": sample_id, "condition": "treated"},
        )
    return expression, pd.DataFrame(metadata_rows)


def _all_constant_tables(
    sample_count: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression_data: dict[str, list[object]] = {
        "gene_id": ["g1", "g2", "g3", "g4"]
    }
    metadata_rows: list[dict[str, str]] = []
    for position in range(sample_count):
        sample_id = f"sample_{position + 1}"
        expression_data[sample_id] = [position + 1] * 4
        metadata_rows.append(
            {
                "sample_id": sample_id,
                "condition": "one_condition",
            }
        )
    return pd.DataFrame(expression_data), pd.DataFrame(metadata_rows)


def _assert_reason(
    expression: object,
    metadata: object,
    expected: CorrelationErrorReason,
) -> CorrelationComputationError:
    with pytest.raises(CorrelationComputationError) as raised:
        compute_sample_correlation(expression, metadata)  # type: ignore[arg-type]
    assert raised.value.reason is expected
    assert str(raised.value)
    assert "Traceback" not in str(raised.value)
    return raised.value


def test_error_reason_values_are_unique() -> None:
    values = [reason.value for reason in CorrelationErrorReason]
    assert len(values) == len(set(values))


def test_core_module_has_no_streamlit_or_session_state_dependency() -> None:
    assert "streamlit" not in correlation_module.__dict__
    assert "st" not in correlation_module.__dict__
    assert "session_state" not in correlation_module.__dict__


@pytest.mark.parametrize(
    ("expression", "metadata"),
    [
        ([], pd.DataFrame()),
        (pd.DataFrame(), []),
    ],
)
def test_non_dataframe_input_has_controlled_reason(
    expression: object,
    metadata: object,
) -> None:
    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.INVALID_TABLE,
    )


def test_demo_matrix_is_six_by_six_and_pearson() -> None:
    candidate = load_demo_candidate()
    assert candidate.tables is not None
    expression, metadata, _ = candidate.tables

    result = compute_sample_correlation(expression, metadata)

    assert result.method == "pearson"
    assert result.gene_count == 120
    assert result.sample_count == 6
    assert result.correlation_matrix.shape == (6, 6)


def test_axes_known_values_symmetry_diagonal_and_bounds(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_correlation(expression, metadata)
    matrix = result.correlation_matrix

    expected_order = ["sample_b", "sample_a", "sample_c", "sample_d"]
    assert list(matrix.index) == expected_order
    assert list(matrix.columns) == expected_order
    assert matrix.loc["sample_b", "sample_a"] == pytest.approx(1.0)
    assert matrix.loc["sample_b", "sample_d"] == pytest.approx(0.8)
    assert matrix.loc["sample_b", "sample_c"] == pytest.approx(-1.0)
    assert matrix.loc["sample_c", "sample_d"] == pytest.approx(-0.8)
    assert np.allclose(
        matrix.to_numpy(),
        matrix.to_numpy().T,
        equal_nan=True,
    )
    assert np.allclose(np.diag(matrix), np.ones(4))
    defined = matrix.to_numpy()[~np.isnan(matrix.to_numpy())]
    assert np.all(defined >= -1.0 - 1e-12)
    assert np.all(defined <= 1.0 + 1e-12)


def test_numeric_strings_use_temporary_copy_without_input_mutation(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    for column in expression.columns[1:]:
        expression[column] = expression[column].map(str)
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)
    expression_dtypes = expression.dtypes.copy(deep=True)
    metadata_dtypes = metadata.dtypes.copy(deep=True)

    result = compute_sample_correlation(expression, metadata)

    assert result.correlation_matrix.loc[
        "sample_b",
        "sample_d",
    ] == pytest.approx(0.8)
    pd.testing.assert_frame_equal(expression, original_expression)
    pd.testing.assert_frame_equal(metadata, original_metadata)
    pd.testing.assert_series_equal(expression.dtypes, expression_dtypes)
    pd.testing.assert_series_equal(metadata.dtypes, metadata_dtypes)


def test_metadata_maps_by_exact_id_and_reports_reordered_input(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_correlation(expression, metadata)

    assert result.metadata_order_matches_expression is False
    assert result.sample_summary["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
        "sample_c",
        "sample_d",
    ]
    assert result.sample_summary["condition"].tolist() == [
        "control",
        "control",
        "treated",
        "treated",
    ]


def test_list_additional_metadata_columns_excludes_required_columns() -> None:
    metadata = pd.DataFrame(
        {"sample_id": ["s1"], "condition": ["control"], "genotype": ["WT"]}
    )

    assert list_additional_metadata_columns(metadata) == ("genotype",)


def test_build_grouped_summaries_match_defaults_for_condition(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_correlation(expression, metadata)

    pd.testing.assert_frame_equal(
        build_grouped_pair_summary(result, metadata, "condition"),
        result.pair_summary,
    )
    pd.testing.assert_frame_equal(
        build_grouped_sample_correlation_summary(result, metadata, "condition"),
        result.sample_summary,
    )
    pd.testing.assert_frame_equal(
        build_grouped_condition_correlation_summary(result, metadata, "condition"),
        result.condition_summary,
    )


def test_build_grouped_summaries_relabel_by_alternate_metadata_column(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    expected_genotype_by_sample = {
        "sample_d": "mutant",
        "sample_a": "WT",
        "sample_c": "mutant",
        "sample_b": "WT",
    }
    metadata = metadata.assign(
        genotype=[
            expected_genotype_by_sample[sample_id]
            for sample_id in metadata["sample_id"]
        ]
    )
    result = compute_sample_correlation(expression, metadata)

    grouped_pairs = build_grouped_pair_summary(result, metadata, "genotype")
    assert "genotype_a" in grouped_pairs.columns
    assert "genotype_b" in grouped_pairs.columns
    assert "condition_a" not in grouped_pairs.columns
    for _, row in grouped_pairs.iterrows():
        assert row["genotype_a"] == expected_genotype_by_sample[row["sample_a"]]
        assert row["genotype_b"] == expected_genotype_by_sample[row["sample_b"]]

    grouped_samples = build_grouped_sample_correlation_summary(
        result, metadata, "genotype"
    )
    assert "genotype" in grouped_samples.columns
    for sample_id, genotype in zip(
        grouped_samples["sample_id"], grouped_samples["genotype"], strict=True
    ):
        assert genotype == expected_genotype_by_sample[sample_id]
    # Correlation statistics are untouched by the relabelling.
    for column in (
        "minimum_correlation",
        "median_correlation",
        "mean_correlation",
        "maximum_correlation",
    ):
        pd.testing.assert_series_equal(
            grouped_samples[column],
            result.sample_summary[column],
            check_names=False,
        )

    grouped_conditions = build_grouped_condition_correlation_summary(
        result, metadata, "genotype"
    )
    assert "genotype_a" in grouped_conditions.columns
    assert "genotype_b" in grouped_conditions.columns
    by_group = grouped_conditions.set_index(["genotype_a", "genotype_b"])
    assert by_group.loc[("WT", "WT"), "sample_count_a"] == 2
    assert by_group.loc[("mutant", "mutant"), "sample_count_a"] == 2


def test_build_grouped_summary_rejects_unknown_column(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_correlation(expression, metadata)

    with pytest.raises(ValueError, match="does not contain column"):
        build_grouped_pair_summary(result, metadata, "tissue")


def test_order_samples_by_group_moves_identical_labels_together() -> None:
    sample_ids = ["s1", "s2", "s3", "s4"]
    group_by_sample = {"s1": "treated", "s2": "control", "s3": "treated", "s4": "control"}

    ordered = order_samples_by_group(sample_ids, group_by_sample)

    assert ordered == ["s2", "s4", "s1", "s3"]


def test_order_samples_by_group_preserves_relative_order_within_a_group() -> None:
    sample_ids = ["s3", "s1", "s4", "s2"]
    group_by_sample = {"s1": "A", "s2": "A", "s3": "B", "s4": "B"}

    ordered = order_samples_by_group(sample_ids, group_by_sample)

    assert ordered == ["s1", "s2", "s3", "s4"]


def test_matching_metadata_order_is_recorded(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.set_index("sample_id").loc[
        expression.columns[1:]
    ]
    metadata.index.name = "sample_id"
    metadata = metadata.reset_index()

    result = compute_sample_correlation(expression, metadata)

    assert result.metadata_order_matches_expression is True


def test_missing_gene_id_column_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    _assert_reason(
        expression.drop(columns="gene_id"),
        metadata,
        CorrelationErrorReason.MISSING_REQUIRED_COLUMN,
    )


def test_duplicate_gene_id_column_is_controlled() -> None:
    expression = pd.DataFrame(
        [["g1", "copy1", 1], ["g2", "copy2", 2]],
        columns=["gene_id", "gene_id", "sample_1"],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_1"], "condition": ["control"]}
    )

    error = _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.DUPLICATE_REQUIRED_COLUMN,
    )

    assert "2 columns named 'gene_id'" in str(error)


@pytest.mark.parametrize("bad_gene_id", ["", "   ", None, pd.NA])
def test_blank_or_missing_gene_id_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    bad_gene_id: object,
) -> None:
    expression, metadata = asymmetric_tables
    expression = expression.copy(deep=True)
    expression.iloc[2, expression.columns.get_loc("gene_id")] = bad_gene_id

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.MISSING_REQUIRED_IDENTIFIER,
    )


def test_duplicate_gene_id_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    expression = expression.copy(deep=True)
    expression.loc[expression.index[1], "gene_id"] = "g1"

    error = _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
    )

    assert "g1" in str(error)


@pytest.mark.parametrize("bad_identifier", ["", "   ", None])
def test_missing_sample_column_identifier_is_controlled(
    bad_identifier: object,
) -> None:
    expression = pd.DataFrame(
        [["g1", 1], ["g2", 2]],
        columns=["gene_id", bad_identifier],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_1"], "condition": ["control"]}
    )

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.MISSING_REQUIRED_IDENTIFIER,
    )


def test_duplicate_sample_column_identifier_is_controlled() -> None:
    expression = pd.DataFrame(
        [["g1", 1, 2], ["g2", 2, 3]],
        columns=["gene_id", "sample_1", "sample_1"],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_1"], "condition": ["control"]}
    )

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
    )


def test_one_constant_sample_has_nan_axis_and_n_minus_one_undefined_pairs() -> None:
    expression, metadata = _constant_tables()

    result = compute_sample_correlation(expression, metadata)

    assert result.constant_samples == ("constant_1",)
    assert result.correlation_matrix.loc["constant_1"].isna().all()
    assert result.correlation_matrix.loc[:, "constant_1"].isna().all()
    assert result.undefined_pair_count == result.sample_count - 1
    undefined = result.pair_summary.loc[~result.pair_summary["defined"]]
    assert len(undefined.index) == result.sample_count - 1
    assert undefined["correlation"].isna().all()


def test_multiple_constant_samples_remain_nan() -> None:
    expression, metadata = _constant_tables(constant_count=2)

    result = compute_sample_correlation(expression, metadata)

    assert result.constant_samples == ("constant_1", "constant_2")
    for sample_id in result.constant_samples:
        assert result.correlation_matrix.loc[sample_id].isna().all()
        assert result.correlation_matrix.loc[:, sample_id].isna().all()


def test_all_samples_constant_contract() -> None:
    expression, metadata = _all_constant_tables(sample_count=4)

    result = compute_sample_correlation(expression, metadata)

    assert result.constant_samples == (
        "sample_1",
        "sample_2",
        "sample_3",
        "sample_4",
    )
    assert result.correlation_matrix.isna().all().all()
    assert len(result.pair_summary.index) == 6
    assert result.undefined_pair_count == 6
    assert (~result.pair_summary["defined"]).all()
    assert result.pair_summary["correlation"].isna().all()
    assert result.sample_summary["defined_pair_count"].eq(0).all()
    assert result.sample_summary["undefined_pair_count"].eq(3).all()
    statistic_columns = SAMPLE_CORRELATION_SUMMARY_COLUMNS[4:]
    assert result.sample_summary.loc[:, statistic_columns].isna().all().all()


@pytest.mark.parametrize(
    ("values", "expected_constant", "expected_diagonal"),
    [
        ([1, 3, 2], (), 1.0),
        ([4, 4, 4], ("only_sample",), float("nan")),
    ],
)
def test_one_sample_contract(
    values: list[int],
    expected_constant: tuple[str, ...],
    expected_diagonal: float,
) -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2", "g3"], "only_sample": values}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["only_sample"], "condition": ["control"]}
    )

    result = compute_sample_correlation(expression, metadata)

    assert result.correlation_matrix.shape == (1, 1)
    if pd.isna(expected_diagonal):
        assert pd.isna(result.correlation_matrix.iat[0, 0])
    else:
        assert result.correlation_matrix.iat[0, 0] == pytest.approx(
            expected_diagonal
        )
    assert result.constant_samples == expected_constant
    assert result.pair_summary.empty
    assert result.undefined_pair_count == 0
    summary = result.sample_summary.iloc[0]
    assert summary["defined_pair_count"] == 0
    assert summary["undefined_pair_count"] == 0
    assert summary.loc[
        list(SAMPLE_CORRELATION_SUMMARY_COLUMNS[4:])
    ].isna().all()


def test_one_gene_is_controlled() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "sample_1": [1], "sample_2": [2]}
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "condition": ["control", "treated"],
        }
    )

    error = _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.INSUFFICIENT_GENE_ROWS,
    )

    assert "1 gene row" in str(error)
    assert "at least 2" in str(error)


def test_zero_genes_is_controlled() -> None:
    expression = pd.DataFrame(columns=["gene_id", "sample_1"])
    metadata = pd.DataFrame(
        {"sample_id": ["sample_1"], "condition": ["control"]}
    )

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.EMPTY_EXPRESSION_TABLE,
    )


def test_no_sample_columns_is_controlled() -> None:
    expression = pd.DataFrame({"gene_id": ["g1", "g2"]})
    metadata = pd.DataFrame(columns=["sample_id", "condition"])

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
    )


def test_multiple_invalid_inputs_follow_fixed_error_priority() -> None:
    valid_metadata = pd.DataFrame(
        {"sample_id": ["sample_1"], "condition": ["control"]}
    )

    _assert_reason(
        pd.DataFrame(),
        [],
        CorrelationErrorReason.INVALID_TABLE,
    )
    _assert_reason(
        pd.DataFrame(),
        pd.DataFrame(),
        CorrelationErrorReason.EMPTY_EXPRESSION_TABLE,
    )
    _assert_reason(
        pd.DataFrame({"sample_1": [1]}),
        valid_metadata,
        CorrelationErrorReason.MISSING_REQUIRED_COLUMN,
    )
    _assert_reason(
        pd.DataFrame({"gene_id": ["g1"]}),
        valid_metadata,
        CorrelationErrorReason.INSUFFICIENT_GENE_ROWS,
    )
    _assert_reason(
        pd.DataFrame({"gene_id": ["", "g2"]}),
        valid_metadata,
        CorrelationErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
    )

    expression_with_invalid_value = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [True, 2]}
    )
    _assert_reason(
        expression_with_invalid_value,
        pd.DataFrame({"sample_id": ["sample_1"]}),
        CorrelationErrorReason.MISSING_REQUIRED_COLUMN,
    )
    mismatch_metadata = pd.DataFrame(
        {"sample_id": ["metadata_only"], "condition": ["control"]}
    )
    expression_with_missing_value = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [None, True]}
    )
    _assert_reason(
        expression_with_missing_value,
        mismatch_metadata,
        CorrelationErrorReason.SAMPLE_MISMATCH,
    )
    missing_error = _assert_reason(
        expression_with_missing_value,
        valid_metadata,
        CorrelationErrorReason.MISSING_EXPRESSION_VALUE,
    )
    assert "1 missing or blank" in str(missing_error)

    boolean_before_non_numeric = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [True, "not-numeric"]}
    )
    boolean_error = _assert_reason(
        boolean_before_non_numeric,
        valid_metadata,
        CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
    )
    assert "1 boolean and 0 complex" in str(boolean_error)


@pytest.mark.parametrize(
    ("bad_value", "reason", "message_part"),
    [
        (
            "not-numeric",
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "1 non-coercible",
        ),
        (
            None,
            CorrelationErrorReason.MISSING_EXPRESSION_VALUE,
            "1 missing or blank",
        ),
        (
            "   ",
            CorrelationErrorReason.MISSING_EXPRESSION_VALUE,
            "1 missing or blank",
        ),
        (
            True,
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "1 boolean and 0 complex",
        ),
        (
            np.bool_(False),
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "1 boolean and 0 complex",
        ),
        (
            1 + 2j,
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "0 boolean and 1 complex",
        ),
        (
            np.complex128(3 + 4j),
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "0 boolean and 1 complex",
        ),
        (
            float("inf"),
            CorrelationErrorReason.INFINITE_EXPRESSION_VALUE,
            "1 positive and 0 negative",
        ),
        (
            float("-inf"),
            CorrelationErrorReason.INFINITE_EXPRESSION_VALUE,
            "0 positive and 1 negative",
        ),
    ],
)
def test_invalid_expression_value_has_specific_controlled_reason(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    bad_value: object,
    reason: CorrelationErrorReason,
    message_part: str,
) -> None:
    expression, metadata = asymmetric_tables
    expression = expression.copy(deep=True)
    expression["sample_b"] = expression["sample_b"].astype("object")
    expression.loc[expression.index[0], "sample_b"] = bad_value

    error = _assert_reason(expression, metadata, reason)

    assert message_part in str(error)


@pytest.mark.parametrize("bad_condition", ["", "   ", None, pd.NA])
def test_missing_or_blank_condition_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    bad_condition: object,
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.copy(deep=True)
    metadata.loc[metadata.index[0], "condition"] = bad_condition

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.MISSING_REQUIRED_VALUE,
    )


@pytest.mark.parametrize("bad_sample_id", ["", "   ", None, pd.NA])
def test_missing_or_blank_metadata_sample_id_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    bad_sample_id: object,
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.copy(deep=True)
    metadata.loc[metadata.index[0], "sample_id"] = bad_sample_id

    _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.MISSING_REQUIRED_IDENTIFIER,
    )


def test_duplicate_metadata_id_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.copy(deep=True)
    metadata.loc[metadata.index[0], "sample_id"] = "sample_a"

    error = _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
    )

    assert "sample_a" in str(error)


@pytest.mark.parametrize("column", ["sample_id", "condition"])
def test_missing_metadata_required_column_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    column: str,
) -> None:
    expression, metadata = asymmetric_tables

    _assert_reason(
        expression,
        metadata.drop(columns=column),
        CorrelationErrorReason.MISSING_REQUIRED_COLUMN,
    )


@pytest.mark.parametrize("column", ["sample_id", "condition"])
def test_duplicate_metadata_required_column_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    column: str,
) -> None:
    expression, metadata = asymmetric_tables
    duplicated = pd.concat([metadata, metadata[[column]]], axis=1)

    _assert_reason(
        expression,
        duplicated,
        CorrelationErrorReason.DUPLICATE_REQUIRED_COLUMN,
    )


def test_metadata_expression_mismatch_has_both_directions(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.copy(deep=True)
    metadata.loc[metadata["sample_id"].eq("sample_d"), "sample_id"] = (
        "metadata_only"
    )

    error = _assert_reason(
        expression,
        metadata,
        CorrelationErrorReason.SAMPLE_MISMATCH,
    )

    assert "missing metadata IDs: sample_d" in str(error)
    assert "metadata-only IDs: metadata_only" in str(error)


def test_pair_summary_columns_count_order_uniqueness_and_labels(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    pairs = compute_sample_correlation(expression, metadata).pair_summary

    assert tuple(pairs.columns) == PAIR_SUMMARY_COLUMNS
    assert len(pairs.index) == 4 * 3 // 2
    assert list(zip(pairs["sample_a"], pairs["sample_b"], strict=True)) == [
        ("sample_b", "sample_a"),
        ("sample_b", "sample_c"),
        ("sample_b", "sample_d"),
        ("sample_a", "sample_c"),
        ("sample_a", "sample_d"),
        ("sample_c", "sample_d"),
    ]
    assert not (pairs["sample_a"] == pairs["sample_b"]).any()
    unordered = pairs.apply(
        lambda row: frozenset((row["sample_a"], row["sample_b"])),
        axis=1,
    )
    assert len(set(unordered)) == len(pairs.index)
    assert pairs["relationship"].tolist() == [
        "within_condition",
        "between_condition",
        "between_condition",
        "between_condition",
        "between_condition",
        "within_condition",
    ]


def test_undefined_pair_count_uses_unique_pair_rows_not_square_cells() -> None:
    expression, metadata = _constant_tables()

    result = compute_sample_correlation(expression, metadata)
    square_nan_count = int(result.correlation_matrix.isna().sum().sum())

    assert result.undefined_pair_count == 3
    assert square_nan_count == 7
    assert result.undefined_pair_count == int(
        (~result.pair_summary["defined"]).sum()
    )


def test_sample_summary_order_self_exclusion_and_defined_statistics(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    summary = compute_sample_correlation(expression, metadata).sample_summary

    assert tuple(summary.columns) == SAMPLE_CORRELATION_SUMMARY_COLUMNS
    assert summary["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
        "sample_c",
        "sample_d",
    ]
    sample_b = summary.loc[summary["sample_id"].eq("sample_b")].iloc[0]
    assert sample_b["defined_pair_count"] == 3
    assert sample_b["undefined_pair_count"] == 0
    assert sample_b["minimum_correlation"] == pytest.approx(-1.0)
    assert sample_b["median_correlation"] == pytest.approx(0.8)
    assert sample_b["mean_correlation"] == pytest.approx(0.8 / 3)
    assert sample_b["maximum_correlation"] == pytest.approx(1.0)


def test_sample_summary_counts_undefined_but_excludes_it_from_statistics() -> None:
    expression, metadata = _constant_tables()

    summary = compute_sample_correlation(expression, metadata).sample_summary

    non_constant = summary.loc[summary["sample_id"].eq("sample_1")].iloc[0]
    assert non_constant["defined_pair_count"] == 2
    assert non_constant["undefined_pair_count"] == 1
    assert non_constant["minimum_correlation"] == pytest.approx(-1.0)
    assert non_constant["maximum_correlation"] == pytest.approx(1.0)
    constant = summary.loc[summary["sample_id"].eq("constant_1")].iloc[0]
    assert constant["defined_pair_count"] == 0
    assert constant["undefined_pair_count"] == 3
    assert constant.loc[
        list(SAMPLE_CORRELATION_SUMMARY_COLUMNS[4:])
    ].isna().all()


def test_condition_summary_order_formulas_counts_and_medians(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    summary = compute_sample_correlation(expression, metadata).condition_summary

    assert tuple(summary.columns) == CONDITION_CORRELATION_SUMMARY_COLUMNS
    assert list(
        zip(summary["condition_a"], summary["condition_b"], strict=True)
    ) == [
        ("control", "control"),
        ("control", "treated"),
        ("treated", "treated"),
    ]
    control = summary.iloc[0]
    assert control["relationship"] == "within_condition"
    assert control["total_pair_count"] == 2 * 1 // 2
    assert control["defined_pair_count"] == 1
    assert control["undefined_pair_count"] == 0
    assert control["median_correlation"] == pytest.approx(1.0)
    between = summary.iloc[1]
    assert between["relationship"] == "between_condition"
    assert between["total_pair_count"] == 2 * 2
    assert between["defined_pair_count"] == 4
    assert between["undefined_pair_count"] == 0
    assert between["median_correlation"] == pytest.approx(-0.1)
    treated = summary.iloc[2]
    assert treated["relationship"] == "within_condition"
    assert treated["median_correlation"] == pytest.approx(-0.8)
    assert (
        summary["defined_pair_count"] + summary["undefined_pair_count"]
    ).equals(summary["total_pair_count"])


def test_condition_median_excludes_undefined_pairs() -> None:
    expression, metadata = _constant_tables()

    summary = compute_sample_correlation(expression, metadata).condition_summary

    treated_within = summary.loc[
        summary["condition_a"].eq("treated")
        & summary["condition_b"].eq("treated")
    ].iloc[0]
    assert treated_within["total_pair_count"] == 1
    assert treated_within["defined_pair_count"] == 0
    assert treated_within["undefined_pair_count"] == 1
    assert pd.isna(treated_within["median_correlation"])
    between = summary.loc[
        summary["relationship"].eq("between_condition")
    ].iloc[0]
    assert between["total_pair_count"] == 4
    assert between["defined_pair_count"] == 2
    assert between["undefined_pair_count"] == 2
    assert between["median_correlation"] == pytest.approx(-1.0)


def test_one_sample_condition_retains_zero_total_within_row(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.copy(deep=True)
    metadata.loc[metadata["sample_id"].eq("sample_d"), "condition"] = "solo"

    summary = compute_sample_correlation(expression, metadata).condition_summary
    solo_within = summary.loc[
        summary["condition_a"].eq("solo")
        & summary["condition_b"].eq("solo")
    ].iloc[0]

    assert solo_within["sample_count_a"] == 1
    assert solo_within["sample_count_b"] == 1
    assert solo_within["total_pair_count"] == 0
    assert solo_within["defined_pair_count"] == 0
    assert solo_within["undefined_pair_count"] == 0
    assert pd.isna(solo_within["median_correlation"])


def test_heatmap_has_n_squared_cells_row_major_order_and_undefined_cells() -> None:
    expression, metadata = _constant_tables()
    result = compute_sample_correlation(expression, metadata)
    original_matrix = result.correlation_matrix.copy(deep=True)

    heatmap = build_heatmap_data(result)

    assert HEATMAP_DATA_COLUMNS == (
        "row_sample",
        "column_sample",
        "correlation",
        "defined",
        "correlation_label",
    )
    assert tuple(heatmap.columns) == HEATMAP_DATA_COLUMNS
    assert len(heatmap.index) == result.sample_count**2
    sample_order = list(result.correlation_matrix.columns)
    assert heatmap["row_sample"].tolist() == [
        row_sample for row_sample in sample_order for _ in sample_order
    ]
    assert heatmap["column_sample"].tolist() == sample_order * len(sample_order)
    undefined = heatmap.loc[~heatmap["defined"]]
    assert len(undefined.index) == 2 * result.sample_count - 1
    assert undefined["correlation"].isna().all()
    assert undefined["correlation_label"].eq("N/A").all()
    pd.testing.assert_frame_equal(
        result.correlation_matrix,
        original_matrix,
    )


def test_repeated_computation_is_equivalent(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    first = compute_sample_correlation(expression, metadata)
    second = compute_sample_correlation(expression, metadata)

    assert first.method == second.method
    assert first.gene_count == second.gene_count
    assert first.sample_count == second.sample_count
    assert first.constant_samples == second.constant_samples
    assert first.undefined_pair_count == second.undefined_pair_count
    pd.testing.assert_frame_equal(
        first.correlation_matrix,
        second.correlation_matrix,
    )
    pd.testing.assert_frame_equal(first.pair_summary, second.pair_summary)
    pd.testing.assert_frame_equal(first.sample_summary, second.sample_summary)
    pd.testing.assert_frame_equal(
        first.condition_summary,
        second.condition_summary,
    )


def test_observations_are_deterministic_deduplicated_and_threshold_free(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_correlation(expression, metadata)
    duplicate_message = "Condition information should be reviewed."
    report = ValidationReport(
        (
            ValidationIssue(
                code=IssueCode.SINGLE_CONDITION,
                severity=Severity.WARNING,
                table="Sample metadata",
                message=duplicate_message,
            ),
            ValidationIssue(
                code=IssueCode.SINGLE_CONDITION,
                severity=Severity.INFORMATION,
                table="Sample metadata",
                message=duplicate_message,
            ),
        )
    )

    observations = build_correlation_observations(result, report)
    joined = " ".join(observations).lower()

    assert isinstance(observations, tuple)
    assert all(isinstance(observation, str) for observation in observations)
    assert observations == build_correlation_observations(result, report)
    assert not any(duplicate_message in item for item in observations)
    assert sum("SINGLE_CONDITION" in item for item in observations) == 1
    assert any("exact sample id" in item.lower() for item in observations)
    assert any(
        "both positive and negative" in item.lower() for item in observations
    )
    for prohibited in (
        "threshold",
        "quality",
        "outlier",
        "exclusion",
        "exclude",
        "ranking",
        "ranked",
        "strong",
        "weak",
        "anomal",
        "passed",
        "failed",
    ):
        assert prohibited not in joined


def test_constant_and_one_sample_observations_are_exact() -> None:
    expression, metadata = _constant_tables()
    constant_result = compute_sample_correlation(expression, metadata)
    one_expression = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "only": [1, 2]}
    )
    one_metadata = pd.DataFrame(
        {"sample_id": ["only"], "condition": ["control"]}
    )
    one_result = compute_sample_correlation(one_expression, one_metadata)

    constant_observations = build_correlation_observations(
        constant_result,
        ValidationReport(),
    )
    one_observations = build_correlation_observations(
        one_result,
        ValidationReport(),
    )

    assert any("constant_1" in item for item in constant_observations)
    assert any("3 unique non-self" in item for item in constant_observations)
    assert any("Only one sample" in item for item in one_observations)


@pytest.mark.parametrize(
    ("second_values", "expected_phrase"),
    [
        ([2, 4, 6], "are positive"),
        ([3, 2, 1], "are negative"),
    ],
)
def test_observations_report_exact_defined_pair_sign(
    second_values: list[int],
    expected_phrase: str,
) -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_1": [1, 2, 3],
            "sample_2": second_values,
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2"],
            "condition": ["control", "treated"],
        }
    )
    result = compute_sample_correlation(expression, metadata)

    observations = build_correlation_observations(
        result,
        ValidationReport(),
    )

    assert any(expected_phrase in item for item in observations)


def test_observations_skip_irrelevant_validation_messages(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_correlation(expression, metadata)
    report = ValidationReport(
        (
            ValidationIssue(
                code=IssueCode.ZERO_ADJUSTED_P_VALUE,
                severity=Severity.INFORMATION,
                table="Differential-expression results",
                message="An unrelated adjusted-p-value observation.",
            ),
        )
    )

    observations = build_correlation_observations(result, report)

    assert not any(
        "adjusted-p-value" in observation for observation in observations
    )


def test_result_schema_has_no_classification_fields(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_correlation(expression, metadata)
    schema_text = " ".join(
        [
            *(field.name for field in fields(SampleCorrelationResult)),
            *result.pair_summary.columns,
            *result.sample_summary.columns,
            *result.condition_summary.columns,
        ]
    ).lower()

    for prohibited in (
        "threshold",
        "quality",
        "outlier",
        "exclusion",
        "pass",
        "fail",
        "rank",
    ):
        assert prohibited not in schema_text


def test_axis_asymmetric_fixture_cannot_confuse_genes_and_samples(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_correlation(expression, metadata)

    assert result.gene_count == 5
    assert result.sample_count == 4
    assert result.correlation_matrix.shape == (4, 4)
    assert len(build_heatmap_data(result).index) == 16
