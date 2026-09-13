"""Tests for descriptive sample principal component analysis."""

import warnings
from dataclasses import fields

import numpy as np
import pandas as pd
import pytest

import plant_expression_explorer.pca as pca_module
from plant_expression_explorer.pca import (
    SCORE_TABLE_FIXED_COLUMNS,
    VARIANCE_TABLE_COLUMNS,
    PcaComputationError,
    PcaErrorReason,
    SamplePcaResult,
    build_grouped_score_plot_data,
    build_pca_observations,
    build_score_plot_data,
    compute_sample_pca,
)
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


@pytest.fixture
def separated_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Three samples, two perfectly correlated genes, well-separated scores.

    Not symmetric, so the largest-absolute-score sample is unambiguous even
    accounting for floating-point rounding (unlike a perfectly symmetric
    fixture, where sub-ULP rounding can make the "tie-break" outcome
    unpredictable by hand).
    """

    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample_a": [1, 2],
            "sample_b": [2, 4],
            "sample_c": [5, 10],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    return expression, metadata


@pytest.fixture
def colinear_symmetric_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Three samples, two perfectly correlated genes, symmetric about the mean.

    Independently checkable by hand: sample_b sits exactly at the mean (score
    0), and sample_a/sample_c are equidistant on opposite sides (antisymmetric
    scores), regardless of which sign convention outcome the tie-break picks
    for "largest magnitude" between them.
    """

    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample_a": [1, 2],
            "sample_b": [2, 4],
            "sample_c": [3, 6],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    return expression, metadata


def _all_constant_tables(
    sample_count: int = 3,
    gene_count: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_ids = [f"g{i + 1}" for i in range(gene_count)]
    data: dict[str, list[object]] = {"gene_id": gene_ids}
    for position in range(sample_count):
        data[f"sample_{position + 1}"] = [10 * (g + 1) for g in range(gene_count)]
    expression = pd.DataFrame(data)
    metadata = pd.DataFrame(
        {
            "sample_id": [f"sample_{position + 1}" for position in range(sample_count)],
            "condition": ["one_condition"] * sample_count,
        }
    )
    return expression, metadata


def _some_constant_gene_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_a": [1, 5, 2],
            "sample_b": [3, 5, 8],
            "sample_c": [2, 5, 4],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    return expression, metadata


def _tied_subspace_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample_1": [1, 0],
            "sample_2": [-1, 0],
            "sample_3": [0, 1],
            "sample_4": [0, -1],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_1", "sample_2", "sample_3", "sample_4"],
            "condition": ["a", "a", "b", "b"],
        }
    )
    return expression, metadata


def _rank_deficient_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "sample_a": [1, 3, 5, 1],
            "sample_b": [2, 6, 1, 2],
            "sample_c": [3, 9, 9, 3],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    return expression, metadata


def _large_finite_value_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Finite but overflow-prone values (1e300-scale), genes non-constant."""

    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_a": [1e300, 2.0, 3.0],
            "sample_b": [-1e300, 4.0, 6.0],
            "sample_c": [0.0, 5.0, 9.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    return expression, metadata


def _tiny_finite_value_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Finite but underflow-prone values (1e-300-scale); every gene is
    non-constant by exact equality, yet the summed ddof=1 variance
    underflows to exactly 0.0 in float64."""

    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_a": [1e-300, 2e-300, 3e-300],
            "sample_b": [-1e-300, -2e-300, -3e-300],
            "sample_c": [0.0, 0.0, 0.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    return expression, metadata


def _assert_reason(
    expression: object,
    metadata: object,
    expected: PcaErrorReason,
) -> PcaComputationError:
    with pytest.raises(PcaComputationError) as raised:
        compute_sample_pca(expression, metadata)  # type: ignore[arg-type]
    assert raised.value.reason is expected
    assert str(raised.value)
    assert "Traceback" not in str(raised.value)
    return raised.value


def test_error_reason_values_are_unique() -> None:
    values = [reason.value for reason in PcaErrorReason]
    assert len(values) == len(set(values))


def test_core_module_has_no_streamlit_or_session_state_dependency() -> None:
    assert "streamlit" not in pca_module.__dict__
    assert "st" not in pca_module.__dict__
    assert "session_state" not in pca_module.__dict__


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
    _assert_reason(expression, metadata, PcaErrorReason.INVALID_TABLE)


def test_colinear_fixture_has_checkable_pca_properties(
    colinear_symmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = colinear_symmetric_tables

    result = compute_sample_pca(expression, metadata)

    assert result.component_count == 2
    scores = result.score_table.set_index("sample_id")["pc1"]
    assert scores["sample_a"] == pytest.approx(-scores["sample_c"], rel=1e-9)
    assert scores["sample_b"] == pytest.approx(0.0, abs=1e-9)
    assert abs(scores["sample_a"]) == pytest.approx(abs(scores["sample_c"]))
    assert abs(scores["sample_a"]) > 0
    variance = result.variance_table
    assert variance.loc[variance["component"] == 1, "explained_variance_ratio"].iloc[0] == pytest.approx(1.0)
    assert variance.loc[variance["component"] == 2, "explained_variance_ratio"].iloc[0] == pytest.approx(0.0)


def test_axis_asymmetric_fixture_cannot_confuse_genes_and_samples(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    assert result.gene_count == 5
    assert result.sample_count == 4
    assert result.component_count == 3
    assert tuple(result.score_table.columns) == (
        *SCORE_TABLE_FIXED_COLUMNS,
        "pc1",
        "pc2",
        "pc3",
    )
    assert len(result.score_table.index) == 4


def test_score_table_row_order_matches_expression_sample_columns(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    assert result.score_table["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
        "sample_c",
        "sample_d",
    ]


def test_metadata_maps_by_exact_id_and_reports_reordered_input(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    assert result.metadata_order_matches_expression is False
    assert result.score_table["condition"].tolist() == [
        "control",
        "control",
        "treated",
        "treated",
    ]


def test_matching_metadata_order_is_recorded(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.set_index("sample_id").loc[expression.columns[1:]]
    metadata.index.name = "sample_id"
    metadata = metadata.reset_index()

    result = compute_sample_pca(expression, metadata)

    assert result.metadata_order_matches_expression is True


def test_variance_table_columns_and_order(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    assert tuple(result.variance_table.columns) == VARIANCE_TABLE_COLUMNS
    assert result.variance_table["component"].tolist() == [1, 2, 3]
    singular_values = result.variance_table["singular_value"].to_numpy()
    assert np.all(np.diff(singular_values) <= 1e-9)


def test_explained_variance_matches_singular_value_formula(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    variance = result.variance_table
    expected = (variance["singular_value"] ** 2) / (result.sample_count - 1)
    assert variance["explained_variance"].to_numpy() == pytest.approx(
        expected.to_numpy()
    )


def test_reported_total_variance_close_to_raw_total_variance(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    sample_columns = [c for c in expression.columns if c != "gene_id"]
    raw_total_variance = float(
        expression.loc[:, sample_columns].astype(float).var(axis=1, ddof=1).sum()
    )

    result = compute_sample_pca(expression, metadata)

    assert result.total_variance == pytest.approx(raw_total_variance, rel=1e-9)


def test_ratios_are_finite_nonnegative_and_match_variance_proportions_exactly(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    variance = result.variance_table
    assert np.all(np.isfinite(variance["explained_variance_ratio"]))
    assert (variance["explained_variance_ratio"] >= 0.0).all()
    assert (variance["explained_variance_ratio"] <= 1.0).all()
    expected_ratio = variance["explained_variance"] / variance["explained_variance"].sum()
    assert variance["explained_variance_ratio"].to_numpy() == pytest.approx(
        expected_ratio.to_numpy()
    )


def test_explained_variance_ratios_sum_to_one_within_tolerance(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    result = compute_sample_pca(expression, metadata)

    ratio_sum = float(result.variance_table["explained_variance_ratio"].sum())
    assert ratio_sum == pytest.approx(1.0, abs=1e-9)


def test_sign_convention_is_deterministic_on_nondegenerate_fixture(
    separated_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = separated_tables

    first = compute_sample_pca(expression, metadata)
    second = compute_sample_pca(expression, metadata)

    pd.testing.assert_frame_equal(first.score_table, second.score_table)
    scores = first.score_table.set_index("sample_id")["pc1"]
    largest_magnitude_sample = scores.abs().idxmax()
    assert scores[largest_magnitude_sample] > 0


def test_tied_singular_value_fixture_preserves_gram_matrix_and_distances() -> None:
    expression, metadata = _tied_subspace_tables()

    result = compute_sample_pca(expression, metadata)

    assert result.component_count == 2
    ratios = result.variance_table["explained_variance_ratio"].to_numpy()
    assert ratios == pytest.approx([0.5, 0.5])

    scores = result.score_table[["pc1", "pc2"]].to_numpy()
    score_gram = scores @ scores.T

    sample_columns = [c for c in expression.columns if c != "gene_id"]
    centered = expression.loc[:, sample_columns].astype(float).to_numpy().T
    centered = centered - centered.mean(axis=0)
    source_gram = centered @ centered.T

    assert score_gram == pytest.approx(source_gram, abs=1e-9)

    score_distances = np.linalg.norm(
        scores[:, None, :] - scores[None, :, :], axis=-1
    )
    source_distances = np.linalg.norm(
        centered[:, None, :] - centered[None, :, :], axis=-1
    )
    assert score_distances == pytest.approx(source_distances, abs=1e-9)


def test_repeated_computation_is_equivalent(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    first = compute_sample_pca(expression, metadata)
    second = compute_sample_pca(expression, metadata)

    assert first.sample_count == second.sample_count
    assert first.gene_count == second.gene_count
    assert first.component_count == second.component_count
    assert first.total_variance == second.total_variance
    pd.testing.assert_frame_equal(first.score_table, second.score_table)
    pd.testing.assert_frame_equal(first.variance_table, second.variance_table)


def test_zero_sample_columns_is_controlled() -> None:
    expression = pd.DataFrame({"gene_id": ["g1", "g2"]})
    metadata = pd.DataFrame(columns=["sample_id", "condition"])

    _assert_reason(
        expression,
        metadata,
        PcaErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
    )


def test_one_sample_is_controlled_insufficient_sample_count() -> None:
    expression = pd.DataFrame({"gene_id": ["g1", "g2"], "sample_a": [1, 2]})
    metadata = pd.DataFrame({"sample_id": ["sample_a"], "condition": ["control"]})
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    error = _assert_reason(
        expression,
        metadata,
        PcaErrorReason.INSUFFICIENT_SAMPLE_COUNT,
    )

    assert "at least two sample columns" in str(error).lower()
    pd.testing.assert_frame_equal(expression, original_expression)
    pd.testing.assert_frame_equal(metadata, original_metadata)


def test_exactly_two_samples_has_one_component() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2", "g3"], "sample_a": [1, 2, 3], "sample_b": [4, 6, 9]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_a", "sample_b"], "condition": ["control", "treated"]}
    )

    result = compute_sample_pca(expression, metadata)

    assert result.component_count == 1
    assert tuple(result.score_table.columns) == (*SCORE_TABLE_FIXED_COLUMNS, "pc1")
    with pytest.raises(ValueError):
        build_score_plot_data(result)


def test_zero_genes_is_controlled() -> None:
    expression = pd.DataFrame(columns=["gene_id", "sample_a", "sample_b"])
    metadata = pd.DataFrame(
        {"sample_id": ["sample_a", "sample_b"], "condition": ["control", "treated"]}
    )

    _assert_reason(expression, metadata, PcaErrorReason.EMPTY_EXPRESSION_TABLE)


def test_one_gene_proceeds_with_observation() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "sample_a": [1], "sample_b": [5], "sample_c": [3]}
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )

    result = compute_sample_pca(expression, metadata)
    observations = build_pca_observations(result, ValidationReport())

    assert result.component_count == 1
    assert result.total_variance > 0
    assert any("only one gene" in item.lower() for item in observations)


def test_constant_genes_retained_with_positive_total_variance_proceeds() -> None:
    expression, metadata = _some_constant_gene_tables()
    original_expression = expression.copy(deep=True)

    result = compute_sample_pca(expression, metadata)
    observations = build_pca_observations(result, ValidationReport())

    assert result.zero_variance_gene_count == 1
    assert result.component_count == 2
    assert result.total_variance > 0
    assert any("1 gene(s)" in item for item in observations)
    pd.testing.assert_frame_equal(expression, original_expression)


def test_all_genes_constant_is_controlled_zero_total_variance() -> None:
    expression, metadata = _all_constant_tables()
    original_expression = expression.copy(deep=True)

    error = _assert_reason(
        expression,
        metadata,
        PcaErrorReason.ZERO_TOTAL_VARIANCE,
    )

    assert "no variation" in str(error).lower()
    assert "invalid" not in str(error).lower()
    assert "low quality" not in str(error).lower()
    pd.testing.assert_frame_equal(expression, original_expression)


def test_large_finite_values_are_numerical_range_error_without_runtime_error() -> None:
    expression, metadata = _large_finite_value_tables()
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        error = _assert_reason(
            expression,
            metadata,
            PcaErrorReason.NUMERICAL_RANGE_ERROR,
        )

    assert "float64" in str(error)
    assert "not modified" in str(error)
    assert "invalid" not in str(error).lower()
    pd.testing.assert_frame_equal(expression, original_expression)
    pd.testing.assert_frame_equal(metadata, original_metadata)


def test_tiny_finite_values_are_numerical_range_error_not_zero_total_variance() -> None:
    expression, metadata = _tiny_finite_value_tables()
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        error = _assert_reason(
            expression,
            metadata,
            PcaErrorReason.NUMERICAL_RANGE_ERROR,
        )

    assert error.reason is not PcaErrorReason.ZERO_TOTAL_VARIANCE
    pd.testing.assert_frame_equal(expression, original_expression)
    pd.testing.assert_frame_equal(metadata, original_metadata)


def test_exactly_constant_ordinary_values_remain_zero_total_variance() -> None:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample_a": [5.0, 10.0],
            "sample_b": [5.0, 10.0],
            "sample_c": [5.0, 10.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )

    _assert_reason(expression, metadata, PcaErrorReason.ZERO_TOTAL_VARIANCE)


def test_normal_finite_data_pca_unaffected_by_numerical_range_guards(
    colinear_symmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = colinear_symmetric_tables

    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        result = compute_sample_pca(expression, metadata)

    scores = result.score_table.set_index("sample_id")["pc1"]
    assert scores["sample_a"] == pytest.approx(-scores["sample_c"], rel=1e-9)
    assert scores["sample_b"] == pytest.approx(0.0, abs=1e-9)
    variance = result.variance_table
    assert variance.loc[
        variance["component"] == 1, "explained_variance_ratio"
    ].iloc[0] == pytest.approx(1.0)


def test_rank_deficient_positive_variance_matrix_proceeds() -> None:
    expression, metadata = _rank_deficient_tables()

    result = compute_sample_pca(expression, metadata)

    assert result.component_count == 2
    ratio_sum = float(result.variance_table["explained_variance_ratio"].sum())
    assert ratio_sum == pytest.approx(1.0, abs=1e-9)


def test_second_component_zero_variance_is_exact_zero_row_not_dropped(
    separated_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = separated_tables

    result = compute_sample_pca(expression, metadata)

    assert len(result.variance_table.index) == 2
    second_row = result.variance_table.loc[result.variance_table["component"] == 2].iloc[0]
    assert second_row["singular_value"] == 0.0
    assert second_row["explained_variance"] == 0.0
    assert second_row["explained_variance_ratio"] == 0.0


def test_missing_gene_id_column_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables

    _assert_reason(
        expression.drop(columns="gene_id"),
        metadata,
        PcaErrorReason.MISSING_REQUIRED_COLUMN,
    )


def test_duplicate_gene_id_column_is_controlled() -> None:
    expression = pd.DataFrame(
        [["g1", "copy1", 1, 2], ["g2", "copy2", 2, 3]],
        columns=["gene_id", "gene_id", "sample_1", "sample_2"],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_1", "sample_2"], "condition": ["control", "treated"]}
    )

    error = _assert_reason(
        expression,
        metadata,
        PcaErrorReason.DUPLICATE_REQUIRED_COLUMN,
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
        PcaErrorReason.MISSING_REQUIRED_IDENTIFIER,
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
        PcaErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
    )

    assert "g1" in str(error)


@pytest.mark.parametrize("bad_identifier", ["", "   ", None])
def test_missing_sample_column_identifier_is_controlled(
    bad_identifier: object,
) -> None:
    expression = pd.DataFrame(
        [["g1", 1, 3], ["g2", 2, 4]],
        columns=["gene_id", bad_identifier, "sample_2"],
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_1", "sample_2"], "condition": ["control", "treated"]}
    )

    _assert_reason(
        expression,
        metadata,
        PcaErrorReason.MISSING_REQUIRED_IDENTIFIER,
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
        PcaErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
    )


@pytest.mark.parametrize("column", ["sample_id", "condition"])
def test_missing_metadata_required_column_is_controlled(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    column: str,
) -> None:
    expression, metadata = asymmetric_tables

    _assert_reason(
        expression,
        metadata.drop(columns=column),
        PcaErrorReason.MISSING_REQUIRED_COLUMN,
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
        PcaErrorReason.DUPLICATE_REQUIRED_COLUMN,
    )


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
        PcaErrorReason.MISSING_REQUIRED_VALUE,
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
        PcaErrorReason.MISSING_REQUIRED_IDENTIFIER,
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
        PcaErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
    )

    assert "sample_a" in str(error)


def test_metadata_expression_mismatch_has_both_directions(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    metadata = metadata.copy(deep=True)
    metadata.loc[metadata["sample_id"].eq("sample_d"), "sample_id"] = "metadata_only"

    error = _assert_reason(
        expression,
        metadata,
        PcaErrorReason.SAMPLE_MISMATCH,
    )

    assert "missing metadata IDs: sample_d" in str(error)
    assert "metadata-only IDs: metadata_only" in str(error)


@pytest.mark.parametrize(
    ("bad_value", "reason", "message_part"),
    [
        (
            "not-numeric",
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "1 non-coercible",
        ),
        (
            None,
            PcaErrorReason.MISSING_EXPRESSION_VALUE,
            "1 missing or blank",
        ),
        (
            "   ",
            PcaErrorReason.MISSING_EXPRESSION_VALUE,
            "1 missing or blank",
        ),
        (
            True,
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "1 boolean and 0 complex",
        ),
        (
            np.bool_(False),
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "1 boolean and 0 complex",
        ),
        (
            1 + 2j,
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "0 boolean and 1 complex",
        ),
        (
            np.complex128(3 + 4j),
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "0 boolean and 1 complex",
        ),
        (
            float("inf"),
            PcaErrorReason.INFINITE_EXPRESSION_VALUE,
            "1 positive and 0 negative",
        ),
        (
            float("-inf"),
            PcaErrorReason.INFINITE_EXPRESSION_VALUE,
            "0 positive and 1 negative",
        ),
    ],
)
def test_invalid_expression_value_has_specific_controlled_reason(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
    bad_value: object,
    reason: PcaErrorReason,
    message_part: str,
) -> None:
    expression, metadata = asymmetric_tables
    expression = expression.copy(deep=True)
    expression["sample_b"] = expression["sample_b"].astype("object")
    expression.loc[expression.index[0], "sample_b"] = bad_value

    error = _assert_reason(expression, metadata, reason)

    assert message_part in str(error)


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

    result = compute_sample_pca(expression, metadata)

    assert result.component_count == 3
    pd.testing.assert_frame_equal(expression, original_expression)
    pd.testing.assert_frame_equal(metadata, original_metadata)
    pd.testing.assert_series_equal(expression.dtypes, expression_dtypes)
    pd.testing.assert_series_equal(metadata.dtypes, metadata_dtypes)


def test_negative_expression_values_are_permitted_and_preserved() -> None:
    # Explicit negative *raw input* values (not merely values that become
    # negative after mean-centering) -- an otherwise valid, non-constant,
    # non-degenerate matrix.
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_a": [-5.0, 2.0, -1.5],
            "sample_b": [3.0, -4.0, 0.5],
            "sample_c": [-2.0, 6.0, 3.5],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    assert (expression[["sample_a", "sample_b", "sample_c"]].to_numpy() < 0).any()
    original_expression = expression.copy(deep=True)
    original_metadata = metadata.copy(deep=True)
    expression_dtypes = expression.dtypes.copy(deep=True)
    metadata_dtypes = metadata.dtypes.copy(deep=True)

    result = compute_sample_pca(expression, metadata)

    assert result.sample_count == 3
    assert result.component_count >= 1

    score_columns = [c for c in result.score_table.columns if c.startswith("pc")]
    scores = result.score_table[score_columns].to_numpy()
    assert np.all(np.isfinite(scores))

    ratios = result.variance_table["explained_variance_ratio"].to_numpy()
    assert np.all(np.isfinite(ratios))
    assert float(ratios.sum()) == pytest.approx(1.0, abs=1e-9)

    pd.testing.assert_frame_equal(expression, original_expression)
    pd.testing.assert_frame_equal(metadata, original_metadata)
    pd.testing.assert_series_equal(expression.dtypes, expression_dtypes)
    pd.testing.assert_series_equal(metadata.dtypes, metadata_dtypes)


def test_multiple_invalid_inputs_follow_fixed_error_priority() -> None:
    valid_metadata = pd.DataFrame(
        {"sample_id": ["sample_1", "sample_2"], "condition": ["control", "treated"]}
    )

    _assert_reason(pd.DataFrame(), [], PcaErrorReason.INVALID_TABLE)
    _assert_reason(pd.DataFrame(), pd.DataFrame(), PcaErrorReason.EMPTY_EXPRESSION_TABLE)
    _assert_reason(
        pd.DataFrame({"sample_1": [1]}),
        valid_metadata,
        PcaErrorReason.MISSING_REQUIRED_COLUMN,
    )
    _assert_reason(
        pd.DataFrame({"gene_id": ["g1"]}),
        valid_metadata,
        PcaErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
    )
    _assert_reason(
        pd.DataFrame({"gene_id": ["g1"], "sample_1": [1]}),
        valid_metadata,
        PcaErrorReason.INSUFFICIENT_SAMPLE_COUNT,
    )

    mismatch_metadata = pd.DataFrame(
        {"sample_id": ["metadata_only", "sample_2"], "condition": ["control", "treated"]}
    )
    expr_blank_gene_and_mismatch = pd.DataFrame(
        {"gene_id": ["", "g2"], "sample_1": [1, 2], "sample_2": [3, 4]}
    )
    _assert_reason(
        expr_blank_gene_and_mismatch,
        mismatch_metadata,
        PcaErrorReason.MISSING_REQUIRED_IDENTIFIER,
    )

    expr_ok = pd.DataFrame({"gene_id": ["g1", "g2"], "sample_1": [1, 2], "sample_2": [3, 4]})
    _assert_reason(expr_ok, mismatch_metadata, PcaErrorReason.SAMPLE_MISMATCH)

    expr_missing_value = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [None, 2], "sample_2": [3, 4]}
    )
    _assert_reason(expr_missing_value, mismatch_metadata, PcaErrorReason.SAMPLE_MISMATCH)

    expr_missing_and_boolean = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [None, True], "sample_2": [3, 4]}
    )
    missing_error = _assert_reason(
        expr_missing_and_boolean,
        valid_metadata,
        PcaErrorReason.MISSING_EXPRESSION_VALUE,
    )
    assert "1 missing or blank" in str(missing_error)

    expr_boolean_and_text = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [True, "not-numeric"], "sample_2": [3, 4]}
    )
    boolean_error = _assert_reason(
        expr_boolean_and_text,
        valid_metadata,
        PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
    )
    assert "1 boolean and 0 complex" in str(boolean_error)

    expr_text_and_inf = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample_1": ["not-numeric", float("inf")],
            "sample_2": [3, 4],
        }
    )
    text_error = _assert_reason(
        expr_text_and_inf,
        valid_metadata,
        PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
    )
    assert "1 non-coercible" in str(text_error)

    expr_inf = pd.DataFrame(
        {"gene_id": ["g1", "g2"], "sample_1": [float("inf"), 2], "sample_2": [3, 4]}
    )
    _assert_reason(expr_inf, valid_metadata, PcaErrorReason.INFINITE_EXPRESSION_VALUE)

    all_constant_expression, all_constant_metadata = _all_constant_tables()
    _assert_reason(
        all_constant_expression,
        all_constant_metadata,
        PcaErrorReason.ZERO_TOTAL_VARIANCE,
    )

    # NUMERICAL_RANGE_ERROR is checked immediately after the exact
    # all-constant (ZERO_TOTAL_VARIANCE) check, and strictly before
    # decomposition is ever attempted: the raw total variance is
    # non-finite/underflowed here despite passing every earlier structural
    # and value check, and no gene is exactly constant.
    large_expression, large_metadata = _large_finite_value_tables()
    _assert_reason(
        large_expression,
        large_metadata,
        PcaErrorReason.NUMERICAL_RANGE_ERROR,
    )
    tiny_expression, tiny_metadata = _tiny_finite_value_tables()
    _assert_reason(
        tiny_expression,
        tiny_metadata,
        PcaErrorReason.NUMERICAL_RANGE_ERROR,
    )


def test_decomposition_failed_wraps_only_linalg_error(
    monkeypatch: pytest.MonkeyPatch,
    separated_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = separated_tables

    def _raise_linalg_error(*args: object, **kwargs: object) -> None:
        raise np.linalg.LinAlgError("forced failure")

    monkeypatch.setattr(np.linalg, "svd", _raise_linalg_error)

    error = _assert_reason(
        expression,
        metadata,
        PcaErrorReason.DECOMPOSITION_FAILED,
    )
    assert "did not converge" in str(error).lower()


def test_non_linalg_exception_from_svd_is_not_wrapped(
    monkeypatch: pytest.MonkeyPatch,
    separated_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = separated_tables

    def _raise_value_error(*args: object, **kwargs: object) -> None:
        raise ValueError("unexpected forced failure")

    monkeypatch.setattr(np.linalg, "svd", _raise_value_error)

    with pytest.raises(ValueError) as raised:
        compute_sample_pca(expression, metadata)
    assert not isinstance(raised.value, PcaComputationError)


def test_observations_are_deterministic_and_threshold_free(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)
    report = ValidationReport()

    observations = build_pca_observations(result, report)
    joined = " ".join(observations).lower()

    assert isinstance(observations, tuple)
    assert observations == build_pca_observations(result, report)
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
        "significant",
        "causation",
    ):
        assert prohibited not in joined


def test_observations_skip_irrelevant_validation_messages(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)
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

    observations = build_pca_observations(result, report)

    assert not any("adjusted-p-value" in observation for observation in observations)


def test_observations_surface_relevant_validation_issue_codes() -> None:
    # A clean, non-degenerate, non-reordered fixture so the *only*
    # observation produced is the one derived from the constructed
    # validation report below (no spontaneous structural observations to
    # disentangle from the validation-report-derived one).
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_a": [1.0, 5.0, 2.0],
            "sample_b": [3.0, 1.0, 8.0],
            "sample_c": [2.0, 9.0, 4.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    result = compute_sample_pca(expression, metadata)
    assert build_pca_observations(result, ValidationReport()) == ()

    negative_value_issue = ValidationIssue(
        code=IssueCode.NEGATIVE_EXPRESSION_VALUE,
        severity=Severity.WARNING,
        table="Expression matrix",
        message="Expression matrix contains negative value(s).",
    )
    single_condition_issue = ValidationIssue(
        code=IssueCode.SINGLE_CONDITION,
        severity=Severity.WARNING,
        table="Sample metadata",
        message="Sample metadata contains only one condition.",
    )
    report = ValidationReport((negative_value_issue, single_condition_issue))
    assert not report.has_errors
    assert all(issue.severity is Severity.WARNING for issue in report.issues)

    observations = build_pca_observations(result, report)

    assert isinstance(observations, tuple)
    assert observations == build_pca_observations(result, report)
    assert len(observations) == 1
    relevant_observation = observations[0]
    assert "2 unique relevant non-blocking message(s)" in relevant_observation
    assert IssueCode.NEGATIVE_EXPRESSION_VALUE.value in relevant_observation
    assert IssueCode.SINGLE_CONDITION.value in relevant_observation

    joined_lower = " ".join(observations).lower()
    for prohibited in (
        "outlier",
        "significant",
        "exclude",
        "exclusion",
        "quality",
        "biological validity",
        "causation",
    ):
        assert prohibited not in joined_lower


def test_observations_report_constant_genes_and_metadata_order() -> None:
    expression, metadata = _some_constant_gene_tables()
    reordered_metadata = metadata.iloc[::-1].reset_index(drop=True)
    result = compute_sample_pca(expression, reordered_metadata)

    observations = build_pca_observations(result, ValidationReport())

    assert any("1 gene(s)" in item for item in observations)
    assert any("metadata sample order differs" in item.lower() for item in observations)


def test_observations_report_zero_variance_second_component(
    separated_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = separated_tables

    result = compute_sample_pca(expression, metadata)
    observations = build_pca_observations(result, ValidationReport())

    assert any("component 2 accounts for 0%" in item.lower() for item in observations)


def test_result_schema_has_no_classification_fields(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)
    schema_text = " ".join(
        [
            *(field.name for field in fields(SamplePcaResult)),
            *result.score_table.columns,
            *result.variance_table.columns,
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
        "significan",
        "cluster",
    ):
        assert prohibited not in schema_text


def test_build_score_plot_data_returns_exact_columns_and_order(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)

    plot_data = build_score_plot_data(result)

    assert list(plot_data.columns) == ["sample_id", "condition", "pc1", "pc2"]
    assert plot_data["sample_id"].tolist() == result.score_table["sample_id"].tolist()
    pd.testing.assert_series_equal(
        plot_data["pc1"], result.score_table["pc1"], check_names=False
    )
    pd.testing.assert_series_equal(
        plot_data["pc2"], result.score_table["pc2"], check_names=False
    )


def test_build_score_plot_data_does_not_mutate_score_table(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)
    original_score_table = result.score_table.copy(deep=True)

    plot_data = build_score_plot_data(result)
    plot_data.loc[0, "pc1"] = 999.0

    pd.testing.assert_frame_equal(result.score_table, original_score_table)


def test_build_score_plot_data_raises_value_error_below_two_components() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2", "g3"], "sample_a": [1, 2, 3], "sample_b": [4, 6, 9]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["sample_a", "sample_b"], "condition": ["control", "treated"]}
    )
    result = compute_sample_pca(expression, metadata)

    with pytest.raises(ValueError, match="requires at least two"):
        build_score_plot_data(result)


def test_list_additional_metadata_columns_excludes_required_columns() -> None:
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a"],
            "condition": ["control"],
            "genotype": ["WT"],
            "batch": ["1"],
        }
    )

    assert list_additional_metadata_columns(metadata) == ("genotype", "batch")


def test_list_additional_metadata_columns_empty_for_two_column_metadata() -> None:
    metadata = pd.DataFrame({"sample_id": ["sample_a"], "condition": ["control"]})

    assert list_additional_metadata_columns(metadata) == ()


def test_list_additional_metadata_columns_returns_empty_for_non_dataframe() -> None:
    assert list_additional_metadata_columns(None) == ()


def test_build_grouped_score_plot_data_matches_condition_plot_by_default(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)

    pd.testing.assert_frame_equal(
        build_grouped_score_plot_data(result, metadata, "condition"),
        build_score_plot_data(result),
    )


def test_build_grouped_score_plot_data_joins_arbitrary_metadata_column_by_sample_id(
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
    result = compute_sample_pca(expression, metadata)

    plot_data = build_grouped_score_plot_data(result, metadata, "genotype")

    assert list(plot_data.columns) == ["sample_id", "genotype", "pc1", "pc2"]
    for sample_id, genotype in zip(
        plot_data["sample_id"], plot_data["genotype"], strict=True
    ):
        assert genotype == expected_genotype_by_sample[sample_id]


def test_build_grouped_score_plot_data_rejects_unknown_column(
    asymmetric_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = asymmetric_tables
    result = compute_sample_pca(expression, metadata)

    with pytest.raises(ValueError, match="does not contain column"):
        build_grouped_score_plot_data(result, metadata, "tissue")


def test_build_grouped_score_plot_data_raises_value_error_below_two_components() -> (
    None
):
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2", "g3"], "sample_a": [1, 2, 3], "sample_b": [4, 6, 9]}
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b"],
            "condition": ["control", "treated"],
            "genotype": ["WT", "WT"],
        }
    )
    result = compute_sample_pca(expression, metadata)

    with pytest.raises(ValueError, match="requires at least two"):
        build_grouped_score_plot_data(result, metadata, "genotype")
