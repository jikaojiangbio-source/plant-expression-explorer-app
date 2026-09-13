"""Pure descriptive sample principal component analysis (PCA).

Axes inside a subspace of equal or numerically indistinguishable singular
values are not uniquely identified: NumPy/BLAS/platform differences may
return a different, equally valid, orthonormal basis within such a subspace.
The sign convention applied here resolves the +/- ambiguity of a single
non-degenerate component; it does not claim to identify individual axes
inside a tied subspace, and no canonical-basis algorithm is applied.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum

import numpy as np
import pandas as pd
from pandas.api.types import is_bool, is_complex

from plant_expression_explorer.qc import count_zero_variance_genes
from plant_expression_explorer.validation import IssueCode, ValidationReport


class PcaErrorReason(StrEnum):
    """Stable reasons for expected PCA input-contract failures."""

    INVALID_TABLE = "INVALID_TABLE"
    EMPTY_EXPRESSION_TABLE = "EMPTY_EXPRESSION_TABLE"
    NO_EXPRESSION_SAMPLE_COLUMNS = "NO_EXPRESSION_SAMPLE_COLUMNS"
    INSUFFICIENT_SAMPLE_COUNT = "INSUFFICIENT_SAMPLE_COUNT"
    MISSING_REQUIRED_COLUMN = "MISSING_REQUIRED_COLUMN"
    DUPLICATE_REQUIRED_COLUMN = "DUPLICATE_REQUIRED_COLUMN"
    MISSING_REQUIRED_IDENTIFIER = "MISSING_REQUIRED_IDENTIFIER"
    DUPLICATE_REQUIRED_IDENTIFIER = "DUPLICATE_REQUIRED_IDENTIFIER"
    MISSING_REQUIRED_VALUE = "MISSING_REQUIRED_VALUE"
    NON_COERCIBLE_EXPRESSION_VALUE = "NON_COERCIBLE_EXPRESSION_VALUE"
    MISSING_EXPRESSION_VALUE = "MISSING_EXPRESSION_VALUE"
    INFINITE_EXPRESSION_VALUE = "INFINITE_EXPRESSION_VALUE"
    SAMPLE_MISMATCH = "SAMPLE_MISMATCH"
    ZERO_TOTAL_VARIANCE = "ZERO_TOTAL_VARIANCE"
    NUMERICAL_RANGE_ERROR = "NUMERICAL_RANGE_ERROR"
    DECOMPOSITION_FAILED = "DECOMPOSITION_FAILED"


class PcaComputationError(ValueError):
    """Expected failure when inputs cannot support principal component analysis."""

    def __init__(self, reason: PcaErrorReason, message: str) -> None:
        """Create an error with a stable reason and user-facing evidence."""

        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class SamplePcaResult:
    """Complete descriptive PCA result for one expression/metadata pair.

    ``frozen=True`` prevents field rebinding. It does not make the newly
    constructed nested pandas DataFrames deeply immutable, so callers should
    treat all result tables as read-only and use copies for display formatting.

    ``total_variance`` is the sum of the reported component explained
    variances (``reported_total_variance``), not the independently
    calculated raw total variance used only to gate :class:`PcaComputationError`
    with reason ``ZERO_TOTAL_VARIANCE``. The two are verified internally to be
    close, but the stored value is always the reported one so that
    ``explained_variance_ratio`` values sum to ``1.0`` by construction.
    """

    sample_count: int
    gene_count: int
    component_count: int
    metadata_order_matches_expression: bool
    total_variance: float
    zero_variance_gene_count: int
    score_table: pd.DataFrame
    variance_table: pd.DataFrame


SCORE_TABLE_FIXED_COLUMNS = ("sample_id", "condition")
VARIANCE_TABLE_COLUMNS = (
    "component",
    "singular_value",
    "explained_variance",
    "explained_variance_ratio",
)

_RATIO_SUM_ABS_TOLERANCE = 1e-9
_TOTAL_VARIANCE_REL_TOLERANCE = 1e-9
_NUMERICAL_RANGE_MESSAGE = (
    "Descriptive PCA could not safely represent the derived centred "
    "variance in float64 for these finite expression values. The source "
    "data were not modified; provide a preprocessed matrix on a "
    "numerically representable scale."
)
_RELEVANT_VALIDATION_OBSERVATIONS = frozenset(
    {
        IssueCode.COERCIBLE_NUMERIC_STRING,
        IssueCode.NEGATIVE_EXPRESSION_VALUE,
        IssueCode.SAMPLE_ORDER_DIFFERS,
        IssueCode.SINGLE_CONDITION,
        IssueCode.CONDITION_WITH_SINGLE_SAMPLE,
    }
)


def compute_sample_pca(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
) -> SamplePcaResult:
    """Calculate non-mutating descriptive PCA sample scores and variance.

    Samples are observations and genes are features. Each gene is
    mean-centred across samples in a temporary numeric copy; genes are not
    scaled to unit variance. Conditions are mapped by exact sample ID and are
    never used to fit the components. Safely coercible numeric strings are
    converted only in the temporary copy. Missing, blank, boolean, complex,
    non-coercible, and infinite values cause a controlled error; no gene or
    sample is filtered, trimmed, imputed, or reordered.

    Expected input failures have a fixed priority: table types; expression
    row, column, and identifier structure; metadata columns, identifiers,
    and conditions; cross-table sample matching; then expression values in
    missing/blank, boolean/complex, non-coercible, and infinite order; then
    the zero-total-variance check; then numerical-range checks (finite
    source values that cannot be safely represented through the float64
    calculation, distinguished from zero variance using the already
    computed exact constant-gene count); then decomposition failure.
    """

    _require_dataframe(expression, "Expression matrix")
    _require_dataframe(metadata, "Sample metadata")
    if len(expression.index) == 0:
        raise PcaComputationError(
            PcaErrorReason.EMPTY_EXPRESSION_TABLE,
            "The expression matrix contains no gene rows.",
        )

    _require_one_column(expression, "gene_id", "Expression matrix")

    sample_columns = [
        column
        for column in expression.columns
        if not _column_name_matches(column, "gene_id")
    ]
    if not sample_columns:
        raise PcaComputationError(
            PcaErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
            "The expression matrix contains no sample columns besides "
            "'gene_id'.",
        )
    if len(sample_columns) == 1:
        raise PcaComputationError(
            PcaErrorReason.INSUFFICIENT_SAMPLE_COUNT,
            "At least two sample columns are required to estimate variance "
            "for principal component analysis; this dataset provides only "
            "one.",
        )

    gene_ids = _series_identifiers(
        expression["gene_id"],
        label="Expression gene IDs",
    )
    sample_ids = _column_identifiers(sample_columns)

    _require_one_column(metadata, "sample_id", "Sample metadata")
    _require_one_column(metadata, "condition", "Sample metadata")
    metadata_sample_ids = _series_identifiers(
        metadata["sample_id"],
        label="Metadata sample IDs",
    )
    conditions = _required_values(
        metadata["condition"],
        label="Metadata conditions",
    )
    _require_matching_samples(sample_ids, metadata_sample_ids)

    numeric_expression = _numeric_expression_copy(
        expression,
        sample_columns,
        sample_ids,
    )
    condition_by_sample = dict(
        zip(metadata_sample_ids, conditions, strict=True)
    )
    sample_conditions = [
        condition_by_sample[sample_id] for sample_id in sample_ids
    ]

    gene_count = len(gene_ids)
    sample_count = len(sample_ids)
    zero_variance_gene_count = count_zero_variance_genes(numeric_expression)

    with np.errstate(over="ignore", invalid="ignore"):
        raw_total_variance = float(
            numeric_expression.var(axis=1, ddof=1).sum()
        )
    _require_finite_scalar(raw_total_variance)
    if raw_total_variance == 0.0:
        if zero_variance_gene_count == gene_count:
            raise PcaComputationError(
                PcaErrorReason.ZERO_TOTAL_VARIANCE,
                "Principal component analysis cannot produce informative "
                "variance components because the temporary centred matrix "
                f"has no variation: every gene is constant across the "
                f"{sample_count} sample(s).",
            )
        # At least one gene is not exactly constant (per the exact-equality
        # check above), yet the ddof=1 variance summed to exactly zero: the
        # true variance underflowed to zero in float64. This is a numerical
        # range failure, not a biologically constant matrix.
        raise PcaComputationError(
            PcaErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        )

    working_matrix = numeric_expression.to_numpy(dtype=float).T
    _require_finite_array(working_matrix)
    with np.errstate(over="ignore", invalid="ignore"):
        gene_means = working_matrix.mean(axis=0)
        centered_matrix = working_matrix - gene_means
    _require_finite_array(gene_means)
    _require_finite_array(centered_matrix)

    try:
        left_vectors, singular_values, _right_vectors = np.linalg.svd(
            centered_matrix,
            full_matrices=False,
        )
    except np.linalg.LinAlgError as exc:
        raise PcaComputationError(
            PcaErrorReason.DECOMPOSITION_FAILED,
            "Singular value decomposition did not converge for the centred "
            "expression matrix.",
        ) from exc
    _require_finite_array(singular_values)

    component_count = min(sample_count - 1, gene_count)
    left_vectors_k = left_vectors[:, :component_count]
    singular_values_k = singular_values[:component_count].copy()

    with np.errstate(over="ignore", invalid="ignore"):
        # Multiply the (small) dimension/epsilon factor first to reduce the
        # risk of unnecessary intermediate overflow when the largest
        # singular value is itself large but still finite.
        dimension_scale = max(sample_count, gene_count) * np.finfo(float).eps
        tolerance = float(singular_values.max()) * dimension_scale
    _require_finite_scalar(tolerance)
    singular_values_snapped = np.where(
        singular_values_k <= tolerance,
        0.0,
        singular_values_k,
    )

    with np.errstate(over="ignore", invalid="ignore"):
        scores = left_vectors_k * singular_values_snapped
    _require_finite_array(scores)
    for component_index in range(component_count):
        if singular_values_snapped[component_index] == 0.0:
            continue
        column = scores[:, component_index]
        pivot = int(np.argmax(np.abs(column)))
        if column[pivot] < 0:
            scores[:, component_index] = -column

    with np.errstate(over="ignore", invalid="ignore"):
        explained_variance = (singular_values_snapped**2) / (sample_count - 1)
    _require_finite_array(explained_variance)

    with np.errstate(over="ignore", invalid="ignore"):
        reported_total_variance = float(explained_variance.sum())
    _require_finite_scalar(reported_total_variance)
    if reported_total_variance <= 0.0:
        raise PcaComputationError(
            PcaErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        )

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        explained_variance_ratio = explained_variance / reported_total_variance
    _require_finite_array(explained_variance_ratio)
    ratio_sum = float(explained_variance_ratio.sum())

    _verify_variance_totals(
        reported_total_variance,
        raw_total_variance,
        ratio_sum,
    )
    _verify_ratio_bounds(explained_variance_ratio)

    score_columns: dict[str, object] = {
        "sample_id": sample_ids,
        "condition": sample_conditions,
    }
    for component_index in range(component_count):
        score_columns[f"pc{component_index + 1}"] = scores[:, component_index]
    score_table = pd.DataFrame(
        score_columns,
        columns=(
            *SCORE_TABLE_FIXED_COLUMNS,
            *(f"pc{i}" for i in range(1, component_count + 1)),
        ),
    )

    variance_table = pd.DataFrame(
        {
            "component": list(range(1, component_count + 1)),
            "singular_value": singular_values_snapped,
            "explained_variance": explained_variance,
            "explained_variance_ratio": explained_variance_ratio,
        },
        columns=VARIANCE_TABLE_COLUMNS,
    )

    return SamplePcaResult(
        sample_count=sample_count,
        gene_count=gene_count,
        component_count=component_count,
        metadata_order_matches_expression=(
            metadata_sample_ids == sample_ids
        ),
        total_variance=reported_total_variance,
        zero_variance_gene_count=zero_variance_gene_count,
        score_table=score_table,
        variance_table=variance_table,
    )


def build_pca_observations(
    result: SamplePcaResult,
    validation_report: ValidationReport,
) -> tuple[str, ...]:
    """Build exact observations without thresholds or quality classification."""

    observations: list[str] = []
    if result.zero_variance_gene_count:
        observations.append(
            f"{result.zero_variance_gene_count} gene(s) have exactly "
            "identical values across all samples and were retained. They "
            "contribute no variance to the components."
        )
    if not result.metadata_order_matches_expression:
        observations.append(
            "Metadata sample order differs from expression-column order. "
            "Conditions were mapped by exact sample ID without modifying "
            "either table."
        )
    if result.sample_count == 2:
        observations.append(
            "Exactly two samples are present, so only one principal "
            "component exists; a PC1-versus-PC2 plot is not available."
        )
    if result.gene_count == 1:
        observations.append(
            "Only one gene is present; the single principal component "
            "reproduces its centred values and no dimensionality reduction "
            "occurs."
        )
    if result.component_count >= 2:
        second_component = result.variance_table.loc[
            result.variance_table["component"] == 2,
            "explained_variance_ratio",
        ]
        if not second_component.empty and float(second_component.iloc[0]) == 0.0:
            observations.append(
                "Component 2 accounts for 0% of total variance; its "
                "position in the plot does not reflect additional "
                "structure."
            )

    relevant_issues = []
    seen_messages: set[str] = set()
    for issue in (*validation_report.warnings, *validation_report.information):
        if issue.code not in _RELEVANT_VALIDATION_OBSERVATIONS:
            continue
        if issue.message in seen_messages:
            continue
        seen_messages.add(issue.message)
        relevant_issues.append(issue)
    if relevant_issues:
        issue_codes = tuple(
            dict.fromkeys(issue.code.value for issue in relevant_issues)
        )
        observations.append(
            f"The validation report contains {len(relevant_issues)} unique "
            "relevant non-blocking message(s) with code(s): "
            + ", ".join(issue_codes)
            + "."
        )
    return tuple(observations)


def build_score_plot_data(result: SamplePcaResult) -> pd.DataFrame:
    """Return a PC1-versus-PC2 plotting copy; requires two or more components.

    This is a programmer-facing precondition, not a data-validation error:
    callers must check ``result.component_count >= 2`` before calling.
    Raises ``ValueError`` (not :class:`PcaComputationError`) otherwise.
    """

    if result.component_count < 2:
        raise ValueError(
            "A PC1-versus-PC2 plot requires at least two principal "
            f"components; this result has {result.component_count}."
        )
    return result.score_table.loc[
        :,
        ["sample_id", "condition", "pc1", "pc2"],
    ].copy(deep=True)


def build_grouped_score_plot_data(
    result: SamplePcaResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return a PC1-versus-PC2 plotting copy coloured by one metadata column.

    ``group_column`` values are joined from ``metadata`` by exact sample ID,
    for display only, exactly like the ``condition`` column already carried
    on ``result.score_table``; they are never used to fit the components.
    Requires two or more components; see :func:`build_score_plot_data`.
    """

    if group_column == "condition":
        return build_score_plot_data(result)
    if result.component_count < 2:
        raise ValueError(
            "A PC1-versus-PC2 plot requires at least two principal "
            f"components; this result has {result.component_count}."
        )
    if not isinstance(metadata, pd.DataFrame) or "sample_id" not in metadata.columns:
        raise ValueError("Sample metadata is missing required column 'sample_id'.")
    if group_column not in metadata.columns:
        raise ValueError(f"Sample metadata does not contain column '{group_column}'.")

    lookup = {
        str(sample_id): value
        for sample_id, value in zip(
            metadata["sample_id"], metadata[group_column], strict=True
        )
    }
    plot_data = result.score_table.loc[:, ["sample_id", "pc1", "pc2"]].copy(deep=True)
    plot_data[group_column] = plot_data["sample_id"].map(lookup)
    return plot_data.loc[:, ["sample_id", group_column, "pc1", "pc2"]]


def _require_finite_scalar(value: float) -> None:
    if not math.isfinite(value):
        raise PcaComputationError(
            PcaErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        )


def _require_finite_array(values: np.ndarray) -> None:
    if not np.all(np.isfinite(values)):
        raise PcaComputationError(
            PcaErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        )


def _verify_variance_totals(
    reported_total_variance: float,
    raw_total_variance: float,
    ratio_sum: float,
) -> None:
    if not math.isclose(
        reported_total_variance,
        raw_total_variance,
        rel_tol=_TOTAL_VARIANCE_REL_TOLERANCE,
    ):
        raise RuntimeError(
            "Reported total variance does not closely match the "
            "independently calculated raw total variance."
        )
    if not math.isclose(
        ratio_sum,
        1.0,
        abs_tol=_RATIO_SUM_ABS_TOLERANCE,
    ):
        raise RuntimeError(
            "Explained-variance ratios do not sum to 1.0 within tolerance."
        )


def _verify_ratio_bounds(explained_variance_ratio: np.ndarray) -> None:
    if not np.all(np.isfinite(explained_variance_ratio)):
        raise RuntimeError("Explained-variance ratios must be finite.")
    if np.any(explained_variance_ratio < 0.0) or np.any(
        explained_variance_ratio > 1.0 + _RATIO_SUM_ABS_TOLERANCE
    ):
        raise RuntimeError(
            "Explained-variance ratios must lie within [0, 1]."
        )


def _require_dataframe(table: object, label: str) -> None:
    if not isinstance(table, pd.DataFrame):
        raise PcaComputationError(
            PcaErrorReason.INVALID_TABLE,
            f"{label} must be a pandas DataFrame.",
        )


def _require_one_column(
    table: pd.DataFrame,
    column: str,
    label: str,
) -> None:
    count = sum(
        _column_name_matches(value, column) for value in table.columns
    )
    if count == 0:
        raise PcaComputationError(
            PcaErrorReason.MISSING_REQUIRED_COLUMN,
            f"{label} is missing required column '{column}'.",
        )
    if count > 1:
        raise PcaComputationError(
            PcaErrorReason.DUPLICATE_REQUIRED_COLUMN,
            f"{label} contains {count} columns named '{column}'.",
        )


def _column_name_matches(value: object, expected: str) -> bool:
    return isinstance(value, str) and value == expected


def _column_identifiers(columns: list[object]) -> list[str]:
    if any(_is_missing_or_blank(value) for value in columns):
        raise PcaComputationError(
            PcaErrorReason.MISSING_REQUIRED_IDENTIFIER,
            "Expression sample-column identifiers must be non-missing and "
            "non-blank.",
        )
    identifiers = [str(value) for value in columns]
    duplicates = _duplicates(identifiers)
    if duplicates:
        raise PcaComputationError(
            PcaErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
            "Expression sample-column identifiers must be unique; duplicate "
            "identifier(s): "
            + ", ".join(duplicates)
            + ".",
        )
    return identifiers


def _series_identifiers(series: pd.Series, *, label: str) -> list[str]:
    values = series.tolist()
    missing_positions = [
        position
        for position, value in enumerate(values)
        if _is_missing_or_blank(value)
    ]
    if missing_positions:
        raise PcaComputationError(
            PcaErrorReason.MISSING_REQUIRED_IDENTIFIER,
            f"{label} contain {len(missing_positions)} missing or blank "
            "value(s) at zero-based row position(s): "
            + ", ".join(str(position) for position in missing_positions[:5])
            + ".",
        )
    identifiers = [str(value) for value in values]
    duplicates = _duplicates(identifiers)
    if duplicates:
        raise PcaComputationError(
            PcaErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
            f"{label} must be unique; duplicate identifier(s): "
            + ", ".join(duplicates)
            + ".",
        )
    return identifiers


def _required_values(series: pd.Series, *, label: str) -> list[str]:
    values = series.tolist()
    missing_positions = [
        position
        for position, value in enumerate(values)
        if _is_missing_or_blank(value)
    ]
    if missing_positions:
        raise PcaComputationError(
            PcaErrorReason.MISSING_REQUIRED_VALUE,
            f"{label} contain {len(missing_positions)} missing or blank "
            "value(s) at zero-based row position(s): "
            + ", ".join(str(position) for position in missing_positions[:5])
            + ".",
        )
    return [str(value) for value in values]


def _require_matching_samples(
    sample_ids: list[str],
    metadata_sample_ids: list[str],
) -> None:
    expression_set = set(sample_ids)
    metadata_set = set(metadata_sample_ids)
    if expression_set == metadata_set:
        return

    missing = [
        sample_id for sample_id in sample_ids if sample_id not in metadata_set
    ]
    unexpected = [
        sample_id
        for sample_id in metadata_sample_ids
        if sample_id not in expression_set
    ]
    details: list[str] = []
    if missing:
        details.append("missing metadata IDs: " + ", ".join(missing))
    if unexpected:
        details.append("metadata-only IDs: " + ", ".join(unexpected))
    raise PcaComputationError(
        PcaErrorReason.SAMPLE_MISMATCH,
        "Expression and metadata sample IDs differ ("
        + "; ".join(details)
        + ").",
    )


def _numeric_expression_copy(
    expression: pd.DataFrame,
    sample_columns: list[object],
    sample_ids: list[str],
) -> pd.DataFrame:
    working = expression.loc[:, sample_columns].copy(deep=True)
    missing_mask = working.apply(
        lambda column: column.map(_is_missing_or_blank)
    )
    missing_count = int(missing_mask.sum().sum())
    if missing_count:
        raise PcaComputationError(
            PcaErrorReason.MISSING_EXPRESSION_VALUE,
            f"The expression matrix contains {missing_count} missing or blank "
            "sample value(s).",
        )

    boolean_mask = working.apply(lambda column: column.map(is_bool))
    complex_mask = working.apply(lambda column: column.map(is_complex))
    boolean_count = int(boolean_mask.sum().sum())
    complex_count = int(complex_mask.sum().sum())
    if boolean_count or complex_count:
        raise PcaComputationError(
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "The expression matrix contains "
            f"{boolean_count} boolean and {complex_count} complex sample "
            "value(s); expression values must be real-valued numbers or "
            "safely coercible numeric strings.",
        )

    numeric = working.apply(
        lambda column: pd.to_numeric(column, errors="coerce")
    )
    non_coercible_count = int(numeric.isna().sum().sum())
    if non_coercible_count:
        raise PcaComputationError(
            PcaErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            f"The expression matrix contains {non_coercible_count} "
            "non-coercible sample value(s).",
        )

    positive_infinity_count = int(numeric.eq(float("inf")).sum().sum())
    negative_infinity_count = int(numeric.eq(float("-inf")).sum().sum())
    if positive_infinity_count or negative_infinity_count:
        raise PcaComputationError(
            PcaErrorReason.INFINITE_EXPRESSION_VALUE,
            "The expression matrix contains "
            f"{positive_infinity_count} positive and "
            f"{negative_infinity_count} negative infinite sample value(s).",
        )

    numeric.columns = sample_ids
    return numeric


def _duplicates(values: list[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: list[str] = []
    for value in values:
        if value in seen and value not in duplicates:
            duplicates.append(value)
        seen.add(value)
    return duplicates


def _is_missing_or_blank(value: object) -> bool:
    if isinstance(value, str):
        return not value.strip()
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
