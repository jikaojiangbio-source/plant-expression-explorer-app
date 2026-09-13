"""Pure descriptive sample-to-sample Pearson correlation calculations."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from typing import Literal

import pandas as pd
from pandas.api.types import is_bool, is_complex

from plant_expression_explorer.qc import find_constant_samples
from plant_expression_explorer.validation import IssueCode, ValidationReport


class CorrelationErrorReason(StrEnum):
    """Stable reasons for expected correlation input-contract failures."""

    INVALID_TABLE = "INVALID_TABLE"
    EMPTY_EXPRESSION_TABLE = "EMPTY_EXPRESSION_TABLE"
    INSUFFICIENT_GENE_ROWS = "INSUFFICIENT_GENE_ROWS"
    NO_EXPRESSION_SAMPLE_COLUMNS = "NO_EXPRESSION_SAMPLE_COLUMNS"
    MISSING_REQUIRED_COLUMN = "MISSING_REQUIRED_COLUMN"
    DUPLICATE_REQUIRED_COLUMN = "DUPLICATE_REQUIRED_COLUMN"
    MISSING_REQUIRED_IDENTIFIER = "MISSING_REQUIRED_IDENTIFIER"
    DUPLICATE_REQUIRED_IDENTIFIER = "DUPLICATE_REQUIRED_IDENTIFIER"
    MISSING_REQUIRED_VALUE = "MISSING_REQUIRED_VALUE"
    NON_COERCIBLE_EXPRESSION_VALUE = "NON_COERCIBLE_EXPRESSION_VALUE"
    MISSING_EXPRESSION_VALUE = "MISSING_EXPRESSION_VALUE"
    INFINITE_EXPRESSION_VALUE = "INFINITE_EXPRESSION_VALUE"
    SAMPLE_MISMATCH = "SAMPLE_MISMATCH"


class CorrelationComputationError(ValueError):
    """Expected failure when inputs cannot support Pearson correlation."""

    def __init__(
        self,
        reason: CorrelationErrorReason,
        message: str,
    ) -> None:
        """Create an error with a stable reason and user-facing evidence."""

        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class SampleCorrelationResult:
    """Complete descriptive Pearson result for one expression/metadata pair.

    ``frozen=True`` prevents field rebinding. It does not make the newly
    constructed nested pandas DataFrames deeply immutable, so callers should
    treat all result tables as read-only and use copies for display formatting.
    """

    method: Literal["pearson"]
    gene_count: int
    sample_count: int
    metadata_order_matches_expression: bool
    correlation_matrix: pd.DataFrame
    pair_summary: pd.DataFrame
    sample_summary: pd.DataFrame
    condition_summary: pd.DataFrame
    constant_samples: tuple[str, ...]
    undefined_pair_count: int


PAIR_SUMMARY_COLUMNS = (
    "sample_a",
    "condition_a",
    "sample_b",
    "condition_b",
    "relationship",
    "correlation",
    "defined",
)
SAMPLE_CORRELATION_SUMMARY_COLUMNS = (
    "sample_id",
    "condition",
    "defined_pair_count",
    "undefined_pair_count",
    "minimum_correlation",
    "median_correlation",
    "mean_correlation",
    "maximum_correlation",
)
CONDITION_CORRELATION_SUMMARY_COLUMNS = (
    "condition_a",
    "sample_count_a",
    "condition_b",
    "sample_count_b",
    "relationship",
    "total_pair_count",
    "defined_pair_count",
    "undefined_pair_count",
    "median_correlation",
)
HEATMAP_DATA_COLUMNS = (
    "row_sample",
    "column_sample",
    "correlation",
    "defined",
    "correlation_label",
)

_MATRIX_TOLERANCE = 1e-12
_RELEVANT_VALIDATION_OBSERVATIONS = frozenset(
    {
        IssueCode.COERCIBLE_NUMERIC_STRING,
        IssueCode.NEGATIVE_EXPRESSION_VALUE,
        IssueCode.SAMPLE_ORDER_DIFFERS,
        IssueCode.SINGLE_CONDITION,
        IssueCode.CONDITION_WITH_SINGLE_SAMPLE,
    }
)


def compute_sample_correlation(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
) -> SampleCorrelationResult:
    """Calculate non-mutating descriptive Pearson sample correlations.

    Safely coercible numeric strings are converted only in a temporary working
    copy. Missing, boolean, complex, non-coercible, and infinite values cause a
    controlled error; no invalid cells, genes, or samples are skipped, imputed,
    or removed.

    Expected input failures have a fixed priority: table types; expression row,
    column, and identifier structure; metadata columns, identifiers, and
    conditions; cross-table sample matching; then expression values in the
    order missing/blank, boolean/complex, non-coercible, and infinite.
    """

    _require_dataframe(expression, "Expression matrix")
    _require_dataframe(metadata, "Sample metadata")
    if len(expression.index) == 0:
        raise CorrelationComputationError(
            CorrelationErrorReason.EMPTY_EXPRESSION_TABLE,
            "The expression matrix contains no gene rows.",
        )

    _require_one_column(expression, "gene_id", "Expression matrix")
    if len(expression.index) < 2:
        raise CorrelationComputationError(
            CorrelationErrorReason.INSUFFICIENT_GENE_ROWS,
            "The expression matrix contains 1 gene row; Pearson correlation "
            "requires at least 2 gene rows.",
        )

    sample_columns = [
        column
        for column in expression.columns
        if not _column_name_matches(column, "gene_id")
    ]
    if not sample_columns:
        raise CorrelationComputationError(
            CorrelationErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
            "The expression matrix contains no sample columns besides "
            "'gene_id'.",
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
    constant_samples = find_constant_samples(numeric_expression)
    correlation_matrix = _pearson_matrix(
        numeric_expression,
        sample_ids,
        constant_samples,
    )
    _verify_correlation_matrix(
        correlation_matrix,
        sample_ids,
        constant_samples,
    )

    pair_summary = build_pair_summary(
        correlation_matrix,
        condition_by_sample,
    )
    sample_summary = build_sample_correlation_summary(
        correlation_matrix,
        condition_by_sample,
    )
    condition_summary = build_condition_correlation_summary(
        pair_summary,
        sample_ids,
        sample_conditions,
    )
    undefined_pair_count = int((~pair_summary["defined"]).sum())

    return SampleCorrelationResult(
        method="pearson",
        gene_count=len(gene_ids),
        sample_count=len(sample_ids),
        metadata_order_matches_expression=(
            metadata_sample_ids == sample_ids
        ),
        correlation_matrix=correlation_matrix,
        pair_summary=pair_summary,
        sample_summary=sample_summary,
        condition_summary=condition_summary,
        constant_samples=constant_samples,
        undefined_pair_count=undefined_pair_count,
    )


def build_pair_summary(
    correlation_matrix: pd.DataFrame,
    condition_by_sample: Mapping[str, str],
) -> pd.DataFrame:
    """Return each unique unordered sample pair in matrix-axis order."""

    sample_ids = [str(sample_id) for sample_id in correlation_matrix.columns]
    records: list[dict[str, object]] = []
    for sample_a_position, sample_a in enumerate(sample_ids):
        for sample_b_position in range(
            sample_a_position + 1,
            len(sample_ids),
        ):
            sample_b = sample_ids[sample_b_position]
            condition_a = condition_by_sample[sample_a]
            condition_b = condition_by_sample[sample_b]
            correlation = correlation_matrix.iat[
                sample_a_position,
                sample_b_position,
            ]
            defined = not pd.isna(correlation)
            records.append(
                {
                    "sample_a": sample_a,
                    "condition_a": condition_a,
                    "sample_b": sample_b,
                    "condition_b": condition_b,
                    "relationship": (
                        "within_condition"
                        if condition_a == condition_b
                        else "between_condition"
                    ),
                    "correlation": (
                        float(correlation) if defined else float("nan")
                    ),
                    "defined": bool(defined),
                }
            )
    return pd.DataFrame.from_records(records, columns=PAIR_SUMMARY_COLUMNS)


def build_sample_correlation_summary(
    correlation_matrix: pd.DataFrame,
    condition_by_sample: Mapping[str, str],
) -> pd.DataFrame:
    """Return ordered non-self descriptive correlation metrics per sample."""

    sample_ids = [str(sample_id) for sample_id in correlation_matrix.columns]
    records: list[dict[str, object]] = []
    for sample_position, sample_id in enumerate(sample_ids):
        non_self = correlation_matrix.iloc[sample_position].drop(
            correlation_matrix.columns[sample_position]
        )
        defined_values = non_self[non_self.notna()]
        defined_pair_count = len(defined_values.index)
        undefined_pair_count = len(non_self.index) - defined_pair_count
        records.append(
            {
                "sample_id": sample_id,
                "condition": condition_by_sample[sample_id],
                "defined_pair_count": defined_pair_count,
                "undefined_pair_count": undefined_pair_count,
                "minimum_correlation": _summary_value(
                    defined_values,
                    "min",
                ),
                "median_correlation": _summary_value(
                    defined_values,
                    "median",
                ),
                "mean_correlation": _summary_value(
                    defined_values,
                    "mean",
                ),
                "maximum_correlation": _summary_value(
                    defined_values,
                    "max",
                ),
            }
        )
    return pd.DataFrame.from_records(
        records,
        columns=SAMPLE_CORRELATION_SUMMARY_COLUMNS,
    )


def build_condition_correlation_summary(
    pair_summary: pd.DataFrame,
    expression_sample_ids: list[str],
    conditions: list[str],
) -> pd.DataFrame:
    """Summarize condition pairs in expression-sample first-appearance order."""

    samples_by_condition: dict[str, list[str]] = {}
    for sample_id, condition in zip(
        expression_sample_ids,
        conditions,
        strict=True,
    ):
        samples_by_condition.setdefault(condition, []).append(sample_id)

    condition_order = list(samples_by_condition)
    records: list[dict[str, object]] = []
    for condition_a_position, condition_a in enumerate(condition_order):
        for condition_b_position in range(
            condition_a_position,
            len(condition_order),
        ):
            condition_b = condition_order[condition_b_position]
            sample_count_a = len(samples_by_condition[condition_a])
            sample_count_b = len(samples_by_condition[condition_b])
            same_condition = condition_a == condition_b
            expected_total = (
                sample_count_a * (sample_count_a - 1) // 2
                if same_condition
                else sample_count_a * sample_count_b
            )
            matching_pairs = _condition_pairs(
                pair_summary,
                condition_a,
                condition_b,
            )
            if len(matching_pairs.index) != expected_total:
                raise RuntimeError(
                    "Condition-pair rows do not match the expected sample-pair "
                    f"count for {condition_a!r} and {condition_b!r}."
                )
            defined_pairs = matching_pairs.loc[matching_pairs["defined"]]
            defined_pair_count = len(defined_pairs.index)
            undefined_pair_count = expected_total - defined_pair_count
            median_correlation = (
                float(defined_pairs["correlation"].median())
                if defined_pair_count
                else float("nan")
            )
            records.append(
                {
                    "condition_a": condition_a,
                    "sample_count_a": sample_count_a,
                    "condition_b": condition_b,
                    "sample_count_b": sample_count_b,
                    "relationship": (
                        "within_condition"
                        if same_condition
                        else "between_condition"
                    ),
                    "total_pair_count": expected_total,
                    "defined_pair_count": defined_pair_count,
                    "undefined_pair_count": undefined_pair_count,
                    "median_correlation": median_correlation,
                }
            )
    return pd.DataFrame.from_records(
        records,
        columns=CONDITION_CORRELATION_SUMMARY_COLUMNS,
    )


def build_correlation_observations(
    result: SampleCorrelationResult,
    validation_report: ValidationReport,
) -> tuple[str, ...]:
    """Build exact observations without thresholds or quality classification."""

    observations: list[str] = []
    if result.constant_samples:
        observations.append(
            f"Constant sample column(s) with undefined Pearson values: "
            + ", ".join(result.constant_samples)
            + "."
        )
    if result.undefined_pair_count:
        observations.append(
            f"{result.undefined_pair_count} unique non-self sample pair(s) "
            "have undefined Pearson correlation."
        )
    if result.sample_count == 1:
        observations.append(
            "Only one sample is present, so pairwise sample similarity cannot "
            "be assessed."
        )
    if not result.metadata_order_matches_expression:
        observations.append(
            "Metadata sample order differs from expression-column order. "
            "Conditions were mapped by exact sample ID without modifying "
            "either table."
        )

    defined_correlations = result.pair_summary.loc[
        result.pair_summary["defined"],
        "correlation",
    ]
    if len(defined_correlations.index):
        if defined_correlations.gt(0).all():
            observations.append(
                "All defined unique-pair Pearson correlations are positive."
            )
        elif defined_correlations.lt(0).all():
            observations.append(
                "All defined unique-pair Pearson correlations are negative."
            )
        elif (
            defined_correlations.gt(0).any()
            and defined_correlations.lt(0).any()
        ):
            observations.append(
                "Defined unique-pair Pearson correlations include both "
                "positive and negative values."
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


def build_grouped_pair_summary(
    result: SampleCorrelationResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return each unique sample pair labelled by one metadata column.

    ``group_column`` values are joined from ``metadata`` by exact sample ID,
    for display only, exactly like 'condition'. Reuses
    :func:`build_pair_summary` on the already-computed correlation matrix; no
    Pearson correlation is recalculated for the new grouping.
    """

    if group_column == "condition":
        return result.pair_summary.copy(deep=True)
    lookup = _metadata_column_lookup(metadata, group_column)
    summary = build_pair_summary(result.correlation_matrix, lookup)
    return summary.rename(
        columns={"condition_a": f"{group_column}_a", "condition_b": f"{group_column}_b"}
    )


def build_grouped_sample_correlation_summary(
    result: SampleCorrelationResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return per-sample correlation metrics labelled by one metadata column.

    Reuses :func:`build_sample_correlation_summary` on the already-computed
    correlation matrix; no Pearson correlation is recalculated.
    """

    if group_column == "condition":
        return result.sample_summary.copy(deep=True)
    lookup = _metadata_column_lookup(metadata, group_column)
    summary = build_sample_correlation_summary(result.correlation_matrix, lookup)
    return summary.rename(columns={"condition": group_column})


def build_grouped_condition_correlation_summary(
    result: SampleCorrelationResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return group-pair correlation summaries for one metadata column.

    Reuses :func:`build_condition_correlation_summary` on the
    already-computed pair summary for that grouping; no Pearson correlation
    is recalculated.
    """

    if group_column == "condition":
        return result.condition_summary.copy(deep=True)
    lookup = _metadata_column_lookup(metadata, group_column)
    sample_ids = list(result.correlation_matrix.columns)
    groups = [lookup[sample_id] for sample_id in sample_ids]
    grouped_pairs = build_pair_summary(result.correlation_matrix, lookup)
    summary = build_condition_correlation_summary(grouped_pairs, sample_ids, groups)
    return summary.rename(
        columns={"condition_a": f"{group_column}_a", "condition_b": f"{group_column}_b"}
    )


def _metadata_column_lookup(
    metadata: pd.DataFrame,
    group_column: str,
) -> dict[str, object]:
    if not isinstance(metadata, pd.DataFrame) or "sample_id" not in metadata.columns:
        raise ValueError("Sample metadata is missing required column 'sample_id'.")
    if group_column not in metadata.columns:
        raise ValueError(f"Sample metadata does not contain column '{group_column}'.")
    return {
        str(sample_id): ("(missing)" if _is_missing_or_blank(value) else value)
        for sample_id, value in zip(
            metadata["sample_id"], metadata[group_column], strict=True
        )
    }


def build_heatmap_data(
    result: SampleCorrelationResult,
) -> pd.DataFrame:
    """Return every matrix cell in deterministic row-major display order."""

    records: list[dict[str, object]] = []
    for row_position, row_sample in enumerate(
        result.correlation_matrix.index
    ):
        for column_position, column_sample in enumerate(
            result.correlation_matrix.columns
        ):
            correlation = result.correlation_matrix.iat[
                row_position,
                column_position,
            ]
            defined = not pd.isna(correlation)
            records.append(
                {
                    "row_sample": str(row_sample),
                    "column_sample": str(column_sample),
                    "correlation": (
                        float(correlation) if defined else float("nan")
                    ),
                    "defined": bool(defined),
                    "correlation_label": (
                        f"{float(correlation):.3f}" if defined else "N/A"
                    ),
                }
            )
    return pd.DataFrame.from_records(records, columns=HEATMAP_DATA_COLUMNS)


def _require_dataframe(table: object, label: str) -> None:
    if not isinstance(table, pd.DataFrame):
        raise CorrelationComputationError(
            CorrelationErrorReason.INVALID_TABLE,
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
        raise CorrelationComputationError(
            CorrelationErrorReason.MISSING_REQUIRED_COLUMN,
            f"{label} is missing required column '{column}'.",
        )
    if count > 1:
        raise CorrelationComputationError(
            CorrelationErrorReason.DUPLICATE_REQUIRED_COLUMN,
            f"{label} contains {count} columns named '{column}'.",
        )


def _column_name_matches(value: object, expected: str) -> bool:
    return isinstance(value, str) and value == expected


def _column_identifiers(columns: list[object]) -> list[str]:
    if any(_is_missing_or_blank(value) for value in columns):
        raise CorrelationComputationError(
            CorrelationErrorReason.MISSING_REQUIRED_IDENTIFIER,
            "Expression sample-column identifiers must be non-missing and "
            "non-blank.",
        )
    identifiers = [str(value) for value in columns]
    duplicates = _duplicates(identifiers)
    if duplicates:
        raise CorrelationComputationError(
            CorrelationErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
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
        raise CorrelationComputationError(
            CorrelationErrorReason.MISSING_REQUIRED_IDENTIFIER,
            f"{label} contain {len(missing_positions)} missing or blank "
            "value(s) at zero-based row position(s): "
            + ", ".join(str(position) for position in missing_positions[:5])
            + ".",
        )
    identifiers = [str(value) for value in values]
    duplicates = _duplicates(identifiers)
    if duplicates:
        raise CorrelationComputationError(
            CorrelationErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
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
        raise CorrelationComputationError(
            CorrelationErrorReason.MISSING_REQUIRED_VALUE,
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
    raise CorrelationComputationError(
        CorrelationErrorReason.SAMPLE_MISMATCH,
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
        raise CorrelationComputationError(
            CorrelationErrorReason.MISSING_EXPRESSION_VALUE,
            f"The expression matrix contains {missing_count} missing or blank "
            "sample value(s).",
        )

    boolean_mask = working.apply(lambda column: column.map(is_bool))
    complex_mask = working.apply(lambda column: column.map(is_complex))
    boolean_count = int(boolean_mask.sum().sum())
    complex_count = int(complex_mask.sum().sum())
    if boolean_count or complex_count:
        raise CorrelationComputationError(
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
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
        raise CorrelationComputationError(
            CorrelationErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            f"The expression matrix contains {non_coercible_count} "
            "non-coercible sample value(s).",
        )

    positive_infinity_count = int(numeric.eq(float("inf")).sum().sum())
    negative_infinity_count = int(numeric.eq(float("-inf")).sum().sum())
    if positive_infinity_count or negative_infinity_count:
        raise CorrelationComputationError(
            CorrelationErrorReason.INFINITE_EXPRESSION_VALUE,
            "The expression matrix contains "
            f"{positive_infinity_count} positive and "
            f"{negative_infinity_count} negative infinite sample value(s).",
        )

    numeric.columns = sample_ids
    return numeric


def _pearson_matrix(
    numeric_expression: pd.DataFrame,
    sample_ids: list[str],
    constant_samples: tuple[str, ...],
) -> pd.DataFrame:
    matrix = numeric_expression.corr(
        method="pearson",
        min_periods=len(numeric_expression.index),
    )
    matrix = matrix.reindex(index=sample_ids, columns=sample_ids).copy(deep=True)
    constant_set = set(constant_samples)
    for sample_id in sample_ids:
        if sample_id in constant_set:
            matrix.loc[sample_id, :] = float("nan")
            matrix.loc[:, sample_id] = float("nan")
        else:
            matrix.loc[sample_id, sample_id] = 1.0
    return matrix


def _verify_correlation_matrix(
    matrix: pd.DataFrame,
    sample_ids: list[str],
    constant_samples: tuple[str, ...],
) -> None:
    expected_shape = (len(sample_ids), len(sample_ids))
    if matrix.shape != expected_shape:
        raise RuntimeError(
            "Pearson matrix shape does not match the expression sample count."
        )
    if list(matrix.index) != sample_ids or list(matrix.columns) != sample_ids:
        raise RuntimeError(
            "Pearson matrix axes do not preserve expression sample order."
        )

    for row_position in range(len(sample_ids)):
        for column_position in range(len(sample_ids)):
            value = matrix.iat[row_position, column_position]
            mirror = matrix.iat[column_position, row_position]
            if pd.isna(value) and pd.isna(mirror):
                continue
            if pd.isna(value) != pd.isna(mirror):
                raise RuntimeError(
                    "Pearson matrix contains asymmetric undefined values."
                )
            if abs(float(value) - float(mirror)) > _MATRIX_TOLERANCE:
                raise RuntimeError("Pearson matrix is not symmetric.")
            if not (
                -1.0 - _MATRIX_TOLERANCE
                <= float(value)
                <= 1.0 + _MATRIX_TOLERANCE
            ):
                raise RuntimeError(
                    "Pearson matrix contains a value outside [-1, 1]."
                )

    constant_set = set(constant_samples)
    for sample_id in sample_ids:
        if sample_id in constant_set:
            if not (
                matrix.loc[sample_id].isna().all()
                and matrix.loc[:, sample_id].isna().all()
            ):
                raise RuntimeError(
                    "A constant sample has a defined Pearson correlation."
                )
        else:
            diagonal = matrix.loc[sample_id, sample_id]
            if pd.isna(diagonal) or (
                abs(float(diagonal) - 1.0) > _MATRIX_TOLERANCE
            ):
                raise RuntimeError(
                    "A non-constant sample has an invalid diagonal "
                    "correlation."
                )


def _summary_value(values: pd.Series, operation: str) -> float:
    if values.empty:
        return float("nan")
    if operation == "min":
        return float(values.min())
    if operation == "median":
        return float(values.median())
    if operation == "mean":
        return float(values.mean())
    if operation == "max":
        return float(values.max())
    raise ValueError(f"Unsupported summary operation: {operation!r}.")


def _condition_pairs(
    pair_summary: pd.DataFrame,
    condition_a: str,
    condition_b: str,
) -> pd.DataFrame:
    if condition_a == condition_b:
        mask = pair_summary["condition_a"].eq(condition_a) & pair_summary[
            "condition_b"
        ].eq(condition_b)
    else:
        forward = pair_summary["condition_a"].eq(
            condition_a
        ) & pair_summary["condition_b"].eq(condition_b)
        reverse = pair_summary["condition_a"].eq(
            condition_b
        ) & pair_summary["condition_b"].eq(condition_a)
        mask = forward | reverse
    return pair_summary.loc[mask].copy(deep=True)


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
