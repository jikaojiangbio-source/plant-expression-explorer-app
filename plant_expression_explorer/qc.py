"""Pure descriptive sample quality-control calculations."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum
from typing import Literal, TypeAlias

import pandas as pd

from plant_expression_explorer.validation import IssueCode, ValidationReport


class QcErrorReason(StrEnum):
    """Stable reasons for expected QC input-contract failures."""

    INVALID_TABLE = "INVALID_TABLE"
    EMPTY_EXPRESSION_TABLE = "EMPTY_EXPRESSION_TABLE"
    NO_EXPRESSION_SAMPLE_COLUMNS = "NO_EXPRESSION_SAMPLE_COLUMNS"
    MISSING_REQUIRED_COLUMN = "MISSING_REQUIRED_COLUMN"
    DUPLICATE_REQUIRED_COLUMN = "DUPLICATE_REQUIRED_COLUMN"
    MISSING_REQUIRED_IDENTIFIER = "MISSING_REQUIRED_IDENTIFIER"
    DUPLICATE_REQUIRED_IDENTIFIER = "DUPLICATE_REQUIRED_IDENTIFIER"
    MISSING_REQUIRED_VALUE = "MISSING_REQUIRED_VALUE"
    NON_COERCIBLE_EXPRESSION_VALUE = "NON_COERCIBLE_EXPRESSION_VALUE"
    EMPTY_EXPRESSION_SAMPLE_COLUMN = "EMPTY_EXPRESSION_SAMPLE_COLUMN"
    INFINITE_EXPRESSION_VALUE = "INFINITE_EXPRESSION_VALUE"
    SAMPLE_MISMATCH = "SAMPLE_MISMATCH"


class QcComputationError(ValueError):
    """Expected failure when QC inputs do not satisfy the numerical contract."""

    def __init__(self, reason: QcErrorReason, message: str) -> None:
        """Create an error with a stable reason and concise user-facing message."""

        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class SampleQcResult:
    """Complete descriptive QC result for one expression/metadata pair.

    ``frozen=True`` prevents field rebinding. It does not make the newly
    constructed nested pandas DataFrames intrinsically immutable, so callers
    should treat those result tables as read-only.
    """

    gene_count: int
    sample_count: int
    expression_cell_count: int
    missing_value_count: int
    non_finite_value_count: int
    zero_value_count: int
    negative_value_count: int
    condition_count: int
    metadata_order_matches_expression: bool
    sample_summary: pd.DataFrame
    condition_summary: pd.DataFrame
    constant_samples: tuple[str, ...]
    zero_variance_gene_count: int


SampleChartMetric: TypeAlias = Literal[
    "median",
    "interquartile_range",
    "zero_value_count",
    "negative_value_count",
]
CountChartMetric: TypeAlias = Literal[
    "zero_value_count",
    "negative_value_count",
]

SAMPLE_SUMMARY_COLUMNS = (
    "sample_id",
    "condition",
    "gene_count",
    "missing_value_count",
    "non_finite_value_count",
    "zero_value_count",
    "negative_value_count",
    "minimum",
    "first_quartile",
    "median",
    "mean",
    "third_quartile",
    "maximum",
    "standard_deviation",
    "interquartile_range",
)
CONDITION_SUMMARY_COLUMNS = ("condition", "sample_count", "sample_ids")
_SAMPLE_CHART_METRICS = frozenset(
    {
        "median",
        "interquartile_range",
        "zero_value_count",
        "negative_value_count",
    }
)
_COUNT_CHART_METRICS = frozenset(
    {"zero_value_count", "negative_value_count"}
)
_RELEVANT_VALIDATION_OBSERVATIONS = frozenset(
    {
        IssueCode.SINGLE_CONDITION,
        IssueCode.CONDITION_WITH_SINGLE_SAMPLE,
    }
)


def compute_sample_qc(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
) -> SampleQcResult:
    """Calculate non-mutating descriptive summaries for one validated dataset.

    Safely coercible numeric strings are converted only in a temporary working
    copy. A missing expression cell is retained as missing (never imputed)
    and excluded from any statistic that requires a value for that cell; a
    sample column that is entirely missing, a non-coercible value, or an
    infinite value each cause a :class:`QcComputationError`. No non-missing
    cell, gene, or sample is skipped.
    """

    _require_dataframe(expression, "Expression matrix")
    _require_dataframe(metadata, "Sample metadata")
    if len(expression.index) == 0:
        raise QcComputationError(
            QcErrorReason.EMPTY_EXPRESSION_TABLE,
            "The expression matrix contains no gene rows.",
        )

    _require_one_column(expression, "gene_id", "Expression matrix")
    _require_one_column(metadata, "sample_id", "Sample metadata")
    _require_one_column(metadata, "condition", "Sample metadata")

    sample_columns = [
        column
        for column in expression.columns
        if not _column_name_matches(column, "gene_id")
    ]
    if not sample_columns:
        raise QcComputationError(
            QcErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
            "The expression matrix contains no sample columns.",
        )

    sample_ids = _column_identifiers(sample_columns)
    _series_identifiers(
        expression["gene_id"],
        label="Expression gene IDs",
    )
    metadata_sample_ids = _series_identifiers(
        metadata["sample_id"],
        label="Metadata sample IDs",
    )
    conditions = _required_values(
        metadata["condition"],
        label="Metadata conditions",
    )

    expression_sample_set = set(sample_ids)
    metadata_sample_set = set(metadata_sample_ids)
    if expression_sample_set != metadata_sample_set:
        missing = [
            sample_id
            for sample_id in sample_ids
            if sample_id not in metadata_sample_set
        ]
        unexpected = [
            sample_id
            for sample_id in metadata_sample_ids
            if sample_id not in expression_sample_set
        ]
        details = []
        if missing:
            details.append("missing metadata IDs: " + ", ".join(missing))
        if unexpected:
            details.append("metadata-only IDs: " + ", ".join(unexpected))
        raise QcComputationError(
            QcErrorReason.SAMPLE_MISMATCH,
            "Expression and metadata sample IDs differ ("
            + "; ".join(details)
            + ").",
        )

    numeric_expression = _numeric_expression_copy(
        expression,
        sample_columns,
        sample_ids,
    )
    condition_by_sample = dict(
        zip(metadata_sample_ids, conditions, strict=True)
    )
    sample_conditions = {
        sample_id: condition_by_sample[sample_id] for sample_id in sample_ids
    }

    sample_summary = build_sample_summary(
        numeric_expression,
        sample_conditions,
    )
    condition_summary = build_condition_summary(
        metadata_sample_ids,
        conditions,
    )
    constant_samples = find_constant_samples(numeric_expression)
    zero_variance_gene_count = count_zero_variance_genes(numeric_expression)

    return SampleQcResult(
        gene_count=len(expression.index),
        sample_count=len(sample_ids),
        expression_cell_count=len(expression.index) * len(sample_ids),
        missing_value_count=int(numeric_expression.isna().sum().sum()),
        non_finite_value_count=0,
        zero_value_count=int(numeric_expression.eq(0).sum().sum()),
        negative_value_count=int(numeric_expression.lt(0).sum().sum()),
        condition_count=len(condition_summary.index),
        metadata_order_matches_expression=(
            metadata_sample_ids == sample_ids
        ),
        sample_summary=sample_summary,
        condition_summary=condition_summary,
        constant_samples=constant_samples,
        zero_variance_gene_count=zero_variance_gene_count,
    )


def build_sample_summary(
    numeric_expression: pd.DataFrame,
    condition_by_sample: Mapping[str, str],
) -> pd.DataFrame:
    """Return one newly constructed descriptive-summary row per sample."""

    records: list[dict[str, object]] = []
    gene_count = len(numeric_expression.index)
    for sample_id in numeric_expression.columns:
        values = numeric_expression[sample_id]
        first_quartile = values.quantile(0.25, interpolation="linear")
        median = values.quantile(0.5, interpolation="linear")
        third_quartile = values.quantile(0.75, interpolation="linear")
        records.append(
            {
                "sample_id": str(sample_id),
                "condition": condition_by_sample[str(sample_id)],
                "gene_count": gene_count,
                "missing_value_count": int(values.isna().sum()),
                "non_finite_value_count": int(
                    (values.eq(float("inf")) | values.eq(float("-inf"))).sum()
                ),
                "zero_value_count": int(values.eq(0).sum()),
                "negative_value_count": int(values.lt(0).sum()),
                "minimum": float(values.min()),
                "first_quartile": float(first_quartile),
                "median": float(median),
                "mean": float(values.mean()),
                "third_quartile": float(third_quartile),
                "maximum": float(values.max()),
                "standard_deviation": float(values.std(ddof=1)),
                "interquartile_range": float(
                    third_quartile - first_quartile
                ),
            }
        )
    return pd.DataFrame.from_records(records, columns=SAMPLE_SUMMARY_COLUMNS)


def build_condition_summary(
    sample_ids: list[str],
    conditions: list[str],
) -> pd.DataFrame:
    """Summarize conditions in first-appearance and metadata row order."""

    samples_by_condition: dict[str, list[str]] = {}
    for sample_id, condition in zip(sample_ids, conditions, strict=True):
        samples_by_condition.setdefault(condition, []).append(sample_id)
    records = [
        {
            "condition": condition,
            "sample_count": len(condition_samples),
            "sample_ids": tuple(condition_samples),
        }
        for condition, condition_samples in samples_by_condition.items()
    ]
    return pd.DataFrame.from_records(
        records,
        columns=CONDITION_SUMMARY_COLUMNS,
    )


def list_additional_metadata_columns(metadata: pd.DataFrame) -> tuple[str, ...]:
    """List metadata columns beyond 'sample_id'/'condition', in file order."""

    if not isinstance(metadata, pd.DataFrame):
        return ()
    return tuple(
        str(column)
        for column in metadata.columns
        if str(column) not in ("sample_id", "condition")
    )


def build_grouped_sample_summary(
    result: SampleQcResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return the per-sample summary labelled by one metadata column.

    ``group_column`` values are joined from ``metadata`` by exact sample ID,
    for display only, exactly like 'condition' on ``result.sample_summary``.
    Every per-sample numeric statistic is copied unchanged from
    :func:`compute_sample_qc`; none is recalculated for the new grouping.
    """

    if group_column == "condition":
        return result.sample_summary.copy(deep=True)
    lookup = _metadata_column_lookup(metadata, group_column)
    display = result.sample_summary.copy(deep=True)
    display[group_column] = display["sample_id"].map(lookup)
    columns = [
        group_column if column == "condition" else column
        for column in SAMPLE_SUMMARY_COLUMNS
    ]
    return display.loc[:, columns]


def build_grouped_condition_summary(
    result: SampleQcResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return group membership counts for one metadata column.

    Reuses :func:`build_condition_summary` with an alternate grouping; no
    per-sample numeric statistic is touched.
    """

    if group_column == "condition":
        return result.condition_summary.copy(deep=True)
    lookup = _metadata_column_lookup(metadata, group_column)
    sample_ids = result.sample_summary["sample_id"].tolist()
    groups = [lookup[sample_id] for sample_id in sample_ids]
    summary = build_condition_summary(sample_ids, groups)
    return summary.rename(columns={"condition": group_column})


def build_grouped_condition_chart_data(
    result: SampleQcResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return ordered group counts for charting under one metadata column."""

    summary = build_grouped_condition_summary(result, metadata, group_column)
    return summary.loc[:, [group_column, "sample_count"]].copy(deep=True)


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


def find_constant_samples(
    numeric_expression: pd.DataFrame,
) -> tuple[str, ...]:
    """Return sample IDs whose non-missing values are identical across genes.

    A sample with no non-missing value at all is not considered constant;
    missing values are excluded from the comparison rather than treated as
    a distinct value.
    """

    unique_counts = numeric_expression.nunique(axis=0, dropna=True)
    return tuple(
        str(sample_id)
        for sample_id in numeric_expression.columns
        if unique_counts[sample_id] == 1
    )


def count_zero_variance_genes(numeric_expression: pd.DataFrame) -> int:
    """Count genes whose non-missing values are identical across samples.

    A gene with no non-missing value at all is not counted; missing values
    are excluded from the comparison rather than treated as a distinct
    value.
    """

    return int(numeric_expression.nunique(axis=1, dropna=True).eq(1).sum())


def build_qc_observations(
    result: SampleQcResult,
    validation_report: ValidationReport,
) -> tuple[str, ...]:
    """Build exact, deterministic observations without quality classification."""

    observations: list[str] = []
    if result.gene_count == 1:
        observations.append(
            "The matrix contains only one gene, so sample-spread and "
            "constant-column diagnostics are not informative."
        )
    elif result.constant_samples:
        observations.append(
            f"{len(result.constant_samples)} constant sample column(s) were "
            "found: "
            + ", ".join(result.constant_samples)
            + ". This exact structural finding may warrant further review; "
            "no samples were removed."
        )

    if result.zero_variance_gene_count:
        observations.append(
            f"{result.zero_variance_gene_count} gene(s) have exactly identical "
            "values across all samples. This may affect later scale-dependent "
            "analyses; no genes were removed."
        )
    if not result.metadata_order_matches_expression:
        observations.append(
            "Metadata sample order differs from expression-column order. "
            "Conditions were mapped by exact sample ID without modifying "
            "either table."
        )
    if result.missing_value_count:
        observations.append(
            f"The matrix contains {result.missing_value_count} missing "
            "value(s), retained as missing and excluded from statistics "
            "that require a value for that cell; per-sample counts are "
            "shown in the table above. No value was imputed."
        )
    if result.zero_value_count:
        observations.append(
            f"The matrix contains {result.zero_value_count} zero value(s). "
            "Their meaning depends on the expression scale and may warrant "
            "further review."
        )
    if result.negative_value_count:
        observations.append(
            f"The matrix contains {result.negative_value_count} negative "
            "value(s). Negative values can be valid on transformed or centred "
            "scales and may warrant scale-aware review."
        )

    seen_messages: set[str] = set()
    for issue in validation_report.issues:
        if issue.code not in _RELEVANT_VALIDATION_OBSERVATIONS:
            continue
        if issue.message in seen_messages:
            continue
        seen_messages.add(issue.message)
        observations.append(
            issue.message + " This may warrant further review."
        )
    return tuple(observations)


def build_sample_chart_data(
    result: SampleQcResult,
    metric: SampleChartMetric,
) -> pd.DataFrame:
    """Return ordered sample IDs and one supported metric for charting."""

    if metric not in _SAMPLE_CHART_METRICS:
        raise ValueError(f"Unsupported sample chart metric: {metric!r}.")
    return result.sample_summary.loc[:, ["sample_id", metric]].copy(deep=True)


def build_condition_chart_data(result: SampleQcResult) -> pd.DataFrame:
    """Return ordered condition counts for charting."""

    return result.condition_summary.loc[
        :,
        ["condition", "sample_count"],
    ].copy(deep=True)


def should_display_count_chart(
    result: SampleQcResult,
    metric: CountChartMetric,
) -> bool:
    """Return whether a zero/negative-count chart has a non-zero bar."""

    if metric not in _COUNT_CHART_METRICS:
        raise ValueError(f"Unsupported count chart metric: {metric!r}.")
    return bool(result.sample_summary[metric].gt(0).any())


def _require_dataframe(table: object, label: str) -> None:
    if not isinstance(table, pd.DataFrame):
        raise QcComputationError(
            QcErrorReason.INVALID_TABLE,
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
        raise QcComputationError(
            QcErrorReason.MISSING_REQUIRED_COLUMN,
            f"{label} is missing required column '{column}'.",
        )
    if count > 1:
        raise QcComputationError(
            QcErrorReason.DUPLICATE_REQUIRED_COLUMN,
            f"{label} contains duplicate column '{column}'.",
        )


def _column_name_matches(value: object, expected: str) -> bool:
    return isinstance(value, str) and value == expected


def _column_identifiers(columns: list[object]) -> list[str]:
    identifiers = [column for column in columns]
    if any(_is_missing_or_blank(value) for value in identifiers):
        raise QcComputationError(
            QcErrorReason.MISSING_REQUIRED_IDENTIFIER,
            "Expression sample IDs must be non-missing and non-blank.",
        )
    as_strings = [str(value) for value in identifiers]
    if len(as_strings) != len(set(as_strings)):
        raise QcComputationError(
            QcErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
            "Expression sample IDs must be unique.",
        )
    return as_strings


def _series_identifiers(series: pd.Series, *, label: str) -> list[str]:
    values = series.tolist()
    if any(_is_missing_or_blank(value) for value in values):
        raise QcComputationError(
            QcErrorReason.MISSING_REQUIRED_IDENTIFIER,
            f"{label} must be non-missing and non-blank.",
        )
    as_strings = [str(value) for value in values]
    if len(as_strings) != len(set(as_strings)):
        raise QcComputationError(
            QcErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
            f"{label} must be unique.",
        )
    return as_strings


def _required_values(series: pd.Series, *, label: str) -> list[str]:
    values = series.tolist()
    if any(_is_missing_or_blank(value) for value in values):
        raise QcComputationError(
            QcErrorReason.MISSING_REQUIRED_VALUE,
            f"{label} must be non-missing and non-blank.",
        )
    return [str(value) for value in values]


def _numeric_expression_copy(
    expression: pd.DataFrame,
    sample_columns: list[object],
    sample_ids: list[str],
) -> pd.DataFrame:
    working = expression.loc[:, sample_columns].copy(deep=True)
    missing_mask = working.apply(
        lambda column: column.map(_is_missing_or_blank)
    )
    empty_columns = [
        sample_id
        for column, sample_id in zip(sample_columns, sample_ids, strict=True)
        if missing_mask[column].all()
    ]
    if empty_columns:
        raise QcComputationError(
            QcErrorReason.EMPTY_EXPRESSION_SAMPLE_COLUMN,
            "Sample column(s) contain only missing or blank values: "
            + ", ".join(empty_columns) + ".",
        )

    numeric = working.apply(
        lambda column: pd.to_numeric(column, errors="coerce")
    )
    non_coercible_count = int((numeric.isna() & ~missing_mask).sum().sum())
    if non_coercible_count:
        raise QcComputationError(
            QcErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            f"The expression matrix contains {non_coercible_count} "
            "non-coercible sample value(s).",
        )

    positive_infinity_count = int(numeric.eq(float("inf")).sum().sum())
    negative_infinity_count = int(numeric.eq(float("-inf")).sum().sum())
    if positive_infinity_count or negative_infinity_count:
        raise QcComputationError(
            QcErrorReason.INFINITE_EXPRESSION_VALUE,
            "The expression matrix contains "
            f"{positive_infinity_count} positive and "
            f"{negative_infinity_count} negative infinite sample value(s).",
        )

    numeric.columns = sample_ids
    return numeric


def _is_missing_or_blank(value: object) -> bool:
    if isinstance(value, str):
        return not value.strip()
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
