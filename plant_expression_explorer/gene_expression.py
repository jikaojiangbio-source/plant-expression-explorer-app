"""Pure descriptive lookup of one supplied gene's expression values."""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import StrEnum

import pandas as pd
from pandas.api.types import is_bool, is_complex

from plant_expression_explorer.validation import IssueCode, ValidationReport


class GeneExpressionErrorReason(StrEnum):
    """Stable reasons for expected gene-expression lookup failures."""

    INVALID_TABLE = "INVALID_TABLE"
    EMPTY_EXPRESSION_TABLE = "EMPTY_EXPRESSION_TABLE"
    NO_EXPRESSION_SAMPLE_COLUMNS = "NO_EXPRESSION_SAMPLE_COLUMNS"
    MISSING_REQUIRED_COLUMN = "MISSING_REQUIRED_COLUMN"
    DUPLICATE_REQUIRED_COLUMN = "DUPLICATE_REQUIRED_COLUMN"
    MISSING_REQUIRED_IDENTIFIER = "MISSING_REQUIRED_IDENTIFIER"
    WHITESPACE_IN_IDENTIFIER = "WHITESPACE_IN_IDENTIFIER"
    DUPLICATE_REQUIRED_IDENTIFIER = "DUPLICATE_REQUIRED_IDENTIFIER"
    MISSING_REQUIRED_VALUE = "MISSING_REQUIRED_VALUE"
    WHITESPACE_IN_REQUIRED_VALUE = "WHITESPACE_IN_REQUIRED_VALUE"
    SAMPLE_MISMATCH = "SAMPLE_MISMATCH"
    INVALID_GENE_ID = "INVALID_GENE_ID"
    UNKNOWN_GENE_ID = "UNKNOWN_GENE_ID"
    MISSING_EXPRESSION_VALUE = "MISSING_EXPRESSION_VALUE"
    NON_COERCIBLE_EXPRESSION_VALUE = "NON_COERCIBLE_EXPRESSION_VALUE"
    INFINITE_EXPRESSION_VALUE = "INFINITE_EXPRESSION_VALUE"
    NUMERICAL_RANGE_ERROR = "NUMERICAL_RANGE_ERROR"


class GeneExpressionComputationError(ValueError):
    """Expected failure when a gene-expression lookup cannot proceed safely."""

    def __init__(
        self,
        reason: GeneExpressionErrorReason,
        message: str,
    ) -> None:
        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class GeneExpressionResult:
    """Descriptive expression values and summaries for one exact gene ID.

    Identifier matching uses exact validated ``str(value)`` representations,
    consistent with the application's existing identifier contract. The source
    DataFrames and their original scalar values are never rewritten. The
    ``expression_value`` column in ``sample_expression`` contains copied source
    scalars, while condition summaries and chart data use separate temporary
    numeric copies.

    ``frozen=True`` prevents field rebinding but does not make the newly
    constructed nested DataFrames deeply immutable. Callers should treat them
    as read-only and use the provided helper for independent chart data.
    """

    gene_id: str
    sample_count: int
    condition_count: int
    metadata_order_matches_expression: bool
    all_values_equal: bool
    sample_expression: pd.DataFrame
    condition_summary: pd.DataFrame


SAMPLE_EXPRESSION_COLUMNS = (
    "sample_id",
    "condition",
    "expression_value",
)
CONDITION_EXPRESSION_SUMMARY_COLUMNS = (
    "condition",
    "sample_count",
    "sample_ids",
    "minimum_expression",
    "median_expression",
    "mean_expression",
    "maximum_expression",
    "standard_deviation",
)
MULTI_GENE_PANEL_COLUMNS = ("gene_id", "sample_id", "condition", "expression_value")
TIME_SERIES_CHART_COLUMNS = ("sample_id", "condition", "time_value", "expression_value")
GENE_EXPRESSION_CHART_COLUMNS = (
    "sample_position",
    "sample_id",
    "condition",
    "expression_value",
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
_NUMERICAL_RANGE_MESSAGE = (
    "Gene-expression summaries could not be represented safely as finite "
    "float64 values. The finite source values were not modified; no partial "
    "summary was returned."
)


def list_gene_ids(expression: pd.DataFrame) -> tuple[str, ...]:
    """Return exact validated gene-ID representations in source row order.

    This function defensively validates the complete expression table relevant
    to lookup, including its sample values. It intentionally does not enforce
    the application's separate two-sample activation rule, because a pure
    single-gene lookup has an explicit one-sample result.
    """

    gene_ids, _sample_columns, _sample_ids, _numeric = _validated_expression(
        expression
    )
    return tuple(gene_ids)


def filter_gene_ids(gene_ids: tuple[str, ...], query: str) -> tuple[str, ...]:
    """Return gene IDs containing ``query`` as a case-insensitive substring.

    An empty or whitespace-only query returns every supplied ID unchanged.
    Matching does not trim, normalize, or reorder identifiers; source row
    order is preserved.
    """

    stripped_query = query.strip()
    if not stripped_query:
        return gene_ids
    lowered_query = stripped_query.casefold()
    return tuple(
        gene_id for gene_id in gene_ids if lowered_query in gene_id.casefold()
    )


def lookup_gene_expression(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_id: object,
) -> GeneExpressionResult:
    """Return a non-mutating descriptive lookup for one exact supplied gene ID.

    The expression and metadata tables are both defensively validated. No
    identifiers are trimmed or case-folded, no samples are intersected, and no
    missing or invalid expression cells are skipped or imputed.
    """

    gene_ids, sample_columns, sample_ids, numeric_expression = (
        _validated_expression(expression)
    )
    metadata_sample_ids, conditions = _validated_metadata(metadata)
    _require_matching_samples(sample_ids, metadata_sample_ids)

    selected_gene_id = _validated_selected_gene_id(gene_id)
    try:
        gene_position = gene_ids.index(selected_gene_id)
    except ValueError as error:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.UNKNOWN_GENE_ID,
            f"Gene ID {selected_gene_id!r} is not present in the expression matrix.",
        ) from error

    condition_by_sample = dict(
        zip(metadata_sample_ids, conditions, strict=True)
    )
    sample_conditions = [condition_by_sample[sample_id] for sample_id in sample_ids]

    source_values = [
        expression.iloc[gene_position][column] for column in sample_columns
    ]
    sample_expression = pd.DataFrame(
        {
            "sample_id": pd.Series(sample_ids, dtype="object"),
            "condition": pd.Series(sample_conditions, dtype="object"),
            "expression_value": pd.Series(source_values, dtype="object"),
        },
        columns=SAMPLE_EXPRESSION_COLUMNS,
    )
    numeric_values = numeric_expression.iloc[gene_position].copy(deep=True)
    condition_summary = _build_condition_summary(
        sample_ids,
        sample_conditions,
        numeric_values,
    )

    result = GeneExpressionResult(
        gene_id=selected_gene_id,
        sample_count=len(sample_ids),
        condition_count=len(condition_summary.index),
        metadata_order_matches_expression=(metadata_sample_ids == sample_ids),
        all_values_equal=bool(numeric_values.nunique(dropna=False) == 1),
        sample_expression=sample_expression,
        condition_summary=condition_summary,
    )
    _assert_result_invariants(result, sample_ids, sample_conditions, source_values)
    return result


def build_multi_gene_panel_data(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_ids: list[object],
) -> pd.DataFrame:
    """Return long-format supplied per-sample values for several exact genes.

    Each gene is looked up independently through :func:`lookup_gene_expression`,
    reusing its full validation and numeric-coercion contract; no additional
    computation, normalization, or cross-gene scaling is introduced. Rows are
    concatenated in the supplied gene order, then each gene's own sample
    order. Different genes may have very different absolute scales; plotting
    them together does not make their magnitudes comparable.
    """

    frames: list[pd.DataFrame] = []
    for gene_id in gene_ids:
        result = lookup_gene_expression(expression, metadata, gene_id)
        gene_frame = result.sample_expression.loc[
            :, ["sample_id", "condition", "expression_value"]
        ].copy(deep=True)
        gene_frame.insert(0, "gene_id", result.gene_id)
        frames.append(gene_frame)
    if not frames:
        return pd.DataFrame(columns=MULTI_GENE_PANEL_COLUMNS)
    combined = pd.concat(frames, ignore_index=True)
    return combined.loc[:, list(MULTI_GENE_PANEL_COLUMNS)]


def build_gene_expression_chart_data(
    result: GeneExpressionResult,
) -> pd.DataFrame:
    """Return independent numeric point-chart data in expression sample order."""

    numeric_values = pd.to_numeric(
        result.sample_expression["expression_value"],
        errors="raise",
    )
    return pd.DataFrame(
        {
            "sample_position": range(result.sample_count),
            "sample_id": result.sample_expression["sample_id"].tolist(),
            "condition": result.sample_expression["condition"].tolist(),
            "expression_value": numeric_values.tolist(),
        },
        columns=GENE_EXPRESSION_CHART_COLUMNS,
    )


def build_grouped_gene_expression_chart_data(
    result: GeneExpressionResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return independent numeric point-chart data labelled by one metadata column.

    ``group_column`` values are joined from ``metadata`` by exact sample ID,
    for display only, exactly like 'condition'. Reuses
    :func:`build_gene_expression_chart_data`'s numeric coercion; the plotted
    values are unchanged.
    """

    if group_column == "condition":
        return build_gene_expression_chart_data(result)
    lookup = _metadata_column_lookup(metadata, group_column)
    numeric_values = pd.to_numeric(
        result.sample_expression["expression_value"],
        errors="raise",
    )
    sample_ids = result.sample_expression["sample_id"].tolist()
    return pd.DataFrame(
        {
            "sample_position": range(result.sample_count),
            "sample_id": sample_ids,
            group_column: [lookup[sample_id] for sample_id in sample_ids],
            "expression_value": numeric_values.tolist(),
        },
        columns=("sample_position", "sample_id", group_column, "expression_value"),
    )


def build_time_series_chart_data(
    result: GeneExpressionResult,
    metadata: pd.DataFrame,
    time_column: str,
) -> pd.DataFrame:
    """Return per-sample values stably sorted by one numeric metadata column.

    ``time_column`` values are joined from ``metadata`` by exact sample ID
    and must be fully numeric (safely coercible via ``pandas.to_numeric``);
    this raises ``ValueError`` naming every sample with a non-numeric or
    missing value rather than silently dropping, imputing, or excluding it.
    Samples are stably sorted by the numeric time value: identical-time
    samples (for example biological replicates at one timepoint) keep their
    original relative order and are never averaged or otherwise combined.
    """

    lookup = _metadata_column_lookup(metadata, time_column)
    sample_ids = result.sample_expression["sample_id"].tolist()
    raw_time_values = [lookup[sample_id] for sample_id in sample_ids]
    numeric_time = pd.to_numeric(pd.Series(raw_time_values), errors="coerce")
    invalid_samples = [
        sample_id
        for sample_id, value in zip(sample_ids, numeric_time, strict=True)
        if pd.isna(value)
    ]
    if invalid_samples:
        raise ValueError(
            f"Metadata column '{time_column}' is not numeric for sample(s): "
            + ", ".join(invalid_samples)
            + "."
        )
    numeric_values = pd.to_numeric(
        result.sample_expression["expression_value"], errors="raise"
    )
    frame = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "condition": result.sample_expression["condition"].tolist(),
            "time_value": numeric_time.tolist(),
            "expression_value": numeric_values.tolist(),
        },
        columns=TIME_SERIES_CHART_COLUMNS,
    )
    return frame.sort_values("time_value", kind="stable", ignore_index=True)


def build_grouped_gene_expression_condition_summary(
    result: GeneExpressionResult,
    metadata: pd.DataFrame,
    group_column: str,
) -> pd.DataFrame:
    """Return the descriptive expression summary for one metadata column.

    Reuses the same descriptive statistics (minimum/median/mean/maximum/
    standard deviation) already used for condition grouping, computed for an
    alternate metadata column chosen by the caller, from the already-copied
    values in ``result.sample_expression``; nothing is re-looked-up or
    re-validated.
    """

    if group_column == "condition":
        return result.condition_summary.copy(deep=True)
    lookup = _metadata_column_lookup(metadata, group_column)
    sample_ids = result.sample_expression["sample_id"].tolist()
    groups = [lookup[sample_id] for sample_id in sample_ids]
    numeric_values = pd.to_numeric(
        result.sample_expression["expression_value"],
        errors="raise",
    )
    summary = _build_condition_summary(sample_ids, groups, numeric_values)
    return summary.rename(columns={"condition": group_column})


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


def build_gene_expression_observations(
    result: GeneExpressionResult,
    validation_report: ValidationReport,
) -> tuple[str, ...]:
    """Build deterministic structural observations without interpretation."""

    observations: list[str] = []
    if not result.metadata_order_matches_expression:
        observations.append(
            "Metadata sample order differs from expression-column order. "
            "Conditions were mapped by exact sample ID without modifying either "
            "source table."
        )
    if result.sample_count == 1:
        observations.append(
            "Only one sample is present, so an across-sample expression pattern "
            "cannot be displayed."
        )
    if result.condition_count == 1:
        observations.append(
            "Only one supplied condition label is present; no between-condition "
            "comparison is available."
        )
    if result.all_values_equal:
        observations.append(
            "The selected gene has exactly identical supplied numeric values "
            "across all samples."
        )

    singleton_conditions = result.condition_summary.loc[
        result.condition_summary["sample_count"].eq(1),
        "condition",
    ].tolist()
    if singleton_conditions:
        observations.append(
            f"{len(singleton_conditions)} condition label(s) contain one sample: "
            + ", ".join(str(value) for value in singleton_conditions)
            + ". Descriptive membership is shown, but replication is not inferred."
        )

    seen_codes: set[IssueCode] = set()
    relevant_codes: list[str] = []
    for issue in (*validation_report.warnings, *validation_report.information):
        if issue.code not in _RELEVANT_VALIDATION_OBSERVATIONS:
            continue
        if issue.code in seen_codes:
            continue
        seen_codes.add(issue.code)
        relevant_codes.append(issue.code.value)
    if relevant_codes:
        observations.append(
            "Relevant non-blocking validation code(s): "
            + ", ".join(relevant_codes)
            + "."
        )
    return tuple(observations)


def _validated_expression(
    expression: object,
) -> tuple[list[str], list[object], list[str], pd.DataFrame]:
    _require_dataframe(expression, "Expression matrix")
    if len(expression.index) == 0:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.EMPTY_EXPRESSION_TABLE,
            "The expression matrix contains no gene rows.",
        )
    _require_one_column(expression, "gene_id", "Expression matrix")

    sample_columns = [
        column
        for column in expression.columns
        if not _column_name_matches(column, "gene_id")
    ]
    if not sample_columns:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.NO_EXPRESSION_SAMPLE_COLUMNS,
            "The expression matrix contains no sample columns besides 'gene_id'.",
        )

    gene_ids = _series_identifiers(
        expression["gene_id"],
        label="Expression gene IDs",
    )
    sample_ids = _identifier_values(
        sample_columns,
        label="Expression sample-column IDs",
    )
    numeric_expression = _numeric_expression_copy(
        expression,
        sample_columns,
        sample_ids,
    )
    return gene_ids, sample_columns, sample_ids, numeric_expression


def _validated_metadata(metadata: object) -> tuple[list[str], list[str]]:
    _require_dataframe(metadata, "Sample metadata")
    _require_one_column(metadata, "sample_id", "Sample metadata")
    _require_one_column(metadata, "condition", "Sample metadata")
    sample_ids = _series_identifiers(
        metadata["sample_id"],
        label="Metadata sample IDs",
    )
    conditions = _required_values(
        metadata["condition"],
        label="Metadata conditions",
    )
    return sample_ids, conditions


def _build_condition_summary(
    sample_ids: list[str],
    conditions: list[str],
    numeric_values: pd.Series,
) -> pd.DataFrame:
    positions_by_condition: dict[str, list[int]] = {}
    for position, condition in enumerate(conditions):
        positions_by_condition.setdefault(condition, []).append(position)

    records: list[dict[str, object]] = []
    for condition, positions in positions_by_condition.items():
        values = numeric_values.iloc[positions]
        try:
            finite_values = [float(value) for value in values.tolist()]
        except (OverflowError, TypeError, ValueError) as error:
            raise GeneExpressionComputationError(
                GeneExpressionErrorReason.NUMERICAL_RANGE_ERROR,
                _NUMERICAL_RANGE_MESSAGE,
            ) from error
        condition_sample_ids = tuple(sample_ids[position] for position in positions)
        minimum = min(finite_values)
        median = _stable_median(finite_values)
        mean = _stable_mean(finite_values)
        maximum = max(finite_values)
        standard_deviation = _stable_sample_standard_deviation(finite_values)
        _require_finite_summary_value(minimum)
        _require_finite_summary_value(median)
        _require_finite_summary_value(mean)
        _require_finite_summary_value(maximum)
        if len(finite_values) == 1:
            if not math.isnan(standard_deviation):
                raise RuntimeError(
                    "One-sample gene-expression standard deviation must be NaN."
                )
        else:
            _require_finite_summary_value(standard_deviation)
        records.append(
            {
                "condition": condition,
                "sample_count": len(positions),
                "sample_ids": condition_sample_ids,
                "minimum_expression": minimum,
                "median_expression": median,
                "mean_expression": mean,
                "maximum_expression": maximum,
                "standard_deviation": standard_deviation,
            }
        )
    return pd.DataFrame.from_records(
        records,
        columns=CONDITION_EXPRESSION_SUMMARY_COLUMNS,
    )


def _stable_mean(values: list[float]) -> float:
    """Return an overflow-resistant arithmetic mean for finite values."""

    scale = max(abs(value) for value in values)
    if scale == 0.0:
        return 0.0
    try:
        scaled_mean = math.fsum(value / scale for value in values) / len(values)
        mean = scale * scaled_mean
    except (OverflowError, ValueError) as error:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        ) from error
    _require_finite_summary_value(mean)
    return mean


def _stable_median(values: list[float]) -> float:
    """Return a median without overflowing the midpoint of two finite values."""

    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]

    lower = ordered[middle - 1]
    upper = ordered[middle]
    if (lower >= 0.0 and upper >= 0.0) or (lower <= 0.0 and upper <= 0.0):
        median = lower + (upper - lower) / 2.0
    else:
        median = (lower + upper) / 2.0
    _require_finite_summary_value(median)
    return median


def _stable_sample_standard_deviation(values: list[float]) -> float:
    """Return scaled sample SD, or NaN for exactly one finite value."""

    if len(values) == 1:
        return float("nan")
    if all(value == values[0] for value in values[1:]):
        return 0.0

    scale = max(abs(value) for value in values)
    if scale == 0.0:
        return 0.0
    try:
        scaled_values = [value / scale for value in values]
        scaled_mean = _stable_mean(scaled_values)
        squared_deviations = math.fsum(
            (value - scaled_mean) ** 2 for value in scaled_values
        )
        scaled_standard_deviation = math.sqrt(
            squared_deviations / (len(values) - 1)
        )
        standard_deviation = scale * scaled_standard_deviation
    except (OverflowError, ValueError) as error:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        ) from error
    _require_finite_summary_value(standard_deviation)
    return standard_deviation


def _require_finite_summary_value(value: float) -> None:
    if not math.isfinite(value):
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.NUMERICAL_RANGE_ERROR,
            _NUMERICAL_RANGE_MESSAGE,
        )


def _require_dataframe(table: object, label: str) -> None:
    if not isinstance(table, pd.DataFrame):
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.INVALID_TABLE,
            f"{label} must be a pandas DataFrame.",
        )


def _require_one_column(
    table: pd.DataFrame,
    column: str,
    label: str,
) -> None:
    count = sum(_column_name_matches(value, column) for value in table.columns)
    if count == 0:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.MISSING_REQUIRED_COLUMN,
            f"{label} is missing required column '{column}'.",
        )
    if count > 1:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.DUPLICATE_REQUIRED_COLUMN,
            f"{label} contains {count} columns named '{column}'.",
        )


def _column_name_matches(value: object, expected: str) -> bool:
    return isinstance(value, str) and value == expected


def _series_identifiers(series: pd.Series, *, label: str) -> list[str]:
    return _identifier_values(series.tolist(), label=label)


def _identifier_values(values: list[object], *, label: str) -> list[str]:
    missing_positions = [
        position
        for position, value in enumerate(values)
        if _is_missing_or_blank(value)
    ]
    if missing_positions:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.MISSING_REQUIRED_IDENTIFIER,
            f"{label} contain {len(missing_positions)} missing or blank value(s) "
            "at zero-based position(s): "
            + ", ".join(str(position) for position in missing_positions[:5])
            + ".",
        )

    identifiers = [str(value) for value in values]
    whitespace = [value for value in identifiers if value != value.strip()]
    if whitespace:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.WHITESPACE_IN_IDENTIFIER,
            f"{label} contain leading or trailing whitespace: "
            + ", ".join(dict.fromkeys(whitespace))
            + ". Validation did not trim them.",
        )
    duplicates = _duplicates(identifiers)
    if duplicates:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.DUPLICATE_REQUIRED_IDENTIFIER,
            f"{label} must be unique by exact string representation; duplicate "
            "identifier(s): "
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
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.MISSING_REQUIRED_VALUE,
            f"{label} contain {len(missing_positions)} missing or blank value(s) "
            "at zero-based row position(s): "
            + ", ".join(str(position) for position in missing_positions[:5])
            + ".",
        )
    representations = [str(value) for value in values]
    whitespace = [value for value in representations if value != value.strip()]
    if whitespace:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.WHITESPACE_IN_REQUIRED_VALUE,
            f"{label} contain leading or trailing whitespace: "
            + ", ".join(dict.fromkeys(whitespace))
            + ". Validation did not trim them.",
        )
    return representations


def _validated_selected_gene_id(gene_id: object) -> str:
    if _is_missing_or_blank(gene_id):
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.INVALID_GENE_ID,
            "The selected gene ID must be non-missing and non-blank.",
        )
    representation = str(gene_id)
    if representation != representation.strip():
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.INVALID_GENE_ID,
            "The selected gene ID must not contain leading or trailing whitespace; "
            "it was not trimmed.",
        )
    return representation


def _require_matching_samples(
    sample_ids: list[str],
    metadata_sample_ids: list[str],
) -> None:
    expression_set = set(sample_ids)
    metadata_set = set(metadata_sample_ids)
    if expression_set == metadata_set:
        return

    missing = [sample_id for sample_id in sample_ids if sample_id not in metadata_set]
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
    raise GeneExpressionComputationError(
        GeneExpressionErrorReason.SAMPLE_MISMATCH,
        "Expression and metadata sample IDs differ ("
        + "; ".join(details)
        + "). No intersection was taken.",
    )


def _numeric_expression_copy(
    expression: pd.DataFrame,
    sample_columns: list[object],
    sample_ids: list[str],
) -> pd.DataFrame:
    working = expression.loc[:, sample_columns].copy(deep=True)
    missing_mask = working.apply(lambda column: column.map(_is_missing_or_blank))
    missing_count = int(missing_mask.sum().sum())
    if missing_count:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.MISSING_EXPRESSION_VALUE,
            f"The expression matrix contains {missing_count} missing or blank "
            "sample value(s). No values were imputed or omitted.",
        )

    boolean_mask = working.apply(lambda column: column.map(is_bool))
    complex_mask = working.apply(lambda column: column.map(is_complex))
    boolean_count = int(boolean_mask.sum().sum())
    complex_count = int(complex_mask.sum().sum())
    if boolean_count or complex_count:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            "The expression matrix contains "
            f"{boolean_count} boolean and {complex_count} complex sample value(s); "
            "expression values must be real-valued numbers or safely coercible "
            "numeric strings.",
        )

    numeric = working.apply(lambda column: pd.to_numeric(column, errors="coerce"))
    non_coercible_count = int(numeric.isna().sum().sum())
    if non_coercible_count:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.NON_COERCIBLE_EXPRESSION_VALUE,
            f"The expression matrix contains {non_coercible_count} non-coercible "
            "sample value(s). No values were omitted.",
        )

    positive_infinity_count = int(numeric.eq(float("inf")).sum().sum())
    negative_infinity_count = int(numeric.eq(float("-inf")).sum().sum())
    if positive_infinity_count or negative_infinity_count:
        raise GeneExpressionComputationError(
            GeneExpressionErrorReason.INFINITE_EXPRESSION_VALUE,
            "The expression matrix contains "
            f"{positive_infinity_count} positive and {negative_infinity_count} "
            "negative infinite sample value(s).",
        )
    numeric.columns = sample_ids
    return numeric


def _assert_result_invariants(
    result: GeneExpressionResult,
    sample_ids: list[str],
    conditions: list[str],
    source_values: list[object],
) -> None:
    if result.sample_count != len(sample_ids):
        raise RuntimeError("Gene-expression sample-count invariant failed.")
    if result.sample_expression.columns.tolist() != list(SAMPLE_EXPRESSION_COLUMNS):
        raise RuntimeError("Gene-expression sample-table schema invariant failed.")
    if result.sample_expression["sample_id"].tolist() != sample_ids:
        raise RuntimeError("Gene-expression sample-order invariant failed.")
    if result.sample_expression["condition"].tolist() != conditions:
        raise RuntimeError("Gene-expression condition-mapping invariant failed.")
    actual_values = result.sample_expression["expression_value"].tolist()
    if len(actual_values) != len(source_values) or any(
        not _scalar_values_equal(actual, expected)
        for actual, expected in zip(actual_values, source_values, strict=True)
    ):
        raise RuntimeError("Gene-expression source-value invariant failed.")
    if int(result.condition_summary["sample_count"].sum()) != result.sample_count:
        raise RuntimeError("Gene-expression condition-count invariant failed.")


def _scalar_values_equal(first: object, second: object) -> bool:
    try:
        if bool(pd.isna(first)) and bool(pd.isna(second)):
            return True
    except (TypeError, ValueError):
        pass
    try:
        return bool(first == second)
    except (TypeError, ValueError):
        return False


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
