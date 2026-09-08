"""Pure exploration of supplied differential-expression results."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from math import isfinite
from numbers import Real

import pandas as pd
from pandas.api.types import is_bool

from plant_expression_explorer.validation import ValidationReport, validate_de_results


STATUS_COLUMN = "exploratory_threshold_status"


class DifferentialExpressionStatus(StrEnum):
    """Stable mutually exclusive exploratory row classifications."""

    NOT_EVALUABLE = "NOT_EVALUABLE"
    POSITIVE_THRESHOLD_MATCH = "POSITIVE_THRESHOLD_MATCH"
    NEGATIVE_THRESHOLD_MATCH = "NEGATIVE_THRESHOLD_MATCH"
    DOES_NOT_MEET_COMBINED_THRESHOLDS = (
        "DOES_NOT_MEET_COMBINED_THRESHOLDS"
    )


class DifferentialExpressionErrorReason(StrEnum):
    """Stable reasons for expected differential-expression exploration failures."""

    BLOCKING_DE_VALIDATION = "BLOCKING_DE_VALIDATION"
    INVALID_ADJUSTED_P_VALUE_THRESHOLD = (
        "INVALID_ADJUSTED_P_VALUE_THRESHOLD"
    )
    INVALID_ABSOLUTE_LOG2_FOLD_CHANGE_THRESHOLD = (
        "INVALID_ABSOLUTE_LOG2_FOLD_CHANGE_THRESHOLD"
    )
    STATUS_COLUMN_CONFLICT = "STATUS_COLUMN_CONFLICT"


class DifferentialExpressionComputationError(ValueError):
    """Expected failure when supplied results cannot be explored safely."""

    def __init__(
        self,
        reason: DifferentialExpressionErrorReason,
        message: str,
        *,
        validation_report: ValidationReport | None = None,
    ) -> None:
        self.reason = reason
        self.validation_report = validation_report
        super().__init__(message)


@dataclass(frozen=True)
class DifferentialExpressionResult:
    """Complete threshold classification for one supplied results table.

    ``frozen=True`` prevents field rebinding but does not make the newly
    constructed annotated DataFrame deeply immutable. Callers should treat it
    as read-only and use :func:`select_rows_by_status` for independent views.
    """

    adjusted_p_value_threshold: float
    absolute_log2_fold_change_threshold: float
    total_row_count: int
    evaluable_row_count: int
    not_evaluable_count: int
    positive_threshold_match_count: int
    negative_threshold_match_count: int
    does_not_meet_combined_thresholds_count: int
    annotated_results: pd.DataFrame


def classify_differential_expression_results(
    de_results: pd.DataFrame,
    *,
    adjusted_p_value_threshold: object = 0.05,
    absolute_log2_fold_change_threshold: object = 1.0,
) -> DifferentialExpressionResult:
    """Classify every accepted source row using two exploratory thresholds.

    The existing validator supplies the complete input contract. Numeric
    coercion is limited to temporary Series used for comparisons; the caller's
    table is never rewritten, filtered, sorted, or otherwise transformed.
    """

    report = validate_de_results(de_results)
    if report.has_errors:
        raise DifferentialExpressionComputationError(
            DifferentialExpressionErrorReason.BLOCKING_DE_VALIDATION,
            "Differential-expression results contain blocking validation Errors.",
            validation_report=report,
        )

    if STATUS_COLUMN in de_results.columns:
        raise DifferentialExpressionComputationError(
            DifferentialExpressionErrorReason.STATUS_COLUMN_CONFLICT,
            f"The supplied table already contains reserved column '{STATUS_COLUMN}'.",
        )

    adjusted_threshold = _validated_threshold(
        adjusted_p_value_threshold,
        reason=DifferentialExpressionErrorReason.INVALID_ADJUSTED_P_VALUE_THRESHOLD,
        label="Adjusted p-value threshold",
        minimum=0.0,
        maximum=1.0,
        minimum_is_inclusive=True,
    )
    fold_change_threshold = _validated_threshold(
        absolute_log2_fold_change_threshold,
        reason=(
            DifferentialExpressionErrorReason.INVALID_ABSOLUTE_LOG2_FOLD_CHANGE_THRESHOLD
        ),
        label="Absolute log2 fold-change threshold",
        minimum=0.0,
        maximum=None,
        minimum_is_inclusive=False,
    )

    numeric_padj = pd.to_numeric(de_results["padj"], errors="coerce")
    numeric_fold_change = pd.to_numeric(
        de_results["log2FoldChange"], errors="coerce"
    )
    not_evaluable = numeric_padj.isna() | numeric_fold_change.isna()
    evaluable = ~not_evaluable
    adjusted_matches = numeric_padj.le(adjusted_threshold).fillna(False)
    positive_matches = (
        evaluable
        & adjusted_matches
        & numeric_fold_change.ge(fold_change_threshold).fillna(False)
    )
    negative_matches = (
        evaluable
        & adjusted_matches
        & numeric_fold_change.le(-fold_change_threshold).fillna(False)
    )
    other_evaluable = evaluable & ~positive_matches & ~negative_matches

    status_values = [
        DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS.value
    ] * len(de_results.index)
    for position in range(len(status_values)):
        if bool(not_evaluable.iloc[position]):
            status_values[position] = DifferentialExpressionStatus.NOT_EVALUABLE.value
        elif bool(positive_matches.iloc[position]):
            status_values[position] = (
                DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH.value
            )
        elif bool(negative_matches.iloc[position]):
            status_values[position] = (
                DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH.value
            )

    annotated = de_results.copy(deep=True)
    annotated[STATUS_COLUMN] = pd.array(status_values, dtype="string")

    counts = {
        DifferentialExpressionStatus.NOT_EVALUABLE: int(not_evaluable.sum()),
        DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH: int(
            positive_matches.sum()
        ),
        DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH: int(
            negative_matches.sum()
        ),
        DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS: int(
            other_evaluable.sum()
        ),
    }
    _assert_result_invariants(de_results, annotated, counts)

    return DifferentialExpressionResult(
        adjusted_p_value_threshold=adjusted_threshold,
        absolute_log2_fold_change_threshold=fold_change_threshold,
        total_row_count=len(de_results.index),
        evaluable_row_count=int(evaluable.sum()),
        not_evaluable_count=counts[DifferentialExpressionStatus.NOT_EVALUABLE],
        positive_threshold_match_count=counts[
            DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH
        ],
        negative_threshold_match_count=counts[
            DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH
        ],
        does_not_meet_combined_thresholds_count=counts[
            DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS
        ],
        annotated_results=annotated,
    )


def build_category_summary(
    result: DifferentialExpressionResult,
) -> pd.DataFrame:
    """Return category counts in stable scientific-contract order."""

    records = (
        {
            "status": DifferentialExpressionStatus.NOT_EVALUABLE.value,
            "row_count": result.not_evaluable_count,
        },
        {
            "status": DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH.value,
            "row_count": result.positive_threshold_match_count,
        },
        {
            "status": DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH.value,
            "row_count": result.negative_threshold_match_count,
        },
        {
            "status": (
                DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS.value
            ),
            "row_count": result.does_not_meet_combined_thresholds_count,
        },
    )
    return pd.DataFrame.from_records(records, columns=("status", "row_count"))


def select_rows_by_status(
    result: DifferentialExpressionResult,
    status: DifferentialExpressionStatus,
) -> pd.DataFrame:
    """Return an independent, order-preserving copy for one status."""

    if not isinstance(status, DifferentialExpressionStatus):
        raise TypeError("status must be a DifferentialExpressionStatus value.")
    mask = result.annotated_results[STATUS_COLUMN].eq(status.value)
    return result.annotated_results.loc[mask].copy(deep=True)


def _validated_threshold(
    value: object,
    *,
    reason: DifferentialExpressionErrorReason,
    label: str,
    minimum: float,
    maximum: float | None,
    minimum_is_inclusive: bool,
) -> float:
    if is_bool(value) or not isinstance(value, Real):
        raise DifferentialExpressionComputationError(
            reason,
            f"{label} must be a finite real number.",
        )
    try:
        threshold = float(value)
    except (OverflowError, TypeError, ValueError) as error:
        raise DifferentialExpressionComputationError(
            reason,
            f"{label} must be a finite real number.",
        ) from error
    if not isfinite(threshold):
        raise DifferentialExpressionComputationError(
            reason,
            f"{label} must be finite.",
        )

    below_minimum = (
        threshold < minimum
        if minimum_is_inclusive
        else threshold <= minimum
    )
    if below_minimum or (maximum is not None and threshold > maximum):
        interval = (
            f"[{minimum:g}, {maximum:g}]"
            if minimum_is_inclusive and maximum is not None
            else f"greater than {minimum:g}"
        )
        raise DifferentialExpressionComputationError(
            reason,
            f"{label} must be {interval}.",
        )
    return threshold


def _assert_result_invariants(
    source: pd.DataFrame,
    annotated: pd.DataFrame,
    counts: dict[DifferentialExpressionStatus, int],
) -> None:
    total = len(source.index)
    if len(annotated.index) != total or sum(counts.values()) != total:
        raise RuntimeError("Differential-expression row-count invariant failed.")
    if not annotated.index.equals(source.index):
        raise RuntimeError("Differential-expression index invariant failed.")
    if annotated.columns[:-1].tolist() != source.columns.tolist():
        raise RuntimeError("Differential-expression column-order invariant failed.")

    preserved_source = annotated.iloc[:, :-1]
    if not preserved_source.equals(source):
        raise RuntimeError("Differential-expression source-value invariant failed.")
    if not preserved_source.dtypes.equals(source.dtypes):
        raise RuntimeError("Differential-expression source-dtype invariant failed.")

    statuses = annotated[STATUS_COLUMN]
    membership = pd.DataFrame(
        {
            status.value: statuses.eq(status.value)
            for status in DifferentialExpressionStatus
        },
        index=annotated.index,
    )
    if not membership.sum(axis=1).eq(1).all():
        raise RuntimeError(
            "Differential-expression status-exclusivity invariant failed."
        )

    status_counts = statuses.value_counts().to_dict()
    expected_counts = {status.value: count for status, count in counts.items()}
    if any(status_counts.get(status, 0) != count for status, count in expected_counts.items()):
        raise RuntimeError("Differential-expression status-count invariant failed.")
    if set(status_counts) - set(expected_counts):
        raise RuntimeError("Differential-expression status-exhaustiveness invariant failed.")
