"""Descriptive exploration of supplied differential-expression results."""

import pandas as pd
import streamlit as st

from plant_expression_explorer.dataset import get_current_dataset
from plant_expression_explorer.differential_expression import (
    DifferentialExpressionComputationError,
    DifferentialExpressionErrorReason,
    DifferentialExpressionStatus,
    build_category_summary,
    classify_differential_expression_results,
    select_rows_by_status,
)
from plant_expression_explorer.validation import Severity, ValidationIssue


_DE_TABLE_NAME = "Differential-expression results"
_CATEGORY_LABELS = {
    DifferentialExpressionStatus.NOT_EVALUABLE.value: (
        "Rows not evaluable with these fields"
    ),
    DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH.value: (
        "Positive fold-change rows meeting both thresholds"
    ),
    DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH.value: (
        "Negative fold-change rows meeting both thresholds"
    ),
    DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS.value: (
        "Evaluable rows that do not meet both thresholds"
    ),
}


def _is_relevant_issue(issue: ValidationIssue) -> bool:
    return issue.table == _DE_TABLE_NAME or issue.related_table == _DE_TABLE_NAME


def _render_issue(issue: ValidationIssue) -> None:
    locations = issue.table
    if issue.related_table is not None:
        locations += f" ↔ {issue.related_table}"
    details = f"**`{issue.code.value}` · {locations}**  \n{issue.message}"
    if issue.severity is Severity.ERROR:
        st.error(details)
    elif issue.severity is Severity.WARNING:
        st.warning(details)
    else:
        st.info(details)


def _show_derived_table(
    title: str,
    table: pd.DataFrame,
    explanation: str,
) -> None:
    st.subheader(title)
    st.write(f"**Rows:** {len(table.index):,}")
    st.dataframe(table, width="stretch")
    st.caption(
        explanation
        + " This is a derived exploratory view in source row order, not a "
        "claim of statistical significance or biological importance."
    )


def _format_applied_threshold(value: float) -> str:
    """Return Python's shortest representation that round-trips to the same float."""

    return repr(float(value))


st.title("🧬 Differential Expression")
st.write(
    "This page applies user-selected exploratory thresholds to supplied, "
    "precomputed differential-expression results. It does not fit a "
    "differential-expression model, calculate or modify p-values, infer the "
    "experimental contrast or reference level, or establish statistical "
    "significance or biological importance."
)
st.caption(
    "Positive and negative labels refer only to the sign of the supplied "
    "log2FoldChange value. Missing adjusted p-values remain missing and are "
    "never replaced with zero."
)

current = get_current_dataset(st.session_state)
if current is None:
    st.info(
        "No validated dataset is currently loaded. Use the Upload Data page "
        "to load the synthetic demo or a complete validated upload."
    )
    st.markdown("Open **Upload Data** from the application navigation.")
    st.stop()

st.header("Current active dataset")
st.write(f"**Source label:** {current.source_label}")
st.write(f"**Source type:** {current.source}")
st.write(f"**Differential-expression rows supplied:** {current.de_row_count:,}")

report = current.validation_report
validation_columns = st.columns(3)
validation_columns[0].metric("Validation Errors", len(report.errors))
validation_columns[1].metric("Validation Warnings", len(report.warnings))
validation_columns[2].metric("Validation Information", len(report.information))

if report.has_errors:
    st.error(
        "Differential-expression exploration is unavailable because the "
        "active dataset contains blocking validation Errors. Return to Upload "
        "Data and correct the source tables."
    )
    for issue in report.errors:
        _render_issue(issue)
    st.stop()

st.header("Relevant validation observations")
relevant_issues = tuple(issue for issue in report.issues if _is_relevant_issue(issue))
if relevant_issues:
    for issue in relevant_issues:
        _render_issue(issue)
else:
    st.write(
        "No differential-expression validation or expression/DEG gene-consistency "
        "issues were reported for the active dataset."
    )

if current.source == "demo":
    st.warning(
        "The demo values are synthetic, its gene IDs are fictional, and its "
        "p-values are constructed. They do not measure uncertainty or support "
        "tomato biological conclusions."
    )

st.header("Exploratory thresholds")
st.write(
    "A row meets a sign-specific threshold only when its supplied adjusted "
    "p-value and fold-change value both meet the inclusive cutoffs. The "
    "cutoffs are user-selected display criteria, not new statistical results."
)
threshold_columns = st.columns(2)
adjusted_p_value_threshold = threshold_columns[0].number_input(
    "Adjusted p-value threshold",
    min_value=0.0,
    max_value=1.0,
    value=0.05,
    step=0.01,
    format="%.17g",
    key="pee_de_adjusted_p_value_threshold",
)
absolute_log2_fold_change_threshold = threshold_columns[1].number_input(
    "Absolute log2 fold-change threshold",
    min_value=0.0,
    value=1.0,
    step=0.1,
    format="%.17g",
    key="pee_de_absolute_log2_fold_change_threshold",
    help="Must be greater than zero so positive and negative categories cannot overlap.",
)

try:
    result = classify_differential_expression_results(
        current.de_results,
        adjusted_p_value_threshold=adjusted_p_value_threshold,
        absolute_log2_fold_change_threshold=absolute_log2_fold_change_threshold,
    )
except DifferentialExpressionComputationError as error:
    st.error(
        "Differential-expression results could not be classified "
        f"({error.reason.value}): {error}"
    )
    if (
        error.reason is DifferentialExpressionErrorReason.BLOCKING_DE_VALIDATION
        and error.validation_report is not None
    ):
        for issue in error.validation_report.errors:
            _render_issue(issue)
    st.stop()

display_adjusted_threshold = _format_applied_threshold(
    result.adjusted_p_value_threshold
)
display_fold_change_threshold = _format_applied_threshold(
    result.absolute_log2_fold_change_threshold
)
display_negative_fold_change_threshold = _format_applied_threshold(
    -result.absolute_log2_fold_change_threshold
)
st.write(
    "**Exact thresholds applied:** "
    f"`padj ≤ {display_adjusted_threshold}`; "
    f"`log2FoldChange ≥ {display_fold_change_threshold}` for positive matches; "
    f"`log2FoldChange ≤ {display_negative_fold_change_threshold}` for negative matches."
)

st.header("Mutually exclusive category counts")
count_columns = st.columns(4)
count_columns[0].metric(
    "Positive threshold matches", result.positive_threshold_match_count
)
count_columns[1].metric(
    "Negative threshold matches", result.negative_threshold_match_count
)
count_columns[2].metric(
    "Other evaluable rows", result.does_not_meet_combined_thresholds_count
)
count_columns[3].metric("Not-evaluable rows", result.not_evaluable_count)
st.write(
    f"**Total supplied rows:** {result.total_row_count:,} · "
    f"**Evaluable rows:** {result.evaluable_row_count:,}"
)

category_summary = build_category_summary(result)
category_summary["category"] = category_summary["status"].map(_CATEGORY_LABELS)
st.dataframe(
    category_summary.loc[:, ["category", "status", "row_count"]],
    hide_index=True,
    width="stretch",
)
st.caption(
    "Each accepted source row appears in exactly one category. Comparisons are "
    "inclusive and use the supplied values without tolerance or imputation."
)

st.header("Complete annotated supplied results")
st.dataframe(result.annotated_results, width="stretch")
st.caption(
    "Every supplied row, index value, source column, value, dtype, and additional "
    "column is preserved in source order. Only the final exploratory status "
    "column is derived; the active source table is not rewritten."
)

positive_rows = select_rows_by_status(
    result, DifferentialExpressionStatus.POSITIVE_THRESHOLD_MATCH
)
negative_rows = select_rows_by_status(
    result, DifferentialExpressionStatus.NEGATIVE_THRESHOLD_MATCH
)
other_rows = select_rows_by_status(
    result, DifferentialExpressionStatus.DOES_NOT_MEET_COMBINED_THRESHOLDS
)
not_evaluable_rows = select_rows_by_status(
    result, DifferentialExpressionStatus.NOT_EVALUABLE
)

_show_derived_table(
    "Positive fold-change rows meeting both thresholds",
    positive_rows,
    (
        f"These rows have padj ≤ {display_adjusted_threshold} and "
        f"log2FoldChange ≥ {display_fold_change_threshold}."
    ),
)
_show_derived_table(
    "Negative fold-change rows meeting both thresholds",
    negative_rows,
    (
        f"These rows have padj ≤ {display_adjusted_threshold} and "
        f"log2FoldChange ≤ {display_negative_fold_change_threshold}."
    ),
)
_show_derived_table(
    "Evaluable rows that do not meet both thresholds",
    other_rows,
    "These rows have usable classification fields but do not meet both current cutoffs.",
)
_show_derived_table(
    "Rows not evaluable with these fields",
    not_evaluable_rows,
    "These rows have a missing or blank padj or log2FoldChange value; missing pvalue alone does not determine this category.",
)

st.header("Scientific and statistical limitations")
st.info(
    "The application uses the supplied padj and log2FoldChange values as-is. It "
    "does not calculate, adjust, replace, or modify p-values, and it does not "
    "run DESeq2 or another differential-expression model."
)
st.info(
    "No contrast direction or reference level is inferred. A positive or "
    "negative supplied fold change cannot be given a condition-specific "
    "interpretation here."
)
st.info(
    "No rows are trimmed, normalised, deduplicated, sorted, ranked, intersected, "
    "aggregated, transformed, imputed, or removed. This phase provides no "
    "volcano plot, gene lookup, or dedicated export workflow."
)
