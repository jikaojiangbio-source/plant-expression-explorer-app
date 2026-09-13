"""Assemble already-computed descriptive result tables into one PDF report."""

from datetime import datetime, timezone

import pandas as pd
import streamlit as st

from plant_expression_explorer.correlation import (
    CorrelationComputationError,
    build_grouped_condition_correlation_summary,
    compute_sample_correlation,
)
from plant_expression_explorer.dataset import (
    ACTIVE_GROUP_COLUMN_KEY,
    ensure_valid_group_column_state,
    get_current_dataset,
)
from plant_expression_explorer.differential_expression import (
    DifferentialExpressionComputationError,
    build_category_summary,
    classify_differential_expression_results,
)
from plant_expression_explorer.pca import PcaComputationError, compute_sample_pca
from plant_expression_explorer.provenance import provenance_display_rows
from plant_expression_explorer.qc import (
    QcComputationError,
    build_grouped_condition_summary,
    compute_sample_qc,
    list_additional_metadata_columns,
)
from plant_expression_explorer.report import (
    PDF_MEDIA_TYPE,
    ReportExportError,
    ReportSection,
    build_report_pdf,
)
from plant_expression_explorer.theme import inject_global_styles

inject_global_styles()

_DE_ADJUSTED_P_KEY = "pee_de_adjusted_p_value_threshold"
_DE_FOLD_CHANGE_KEY = "pee_de_absolute_log2_fold_change_threshold"
_DEFAULT_ADJUSTED_P_THRESHOLD = 0.05
_DEFAULT_FOLD_CHANGE_THRESHOLD = 1.0
_REPORT_TITLE = "Plant Expression Explorer report"
_REPORT_FILENAME = "plant-expression-explorer-report.pdf"


st.title("📄 Report Export")
st.write(
    "Assemble the descriptive result tables already shown on other pages "
    "into one downloadable PDF snapshot of the currently active dataset."
)
st.caption(
    "This page computes nothing new: every table below is produced by the "
    "same functions used on the Sample Quality Control, PCA, Sample "
    "Correlation, and Differential Expression pages, using the same shared "
    "'group by' choice and (for Differential Expression) the same "
    "exploratory thresholds."
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
de_description = (
    f"{current.de_row_count:,} differential-expression rows"
    if current.has_de_results
    else "no differential-expression results supplied"
)
st.write(
    f"**Bundle contents:** {current.gene_count:,} genes, "
    f"{current.sample_count:,} samples, and {de_description}."
)

additional_columns = list_additional_metadata_columns(current.metadata)
group_column = "condition"
if additional_columns:
    ensure_valid_group_column_state(st.session_state, additional_columns)
    group_column = st.selectbox(
        "Group report tables by",
        options=("condition", *additional_columns),
        key=ACTIVE_GROUP_COLUMN_KEY,
        help=(
            "Any column present in the uploaded sample metadata beyond "
            "'sample_id' and 'condition' can relabel the Sample Quality "
            "Control and Sample Correlation tables below. This choice is "
            "shared with the PCA, Sample Quality Control, Sample "
            "Correlation, and Gene Expression pages."
        ),
    )

adjusted_p_value_threshold = st.session_state.get(
    _DE_ADJUSTED_P_KEY, _DEFAULT_ADJUSTED_P_THRESHOLD
)
absolute_log2_fold_change_threshold = st.session_state.get(
    _DE_FOLD_CHANGE_KEY, _DEFAULT_FOLD_CHANGE_THRESHOLD
)


def _dataset_summary_table() -> pd.DataFrame:
    rows = [
        ("Source label", current.source_label),
        ("Source type", current.source),
        ("Genes", f"{current.gene_count:,}"),
        ("Samples", f"{current.sample_count:,}"),
        ("Differential-expression rows", de_description),
    ]
    return pd.DataFrame(rows, columns=["Field", "Value"])


def _dataset_context_table() -> pd.DataFrame:
    rows = provenance_display_rows(current.provenance)
    return pd.DataFrame(rows, columns=["Field", "Value"])


def _qc_section() -> ReportSection:
    title = "Sample quality control — condition summary"
    try:
        result = compute_sample_qc(current.expression, current.metadata)
    except QcComputationError as error:
        return ReportSection(
            title=title,
            note=f"Not included ({error.reason.value}): {error}",
        )
    table = build_grouped_condition_summary(result, current.metadata, group_column)
    return ReportSection(title=title, table=table)


def _pca_section() -> ReportSection:
    title = "PCA — explained variance"
    try:
        result = compute_sample_pca(current.expression, current.metadata)
    except PcaComputationError as error:
        return ReportSection(
            title=title,
            note=f"Not included ({error.reason.value}): {error}",
        )
    return ReportSection(title=title, table=result.variance_table)


def _correlation_section() -> ReportSection:
    title = "Sample correlation — condition-pair summary"
    try:
        result = compute_sample_correlation(current.expression, current.metadata)
    except CorrelationComputationError as error:
        return ReportSection(
            title=title,
            note=f"Not included ({error.reason.value}): {error}",
        )
    table = build_grouped_condition_correlation_summary(
        result, current.metadata, group_column
    )
    return ReportSection(title=title, table=table)


def _differential_expression_section() -> ReportSection:
    title = "Differential expression — category summary"
    if not current.has_de_results:
        return ReportSection(
            title=title,
            note="Not included: no differential-expression results were supplied "
            "for this dataset.",
        )
    try:
        result = classify_differential_expression_results(
            current.de_results,
            adjusted_p_value_threshold=adjusted_p_value_threshold,
            absolute_log2_fold_change_threshold=absolute_log2_fold_change_threshold,
        )
    except DifferentialExpressionComputationError as error:
        return ReportSection(
            title=title,
            note=f"Not included ({error.reason.value}): {error}",
        )
    return ReportSection(title=title, table=build_category_summary(result))


sections = (
    ReportSection(title="Dataset summary", table=_dataset_summary_table()),
    ReportSection(title="Dataset context", table=_dataset_context_table()),
    _qc_section(),
    _pca_section(),
    _correlation_section(),
    _differential_expression_section(),
)

grouping_line = (
    "Sample Quality Control and Sample Correlation tables are grouped by "
    f"metadata column '{group_column}'."
    if group_column != "condition"
    else "Sample Quality Control and Sample Correlation tables are grouped "
    "by 'condition'."
)
de_threshold_line = (
    "Differential-expression thresholds applied: "
    f"padj <= {adjusted_p_value_threshold}; "
    f"|log2FoldChange| >= {absolute_log2_fold_change_threshold}."
    if current.has_de_results
    else "No differential-expression results were supplied, so no threshold "
    "was applied."
)
intro_lines = (
    "Every table below reproduces an already-computed, already-disclosed "
    "descriptive result table shown elsewhere in this application; this "
    "page performs no new calculation.",
    grouping_line,
    de_threshold_line,
    "This report does not include the Gene Expression page's per-gene "
    "lookup or its multi-gene/time-series charts: those depend on a gene "
    "and chart choice made within that page's own session and are not "
    "part of this dataset-level snapshot.",
    "This report contains no chart images. No normalisation, "
    "transformation, filtering, imputation, clustering, or "
    "differential-expression modelling is performed by this application.",
)

st.header("Report preview")
st.write(
    "This is the exact content the PDF below will contain, laid out for "
    "on-screen review before download."
)
for line in intro_lines:
    st.caption(line)
for section in sections:
    st.subheader(section.title)
    if section.table is not None:
        st.dataframe(section.table, hide_index=True, width="stretch")
    else:
        st.info(section.note)

try:
    artifact = build_report_pdf(
        title=_REPORT_TITLE,
        generated_at=datetime.now(timezone.utc),
        intro_lines=intro_lines,
        sections=sections,
        filename=_REPORT_FILENAME,
    )
except ReportExportError as error:
    st.error(
        "The PDF report could not be assembled "
        f"({error.reason.value}): {error} No download button was shown."
    )
else:
    if current.source == "demo":
        st.warning(
            "The demo values are synthetic, its gene IDs are fictional, and "
            "its p-values are constructed. This report does not support real "
            "tomato nitrate-response conclusions."
        )

    with st.expander("Scientific and statistical limitations"):
        st.info(
            "This report is a formatting convenience over already-shown "
            "descriptive tables. It introduces no new statistic, fits no "
            "model, and does not establish biological validity, replicate "
            "adequacy, or differential-expression significance."
        )
        st.info(
            "A section is omitted with an explicit note, not silently "
            "dropped, whenever its underlying computation is not currently "
            "possible for the active dataset (for example, too few samples "
            "for PCA, or no differential-expression results supplied)."
        )

    st.header("Download report")
    st.download_button(
        "Download descriptive report (PDF)",
        data=artifact.data,
        file_name=artifact.filename,
        mime=PDF_MEDIA_TYPE,
        on_click="ignore",
    )
    st.caption(
        "Regenerating the report from the same active dataset, group-by "
        "choice, and differential-expression thresholds produces the same "
        "table contents; only the 'Generated' timestamp changes."
    )
