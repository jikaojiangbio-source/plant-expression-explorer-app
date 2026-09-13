"""Descriptive Pearson sample-to-sample correlation summaries."""

import pandas as pd
import streamlit as st

from plant_expression_explorer.correlation import (
    CorrelationComputationError,
    build_correlation_observations,
    build_grouped_condition_correlation_summary,
    build_grouped_pair_summary,
    build_grouped_sample_correlation_summary,
    build_heatmap_data,
    compute_sample_correlation,
)
from plant_expression_explorer.dataset import get_current_dataset
from plant_expression_explorer.exports import CsvExportError, build_csv_export
from plant_expression_explorer.provenance import (
    DatasetProvenance,
    provenance_display_rows,
)
from plant_expression_explorer.qc import list_additional_metadata_columns


def _display_matrix(correlation_matrix: pd.DataFrame) -> pd.DataFrame:
    """Return a rounded display copy with undefined values shown as N/A."""

    display = correlation_matrix.copy(deep=True).round(3).astype("object")
    return display.where(display.notna(), "N/A")


def _display_summary(
    summary: pd.DataFrame,
    correlation_columns: tuple[str, ...],
) -> pd.DataFrame:
    """Return a display copy without changing computational NaN values."""

    display = summary.copy(deep=True)
    for column in correlation_columns:
        values = display[column].round(3).astype("object")
        display[column] = values.where(values.notna(), "N/A")
    return display


def _render_dataset_context(provenance: DatasetProvenance | None) -> None:
    with st.expander("Dataset context (descriptive only)"):
        st.caption(
            "Context is displayed verbatim and is not scientifically verified, "
            "parsed, or used in this calculation."
        )
        for label, value in provenance_display_rows(provenance):
            st.caption(label)
            st.code(value, language=None)


def _render_csv_downloads(
    downloads: tuple[tuple[str, pd.DataFrame, str], ...],
) -> None:
    prepared_downloads = []
    try:
        for label, table, filename in downloads:
            artifact = build_csv_export(table, filename=filename)
            prepared_downloads.append((label, artifact))
    except CsvExportError as error:
        st.error(
            "CSV downloads are unavailable "
            f"({error.reason.value}): {error} No download buttons were shown."
        )
        return
    for label, artifact in prepared_downloads:
        st.download_button(
            label,
            data=artifact.data,
            file_name=artifact.filename,
            mime=artifact.media_type,
            on_click="ignore",
        )


st.title("🔥 Sample Correlation")
st.write(
    "Pearson sample-to-sample correlations are descriptive summaries of the "
    "active expression matrix. Correlation does not establish biological "
    "validity and does not imply causation."
)
st.caption(
    "Unusual relationships may warrant further review, but they are not "
    "automatic exclusion criteria."
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
st.write(
    f"**Expression contents:** {current.gene_count:,} genes and "
    f"{current.sample_count:,} samples."
)
_render_dataset_context(current.provenance)

report = current.validation_report
validation_columns = st.columns(3)
validation_columns[0].metric("Validation Errors", len(report.errors))
validation_columns[1].metric("Validation Warnings", len(report.warnings))
validation_columns[2].metric(
    "Validation Information",
    len(report.information),
)

with st.expander("Method"):
    st.write(
        "Pearson correlation is calculated between every sample pair across all "
        "gene rows using complete observations. Sample order follows the "
        "expression-matrix columns; no similarity-based reordering is performed."
    )

try:
    result = compute_sample_correlation(current.expression, current.metadata)
except CorrelationComputationError as error:
    st.error(
        "Descriptive Pearson correlations could not be calculated "
        f"({error.reason.value}): {error}"
    )
    st.stop()

st.header("Correlation matrix")
st.dataframe(
    _display_matrix(result.correlation_matrix),
    width="stretch",
    height=min(max(240, result.sample_count * 42), 700),
)
st.caption(
    "Rows and columns preserve expression sample order. Values are rounded to "
    "three decimals for display only; N/A denotes undefined Pearson "
    "correlation."
)

st.subheader("Correlation heatmap")
heatmap_data = build_heatmap_data(result)
sample_order = list(result.correlation_matrix.columns)
heatmap_spec = {
    "mark": {
        "type": "rect",
        "stroke": "white",
        "strokeWidth": 0.5,
    },
    "encoding": {
        "x": {
            "field": "column_sample",
            "type": "nominal",
            "sort": sample_order,
            "title": "Sample",
            "axis": {"labelAngle": -45},
        },
        "y": {
            "field": "row_sample",
            "type": "nominal",
            "sort": sample_order,
            "title": "Sample",
        },
        "color": {
            "condition": {
                "test": "datum.defined === true",
                "field": "correlation",
                "type": "quantitative",
                "scale": {
                    "domain": [-1, 1],
                    "scheme": "redblue",
                },
                "legend": {"title": "Pearson r"},
            },
            "value": "#b8b8b8",
        },
        "tooltip": [
            {
                "field": "row_sample",
                "type": "nominal",
                "title": "Row sample",
            },
            {
                "field": "column_sample",
                "type": "nominal",
                "title": "Column sample",
            },
            {
                "field": "correlation_label",
                "type": "nominal",
                "title": "Pearson r",
            },
        ],
    },
}
st.vega_lite_chart(
    heatmap_data,
    spec=heatmap_spec,
    width="stretch",
    height=min(max(320, result.sample_count * 48), 900),
)
st.caption(
    "The colour domain is fixed at -1 to 1. Neutral grey cells are undefined; "
    "the heatmap is not clustered and does not imply statistical significance."
)
st.caption(
    "No dedicated application export workflow is implemented. Streamlit "
    "components may expose framework-provided table or chart actions."
)

st.header("Constant samples and undefined correlations")
if result.constant_samples:
    st.info(
        f"{len(result.constant_samples)} constant sample column(s) produce "
        "undefined Pearson correlations: "
        + ", ".join(result.constant_samples)
        + ". They remain in every summary."
    )
else:
    st.write("No constant sample columns were found.")
st.write(
    f"**Undefined unique non-self sample pairs:** "
    f"{result.undefined_pair_count:,}."
)
st.caption(
    "Constant samples have identical expression values across all genes. "
    "Undefined correlations remain N/A and are never replaced with zero or one."
)

additional_columns = list_additional_metadata_columns(current.metadata)
group_column = "condition"
if additional_columns:
    group_column = st.selectbox(
        "Group summaries by",
        options=("condition", *additional_columns),
        help=(
            "Any column present in the uploaded sample metadata beyond "
            "'sample_id' and 'condition' can relabel the tables below. No "
            "Pearson correlation is recalculated for the new grouping; only "
            "the group label and within/between-group membership change."
        ),
    )
grouped_sample_summary = build_grouped_sample_correlation_summary(
    result, current.metadata, group_column
)
grouped_pair_summary = build_grouped_pair_summary(
    result, current.metadata, group_column
)
grouped_condition_summary = build_grouped_condition_correlation_summary(
    result, current.metadata, group_column
)

st.header("Per-sample correlation summary")
if group_column != "condition":
    st.caption(f"Grouped by metadata column '{group_column}', not 'condition'.")
st.dataframe(
    _display_summary(
        grouped_sample_summary,
        (
            "minimum_correlation",
            "median_correlation",
            "mean_correlation",
            "maximum_correlation",
        ),
    ),
    hide_index=True,
    width="stretch",
)
st.caption(
    "Self-correlation is excluded. Statistics use defined non-self pairs only; "
    "samples are not ranked or classified."
)

st.header("Unique sample-pair summary")
st.dataframe(
    _display_summary(grouped_pair_summary, ("correlation",)),
    hide_index=True,
    width="stretch",
)
st.caption(
    "Each unordered pair appears once in expression-column combination order. "
    "Undefined pairs are retained."
)

st.header("Condition-pair descriptive summary")
st.dataframe(
    _display_summary(
        grouped_condition_summary,
        ("median_correlation",),
    ),
    hide_index=True,
    width="stretch",
)
st.caption(
    "Condition order follows first appearance across expression sample "
    "columns. Medians use defined pairs only; condition labels do not validate "
    "the experimental design, and no comparison test is performed."
)

st.header("Deterministic structural observations")
observations = build_correlation_observations(result, report)
if observations:
    for observation in observations:
        st.markdown(f"- {observation}")
else:
    st.write(
        "No additional exact structural observations were generated by these "
        "descriptive summaries."
    )

if current.source == "demo":
    st.warning(
        "The demo values are synthetic, its gene IDs are fictional, and its "
        "p-values are constructed. It does not define real correlation "
        "thresholds or support tomato biological conclusions."
    )

with st.expander("Scientific and statistical limitations"):
    st.info(
        "Input scale and transformation affect Pearson correlation, and gene "
        "filtering or selection can change every displayed value. Broad expression "
        "distributions may dominate these summaries."
    )
    st.info(
        "High correlation does not prove replicate validity. A low or negative "
        "correlation alone does not prove that a sample is unsuitable. Constant "
        "samples produce undefined Pearson correlations."
    )
    st.info(
        "No samples or genes are modified or removed. This page performs no PCA, "
        "clustering, distance analysis, batch correction, hypothesis testing, "
        "correlation p-value calculation, or differential-expression inference."
    )

st.header("Download descriptive results")
st.write(
    "Downloads are UTF-8 CSV copies of the underlying unrounded result tables. "
    "They preserve result row and column order, omit pandas index and dtype "
    "metadata, and represent undefined correlations as empty fields."
)
st.warning(
    "CSV text is not prefixed or rewritten. Spreadsheet software may interpret "
    "formula-like leading characters in untrusted text; review such data and "
    "import it as plain text when needed."
)
matrix_long = heatmap_data.loc[
    :, ["row_sample", "column_sample", "correlation", "defined"]
].copy(deep=True)
_group_file_suffix = "" if group_column == "condition" else f"-by-{group_column}"
_render_csv_downloads(
    (
        (
            "Download correlation matrix in long form (CSV)",
            matrix_long,
            "sample-correlation-matrix-long.csv",
        ),
        (
            "Download per-sample correlation summary (CSV)",
            grouped_sample_summary,
            f"sample-correlation-sample-summary{_group_file_suffix}.csv",
        ),
        (
            "Download unique sample pairs (CSV)",
            grouped_pair_summary,
            f"sample-correlation-unique-pairs{_group_file_suffix}.csv",
        ),
        (
            "Download condition-pair summary (CSV)",
            grouped_condition_summary,
            f"sample-correlation-condition-pairs{_group_file_suffix}.csv",
        ),
    )
)
st.caption(
    "The long-form matrix includes every ordered row/column sample position and "
    "a defined flag. These exports do not classify correlation strength, sample "
    "quality, or biological validity."
)
