"""Descriptive sample quality-control summaries for the active dataset."""

import pandas as pd
import streamlit as st

from plant_expression_explorer.dataset import get_current_dataset
from plant_expression_explorer.exports import CsvExportError, build_csv_export
from plant_expression_explorer.qc import (
    QcComputationError,
    build_condition_chart_data,
    build_qc_observations,
    build_sample_chart_data,
    compute_sample_qc,
    should_display_count_chart,
)


st.html(
    """
    <style>
    button[aria-label="Download as CSV"],
    button[aria-label="Download as PNG"],
    button[aria-label="Copy Vega-Lite spec"] {
        display: none;
    }
    </style>
    """
)


def _display_sample_summary(sample_summary: pd.DataFrame) -> pd.DataFrame:
    """Return a display copy that shows undefined standard deviation as N/A."""

    display = sample_summary.copy(deep=True)
    standard_deviation = display["standard_deviation"].astype("object")
    display["standard_deviation"] = standard_deviation.where(
        standard_deviation.notna(),
        "N/A",
    )
    return display


def _render_csv_downloads(
    downloads: tuple[
        tuple[str, pd.DataFrame, str, tuple[str, ...]],
        ...,
    ],
) -> None:
    prepared_downloads = []
    try:
        for label, table, filename, json_sequence_columns in downloads:
            artifact = build_csv_export(
                table,
                filename=filename,
                json_sequence_columns=json_sequence_columns,
            )
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


st.title("🧪 Sample Quality Control")
st.write(
    "These descriptive summaries show the structure and value distributions "
    "of the active expression matrix. They do not establish biological validity."
)
st.caption(
    "This page does not replace FASTQ, mapping, library-quality, experimental-"
    "design, or upstream processing assessment."
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
    f"**Bundle contents:** {current.gene_count:,} genes, "
    f"{current.sample_count:,} samples, and "
    f"{current.de_row_count:,} differential-expression rows."
)

report = current.validation_report
validation_columns = st.columns(3)
validation_columns[0].metric("Validation Errors", len(report.errors))
validation_columns[1].metric("Validation Warnings", len(report.warnings))
validation_columns[2].metric(
    "Validation Information",
    len(report.information),
)

try:
    result = compute_sample_qc(current.expression, current.metadata)
except QcComputationError as error:
    st.error(
        "Descriptive QC summaries could not be calculated "
        f"({error.reason.value}): {error}"
    )
    st.stop()

st.header("Matrix overview")
overview_rows = (
    (
        ("Genes", result.gene_count),
        ("Samples", result.sample_count),
        ("Expression cells", result.expression_cell_count),
    ),
    (
        ("Missing values", result.missing_value_count),
        ("Non-finite values", result.non_finite_value_count),
        ("Conditions", result.condition_count),
    ),
    (
        ("Zero values", result.zero_value_count),
        ("Negative values", result.negative_value_count),
        ("Constant samples", len(result.constant_samples)),
    ),
)
for row in overview_rows:
    columns = st.columns(3)
    for column, (label, value) in zip(columns, row, strict=True):
        column.metric(label, value)
st.write(
    "**Metadata order relative to expression columns:** "
    + (
        "corresponds exactly."
        if result.metadata_order_matches_expression
        else "differs; conditions were mapped by exact sample ID."
    )
)

st.header("Condition and replicate summary")
st.dataframe(
    result.condition_summary,
    hide_index=True,
    width="stretch",
)
st.caption(
    "Counts and sample membership are descriptive. Replicate counts do not "
    "determine biological validity."
)

st.header("Per-sample statistics")
st.dataframe(
    _display_sample_summary(result.sample_summary),
    hide_index=True,
    width="stretch",
)
st.caption(
    "Quartiles use pandas linear interpolation. Standard deviation uses "
    "ddof=1; it is shown as N/A when undefined."
)

st.header("Descriptive charts")
st.subheader("Median expression by sample")
st.bar_chart(
    build_sample_chart_data(result, "median"),
    x="sample_id",
    y="median",
    x_label="Sample",
    y_label="Median expression",
    sort=False,
)

st.subheader("Interquartile range by sample")
st.bar_chart(
    build_sample_chart_data(result, "interquartile_range"),
    x="sample_id",
    y="interquartile_range",
    x_label="Sample",
    y_label="Interquartile range",
    sort=False,
)

st.subheader("Zero-value count by sample")
if should_display_count_chart(result, "zero_value_count"):
    st.bar_chart(
        build_sample_chart_data(result, "zero_value_count"),
        x="sample_id",
        y="zero_value_count",
        x_label="Sample",
        y_label="Zero-value count",
        sort=False,
    )
else:
    st.info("Every sample has a zero-value count of 0; the chart is omitted.")

st.subheader("Negative-value count by sample")
if should_display_count_chart(result, "negative_value_count"):
    st.bar_chart(
        build_sample_chart_data(result, "negative_value_count"),
        x="sample_id",
        y="negative_value_count",
        x_label="Sample",
        y_label="Negative-value count",
        sort=False,
    )
else:
    st.info(
        "Every sample has a negative-value count of 0; the chart is omitted."
    )

st.subheader("Sample count by condition")
st.bar_chart(
    build_condition_chart_data(result),
    x="condition",
    y="sample_count",
    x_label="Condition",
    y_label="Sample count",
    sort=False,
)

st.header("Constant-sample and zero-variance findings")
if result.gene_count == 1:
    st.info(
        "With one gene, every sample column meets the exact constant-column "
        "definition, so this diagnostic is not informative."
    )
elif result.constant_samples:
    st.info(
        f"{len(result.constant_samples)} constant sample column(s): "
        + ", ".join(result.constant_samples)
        + ". They were retained."
    )
else:
    st.write("No constant sample columns were found.")
st.write(
    f"**Zero-variance genes:** {result.zero_variance_gene_count:,}. "
    "These genes were retained."
)
st.caption(
    "Constant samples and zero-variance genes use exact equality. Such findings "
    "may affect later scale-dependent analyses but are not automatic exclusion "
    "criteria."
)

st.header("Observations that may warrant further review")
observations = build_qc_observations(result, report)
if observations:
    for observation in observations:
        st.markdown(f"- {observation}")
else:
    st.write(
        "No additional exact structural observations were generated by these "
        "descriptive checks."
    )

if current.source == "demo":
    st.warning(
        "The demo values are synthetic, its gene IDs are fictional, and its "
        "p-values are constructed. It does not establish real QC thresholds "
        "or support tomato nitrate-response conclusions."
    )

st.header("Scientific and statistical limitations")
st.info(
    "Different raw, normalised, transformed, and centred expression scales "
    "require different interpretations. Zero and negative values therefore "
    "have scale-dependent meanings. Sample medians, spreads, and unusual "
    "descriptive metrics do not establish biological validity and are not "
    "automatic exclusion criteria."
)
st.info(
    "No samples or genes are modified or removed. This page performs no "
    "normalisation, PCA, correlation, clustering, hypothesis testing, or "
    "differential-expression inference."
)

st.header("Download descriptive results")
st.write(
    "Downloads are UTF-8 CSV copies of the underlying unrounded result tables. "
    "They preserve result row and column order, do not include the pandas index "
    "or dtype metadata, and represent missing values as empty fields."
)
st.warning(
    "CSV text is not prefixed or rewritten. Spreadsheet software may interpret "
    "formula-like leading characters in untrusted text; review such data and "
    "import it as plain text when needed."
)
_render_csv_downloads(
    (
        (
            "Download condition and sample membership (CSV)",
            result.condition_summary,
            "sample-qc-condition-membership.csv",
            ("sample_ids",),
        ),
        (
            "Download per-sample statistics (CSV)",
            result.sample_summary,
            "sample-qc-sample-statistics.csv",
            (),
        ),
    )
)
st.caption(
    "The sample_ids field in the condition-membership file is a JSON array so "
    "every exact sample identifier remains distinguishable. These files contain "
    "descriptive summaries, not quality classifications or exclusion decisions."
)
