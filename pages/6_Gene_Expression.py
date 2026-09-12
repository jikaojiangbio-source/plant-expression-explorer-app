"""Descriptive expression lookup for one exact supplied gene identifier."""

import pandas as pd
import streamlit as st

from plant_expression_explorer.dataset import get_current_dataset
from plant_expression_explorer.exports import CsvExportError, build_csv_export
from plant_expression_explorer.gene_expression import (
    GeneExpressionComputationError,
    build_gene_expression_chart_data,
    build_gene_expression_observations,
    list_gene_ids,
    lookup_gene_expression,
)
from plant_expression_explorer.validation import Severity, ValidationIssue


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


def _is_expression_metadata_issue(issue: ValidationIssue) -> bool:
    relevant_tables = {"Expression matrix", "Sample metadata"}
    issue_tables = {issue.table}
    if issue.related_table is not None:
        issue_tables.add(issue.related_table)
    return issue_tables.issubset(relevant_tables)


def _display_condition_summary(summary: pd.DataFrame) -> pd.DataFrame:
    """Return a display copy with undefined sample SD shown as N/A."""

    display = summary.copy(deep=True)
    display["sample_ids"] = display["sample_ids"].map(", ".join)
    display["standard_deviation"] = display["standard_deviation"].map(
        lambda value: "N/A" if pd.isna(value) else repr(float(value))
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


st.title("🌿 Gene Expression")
st.write(
    "This page displays the supplied preprocessed expression values for one "
    "exact `gene_id` across samples. It does not normalize, transform, impute, "
    "or model the data, and it does not test condition effects or infer "
    "differential expression, statistical significance, contrast direction, "
    "reference levels, or biological importance."
)
st.caption(
    "Results remain on the scale supplied by the user. Any normalization, "
    "transformation, or filtering was performed upstream and is not inferred "
    "or changed here."
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

report = current.validation_report
validation_columns = st.columns(3)
validation_columns[0].metric("Validation Errors", len(report.errors))
validation_columns[1].metric("Validation Warnings", len(report.warnings))
validation_columns[2].metric(
    "Validation Information",
    len(report.information),
)

if report.has_errors:
    st.error(
        "Gene-expression lookup is unavailable because the active dataset "
        "contains blocking validation Errors. Return to Upload Data and "
        "correct the source tables."
    )
    for issue in report.errors:
        _render_issue(issue)
    st.stop()

st.header("Relevant validation observations")
relevant_issues = tuple(
    issue for issue in report.issues if _is_expression_metadata_issue(issue)
)
if relevant_issues:
    for issue in relevant_issues:
        _render_issue(issue)
else:
    st.write(
        "No expression-matrix or sample-metadata validation issues were "
        "reported for the active dataset."
    )

if current.source == "demo":
    st.warning(
        "These expression values are synthetic and the gene IDs are fictional. "
        "The display demonstrates software behaviour only and does not support "
        "conclusions about tomato biology or nitrate response."
    )

try:
    gene_ids = list_gene_ids(current.expression)
except GeneExpressionComputationError as error:
    st.error(
        "Gene options could not be built "
        f"({error.reason.value}): {error}"
    )
    st.stop()

st.header("Select one supplied gene")
st.write(
    "Options follow expression-matrix row order. Matching uses the exact "
    "validated `str(value)` representation of each supplied identifier: no "
    "identifier is trimmed, case-folded, normalized, or searched by alias, and "
    "the source values remain unchanged."
)
selected_gene_id = st.selectbox(
    "Exact gene ID",
    gene_ids,
    index=None,
    placeholder="Select a supplied gene ID",
    help="Typing filters the supplied options; it does not create a new identifier.",
)
if selected_gene_id is None:
    st.info("Select one exact supplied gene ID to display its expression values.")
    st.stop()

try:
    result = lookup_gene_expression(
        current.expression,
        current.metadata,
        selected_gene_id,
    )
except GeneExpressionComputationError as error:
    st.error(
        "Gene-expression values could not be displayed "
        f"({error.reason.value}): {error}"
    )
    st.stop()

st.header(f"Expression values for {result.gene_id}")
summary_columns = st.columns(3)
summary_columns[0].metric("Samples displayed", result.sample_count)
summary_columns[1].metric("Condition labels", result.condition_count)
summary_columns[2].metric(
    "All values exactly equal",
    "Yes" if result.all_values_equal else "No",
)
st.write(
    "**Metadata order relative to expression columns:** "
    + (
        "corresponds exactly."
        if result.metadata_order_matches_expression
        else "differs; conditions were mapped by exact sample ID."
    )
)

st.subheader("Supplied per-sample values")
st.dataframe(
    result.sample_expression,
    hide_index=True,
    width="stretch",
)
st.caption(
    "Every expression sample column appears exactly once and in source column "
    "order. Expression values are copied from the selected source row without "
    "rounding or replacement; neither source table is reordered or rewritten."
)

st.subheader("Per-sample expression plot")
if result.sample_count >= 2:
    chart_data = build_gene_expression_chart_data(result)
    sample_order = result.sample_expression["sample_id"].tolist()
    condition_order = list(dict.fromkeys(result.sample_expression["condition"]))
    point_spec = {
        "mark": {"type": "point", "filled": True, "size": 110},
        "encoding": {
            "x": {
                "field": "sample_id",
                "type": "nominal",
                "sort": sample_order,
                "title": "Sample",
                "axis": {"labelAngle": -45},
            },
            "y": {
                "field": "expression_value",
                "type": "quantitative",
                "title": "Supplied preprocessed expression value",
                "scale": {"zero": True},
            },
            "color": {
                "field": "condition",
                "type": "nominal",
                "sort": condition_order,
                "legend": {"title": "Condition"},
            },
            "order": {
                "field": "sample_position",
                "type": "quantitative",
            },
            "tooltip": [
                {"field": "sample_id", "type": "nominal", "title": "Sample"},
                {
                    "field": "condition",
                    "type": "nominal",
                    "title": "Condition",
                },
                {
                    "field": "expression_value",
                    "type": "quantitative",
                    "title": "Expression value",
                },
            ],
        },
    }
    st.vega_lite_chart(
        chart_data,
        spec=point_spec,
        width="stretch",
        height=420,
    )
    st.caption(
        "Each point represents one supplied sample. Colour shows the exact "
        "metadata condition label for context only. Points are not connected, "
        "the y-axis includes zero on the supplied scale, and no group estimate "
        "or statistical comparison is shown."
    )
else:
    st.info(
        "Only one sample is available, so an across-sample expression plot is "
        "not shown. The supplied value remains visible in the table."
    )

st.header("Condition-grouped descriptive summary")
st.dataframe(
    _display_condition_summary(result.condition_summary),
    hide_index=True,
    width="stretch",
)
st.caption(
    "These are arithmetic summaries of the displayed sample values on the "
    "supplied scale. Samples are grouped only by exact condition label. "
    "Standard deviation uses ddof=1 and is N/A for a one-sample group. These "
    "summaries do not establish biological replication, a condition effect, or "
    "statistical significance."
)

st.header("Descriptive structural observations")
observations = build_gene_expression_observations(result, report)
if observations:
    for observation in observations:
        st.markdown(f"- {observation}")
else:
    st.write(
        "No additional exact structural observations were generated for this "
        "selected gene."
    )

st.header("Scientific and statistical limitations")
st.info(
    "Condition labels provide display context only. This page does not infer "
    "a reference level, contrast direction, regulation, a condition effect, "
    "statistical significance, or biological importance."
)
st.info(
    "No genes, samples, identifiers, conditions, or expression values are "
    "trimmed, normalized, transformed, deduplicated, aggregated into the source "
    "data, imputed, intersected, reordered, or discarded."
)
st.info(
    "This page performs no FASTQ processing, batch correction, hypothesis "
    "testing, differential-expression inference, or volcano plotting."
)

st.header("Download descriptive results")
st.write(
    "Downloads are UTF-8 CSV copies of the underlying unrounded result tables. "
    "The selected exact gene ID is added as an explicit context column. Result "
    "row and column order is preserved; pandas index and dtype metadata are not "
    "included, and missing values are empty fields."
)
st.warning(
    "CSV text is not prefixed or rewritten. Spreadsheet software may interpret "
    "formula-like leading characters in untrusted text; review such data and "
    "import it as plain text when needed."
)
sample_export = result.sample_expression.copy(deep=True)
sample_export.insert(0, "gene_id", result.gene_id)
condition_export = result.condition_summary.copy(deep=True)
condition_export.insert(0, "gene_id", result.gene_id)
_render_csv_downloads(
    (
        (
            "Download per-sample expression values (CSV)",
            sample_export,
            "gene-expression-sample-values.csv",
            (),
        ),
        (
            "Download condition-grouped summary (CSV)",
            condition_export,
            "gene-expression-condition-summary.csv",
            ("sample_ids",),
        ),
    )
)
st.caption(
    "The sample_ids field in the condition summary is a JSON array. These files "
    "contain supplied values and disclosed descriptive aggregations only; they "
    "do not report condition effects, significance, regulation, or importance."
)
