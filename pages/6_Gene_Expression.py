"""Descriptive expression lookup for one exact supplied gene identifier."""

import pandas as pd
import plotly.express as px
import streamlit as st

from plant_expression_explorer.annotations import (
    identifier_format_hint,
    list_supported_species,
    lookup_gene_annotation,
)
from plant_expression_explorer.dataset import (
    ACTIVE_GROUP_COLUMN_KEY,
    ensure_valid_group_column_state,
    get_current_dataset,
)
from plant_expression_explorer.exports import CsvExportError, build_csv_export
from plant_expression_explorer.gene_expression import (
    GeneExpressionComputationError,
    build_gene_expression_chart_data,
    build_gene_expression_observations,
    build_grouped_gene_expression_chart_data,
    build_grouped_gene_expression_condition_summary,
    build_multi_gene_panel_data,
    build_time_series_chart_data,
    filter_gene_ids,
    list_gene_ids,
    lookup_gene_expression,
)
from plant_expression_explorer.provenance import (
    DatasetProvenance,
    provenance_display_rows,
)
from plant_expression_explorer.qc import list_additional_metadata_columns
from plant_expression_explorer.theme import inject_global_styles
from plant_expression_explorer.validation import Severity, ValidationIssue

_SEARCH_BOX_GENE_COUNT_THRESHOLD = 200
_SEARCH_RESULT_DISPLAY_LIMIT = 500


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


def _render_dataset_context(provenance: DatasetProvenance | None) -> None:
    with st.expander("Dataset context (descriptive only)"):
        st.caption(
            "Context is displayed verbatim and is not scientifically verified, "
            "parsed, or used in this calculation."
        )
        for label, value in provenance_display_rows(provenance):
            st.caption(label)
            st.code(value, language=None)


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


inject_global_styles()
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
_render_dataset_context(current.provenance)

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

selectable_gene_ids = gene_ids
if len(gene_ids) > _SEARCH_BOX_GENE_COUNT_THRESHOLD:
    search_query = st.text_input(
        "Search gene IDs",
        placeholder="Type part of a gene ID to narrow the list below",
        help=(
            "Case-insensitive substring match against the exact supplied gene "
            "IDs. This dataset has "
            f"{len(gene_ids):,} genes; searching keeps the selector responsive."
        ),
    )
    selectable_gene_ids = filter_gene_ids(gene_ids, search_query)
    if len(selectable_gene_ids) > _SEARCH_RESULT_DISPLAY_LIMIT:
        st.info(
            f"{len(selectable_gene_ids):,} gene IDs match; showing the first "
            f"{_SEARCH_RESULT_DISPLAY_LIMIT:,} in matrix row order. Narrow your "
            "search to see others."
        )
        selectable_gene_ids = selectable_gene_ids[:_SEARCH_RESULT_DISPLAY_LIMIT]
    elif not selectable_gene_ids:
        st.info("No supplied gene ID contains that text.")

selected_gene_id = st.selectbox(
    "Exact gene ID",
    selectable_gene_ids,
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

st.header("Multi-gene panel (optional)")
st.write(
    "Optionally compare the gene selected above with additional exact gene "
    "IDs on one combined chart and table. Each gene is looked up "
    "independently and shown on its own supplied scale; no normalization or "
    "cross-gene scaling is applied, so different genes' absolute magnitudes "
    "are not directly comparable."
)
additional_gene_ids = st.multiselect(
    "Compare with additional gene IDs",
    [gene_id for gene_id in selectable_gene_ids if gene_id != result.gene_id],
    default=[],
    help=(
        "Selecting one or more genes here adds a combined panel below; it "
        "does not replace the detailed single-gene view further down."
    ),
)
if additional_gene_ids:
    panel_gene_ids = [result.gene_id, *additional_gene_ids]
    try:
        panel_data = build_multi_gene_panel_data(
            current.expression, current.metadata, panel_gene_ids
        )
    except GeneExpressionComputationError as error:
        st.error(
            "The multi-gene panel could not be displayed "
            f"({error.reason.value}): {error}"
        )
    else:
        panel_sample_order = result.sample_expression["sample_id"].tolist()
        panel_gene_order = list(dict.fromkeys(panel_data["gene_id"]))
        panel_figure = px.line(
            panel_data,
            x="sample_id",
            y="expression_value",
            color="gene_id",
            markers=True,
            category_orders={
                "sample_id": panel_sample_order,
                "gene_id": panel_gene_order,
            },
            labels={
                "sample_id": "Sample",
                "expression_value": "Supplied preprocessed expression value",
                "gene_id": "Gene ID",
            },
        )
        panel_figure.update_traces(marker=dict(size=9), line=dict(width=1.5))
        panel_figure.update_layout(
            height=440,
            margin=dict(l=10, r=10, t=10, b=10),
            legend_title_text="Gene ID",
            xaxis=dict(tickangle=-45),
        )
        st.plotly_chart(panel_figure, width="stretch")
        st.caption(
            "Lines connect one gene's own points across samples in "
            "expression-column order only, to make each gene's series "
            "easier to trace; they are not a fitted trend or a claim about "
            "intermediate values. Colour identifies the gene, not a "
            "condition or group."
        )

        panel_wide = panel_data.pivot(
            index=["sample_id", "condition"],
            columns="gene_id",
            values="expression_value",
        ).reset_index()
        panel_wide.columns.name = None
        st.dataframe(panel_wide, hide_index=True, width="stretch")
        st.caption(
            "One row per sample; one column per selected gene, each holding "
            "that gene's exact supplied value for that sample."
        )

        panel_export = panel_data.copy(deep=True)
        _render_csv_downloads(
            (
                (
                    "Download multi-gene panel values (CSV)",
                    panel_export,
                    "gene-expression-multi-gene-panel.csv",
                    (),
                ),
            )
        )

st.header(f"Expression values for {result.gene_id}")

species_label = st.selectbox(
    "Species reference (optional)",
    options=("Not selected", *list_supported_species()),
    help=(
        "Matches the exact selected gene ID against a small, hand-curated "
        "list of well-known reference genes for one species (see "
        "data/annotations/README.md). This is not a genome annotation, is "
        "not used in any calculation, and most real gene IDs will not "
        "have an entry."
    ),
)
if species_label != "Not selected":
    annotation = lookup_gene_annotation(species_label, result.gene_id)
    if annotation is None:
        st.caption(
            f"No entry for '{result.gene_id}' in the small curated reference "
            f"list for {species_label}. This is expected for most gene IDs; "
            "it does not indicate a problem with the gene ID."
        )
        format_hint = identifier_format_hint(species_label, result.gene_id)
        if format_hint is not None:
            st.caption(
                f"⚠️ {format_hint} This is a descriptive note about "
                "identifier shape, not a lookup: the gene ID is never "
                "rewritten or searched under any other form."
            )
    else:
        st.info(
            f"**{annotation.symbol}** — {annotation.description}  \n"
            f"Source: {annotation.source}."
        )

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

additional_columns = list_additional_metadata_columns(current.metadata)
group_column = "condition"
if additional_columns:
    ensure_valid_group_column_state(st.session_state, additional_columns)
    group_column = st.selectbox(
        "Group plot and summary by",
        options=("condition", *additional_columns),
        key=ACTIVE_GROUP_COLUMN_KEY,
        help=(
            "Any column present in the uploaded sample metadata beyond "
            "'sample_id' and 'condition' can relabel the plot and summary "
            "below. The same descriptive statistics are recomputed for the "
            "new grouping; the underlying per-sample values are unchanged. "
            "This choice is shared with the PCA, Sample Quality Control, and "
            "Sample Correlation pages."
        ),
    )
group_title = group_column.replace("_", " ").capitalize()

time_column = None
if additional_columns:
    axis_options = (
        "Sample (upload order)",
        *[f"Numeric time/order: {col}" for col in additional_columns],
    )
    axis_choice = st.selectbox(
        "Chart x-axis",
        options=axis_options,
        help=(
            "'Sample (upload order)' plots samples as categories in "
            "expression-column order, coloured by the grouping above. A "
            "numeric option instead connects points in ascending order of "
            "that metadata column's exact supplied value; every sample "
            "must have a numeric value in that column, or the chart shows "
            "a controlled error instead of skipping samples."
        ),
    )
    if axis_choice != "Sample (upload order)":
        time_column = axis_choice.removeprefix("Numeric time/order: ")

st.subheader("Per-sample expression plot")
if result.sample_count < 2:
    st.info(
        "Only one sample is available, so an across-sample expression plot is "
        "not shown. The supplied value remains visible in the table."
    )
elif time_column is not None:
    try:
        time_series = build_time_series_chart_data(
            result, current.metadata, time_column
        )
    except ValueError as error:
        st.error(
            f"'{time_column}' cannot be used as a numeric time axis: {error}"
        )
    else:
        time_figure = px.line(
            time_series,
            x="time_value",
            y="expression_value",
            markers=True,
            labels={
                "time_value": time_column,
                "expression_value": "Supplied preprocessed expression value",
            },
            hover_data={"sample_id": True, "condition": True},
        )
        time_figure.update_traces(marker=dict(size=10), line=dict(width=1.5))
        time_figure.update_layout(
            height=420,
            margin=dict(l=10, r=10, t=10, b=10),
        )
        st.plotly_chart(time_figure, width="stretch")
        st.caption(
            f"Points are ordered by the exact supplied '{time_column}' "
            "value and connected by a line only to make the sequence "
            "easier to trace; this is not a fitted trend, interpolation, "
            "or a claim about values between samples. Samples sharing the "
            "same time value (for example replicates at one timepoint) "
            "are shown individually, not averaged."
        )
if result.sample_count >= 2 and time_column is None:
    chart_data = build_grouped_gene_expression_chart_data(
        result, current.metadata, group_column
    )
    sample_order = result.sample_expression["sample_id"].tolist()
    group_order = list(dict.fromkeys(chart_data[group_column]))
    point_figure = px.scatter(
        chart_data,
        x="sample_id",
        y="expression_value",
        color=group_column,
        category_orders={"sample_id": sample_order, group_column: group_order},
        labels={
            "sample_id": "Sample",
            "expression_value": "Supplied preprocessed expression value",
            group_column: group_title,
        },
    )
    point_figure.update_traces(marker=dict(size=12, symbol="circle", line=dict(width=0)))
    point_figure.update_layout(
        height=420,
        margin=dict(l=10, r=10, t=10, b=10),
        legend_title_text=group_title,
        xaxis=dict(tickangle=-45),
        yaxis=dict(rangemode="tozero"),
    )
    st.plotly_chart(point_figure, width="stretch")
    st.caption(
        "Each point represents one supplied sample. Colour shows the selected "
        "metadata column for context only. Points are not connected, the "
        "y-axis includes zero on the supplied scale, and no group estimate or "
        "statistical comparison is shown. Zoom, pan, and hover are Plotly's "
        "built-in interactions and do not change the underlying values."
    )

st.header("Condition-grouped descriptive summary")
if group_column != "condition":
    st.caption(f"Grouped by metadata column '{group_column}', not 'condition'.")
grouped_condition_summary = build_grouped_gene_expression_condition_summary(
    result, current.metadata, group_column
)
st.dataframe(
    _display_condition_summary(grouped_condition_summary),
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

with st.expander("Scientific and statistical limitations"):
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
_group_file_suffix = "" if group_column == "condition" else f"-by-{group_column}"
sample_export = result.sample_expression.copy(deep=True)
sample_export.insert(0, "gene_id", result.gene_id)
condition_export = grouped_condition_summary.copy(deep=True)
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
            f"gene-expression-condition-summary{_group_file_suffix}.csv",
            ("sample_ids",),
        ),
    )
)
st.caption(
    "The sample_ids field in the condition summary is a JSON array. These files "
    "contain supplied values and disclosed descriptive aggregations only; they "
    "do not report condition effects, significance, regulation, or importance."
)
