"""Descriptive sample principal component analysis (PCA)."""

import pandas as pd
import plotly.express as px
import streamlit as st

from plant_expression_explorer.dataset import (
    ACTIVE_GROUP_COLUMN_KEY,
    ensure_valid_group_column_state,
    get_current_dataset,
)
from plant_expression_explorer.exports import CsvExportError, build_csv_export
from plant_expression_explorer.pca import (
    PcaComputationError,
    build_grouped_score_plot_data,
    build_pca_observations,
    build_score_plot_data,
    compute_sample_pca,
)
from plant_expression_explorer.provenance import (
    DatasetProvenance,
    provenance_display_rows,
)
from plant_expression_explorer.qc import list_additional_metadata_columns
from plant_expression_explorer.theme import inject_global_styles


def _display_variance_table(variance_table: pd.DataFrame) -> pd.DataFrame:
    """Return a rounded display copy; computational values are unchanged."""

    display = variance_table.copy(deep=True)
    display["component"] = display["component"].map(lambda value: f"PC{value}")
    display["singular_value"] = display["singular_value"].round(3)
    display["explained_variance"] = display["explained_variance"].round(3)
    display["explained_variance_ratio"] = display["explained_variance_ratio"].round(4)
    return display


def _display_score_table(score_table: pd.DataFrame) -> pd.DataFrame:
    """Return a rounded display copy; computational values are unchanged."""

    display = score_table.copy(deep=True)
    for column in display.columns:
        if column.startswith("pc"):
            display[column] = display[column].round(3)
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


inject_global_styles()
st.title("📊 PCA")
st.write(
    "Principal component analysis (PCA) is a descriptive, unsupervised summary "
    "of the active expression matrix. It does not use sample conditions to fit "
    "the components; condition labels are shown only as visual colour coding. "
    "Proximity in a PCA plot does not prove biological similarity or replicate "
    "validity, and separation between groups does not prove a condition effect "
    "or establish statistical significance. PCA does not identify failed or "
    "low-quality samples."
)
st.caption(
    "Results depend on whatever normalization, transformation, filtering, and "
    "gene selection were applied upstream."
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
        "Each gene is mean-centred across samples in a temporary computation copy "
        "before singular value decomposition; the active expression matrix is "
        "never rewritten. Genes are not scaled to unit variance. This preserves "
        "the relative variance structure of the supplied preprocessed matrix and "
        "avoids adding a separate standardization step. Genes with larger "
        "variance in the supplied matrix consequently contribute more strongly to "
        "the components; this is a disclosed analysis policy, not a claim that "
        "such genes are more biologically informative, and not a claim that this "
        "is a universally superior PCA method."
    )
    st.caption(
        "Finite expression values that cannot be safely represented through the "
        "float64 PCA calculation (for example, extremely large or extremely "
        "small magnitudes) produce a controlled error; they are never "
        "automatically rescaled or rewritten."
    )

try:
    with st.spinner("Computing PCA…"):
        result = compute_sample_pca(current.expression, current.metadata)
except PcaComputationError as error:
    st.error(
        "Descriptive PCA could not be calculated "
        f"({error.reason.value}): {error}"
    )
    st.stop()

st.header("Dataset summary for PCA")
summary_columns = st.columns(4)
summary_columns[0].metric("Samples", result.sample_count)
summary_columns[1].metric("Genes", result.gene_count)
summary_columns[2].metric("Components", result.component_count)
summary_columns[3].metric(
    "Constant genes",
    result.zero_variance_gene_count,
    help=(
        "Genes with zero variance across samples. They are retained in the "
        "matrix, not removed, and contribute no signal to the PCA."
    ),
)
st.write(
    "**Metadata order relative to expression columns:** "
    + (
        "corresponds exactly."
        if result.metadata_order_matches_expression
        else "differs; conditions were mapped by exact sample ID."
    )
)

st.header("Explained variance")
st.dataframe(
    _display_variance_table(result.variance_table),
    hide_index=True,
    width="stretch",
)
st.caption(
    "Singular values at or below a numerical-noise tolerance are reported as "
    "exactly 0, along with their explained variance; this is treatment of "
    "derived decomposition output, not a change to the source or centred "
    "expression matrix. Ratios sum to 1.0 within a small floating-point "
    "tolerance, not exactly."
)
variance_chart_data = result.variance_table.assign(
    component_label=result.variance_table["component"].map(
        lambda value: f"PC{value}"
    )
)
st.bar_chart(
    variance_chart_data,
    x="component_label",
    y="explained_variance_ratio",
    x_label="Component",
    y_label="Explained variance ratio",
    sort=False,
)

st.header("Sample scores")
st.dataframe(
    _display_score_table(result.score_table),
    hide_index=True,
    width="stretch",
)
st.caption(
    "Rows preserve expression sample order. Conditions are shown for visual "
    "grouping only and were not used to fit the components."
)

if result.component_count >= 2:
    st.subheader("PC1-versus-PC2 sample plot")
    additional_columns = list_additional_metadata_columns(current.metadata)
    group_column = "condition"
    if additional_columns:
        ensure_valid_group_column_state(st.session_state, additional_columns)
        group_column = st.selectbox(
            "Colour points by",
            options=("condition", *additional_columns),
            key=ACTIVE_GROUP_COLUMN_KEY,
            help=(
                "Any column present in the uploaded sample metadata beyond "
                "'sample_id' and 'condition' can be used for visual grouping "
                "only; the choice never changes how components were fit. This "
                "choice is shared with the Sample Quality Control, Sample "
                "Correlation, and Gene Expression pages."
            ),
        )
    plot_data = build_grouped_score_plot_data(result, current.metadata, group_column)
    ratio_by_component = result.variance_table.set_index("component")[
        "explained_variance_ratio"
    ]
    pc1_ratio = float(ratio_by_component.loc[1])
    pc2_ratio = float(ratio_by_component.loc[2])
    group_order = list(dict.fromkeys(plot_data[group_column]))
    group_title = group_column.replace("_", " ").capitalize()
    scatter_figure = px.scatter(
        plot_data,
        x="pc1",
        y="pc2",
        color=group_column,
        category_orders={group_column: group_order},
        hover_name="sample_id",
        hover_data={"pc1": ":.4f", "pc2": ":.4f", group_column: True},
        labels={
            "pc1": f"PC1 ({pc1_ratio:.1%} explained variance)",
            "pc2": f"PC2 ({pc2_ratio:.1%} explained variance)",
            group_column: group_title,
        },
    )
    scatter_figure.update_traces(marker=dict(size=12, line=dict(width=0)))
    scatter_figure.update_layout(
        height=420,
        margin=dict(l=10, r=10, t=10, b=10),
        legend_title_text=group_title,
    )
    st.plotly_chart(scatter_figure, width="stretch")
    st.caption(
        "Colour reflects the selected metadata column for visual grouping "
        "only; it was not used to fit the components. Proximity does not "
        "prove biological similarity, and separation does not prove a group "
        "effect. Zoom, pan, and hover are Plotly's built-in interactions and "
        "do not change the underlying values."
    )
else:
    st.info(
        "Fewer than two principal components are available, so a "
        "PC1-versus-PC2 plot is not shown."
    )

st.header("Descriptive structural observations")
observations = build_pca_observations(result, report)
if observations:
    for observation in observations:
        st.markdown(f"- {observation}")
else:
    st.write(
        "No additional exact structural observations were generated by this "
        "descriptive summary."
    )

if current.source == "demo":
    st.warning(
        "The demo values are synthetic, its gene IDs are fictional, and its "
        "p-values are constructed. It does not define real PCA structure or "
        "support tomato biological conclusions."
    )

with st.expander("Scientific and statistical limitations"):
    st.info(
        "PCA is descriptive and unsupervised. It does not calculate a hypothesis "
        "test, a confidence interval, or a correlation p-value, and it does not "
        "perform differential-expression inference."
    )
    st.info(
        "This page does not cluster samples, assign quality scores, classify or "
        "rank samples, or recommend sample exclusion. Axes within a subspace of "
        "equal or near-equal explained variance are not uniquely identified and "
        "may differ between runs on different NumPy, BLAS, or platform "
        "combinations."
    )
    st.info(
        "No samples or genes are modified or removed. Upstream normalization, "
        "transformation, filtering, and gene selection materially affect every "
        "displayed value."
    )

st.header("Download descriptive results")
st.write(
    "Downloads are UTF-8 CSV copies of the underlying full-precision PCA result "
    "tables, not the rounded display copies. They preserve result row and column "
    "order, omit the pandas index and dtype metadata, and use empty fields for "
    "missing values."
)
st.warning(
    "CSV text is not prefixed or rewritten. Spreadsheet software may interpret "
    "formula-like leading characters in untrusted text; review such data and "
    "import it as plain text when needed."
)
_render_csv_downloads(
    (
        (
            "Download explained variance (CSV)",
            result.variance_table,
            "pca-explained-variance.csv",
        ),
        (
            "Download sample scores (CSV)",
            result.score_table,
            "pca-sample-scores.csv",
        ),
    )
)
st.caption(
    "These files contain descriptive PCA outputs only. They do not contain "
    "sample rankings, quality labels, condition-effect tests, or significance claims."
)
