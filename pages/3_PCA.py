"""Descriptive sample principal component analysis (PCA)."""

import pandas as pd
import streamlit as st

from plant_expression_explorer.dataset import get_current_dataset
from plant_expression_explorer.pca import (
    PcaComputationError,
    build_pca_observations,
    build_score_plot_data,
    compute_sample_pca,
)


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

report = current.validation_report
validation_columns = st.columns(3)
validation_columns[0].metric("Validation Errors", len(report.errors))
validation_columns[1].metric("Validation Warnings", len(report.warnings))
validation_columns[2].metric(
    "Validation Information",
    len(report.information),
)

st.header("Method")
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
    "Constant genes retained",
    result.zero_variance_gene_count,
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
    plot_data = build_score_plot_data(result)
    ratio_by_component = result.variance_table.set_index("component")[
        "explained_variance_ratio"
    ]
    pc1_ratio = float(ratio_by_component.loc[1])
    pc2_ratio = float(ratio_by_component.loc[2])
    condition_order = list(dict.fromkeys(plot_data["condition"]))
    scatter_spec = {
        "mark": {"type": "circle", "size": 120},
        "encoding": {
            "x": {
                "field": "pc1",
                "type": "quantitative",
                "title": f"PC1 ({pc1_ratio:.1%} explained variance)",
            },
            "y": {
                "field": "pc2",
                "type": "quantitative",
                "title": f"PC2 ({pc2_ratio:.1%} explained variance)",
            },
            "color": {
                "field": "condition",
                "type": "nominal",
                "sort": condition_order,
                "legend": {"title": "Condition"},
            },
            "tooltip": [
                {"field": "sample_id", "type": "nominal", "title": "Sample"},
                {"field": "condition", "type": "nominal", "title": "Condition"},
                {"field": "pc1", "type": "quantitative", "title": "PC1"},
                {"field": "pc2", "type": "quantitative", "title": "PC2"},
            ],
        },
    }
    st.vega_lite_chart(
        plot_data,
        spec=scatter_spec,
        width="stretch",
        height=420,
    )
    st.caption(
        "Colour reflects sample condition for visual grouping only; "
        "conditions were not used to fit the components. Proximity does not "
        "prove biological similarity, and separation does not prove a "
        "condition effect."
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

st.header("Scientific and statistical limitations")
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
