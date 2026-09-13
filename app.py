"""Home page for Plant Expression Explorer."""

import streamlit as st

from plant_expression_explorer.theme import (
    inject_global_styles,
    render_chip_row,
    render_hero_illustration,
)

st.set_page_config(
    page_title="Plant Expression Explorer",
    page_icon="🌱",
    layout="wide",
)
inject_global_styles()

hero_column, illustration_column = st.columns([1.15, 1], gap="large")

with hero_column:
    st.html('<span class="pee-eyebrow">Descriptive · Non-inferential</span>')
    st.title("🌱 Plant Expression Explorer")
    st.subheader("Explore preprocessed plant transcriptomics data")

    st.write(
        "Plant Expression Explorer analyses preprocessed expression matrices and "
        "sample metadata, with optional precomputed differential-expression results."
    )

    st.page_link(
        "pages/1_Upload_Data.py", label="Start with Upload Data", icon="📤"
    )

with illustration_column:
    render_hero_illustration()
    st.caption(
        "A stylised preview of this app's descriptive chart types; not real data."
    )

st.info(
    "This application does not process FASTQ files, normalize raw counts, "
    "or independently perform differential-expression inference."
)

st.markdown(
    """
<div class="pee-dark-panel">
<h2>How it works</h2>
<ol>
<li><strong>Upload data</strong> and check that the required expression and metadata tables are compatible.</li>
<li>Review <strong>sample quality-control summaries</strong>.</li>
<li>Explore descriptive <strong>PCA</strong> sample scores and explained variance.</li>
<li>Inspect descriptive <strong>sample correlation</strong> summaries.</li>
<li>Explore supplied <strong>differential-expression results</strong> with descriptive thresholds.</li>
<li>Look up one exact supplied gene in descriptive <strong>gene expression</strong> views.</li>
<li>Download the current descriptive result tables as dedicated UTF-8 CSV files.</li>
</ol>
</div>
""",
    unsafe_allow_html=True,
)

st.header("Phase 12")
st.write(
    "Validated synthetic-demo loading, expression-plus-metadata CSV upload, and "
    "session Reset are available, together with descriptive Sample Quality "
    "Control, Pearson Sample Correlation, and PCA summaries for the active "
    "dataset. Precomputed differential-expression results are now optional, and "
    "descriptive dataset context can record the supplied scale, organism, "
    "reference annotation, feature level, and analysis provenance without "
    "changing calculations."
)
render_chip_row(
    [
        "Sample QC",
        "PCA",
        "Sample correlation",
        "DE threshold exploration",
        "Volcano plot",
        "Gene lookup",
        "Dataset provenance",
        "CSV exports",
    ]
)
st.info(
    "Clustering and differential-expression modelling remain unimplemented. "
    "A descriptive volcano plot of already-supplied, already-classified "
    "differential-expression results is available on the Differential "
    "Expression page; it fits no model and calculates no statistic."
)

st.header("Browse pages")
link_columns = st.columns(3)
link_columns[0].page_link(
    "pages/2_Sample_Quality_Control.py",
    label="Review Sample Quality Control",
    icon="🧪",
)
link_columns[1].page_link(
    "pages/3_PCA.py",
    label="Review PCA",
    icon="📊",
)
link_columns[2].page_link(
    "pages/4_Sample_Correlation.py",
    label="Review Sample Correlation",
    icon="🔥",
)
link_columns2 = st.columns(3)
link_columns2[0].page_link(
    "pages/5_Differential_Expression.py",
    label="Explore Differential Expression",
    icon="🧬",
)
link_columns2[1].page_link(
    "pages/6_Gene_Expression.py",
    label="Explore Gene Expression",
    icon="🌿",
)
