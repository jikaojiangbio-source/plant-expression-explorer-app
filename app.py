"""Home page for Plant Expression Explorer."""

import streamlit as st


st.set_page_config(
    page_title="Plant Expression Explorer",
    page_icon="🌱",
    layout="wide",
)

st.title("🌱 Plant Expression Explorer")
st.subheader("Explore preprocessed plant transcriptomics data")

st.write(
    "Plant Expression Explorer analyses preprocessed expression matrices and "
    "sample metadata, with optional precomputed differential-expression results."
)

st.info(
    "This application does not process FASTQ files, normalize raw counts, "
    "or independently perform differential-expression inference."
)

st.header("Planned workflow")
st.markdown(
    """
1. **Upload data** and check that the required expression and metadata tables are compatible.
2. Review **sample quality-control summaries**.
3. Explore descriptive **PCA** sample scores and explained variance.
4. Inspect descriptive **sample correlation** summaries.
5. Explore supplied **differential-expression results** with descriptive thresholds.
6. Look up one exact supplied gene in descriptive **gene expression** views.
7. Download the current descriptive result tables as dedicated UTF-8 CSV files.
"""
)

st.header("Phase 12")
st.write(
    "Validated synthetic-demo loading, expression-plus-metadata CSV upload, and "
    "session Reset "
    "are available. Descriptive Sample Quality Control summaries are also "
    "available, together with descriptive Pearson Sample Correlation summaries "
    "and descriptive PCA sample scores and explained variance, for the active "
    "dataset. Descriptive threshold exploration of supplied, precomputed "
    "differential-expression results is also available. Exact, descriptive "
    "single-gene expression lookup is available for the active dataset. "
    "Dedicated in-memory CSV downloads of the current unrounded descriptive "
    "result tables are available on each analysis page. Precomputed "
    "differential-expression results are now optional, and descriptive dataset "
    "context can record the supplied scale, organism, reference annotation, "
    "feature level, and analysis provenance without changing calculations."
)
st.info(
    "Clustering and differential-expression modelling remain unimplemented. "
    "A descriptive volcano plot of already-supplied, already-classified "
    "differential-expression results is available on the Differential "
    "Expression page; it fits no model and calculates no statistic."
)
st.page_link("pages/1_Upload_Data.py", label="Start with Upload Data", icon="📤")
st.page_link(
    "pages/2_Sample_Quality_Control.py",
    label="Review Sample Quality Control",
    icon="🧪",
)
st.page_link(
    "pages/3_PCA.py",
    label="Review PCA",
    icon="📊",
)
st.page_link(
    "pages/4_Sample_Correlation.py",
    label="Review Sample Correlation",
    icon="🔥",
)
st.page_link(
    "pages/5_Differential_Expression.py",
    label="Explore Differential Expression",
    icon="🧬",
)
st.page_link(
    "pages/6_Gene_Expression.py",
    label="Explore Gene Expression",
    icon="🌿",
)
