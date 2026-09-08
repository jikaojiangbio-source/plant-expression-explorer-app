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
    "Plant Expression Explorer analyses normalized expression matrices, "
    "sample metadata, and precomputed differential-expression results."
)

st.info(
    "This application does not process FASTQ files, normalize raw counts, "
    "or independently perform differential-expression inference."
)

st.header("Planned workflow")
st.markdown(
    """
1. **Upload data** and check that the three input tables are compatible.
2. Review **sample quality-control summaries**.
3. Explore descriptive **PCA** sample scores and explained variance.
4. Inspect descriptive **sample correlation** summaries.
5. Explore supplied **differential-expression results** with descriptive thresholds.
6. Inspect planned **gene expression** views and dedicated exports.
"""
)

st.header("Phase 8")
st.write(
    "Validated synthetic-demo loading, three-file CSV upload, and session Reset "
    "are available. Descriptive Sample Quality Control summaries are also "
    "available, together with descriptive Pearson Sample Correlation summaries "
    "and descriptive PCA sample scores and explained variance, for the active "
    "dataset. Descriptive threshold exploration of supplied, precomputed "
    "differential-expression results is also available."
)
st.info(
    "Clustering, volcano plots, gene lookup, and dedicated exports remain "
    "unimplemented."
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
