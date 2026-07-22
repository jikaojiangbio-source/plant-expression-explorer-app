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
3. Explore **PCA** and **sample correlation**.
4. Filter supplied **differential-expression results**.
5. Inspect **gene expression** and export selected results.
"""
)

st.header("Phase 1")
st.write(
    "CSV upload and validation are available now. The analysis pages are "
    "visible as placeholders so the intended MVP structure is clear."
)
st.page_link("pages/1_Upload_Data.py", label="Start with Upload Data", icon="📤")

