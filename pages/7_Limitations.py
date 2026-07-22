"""Scientific scope and limitations."""

import streamlit as st


st.title("⚠️ Scientific Limitations")
st.warning(
    "Plant Expression Explorer is an exploratory visualization tool for "
    "preprocessed data. It does not process FASTQ files, normalize raw counts, "
    "correct batch effects, or calculate differential expression."
)
st.write(
    "Interpretation depends on the upstream experimental design, normalization, "
    "statistical model, contrasts, covariates, and multiple-testing correction."
)

