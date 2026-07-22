"""Upload and validate the three Phase 1 input tables."""

import streamlit as st

from plant_expression_explorer.data import (
    DataValidationError,
    align_expression_and_metadata,
    read_csv,
    validate_de_results,
)


st.title("📤 Upload Data")
st.write("Upload the three preprocessed CSV files required by the explorer.")

expression_file = st.file_uploader(
    "Normalized expression matrix",
    type="csv",
    help="Required: gene_id and one numeric column per sample.",
)
metadata_file = st.file_uploader(
    "Sample metadata",
    type="csv",
    help="Required: sample_id and condition.",
)
de_file = st.file_uploader(
    "Precomputed differential-expression results",
    type="csv",
    help="Required: gene_id, log2FoldChange, pvalue, and padj.",
)

if expression_file and metadata_file and de_file:
    try:
        expression = read_csv(expression_file)
        metadata = read_csv(metadata_file)
        de_results = read_csv(de_file)
        aligned_metadata = align_expression_and_metadata(expression, metadata)
        validate_de_results(de_results)
    except DataValidationError as exc:
        st.error(str(exc))
    else:
        st.session_state["expression"] = expression
        st.session_state["metadata"] = aligned_metadata
        st.session_state["de_results"] = de_results
        st.success("All three files are valid and their samples are aligned.")
        st.write(
            f"Loaded {len(expression):,} genes, "
            f"{len(aligned_metadata):,} samples, and "
            f"{len(de_results):,} differential-expression rows."
        )
else:
    st.info("Choose all three CSV files to run validation.")
