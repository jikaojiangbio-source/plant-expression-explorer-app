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
st.info(
    "Dataset context is descriptive text supplied with the active dataset. It is "
    "displayed verbatim but is not scientifically verified, parsed, or used to "
    "alter calculations. Missing context remains explicitly marked as not supplied."
)

st.header("What each analysis can and cannot show")
st.markdown(
    """
- **Sample Quality Control** describes the supplied matrix and sample-level value distributions. It does not assess FASTQ, alignment, library, or experimental quality and does not decide whether a sample should be removed.
- **PCA** summarizes variance in the supplied, mean-centred and unscaled matrix. Separation does not prove a condition effect, biological similarity, or statistical significance.
- **Sample Correlation** reports descriptive Pearson relationships in expression-column order. High correlation does not prove replicate validity, and unusual correlation is not an automatic exclusion rule.
- **Differential Expression** applies display thresholds only to an optional, supplied precomputed results table. It does not fit a model, determine a contrast or reference level, or establish significance, regulation, or biological importance.
- **Gene Expression** shows one exact supplied gene identifier across samples and descriptive condition groupings. It does not test condition effects or infer differential expression.
"""
)

st.header("Plant-experiment context not assessed")
st.write(
    "The application does not determine whether genotype, cultivar, tissue, "
    "developmental stage, treatment dose, timepoint, circadian timing, growth "
    "environment, block, batch, or biological replication are appropriate or "
    "confounded. Preserve these as separate metadata columns when available and "
    "interpret every view against the original experimental design."
)
st.write(
    "Gene identifiers are matched exactly. Annotation releases, gene-versus-"
    "transcript identifiers, isoform or version suffixes, and alias conventions "
    "can cause mismatches; the application does not rewrite or resolve them."
)

st.header("Outside this application's scope")
st.write(
    "Read-level quality control, alignment or pseudoalignment, quantification, "
    "filtering and normalization, power and design assessment, batch or block "
    "modelling, differential-expression fitting, multiple-testing correction, "
    "annotation mapping, enrichment, and biological validation must be performed "
    "and documented in an appropriate upstream or downstream workflow."
)
