# Plant Expression Explorer contributor guide

## Project scope

This is a beginner-friendly Streamlit application for exploring preprocessed
plant transcriptomics data. It does not process FASTQ files, normalize raw
counts, correct batch effects, or perform differential-expression inference.

## Development guidelines

- Keep Streamlit presentation code in `app.py` and `pages/`.
- Keep reusable data and analysis logic in `plant_expression_explorer/` as
  small, testable functions.
- Prefer clear pandas operations over clever abstractions.
- Validate user-supplied data before analysis and return actionable errors.
- Do not silently impute, transform, aggregate, or discard biological data.
- Label descriptive results clearly and document scientific assumptions.
- Add or update pytest tests whenever behavior changes.
- Run `python -m pytest` before handing changes back for review.

## Phase 1 boundary (historical)

Phase 1 includes project setup, CSV upload and validation foundations, and
placeholder navigation. PCA, correlation analysis, differential-expression
exploration, volcano plots, gene lookup, and exports are intentionally deferred.

## Current implementation boundary

As of Phase 12, the application implements: project setup; a validated CSV
upload and data-loading workflow requiring an expression matrix and sample
metadata, with precomputed differential-expression results now optional; a
bundled reproducible synthetic demo dataset; descriptive, non-mutating sample
quality-control summaries; descriptive, non-mutating Pearson sample-correlation
summaries; and descriptive, non-mutating sample PCA (mean-centred, unscaled,
SVD-based sample scores and explained variance); with a shared, session-wide
"group by" selector, on PCA, Sample Quality Control, Sample Correlation, and
Gene Expression, that relabels or recomputes only already-disclosed
descriptive statistics for any additional supplied metadata column, never a
PCA component or a Pearson correlation; and descriptive, non-mutating
threshold exploration of supplied, precomputed differential-expression
results when supplied; and descriptive, non-mutating lookup of one exact
supplied gene's expression values across samples with sample-condition
context, with a substring search box once a dataset supplies more than 200
genes; and deterministic, non-mutating UTF-8 CSV downloads of the current
full-precision descriptive result tables on each analysis page; and optional,
verbatim, non-scientifically-verified dataset context (title, organism,
expression scale, upstream normalization, reference annotation, feature
level, DE contrast description, notes) shown on every analysis page; and
example CSV templates, an at-a-glance dataset overview card, and clearer
error messages for common real-world upload mistakes (delimiter mismatch,
blank spreadsheet-export columns, ad hoc missing-value placeholders); and
GitHub release-readiness foundations comprising an exact Python 3.13.9
version pin, pinned Python dependencies, least-privilege continuous
integration, contributor guidance, and repository-readiness regression
tests. Differential-expression modelling and volcano plots remain
intentionally deferred.
