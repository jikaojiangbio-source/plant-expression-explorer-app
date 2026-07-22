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

## Phase 1 boundary

Phase 1 includes project setup, CSV upload and validation foundations, and
placeholder navigation. PCA, correlation analysis, differential-expression
exploration, volcano plots, gene lookup, and exports are intentionally deferred.

