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
a descriptive volcano plot on the Differential Expression page, plotting
the already-supplied, already-classified log2FoldChange/padj values with
no model fit and no calculated statistic; a "sort by group" checkbox on
the Sample Correlation heatmap that groups identical labels together
without computing any similarity or distance (not a clustering or
dendrogram-based reordering); Plotly-rendered interactive charts (zoom,
pan, box/lasso-select, hover) for the PCA sample plot, correlation
heatmap, gene-expression plot, and volcano plot; and GitHub
release-readiness foundations comprising an exact Python 3.13.9 version
pin, pinned Python dependencies, least-privilege continuous integration,
contributor guidance, and repository-readiness regression tests; an
optional, display-only "Species reference" lookup on the Gene Expression
page against a small, hand-curated, individually source-verified
per-species gene table (`data/annotations/`), matched by the same exact
`str(value)` identifier equality as everywhere else in the application; an
optional multi-gene panel on the Gene Expression page comparing several
independently looked-up exact gene IDs in one chart and wide table without
averaging or ranking them against each other; and an optional "Chart
x-axis" selector on the Gene Expression page that re-renders the per-sample
plot ordered by one numeric supplied metadata column instead of upload
order, with a controlled error (not silent coercion) when that column is
not numeric for every sample; and a Report Export page that assembles a
single PDF from already-computed descriptive result tables (dataset
summary and context, Sample Quality Control, PCA, Sample Correlation, and,
when supplied, Differential Expression), previewing every section on
screen before download and noting, rather than silently omitting, a
section whose computation is not currently possible for the active
dataset.
Differential-expression modelling and clustering remain intentionally
deferred.
