# Plant Expression Explorer

Plant Expression Explorer is a beginner-friendly Streamlit application for
exploring **preprocessed plant transcriptomics data**.

The application accepts three CSV inputs:

1. A normalized expression matrix with `gene_id` and one numeric column per
   sample.
2. Sample metadata with `sample_id` and `condition`.
3. Precomputed differential-expression results with `gene_id`,
   `log2FoldChange`, `pvalue`, and `padj`.

## Current status

Phases 1–5 provide:

- the Python project structure;
- a Streamlit home page and navigation;
- a validated three-file CSV upload workflow;
- a bundled reproducible synthetic tomato-style demonstration dataset;
- a **Load synthetic demo data** route;
- one complete validated dataset bundle in the active Streamlit session;
- explicit validation feedback, table previews, source status, and Reset Data;
- descriptive, non-mutating sample quality-control summaries for the active
  dataset;
- placeholders for the later analysis pages;
- automated tests.

PCA, sample correlation, clustering, differential-expression filtering,
volcano plots, gene lookup, and exports are not implemented yet.

## Setup

Python 3.11 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the application:

```bash
streamlit run app.py
```

Run the tests:

```bash
python -m pytest
```

## Bundled synthetic demonstration data

Phase 3 provides three reproducible CSV files in
[`data/demo/`](data/demo/README.md). They contain fictional gene identifiers
and synthetic log2-normalised expression-like values for three `Control` and
three `High_nitrate` samples.

Regenerate the files with the fixed default seed:

```bash
.venv/bin/python scripts/generate_demo_data.py
```

Select **Load synthetic demo data** on the Upload Data page to read the three
committed files through the same CSV reader and aggregate validation workflow
used for uploaded data. Page loading never regenerates or rewrites these files.

All values are synthetic and are for software testing and demonstration only.
The gene identifiers are fictional, the expression values are not raw counts,
and the p-values are directly constructed demonstration values. No personal or
unpublished experimental data were used, no real tomato nitrate-response
claims are made, and DESeq2, a t-test, or another fitted RNA-seq model was not
run.

## Phase 4 data-loading workflow

The Upload Data page provides two routes:

1. Load the three bundled synthetic demo CSV files.
2. Supply an expression matrix, sample metadata, and precomputed
   differential-expression results.

Uploaded files are not fully validated until all three required sources are
present. Both routes use the existing CSV reader and the same Phase 2 aggregate
validation API. CSV and header failures are presented as structured Errors
rather than raw exception tracebacks.

After successful validation, the application stores one complete current
dataset bundle containing the three DataFrames, their source and source label,
and the complete validation report. A frozen dataclass prevents bundle fields
from being rebound, but it does not make nested pandas DataFrames intrinsically
immutable. Application helpers and rendering therefore do not modify these
DataFrames in place.

Validation outcomes control activation as follows:

- Error — the candidate does not become the current dataset.
- Warning — shown to the user but does not block activation.
- Information — shown to the user and does not block activation.

If a valid dataset is already active and a replacement candidate fails, the
previous successfully validated bundle remains active. The failed candidate is
reported separately so its filenames or source are never presented as the
active dataset.

**Reset Data** removes the current bundle, failed-candidate feedback, uploader
widget values, and legacy application data keys from the active session. It
does not delete user files or bundled demo files from disk. Uploaded tables are
retained in the active Streamlit session for application use; the application
does not intentionally write uploaded tables to project files.

Phase 4 performs data loading and validation only. Phase 5 consumes its single
active validated dataset bundle without creating another session-state copy.

## Phase 5 descriptive sample quality control

The Sample Quality Control page describes the active expression matrix and its
sample metadata. It reports gene, sample, expression-cell, condition, missing,
non-finite, zero, and negative-value counts; whether metadata and expression
sample orders correspond; constant sample columns; and zero-variance genes.

For every sample, in expression-column order, it reports the mapped condition,
gene count, missing/non-finite/zero/negative counts, minimum, first quartile,
median, arithmetic mean, third quartile, maximum, standard deviation, and
interquartile range. Quartiles use pandas linear interpolation. Standard
deviation uses pandas sample standard deviation with `ddof=1`, and IQR is
`Q3 - Q1`. For a one-gene matrix, the `ddof=1` standard deviation is undefined:
the computational result retains `NaN`, while the page displays `N/A` from a
separate display copy.

A constant sample is a sample column whose expression values are exactly
identical across all genes. A zero-variance gene is a gene row whose expression
values are exactly identical across all samples. Both use exact equality. They
are reported but never removed.

Safely coercible numeric strings are converted only in a temporary numeric
copy used for calculation. The active expression and metadata DataFrames,
including their values, dtypes, indices, row order, and column order, are not
modified. Missing, non-coercible, and infinite expression cells produce a
controlled QC error; cells, genes, and samples are never silently skipped,
imputed, trimmed, transformed, or discarded.

The result is held in a frozen dataclass. Freezing prevents result fields from
being rebound, but it does not make nested pandas DataFrames intrinsically
immutable.

These summaries are descriptive and do not establish biological validity.
Mean, standard deviation, median, IQR, zero values, and negative values have
scale-dependent interpretations, and uploaded matrices may already have been
filtered upstream. The page does not assign quality scores, classify unusual
samples, recommend exclusion, or automatically filter data. It performs no
normalisation, PCA, correlation, clustering, batch correction, hypothesis
testing, or differential-expression inference.

## Input validation

Phase 2 validates three preprocessed CSV tables before making a validated
dataset bundle available to downstream pages.

### Expression matrix

- Required column: `gene_id`.
- At least two sample columns are required in addition to `gene_id`.
- Gene identifiers must be non-empty and unique.
- Expression values must be numeric or safely coercible to numeric values,
  finite, and non-missing.
- Negative expression values are permitted because transformed or centred
  matrices may legitimately contain them. They produce a Warning, remain
  unchanged, and are not treated as an Error.

### Sample metadata

- Required columns: `sample_id` and `condition`.
- Sample identifiers must be non-empty and unique.
- Condition values must be non-empty.
- Additional columns are preserved but are not scientifically interpreted by
  Phase 2 validation.

### Precomputed differential-expression results

- Required columns: `gene_id`, `log2FoldChange`, `pvalue`, and `padj`.
- Gene identifiers must be non-empty and unique.
- Non-missing statistic values must be numeric or safely coercible to numeric
  values and finite.
- Non-missing `pvalue` and `padj` values must be within the inclusive range
  `[0, 1]`.
- Partially missing statistics are preserved and reported as a Warning.
- Missing `padj` values are never replaced with zero.
- If every value in a required statistic column is missing, validation reports
  `NO_USABLE_VALUES` as an Error. Non-numeric or infinite values are blocked by
  their corresponding Errors.
- Phase 2 does not classify genes as significant or non-significant.

### Validation outcomes

- Error — the validated dataset bundle is not made available to downstream
  pages.
- Warning — non-blocking, but requires user attention.
- Information — non-blocking observation.

Duplicate expression-matrix `gene_id`, metadata `sample_id`, and
differential-expression `gene_id` values are Errors. Duplicate rows are
reported but are not silently removed or aggregated.

Safely coercible numeric strings are used only to determine validation
outcomes. Validation does not silently rewrite the original values or change
their stored data types.

### Cross-file consistency

Partial expression/metadata sample mismatches are reported in both directions;
no sample overlap is reported as one blocking Error. Partial expression/DEG
gene-coverage mismatches are reported in both directions as Warnings; no gene
overlap is reported as one blocking Error.

Cross-file checks compare identifier sets without changing the source tables.
Validation does not trim identifiers, retain only intersections, remove rows,
impute values, normalize values, or reorder either table.

Phase 2 validation checks table structure, identifier consistency, and basic
numeric constraints. Passing validation does not prove that the experimental
design, statistical analysis, or biological data quality is appropriate.

## Scientific limitations

Plant Expression Explorer is a visualization and exploration tool. It does not:

- process FASTQ reads or assess sequencing/read quality;
- normalize raw counts;
- detect or correct batch effects;
- fit statistical models or calculate differential expression;
- establish causal or biological significance.

Users remain responsible for the upstream normalization, experimental design,
statistical model, contrasts, covariates, and multiple-testing correction used
to create the supplied files.
