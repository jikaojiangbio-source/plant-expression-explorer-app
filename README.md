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

Phases 1–3 provide:

- the Python project structure;
- a Streamlit home page and navigation;
- CSV upload and foundational validation;
- a bundled reproducible synthetic tomato-style demonstration dataset;
- placeholders for the later analysis pages;
- automated tests.

PCA, sample correlation, differential-expression filtering, volcano plots,
gene lookup, a Load Demo button, and exports are not implemented yet.

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

All values are synthetic and are for software testing and demonstration only.
The expression values are not raw counts, and the p-values are directly
constructed demonstration values. No personal or unpublished experimental data
were used, no real tomato nitrate-response claims are made, and DESeq2, a
t-test, or another fitted RNA-seq model was not run. A Load Demo button is not
implemented yet; the files can be uploaded manually through the existing input
workflow.

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
