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

Phase 1 provides:

- the Python project structure;
- a Streamlit home page and navigation;
- CSV upload and foundational validation;
- placeholders for the later analysis pages;
- automated tests.

PCA, sample correlation, differential-expression filtering, volcano plots,
gene lookup, demo data, and exports are not implemented yet.

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

## Input validation

Phase 1 checks required columns, unique and non-empty identifiers, numeric and
finite expression values, matching expression/metadata samples, and valid
numeric differential-expression fields. Missing differential-expression
statistics are accepted because they occur in real precomputed results, but
later analyses must report and exclude them where necessary.

Normalized expression values are not restricted to non-negative numbers,
because some valid transformed expression matrices contain negative values.

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

