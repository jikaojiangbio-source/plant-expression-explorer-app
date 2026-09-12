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

Phases 1–10 provide:

- the Python project structure;
- a Streamlit home page and navigation;
- a validated three-file CSV upload workflow;
- a bundled reproducible synthetic tomato-style demonstration dataset;
- a **Load synthetic demo data** route;
- one complete validated dataset bundle in the active Streamlit session;
- explicit validation feedback, table previews, source status, and Reset Data;
- descriptive, non-mutating sample quality-control summaries for the active
  dataset;
- descriptive, non-mutating Pearson sample-correlation summaries for the active
  dataset;
- descriptive, non-mutating sample PCA (mean-centred, unscaled, SVD-based
  sample scores and explained variance) for the active dataset;
- descriptive, non-mutating threshold exploration of supplied, precomputed
  differential-expression results;
- descriptive, non-mutating lookup of one exact supplied gene's expression
  values across samples with sample-condition context;
- deterministic, non-mutating UTF-8 CSV downloads of current full-precision
  descriptive result tables from each analysis page;
- placeholders for the later analysis pages;
- automated tests.

Clustering, differential-expression modelling, and volcano plots are not
implemented.

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

## Phase 6 descriptive sample correlation

The Sample Correlation page calculates Pearson correlation between every pair
of sample columns across all gene rows in the active expression matrix. Phase 6
supports Pearson only: it does not silently select another correlation method
and does not calculate correlation p-values or confidence intervals.

For `n` expression samples, the correlation matrix is `n × n`. Both axes follow
the original expression sample-column order; the matrix is never clustered or
reordered by similarity. Non-constant sample diagonal values are `1.0`.
A constant sample has exactly identical expression values across every gene,
so its Pearson correlations, including its diagonal, are undefined and remain
`NaN` in the computational result. Undefined correlations are never replaced
with zero or one.

The unique sample-pair summary excludes the diagonal and contains each unordered
pair exactly once, for `n(n-1)/2` rows. Its deterministic order follows
expression-column combinations: `(0,1), (0,2), …, (0,n-1), (1,2), …`.
Undefined pairs remain present and are marked as undefined.

The per-sample summary excludes self-correlation. It reports defined and
undefined non-self pair counts plus minimum, median, mean, and maximum
correlation calculated only from defined pairs. When no defined pair exists,
the numerical summaries remain `NaN`.

Condition-pair summaries follow condition first-appearance order while walking
the expression-matrix sample columns. Metadata row order does not control this
order. A condition with `n` samples has `n(n-1)/2` possible within-condition
pairs; two different conditions of sizes `n_a` and `n_b` have `n_a * n_b`
possible between-condition pairs. Defined and undefined pairs are counted
separately, and medians use defined correlations only. These summaries do not
perform a condition-comparison test or validate the experimental design.

A helper called directly with one non-constant sample and at least two gene
rows returns a `1 × 1` matrix containing `1.0`; a single constant sample returns
a `1 × 1` matrix containing `NaN`. In either case there are no unique non-self
pairs and the per-sample correlation statistics remain undefined. Although the
validated application input contract normally requires more than one sample,
these behaviours make the pure computation API explicit. A one-gene matrix
cannot support Pearson correlation and produces a controlled
`INSUFFICIENT_GENE_ROWS` error.

Safely coercible numeric strings are converted only in a temporary deep copy.
Missing, blank, boolean, complex, non-coercible, or infinite expression cells
produce a controlled correlation error; booleans are not treated as the numbers
zero or one. No pairwise deletion, imputation, normalisation, gene filtering,
sample filtering, or transformation is performed. Conditions are mapped by
exact sample ID even when metadata row order differs from expression column
order.

Expected input failures have a fixed priority: table types; expression row,
column, and identifier structure; metadata columns, identifiers, and conditions;
cross-table sample matching; then expression values in missing/blank,
boolean/complex, non-coercible, and infinite order. This keeps the reported
`CorrelationErrorReason` deterministic when more than one defect is present.

Heatmap data contain exactly `row_sample`, `column_sample`, `correlation`,
`defined`, and `correlation_label`. Observations are deterministic structural
facts only; they do not assign strength, quality, anomaly, rank, or exclusion
labels.

The active expression and metadata DataFrames retain their values, dtypes,
indices, and order. Result tables are newly constructed. The frozen result
dataclass prevents field rebinding but does not make nested pandas DataFrames
deeply immutable. Matrix and summary rounding, plus display of undefined values
as `N/A`, occur only on separate presentation copies.

Pearson correlations are descriptive. Correlation does not establish biological
validity and does not imply causation. Input scale, transformation, gene
filtering, and gene selection affect the result. High correlation does not
prove replicate validity, while low or negative correlation alone does not
prove that a sample is unsuitable. Unusual relationships may warrant review but
are not automatic exclusion criteria. No samples or genes are modified or
removed, and Phase 6 performs no PCA, clustering, distance analysis, batch
correction, hypothesis testing, or differential-expression inference.

No dedicated application export workflow is implemented. Streamlit components
may expose framework-provided table or chart actions.

The synthetic demo values, fictional gene IDs, and constructed p-values do not
define real correlation thresholds or support tomato biological conclusions.

## Phase 7 descriptive sample PCA

The PCA page calculates a descriptive, unsupervised principal component
analysis of the active expression matrix. Samples are PCA observations and
genes are PCA features. Output sample order follows the original expression
sample-column order; the matrix is never clustered or reordered. Sample
conditions are mapped by exact `sample_id` and are shown only as visual
colour coding; they are never used to fit the components.

Each gene is mean-centred across samples in a temporary numeric copy before
singular value decomposition (SVD); the active expression matrix is never
rewritten. Genes are not scaled to unit variance. This preserves the relative
variance structure of the supplied preprocessed matrix and avoids adding a
separate standardization step. Genes with larger variance in the supplied
matrix consequently contribute more strongly to the components; this is a
disclosed analysis policy, not a claim that such genes are more biologically
informative, and not a claim that this is a universally superior PCA method.
Upstream normalization, transformation, filtering, and gene selection
materially affect the result. No log transformation, normalisation,
filtering, trimming, aggregation, reordering, or imputation is performed. All
gene features are retained, including constant genes.

PCA uses NumPy's singular value decomposition directly (`numpy.linalg.svd`).
The component count is `min(sample_count - 1, gene_count)`. At least two
sample columns are required to estimate variance; exactly one sample column
produces a controlled `INSUFFICIENT_SAMPLE_COUNT` error rather than a
degenerate result. Exactly two samples produce exactly one component, and no
PC1-versus-PC2 plot is shown. One non-constant gene with at least two samples
is permitted and proceeds with a descriptive observation that no
dimensionality reduction occurs. Some constant genes are permitted and
retained provided total variance remains positive, with a descriptive
observation reporting their count; if every gene is constant across samples,
the temporary centred matrix has no variation and PCA produces a controlled
`ZERO_TOTAL_VARIANCE` error rather than an apparently successful result or
arbitrary components. Rank-deficient matrices with positive total variance
remain permitted. A trailing component with zero explained variance remains
in the variance table rather than being dropped, and is distinct from the
all-genes-constant case.

Finite expression values are otherwise accepted regardless of magnitude, but
if they cannot be safely represented through the float64 PCA calculation
(for example, extremely large or extremely small magnitudes that overflow or
underflow during variance, decomposition, or ratio calculations), PCA
produces a controlled `NUMERICAL_RANGE_ERROR` rather than an uncontrolled
failure, a misclassified `ZERO_TOTAL_VARIANCE` result, or non-finite scores
or ratios. Values are never automatically rescaled, clipped, or rewritten;
the source and temporary centred matrices are unchanged.

Explained variance for component `j` is `singular_value_j**2 / (sample_count
- 1)`. Explained-variance ratios use the sum of the reported components'
explained variances as the denominator, so they sum to `1.0` within a small
floating-point tolerance, not exactly. Singular values at or below a
numerical-noise tolerance (scaled to the largest singular value, matrix
dimensions, and floating-point precision) are reported as exactly `0.0` in
the derived scores and explained variance; this is treatment of derived
decomposition output only, never a change to the source or temporary centred
expression matrix.

Component ordering follows NumPy's SVD convention: singular values are
non-increasing. When singular values are distinct, this defines component
order unambiguously. When singular values are equal or numerically
indistinguishable, the individual axes within that tied subspace are not
uniquely identified: different NumPy, BLAS, or platform combinations may
return a different, equally valid, orthonormal basis for the same tied
subspace, even though the components' explained variance and the
sample-to-sample distances they represent remain the same. A deterministic
sign convention (the sample with the largest-magnitude score on a component
is positive, ties broken by expression sample-column order) resolves the
positive/negative ambiguity for a single non-degenerate component only; it
does not identify axes inside a tied subspace, and pinning the NumPy version
does not guarantee identical floating-point coordinates across platforms or
BLAS implementations.

The active expression and metadata DataFrames retain their values, dtypes,
indices, and order. Result tables are newly constructed. The frozen result
dataclass prevents field rebinding but does not make nested pandas DataFrames
deeply immutable. Score and variance-table rounding for display occurs only
on separate presentation copies.

PCA is descriptive and unsupervised. Proximity in a PCA plot does not prove
biological similarity or replicate validity, and separation between groups
does not prove a condition effect or establish statistical significance.
PCA does not identify failed or low-quality samples, does not calculate a
hypothesis test, a confidence interval, or a correlation p-value, and does
not perform clustering or differential-expression inference. No samples or
genes are modified or removed. The synthetic demo values and fictional gene
IDs do not define real PCA structure or support tomato biological
conclusions.

## Phase 8 descriptive differential-expression exploration

The Differential Expression page displays the complete supplied, precomputed
results table and assigns one exploratory status to every accepted row. Users
choose an adjusted-p-value threshold in the inclusive range `[0, 1]` and an
absolute `log2FoldChange` threshold strictly greater than zero. Comparisons are
inclusive: positive matches require `padj` at or below the selected threshold
and `log2FoldChange` at or above the positive threshold; negative matches use
the corresponding negative boundary. Every other row with usable values is an
other-evaluable row.

A row is not evaluable when its supplied `padj` or `log2FoldChange` is missing
or blank. Missing `pvalue` alone does not affect classification when the two
classification fields are usable. Missing values remain missing, including
`padj`; zero adjusted p-values remain valid zeros. Safely coercible numeric
strings are converted only in temporary calculation Series. The original DEG
table, its index, row and column order, values, dtypes, and additional columns
are not rewritten. No identifier is trimmed, normalised, deduplicated, or
silently intersected with expression data.

These categories mean only that rows meet the current user-selected exploratory
thresholds. Positive and negative labels refer solely to the sign of the
supplied fold-change value because the application does not know the contrast
direction or reference level. Threshold matches do not establish statistical
significance or biological importance. Phase 8 does not calculate or modify
p-values, run DESeq2 or another differential-expression model, process FASTQ
files, create a volcano plot, perform gene lookup, or provide a dedicated
export workflow.

## Phase 9 descriptive gene-expression lookup

The Gene Expression page lets users select one exact `gene_id` supplied in the
active expression matrix and inspect its preprocessed expression values across
all samples. Gene options and result rows preserve expression-matrix order. No
gene is selected by default, ranked, recommended, or looked up through aliases.

Identifier matching follows the application's established contract: each
validated source identifier is represented by exact `str(value)` for comparison
and display. These representations are not trimmed, case-folded, normalized, or
otherwise rewritten. The source DataFrames and their original scalar values,
dtypes, indices, row order, and column order remain unchanged. The per-sample
result table copies the original selected-row expression scalars; safely
coercible numeric strings are converted only in independent numeric working
data used for charting and descriptive summaries.

Each result contains one row per expression sample column, in expression-column
order. Conditions are mapped by exact sample ID even when metadata row order
differs. A sample-set mismatch is a controlled error in both directions; Phase
9 does not take an intersection, fill missing metadata, or reorder either source
table.

The per-sample chart contains unconnected points in explicit expression sample
order. Condition colour provides context only. The condition-grouped table
reports sample membership, count, minimum, median, arithmetic mean, maximum,
and sample standard deviation (`ddof=1`) on the supplied scale. Standard
deviation is undefined for a one-sample condition and is retained as `NaN` in
the computational result while a separate display copy shows `N/A`. Derived
statistics use scaled calculations to avoid overflowing intermediate sums. All
returned statistics must be finite except that defined one-sample `NaN`; if a
derived statistic cannot be represented safely as a finite float64 value, the
lookup returns the controlled `NUMERICAL_RANGE_ERROR` and no partial summary.
These are disclosed descriptive aggregations of the displayed values, not
effect estimates, tests, confidence intervals, or evidence of biological
replication.

The pure lookup API defines a one-sample result, although the application's
validated upload contract continues to require at least two expression sample
columns. With one defensively supplied sample, the value and one-row condition
summary remain available but an across-sample plot is omitted. One-condition
data remain descriptive and do not create a between-condition comparison.

Phase 9 performs no FASTQ processing, normalization, transformation, filtering,
batch correction, imputation, statistical testing, differential-expression
inference, contrast or reference-level interpretation, biological-importance
classification, volcano plotting, or dedicated export workflow.

## Phase 10 descriptive result exports

Each implemented analysis page provides dedicated downloads for its current
successfully calculated descriptive result tables:

- Sample Quality Control exports condition membership and per-sample statistics.
- PCA exports full-precision explained variance and sample scores.
- Sample Correlation exports the full matrix in ordered long form, plus
  per-sample, unique-pair, and condition-pair summaries.
- Differential Expression exports category counts and the complete supplied
  results with the current exploratory status column. Both files include the
  exact adjusted-p-value and absolute log2-fold-change thresholds applied to
  every exported classification.
- Gene Expression exports the selected gene's per-sample values and
  condition-grouped descriptive summary; the exact selected `gene_id` is an
  explicit column in both files.

Downloads are built in memory and do not write files on the application server.
They are generated from the underlying result DataFrames, not from rounded or
`N/A`-formatted display copies. Result row and column order is retained and the
active dataset and result objects are not mutated. Repeated export of the same
result produces the same UTF-8 bytes with `\n` line endings.

CSV is a serialized interchange format, not a lossless pandas archive. Exports
do not include the pandas index or dtype metadata. Missing scalar values are
written as empty fields, so CSV alone cannot distinguish a missing value from a
source empty string. Tuple/list sample memberships are exported as compact JSON
arrays in the `sample_ids` field, making exact membership boundaries explicit.
No identifiers are placed unvalidated into filenames. Every download for a page
is built before any download button is shown; unsupported objects, infinite
values, invalid or conflicting columns, unsafe filenames, or serialization
failures therefore produce one controlled error and no partial set of download
buttons.

CSV text cells are preserved without spreadsheet-specific prefixing. Some
spreadsheet programs may interpret formula-like leading characters in untrusted
text as formulas, so each download section tells users to review such content
and import it as plain text when needed. The application does not silently alter
identifiers or annotations to mitigate consumer-specific behaviour.

These downloads do not recreate the original uploaded CSV bytes and are not
described as backups of the source files. They do not add normalization,
transformation, filtering, imputation, ranking, hypothesis testing,
differential-expression inference, reference-level or contrast interpretation,
biological-importance claims, chart-image export, Excel/ZIP packaging, or a
combined report.

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
