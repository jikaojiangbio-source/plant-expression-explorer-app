# User guide

Detailed reference for each page and input contract, moved verbatim from the
earlier README. Start with the [README](../README.md) for an overview.

## Bundled synthetic demonstration data

Phase 3 provides three reproducible CSV files in
[`data/demo/`](../data/demo/README.md). They contain fictional gene identifiers
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

## Volcano plot

The Differential Expression page includes a descriptive volcano plot: the
supplied `log2FoldChange` against `-log10(padj)` for evaluable rows, coloured
by the same exploratory threshold status already shown in the category table,
with dashed guides at the exact applied thresholds. It introduces no new
statistic; it plots values already computed by threshold classification. A
supplied `padj` of exactly 0 has no finite `-log10` value, so such rows are
excluded from the plot only (never from any table or download) and counted in
an explicit disclosure message.

## MA plot

The Differential Expression page also includes a descriptive MA plot: each
evaluable gene's `mean_expression` (the arithmetic mean of that gene's own
supplied per-sample expression values, missing values excluded and never
imputed) against the supplied `log2FoldChange`, coloured by the same
exploratory threshold status. `mean_expression` is computed from the active
expression matrix already used elsewhere in this application; it is not a
library-size-normalized `baseMean` (DESeq2) or `AveExpr` (limma) statistic
from any specific external tool, and no statistic is calculated, adjusted,
or inferred here. An evaluable row whose gene_id has no exact, unambiguous
match in the expression matrix, or whose every supplied expression value
for that gene is missing, is excluded from the plot only (never from any
table or download) and counted in an explicit disclosure message.

## Visual design system

`plant_expression_explorer/theme.py` defines a shared, presentation-only
stylesheet (`inject_global_styles()`) called near the top of every page. It
restyles existing Streamlit elements only — alerts as quieter cards with a
coloured left accent instead of a solid pastel fill, `st.metric` as a small
card with a brand-coloured top border, buttons and `st.page_link` rows as
rounded cards with a hover lift, a sidebar-link hover state, and the Inter
typeface — and hides the Streamlit Community Cloud "Deploy" affordance, since
this is a finished, purpose-built app rather than a work-in-progress
template. No rule changes any element's text, order, or presence, and no
rule affects computation.

The home page adds a two-column hero (title/scope/primary call to action
beside a hand-authored, clearly-labelled illustrative SVG preview of the
app's chart types — abstract shapes, not real data), a dark rounded "How it
works" panel presenting the same seven workflow steps as numbered cards, and
a row of capability chips. These are additive presentational elements built
from the same underlying text already required elsewhere; no scientific
wording changed.

## Species reference (optional, display-only)

The Gene Expression page offers an optional "Species reference" selector.
When a supported species is chosen and the currently selected gene ID
exactly matches an entry in `data/annotations/<species>.csv`, its published
gene symbol, a short description, and its data source are shown. Matching is
exact `str(value)` equality only, identical to every other identifier match
in this application: no case-folding, trimming, or alias resolution, and no
computation ever uses this lookup.

This is a small, hand-curated list of well-known reference/marker genes per
species (currently Arabidopsis thaliana, rice, maize, and soybean; not yet
tomato — see `data/annotations/README.md` for why), not a genome annotation.
Most real gene IDs will not have an entry; that is expected, not an error.
Every bundled `gene_id`/`symbol` pair was individually verified against the
Ensembl Plants REST API in the session that added it; rows whose description
was hand-written from established literature (because the API returned
none) are labelled as such in the `source` column rather than attributed to
Ensembl.

An unmatched gene ID is often not because the gene is unlisted, but because
several of these species have two or more non-interchangeable identifier
systems for the same genome (for example rice's RAP-DB vs MSU/TIGR locus
IDs, maize's Ensembl Plants Zm00001eb-style IDs vs older GRMZM/Zm00001d
IDs from earlier assembly versions, soybean's underscore- vs dot-notation
Wm82.a2 IDs, or an Arabidopsis TAIR locus ID with a transcript/splice-
variant version suffix). When the supplied gene ID's shape matches one of
these well-documented alternate systems, an additional caption names which
system it looks like and which system the bundled list actually uses. This
is a descriptive note only: `plant_expression_explorer/annotations.py`
never rewrites, strips, or looks up the identifier under any other form,
consistent with the exact-match-only contract above.

## Custom annotation upload (optional, display-only)

The bundled species reference above covers only a few dozen marker genes
across four species, so the Gene Expression page also offers an
independent "Upload gene annotation CSV" file uploader. A user-supplied
CSV with `gene_id`, `symbol`, and `description` columns (an optional
`source` column is shown verbatim when present and non-blank; otherwise a
row's source is reported as an unverified user upload) is indexed by exact
`gene_id` and matched against the currently selected gene, using the same
exact `str(value)` equality as every other identifier lookup in this
application. Unlike the bundled per-species lists, an uploaded annotation
file's accuracy is never independently verified by this application: the
uploader is responsible for its contents, exactly as for the expression,
metadata, and differential-expression tables. Like the bundled species
reference, this is never a genome annotation and is never used in any
calculation, and a malformed file (missing a required column, a blank
`gene_id`/`symbol`, or a duplicate `gene_id`) produces a controlled error
naming the problem rather than a partial or guessed table.

## Multi-gene panel (optional)

The Gene Expression page offers an optional "Compare with additional gene
IDs" multiselect, independent of the primary "Exact gene ID" selector used
throughout the rest of the page. Choosing one or more additional exact gene
IDs adds a combined panel above the single-gene detail view: a Plotly line
chart (one coloured series per gene, samples in expression-column order) and
a wide-format table (`sample_id`, `condition`, one column per selected gene)
with its own CSV download. Each gene's values come from an independent call
to the same per-sample lookup used by the single-gene view; genes are never
averaged, combined into a score, or ranked against one another, and adding
this panel does not change the primary gene's detail view below it.

## Time-series chart x-axis (optional)

When the active dataset's sample metadata has at least one column beyond
`sample_id` and `condition`, the Gene Expression page offers a "Chart
x-axis" selector. The default, "Sample (upload order)", is the existing
categorical per-sample plot. Choosing "Numeric time/order: `<column>`"
instead re-renders the per-sample plot as a Plotly line chart with that
column's exact supplied value, converted to numeric, on the x-axis; points
are connected only to make the sample sequence easier to trace, not as a
fitted trend or interpolation, and replicates sharing one time value are
plotted individually rather than averaged. If the chosen column is not
numeric for every sample, the page shows a controlled error naming every
offending sample instead of skipping, coercing, or guessing values.

## Report export (optional, PDF)

The Report Export page assembles a single downloadable PDF from
already-computed descriptive result tables: a dataset summary, dataset
context, the Sample Quality Control condition summary, the PCA explained-
variance table, the Sample Correlation condition-pair summary, and (when
supplied) the Differential Expression category summary. It performs no new
calculation; every table is produced by the same functions used on their
respective pages, using the same shared "group by" column and, for
Differential Expression, the same exploratory thresholds set on that page
(or the same defaults it uses when unvisited). The page also renders every
section on screen before offering the download, so the PDF never contains
content the user has not already seen. A section whose underlying
computation is not currently possible for the active dataset (for example,
no differential-expression results supplied) is included as an explicit
note rather than silently dropped.

The report does not include the Gene Expression page's per-gene lookup or
its multi-gene/time-series charts, since those depend on a page-local gene
and chart choice rather than dataset-level state; it also contains no chart
images. `plant_expression_explorer/report.py` builds the PDF with
`reportlab`; a table cell, title, or note containing a character its
built-in font cannot render (outside Windows-1252) produces a controlled
error naming the offending section instead of silently dropping or
mangling that character.

The Dataset summary table also includes a SHA-256 checksum for each
uploaded/demo input file, computed from the exact source CSV bytes before
any parsing (`plant_expression_explorer/dataset.py`'s `DatasetChecksums`),
so a later reader can confirm a specific report was generated from a
specific, byte-identical input file. Two uploads of the same file always
produce the same digest, and this is purely a disclosure: it is never
compared, validated, or used to gate anything. A dataset built without
supplying checksums (as some older or directly-constructed bundles are)
shows an explicit "not available" note instead of a fabricated value.

## Input validation

Phase 2 validates three preprocessed CSV tables before making a validated
dataset bundle available to downstream pages.

### Expression matrix

- Required column: `gene_id`.
- At least two sample columns are required in addition to `gene_id`.
- Gene identifiers must be non-empty and unique.
- Expression values must be numeric or safely coercible to numeric values,
  and finite.
- A missing expression value is permitted (a Warning, not an Error): it is
  retained as missing and never imputed. Each descriptive calculation
  documents its own handling — Sample Quality Control's per-sample
  statistics and PCA's and Sample Correlation's constant-sample check
  exclude missing values gene-wise; PCA excludes any gene with a missing
  value from that calculation entirely (disclosed exact count); Sample
  Correlation computes each sample pair from only the gene rows where both
  samples have a value ("pairwise complete"), reporting a pair as undefined
  when fewer than 2 such rows exist; Gene Expression excludes a missing
  sample from that gene's condition summary statistics. A sample column
  with no non-missing value at all is still a blocking Error, since there is
  nothing to summarize.
- Negative expression values are permitted because transformed or centred
  matrices may legitimately contain them. They produce a Warning, remain
  unchanged, and are not treated as an Error.

### Sample metadata

- Required columns: `sample_id` and `condition`.
- Sample identifiers must be non-empty and unique.
- Condition values must be non-empty.
- Additional columns are preserved but are not scientifically interpreted by
  Phase 2 validation. Phase 12 lets the PCA page use one such column, chosen
  by the user, for display-only sample-plot colouring (see "Phase 12 real-world
  usability").

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
