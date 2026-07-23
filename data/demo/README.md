# Synthetic tomato-style demonstration data

All values in this directory are synthetic and are provided only for software
testing and demonstration. All gene identifiers are fictional. No personal,
unpublished, or real experimental data were used, and the data make no claims
about real tomato nitrate responses.

## Files and contracts

`expression_matrix.csv` has 120 rows and 7 columns:

1. `gene_id`
2. `Control_1`
3. `Control_2`
4. `Control_3`
5. `High_nitrate_1`
6. `High_nitrate_2`
7. `High_nitrate_3`

Its values are synthetic log2-normalised expression-like demonstration values,
not raw counts. The biological-replicate labels are part of the software demo;
the generated variation is not real experimental replication or statistical
evidence.

`metadata.csv` has 6 rows and 2 columns, `sample_id` and `condition`. Its rows
follow the same sample order as the expression matrix: three `Control` samples
followed by three `High_nitrate` samples.

`deg_results.csv` has 120 rows and 4 columns:

1. `gene_id`
2. `log2FoldChange`
3. `pvalue`
4. `padj`

`log2FoldChange` is the difference between the arithmetic mean of the three
quantised `High_nitrate` values and the arithmetic mean of the three quantised
`Control` values. It is a demonstration mean difference, not a fitted
coefficient or statistical estimate.

## Synthetic construction

The fixed default seed is `20260722`. Gene identifiers run from
`SYN_Solyc_0001` through `SYN_Solyc_0120`, with these internal synthetic
effect assignments:

- `SYN_Solyc_0001`–`SYN_Solyc_0020`: upregulated effects of approximately
  `+1.2` to `+2.5`;
- `SYN_Solyc_0021`–`SYN_Solyc_0040`: downregulated effects of approximately
  `-2.5` to `-1.2`;
- `SYN_Solyc_0041`–`SYN_Solyc_0120`: unchanged effects of approximately
  `-0.25` to `+0.25`.

Condition-specific, condition-centred synthetic noise makes the sample columns
non-identical while retaining a useful group signal. All numeric values are
quantised to 8 decimal places before being returned or written.

The p-values are constructed demonstration values. Stronger preset effects
generally receive smaller values, while unchanged genes generally receive
larger values. `padj` is calculated across all 120 constructed p-values using
the Benjamini–Hochberg procedure, including its reverse cumulative-minimum
step.

DESeq2 was not run. No t-test, ANOVA, or fitted RNA-seq model was run. The
p-values do not measure uncertainty and neither they nor the adjusted values
support biological significance claims.

## Regeneration

From the repository root, run:

```bash
.venv/bin/python scripts/generate_demo_data.py
```

The script regenerates only `expression_matrix.csv`, `metadata.csv`, and
`deg_results.csv` in this directory. Use `--seed` or `--output-dir` to generate
an alternative deterministic dataset elsewhere.
