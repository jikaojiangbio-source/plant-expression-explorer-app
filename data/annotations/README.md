# Small, hand-curated species gene-annotation reference

These files back the optional "Species reference" lookup on the Gene
Expression page. They are **not** a genome annotation. Each file lists a
small, hand-picked set of well-known reference/marker genes for one species
(a few to a dozen genes), so that a user who selects a matching species and
looks up one of these specific genes sees its published gene symbol and a
short functional description alongside the exact supplied gene ID. Any gene
ID not in these files simply shows no reference entry; this is expected, not
an error, since real datasets contain thousands of genes and these files
cover only a handful of especially well-known ones per species.

## Files and contract

One CSV per species, named by its Ensembl-style species key:

- `arabidopsis_thaliana.csv` (12 genes)
- `oryza_sativa.csv` (5 genes)
- `zea_mays.csv` (7 genes)
- `glycine_max.csv` (4 genes)

Each file has exactly four columns: `gene_id`, `symbol`, `description`,
`source`. `gene_id` uses each species' own genome-consortium locus identifier
convention (for example `AT1G65480` for Arabidopsis, `Os06g0157700` for
rice), matching what a real preprocessed expression matrix from that species
would typically use. Lookup in the application matches `gene_id` by exact
`str(value)` equality only, consistent with every other identifier match in
this application: no case-folding, trimming, or alias resolution.

## Provenance

Every `gene_id` and `symbol` pair was individually verified in this
session (2026-09-13) against the Ensembl Plants REST API
(`https://rest.ensembl.org/xrefs/symbol/...` and `.../lookup/id/...`), a
public database. The `source` column records, per row, whether its
`description` text was returned directly by that API call or, when the API
returned no description, hand-written from well-established literature for
an iconic, frequently published gene (for example rice `Hd1`/`Hd3a`, the
"Green Revolution" gene `SD1`, or maize `vp1`/`su1`/`y1`); such rows are
marked `hand-curated from established literature` rather than attributed to
Ensembl. No value in these files was invented without a stated source.

**Solanum lycopersicum (tomato) is intentionally not included yet.** Despite
attempting verification through Ensembl's REST API, UniProt, NCBI Gene
(eutils), and SGN/ITAG, this session could not confirm an exact
`Solyc##g######`-style locus identifier against a live authoritative source
for any candidate tomato marker gene. Rather than ship an unverified gene ID
in a scientific tool, tomato support was left out of this first version; it
can be added once a specific gene ID has been independently verified.

## Extending this list

To add a species or gene, verify the exact `gene_id` against a live,
authoritative source first (a genome browser or REST API for that species'
reference annotation), record the source, and add a row. Do not add a
`gene_id` from memory alone.
