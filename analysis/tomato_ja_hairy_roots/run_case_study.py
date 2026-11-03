"""Reproduce the tomato hairy-root jasmonate case study end to end.

Usage (from the repository root, in the analysis environment)::

    python -m pip install -r analysis/tomato_ja_hairy_roots/requirements.txt
    python analysis/tomato_ja_hairy_roots/run_case_study.py

Steps
-----
1. Download public inputs into ``analysis/tomato_ja_hairy_roots/.cache/``
   and record their SHA-256 checksums (Expression Atlas E-MTAB-11876 raw
   counts and design, the Ensembl Plants release 53 tomato proteome, and five
   Pfam profile HMMs from InterPro).
2. Annotate terpene synthase (TPS), MYB and JAZ genes by searching the
   proteome with the Pfam HMMs at their curated gathering thresholds.
3. Run DESeq2 (pyDESeq2) for jasmonate vs mock, twice: the standard Wald
   test (H0: log2FC = 0), which feeds the app, and a thresholded test
   (H0: |log2FC| <= 1), which defines the "robustly responsive" gene lists.
4. Cross-check fold changes and calls against Expression Atlas's own,
   independent DESeq2 analysis of the same counts.
5. Write app-ready CSVs to ``data/case_studies/tomato_ja_hairy_roots/``,
   full results to ``results/`` and figures to ``figures/``.

The run is deterministic: re-running with the same inputs reproduces every
output byte for byte, except the ``run_utc`` timestamp in provenance.json.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import platform
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
sys.path.insert(0, str(HERE))

import case_study_lib as lib  # noqa: E402

CACHE = HERE / ".cache"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
APP_DATA = REPOSITORY / "data" / "case_studies" / "tomato_ja_hairy_roots"

ATLAS = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/E-MTAB-11876"
ENSEMBL_PEP = (
    "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-53/fasta/"
    "solanum_lycopersicum/pep/Solanum_lycopersicum.SL3.0.pep.all.fa.gz"
)
PFAM_HMM = "https://www.ebi.ac.uk/interpro/api/entry/pfam/{accession}?annotation=hmm"

INPUTS = {
    "raw_counts": f"{ATLAS}/E-MTAB-11876-raw-counts.tsv",
    "configuration": f"{ATLAS}/E-MTAB-11876-configuration.xml",
    "atlas_analytics": f"{ATLAS}/E-MTAB-11876-analytics.tsv",
    "sdrf": f"{ATLAS}/E-MTAB-11876.condensed-sdrf.tsv",
    "idf": f"{ATLAS}/E-MTAB-11876.idf.txt",
    "proteome": ENSEMBL_PEP,
    **{
        f"pfam_{accession}": PFAM_HMM.format(accession=accession)
        for accession in (
            lib.PFAM_TERPENE_SYNTH_N,
            lib.PFAM_TERPENE_SYNTH_C,
            lib.PFAM_MYB,
            lib.PFAM_TIFY,
            lib.PFAM_JAS,
        )
    },
}

MIN_COUNT = 10
MIN_SAMPLES = 3
ALPHA = 0.05
LFC_THRESHOLD = 1.0


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def fetch(name: str, url: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    suffix = ".hmm.gz" if name.startswith("pfam_") else Path(url.split("?")[0]).suffix
    if url.endswith(".fa.gz"):
        suffix = ".fa.gz"
    path = CACHE / f"{name}{suffix}"
    if not path.exists():
        print(f"  downloading {name}")
        request = urllib.request.Request(url, headers={"User-Agent": "plant-expression-explorer case study"})
        with urllib.request.urlopen(request, timeout=300) as response:
            path.write_bytes(response.read())
    return path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# Gene families
# ---------------------------------------------------------------------------


def search_domains(proteome: Path, hmm_paths: dict[str, Path]) -> tuple[pd.DataFrame, dict[str, str]]:
    """hmmsearch every Pfam HMM against the proteome at gathering cutoffs."""

    import pyhmmer

    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.easel.SequenceFile(str(proteome), digital=True, alphabet=alphabet) as handle:
        sequences = handle.read_block()

    protein_to_gene: dict[str, str] = {}
    with gzip.open(proteome, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].split()
                gene = next(f.split(":", 1)[1] for f in fields if f.startswith("gene:"))
                protein_to_gene[fields[0]] = gene

    rows: list[dict[str, object]] = []
    versions: dict[str, str] = {}
    for accession, path in hmm_paths.items():
        with gzip.open(path, "rb") as raw, pyhmmer.plan7.HMMFile(raw) as hmm_file:
            hmm = hmm_file.read()
        versions[accession] = hmm.accession.decode() if isinstance(hmm.accession, bytes) else str(hmm.accession)
        pipeline = pyhmmer.plan7.Pipeline(alphabet, bit_cutoffs="gathering")
        hits = pipeline.search_hmm(hmm, sequences)
        for hit in hits:
            if not hit.included:
                continue
            name = hit.name.decode() if isinstance(hit.name, bytes) else str(hit.name)
            for domain in hit.domains:
                if not domain.included:
                    continue
                rows.append(
                    {
                        "gene_id": protein_to_gene[name],
                        "protein_id": name,
                        "pfam": accession,
                        "env_from": domain.env_from,
                        "env_to": domain.env_to,
                        "score": round(float(domain.score), 2),
                    }
                )
    hits_table = pd.DataFrame(rows).sort_values(["gene_id", "protein_id", "pfam", "env_from"])
    return hits_table.reset_index(drop=True), versions


# ---------------------------------------------------------------------------
# Differential expression
# ---------------------------------------------------------------------------


def run_deseq2(counts: pd.DataFrame, metadata: pd.DataFrame) -> dict[str, object]:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    design = metadata.set_index("sample_id")[["condition"]].copy()
    design["condition"] = pd.Categorical(
        design["condition"], categories=[lib.REFERENCE_CONDITION, lib.TEST_CONDITION]
    )
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts.T,
        metadata=design,
        design="~condition",
        refit_cooks=True,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    coefficients = list(dds.obsm["design_matrix"].columns)
    contrast = ["condition", lib.TEST_CONDITION, lib.REFERENCE_CONDITION]

    standard = DeseqStats(dds, contrast=contrast, alpha=ALPHA, inference=inference, quiet=True)
    standard.summary()
    results = standard.results_df.copy()

    shrink_coefficient = f"condition[T.{lib.TEST_CONDITION}]"
    if shrink_coefficient not in coefficients:
        raise RuntimeError(f"Unexpected design columns {coefficients}.")
    standard.lfc_shrink(coeff=shrink_coefficient)
    results["log2FoldChange_shrunk"] = standard.results_df["log2FoldChange"]

    thresholded = DeseqStats(
        dds,
        contrast=contrast,
        alpha=ALPHA,
        lfc_null=LFC_THRESHOLD,
        alt_hypothesis="greaterAbs",
        inference=inference,
        quiet=True,
    )
    thresholded.summary()
    results["padj_lfc1"] = thresholded.results_df["padj"]
    results["pvalue_lfc1"] = thresholded.results_df["pvalue"]

    size_factors = pd.Series(np.asarray(dds.obs["size_factors"]), index=dds.obs_names)
    return {"results": results, "size_factors": size_factors, "coefficients": coefficients}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args(argv)

    print("1/5 Fetching inputs")
    paths = {name: fetch(name, url) for name, url in INPUTS.items()}
    input_records = {
        name: {"url": INPUTS[name], "sha256": sha256(path), "bytes": path.stat().st_size}
        for name, path in paths.items()
    }

    print("2/5 Annotating gene families with Pfam HMMs")
    hmm_paths = {name.removeprefix("pfam_"): path for name, path in paths.items() if name.startswith("pfam_")}
    domain_hits, pfam_versions = search_domains(paths["proteome"], hmm_paths)
    families = lib.classify_families(domain_hits)
    family_sizes = families["family"].value_counts().to_dict()
    print("  family sizes:", {k: family_sizes[k] for k in sorted(family_sizes)})

    print("3/5 Differential expression (pyDESeq2)")
    contrast = lib.parse_atlas_configuration(paths["configuration"].read_text(encoding="utf-8"))
    metadata = lib.build_sample_metadata(contrast)
    raw = pd.read_csv(paths["raw_counts"], sep="\t")
    counts_all = lib.prepare_counts(raw, metadata)
    counts = lib.filter_low_counts(counts_all, MIN_COUNT, MIN_SAMPLES)
    print(f"  {len(counts_all)} gene models, {len(counts)} pass the expression filter")

    fitted = run_deseq2(counts, metadata)
    results: pd.DataFrame = fitted["results"]
    own_factors = lib.size_factors_median_of_ratios(counts)
    if not np.allclose(own_factors.to_numpy(), fitted["size_factors"].reindex(own_factors.index).to_numpy(), rtol=1e-6):
        raise RuntimeError("Size factors disagree with pyDESeq2; refusing to write outputs.")

    results["regulation_lfc1"] = lib.classify_regulation(results["padj_lfc1"], results["log2FoldChange"], ALPHA)
    names = raw.set_index("Gene ID")["Gene Name"]
    results["atlas_gene_name"] = names.reindex(results.index)
    family_wide = families.groupby("gene_id")["family"].agg(lambda f: "; ".join(sorted(f)))
    results["family"] = family_wide.reindex(results.index)
    results.index.name = "gene_id"

    print("4/5 Cross-checking against the Expression Atlas DESeq2 analysis")
    atlas = pd.read_csv(paths["atlas_analytics"], sep="\t")
    comparison = lib.compare_with_atlas(results[["log2FoldChange", "padj"]], atlas, ALPHA, LFC_THRESHOLD)
    print("  ", comparison)
    pvalue_scale = lib.atlas_pvalue_scale_check(results[["pvalue", "padj"]], atlas)
    print("  ", pvalue_scale)

    print("5/5 Writing outputs")
    RESULTS.mkdir(exist_ok=True)
    APP_DATA.mkdir(parents=True, exist_ok=True)

    tested_regulation = results["regulation_lfc1"]
    summary_table = lib.family_response_summary(tested_regulation, families, family_sizes)
    all_tested = int((tested_regulation != lib.REGULATION_NOT_TESTED).sum())
    all_up = int((tested_regulation == lib.REGULATION_UP).sum())
    all_down = int((tested_regulation == lib.REGULATION_DOWN).sum())
    enrichment = {}
    for row in summary_table.to_dict("records"):
        up, down = int(row[lib.REGULATION_UP]), int(row[lib.REGULATION_DOWN])
        family_tested = up + down + int(row[lib.REGULATION_NS])
        enrichment[row["family"]] = {
            "tested": family_tested,
            "up": up,
            "down": down,
            "p_enriched_up": lib.hypergeometric_enrichment(up, family_tested, all_up, all_tested),
            "p_enriched_down": lib.hypergeometric_enrichment(down, family_tested, all_down, all_tested),
        }

    jaz_reference = pd.read_csv(HERE / "reference" / "chini2017_tomato_jaz.csv")
    jaz_domain_ids = families.loc[families["family"] == lib.FAMILY_JAZ, "gene_id"]
    positive_control = lib.reference_response_table(
        jaz_reference, counts_all, metadata, results, jaz_domain_ids
    )
    positive_control.to_csv(RESULTS / "positive_control_jaz.csv", index=False, float_format="%.6g")

    tested_jaz = positive_control[positive_control["status"].isin([lib.REGULATION_UP, lib.REGULATION_DOWN, lib.REGULATION_NS])]
    jaz_up = int((tested_jaz["status"] == lib.REGULATION_UP).sum())
    jaz_down = int((tested_jaz["status"] == lib.REGULATION_DOWN).sum())
    enrichment["JAZ (Chini et al. 2017 list)"] = {
        "tested": int(len(tested_jaz)), "up": jaz_up, "down": jaz_down,
        "p_enriched_up": lib.hypergeometric_enrichment(jaz_up, len(tested_jaz), all_up, all_tested),
        "p_enriched_down": lib.hypergeometric_enrichment(jaz_down, len(tested_jaz), all_down, all_tested),
    }
    family_rows = [("all tested genes", all_up, all_down, all_tested),
                   ("JAZ (literature list)", jaz_up, jaz_down, int(len(tested_jaz)))]
    for family in (lib.FAMILY_TPS, lib.FAMILY_R2R3_MYB):
        record = enrichment[family]
        family_rows.append((family, record["up"], record["down"], record["tested"]))

    full = results.reset_index()[
        [
            "gene_id", "atlas_gene_name", "family", "baseMean", "log2FoldChange", "lfcSE",
            "stat", "pvalue", "padj", "log2FoldChange_shrunk", "pvalue_lfc1", "padj_lfc1", "regulation_lfc1",
        ]
    ]
    full.to_csv(RESULTS / "deseq2_jasmonate_vs_mock.csv", index=False, float_format="%.6g")
    families.to_csv(RESULTS / "gene_families.csv", index=False)
    domain_hits.to_csv(RESULTS / "pfam_domain_hits.csv", index=False)
    summary_table.to_csv(RESULTS / "family_response_summary.csv", index=False)

    # App-ready tables. Gene order follows the expression filter (genome order).
    expression = lib.log2_normalised_expression(counts, own_factors)
    expression.to_csv(APP_DATA / "expression_matrix.csv", index=False)
    metadata.to_csv(APP_DATA / "metadata.csv", index=False)
    lib.app_de_table(results.loc[counts.index]).to_csv(APP_DATA / "deg_results.csv", index=False)
    annotation = build_annotation(families, names, set(counts.index))
    annotation.to_csv(APP_DATA / "gene_annotation.csv", index=False)

    summary = {
        "gene_models": int(len(counts_all)),
        "genes_tested": int(len(counts)),
        "size_factors": {k: round(float(v), 6) for k, v in own_factors.items()},
        "library_sizes": {k: int(v) for k, v in counts_all.sum().items()},
        "standard_test": {
            "padj_lt_0.05": int((results["padj"] < ALPHA).sum()),
            "padj_missing": int(results["padj"].isna().sum()),
            "pvalue_missing": int(results["pvalue"].isna().sum()),
        },
        "lfc1_test": {"up": all_up, "down": all_down, "tested": all_tested},
        "family_enrichment_lfc1": enrichment,
        "atlas_comparison": comparison,
        "atlas_pvalue_scale_check": pvalue_scale,
        "positive_control_jaz": {
            "literature_genes": int(len(positive_control)),
            "recovered_by_domain_rule": int(positive_control["in_domain_call"].sum()),
            "domain_rule_calls": int(len(jaz_domain_ids)),
            "status_counts": positive_control["status"].value_counts().to_dict(),
        },
        "design_coefficients": fitted["coefficients"],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    provenance = {
        "run_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "python": platform.python_version(),
        "packages": {
            package: importlib.metadata.version(package)
            for package in ("pydeseq2", "pyhmmer", "pandas", "numpy", "scipy", "matplotlib")
        },
        "parameters": {
            "min_count": MIN_COUNT, "min_samples": MIN_SAMPLES, "alpha": ALPHA,
            "lfc_threshold": LFC_THRESHOLD, "design": "~condition",
            "reference_level": lib.REFERENCE_CONDITION,
            "hmm_cutoff": "Pfam gathering (GA) thresholds",
        },
        "pfam_versions": pfam_versions,
        "inputs": input_records,
    }
    (RESULTS / "provenance.json").write_text(json.dumps(provenance, indent=2, sort_keys=True) + "\n")

    if not args.skip_figures:
        import figures

        figures.make_all(
            results, expression, metadata, families, atlas, FIGURES,
            jaz_ids=set(positive_control["gene_id"]) & set(results.index),
            family_rows=family_rows,
        )
    print("done")
    return 0


def build_annotation(families: pd.DataFrame, names: pd.Series, expressed: set[str]) -> pd.DataFrame:
    """Gene-annotation CSV for the app's 'Upload gene annotation CSV' box."""

    source = (
        "Pfam HMM search (pyhmmer, gathering thresholds) of the Ensembl Plants "
        "release 53 tomato proteome; gene name from Expression Atlas where present"
    )
    by_gene = families.groupby("gene_id").agg(
        family=("family", lambda f: "; ".join(sorted(f))),
        evidence=("evidence", lambda e: "; ".join(e)),
    )
    rows = []
    for gene_id, row in by_gene.iterrows():
        if gene_id not in expressed:
            continue
        name = lib.first_present([names.get(gene_id)])
        rows.append(
            {
                "gene_id": gene_id,
                "symbol": name or row["family"],
                "description": f"{row['family']}: {row['evidence']}",
                "source": source,
            }
        )
    named = names[names.notna() & names.index.isin(expressed) & ~names.index.isin(by_gene.index)]
    for gene_id, name in named.items():
        rows.append(
            {
                "gene_id": gene_id,
                "symbol": str(name),
                "description": "Gene name from Expression Atlas (Ensembl Plants 53); no family call",
                "source": "Expression Atlas E-MTAB-11876 raw-counts table",
            }
        )
    return pd.DataFrame(rows).sort_values("gene_id").reset_index(drop=True)


if __name__ == "__main__":
    raise SystemExit(main())
