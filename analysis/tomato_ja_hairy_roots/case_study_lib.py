"""Pure, testable helpers for the tomato jasmonate hairy-root case study.

This module deliberately imports only the Python standard library, numpy and
pandas, so its unit tests run in the application's own pinned CI environment.
Everything that needs the heavier analysis stack (pyDESeq2, pyhmmer,
matplotlib) lives in ``run_case_study.py``.

Conventions used throughout:

* ``counts`` tables are genes x samples, indexed by ``gene_id``, integers.
* Condition labels are ``"mock"`` (reference) and ``"jasmonate"`` (test).
* A positive ``log2FoldChange`` means higher in jasmonate than in mock.
"""

from __future__ import annotations

import math
import xml.etree.ElementTree as ElementTree
from collections.abc import Iterable, Mapping
from dataclasses import dataclass

import numpy as np
import pandas as pd

REFERENCE_CONDITION = "mock"
TEST_CONDITION = "jasmonate"

# Pfam accessions (without version) used for the gene-family calls.
PFAM_TERPENE_SYNTH_N = "PF01397"  # Terpene_synth: terpene synthase, N-terminal domain
PFAM_TERPENE_SYNTH_C = "PF03936"  # Terpene_synth_C: terpene synthase, metal-binding domain
PFAM_MYB = "PF00249"  # Myb_DNA-binding: one MYB repeat
PFAM_TIFY = "PF06200"  # tify domain
PFAM_JAS = "PF09425"  # Jas motif (JAZ degron / MYC-interaction motif)

FAMILY_TPS = "TPS"
FAMILY_R2R3_MYB = "R2R3-MYB"
FAMILY_JAZ = "JAZ"

REGULATION_UP = "up"
REGULATION_DOWN = "down"
REGULATION_NS = "not significant"
REGULATION_NOT_TESTED = "not tested"


# ---------------------------------------------------------------------------
# Experimental design
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AtlasContrast:
    """One two-group contrast as declared in an Expression Atlas configuration."""

    reference_label: str
    test_label: str
    reference_assays: tuple[str, ...]
    test_assays: tuple[str, ...]


def parse_atlas_configuration(xml_text: str) -> AtlasContrast:
    """Parse the single two-group contrast from an Atlas configuration XML.

    Raises ``ValueError`` if the file does not declare exactly one contrast or
    if an assay appears in both groups.
    """

    root = ElementTree.fromstring(xml_text)
    groups: dict[str, tuple[str, tuple[str, ...]]] = {}
    for group in root.iter("assay_group"):
        assays = tuple(
            (assay.text or "").strip() for assay in group.iter("assay")
        )
        groups[group.attrib["id"]] = (group.attrib.get("label", ""), assays)

    contrasts = list(root.iter("contrast"))
    if len(contrasts) != 1:
        raise ValueError(
            f"Expected exactly one contrast, found {len(contrasts)}."
        )
    contrast = contrasts[0]
    reference_id = contrast.findtext("reference_assay_group", "").strip()
    test_id = contrast.findtext("test_assay_group", "").strip()
    if reference_id not in groups or test_id not in groups:
        raise ValueError("Contrast refers to an undeclared assay group.")

    reference_label, reference_assays = groups[reference_id]
    test_label, test_assays = groups[test_id]
    overlap = set(reference_assays) & set(test_assays)
    if overlap:
        raise ValueError(f"Assays appear in both groups: {sorted(overlap)}.")
    return AtlasContrast(
        reference_label=reference_label,
        test_label=test_label,
        reference_assays=tuple(sorted(reference_assays)),
        test_assays=tuple(sorted(test_assays)),
    )


def build_sample_metadata(contrast: AtlasContrast) -> pd.DataFrame:
    """Return app-ready sample metadata, reference group first.

    ``sample_id`` values are human-readable (``mock_1`` ...); the original ENA
    run accession is kept in ``run_accession`` so every sample stays traceable.
    """

    rows = []
    for condition, assays in (
        (REFERENCE_CONDITION, contrast.reference_assays),
        (TEST_CONDITION, contrast.test_assays),
    ):
        for replicate, accession in enumerate(assays, start=1):
            rows.append(
                {
                    "sample_id": f"{condition}_{replicate}",
                    "condition": condition,
                    "run_accession": accession,
                }
            )
    return pd.DataFrame(rows, columns=["sample_id", "condition", "run_accession"])


# ---------------------------------------------------------------------------
# Counts
# ---------------------------------------------------------------------------


def prepare_counts(
    raw_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_id_prefix: str = "Solyc",
) -> pd.DataFrame:
    """Select nuclear gene models and rename run accessions to sample IDs.

    ``raw_counts`` is the Atlas ``*-raw-counts.tsv`` table (``Gene ID``,
    ``Gene Name`` and one integer column per run). Only rows whose ID starts
    with ``gene_id_prefix`` are kept: the Atlas table also contains ~1,200
    Ensembl ncRNA features (tRNA, snoRNA, miRNA) that a poly(A)-selected
    library does not measure meaningfully and that have no protein for the
    family annotation.
    """

    missing = [
        accession
        for accession in metadata["run_accession"]
        if accession not in raw_counts.columns
    ]
    if missing:
        raise ValueError(f"Run accessions missing from counts: {missing}.")

    selected = raw_counts[
        raw_counts["Gene ID"].astype(str).str.startswith(gene_id_prefix)
    ]
    if selected["Gene ID"].duplicated().any():
        raise ValueError("Duplicate gene IDs in the raw counts table.")

    counts = selected.set_index("Gene ID")[list(metadata["run_accession"])]
    counts.columns = list(metadata["sample_id"])
    counts.index.name = "gene_id"

    values = counts.to_numpy()
    if not np.issubdtype(values.dtype, np.integer):
        if not np.all(np.isfinite(values)) or not np.all(values == np.round(values)):
            raise ValueError("Raw counts must be finite integers.")
        counts = counts.astype("int64")
    if (counts < 0).to_numpy().any():
        raise ValueError("Raw counts must be non-negative.")
    return counts


def filter_low_counts(
    counts: pd.DataFrame,
    min_count: int = 10,
    min_samples: int = 3,
) -> pd.DataFrame:
    """Keep genes with at least ``min_count`` reads in ``min_samples`` samples.

    ``min_samples`` defaults to the smallest group size (3), the rule of thumb
    from the DESeq2 vignette, so a gene expressed in only one condition is
    still retained.
    """

    if min_count < 0 or min_samples < 1:
        raise ValueError("min_count must be >= 0 and min_samples >= 1.")
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    return counts.loc[keep]


def size_factors_median_of_ratios(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 median-of-ratios size factors (genes with any zero excluded).

    Implemented here (rather than taken from pyDESeq2) so the normalised
    expression matrix shipped to the app can be recomputed and unit-tested
    without the analysis stack. ``run_case_study.py`` asserts that these agree
    with pyDESeq2's own size factors.
    """

    positive = counts.loc[(counts > 0).all(axis=1)]
    if positive.empty:
        raise ValueError("No gene has a positive count in every sample.")
    log_counts = np.log(positive.to_numpy(dtype=float))
    log_geo_means = log_counts.mean(axis=1, keepdims=True)
    factors = np.exp(np.median(log_counts - log_geo_means, axis=0))
    return pd.Series(factors, index=counts.columns, name="size_factor")


def log2_normalised_expression(
    counts: pd.DataFrame,
    size_factors: pd.Series,
    decimals: int = 6,
) -> pd.DataFrame:
    """Return ``log2(count / size_factor + 1)`` with a ``gene_id`` column."""

    aligned = size_factors.reindex(counts.columns)
    if aligned.isna().any():
        raise ValueError("Every sample needs a size factor.")
    normalised = counts.divide(aligned, axis=1)
    expression = np.log2(normalised + 1.0).round(decimals)
    return expression.reset_index()


# ---------------------------------------------------------------------------
# Gene families from Pfam domain hits
# ---------------------------------------------------------------------------


def myb_repeat_class(repeat_count: int) -> str:
    """Name the MYB class implied by the number of MYB (PF00249) repeats."""

    if repeat_count < 1:
        raise ValueError("A MYB protein needs at least one repeat.")
    if repeat_count == 1:
        return "1R-MYB (MYB-related)"
    if repeat_count == 2:
        return FAMILY_R2R3_MYB
    if repeat_count == 3:
        return "3R-MYB"
    return f"{repeat_count}R-MYB"


def classify_families(domain_hits: pd.DataFrame) -> pd.DataFrame:
    """Call TPS, R2R3-MYB and JAZ families from per-protein domain hits.

    ``domain_hits`` has one row per significant domain hit with columns
    ``gene_id``, ``protein_id`` and ``pfam`` (accession without version).
    Calls are made per protein isoform and a gene inherits the call of its
    best-supported isoform:

    * TPS: at least one terpene-synthase domain (N-terminal PF01397 or
      metal-binding PF03936). ``evidence`` records which were found.
    * MYB class: number of PF00249 repeats on the isoform with the most
      repeats; only exactly two repeats is called R2R3-MYB.
    * JAZ: tify (PF06200) and Jas (PF09425) on the same isoform. tify-only
      proteins (PPD, ZIM/ZML, TIFY8) are deliberately not called JAZ.

    Returns one row per (gene_id, family) with columns ``gene_id``,
    ``family``, ``evidence``.
    """

    required = {"gene_id", "protein_id", "pfam"}
    if not required.issubset(domain_hits.columns):
        raise ValueError(f"domain_hits needs columns {sorted(required)}.")

    per_protein = (
        domain_hits.groupby(["gene_id", "protein_id", "pfam"]).size().unstack(fill_value=0)
    )
    for accession in (
        PFAM_TERPENE_SYNTH_N,
        PFAM_TERPENE_SYNTH_C,
        PFAM_MYB,
        PFAM_TIFY,
        PFAM_JAS,
    ):
        if accession not in per_protein.columns:
            per_protein[accession] = 0
    per_protein = per_protein.reset_index()

    records: list[dict[str, str]] = []

    tps = per_protein[
        (per_protein[PFAM_TERPENE_SYNTH_N] > 0) | (per_protein[PFAM_TERPENE_SYNTH_C] > 0)
    ]
    for gene_id, group in tps.groupby("gene_id"):
        has_n = bool((group[PFAM_TERPENE_SYNTH_N] > 0).any())
        has_c = bool((group[PFAM_TERPENE_SYNTH_C] > 0).any())
        both = bool(
            ((group[PFAM_TERPENE_SYNTH_N] > 0) & (group[PFAM_TERPENE_SYNTH_C] > 0)).any()
        )
        if both:
            evidence = "N-terminal + metal-binding domains (PF01397, PF03936)"
        elif has_n:
            evidence = "N-terminal domain only (PF01397); possibly partial"
        else:
            evidence = "metal-binding domain only (PF03936); possibly partial"
        records.append({"gene_id": gene_id, "family": FAMILY_TPS, "evidence": evidence})

    myb = per_protein[per_protein[PFAM_MYB] > 0]
    for gene_id, group in myb.groupby("gene_id"):
        repeats = int(group[PFAM_MYB].max())
        records.append(
            {
                "gene_id": gene_id,
                "family": myb_repeat_class(repeats),
                "evidence": f"{repeats} MYB repeat(s) (PF00249) on best isoform",
            }
        )

    jaz = per_protein[(per_protein[PFAM_TIFY] > 0) & (per_protein[PFAM_JAS] > 0)]
    for gene_id in sorted(jaz["gene_id"].unique()):
        records.append(
            {
                "gene_id": gene_id,
                "family": FAMILY_JAZ,
                "evidence": "tify (PF06200) + Jas (PF09425) on one isoform",
            }
        )

    return (
        pd.DataFrame(records, columns=["gene_id", "family", "evidence"])
        .sort_values(["family", "gene_id"])
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Differential-expression results
# ---------------------------------------------------------------------------


def app_de_table(results: pd.DataFrame, decimals: int = 8) -> pd.DataFrame:
    """Shape DESeq2-style results into the app's four-column contract.

    ``results`` is indexed by gene ID and has ``log2FoldChange``, ``pvalue``
    and ``padj``. Missing p-values (Cook's outliers, independent filtering)
    are kept as missing, never filled. Values are rounded only for file size;
    rounding never moves a p-value outside [0, 1].
    """

    table = results[["log2FoldChange", "pvalue", "padj"]].copy()
    table["log2FoldChange"] = table["log2FoldChange"].round(decimals)
    for column in ("pvalue", "padj"):
        # Keep tiny p-values informative: use significant figures, not decimals.
        table[column] = table[column].map(_round_significant)
    table.index.name = "gene_id"
    return table.reset_index()


def _round_significant(value: float, figures: int = 6) -> float:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return float("nan")
    if value == 0:
        return 0.0
    rounded = float(f"{value:.{figures}g}")
    return min(max(rounded, 0.0), 1.0)


def classify_regulation(
    padj: pd.Series,
    log2_fold_change: pd.Series,
    alpha: float = 0.05,
) -> pd.Series:
    """Label each gene up / down / not significant / not tested.

    ``padj`` here should come from a test that already includes any effect
    size threshold (for example DESeq2's ``lfcThreshold``), so this function
    applies only the FDR cutoff and the sign; it never re-thresholds the fold
    change post hoc.
    """

    labels = pd.Series(REGULATION_NS, index=padj.index, dtype="object")
    tested = padj.notna()
    significant = tested & (padj < alpha)
    labels[~tested] = REGULATION_NOT_TESTED
    labels[significant & (log2_fold_change > 0)] = REGULATION_UP
    labels[significant & (log2_fold_change < 0)] = REGULATION_DOWN
    return labels


def family_response_summary(
    regulation: pd.Series,
    families: pd.DataFrame,
    genes_in_genome: Mapping[str, int],
) -> pd.DataFrame:
    """Count family members by regulation class.

    ``regulation`` is indexed by gene ID for the genes that entered the test.
    ``genes_in_genome`` gives each family's size in the whole proteome, so
    the table also shows how many members were too lowly expressed to test.
    """

    rows = []
    for family, group in families.groupby("family"):
        members = set(group["gene_id"])
        tested = regulation.reindex(sorted(members)).dropna()
        counts = tested.value_counts()
        rows.append(
            {
                "family": family,
                "in_proteome": int(genes_in_genome.get(family, len(members))),
                "expressed_and_tested": int(len(tested)),
                REGULATION_UP: int(counts.get(REGULATION_UP, 0)),
                REGULATION_DOWN: int(counts.get(REGULATION_DOWN, 0)),
                REGULATION_NS: int(counts.get(REGULATION_NS, 0)),
                REGULATION_NOT_TESTED: int(counts.get(REGULATION_NOT_TESTED, 0)),
            }
        )
    return pd.DataFrame(rows).sort_values("family").reset_index(drop=True)


def hypergeometric_enrichment(
    family_hits: int,
    family_tested: int,
    all_hits: int,
    all_tested: int,
) -> float:
    """One-sided P(X >= family_hits) for a family among tested genes.

    Exact hypergeometric upper tail computed with log-gamma so the module
    needs no SciPy. Used to ask whether e.g. JA-induced genes are
    over-represented among expressed TPS genes relative to all expressed
    genes.
    """

    if not (0 <= family_hits <= family_tested <= all_tested):
        raise ValueError("Require 0 <= family_hits <= family_tested <= all_tested.")
    if not (family_hits <= all_hits <= all_tested):
        raise ValueError("Require family_hits <= all_hits <= all_tested.")

    def log_comb(n: int, k: int) -> float:
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    upper = min(family_tested, all_hits)
    log_total = log_comb(all_tested, family_tested)
    terms = [
        log_comb(all_hits, k)
        + log_comb(all_tested - all_hits, family_tested - k)
        - log_total
        for k in range(family_hits, upper + 1)
        if family_tested - k <= all_tested - all_hits
    ]
    if not terms:
        return 0.0
    peak = max(terms)
    return float(min(1.0, math.exp(peak) * sum(math.exp(t - peak) for t in terms)))


def strip_version(gene_id: str) -> str:
    """``Solyc12g009220.2`` -> ``Solyc12g009220`` (literature IDs are unversioned)."""

    head, dot, tail = str(gene_id).rpartition(".")
    return head if dot and tail.isdigit() else str(gene_id)


def match_unversioned_ids(
    unversioned: Iterable[str],
    versioned: Iterable[str],
) -> dict[str, str | None]:
    """Map each literature (unversioned) ID to the one annotation ID it denotes.

    Raises ``ValueError`` if an unversioned ID matches more than one versioned
    ID, because silently picking one would hide an annotation ambiguity.
    """

    by_base: dict[str, list[str]] = {}
    for gene_id in versioned:
        by_base.setdefault(strip_version(gene_id), []).append(str(gene_id))
    mapping: dict[str, str | None] = {}
    for gene_id in unversioned:
        matches = by_base.get(str(gene_id), [])
        if len(matches) > 1:
            raise ValueError(f"{gene_id} matches several gene models: {matches}.")
        mapping[str(gene_id)] = matches[0] if matches else None
    return mapping


def reference_response_table(
    reference: pd.DataFrame,
    counts_all: pd.DataFrame,
    metadata: pd.DataFrame,
    results: pd.DataFrame,
    domain_family_ids: Iterable[str] = (),
) -> pd.DataFrame:
    """Per-gene response of a literature-defined gene set (e.g. positive controls).

    ``reference`` needs ``symbol`` and unversioned ``gene_id``. Genes that were
    removed by the expression filter are reported with their raw mean counts
    and ``status = "below expression filter"`` instead of disappearing.
    """

    mapping = match_unversioned_ids(reference["gene_id"], counts_all.index)
    domain_family_ids = set(domain_family_ids)
    mock = list(metadata.loc[metadata["condition"] == REFERENCE_CONDITION, "sample_id"])
    test = list(metadata.loc[metadata["condition"] == TEST_CONDITION, "sample_id"])
    rows = []
    for record in reference.to_dict("records"):
        versioned = mapping[str(record["gene_id"])]
        row = {
            "symbol": record["symbol"],
            "gene_id": versioned or record["gene_id"],
            "in_domain_call": bool(versioned in domain_family_ids) if versioned else False,
            "mean_count_mock": np.nan,
            "mean_count_jasmonate": np.nan,
            "log2FoldChange": np.nan,
            "padj_lfc1": np.nan,
            "status": "absent from annotation",
        }
        if versioned is not None:
            row["mean_count_mock"] = round(float(counts_all.loc[versioned, mock].mean()), 2)
            row["mean_count_jasmonate"] = round(float(counts_all.loc[versioned, test].mean()), 2)
            if versioned in results.index:
                row["log2FoldChange"] = float(results.loc[versioned, "log2FoldChange"])
                row["padj_lfc1"] = float(results.loc[versioned, "padj_lfc1"])
                row["status"] = str(results.loc[versioned, "regulation_lfc1"])
            else:
                row["status"] = "below expression filter"
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Cross-check against the Expression Atlas DESeq2 analysis
# ---------------------------------------------------------------------------


def compare_with_atlas(
    ours: pd.DataFrame,
    atlas: pd.DataFrame,
    alpha: float = 0.05,
    min_abs_lfc: float = 1.0,
) -> dict[str, float | int]:
    """Agreement between this analysis and Atlas's independent DESeq2 run.

    ``ours`` is indexed by gene ID with ``log2FoldChange`` and ``padj``.
    ``atlas`` is the Atlas analytics table with ``Gene ID``, a
    ``*.p-value`` column (treated as Benjamini-Hochberg adjusted; verified
    empirically by :func:`atlas_pvalue_scale_check`)
    and a ``*.log2foldchange`` column rounded to 0.1. Both sides are
    classified with the same rule (padj < alpha and |log2FC| >= min_abs_lfc)
    so the comparison is like-for-like.
    """

    pcol = _single_column(atlas, ".p-value")
    lcol = _single_column(atlas, ".log2foldchange")
    theirs = atlas.set_index("Gene ID")[[lcol, pcol]].rename(
        columns={lcol: "log2FoldChange", pcol: "padj"}
    )
    shared = ours.index.intersection(theirs.index)
    a = ours.loc[shared]
    b = theirs.loc[shared]
    both_lfc = a["log2FoldChange"].notna() & b["log2FoldChange"].notna()
    spearman = float(
        a.loc[both_lfc, "log2FoldChange"].rank().corr(
            b.loc[both_lfc, "log2FoldChange"].rank()
        )
    )
    pearson = float(
        a.loc[both_lfc, "log2FoldChange"].corr(b.loc[both_lfc, "log2FoldChange"])
    )

    def calls(frame: pd.DataFrame) -> set[str]:
        hit = (frame["padj"] < alpha) & (frame["log2FoldChange"].abs() >= min_abs_lfc)
        return set(frame.index[hit.fillna(False)])

    ours_calls = calls(a)
    atlas_calls = calls(b)
    union = ours_calls | atlas_calls
    return {
        "shared_genes": int(len(shared)),
        "genes_with_both_lfc": int(both_lfc.sum()),
        "lfc_pearson": round(pearson, 4),
        "lfc_spearman": round(spearman, 4),
        "deg_ours": int(len(ours_calls)),
        "deg_atlas": int(len(atlas_calls)),
        "deg_both": int(len(ours_calls & atlas_calls)),
        "deg_jaccard": round(len(ours_calls & atlas_calls) / len(union), 4)
        if union
        else 1.0,
    }


def atlas_pvalue_scale_check(ours: pd.DataFrame, atlas: pd.DataFrame) -> dict[str, float | int]:
    """Is the Atlas ``*.p-value`` column raw or adjusted? Check empirically.

    Compares it, on a log10 scale, with this analysis's raw ``pvalue`` and its
    adjusted ``padj``. An adjusted column sits closer to ``padj`` and is
    almost always >= the raw p-value.
    """

    pcol = _single_column(atlas, ".p-value")
    theirs = atlas.set_index("Gene ID")[pcol].rename("atlas")
    joined = ours[["pvalue", "padj"]].join(theirs, how="inner").dropna()
    joined = joined[(joined > 0).all(axis=1)]
    log_atlas = np.log10(joined["atlas"])
    return {
        "genes": int(len(joined)),
        "median_abs_log10_diff_vs_pvalue": round(float(np.median(np.abs(log_atlas - np.log10(joined["pvalue"])))), 4),
        "median_abs_log10_diff_vs_padj": round(float(np.median(np.abs(log_atlas - np.log10(joined["padj"])))), 4),
        "fraction_atlas_ge_pvalue": round(float((joined["atlas"] >= joined["pvalue"]).mean()), 4),
    }


def _single_column(table: pd.DataFrame, suffix: str) -> str:
    matches = [column for column in table.columns if column.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"Expected one column ending in {suffix!r}, got {matches}.")
    return matches[0]


def first_present(values: Iterable[object]) -> str | None:
    """Return the first non-missing, non-blank string in ``values``."""

    for value in values:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            continue
        text = str(value).strip()
        if text:
            return text
    return None
