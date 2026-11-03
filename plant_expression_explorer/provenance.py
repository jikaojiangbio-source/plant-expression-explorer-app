"""Descriptive dataset context that is never used in computation."""

from __future__ import annotations

from dataclasses import dataclass, fields


NOT_SUPPLIED = "Not supplied"

PROVENANCE_FIELD_LABELS: tuple[tuple[str, str], ...] = (
    ("dataset_title", "Dataset title"),
    ("organism", "Organism / taxon"),
    (
        "expression_scale_description",
        "Expression scale or preprocessing description",
    ),
    ("upstream_normalization_method", "Upstream normalization method"),
    (
        "reference_genome_annotation",
        "Reference genome / annotation release",
    ),
    ("feature_level", "Feature level"),
    (
        "de_contrast_description",
        "Differential-expression contrast description",
    ),
    ("notes", "Notes"),
)


@dataclass(frozen=True)
class DatasetProvenance:
    """Verbatim descriptive context supplied with a dataset.

    Values are preserved exactly as provided, including empty strings,
    whitespace, punctuation, and Unicode. They are not parsed, normalized,
    scientifically verified, or used by any numerical computation.
    """

    dataset_title: str | None = None
    organism: str | None = None
    expression_scale_description: str | None = None
    upstream_normalization_method: str | None = None
    reference_genome_annotation: str | None = None
    feature_level: str | None = None
    de_contrast_description: str | None = None
    notes: str | None = None

    def __post_init__(self) -> None:
        for field in fields(self):
            value = getattr(self, field.name)
            if value is not None and not isinstance(value, str):
                raise TypeError(
                    f"Dataset provenance field '{field.name}' must be a string "
                    "or None."
                )


def provenance_display_rows(
    provenance: DatasetProvenance | None,
) -> tuple[tuple[str, str], ...]:
    """Return stable label/value rows without rewriting supplied text."""

    if provenance is not None and not isinstance(provenance, DatasetProvenance):
        raise TypeError("provenance must be a DatasetProvenance instance or None.")

    return tuple(
        (
            label,
            NOT_SUPPLIED
            if provenance is None or getattr(provenance, attribute) is None
            else getattr(provenance, attribute),
        )
        for attribute, label in PROVENANCE_FIELD_LABELS
    )


DEMO_PROVENANCE = DatasetProvenance(
    dataset_title="Synthetic tomato-style demonstration data",
    organism=(
        "Synthetic tomato-style example; not observations from a real organism"
    ),
    expression_scale_description=(
        "Synthetic log2-normalised expression-like values; not raw counts"
    ),
    upstream_normalization_method=(
        "Not applicable: values were constructed directly; raw counts were not normalized"
    ),
    reference_genome_annotation=(
        "Not applicable: all supplied gene identifiers are fictional"
    ),
    feature_level="Fictional gene-level identifiers for software demonstration",
    de_contrast_description=(
        "Constructed High_nitrate minus Control arithmetic-mean difference; "
        "not a fitted coefficient"
    ),
    notes=(
        "Constructed p-values do not measure uncertainty and support no tomato "
        "biological conclusions"
    ),
)


CASE_STUDY_PROVENANCE = DatasetProvenance(
    dataset_title=(
        "Tomato (cv. Moneymaker) hairy roots expressing a GUS control construct, "
        "50 µM jasmonate vs mock for 24 h, 3 biological replicates each "
        "(Expression Atlas / ArrayExpress E-MTAB-11876)"
    ),
    organism="Solanum lycopersicum",
    expression_scale_description=(
        "log2(DESeq2 median-of-ratios normalised count + 1), computed offline "
        "from Expression Atlas raw counts"
    ),
    upstream_normalization_method=(
        "Expression Atlas iRAP pipeline (TopHat2 alignment, HTSeq counts); "
        "DESeq2 median-of-ratios size factors in analysis/tomato_ja_hairy_roots"
    ),
    reference_genome_annotation="Ensembl Plants release 53, tomato assembly SL3.0",
    feature_level=(
        "Gene level; nuclear Solyc gene models with >= 10 reads in >= 3 samples "
        "(19,522 of 34,658)"
    ),
    de_contrast_description=(
        "pyDESeq2 0.5.4 Wald test, design ~condition, jasmonate vs mock "
        "(reference); unshrunk maximum-likelihood log2 fold changes; "
        "Benjamini-Hochberg padj with DESeq2 independent filtering"
    ),
    notes=(
        "Real public data. Statistics were computed offline, not by this "
        "application; see analysis/tomato_ja_hairy_roots/README.md for methods, "
        "cross-checks against Expression Atlas and limitations. Gene-family "
        "labels: upload data/case_studies/tomato_ja_hairy_roots/gene_annotation.csv "
        "on the Gene Expression page."
    ),
)
