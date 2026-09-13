"""Small, hand-curated species gene-annotation reference lookups.

See ``data/annotations/README.md`` for scope, provenance, and the exact
verification method used for every bundled entry. This is never a genome
annotation and is never used in any calculation; it only supplies optional,
display-only context for an exact matching ``gene_id``.
"""

from __future__ import annotations

import csv
import functools
import re
from dataclasses import dataclass
from pathlib import Path

ANNOTATIONS_DIRECTORY = Path(__file__).resolve().parents[1] / "data" / "annotations"

SPECIES_LABELS: tuple[tuple[str, str], ...] = (
    ("arabidopsis_thaliana", "Arabidopsis thaliana"),
    ("oryza_sativa", "Oryza sativa (rice)"),
    ("zea_mays", "Zea mays (maize)"),
    ("glycine_max", "Glycine max (soybean)"),
)
_LABEL_TO_KEY = {label: key for key, label in SPECIES_LABELS}

# Each entry recognizes one *other*, non-interchangeable identifier system
# documented for that species, distinct from the system used by this
# application's bundled reference list (see data/annotations/README.md).
# These patterns are used only to write a more specific "no match" message;
# they never feed a lookup, so no fuzzy or cross-system matching occurs.
_ALTERNATE_IDENTIFIER_SYSTEMS: dict[str, tuple[tuple[re.Pattern[str], str], ...]] = {
    "arabidopsis_thaliana": (
        (
            re.compile(r"^AT[1-5MC]G\d{5}\.\d+$", re.IGNORECASE),
            "This looks like a TAIR locus ID with a transcript/splice-variant "
            "version suffix (e.g. '.1'). The bundled reference list is keyed "
            "by the bare locus ID without a suffix; that suffix is not "
            "stripped automatically.",
        ),
    ),
    "oryza_sativa": (
        (
            re.compile(r"^LOC_Os\d{2}g\d+$", re.IGNORECASE),
            "This looks like an MSU/TIGR locus ID (LOC_OsNNgNNNNN). The "
            "bundled reference list is keyed by RAP-DB IDs (e.g. "
            "Os01g0883800) instead; the two naming systems both describe "
            "the Nipponbare genome but are not interchangeable strings.",
        ),
    ),
    "zea_mays": (
        (
            re.compile(r"^GRMZM\d+G\d+$", re.IGNORECASE),
            "This looks like an older GRMZM-style (AGPv2/v3) maize gene ID. "
            "The bundled reference list is keyed by Ensembl Plants "
            "Zm00001eb-style IDs (B73 RefGen_v5) instead; the same gene can "
            "have different identifiers across genome-assembly versions.",
        ),
        (
            re.compile(r"^Zm00001d\d+$", re.IGNORECASE),
            "This looks like a Zm00001d-style (B73 RefGen_v4) maize gene "
            "ID. The bundled reference list is keyed by the newer "
            "Zm00001eb-style IDs (B73 RefGen_v5) instead; the same gene can "
            "have different identifiers across genome-assembly versions.",
        ),
    ),
    "glycine_max": (
        (
            re.compile(r"^Glyma\.\d+G\d+$", re.IGNORECASE),
            "This looks like a Wm82.a2 dot-notation soybean gene ID "
            "(Glyma.10G246300 style). The bundled reference list is keyed "
            "by the underscore-notation form (GLYMA_10G246300 style) "
            "instead; both describe the same assembly but are different "
            "strings.",
        ),
    ),
}


def identifier_format_hint(species_label: str, gene_id: object) -> str | None:
    """Return a note about a *different* documented ID system, or ``None``.

    This never performs a lookup with a rewritten identifier: it only
    recognizes the shape of a well-documented alternate identifier system
    for the species, to explain a likely reason an exact match failed.
    Returns ``None`` when the species is unrecognised or the identifier's
    shape does not match any recognized alternate system.
    """

    patterns = _ALTERNATE_IDENTIFIER_SYSTEMS.get(_LABEL_TO_KEY.get(species_label, ""))
    if patterns is None:
        return None
    text = str(gene_id)
    for pattern, hint in patterns:
        if pattern.match(text):
            return hint
    return None


@dataclass(frozen=True)
class GeneAnnotation:
    """One verified gene_id/symbol/description row for a supported species."""

    gene_id: str
    symbol: str
    description: str
    source: str


def list_supported_species() -> tuple[str, ...]:
    """Return supported species display labels in a stable, fixed order."""

    return tuple(label for _key, label in SPECIES_LABELS)


def lookup_gene_annotation(
    species_label: str,
    gene_id: object,
) -> GeneAnnotation | None:
    """Return the bundled annotation for one exact gene ID, or ``None``.

    Matching is exact ``str(value)`` equality only, consistent with every
    other identifier match in this application: no case-folding, trimming,
    or alias resolution. Returns ``None`` when the species is unrecognised
    or the gene ID has no entry in that species' small reference list.
    """

    species_key = _LABEL_TO_KEY.get(species_label)
    if species_key is None:
        return None
    table = _load_species_table(species_key)
    return table.get(str(gene_id))


@functools.lru_cache(maxsize=None)
def _load_species_table(species_key: str) -> dict[str, GeneAnnotation]:
    path = ANNOTATIONS_DIRECTORY / f"{species_key}.csv"
    table: dict[str, GeneAnnotation] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            table[row["gene_id"]] = GeneAnnotation(
                gene_id=row["gene_id"],
                symbol=row["symbol"],
                description=row["description"],
                source=row["source"],
            )
    return table
