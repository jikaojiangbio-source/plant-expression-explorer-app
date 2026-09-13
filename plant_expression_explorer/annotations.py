"""Small, hand-curated species gene-annotation reference lookups.

See ``data/annotations/README.md`` for scope, provenance, and the exact
verification method used for every bundled entry. This is never a genome
annotation and is never used in any calculation; it only supplies optional,
display-only context for an exact matching ``gene_id``.
"""

from __future__ import annotations

import csv
import functools
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
