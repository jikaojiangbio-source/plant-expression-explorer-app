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
from enum import StrEnum
from pathlib import Path

import pandas as pd

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


class CustomAnnotationErrorReason(StrEnum):
    """Stable reasons a user-supplied annotation table cannot be indexed."""

    INVALID_TABLE = "INVALID_TABLE"
    EMPTY_TABLE = "EMPTY_TABLE"
    MISSING_REQUIRED_COLUMN = "MISSING_REQUIRED_COLUMN"
    MISSING_REQUIRED_VALUE = "MISSING_REQUIRED_VALUE"
    DUPLICATE_GENE_ID = "DUPLICATE_GENE_ID"


class CustomAnnotationError(ValueError):
    """Expected failure when a user-supplied annotation table is malformed."""

    def __init__(self, reason: CustomAnnotationErrorReason, message: str) -> None:
        self.reason = reason
        super().__init__(message)


CUSTOM_ANNOTATION_REQUIRED_COLUMNS = ("gene_id", "symbol", "description")
CUSTOM_ANNOTATION_SOURCE_LABEL = (
    "user-uploaded annotation file (not independently verified by this "
    "application)"
)


def parse_custom_annotation_table(table: object) -> dict[str, GeneAnnotation]:
    """Validate and index a user-supplied gene-annotation table by exact ID.

    Required columns: ``gene_id``, ``symbol``, ``description``. An optional
    ``source`` column is carried verbatim when present and non-blank for a
    row; otherwise that row's source is reported as
    :data:`CUSTOM_ANNOTATION_SOURCE_LABEL`. Gene identifiers must be
    non-empty and unique, matched later by the same exact ``str(value)``
    equality as every other identifier lookup in this application.

    Unlike the bundled per-species reference lists, this application never
    independently verifies a user-supplied annotation's accuracy; the
    uploader is responsible for its contents, exactly as for the expression,
    metadata, and differential-expression tables.
    """

    if not isinstance(table, pd.DataFrame):
        raise CustomAnnotationError(
            CustomAnnotationErrorReason.INVALID_TABLE,
            "The annotation table must be a pandas DataFrame.",
        )
    if len(table.index) == 0:
        raise CustomAnnotationError(
            CustomAnnotationErrorReason.EMPTY_TABLE,
            "The annotation table contains no rows.",
        )
    missing_columns = [
        column
        for column in CUSTOM_ANNOTATION_REQUIRED_COLUMNS
        if column not in table.columns
    ]
    if missing_columns:
        raise CustomAnnotationError(
            CustomAnnotationErrorReason.MISSING_REQUIRED_COLUMN,
            "The annotation table is missing required column(s): "
            + ", ".join(missing_columns) + ".",
        )

    has_source_column = "source" in table.columns
    result: dict[str, GeneAnnotation] = {}
    seen_ids: set[str] = set()
    duplicate_ids: list[str] = []
    for _, row in table.iterrows():
        gene_id = row["gene_id"]
        symbol = row["symbol"]
        description = row["description"]
        if _is_missing_or_blank(gene_id) or _is_missing_or_blank(symbol):
            raise CustomAnnotationError(
                CustomAnnotationErrorReason.MISSING_REQUIRED_VALUE,
                "Every row must supply a non-blank 'gene_id' and 'symbol'.",
            )
        gene_id_text = str(gene_id)
        if gene_id_text in seen_ids:
            duplicate_ids.append(gene_id_text)
            continue
        seen_ids.add(gene_id_text)
        source_value = row["source"] if has_source_column else None
        source_text = (
            CUSTOM_ANNOTATION_SOURCE_LABEL
            if _is_missing_or_blank(source_value)
            else str(source_value)
        )
        result[gene_id_text] = GeneAnnotation(
            gene_id=gene_id_text,
            symbol=str(symbol),
            description=(
                "" if _is_missing_or_blank(description) else str(description)
            ),
            source=source_text,
        )
    if duplicate_ids:
        raise CustomAnnotationError(
            CustomAnnotationErrorReason.DUPLICATE_GENE_ID,
            "The annotation table contains duplicate gene_id value(s): "
            + ", ".join(sorted(set(duplicate_ids))) + ".",
        )
    return result


def _is_missing_or_blank(value: object) -> bool:
    if isinstance(value, str):
        return not value.strip()
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False
