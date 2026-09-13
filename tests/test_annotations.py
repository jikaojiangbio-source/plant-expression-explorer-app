"""Tests for the small, hand-curated species gene-annotation reference."""

import csv

from plant_expression_explorer.annotations import (
    ANNOTATIONS_DIRECTORY,
    GeneAnnotation,
    list_supported_species,
    lookup_gene_annotation,
)


def test_list_supported_species_returns_fixed_order() -> None:
    assert list_supported_species() == (
        "Arabidopsis thaliana",
        "Oryza sativa (rice)",
        "Zea mays (maize)",
        "Glycine max (soybean)",
    )


def test_lookup_returns_a_known_arabidopsis_gene() -> None:
    result = lookup_gene_annotation("Arabidopsis thaliana", "AT1G65480")

    assert result == GeneAnnotation(
        gene_id="AT1G65480",
        symbol="FT",
        description="PEBP (phosphatidylethanolamine-binding protein) family protein",
        source="Ensembl Plants REST API (rest.ensembl.org)",
    )


def test_lookup_matches_exact_string_only() -> None:
    assert lookup_gene_annotation("Arabidopsis thaliana", "at1g65480") is None
    assert lookup_gene_annotation("Arabidopsis thaliana", " AT1G65480 ") is None
    assert lookup_gene_annotation("Arabidopsis thaliana", "AT1G65480") is not None


def test_lookup_returns_none_for_an_unrecognised_gene_id() -> None:
    assert lookup_gene_annotation("Arabidopsis thaliana", "NOT_A_REAL_GENE") is None


def test_lookup_returns_none_for_an_unrecognised_species() -> None:
    assert lookup_gene_annotation("Solanum lycopersicum", "Solyc05g012020") is None
    assert lookup_gene_annotation("", "AT1G65480") is None


def test_lookup_accepts_non_string_gene_id_values() -> None:
    assert lookup_gene_annotation("Arabidopsis thaliana", 12345) is None


def test_every_bundled_species_file_matches_the_documented_contract() -> None:
    for path in sorted(ANNOTATIONS_DIRECTORY.glob("*.csv")):
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            assert reader.fieldnames == ["gene_id", "symbol", "description", "source"]
            rows = list(reader)
        assert rows, f"{path} must contain at least one row"
        seen_ids = set()
        for row in rows:
            assert row["gene_id"].strip() == row["gene_id"]
            assert row["gene_id"] not in seen_ids, f"duplicate gene_id in {path}"
            seen_ids.add(row["gene_id"])
            assert row["symbol"].strip()
            assert row["source"].strip() in (
                "Ensembl Plants REST API (rest.ensembl.org)",
                "hand-curated from established literature",
            )


def test_every_bundled_species_key_has_a_documented_label() -> None:
    from plant_expression_explorer.annotations import SPECIES_LABELS

    bundled_keys = {path.stem for path in ANNOTATIONS_DIRECTORY.glob("*.csv")}
    documented_keys = {key for key, _label in SPECIES_LABELS}
    assert bundled_keys == documented_keys
