"""Tests for the small, hand-curated species gene-annotation reference."""

import csv

from plant_expression_explorer.annotations import (
    ANNOTATIONS_DIRECTORY,
    SPECIES_LABELS,
    GeneAnnotation,
    identifier_format_hint,
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


def test_identifier_format_hint_recognizes_versioned_arabidopsis_ids() -> None:
    hint = identifier_format_hint("Arabidopsis thaliana", "AT1G65480.1")
    assert hint is not None
    assert "version suffix" in hint


def test_identifier_format_hint_recognizes_msu_rice_ids() -> None:
    hint = identifier_format_hint("Oryza sativa (rice)", "LOC_Os01g67980")
    assert hint is not None
    assert "MSU/TIGR" in hint
    assert "RAP-DB" in hint


def test_identifier_format_hint_recognizes_older_maize_ids() -> None:
    grmzm_hint = identifier_format_hint("Zea mays (maize)", "GRMZM2G017053")
    assert grmzm_hint is not None
    assert "GRMZM" in grmzm_hint

    zm00001d_hint = identifier_format_hint("Zea mays (maize)", "Zm00001d048373")
    assert zm00001d_hint is not None
    assert "RefGen_v4" in zm00001d_hint


def test_identifier_format_hint_recognizes_dot_notation_soybean_ids() -> None:
    hint = identifier_format_hint("Glycine max (soybean)", "Glyma.10G246300")
    assert hint is not None
    assert "dot-notation" in hint


def test_identifier_format_hint_returns_none_for_an_unrecognized_shape() -> None:
    assert identifier_format_hint("Arabidopsis thaliana", "AT1G65480") is None
    assert identifier_format_hint("Arabidopsis thaliana", "totally-unrelated") is None


def test_identifier_format_hint_returns_none_for_an_unrecognised_species() -> None:
    assert identifier_format_hint("Solanum lycopersicum", "Solyc05g012020") is None
    assert identifier_format_hint("", "AT1G65480.1") is None


def test_identifier_format_hint_accepts_non_string_gene_id_values() -> None:
    assert identifier_format_hint("Arabidopsis thaliana", 12345) is None


def test_identifier_format_hint_never_matches_a_bundled_gene_id() -> None:
    # A format hint should only ever fire when the exact lookup already
    # failed; bundled IDs themselves must never also match an "alternate
    # system" pattern, which would make the two messages contradict.
    for key, label in SPECIES_LABELS:
        path = ANNOTATIONS_DIRECTORY / f"{key}.csv"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                assert identifier_format_hint(label, row["gene_id"]) is None


def test_every_bundled_species_key_has_a_documented_label() -> None:
    bundled_keys = {path.stem for path in ANNOTATIONS_DIRECTORY.glob("*.csv")}
    documented_keys = {key for key, _label in SPECIES_LABELS}
    assert bundled_keys == documented_keys
