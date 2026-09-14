"""Regression tests for descriptive dataset provenance."""

from dataclasses import FrozenInstanceError

import pytest

from plant_expression_explorer.provenance import (
    DEMO_PROVENANCE,
    NOT_SUPPLIED,
    PROVENANCE_FIELD_LABELS,
    DatasetProvenance,
    provenance_display_rows,
)


def test_provenance_preserves_every_supplied_string_exactly() -> None:
    provenance = DatasetProvenance(
        dataset_title="  盐胁迫 RNA-seq  ",
        organism="Solanum lycopersicum 🍅",
        expression_scale_description="log2(TPM + 1)\nrelease value",
        upstream_normalization_method="",
        reference_genome_annotation="SL4.0 / ITAG4.1",
        feature_level="gene   model",
        de_contrast_description="treated - control",
        notes="\tkeep leading tab and trailing spaces  ",
    )

    rows = dict(provenance_display_rows(provenance))

    for attribute, label in PROVENANCE_FIELD_LABELS:
        assert rows[label] == getattr(provenance, attribute)


def test_missing_provenance_displays_not_supplied_for_every_field() -> None:
    rows = provenance_display_rows(None)

    assert tuple(label for _, label in PROVENANCE_FIELD_LABELS) == tuple(
        label for label, _ in rows
    )
    assert {value for _, value in rows} == {NOT_SUPPLIED}


@pytest.mark.parametrize("value", [True, 1, 1.5, object()])
def test_provenance_rejects_non_string_values_without_coercion(value: object) -> None:
    with pytest.raises(TypeError, match="organism"):
        DatasetProvenance(organism=value)  # type: ignore[arg-type]


def test_provenance_is_frozen() -> None:
    provenance = DatasetProvenance(dataset_title="unchanged")

    with pytest.raises(FrozenInstanceError):
        provenance.dataset_title = "changed"  # type: ignore[misc]


def test_display_rows_reject_an_invalid_container() -> None:
    with pytest.raises(TypeError, match="DatasetProvenance"):
        provenance_display_rows({"organism": "tomato"})  # type: ignore[arg-type]


def test_demo_provenance_discloses_synthetic_scale_and_contrast() -> None:
    assert "Synthetic" in (DEMO_PROVENANCE.dataset_title or "")
    assert "not raw counts" in (
        DEMO_PROVENANCE.expression_scale_description or ""
    )
    assert "not a fitted coefficient" in (
        DEMO_PROVENANCE.de_contrast_description or ""
    )
    assert "fictional" in (DEMO_PROVENANCE.reference_genome_annotation or "")
