"""Regression tests for deterministic descriptive PDF report assembly."""

from datetime import datetime, timezone

import pandas as pd
import pytest
from pypdf import PdfReader

import plant_expression_explorer as package
from plant_expression_explorer.report import (
    PDF_MEDIA_TYPE,
    ReportArtifact,
    ReportExportError,
    ReportExportErrorReason,
    ReportSection,
    build_report_pdf,
)

_GENERATED_AT = datetime(2026, 9, 13, 12, 0, 0, tzinfo=timezone.utc)


def _capture_error(**overrides: object) -> ReportExportError:
    kwargs = {
        "title": "Plant Expression Explorer report",
        "generated_at": _GENERATED_AT,
        "intro_lines": ("An intro line.",),
        "sections": (
            ReportSection(title="Section one", table=pd.DataFrame({"a": [1]})),
        ),
        "filename": "report.pdf",
    }
    kwargs.update(overrides)
    with pytest.raises(ReportExportError) as captured:
        build_report_pdf(**kwargs)
    return captured.value


def _extract_text(data: bytes) -> str:
    reader = PdfReader(__import__("io").BytesIO(data))
    return "\n".join(page.extract_text() for page in reader.pages)


def test_package_exports_the_report_public_api() -> None:
    assert package.PDF_MEDIA_TYPE == PDF_MEDIA_TYPE
    assert package.ReportArtifact is ReportArtifact
    assert package.ReportExportError is ReportExportError
    assert package.ReportExportErrorReason is ReportExportErrorReason
    assert package.ReportSection is ReportSection
    assert package.build_report_pdf is build_report_pdf


def test_error_reasons_are_stable_and_unique() -> None:
    values = [reason.value for reason in ReportExportErrorReason]
    assert len(values) == len(set(values))


def test_build_report_pdf_returns_a_well_formed_pdf_with_expected_sections() -> None:
    qc_table = pd.DataFrame(
        {
            "condition": ["control", "treated"],
            "sample_count": [3, 3],
            "sample_ids": [("s1", "s2", "s3"), ("s4", "s5", "s6")],
        }
    )
    artifact = build_report_pdf(
        title="Plant Expression Explorer report",
        generated_at=_GENERATED_AT,
        intro_lines=(
            "This report mirrors descriptive result tables already shown "
            "in the application.",
        ),
        sections=(
            ReportSection(title="Dataset context", table=pd.DataFrame({"Field": ["Title"], "Value": ["Demo"]})),
            ReportSection(title="Sample quality control", table=qc_table),
            ReportSection(
                title="Differential expression",
                note="Not included: no differential-expression results were supplied.",
            ),
        ),
        filename="pee-report.pdf",
    )

    assert isinstance(artifact, ReportArtifact)
    assert artifact.filename == "pee-report.pdf"
    assert artifact.media_type == PDF_MEDIA_TYPE
    assert artifact.section_titles == (
        "Dataset context",
        "Sample quality control",
        "Differential expression",
    )
    assert artifact.data.startswith(b"%PDF-")
    assert artifact.data.rstrip().endswith(b"%%EOF")

    text = _extract_text(artifact.data)
    assert "Plant Expression Explorer report" in text
    assert "Dataset context" in text
    assert "Sample quality control" in text
    assert "control" in text
    assert "s1, s2, s3" in text
    assert "Differential expression" in text
    assert "Not included: no differential-expression results were supplied." in text


def test_build_report_pdf_renders_missing_values_as_empty_and_preserves_precision() -> None:
    table = pd.DataFrame({"gene_id": ["g1"], "value": [1.0 / 3.0], "flag": [pd.NA]})
    artifact = build_report_pdf(
        title="Report",
        generated_at=_GENERATED_AT,
        intro_lines=(),
        sections=(ReportSection(title="Values", table=table),),
    )
    text = _extract_text(artifact.data)
    assert str(1.0 / 3.0) in text


def test_build_report_pdf_does_not_mutate_the_supplied_table() -> None:
    table = pd.DataFrame({"a": [1, 2], "b": ["x", "y"]})
    original = table.copy(deep=True)
    build_report_pdf(
        title="Report",
        generated_at=_GENERATED_AT,
        intro_lines=(),
        sections=(ReportSection(title="Section", table=table),),
    )
    pd.testing.assert_frame_equal(table, original)


def test_rejects_a_blank_title() -> None:
    error = _capture_error(title="   ")
    assert error.reason is ReportExportErrorReason.INVALID_TITLE


def test_rejects_a_non_string_title() -> None:
    error = _capture_error(title=123)
    assert error.reason is ReportExportErrorReason.INVALID_TITLE


def test_rejects_an_unsafe_filename() -> None:
    error = _capture_error(filename="../report.pdf")
    assert error.reason is ReportExportErrorReason.INVALID_FILENAME


def test_rejects_a_filename_without_pdf_extension() -> None:
    error = _capture_error(filename="report.csv")
    assert error.reason is ReportExportErrorReason.INVALID_FILENAME


def test_rejects_an_empty_sections_tuple() -> None:
    error = _capture_error(sections=())
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_a_non_reportsection_item() -> None:
    error = _capture_error(sections=({"title": "x"},))
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_a_section_with_a_blank_title() -> None:
    error = _capture_error(
        sections=(ReportSection(title="  ", table=pd.DataFrame({"a": [1]})),)
    )
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_duplicate_section_titles() -> None:
    error = _capture_error(
        sections=(
            ReportSection(title="Dup", table=pd.DataFrame({"a": [1]})),
            ReportSection(title="Dup", table=pd.DataFrame({"a": [2]})),
        )
    )
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_a_section_with_both_table_and_note() -> None:
    error = _capture_error(
        sections=(
            ReportSection(
                title="Both", table=pd.DataFrame({"a": [1]}), note="also a note"
            ),
        )
    )
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_a_section_with_neither_table_nor_note() -> None:
    error = _capture_error(sections=(ReportSection(title="Neither"),))
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_a_section_table_that_is_not_a_dataframe() -> None:
    error = _capture_error(
        sections=(ReportSection(title="Bad table", table="not a dataframe"),)
    )
    assert error.reason is ReportExportErrorReason.INVALID_SECTIONS


def test_rejects_an_unencodable_character_in_a_table_cell() -> None:
    error = _capture_error(
        sections=(
            ReportSection(
                title="Unicode", table=pd.DataFrame({"gene_id": ["基因1"]})
            ),
        )
    )
    assert error.reason is ReportExportErrorReason.UNSUPPORTED_CHARACTER


def test_rejects_an_unencodable_character_in_the_title() -> None:
    error = _capture_error(title="基因报告")
    assert error.reason is ReportExportErrorReason.UNSUPPORTED_CHARACTER


def test_rejects_an_unencodable_character_in_an_intro_line() -> None:
    error = _capture_error(intro_lines=("基因报告说明",))
    assert error.reason is ReportExportErrorReason.UNSUPPORTED_CHARACTER


def test_rejects_an_unencodable_character_in_a_note() -> None:
    error = _capture_error(
        sections=(ReportSection(title="Note section", note="基因未提供"),)
    )
    assert error.reason is ReportExportErrorReason.UNSUPPORTED_CHARACTER


def test_report_section_and_artifact_are_frozen() -> None:
    from dataclasses import FrozenInstanceError

    section = ReportSection(title="x", table=pd.DataFrame({"a": [1]}))
    with pytest.raises(FrozenInstanceError):
        section.title = "y"  # type: ignore[misc]

    artifact = build_report_pdf(
        title="Report",
        generated_at=_GENERATED_AT,
        intro_lines=(),
        sections=(ReportSection(title="Section", table=pd.DataFrame({"a": [1]})),),
    )
    with pytest.raises(FrozenInstanceError):
        artifact.filename = "other.pdf"  # type: ignore[misc]
