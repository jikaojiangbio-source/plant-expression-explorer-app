"""Deterministic, non-mutating PDF report assembly for descriptive result tables.

This module performs no scientific calculation of its own. Callers supply
already-computed descriptive tables (the same ones already shown and
exported elsewhere in the application); this module only lays them out as a
single downloadable PDF. It never fits a model, computes a statistic, or
reorders/renames anything it is given.
"""

from __future__ import annotations

import io
import re
from dataclasses import dataclass
from datetime import datetime
from enum import StrEnum
from xml.sax.saxutils import escape

import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


class ReportExportErrorReason(StrEnum):
    """Stable reasons for controlled report-assembly failures."""

    INVALID_TITLE = "INVALID_TITLE"
    INVALID_SECTIONS = "INVALID_SECTIONS"
    INVALID_FILENAME = "INVALID_FILENAME"
    UNSUPPORTED_CHARACTER = "UNSUPPORTED_CHARACTER"


class ReportExportError(ValueError):
    """Expected failure when a descriptive PDF report cannot be assembled safely."""

    def __init__(self, reason: ReportExportErrorReason, message: str) -> None:
        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class ReportSection:
    """One titled block of the generated report: a descriptive table, or a note.

    Exactly one of ``table`` and ``note`` must be supplied. ``table`` is
    rendered as a grid, every cell rendered by ``str(value)`` with no
    rounding; ``note`` is rendered as a short paragraph, used when the caller
    intentionally omits a section it could not compute (for example, no
    differential-expression results were supplied) and wants to disclose why.
    """

    title: str
    table: pd.DataFrame | None = None
    note: str | None = None


@dataclass(frozen=True)
class ReportArtifact:
    """One in-memory PDF report and its stable descriptive metadata."""

    filename: str
    media_type: str
    data: bytes
    section_titles: tuple[str, ...]


PDF_MEDIA_TYPE = "application/pdf"
_SAFE_FILENAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.pdf")
_TEXT_ENCODING = "cp1252"


def build_report_pdf(
    *,
    title: str,
    generated_at: datetime,
    intro_lines: tuple[str, ...],
    sections: tuple[ReportSection, ...],
    filename: str = "plant-expression-explorer-report.pdf",
) -> ReportArtifact:
    """Assemble a descriptive PDF report from already-computed inputs.

    Row and column order in each supplied table are preserved exactly. A
    tuple/list cell (for example a sample-ID membership column) is rendered
    as a comma-separated list; a missing scalar value (as detected by
    ``pandas.isna``) is rendered as an empty cell, matching this
    application's CSV export convention. The built-in report font supports
    Windows-1252 (a superset of Latin-1); a title, intro line, section
    title, note, or table cell containing any other character raises a
    controlled error naming the offending section instead of silently
    dropping or substituting that character.
    """

    safe_title = _validated_title(title)
    safe_filename = _validated_filename(filename)
    safe_intro_lines = tuple(intro_lines)
    for index, line in enumerate(safe_intro_lines):
        _ensure_encodable(line, f"intro line {index + 1}")
    safe_sections = _validated_sections(sections)

    buffer = io.BytesIO()
    _render_pdf(buffer, safe_title, generated_at, safe_intro_lines, safe_sections)

    return ReportArtifact(
        filename=safe_filename,
        media_type=PDF_MEDIA_TYPE,
        data=buffer.getvalue(),
        section_titles=tuple(section.title for section in safe_sections),
    )


def _validated_title(title: object) -> str:
    if not isinstance(title, str) or not title.strip():
        raise ReportExportError(
            ReportExportErrorReason.INVALID_TITLE,
            "The report title must be a non-empty string.",
        )
    _ensure_encodable(title, "the report title")
    return title


def _validated_filename(filename: object) -> str:
    if not isinstance(filename, str) or _SAFE_FILENAME.fullmatch(filename) is None:
        raise ReportExportError(
            ReportExportErrorReason.INVALID_FILENAME,
            "The download filename must be a simple .pdf name containing only "
            "letters, numbers, periods, underscores, and hyphens.",
        )
    return filename


def _validated_sections(
    sections: object,
) -> tuple[ReportSection, ...]:
    if not isinstance(sections, (tuple, list)) or len(sections) == 0:
        raise ReportExportError(
            ReportExportErrorReason.INVALID_SECTIONS,
            "At least one report section is required.",
        )
    validated: list[ReportSection] = []
    seen_titles: set[str] = set()
    for section in sections:
        if not isinstance(section, ReportSection):
            raise ReportExportError(
                ReportExportErrorReason.INVALID_SECTIONS,
                "Every report section must be a ReportSection instance.",
            )
        if not isinstance(section.title, str) or not section.title.strip():
            raise ReportExportError(
                ReportExportErrorReason.INVALID_SECTIONS,
                "Every report section must have a non-empty title.",
            )
        if section.title in seen_titles:
            raise ReportExportError(
                ReportExportErrorReason.INVALID_SECTIONS,
                f"Report section title '{section.title}' is repeated.",
            )
        seen_titles.add(section.title)
        has_table = section.table is not None
        has_note = section.note is not None
        if has_table == has_note:
            raise ReportExportError(
                ReportExportErrorReason.INVALID_SECTIONS,
                f"Report section '{section.title}' must supply exactly one of "
                "'table' or 'note'.",
            )
        if has_table and not isinstance(section.table, pd.DataFrame):
            raise ReportExportError(
                ReportExportErrorReason.INVALID_SECTIONS,
                f"Report section '{section.title}' table must be a pandas "
                "DataFrame.",
            )
        _ensure_encodable(section.title, f"the title of section '{section.title}'")
        if has_note:
            _ensure_encodable(section.note, f"the note in section '{section.title}'")
        validated.append(section)
    return tuple(validated)


def _ensure_encodable(text: str, context: str) -> None:
    try:
        text.encode(_TEXT_ENCODING)
    except UnicodeEncodeError as error:
        raise ReportExportError(
            ReportExportErrorReason.UNSUPPORTED_CHARACTER,
            f"{context.capitalize()} contains a character "
            f"('{error.object[error.start:error.end]}') that the PDF report's "
            "supported font cannot render.",
        ) from error


def _cell_to_text(value: object, context: str) -> str:
    if isinstance(value, (tuple, list)):
        text = ", ".join(str(item) for item in value)
    else:
        try:
            is_missing = bool(pd.isna(value))
        except (TypeError, ValueError):
            is_missing = False
        text = "" if is_missing else str(value)
    _ensure_encodable(text, context)
    return text


def _render_pdf(
    buffer: io.BytesIO,
    title: str,
    generated_at: datetime,
    intro_lines: tuple[str, ...],
    sections: tuple[ReportSection, ...],
) -> None:
    styles = getSampleStyleSheet()
    body_style = styles["BodyText"]
    heading_style = styles["Heading2"]
    note_style = ParagraphStyle(
        "ReportNote",
        parent=body_style,
        textColor=colors.HexColor("#5b5b5b"),
        leading=14,
    )

    doc = SimpleDocTemplate(
        buffer,
        pagesize=LETTER,
        leftMargin=0.75 * inch,
        rightMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
        title=title,
    )

    flowables: list[object] = [
        Paragraph(escape(title), styles["Title"]),
        Paragraph(
            f"Generated {generated_at.isoformat(timespec='seconds')} from the "
            "currently active dataset in this session.",
            body_style,
        ),
        Spacer(1, 0.15 * inch),
    ]
    for line in intro_lines:
        flowables.append(Paragraph(escape(line), body_style))
    flowables.append(Spacer(1, 0.2 * inch))

    for section in sections:
        flowables.append(Paragraph(escape(section.title), heading_style))
        if section.table is not None:
            flowables.append(_build_table_flowable(section.table, section.title))
        else:
            flowables.append(Paragraph(escape(section.note), note_style))
        flowables.append(Spacer(1, 0.25 * inch))

    doc.build(flowables)


def _build_table_flowable(table: pd.DataFrame, section_title: str) -> Table:
    header = [
        _cell_to_text(column, f"a column header in section '{section_title}'")
        for column in table.columns
    ]
    rows = [
        [
            _cell_to_text(value, f"a table cell in section '{section_title}'")
            for value in row
        ]
        for row in table.itertuples(index=False, name=None)
    ]
    data = [header, *rows]
    grid = Table(data, repeatRows=1, hAlign="LEFT")
    grid.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#eef3ee")),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
                ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#c9c9c9")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f7f9f7")]),
            ]
        )
    )
    return grid
