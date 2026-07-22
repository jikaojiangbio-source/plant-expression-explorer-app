"""CSV reading and raw-header integrity checks."""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path
from typing import IO, Any

import pandas as pd

from plant_expression_explorer.validation import IssueCode


class CsvReadError(ValueError):
    """Raised when a source cannot be represented safely as a CSV table."""

    def __init__(self, message: str, *, code: IssueCode | None = None) -> None:
        super().__init__(message)
        self.code = code


def read_csv(source: str | Path | IO[Any]) -> pd.DataFrame:
    """Read CSV data after checking the raw header for duplicate names."""

    try:
        header = _read_raw_header(source)
        _validate_raw_header(header)
        if hasattr(source, "seek"):
            source.seek(0)
        return pd.read_csv(source)
    except CsvReadError:
        raise
    except (
        csv.Error,
        pd.errors.EmptyDataError,
        pd.errors.ParserError,
        UnicodeDecodeError,
    ) as exc:
        raise CsvReadError(f"Could not read this file as CSV: {exc}") from exc


def _read_raw_header(source: str | Path | IO[Any]) -> list[str]:
    if isinstance(source, (str, Path)):
        try:
            with Path(source).open("r", encoding="utf-8-sig", newline="") as handle:
                return next(csv.reader(handle))
        except StopIteration as exc:
            raise CsvReadError("Could not read this file as CSV: the file is empty.") from exc

    if not hasattr(source, "read"):
        raise CsvReadError("Could not read this file as CSV: unsupported source.")

    source.seek(0)
    content = source.read()
    source.seek(0)
    if isinstance(content, bytes):
        text = content.decode("utf-8-sig")
    else:
        text = str(content).lstrip("\ufeff")
    try:
        return next(csv.reader(StringIO(text)))
    except StopIteration as exc:
        raise CsvReadError("Could not read this file as CSV: the file is empty.") from exc


def _validate_raw_header(header: list[str]) -> None:
    duplicate_names = list(
        dict.fromkeys(name for name in header if header.count(name) > 1)
    )
    if duplicate_names:
        raise CsvReadError(
            "CSV contains duplicate column name(s): "
            + ", ".join(duplicate_names)
            + ". Rename them explicitly before continuing.",
            code=IssueCode.DUPLICATE_COLUMN_NAME,
        )
