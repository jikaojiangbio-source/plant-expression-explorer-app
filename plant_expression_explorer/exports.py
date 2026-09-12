"""Deterministic, non-mutating CSV exports for descriptive result tables."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from enum import StrEnum
from numbers import Integral, Real
from typing import Collection

import pandas as pd
from pandas.api.types import is_bool


class CsvExportErrorReason(StrEnum):
    """Stable reasons for controlled result-export failures."""

    INVALID_TABLE = "INVALID_TABLE"
    INVALID_FILENAME = "INVALID_FILENAME"
    INVALID_COLUMNS = "INVALID_COLUMNS"
    DUPLICATE_COLUMNS = "DUPLICATE_COLUMNS"
    INVALID_SEQUENCE_COLUMN = "INVALID_SEQUENCE_COLUMN"
    UNSUPPORTED_VALUE = "UNSUPPORTED_VALUE"
    INFINITE_VALUE = "INFINITE_VALUE"
    SERIALIZATION_ERROR = "SERIALIZATION_ERROR"


class CsvExportError(ValueError):
    """Expected failure when a result table cannot be exported safely."""

    def __init__(self, reason: CsvExportErrorReason, message: str) -> None:
        self.reason = reason
        super().__init__(message)


@dataclass(frozen=True)
class CsvExportArtifact:
    """One in-memory CSV download and its stable descriptive metadata."""

    filename: str
    media_type: str
    data: bytes
    row_count: int
    column_count: int
    column_names: tuple[str, ...]


CSV_MEDIA_TYPE = "text/csv; charset=utf-8"
_SAFE_FILENAME = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.csv")


def build_csv_export(
    table: pd.DataFrame,
    *,
    filename: str,
    json_sequence_columns: Collection[str] = (),
) -> CsvExportArtifact:
    """Serialize a descriptive result table without modifying the input.

    Row and column order are retained, the pandas index is intentionally not
    exported, and no display rounding is applied. Missing scalar values are
    represented by empty CSV fields. Columns named in ``json_sequence_columns``
    must contain tuples or lists of strings; each sequence is encoded as a
    compact JSON array in an independent working copy so membership remains
    explicit and reversible.
    """

    _require_dataframe(table)
    safe_filename = _validated_filename(filename)
    column_names = _validated_columns(table)
    sequence_columns = _validated_sequence_columns(
        json_sequence_columns,
        column_names,
    )
    export_table = table.copy(deep=True)
    _validate_and_encode_values(export_table, sequence_columns)

    try:
        csv_text = export_table.to_csv(
            index=False,
            na_rep="",
            lineterminator="\n",
        )
        data = csv_text.encode("utf-8")
    except (TypeError, UnicodeError, ValueError) as error:
        raise CsvExportError(
            CsvExportErrorReason.SERIALIZATION_ERROR,
            "The descriptive result table could not be serialized as UTF-8 CSV.",
        ) from error

    return CsvExportArtifact(
        filename=safe_filename,
        media_type=CSV_MEDIA_TYPE,
        data=data,
        row_count=len(table.index),
        column_count=len(column_names),
        column_names=column_names,
    )


def _require_dataframe(table: object) -> None:
    if not isinstance(table, pd.DataFrame):
        raise CsvExportError(
            CsvExportErrorReason.INVALID_TABLE,
            "The export source must be a pandas DataFrame.",
        )


def _validated_filename(filename: object) -> str:
    if not isinstance(filename, str) or _SAFE_FILENAME.fullmatch(filename) is None:
        raise CsvExportError(
            CsvExportErrorReason.INVALID_FILENAME,
            "The download filename must be a simple .csv name containing only "
            "letters, numbers, periods, underscores, and hyphens.",
        )
    return filename


def _validated_columns(table: pd.DataFrame) -> tuple[str, ...]:
    if len(table.columns) == 0:
        raise CsvExportError(
            CsvExportErrorReason.INVALID_COLUMNS,
            "The export table must contain at least one named column.",
        )
    if any(not isinstance(column, str) or not column for column in table.columns):
        raise CsvExportError(
            CsvExportErrorReason.INVALID_COLUMNS,
            "Every exported column must have a non-empty string name.",
        )
    column_names = tuple(table.columns)
    if len(set(column_names)) != len(column_names):
        raise CsvExportError(
            CsvExportErrorReason.DUPLICATE_COLUMNS,
            "Exported column names must be unique; no columns were renamed or removed.",
        )
    return column_names


def _validated_sequence_columns(
    columns: Collection[str],
    table_columns: tuple[str, ...],
) -> frozenset[str]:
    try:
        sequence_columns = tuple(columns)
    except TypeError as error:
        raise CsvExportError(
            CsvExportErrorReason.INVALID_SEQUENCE_COLUMN,
            "JSON sequence columns must be supplied as a collection of column names.",
        ) from error
    if any(not isinstance(column, str) for column in sequence_columns):
        raise CsvExportError(
            CsvExportErrorReason.INVALID_SEQUENCE_COLUMN,
            "Every JSON sequence column name must be a string.",
        )
    if len(set(sequence_columns)) != len(sequence_columns):
        raise CsvExportError(
            CsvExportErrorReason.INVALID_SEQUENCE_COLUMN,
            "JSON sequence column names must be unique.",
        )
    missing = [column for column in sequence_columns if column not in table_columns]
    if missing:
        raise CsvExportError(
            CsvExportErrorReason.INVALID_SEQUENCE_COLUMN,
            "JSON sequence column(s) are absent from the export table: "
            + ", ".join(missing)
            + ".",
        )
    return frozenset(sequence_columns)


def _validate_and_encode_values(
    table: pd.DataFrame,
    sequence_columns: frozenset[str],
) -> None:
    for column in table.columns:
        for row_position, value in enumerate(table[column].tolist()):
            if column in sequence_columns:
                table.iat[row_position, table.columns.get_loc(column)] = (
                    _encoded_string_sequence(value, column, row_position)
                )
                continue
            _validate_scalar(value, column, row_position)


def _encoded_string_sequence(
    value: object,
    column: str,
    row_position: int,
) -> str:
    if not isinstance(value, (tuple, list)) or any(
        not isinstance(item, str) for item in value
    ):
        raise CsvExportError(
            CsvExportErrorReason.UNSUPPORTED_VALUE,
            f"Column {column!r} contains a non-string sequence at zero-based "
            f"row position {row_position}.",
        )
    return json.dumps(list(value), ensure_ascii=False, separators=(",", ":"))


def _validate_scalar(value: object, column: str, row_position: int) -> None:
    if _is_missing(value):
        return
    if isinstance(value, str) or is_bool(value) or isinstance(value, Integral):
        return
    if isinstance(value, Real):
        if not math.isfinite(float(value)):
            raise CsvExportError(
                CsvExportErrorReason.INFINITE_VALUE,
                f"Column {column!r} contains an infinite value at zero-based "
                f"row position {row_position}; no partial export was produced.",
            )
        return
    raise CsvExportError(
        CsvExportErrorReason.UNSUPPORTED_VALUE,
        f"Column {column!r} contains unsupported value type "
        f"{type(value).__name__!r} at zero-based row position {row_position}.",
    )


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return isinstance(missing, bool) and missing
