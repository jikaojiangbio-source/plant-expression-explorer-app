"""CSV loading and validation for the application's three input tables."""

from __future__ import annotations

from math import isfinite
from pathlib import Path
from typing import IO, Any

import pandas as pd


class DataValidationError(ValueError):
    """Raised when an input table does not satisfy the documented contract."""


def read_csv(source: str | Path | IO[Any]) -> pd.DataFrame:
    """Read a CSV source and turn parser failures into a user-facing error."""

    if hasattr(source, "seek"):
        source.seek(0)

    try:
        return pd.read_csv(source)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError) as exc:
        raise DataValidationError(f"Could not read this file as CSV: {exc}") from exc


def validate_expression_matrix(expression: pd.DataFrame) -> None:
    """Validate a gene-by-sample normalized expression matrix."""

    _validate_table(expression, "Expression matrix")
    _require_columns(expression, ("gene_id",), "Expression matrix")
    _validate_identifier(expression["gene_id"], "gene_id", "Expression matrix")

    sample_columns = [column for column in expression.columns if column != "gene_id"]
    if not sample_columns:
        raise DataValidationError(
            "Expression matrix must contain at least one numeric sample column."
        )
    if any(not str(column).strip() for column in sample_columns):
        raise DataValidationError("Expression matrix contains an empty sample name.")

    for column in sample_columns:
        numeric = _as_numeric(expression[column], str(column), allow_missing=False)
        if not numeric.map(isfinite).all():
            raise DataValidationError(
                f"Expression sample column '{column}' contains an infinite value."
            )


def validate_metadata(metadata: pd.DataFrame) -> None:
    """Validate sample metadata containing sample identifiers and conditions."""

    _validate_table(metadata, "Sample metadata")
    _require_columns(metadata, ("sample_id", "condition"), "Sample metadata")
    _validate_identifier(metadata["sample_id"], "sample_id", "Sample metadata")

    condition = metadata["condition"]
    if condition.isna().any() or condition.astype("string").str.strip().eq("").any():
        raise DataValidationError(
            "Sample metadata column 'condition' contains a missing or empty value."
        )


def validate_de_results(de_results: pd.DataFrame) -> None:
    """Validate precomputed gene-level differential-expression results."""

    _validate_table(de_results, "Differential-expression results")
    required = ("gene_id", "log2FoldChange", "pvalue", "padj")
    _require_columns(de_results, required, "Differential-expression results")
    _validate_identifier(
        de_results["gene_id"], "gene_id", "Differential-expression results"
    )

    for column in ("log2FoldChange", "pvalue", "padj"):
        numeric = _as_numeric(de_results[column], column, allow_missing=True)
        finite_values = numeric.dropna()
        if not finite_values.map(isfinite).all():
            raise DataValidationError(
                f"Differential-expression column '{column}' contains an infinite value."
            )

    for column in ("pvalue", "padj"):
        numeric = pd.to_numeric(de_results[column], errors="coerce").dropna()
        if not numeric.between(0, 1, inclusive="both").all():
            raise DataValidationError(
                f"Differential-expression column '{column}' must be between 0 and 1."
            )


def align_expression_and_metadata(
    expression: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """Validate both tables and return metadata in expression sample order."""

    validate_expression_matrix(expression)
    validate_metadata(metadata)

    expression_samples = [str(column) for column in expression.columns if column != "gene_id"]
    metadata_samples = metadata["sample_id"].astype("string").tolist()

    missing_metadata = sorted(set(expression_samples) - set(metadata_samples))
    unexpected_metadata = sorted(set(metadata_samples) - set(expression_samples))
    if missing_metadata or unexpected_metadata:
        details: list[str] = []
        if missing_metadata:
            details.append("missing metadata for: " + ", ".join(missing_metadata))
        if unexpected_metadata:
            details.append("metadata without expression columns: " + ", ".join(unexpected_metadata))
        raise DataValidationError("Sample identifiers do not match; " + "; ".join(details) + ".")

    aligned = metadata.copy()
    aligned["sample_id"] = aligned["sample_id"].astype("string")
    return aligned.set_index("sample_id").loc[expression_samples].reset_index()


def _validate_table(table: pd.DataFrame, table_name: str) -> None:
    if not isinstance(table, pd.DataFrame):
        raise DataValidationError(f"{table_name} must be a table.")
    if table.empty:
        raise DataValidationError(f"{table_name} must contain at least one row.")
    duplicate_columns = table.columns[table.columns.duplicated()].tolist()
    if duplicate_columns:
        raise DataValidationError(
            f"{table_name} contains duplicate columns: {', '.join(map(str, duplicate_columns))}."
        )


def _require_columns(table: pd.DataFrame, required: tuple[str, ...], table_name: str) -> None:
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise DataValidationError(
            f"{table_name} is missing required column(s): {', '.join(missing)}."
        )


def _validate_identifier(series: pd.Series, column: str, table_name: str) -> None:
    normalized = series.astype("string").str.strip()
    if normalized.isna().any() or normalized.eq("").any():
        raise DataValidationError(
            f"{table_name} column '{column}' contains a missing or empty identifier."
        )
    if normalized.duplicated().any():
        duplicates = normalized[normalized.duplicated(keep=False)].unique().tolist()
        raise DataValidationError(
            f"{table_name} column '{column}' contains duplicate identifiers: "
            + ", ".join(map(str, duplicates))
            + "."
        )


def _as_numeric(series: pd.Series, column: str, *, allow_missing: bool) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    blank = series.astype("string").str.strip().eq("")
    invalid = numeric.isna() & series.notna() & ~blank
    if invalid.any():
        raise DataValidationError(f"Column '{column}' must contain numeric values.")
    if not allow_missing and numeric.isna().any():
        raise DataValidationError(f"Column '{column}' contains a missing value.")
    return numeric
