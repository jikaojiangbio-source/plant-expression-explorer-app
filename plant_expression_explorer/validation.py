"""Pure, non-mutating validation for the biological input tables."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum
from math import isfinite
from typing import Iterable

import pandas as pd
from pandas.api.types import is_bool, is_complex


class Severity(StrEnum):
    """Severity of a validation issue."""

    ERROR = "error"
    WARNING = "warning"
    INFORMATION = "information"


class IssueCode(StrEnum):
    """Stable machine-readable identifiers for loading and validation issues."""

    CSV_READ_ERROR = "CSV_READ_ERROR"
    POSSIBLE_DELIMITER_MISMATCH = "POSSIBLE_DELIMITER_MISMATCH"
    DEMO_FILE_MISSING = "DEMO_FILE_MISSING"
    DE_RESULTS_NOT_SUPPLIED = "DE_RESULTS_NOT_SUPPLIED"
    NOT_A_TABLE = "NOT_A_TABLE"
    EMPTY_TABLE = "EMPTY_TABLE"
    MISSING_REQUIRED_COLUMN = "MISSING_REQUIRED_COLUMN"
    EMPTY_COLUMN_NAME = "EMPTY_COLUMN_NAME"
    DUPLICATE_COLUMN_NAME = "DUPLICATE_COLUMN_NAME"
    NO_SAMPLE_COLUMNS = "NO_SAMPLE_COLUMNS"
    LOW_SAMPLE_COUNT = "LOW_SAMPLE_COUNT"
    MISSING_IDENTIFIER = "MISSING_IDENTIFIER"
    DUPLICATE_IDENTIFIER = "DUPLICATE_IDENTIFIER"
    WHITESPACE_IN_IDENTIFIER = "WHITESPACE_IN_IDENTIFIER"
    WHITESPACE_IN_SAMPLE_NAME = "WHITESPACE_IN_SAMPLE_NAME"
    WHITESPACE_IN_CATEGORY_VALUE = "WHITESPACE_IN_CATEGORY_VALUE"
    MISSING_REQUIRED_VALUE = "MISSING_REQUIRED_VALUE"
    SINGLE_CONDITION = "SINGLE_CONDITION"
    CONDITION_WITH_SINGLE_SAMPLE = "CONDITION_WITH_SINGLE_SAMPLE"
    NON_NUMERIC_VALUE = "NON_NUMERIC_VALUE"
    COERCIBLE_NUMERIC_STRING = "COERCIBLE_NUMERIC_STRING"
    MISSING_EXPRESSION_VALUE = "MISSING_EXPRESSION_VALUE"
    NON_FINITE_VALUE = "NON_FINITE_VALUE"
    NEGATIVE_EXPRESSION_VALUE = "NEGATIVE_EXPRESSION_VALUE"
    MISSING_DE_STATISTIC = "MISSING_DE_STATISTIC"
    NO_USABLE_VALUES = "NO_USABLE_VALUES"
    VALUE_OUT_OF_RANGE = "VALUE_OUT_OF_RANGE"
    ZERO_ADJUSTED_P_VALUE = "ZERO_ADJUSTED_P_VALUE"
    MISSING_METADATA_SAMPLE = "MISSING_METADATA_SAMPLE"
    UNEXPECTED_METADATA_SAMPLE = "UNEXPECTED_METADATA_SAMPLE"
    NO_SAMPLE_OVERLAP = "NO_SAMPLE_OVERLAP"
    SAMPLE_ORDER_DIFFERS = "SAMPLE_ORDER_DIFFERS"
    DE_GENE_NOT_IN_EXPRESSION = "DE_GENE_NOT_IN_EXPRESSION"
    EXPRESSION_GENE_NOT_IN_DE = "EXPRESSION_GENE_NOT_IN_DE"
    NO_GENE_OVERLAP = "NO_GENE_OVERLAP"


_UNNAMED_COLUMN_PATTERN = re.compile(r"Unnamed: \d+")

_LIKELY_MISSING_VALUE_MARKERS = frozenset(
    {"-", "--", ".", "na", "n.a.", "?", "nd", "n.d.", "nr", "n.r.", "tbd"}
)


@dataclass(frozen=True)
class ValidationIssue:
    """One actionable problem or informational observation."""

    code: IssueCode
    severity: Severity
    table: str
    message: str
    column: str | None = None
    count: int = 1
    row_positions: tuple[int, ...] = ()
    example_values: tuple[str, ...] = ()
    related_table: str | None = None


@dataclass(frozen=True)
class ValidationReport:
    """Complete validation outcome, including non-blocking observations."""

    issues: tuple[ValidationIssue, ...] = ()

    @property
    def errors(self) -> tuple[ValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.severity is Severity.ERROR)

    @property
    def warnings(self) -> tuple[ValidationIssue, ...]:
        return tuple(issue for issue in self.issues if issue.severity is Severity.WARNING)

    @property
    def information(self) -> tuple[ValidationIssue, ...]:
        return tuple(
            issue for issue in self.issues if issue.severity is Severity.INFORMATION
        )

    @property
    def has_errors(self) -> bool:
        return bool(self.errors)


def validate_expression_matrix(expression: pd.DataFrame) -> ValidationReport:
    """Validate one gene-by-sample expression matrix without changing it."""

    table_name = "Expression matrix"
    issues = _table_issues(expression, table_name)
    if not isinstance(expression, pd.DataFrame):
        return ValidationReport(tuple(issues))

    if "gene_id" not in expression.columns:
        issues.append(
            _issue(
                IssueCode.MISSING_REQUIRED_COLUMN,
                Severity.ERROR,
                table_name,
                "Expression matrix is missing required column 'gene_id'.",
                column="gene_id",
            )
        )
    elif _column_is_unique(expression, "gene_id"):
        issues.extend(_identifier_issues(expression["gene_id"], "gene_id", table_name))

    sample_columns = [column for column in expression.columns if column != "gene_id"]
    if not sample_columns:
        issues.append(
            _issue(
                IssueCode.NO_SAMPLE_COLUMNS,
                Severity.ERROR,
                table_name,
                "Expression matrix must contain at least one sample column in addition to 'gene_id'.",
            )
        )
    elif len(sample_columns) == 1:
        issues.append(
            _issue(
                IssueCode.LOW_SAMPLE_COUNT,
                Severity.ERROR,
                table_name,
                "Expression matrix contains only one sample column; this application's input contract requires at least two sample columns before downstream exploration. This does not assess biological validity.",
            )
        )

    whitespace_samples = [
        str(column)
        for column in sample_columns
        if str(column) != str(column).strip()
    ]
    if whitespace_samples:
        issues.append(
            _issue(
                IssueCode.WHITESPACE_IN_SAMPLE_NAME,
                Severity.ERROR,
                table_name,
                "Expression matrix contains sample name(s) with leading or trailing whitespace: "
                + _format_examples(whitespace_samples)
                + ". Correct them explicitly; validation did not trim them.",
                count=len(whitespace_samples),
                examples=whitespace_samples,
            )
        )

    for column in sample_columns:
        if not _column_is_unique(expression, column):
            continue
        issues.extend(_expression_value_issues(expression[column], str(column)))

    return ValidationReport(tuple(issues))


def validate_sample_metadata(metadata: pd.DataFrame) -> ValidationReport:
    """Validate sample identifiers and conditions without changing metadata."""

    table_name = "Sample metadata"
    issues = _table_issues(metadata, table_name)
    if not isinstance(metadata, pd.DataFrame):
        return ValidationReport(tuple(issues))

    for column in ("sample_id", "condition"):
        if column not in metadata.columns:
            issues.append(
                _issue(
                    IssueCode.MISSING_REQUIRED_COLUMN,
                    Severity.ERROR,
                    table_name,
                    f"Sample metadata is missing required column '{column}'.",
                    column=column,
                )
            )

    if "sample_id" in metadata.columns and _column_is_unique(metadata, "sample_id"):
        issues.extend(_identifier_issues(metadata["sample_id"], "sample_id", table_name))

    if "condition" in metadata.columns and _column_is_unique(metadata, "condition"):
        condition = metadata["condition"]
        missing = _missing_mask(condition)
        if missing.any():
            issues.append(
                _issue_from_mask(
                    IssueCode.MISSING_REQUIRED_VALUE,
                    Severity.ERROR,
                    table_name,
                    "Sample metadata column 'condition' contains {count} missing or empty value(s).",
                    condition,
                    missing,
                    column="condition",
                )
            )

        whitespace = _whitespace_mask(condition) & ~missing
        if whitespace.any():
            issues.append(
                _issue_from_mask(
                    IssueCode.WHITESPACE_IN_CATEGORY_VALUE,
                    Severity.ERROR,
                    table_name,
                    "Sample metadata column 'condition' contains value(s) with leading or trailing whitespace: {examples}. Correct them explicitly; validation did not trim them.",
                    condition,
                    whitespace,
                    column="condition",
                )
            )

        usable_conditions = condition[~missing].astype("string")
        counts = usable_conditions.value_counts(dropna=False)
        if len(counts) == 1:
            issues.append(
                _issue(
                    IssueCode.SINGLE_CONDITION,
                    Severity.WARNING,
                    table_name,
                    "Sample metadata contains only one condition; descriptive analyses can continue, but between-condition views will be unavailable.",
                    column="condition",
                )
            )

        singleton_conditions = sorted(str(value) for value in counts[counts == 1].index)
        if singleton_conditions:
            issues.append(
                _issue(
                    IssueCode.CONDITION_WITH_SINGLE_SAMPLE,
                    Severity.WARNING,
                    table_name,
                    "Condition(s) "
                    + _format_examples(singleton_conditions)
                    + " contain only one sample; descriptive display can continue, but replication must not be inferred.",
                    column="condition",
                    count=len(singleton_conditions),
                    examples=singleton_conditions,
                )
            )

    return ValidationReport(tuple(issues))


def validate_de_results(de_results: pd.DataFrame) -> ValidationReport:
    """Validate precomputed differential-expression results without changing them."""

    table_name = "Differential-expression results"
    issues = _table_issues(de_results, table_name)
    if not isinstance(de_results, pd.DataFrame):
        return ValidationReport(tuple(issues))

    required = ("gene_id", "log2FoldChange", "pvalue", "padj")
    for column in required:
        if column not in de_results.columns:
            suffix = (
                " Adjusted p-values will not be calculated or inferred by this application."
                if column == "padj"
                else ""
            )
            issues.append(
                _issue(
                    IssueCode.MISSING_REQUIRED_COLUMN,
                    Severity.ERROR,
                    table_name,
                    f"Differential-expression results are missing required column '{column}'.{suffix}",
                    column=column,
                )
            )

    if "gene_id" in de_results.columns and _column_is_unique(de_results, "gene_id"):
        issues.extend(_identifier_issues(de_results["gene_id"], "gene_id", table_name))

    for column in ("log2FoldChange", "pvalue", "padj"):
        if column not in de_results.columns or not _column_is_unique(de_results, column):
            continue
        issues.extend(_de_value_issues(de_results[column], column))

    return ValidationReport(tuple(issues))


def _table_issues(table: object, table_name: str) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    if not isinstance(table, pd.DataFrame):
        return [
            _issue(
                IssueCode.NOT_A_TABLE,
                Severity.ERROR,
                table_name,
                f"{table_name} must be provided as a table.",
            )
        ]

    if table.empty:
        issues.append(
            _issue(
                IssueCode.EMPTY_TABLE,
                Severity.ERROR,
                table_name,
                f"{table_name} must contain at least one data row.",
            )
        )

    empty_columns = [
        str(column)
        for column in table.columns
        if not str(column).strip() or _UNNAMED_COLUMN_PATTERN.fullmatch(str(column))
    ]
    if empty_columns:
        issues.append(
            _issue(
                IssueCode.EMPTY_COLUMN_NAME,
                Severity.ERROR,
                table_name,
                f"{table_name} contains {len(empty_columns)} column(s) with a blank "
                "header: "
                + _format_examples(empty_columns)
                + ". A name shown as 'Unnamed: N' means the original header cell "
                "was blank, which commonly happens when a spreadsheet export "
                "leaves a trailing empty column; remove that column and "
                "re-upload.",
                count=len(empty_columns),
                examples=empty_columns,
            )
        )

    duplicate_columns = table.columns[table.columns.duplicated(keep=False)].tolist()
    if duplicate_columns:
        unique_duplicates = list(dict.fromkeys(str(value) for value in duplicate_columns))
        issues.append(
            _issue(
                IssueCode.DUPLICATE_COLUMN_NAME,
                Severity.ERROR,
                table_name,
                f"{table_name} contains duplicate column name(s): "
                + _format_examples(unique_duplicates)
                + ".",
                count=len(unique_duplicates),
                examples=unique_duplicates,
            )
        )
    return issues


def _identifier_issues(
    series: pd.Series, column: str, table_name: str
) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    missing = _missing_mask(series)
    if missing.any():
        issues.append(
            _issue_from_mask(
                IssueCode.MISSING_IDENTIFIER,
                Severity.ERROR,
                table_name,
                f"{table_name} column '{column}' contains {{count}} missing or empty identifier(s).",
                series,
                missing,
                column=column,
            )
        )

    whitespace = _whitespace_mask(series) & ~missing
    if whitespace.any():
        issues.append(
            _issue_from_mask(
                IssueCode.WHITESPACE_IN_IDENTIFIER,
                Severity.ERROR,
                table_name,
                f"{table_name} column '{column}' contains identifier(s) with leading or trailing whitespace: {{examples}}. Correct them explicitly; validation did not trim them.",
                series,
                whitespace,
                column=column,
            )
        )

    comparable = series.astype("string")
    duplicate = comparable.duplicated(keep=False) & ~missing
    if duplicate.any():
        issues.append(
            _issue_from_mask(
                IssueCode.DUPLICATE_IDENTIFIER,
                Severity.ERROR,
                table_name,
                f"{table_name} column '{column}' contains duplicate identifier(s): {{examples}}. Resolve them explicitly; no rows were removed or aggregated.",
                series,
                duplicate,
                column=column,
            )
        )
    return issues


def _missing_value_marker_hint(series: pd.Series, invalid: pd.Series) -> str:
    invalid_values = series[invalid]
    if invalid_values.map(_looks_like_missing_marker).any():
        return (
            " If any of these represent a missing measurement, leave the cell "
            "blank rather than using a placeholder such as '-' or '.'; this "
            "application does not treat such placeholders as missing values."
        )
    return ""


def _looks_like_missing_marker(value: object) -> bool:
    return isinstance(value, str) and value.strip().casefold() in _LIKELY_MISSING_VALUE_MARKERS


def _expression_value_issues(series: pd.Series, column: str) -> list[ValidationIssue]:
    table_name = "Expression matrix"
    numeric, missing, invalid, coercible, non_finite = _numeric_masks(series)
    issues: list[ValidationIssue] = []

    if invalid.any():
        issues.append(
            _issue_from_mask(
                IssueCode.NON_NUMERIC_VALUE,
                Severity.ERROR,
                table_name,
                f"Expression column '{column}' contains {{count}} non-numeric value(s): {{examples}}."
                + _missing_value_marker_hint(series, invalid),
                series,
                invalid,
                column=column,
            )
        )
    if coercible.any():
        issues.append(
            _issue_from_mask(
                IssueCode.COERCIBLE_NUMERIC_STRING,
                Severity.INFORMATION,
                table_name,
                f"Expression column '{column}' contains {{count}} numeric value(s) stored as text; they are safely coercible and the uploaded table has not been modified.",
                series,
                coercible,
                column=column,
            )
        )
    if missing.any():
        issues.append(
            _issue_from_mask(
                IssueCode.MISSING_EXPRESSION_VALUE,
                Severity.ERROR,
                table_name,
                f"Expression column '{column}' contains {{count}} missing value(s). Provide complete values; validation did not impute or remove data.",
                series,
                missing,
                column=column,
            )
        )
    if non_finite.any():
        issues.append(
            _issue_from_mask(
                IssueCode.NON_FINITE_VALUE,
                Severity.ERROR,
                table_name,
                f"Expression column '{column}' contains {{count}} infinite value(s).",
                series,
                non_finite,
                column=column,
            )
        )

    usable = ~(missing | invalid | non_finite)
    negative = usable & numeric.lt(0)
    if negative.any():
        issues.append(
            _issue_from_mask(
                IssueCode.NEGATIVE_EXPRESSION_VALUE,
                Severity.WARNING,
                table_name,
                "Expression matrix contains {count} negative value(s). Negative values may be valid for transformed or centred expression data and have been retained; confirm that the matrix data type is appropriate.",
                series,
                negative,
                column=column,
            )
        )
    return issues


def _de_value_issues(series: pd.Series, column: str) -> list[ValidationIssue]:
    table_name = "Differential-expression results"
    numeric, missing, invalid, coercible, non_finite = _numeric_masks(series)
    issues: list[ValidationIssue] = []

    if invalid.any():
        issues.append(
            _issue_from_mask(
                IssueCode.NON_NUMERIC_VALUE,
                Severity.ERROR,
                table_name,
                f"Differential-expression column '{column}' contains {{count}} non-numeric value(s): {{examples}}."
                + _missing_value_marker_hint(series, invalid),
                series,
                invalid,
                column=column,
            )
        )
    if coercible.any():
        issues.append(
            _issue_from_mask(
                IssueCode.COERCIBLE_NUMERIC_STRING,
                Severity.INFORMATION,
                table_name,
                f"Differential-expression column '{column}' contains {{count}} numeric value(s) stored as text; they are safely coercible and the uploaded table has not been modified.",
                series,
                coercible,
                column=column,
            )
        )

    if len(series) and missing.all():
        issues.append(
            _issue(
                IssueCode.NO_USABLE_VALUES,
                Severity.ERROR,
                table_name,
                f"Differential-expression column '{column}' contains no usable values.",
                column=column,
                count=int(missing.sum()),
            )
        )
    elif missing.any():
        suffix = (
            " They have been retained and must not be treated as zero."
            if column == "padj"
            else " They have been retained and will not be imputed."
        )
        issues.append(
            _issue_from_mask(
                IssueCode.MISSING_DE_STATISTIC,
                Severity.WARNING,
                table_name,
                f"Differential-expression column '{column}' contains {{count}} missing value(s).{suffix}",
                series,
                missing,
                column=column,
            )
        )

    if non_finite.any():
        issues.append(
            _issue_from_mask(
                IssueCode.NON_FINITE_VALUE,
                Severity.ERROR,
                table_name,
                f"Differential-expression column '{column}' contains {{count}} infinite value(s).",
                series,
                non_finite,
                column=column,
            )
        )

    usable = ~(missing | invalid | non_finite)
    if column in ("pvalue", "padj"):
        outside = usable & ~numeric.between(0, 1, inclusive="both")
        if outside.any():
            issues.append(
                _issue_from_mask(
                    IssueCode.VALUE_OUT_OF_RANGE,
                    Severity.ERROR,
                    table_name,
                    f"Differential-expression column '{column}' contains {{count}} value(s) outside the inclusive range 0 to 1: {{examples}}.",
                    series,
                    outside,
                    column=column,
                )
            )

    if column == "padj":
        zero = usable & numeric.eq(0)
        if zero.any():
            issues.append(
                _issue_from_mask(
                    IssueCode.ZERO_ADJUSTED_P_VALUE,
                    Severity.INFORMATION,
                    table_name,
                    "Differential-expression results contain {count} adjusted p-value(s) equal to 0. Zero is a valid boundary value and has been retained; no significance classification was applied.",
                    series,
                    zero,
                    column=column,
                )
            )
    return issues


def _numeric_masks(
    series: pd.Series,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series, pd.Series]:
    missing = _missing_mask(series)
    complex_value = series.map(is_complex)
    numeric_source = series.astype("object").mask(complex_value, other=None)
    numeric = pd.to_numeric(numeric_source, errors="coerce")
    boolean = series.map(is_bool)
    invalid = (numeric.isna() & ~missing) | boolean | complex_value
    finite = numeric.map(lambda value: isfinite(value) if pd.notna(value) else True)
    non_finite = numeric.notna() & ~invalid & ~finite
    is_text = series.map(lambda value: isinstance(value, str))
    coercible = is_text & ~missing & ~invalid & ~non_finite
    return numeric, missing, invalid, coercible, non_finite


def _missing_mask(series: pd.Series) -> pd.Series:
    blank = series.astype("string").str.strip().eq("").fillna(False)
    return series.isna() | blank


def _whitespace_mask(series: pd.Series) -> pd.Series:
    text = series.astype("string")
    return (text != text.str.strip()).fillna(False)


def _column_is_unique(table: pd.DataFrame, column: object) -> bool:
    return sum(existing == column for existing in table.columns) == 1


def _issue(
    code: IssueCode,
    severity: Severity,
    table: str,
    message: str,
    *,
    column: str | None = None,
    count: int = 1,
    rows: Iterable[int] = (),
    examples: Iterable[object] = (),
    related_table: str | None = None,
) -> ValidationIssue:
    return ValidationIssue(
        code=code,
        severity=severity,
        table=table,
        message=message,
        column=column,
        count=count,
        row_positions=tuple(rows)[:5],
        example_values=tuple(str(value) for value in examples)[:5],
        related_table=related_table,
    )


def _issue_from_mask(
    code: IssueCode,
    severity: Severity,
    table: str,
    message_template: str,
    series: pd.Series,
    mask: pd.Series,
    *,
    column: str,
) -> ValidationIssue:
    count = int(mask.sum())
    values = [str(value) for value in series[mask].tolist()]
    examples = tuple(dict.fromkeys(values))[:5]
    positions = tuple(position for position, value in enumerate(mask.tolist(), 1) if value)[:5]
    message = message_template.format(
        count=count,
        examples=_format_examples(examples),
    )
    return _issue(
        code,
        severity,
        table,
        message,
        column=column,
        count=count,
        rows=positions,
        examples=examples,
    )


def _format_examples(values: Iterable[object]) -> str:
    return ", ".join(str(value) for value in list(values)[:5])
