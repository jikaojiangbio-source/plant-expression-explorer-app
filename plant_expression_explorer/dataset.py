"""Dataset loading, activation, and session-state helpers."""

from __future__ import annotations

from collections.abc import Mapping, MutableMapping
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, Literal, TypeAlias

import pandas as pd

from plant_expression_explorer.consistency import validate_input_tables
from plant_expression_explorer.data import CsvReadError, read_csv
from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
)

DatasetSource: TypeAlias = Literal["demo", "uploaded"]
CandidateStatus: TypeAlias = Literal["incomplete", "invalid", "valid"]
CsvSource: TypeAlias = str | Path | IO[Any]
DatasetTables: TypeAlias = tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]

CURRENT_DATASET_KEY = "pee_current_dataset"
CANDIDATE_REPORT_KEY = "pee_candidate_report"
CANDIDATE_SOURCE_KEY = "pee_candidate_source"
CANDIDATE_LABEL_KEY = "pee_candidate_label"
EXPRESSION_UPLOAD_KEY = "pee_upload_expression"
METADATA_UPLOAD_KEY = "pee_upload_metadata"
DE_RESULTS_UPLOAD_KEY = "pee_upload_de_results"

CANDIDATE_FEEDBACK_KEYS = (
    CANDIDATE_REPORT_KEY,
    CANDIDATE_SOURCE_KEY,
    CANDIDATE_LABEL_KEY,
)
UPLOADER_KEYS = (
    EXPRESSION_UPLOAD_KEY,
    METADATA_UPLOAD_KEY,
    DE_RESULTS_UPLOAD_KEY,
)
LEGACY_DATA_KEYS = ("expression", "metadata", "de_results")
APPLICATION_DATA_KEYS = (
    CURRENT_DATASET_KEY,
    *CANDIDATE_FEEDBACK_KEYS,
    *UPLOADER_KEYS,
    *LEGACY_DATA_KEYS,
)

DEMO_DIRECTORY = Path(__file__).resolve().parents[1] / "data" / "demo"
DEMO_SOURCE_LABEL = "Bundled synthetic demonstration data"

_TABLE_SPECS = (
    ("expression", "Expression matrix", "expression_matrix.csv"),
    ("metadata", "Sample metadata", "metadata.csv"),
    (
        "de_results",
        "Differential-expression results",
        "deg_results.csv",
    ),
)


@dataclass(frozen=True)
class DatasetBundle:
    """One complete validated dataset available to later pages.

    ``frozen=True`` prevents field rebinding. It does not make the nested pandas
    DataFrames intrinsically immutable, so callers must not modify them in place.
    """

    expression: pd.DataFrame
    metadata: pd.DataFrame
    de_results: pd.DataFrame
    source: DatasetSource
    source_label: str
    validation_report: ValidationReport

    @property
    def gene_count(self) -> int:
        return len(self.expression)

    @property
    def sample_count(self) -> int:
        return len(self.metadata)

    @property
    def de_row_count(self) -> int:
        return len(self.de_results)


@dataclass(frozen=True)
class CandidateResult:
    """Outcome of evaluating one complete or incomplete dataset candidate."""

    status: CandidateStatus
    tables: DatasetTables | None
    report: ValidationReport

    def __post_init__(self) -> None:
        if self.status == "incomplete":
            if self.tables is not None or self.report.has_errors:
                raise ValueError(
                    "An incomplete candidate must have no tables and no Errors."
                )
            return

        if self.status == "invalid":
            if self.tables is not None or not self.report.has_errors:
                raise ValueError(
                    "An invalid candidate must hide its tables and contain Errors."
                )
            return

        if self.status == "valid":
            if self.tables is None or self.report.has_errors:
                raise ValueError(
                    "A valid candidate must contain tables and have no Errors."
                )
            if len(self.tables) != 3 or not all(
                isinstance(table, pd.DataFrame) for table in self.tables
            ):
                raise ValueError(
                    "A valid candidate must contain three pandas DataFrames."
                )
            return

        raise ValueError(f"Unsupported candidate status: {self.status!r}.")


def load_demo_candidate(demo_dir: Path | None = None) -> CandidateResult:
    """Load and validate the three committed synthetic demonstration tables."""

    directory = DEMO_DIRECTORY if demo_dir is None else Path(demo_dir)
    sources = tuple(
        (key, table_name, directory / filename)
        for key, table_name, filename in _TABLE_SPECS
    )

    missing_issues = tuple(
        ValidationIssue(
            code=IssueCode.DEMO_FILE_MISSING,
            severity=Severity.ERROR,
            table=table_name,
            message=(
                f"Required synthetic demo file '{path.name}' is missing. "
                "Reinstall or restore the committed demo data."
            ),
            example_values=(path.name,),
        )
        for _, table_name, path in sources
        if not path.is_file()
    )
    if missing_issues:
        return CandidateResult(
            status="invalid",
            tables=None,
            report=ValidationReport(missing_issues),
        )

    return _load_and_validate_sources(
        sources,
        missing_code=IssueCode.DEMO_FILE_MISSING,
    )


def load_uploaded_candidate(
    expression_source: CsvSource | None,
    metadata_source: CsvSource | None,
    de_results_source: CsvSource | None,
) -> CandidateResult:
    """Load and validate uploaded sources, or report an incomplete selection."""

    supplied_sources = (
        expression_source,
        metadata_source,
        de_results_source,
    )
    if any(source is None for source in supplied_sources):
        return CandidateResult(
            status="incomplete",
            tables=None,
            report=ValidationReport(),
        )

    sources = tuple(
        (key, table_name, source)
        for (key, table_name, _), source in zip(
            _TABLE_SPECS,
            supplied_sources,
            strict=True,
        )
    )
    return _load_and_validate_sources(
        sources,
        missing_code=IssueCode.CSV_READ_ERROR,
    )


def build_dataset_bundle(
    tables: DatasetTables,
    *,
    source: DatasetSource,
    source_label: str,
    report: ValidationReport,
) -> DatasetBundle:
    """Build a complete bundle only from a non-blocking validation report."""

    if report.has_errors:
        raise ValueError("Cannot build a dataset bundle from a report with Errors.")
    if source not in ("demo", "uploaded"):
        raise ValueError(f"Unsupported dataset source: {source!r}.")
    if len(tables) != 3 or not all(
        isinstance(table, pd.DataFrame) for table in tables
    ):
        raise ValueError("A dataset bundle requires three pandas DataFrames.")
    if not source_label.strip():
        raise ValueError("A dataset bundle requires a source label.")

    expression, metadata, de_results = tables
    return DatasetBundle(
        expression=expression,
        metadata=metadata,
        de_results=de_results,
        source=source,
        source_label=source_label,
        validation_report=report,
    )


def uploaded_source_label(
    expression_source: CsvSource,
    metadata_source: CsvSource,
    de_results_source: CsvSource,
) -> str:
    """Return a human-readable label without exposing local directory paths."""

    names = [
        _display_source_name(source)
        for source in (
            expression_source,
            metadata_source,
            de_results_source,
        )
    ]
    visible_names = [name for name in names if name is not None]
    if not visible_names:
        return "User-uploaded CSV tables"
    return "User-uploaded CSV tables (" + ", ".join(visible_names) + ")"


def get_current_dataset(state: Mapping[str, Any]) -> DatasetBundle | None:
    """Return the current bundle, or ``None`` when no data are active."""

    bundle = state.get(CURRENT_DATASET_KEY)
    return bundle if isinstance(bundle, DatasetBundle) else None


def set_current_dataset(
    state: MutableMapping[str, Any],
    bundle: DatasetBundle,
) -> None:
    """Replace the single current bundle and clear failed-candidate feedback."""

    if bundle.validation_report.has_errors:
        raise ValueError("Cannot activate a dataset bundle with Errors.")
    state[CURRENT_DATASET_KEY] = bundle
    clear_candidate_feedback(state)


def record_candidate_failure(
    state: MutableMapping[str, Any],
    *,
    source: DatasetSource,
    source_label: str,
    report: ValidationReport,
) -> None:
    """Record a failed candidate without disturbing the current valid bundle."""

    if not report.has_errors:
        raise ValueError("Candidate failure feedback requires at least one Error.")
    if source not in ("demo", "uploaded"):
        raise ValueError(f"Unsupported candidate source: {source!r}.")
    if not source_label.strip():
        raise ValueError("Candidate failure feedback requires a source label.")

    state[CANDIDATE_REPORT_KEY] = report
    state[CANDIDATE_SOURCE_KEY] = source
    state[CANDIDATE_LABEL_KEY] = source_label


def clear_candidate_feedback(state: MutableMapping[str, Any]) -> None:
    """Remove only the failed-candidate report and its source metadata."""

    for key in CANDIDATE_FEEDBACK_KEYS:
        state.pop(key, None)


def clear_uploader_state(state: MutableMapping[str, Any]) -> None:
    """Remove the three uploader widget values before widget instantiation."""

    for key in UPLOADER_KEYS:
        state.pop(key, None)


def clear_legacy_data_state(state: MutableMapping[str, Any]) -> None:
    """Remove the three independent session keys used before Phase 4."""

    for key in LEGACY_DATA_KEYS:
        state.pop(key, None)


def reset_data_state(state: MutableMapping[str, Any]) -> None:
    """Delete only the explicit application-owned data and uploader keys."""

    for key in APPLICATION_DATA_KEYS:
        state.pop(key, None)


def table_preview(table: pd.DataFrame, rows: int = 5) -> pd.DataFrame:
    """Return an order-preserving copy of the first requested rows."""

    if rows <= 0:
        return table.iloc[0:0].copy(deep=True)
    return table.iloc[:rows].copy(deep=True)


def _load_and_validate_sources(
    sources: tuple[tuple[str, str, CsvSource], ...],
    *,
    missing_code: IssueCode,
) -> CandidateResult:
    tables: dict[str, pd.DataFrame] = {}
    issues: list[ValidationIssue] = []

    for key, table_name, source in sources:
        cursor_position = _cursor_position(source)
        try:
            tables[key] = read_csv(source)
        except CsvReadError as exc:
            issues.append(
                _read_issue(
                    table_name,
                    source,
                    exc.code or IssueCode.CSV_READ_ERROR,
                    str(exc),
                )
            )
        except FileNotFoundError:
            issues.append(
                _read_issue(
                    table_name,
                    source,
                    missing_code,
                    "The required CSV file could not be found.",
                )
            )
        except (PermissionError, IsADirectoryError) as exc:
            issues.append(
                _read_issue(
                    table_name,
                    source,
                    IssueCode.CSV_READ_ERROR,
                    f"The CSV file could not be read: {exc.strerror or str(exc)}",
                )
            )
        finally:
            _restore_cursor_position(source, cursor_position)

    if issues:
        return CandidateResult(
            status="invalid",
            tables=None,
            report=ValidationReport(tuple(issues)),
        )

    dataset_tables = (
        tables["expression"],
        tables["metadata"],
        tables["de_results"],
    )
    report = validate_input_tables(*dataset_tables)
    if report.has_errors:
        return CandidateResult(status="invalid", tables=None, report=report)
    return CandidateResult(status="valid", tables=dataset_tables, report=report)


def _read_issue(
    table_name: str,
    source: CsvSource,
    code: IssueCode,
    detail: str,
) -> ValidationIssue:
    source_name = _display_source_name(source)
    prefix = f"'{source_name}': " if source_name is not None else ""
    return ValidationIssue(
        code=code,
        severity=Severity.ERROR,
        table=table_name,
        message=prefix + detail,
        example_values=(source_name,) if source_name is not None else (),
    )


def _display_source_name(source: CsvSource) -> str | None:
    name = getattr(source, "name", None)
    if name is None and isinstance(source, (str, Path)):
        name = str(source)
    if name is None:
        return None
    return Path(str(name)).name


def _cursor_position(source: CsvSource) -> int | None:
    if not hasattr(source, "tell"):
        return None
    try:
        return source.tell()
    except (OSError, ValueError):
        return None


def _restore_cursor_position(source: CsvSource, position: int | None) -> None:
    if position is None or not hasattr(source, "seek"):
        return
    try:
        source.seek(position)
    except (OSError, ValueError):
        pass
