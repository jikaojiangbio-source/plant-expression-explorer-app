"""Regression tests for deterministic descriptive CSV exports."""

import csv
import io
from dataclasses import FrozenInstanceError

import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import plant_expression_explorer as package
from plant_expression_explorer.exports import (
    CSV_MEDIA_TYPE,
    CsvExportArtifact,
    CsvExportError,
    CsvExportErrorReason,
    build_csv_export,
)


def _capture_error(
    table: object,
    *,
    filename: object = "result.csv",
    json_sequence_columns: object = (),
) -> CsvExportError:
    with pytest.raises(CsvExportError) as captured:
        build_csv_export(
            table,
            filename=filename,
            json_sequence_columns=json_sequence_columns,
        )
    return captured.value


def test_package_exports_the_phase_10_public_api() -> None:
    assert package.CSV_MEDIA_TYPE == CSV_MEDIA_TYPE
    assert package.CsvExportArtifact is CsvExportArtifact
    assert package.CsvExportError is CsvExportError
    assert package.CsvExportErrorReason is CsvExportErrorReason
    assert package.build_csv_export is build_csv_export


def test_error_reasons_are_stable_and_unique() -> None:
    values = [reason.value for reason in CsvExportErrorReason]
    assert len(values) == len(set(values))


def test_export_is_utf8_deterministic_unrounded_and_csv_safe() -> None:
    table = pd.DataFrame(
        {
            "sample_id": ["样本,1", 'sample"2', "sample\n3"],
            "value": [1e308, np.nextafter(1.0, 2.0), np.nan],
            "defined": [True, False, pd.NA],
        },
        index=pd.Index([9, 4, 1], name="source_index"),
    )
    original = table.copy(deep=True)

    first = build_csv_export(table, filename="sample-results.csv")
    second = build_csv_export(table, filename="sample-results.csv")

    assert first == second
    assert first.filename == "sample-results.csv"
    assert first.media_type == "text/csv; charset=utf-8"
    assert first.row_count == 3
    assert first.column_count == 3
    assert first.column_names == ("sample_id", "value", "defined")
    assert first.data.startswith(b"sample_id,value,defined\n")
    assert b"source_index" not in first.data
    rows = list(csv.reader(io.StringIO(first.data.decode("utf-8"))))
    assert rows[1] == ["样本,1", "1e+308", "True"]
    assert rows[2][0] == 'sample"2'
    assert float(rows[2][1]) == np.nextafter(1.0, 2.0)
    assert rows[2][2] == "False"
    assert rows[3] == ["sample\n3", "", ""]
    assert_frame_equal(table, original, check_exact=True)


def test_sequence_columns_are_explicit_compact_json_in_a_copy() -> None:
    table = pd.DataFrame(
        {
            "condition": ["A", "B"],
            "sample_ids": [("s1", "s,2"), ["样本3"]],
            "sample_count": [2, 1],
        }
    )
    original = table.copy(deep=True)

    artifact = build_csv_export(
        table,
        filename="condition-summary.csv",
        json_sequence_columns=("sample_ids",),
    )

    rows = list(csv.DictReader(io.StringIO(artifact.data.decode("utf-8"))))
    assert rows[0]["sample_ids"] == '["s1","s,2"]'
    assert rows[1]["sample_ids"] == '["样本3"]'
    assert_frame_equal(table, original, check_exact=True)


def test_formula_like_text_is_preserved_without_silent_prefixing() -> None:
    values = ["=1+1", "+command", "-identifier", "@function", "\t=1+1"]
    table = pd.DataFrame({"annotation": values})

    artifact = build_csv_export(table, filename="exact-text.csv")

    rows = list(csv.DictReader(io.StringIO(artifact.data.decode("utf-8"))))
    assert [row["annotation"] for row in rows] == values


def test_empty_table_exports_its_header_without_inventing_rows() -> None:
    table = pd.DataFrame(columns=["gene_id", "status"])

    artifact = build_csv_export(table, filename="empty-result.csv")

    assert artifact.data == b"gene_id,status\n"
    assert artifact.row_count == 0


@pytest.mark.parametrize("table", [None, [], {"a": [1]}])
def test_non_dataframe_inputs_are_controlled(table: object) -> None:
    assert _capture_error(table).reason is CsvExportErrorReason.INVALID_TABLE


@pytest.mark.parametrize(
    "filename",
    [None, "", "result", "../result.csv", "/result.csv", "result/a.csv", True],
)
def test_unsafe_or_invalid_filenames_are_controlled(filename: object) -> None:
    error = _capture_error(pd.DataFrame({"a": [1]}), filename=filename)
    assert error.reason is CsvExportErrorReason.INVALID_FILENAME


@pytest.mark.parametrize(
    ("columns", "reason"),
    [
        (["", "value"], CsvExportErrorReason.INVALID_COLUMNS),
        ([1, "value"], CsvExportErrorReason.INVALID_COLUMNS),
        (["value", "value"], CsvExportErrorReason.DUPLICATE_COLUMNS),
    ],
)
def test_invalid_columns_are_not_renamed_or_removed(
    columns: list[object],
    reason: CsvExportErrorReason,
) -> None:
    table = pd.DataFrame([[1, 2]], columns=columns)
    original = table.copy(deep=True)

    assert _capture_error(table).reason is reason
    assert_frame_equal(table, original, check_exact=True)


def test_table_without_columns_is_controlled() -> None:
    table = pd.DataFrame(index=[0])

    error = _capture_error(table)

    assert error.reason is CsvExportErrorReason.INVALID_COLUMNS


@pytest.mark.parametrize(
    "sequence_columns",
    [None, ("missing",), ("sample_ids", "sample_ids"), (1,)],
)
def test_invalid_sequence_column_requests_are_controlled(
    sequence_columns: object,
) -> None:
    table = pd.DataFrame({"sample_ids": [("s1",)]})
    error = _capture_error(
        table,
        json_sequence_columns=sequence_columns,
    )
    assert error.reason is CsvExportErrorReason.INVALID_SEQUENCE_COLUMN


@pytest.mark.parametrize("value", ["s1", ("s1", 2), {"s1", "s2"}])
def test_invalid_sequence_values_are_controlled_without_mutation(value: object) -> None:
    table = pd.DataFrame({"sample_ids": [value]})
    original = table.copy(deep=True)

    error = _capture_error(
        table,
        json_sequence_columns=("sample_ids",),
    )

    assert error.reason is CsvExportErrorReason.UNSUPPORTED_VALUE
    assert_frame_equal(table, original, check_exact=True)


@pytest.mark.parametrize(
    ("value", "reason"),
    [
        (object(), CsvExportErrorReason.UNSUPPORTED_VALUE),
        (1 + 2j, CsvExportErrorReason.UNSUPPORTED_VALUE),
        (float("inf"), CsvExportErrorReason.INFINITE_VALUE),
        (float("-inf"), CsvExportErrorReason.INFINITE_VALUE),
    ],
)
def test_unsupported_and_infinite_values_produce_no_partial_export(
    value: object,
    reason: CsvExportErrorReason,
) -> None:
    table = pd.DataFrame({"value": [1, value]}, dtype="object")
    original = table.copy(deep=True)

    error = _capture_error(table)

    assert error.reason is reason
    assert_frame_equal(table, original, check_exact=True)


@pytest.mark.parametrize("value", [True, False, np.bool_(True), np.bool_(False)])
def test_python_and_numpy_boolean_values_are_exported_without_coercing_input(
    value: object,
) -> None:
    table = pd.DataFrame({"defined": [value]}, dtype="object")
    original = table.copy(deep=True)

    artifact = build_csv_export(table, filename="boolean.csv")

    assert artifact.data.decode("utf-8").splitlines()[1] == str(bool(value))
    assert_frame_equal(table, original, check_exact=True)


def test_serialization_failure_is_controlled(monkeypatch: pytest.MonkeyPatch) -> None:
    table = pd.DataFrame({"value": [1]})

    def fail_serialization(*args: object, **kwargs: object) -> str:
        raise UnicodeError("test failure")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_serialization)

    error = _capture_error(table)

    assert error.reason is CsvExportErrorReason.SERIALIZATION_ERROR


def test_artifact_fields_cannot_be_rebound() -> None:
    artifact = build_csv_export(pd.DataFrame({"value": [1]}), filename="x.csv")

    with pytest.raises(FrozenInstanceError):
        artifact.filename = "changed.csv"  # type: ignore[misc]
