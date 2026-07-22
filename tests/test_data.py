"""Tests for CSV reading and raw-header integrity."""

from io import StringIO

import pytest

from plant_expression_explorer.data import CsvReadError, read_csv
from plant_expression_explorer.validation import IssueCode


def test_read_csv_reads_valid_table() -> None:
    table = read_csv(StringIO("gene_id,s1\ng1,1.0\n"))

    assert table.columns.tolist() == ["gene_id", "s1"]
    assert table.loc[0, "gene_id"] == "g1"


def test_read_csv_detects_duplicate_headers_before_pandas_mangling() -> None:
    source = StringIO("gene_id,s1,s1\ng1,1.0,2.0\n")

    with pytest.raises(CsvReadError) as caught:
        read_csv(source)

    assert caught.value.code is IssueCode.DUPLICATE_COLUMN_NAME
    assert "s1" in str(caught.value)


def test_read_csv_reports_empty_file() -> None:
    with pytest.raises(CsvReadError, match="empty"):
        read_csv(StringIO(""))


def test_read_csv_reports_parser_failure() -> None:
    with pytest.raises(CsvReadError, match="Could not read"):
        read_csv(StringIO('a,b\n"unterminated,1\n'))
