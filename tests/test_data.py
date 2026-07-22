"""Tests for the Phase 1 input contracts."""

from io import StringIO

import pandas as pd
import pytest

from plant_expression_explorer.data import (
    DataValidationError,
    align_expression_and_metadata,
    read_csv,
    validate_de_results,
    validate_expression_matrix,
)


def test_valid_tables_are_read_validated_and_aligned() -> None:
    expression = read_csv(StringIO("gene_id,s2,s1\nSolyc01g1,2.0,1.0\nSolyc01g2,4.0,3.0\n"))
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["control", "treated"]}
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["Solyc01g1", "Solyc01g2"],
            "log2FoldChange": [1.2, -0.5],
            "pvalue": [0.001, None],
            "padj": [0.01, None],
        }
    )

    aligned = align_expression_and_metadata(expression, metadata)
    validate_de_results(de_results)

    assert aligned["sample_id"].tolist() == ["s2", "s1"]


def test_expression_rejects_non_numeric_values() -> None:
    expression = pd.DataFrame({"gene_id": ["gene1"], "sample1": ["high"]})

    with pytest.raises(DataValidationError, match="numeric"):
        validate_expression_matrix(expression)


def test_expression_rejects_duplicate_gene_ids() -> None:
    expression = pd.DataFrame(
        {"gene_id": ["gene1", "gene1"], "sample1": [1.0, 2.0]}
    )

    with pytest.raises(DataValidationError, match="duplicate identifiers"):
        validate_expression_matrix(expression)


def test_alignment_rejects_mismatched_sample_ids() -> None:
    expression = pd.DataFrame({"gene_id": ["gene1"], "sample1": [1.0]})
    metadata = pd.DataFrame(
        {"sample_id": ["different_sample"], "condition": ["control"]}
    )

    with pytest.raises(DataValidationError, match="do not match"):
        align_expression_and_metadata(expression, metadata)


@pytest.mark.parametrize("column", ["pvalue", "padj"])
def test_de_results_reject_out_of_range_probabilities(column: str) -> None:
    de_results = pd.DataFrame(
        {
            "gene_id": ["gene1"],
            "log2FoldChange": [1.0],
            "pvalue": [0.1],
            "padj": [0.2],
        }
    )
    de_results.loc[0, column] = 1.5

    with pytest.raises(DataValidationError, match="between 0 and 1"):
        validate_de_results(de_results)
