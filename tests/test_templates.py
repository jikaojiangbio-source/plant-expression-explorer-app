"""Tests that example CSV templates are valid, loadable, and match the app contract."""

from io import StringIO

from plant_expression_explorer.data import read_csv
from plant_expression_explorer.templates import (
    build_example_de_results_template,
    build_example_expression_template,
    build_example_metadata_template,
)
from plant_expression_explorer.validation import (
    validate_de_results,
    validate_expression_matrix,
    validate_sample_metadata,
)


def _round_trip_csv(builder) -> None:
    csv_text = builder().to_csv(index=False)
    return read_csv(StringIO(csv_text))


def test_expression_template_is_a_valid_expression_matrix() -> None:
    table = _round_trip_csv(build_example_expression_template)

    report = validate_expression_matrix(table)

    assert not report.has_errors


def test_metadata_template_is_valid_sample_metadata() -> None:
    table = _round_trip_csv(build_example_metadata_template)

    report = validate_sample_metadata(table)

    assert not report.has_errors


def test_de_results_template_is_valid() -> None:
    table = _round_trip_csv(build_example_de_results_template)

    report = validate_de_results(table)

    assert not report.has_errors


def test_expression_and_metadata_templates_share_sample_ids() -> None:
    expression_columns = set(build_example_expression_template().columns) - {"gene_id"}
    metadata_sample_ids = set(build_example_metadata_template()["sample_id"])

    assert expression_columns == metadata_sample_ids


def test_expression_and_de_templates_share_gene_ids() -> None:
    expression_genes = set(build_example_expression_template()["gene_id"])
    de_genes = set(build_example_de_results_template()["gene_id"])

    assert expression_genes == de_genes
