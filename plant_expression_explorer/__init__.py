"""Core validation functionality for Plant Expression Explorer."""

from plant_expression_explorer.consistency import (
    validate_expression_de_consistency,
    validate_expression_metadata_consistency,
    validate_input_tables,
)
from plant_expression_explorer.data import CsvReadError, read_csv
from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
    validate_de_results,
    validate_expression_matrix,
    validate_sample_metadata,
)

__all__ = [
    "CsvReadError",
    "IssueCode",
    "Severity",
    "ValidationIssue",
    "ValidationReport",
    "read_csv",
    "validate_de_results",
    "validate_expression_de_consistency",
    "validate_expression_matrix",
    "validate_expression_metadata_consistency",
    "validate_input_tables",
    "validate_sample_metadata",
]
