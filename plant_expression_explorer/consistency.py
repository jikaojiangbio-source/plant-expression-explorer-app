"""Non-mutating consistency checks across the biological input tables."""

from __future__ import annotations

import pandas as pd

from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
    validate_de_results,
    validate_expression_matrix,
    validate_sample_metadata,
)


def validate_expression_metadata_consistency(
    expression: pd.DataFrame, metadata: pd.DataFrame
) -> ValidationReport:
    """Compare expression sample columns with metadata sample identifiers."""

    if not isinstance(expression, pd.DataFrame) or not isinstance(metadata, pd.DataFrame):
        return ValidationReport()
    if not _has_one_column(expression, "gene_id") or not _has_one_column(
        metadata, "sample_id"
    ):
        return ValidationReport()
    if expression.columns.duplicated().any():
        return ValidationReport()

    expression_samples = [str(column) for column in expression.columns if column != "gene_id"]
    metadata_samples = [str(value) for value in metadata["sample_id"].dropna().tolist()]
    expression_set = set(expression_samples)
    metadata_set = set(metadata_samples)

    if expression_set and metadata_set and expression_set.isdisjoint(metadata_set):
        return ValidationReport(
            (
                ValidationIssue(
                    code=IssueCode.NO_SAMPLE_OVERLAP,
                    severity=Severity.ERROR,
                    table="Expression matrix",
                    related_table="Sample metadata",
                    message="The expression matrix and sample metadata share no sample identifiers.",
                ),
            )
        )

    issues: list[ValidationIssue] = []
    missing_metadata = sorted(expression_set - metadata_set)
    unexpected_metadata = sorted(metadata_set - expression_set)
    if missing_metadata:
        issues.append(
            _set_issue(
                IssueCode.MISSING_METADATA_SAMPLE,
                Severity.ERROR,
                "Expression matrix",
                "Sample metadata",
                f"Sample metadata is missing {len(missing_metadata)} expression sample(s): "
                + _format_examples(missing_metadata)
                + ".",
                missing_metadata,
            )
        )
    if unexpected_metadata:
        issues.append(
            _set_issue(
                IssueCode.UNEXPECTED_METADATA_SAMPLE,
                Severity.ERROR,
                "Sample metadata",
                "Expression matrix",
                f"Sample metadata contains {len(unexpected_metadata)} sample(s) absent from the expression matrix: "
                + _format_examples(unexpected_metadata)
                + ".",
                unexpected_metadata,
            )
        )

    identifiers_are_unique = (
        len(expression_samples) == len(expression_set)
        and len(metadata_samples) == len(metadata_set)
    )
    if (
        not missing_metadata
        and not unexpected_metadata
        and identifiers_are_unique
        and expression_samples != metadata_samples
    ):
        issues.append(
            ValidationIssue(
                code=IssueCode.SAMPLE_ORDER_DIFFERS,
                severity=Severity.INFORMATION,
                table="Expression matrix",
                related_table="Sample metadata",
                message="Expression samples and metadata samples contain the same identifiers in a different order; neither table was reordered.",
            )
        )
    return ValidationReport(tuple(issues))


def validate_expression_de_consistency(
    expression: pd.DataFrame, de_results: pd.DataFrame
) -> ValidationReport:
    """Compare expression and differential-expression gene identifiers."""

    if not isinstance(expression, pd.DataFrame) or not isinstance(de_results, pd.DataFrame):
        return ValidationReport()
    if not _has_one_column(expression, "gene_id") or not _has_one_column(
        de_results, "gene_id"
    ):
        return ValidationReport()

    expression_genes = _identifier_set(expression["gene_id"])
    de_genes = _identifier_set(de_results["gene_id"])
    if expression_genes and de_genes and expression_genes.isdisjoint(de_genes):
        return ValidationReport(
            (
                ValidationIssue(
                    code=IssueCode.NO_GENE_OVERLAP,
                    severity=Severity.ERROR,
                    table="Expression matrix",
                    related_table="Differential-expression results",
                    message=(
                        "The expression matrix and differential-expression results "
                        "share no gene identifiers. Common plant-genomics causes "
                        "include different annotation releases, gene-versus-transcript "
                        "identifiers, isoform or version suffixes, and alias conventions. "
                        "Validation did not rewrite any identifier."
                    ),
                ),
            )
        )

    issues: list[ValidationIssue] = []
    de_only = sorted(de_genes - expression_genes)
    expression_only = sorted(expression_genes - de_genes)
    if de_only:
        issues.append(
            _set_issue(
                IssueCode.DE_GENE_NOT_IN_EXPRESSION,
                Severity.WARNING,
                "Differential-expression results",
                "Expression matrix",
                f"{len(de_only)} differential-expression gene identifier(s) are absent from the expression matrix: "
                + _format_examples(de_only)
                + ". The rows have been retained. Common plant-genomics causes "
                "include different annotation releases, gene-versus-transcript "
                "identifiers, isoform or version suffixes, and alias conventions; "
                "validation did not rewrite any identifier.",
                de_only,
            )
        )
    if expression_only:
        issues.append(
            _set_issue(
                IssueCode.EXPRESSION_GENE_NOT_IN_DE,
                Severity.WARNING,
                "Expression matrix",
                "Differential-expression results",
                f"{len(expression_only)} expression gene identifier(s) are absent from the differential-expression results: "
                + _format_examples(expression_only)
                + ". No expression rows were removed. Common plant-genomics causes "
                "include different annotation releases, gene-versus-transcript "
                "identifiers, isoform or version suffixes, and alias conventions; "
                "validation did not rewrite any identifier.",
                expression_only,
            )
        )
    return ValidationReport(tuple(issues))


def validate_input_tables(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    de_results: pd.DataFrame | None,
) -> ValidationReport:
    """Combine table and applicable cross-file validation reports."""

    reports = [
        validate_expression_matrix(expression),
        validate_sample_metadata(metadata),
        validate_expression_metadata_consistency(expression, metadata),
    ]
    if de_results is None:
        reports.append(
            ValidationReport(
                (
                    ValidationIssue(
                        code=IssueCode.DE_RESULTS_NOT_SUPPLIED,
                        severity=Severity.INFORMATION,
                        table="Differential-expression results",
                        message=(
                            "No precomputed differential-expression results were "
                            "supplied. Expression- and metadata-based descriptive "
                            "pages remain available; Differential Expression is "
                            "unavailable for this dataset."
                        ),
                    ),
                )
            )
        )
    else:
        reports.extend(
            (
                validate_de_results(de_results),
                validate_expression_de_consistency(expression, de_results),
            )
        )
    return ValidationReport(tuple(issue for report in reports for issue in report.issues))


def _identifier_set(series: pd.Series) -> set[str]:
    return {
        str(value)
        for value in series.tolist()
        if pd.notna(value) and str(value).strip()
    }


def _has_one_column(table: pd.DataFrame, column: str) -> bool:
    return sum(existing == column for existing in table.columns) == 1


def _set_issue(
    code: IssueCode,
    severity: Severity,
    table: str,
    related_table: str,
    message: str,
    values: list[str],
) -> ValidationIssue:
    return ValidationIssue(
        code=code,
        severity=severity,
        table=table,
        related_table=related_table,
        message=message,
        count=len(values),
        example_values=tuple(values[:5]),
    )


def _format_examples(values: list[str]) -> str:
    return ", ".join(values[:5])
