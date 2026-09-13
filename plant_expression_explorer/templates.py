"""Illustrative example CSV templates for the Upload Data page.

Every value is a fabricated placeholder that exists only to demonstrate the
required column names and layout; templates are never loaded as a dataset.
"""

from __future__ import annotations

import pandas as pd


def build_example_expression_template() -> pd.DataFrame:
    """Return a small illustrative expression-matrix template."""

    return pd.DataFrame(
        {
            "gene_id": ["GENE_EXAMPLE_1", "GENE_EXAMPLE_2", "GENE_EXAMPLE_3"],
            "Sample_1": [10.2, 3.4, 6.1],
            "Sample_2": [9.8, 3.1, 6.4],
            "Sample_3": [15.1, 3.6, 5.9],
            "Sample_4": [14.7, 3.9, 6.0],
        }
    )


def build_example_metadata_template() -> pd.DataFrame:
    """Return a small illustrative sample-metadata template."""

    return pd.DataFrame(
        {
            "sample_id": ["Sample_1", "Sample_2", "Sample_3", "Sample_4"],
            "condition": ["Control", "Control", "Treatment", "Treatment"],
            "genotype": ["WT", "WT", "WT", "WT"],
            "biological_replicate": [1, 2, 1, 2],
        }
    )


def build_example_de_results_template() -> pd.DataFrame:
    """Return a small illustrative differential-expression-results template."""

    return pd.DataFrame(
        {
            "gene_id": ["GENE_EXAMPLE_1", "GENE_EXAMPLE_2", "GENE_EXAMPLE_3"],
            "log2FoldChange": [1.25, -0.42, 0.03],
            "pvalue": [0.01, 0.35, 0.98],
            "padj": [0.03, 0.51, 0.98],
        }
    )
