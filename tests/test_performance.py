"""Regression guard against accidental non-linear blow-ups on large matrices.

These are coarse smoke tests, not benchmarks: the time budget is deliberately
generous so ordinary machine and CI variance never causes a false failure.
Their purpose is only to catch an accidentally quadratic-or-worse
implementation change on inputs sized like a realistic plant genome.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd
import pytest

from plant_expression_explorer.correlation import compute_sample_correlation
from plant_expression_explorer.pca import compute_sample_pca
from plant_expression_explorer.qc import compute_sample_qc

_GENE_COUNT = 20_000
_SAMPLE_COUNT = 80
_TIME_BUDGET_SECONDS = 20.0


@pytest.fixture(scope="module")
def large_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(0)
    sample_ids = [f"sample_{i}" for i in range(_SAMPLE_COUNT)]
    expression = pd.DataFrame(
        rng.normal(size=(_GENE_COUNT, _SAMPLE_COUNT)),
        columns=sample_ids,
    )
    expression.insert(0, "gene_id", [f"gene_{i}" for i in range(_GENE_COUNT)])
    metadata = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "condition": ["control", "treated"] * (_SAMPLE_COUNT // 2),
        }
    )
    return expression, metadata


def test_pca_completes_within_budget_on_a_genome_scale_matrix(
    large_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = large_tables

    started = time.perf_counter()
    result = compute_sample_pca(expression, metadata)
    elapsed = time.perf_counter() - started

    assert result.gene_count == _GENE_COUNT
    assert elapsed < _TIME_BUDGET_SECONDS


def test_correlation_completes_within_budget_on_a_genome_scale_matrix(
    large_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = large_tables

    started = time.perf_counter()
    result = compute_sample_correlation(expression, metadata)
    elapsed = time.perf_counter() - started

    assert len(result.correlation_matrix) == _SAMPLE_COUNT
    assert elapsed < _TIME_BUDGET_SECONDS


def test_qc_completes_within_budget_on_a_genome_scale_matrix(
    large_tables: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata = large_tables

    started = time.perf_counter()
    result = compute_sample_qc(expression, metadata)
    elapsed = time.perf_counter() - started

    assert len(result.sample_summary) == _SAMPLE_COUNT
    assert elapsed < _TIME_BUDGET_SECONDS
