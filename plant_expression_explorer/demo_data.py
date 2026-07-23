"""Deterministic synthetic data for software demonstration and testing."""

from __future__ import annotations

import math
import random
from pathlib import Path
from typing import Sequence

import pandas as pd

DEMO_SEED = 20260722
DECIMAL_PLACES = 8

CONTROL_SAMPLES = ("Control_1", "Control_2", "Control_3")
HIGH_NITRATE_SAMPLES = (
    "High_nitrate_1",
    "High_nitrate_2",
    "High_nitrate_3",
)
SAMPLE_IDS = CONTROL_SAMPLES + HIGH_NITRATE_SAMPLES
GENE_IDS = tuple(f"SYN_Solyc_{index:04d}" for index in range(1, 121))

EXPRESSION_COLUMNS = ("gene_id", *SAMPLE_IDS)
METADATA_COLUMNS = ("sample_id", "condition")
DEG_COLUMNS = ("gene_id", "log2FoldChange", "pvalue", "padj")
OUTPUT_FILENAMES = (
    "expression_matrix.csv",
    "metadata.csv",
    "deg_results.csv",
)


def generate_demo_data(
    seed: int = DEMO_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Return expression, metadata, and DEG tables in that order."""

    rng = random.Random(seed)
    expression = _generate_expression_matrix(rng)
    metadata = _build_metadata()
    deg_results = _build_deg_results(expression)
    return expression, metadata, deg_results


def write_demo_data(
    output_dir: Path,
    seed: int = DEMO_SEED,
) -> tuple[Path, Path, Path]:
    """Generate and write the three demo CSV files to ``output_dir``."""

    resolved_output_dir = Path(output_dir).expanduser().resolve()
    resolved_output_dir.mkdir(parents=True, exist_ok=True)

    paths = tuple(resolved_output_dir / filename for filename in OUTPUT_FILENAMES)
    for path in paths:
        if path.is_symlink():
            raise ValueError(f"Refusing to replace symbolic link: {path}")
        if path.exists() and not path.is_file():
            raise ValueError(f"Refusing to replace non-file path: {path}")

    tables = generate_demo_data(seed=seed)
    for table, path in zip(tables, paths, strict=True):
        table.to_csv(
            path,
            index=False,
            float_format=f"%.{DECIMAL_PLACES}f",
            lineterminator="\n",
        )
    return paths


def _generate_expression_matrix(rng: random.Random) -> pd.DataFrame:
    rows: list[dict[str, str | float]] = []
    for position, gene_id in enumerate(GENE_IDS, start=1):
        baseline = rng.uniform(6.0, 12.0)
        if position <= 20:
            effect = rng.uniform(1.2, 2.5)
        elif position <= 40:
            effect = -rng.uniform(1.2, 2.5)
        else:
            effect = rng.uniform(-0.25, 0.25)

        control_noise = _condition_centred_noise(rng)
        high_nitrate_noise = _condition_centred_noise(rng)
        values = (
            *(baseline + noise for noise in control_noise),
            *(baseline + effect + noise for noise in high_nitrate_noise),
        )
        row: dict[str, str | float] = {"gene_id": gene_id}
        row.update(
            {
                sample_id: _quantise(value)
                for sample_id, value in zip(SAMPLE_IDS, values, strict=True)
            }
        )
        rows.append(row)

    return pd.DataFrame(rows, columns=EXPRESSION_COLUMNS)


def _condition_centred_noise(rng: random.Random) -> tuple[float, float, float]:
    noise = [rng.uniform(-0.4, 0.4) for _ in range(3)]
    centre = sum(noise) / len(noise)
    return (
        noise[0] - centre,
        noise[1] - centre,
        noise[2] - centre,
    )


def _build_metadata() -> pd.DataFrame:
    rows = [
        {"sample_id": sample_id, "condition": "Control"}
        for sample_id in CONTROL_SAMPLES
    ]
    rows.extend(
        {"sample_id": sample_id, "condition": "High_nitrate"}
        for sample_id in HIGH_NITRATE_SAMPLES
    )
    return pd.DataFrame(rows, columns=METADATA_COLUMNS)


def _build_deg_results(expression: pd.DataFrame) -> pd.DataFrame:
    control_means = expression.loc[:, CONTROL_SAMPLES].mean(axis=1)
    high_nitrate_means = expression.loc[:, HIGH_NITRATE_SAMPLES].mean(axis=1)
    fold_changes = [
        _quantise(value)
        for value in (high_nitrate_means - control_means).tolist()
    ]

    pvalues = [
        _quantise(_synthetic_pvalue(position, fold_change))
        for position, fold_change in enumerate(fold_changes, start=1)
    ]
    adjusted_pvalues = [
        _quantise(value) for value in _benjamini_hochberg(pvalues)
    ]

    return pd.DataFrame(
        {
            "gene_id": expression["gene_id"].tolist(),
            "log2FoldChange": fold_changes,
            "pvalue": pvalues,
            "padj": adjusted_pvalues,
        },
        columns=DEG_COLUMNS,
    )


def _synthetic_pvalue(position: int, fold_change: float) -> float:
    magnitude = abs(fold_change)
    if position <= 40:
        scaled_effect = _bounded((magnitude - 1.2) / (2.5 - 1.2))
        return 0.012 * math.exp(-3.7 * scaled_effect)

    scaled_effect = _bounded(magnitude / 0.25)
    return 0.95 - (0.70 * scaled_effect)


def _benjamini_hochberg(pvalues: Sequence[float]) -> list[float]:
    """Return Benjamini-Hochberg adjusted p-values in original order."""

    values = [float(value) for value in pvalues]
    if any(not math.isfinite(value) or not 0.0 <= value <= 1.0 for value in values):
        raise ValueError("p-values must be finite values within [0, 1].")
    if not values:
        return []

    count = len(values)
    order = sorted(range(count), key=lambda index: (values[index], index))
    adjusted_sorted = [0.0] * count
    running_minimum = 1.0

    for sorted_index in range(count - 1, -1, -1):
        original_index = order[sorted_index]
        rank = sorted_index + 1
        candidate = values[original_index] * count / rank
        running_minimum = min(running_minimum, candidate)
        adjusted_sorted[sorted_index] = max(0.0, min(1.0, running_minimum))

    adjusted = [0.0] * count
    for sorted_index, original_index in enumerate(order):
        adjusted[original_index] = adjusted_sorted[sorted_index]
    return adjusted


def _bounded(value: float) -> float:
    return max(0.0, min(1.0, value))


def _quantise(value: float) -> float:
    return float(f"{value:.{DECIMAL_PLACES}f}")
