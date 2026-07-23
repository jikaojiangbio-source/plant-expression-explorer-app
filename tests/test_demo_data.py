from __future__ import annotations

import math
import sys
from itertools import combinations
from pathlib import Path

import pandas as pd
import pytest

from scripts import generate_demo_data as demo_data_cli
from plant_expression_explorer.consistency import validate_input_tables
from plant_expression_explorer.data import read_csv
from plant_expression_explorer.demo_data import (
    DEMO_SEED,
    _benjamini_hochberg,
    generate_demo_data,
    write_demo_data,
)
from plant_expression_explorer.validation import (
    validate_de_results,
    validate_expression_matrix,
    validate_sample_metadata,
)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEMO_DIRECTORY = PROJECT_ROOT / "data" / "demo"
PRECISION = 8

EXPECTED_GENE_IDS = [f"SYN_Solyc_{index:04d}" for index in range(1, 121)]
CONTROL_SAMPLES = ["Control_1", "Control_2", "Control_3"]
HIGH_NITRATE_SAMPLES = [
    "High_nitrate_1",
    "High_nitrate_2",
    "High_nitrate_3",
]
SAMPLE_IDS = CONTROL_SAMPLES + HIGH_NITRATE_SAMPLES
EXPECTED_FILENAMES = (
    "expression_matrix.csv",
    "metadata.csv",
    "deg_results.csv",
)


@pytest.fixture(scope="module")
def demo_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return generate_demo_data()


@pytest.fixture(scope="module")
def bundled_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    return tuple(
        read_csv(DEMO_DIRECTORY / filename) for filename in EXPECTED_FILENAMES
    )


def test_exact_table_dimensions_and_column_order(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata, deg_results = demo_tables

    assert expression.shape == (120, 7)
    assert expression.columns.tolist() == ["gene_id", *SAMPLE_IDS]
    assert metadata.shape == (6, 2)
    assert metadata.columns.tolist() == ["sample_id", "condition"]
    assert deg_results.shape == (120, 4)
    assert deg_results.columns.tolist() == [
        "gene_id",
        "log2FoldChange",
        "pvalue",
        "padj",
    ]


def test_exact_sample_metadata_contract(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    _, metadata, _ = demo_tables
    expected = pd.DataFrame(
        {
            "sample_id": SAMPLE_IDS,
            "condition": ["Control"] * 3 + ["High_nitrate"] * 3,
        }
    )

    pd.testing.assert_frame_equal(metadata, expected, check_exact=True)


def test_exact_gene_identifiers_are_ordered_unique_and_complete(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    expression, _, deg_results = demo_tables

    for identifiers in (expression["gene_id"], deg_results["gene_id"]):
        assert identifiers.tolist() == EXPECTED_GENE_IDS
        assert identifiers.is_unique
        assert identifiers.notna().all()
        assert identifiers.str.strip().ne("").all()


def test_generation_is_reproducible_for_the_same_seed() -> None:
    first = generate_demo_data(seed=DEMO_SEED)
    second = generate_demo_data(seed=DEMO_SEED)

    for first_table, second_table in zip(first, second, strict=True):
        pd.testing.assert_frame_equal(first_table, second_table, check_exact=True)


def test_different_seeds_change_expression_and_deg_values() -> None:
    default_expression, default_metadata, default_deg = generate_demo_data(
        seed=DEMO_SEED
    )
    other_expression, other_metadata, other_deg = generate_demo_data(
        seed=DEMO_SEED + 1
    )

    assert not default_expression[SAMPLE_IDS].equals(other_expression[SAMPLE_IDS])
    assert not default_deg.drop(columns="gene_id").equals(
        other_deg.drop(columns="gene_id")
    )
    pd.testing.assert_frame_equal(default_metadata, other_metadata, check_exact=True)


def test_biological_replicate_columns_are_not_identical(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    expression, _, _ = demo_tables

    for group in (CONTROL_SAMPLES, HIGH_NITRATE_SAMPLES):
        for first, second in combinations(group, 2):
            assert not expression[first].equals(expression[second])


def test_numeric_values_are_finite_complete_and_expression_is_positive(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    expression, _, deg_results = demo_tables
    numeric_groups = (
        expression[SAMPLE_IDS],
        deg_results[["log2FoldChange", "pvalue", "padj"]],
    )

    for values in numeric_groups:
        assert all(pd.api.types.is_numeric_dtype(values[column]) for column in values)
        assert values.notna().all().all()
        assert values.map(math.isfinite).all().all()
    assert expression[SAMPLE_IDS].gt(0).all().all()


def test_fixed_gene_ranges_have_exact_synthetic_effect_contract(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    _, _, deg_results = demo_tables
    upregulated = deg_results.iloc[:20]
    downregulated = deg_results.iloc[20:40]
    unchanged = deg_results.iloc[40:]

    assert len(upregulated) == 20
    assert len(downregulated) == 20
    assert len(unchanged) == 80
    assert upregulated["gene_id"].tolist() == EXPECTED_GENE_IDS[:20]
    assert downregulated["gene_id"].tolist() == EXPECTED_GENE_IDS[20:40]
    assert unchanged["gene_id"].tolist() == EXPECTED_GENE_IDS[40:]
    assert upregulated["log2FoldChange"].between(1.2, 2.5).all()
    assert downregulated["log2FoldChange"].between(-2.5, -1.2).all()
    assert unchanged["log2FoldChange"].between(-0.25, 0.25).all()


def test_log2_fold_change_matches_quantised_expression_means_exactly(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    expression, _, deg_results = demo_tables
    derived = (
        expression[HIGH_NITRATE_SAMPLES].mean(axis=1)
        - expression[CONTROL_SAMPLES].mean(axis=1)
    ).map(lambda value: float(f"{value:.{PRECISION}f}"))

    pd.testing.assert_series_equal(
        deg_results["log2FoldChange"],
        derived.rename("log2FoldChange"),
        check_exact=True,
    )


def test_pvalues_and_adjusted_pvalues_are_strictly_positive_and_bounded(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    _, _, deg_results = demo_tables

    for column in ("pvalue", "padj"):
        assert deg_results[column].gt(0).all()
        assert deg_results[column].le(1).all()


def test_benjamini_hochberg_known_vector_handles_ties_and_restores_order() -> None:
    pvalues = [0.04, 0.01, 0.01, 0.20]

    assert _benjamini_hochberg(pvalues) == pytest.approx(
        [0.05333333333333334, 0.02, 0.02, 0.20]
    )


def test_adjusted_pvalues_are_monotonic_when_sorted_by_pvalue(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    _, _, deg_results = demo_tables
    ordered = deg_results.sort_values("pvalue", kind="stable")
    adjusted = ordered["padj"].tolist()

    assert all(first <= second for first, second in zip(adjusted, adjusted[1:]))


def test_default_dataset_has_useful_positive_and_negative_demo_points(
    demo_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    _, _, deg_results = demo_tables
    passes_thresholds = (
        deg_results["log2FoldChange"].abs().ge(1)
        & deg_results["padj"].le(0.05)
    )

    assert (passes_thresholds & deg_results["log2FoldChange"].gt(0)).sum() >= 10
    assert (passes_thresholds & deg_results["log2FoldChange"].lt(0)).sum() >= 10


def test_bundled_tables_pass_all_phase_2_single_table_validation(
    bundled_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    expression, metadata, deg_results = bundled_tables
    reports = (
        validate_expression_matrix(expression),
        validate_sample_metadata(metadata),
        validate_de_results(deg_results),
    )

    assert all(report.issues == () for report in reports)


def test_bundled_tables_pass_phase_2_aggregate_validation_without_issues(
    bundled_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    report = validate_input_tables(*bundled_tables)

    assert report.issues == ()
    assert report.errors == ()
    assert report.warnings == ()
    assert report.information == ()


def test_bundled_csv_files_match_fresh_default_generation(
    bundled_tables: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame],
) -> None:
    regenerated = generate_demo_data()

    for bundled, fresh in zip(bundled_tables, regenerated, strict=True):
        pd.testing.assert_frame_equal(bundled, fresh, check_exact=True)


def test_repeated_default_writes_are_byte_identical(tmp_path: Path) -> None:
    first_paths = write_demo_data(tmp_path / "first")
    second_paths = write_demo_data(tmp_path / "second")

    for first, second in zip(first_paths, second_paths, strict=True):
        assert first.read_bytes() == second.read_bytes()


def test_writer_uses_supplied_seed_and_returns_only_intended_paths(
    tmp_path: Path,
) -> None:
    output_directory = tmp_path / "custom"
    paths = write_demo_data(output_directory, seed=12345)
    written = tuple(read_csv(path) for path in paths)
    expected = generate_demo_data(seed=12345)

    assert tuple(path.name for path in paths) == EXPECTED_FILENAMES
    assert all(path.parent == output_directory.resolve() for path in paths)
    for written_table, expected_table in zip(written, expected, strict=True):
        pd.testing.assert_frame_equal(
            written_table, expected_table, check_exact=True
        )


def test_writer_preserves_unrelated_files_and_manages_only_three_csvs(
    tmp_path: Path,
) -> None:
    output_directory = tmp_path / "safe-output"
    output_directory.mkdir()
    sentinel = output_directory / "README.md"
    unrelated = output_directory / "notes.txt"
    sentinel.write_text("keep this documentation\n", encoding="utf-8")
    unrelated.write_text("keep this note\n", encoding="utf-8")
    before = {
        sentinel.name: sentinel.read_bytes(),
        unrelated.name: unrelated.read_bytes(),
    }

    paths = write_demo_data(output_directory)

    assert sentinel.read_bytes() == before[sentinel.name]
    assert unrelated.read_bytes() == before[unrelated.name]
    assert {path.name for path in output_directory.iterdir()} == {
        *EXPECTED_FILENAMES,
        sentinel.name,
        unrelated.name,
    }
    assert {path.name for path in paths} == set(EXPECTED_FILENAMES)


def test_cli_default_output_directory_is_repository_relative(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    external_working_directory = tmp_path / "outside-repository"
    external_working_directory.mkdir()
    assert not external_working_directory.is_relative_to(PROJECT_ROOT)
    calls: list[tuple[Path, int]] = []

    def record_write(
        output_dir: Path,
        seed: int = DEMO_SEED,
    ) -> tuple[Path, Path, Path]:
        resolved_output_dir = Path(output_dir).resolve()
        calls.append((resolved_output_dir, seed))
        return tuple(
            resolved_output_dir / filename for filename in EXPECTED_FILENAMES
        )

    monkeypatch.chdir(external_working_directory)
    monkeypatch.setattr(demo_data_cli, "write_demo_data", record_write)
    monkeypatch.setattr(sys, "argv", ["generate_demo_data.py"])

    demo_data_cli.main()

    expected_output_directory = (
        Path(demo_data_cli.__file__).resolve().parents[1] / "data" / "demo"
    )
    assert expected_output_directory == PROJECT_ROOT / "data" / "demo"
    assert calls == [(expected_output_directory, DEMO_SEED)]
    assert not calls[0][0].is_relative_to(external_working_directory)


def test_cli_explicit_output_directory_and_seed_override_defaults(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    external_working_directory = tmp_path / "outside-repository"
    external_working_directory.mkdir()
    custom_output_directory = tmp_path / "custom-output"
    calls: list[tuple[Path, int]] = []

    def record_write(
        output_dir: Path,
        seed: int = DEMO_SEED,
    ) -> tuple[Path, Path, Path]:
        resolved_output_dir = Path(output_dir).resolve()
        calls.append((resolved_output_dir, seed))
        return tuple(
            resolved_output_dir / filename for filename in EXPECTED_FILENAMES
        )

    monkeypatch.chdir(external_working_directory)
    monkeypatch.setattr(demo_data_cli, "write_demo_data", record_write)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "generate_demo_data.py",
            "--output-dir",
            str(custom_output_directory),
            "--seed",
            "12345",
        ],
    )

    demo_data_cli.main()

    assert calls == [(custom_output_directory.resolve(), 12345)]
    assert calls[0][0] != demo_data_cli.DEFAULT_OUTPUT_DIR
    assert not custom_output_directory.exists()
