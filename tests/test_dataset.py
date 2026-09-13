"""Tests for Phase 4 dataset loading and session-state contracts."""

from __future__ import annotations

import shutil
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

import plant_expression_explorer.dataset as dataset_module
from plant_expression_explorer.consistency import validate_input_tables
from plant_expression_explorer.dataset import (
    ACTIVE_GROUP_COLUMN_KEY,
    APPLICATION_DATA_KEYS,
    CANDIDATE_LABEL_KEY,
    CANDIDATE_REPORT_KEY,
    CANDIDATE_SOURCE_KEY,
    DATASET_CONTEXT_KEYS,
    CURRENT_DATASET_KEY,
    DE_RESULTS_UPLOAD_KEY,
    DEMO_DIRECTORY,
    DEMO_SOURCE_LABEL,
    EXPRESSION_UPLOAD_KEY,
    LEGACY_DATA_KEYS,
    METADATA_UPLOAD_KEY,
    CandidateResult,
    DatasetBundle,
    build_dataset_bundle,
    clear_dataset_context_state,
    clear_legacy_data_state,
    clear_uploader_state,
    ensure_valid_group_column_state,
    get_current_dataset,
    load_demo_candidate,
    load_uploaded_candidate,
    record_candidate_failure,
    reset_data_state,
    set_current_dataset,
    table_preview,
    uploaded_source_label,
)
from plant_expression_explorer.provenance import DatasetProvenance
from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
)


class NamedStringIO(StringIO):
    """Small seekable uploaded-file stand-in with a user-facing name."""

    def __init__(self, content: str, name: str) -> None:
        super().__init__(content)
        self.name = name


def _valid_uploaded_sources() -> tuple[NamedStringIO, NamedStringIO, NamedStringIO]:
    expression = NamedStringIO(
        "gene_id,s1,s2,s3,s4\n"
        "g1,1.0,1.1,2.0,2.1\n"
        "g2,3.0,3.1,4.0,4.1\n",
        "expression.csv",
    )
    metadata = NamedStringIO(
        "sample_id,condition\n"
        "s1,control\n"
        "s2,control\n"
        "s3,treated\n"
        "s4,treated\n",
        "metadata.csv",
    )
    de_results = NamedStringIO(
        "gene_id,log2FoldChange,pvalue,padj\n"
        "g1,1.0,0.01,0.02\n"
        "g2,-1.0,0.02,0.03\n",
        "deg.csv",
    )
    return expression, metadata, de_results


def _warning_uploaded_sources() -> tuple[NamedStringIO, NamedStringIO, NamedStringIO]:
    sources = _valid_uploaded_sources()
    sources[0].seek(0)
    sources = (
        NamedStringIO(
            "gene_id,s1,s2,s3,s4\n"
            "g1,-1.0,1.1,2.0,2.1\n"
            "g2,3.0,3.1,4.0,4.1\n",
            "expression-warning.csv",
        ),
        sources[1],
        sources[2],
    )
    return sources


def _information_uploaded_sources(
) -> tuple[NamedStringIO, NamedStringIO, NamedStringIO]:
    expression, _, de_results = _valid_uploaded_sources()
    metadata = NamedStringIO(
        "sample_id,condition\n"
        "s2,control\n"
        "s1,control\n"
        "s3,treated\n"
        "s4,treated\n",
        "metadata-reordered.csv",
    )
    return expression, metadata, de_results


def _invalid_uploaded_sources() -> tuple[NamedStringIO, NamedStringIO, NamedStringIO]:
    expression, metadata, de_results = _valid_uploaded_sources()
    expression = NamedStringIO(
        "wrong_id,s1,s2,s3,s4\n"
        "g1,1.0,1.1,2.0,2.1\n",
        "invalid-expression.csv",
    )
    return expression, metadata, de_results


def _bundle_from_candidate(
    candidate: CandidateResult,
    *,
    source: str,
    source_label: str,
) -> DatasetBundle:
    assert candidate.status == "valid"
    assert candidate.tables is not None
    return build_dataset_bundle(
        candidate.tables,
        source=source,
        source_label=source_label,
        report=candidate.report,
    )


def _copy_demo_directory(destination: Path) -> Path:
    copied = destination / "demo"
    shutil.copytree(DEMO_DIRECTORY, copied)
    return copied


def test_default_demo_loading_returns_clean_valid_candidate() -> None:
    candidate = load_demo_candidate()

    assert candidate.status == "valid"
    assert candidate.tables is not None
    assert [table.shape for table in candidate.tables] == [
        (120, 7),
        (6, 2),
        (120, 4),
    ]
    assert candidate.report.issues == ()


def test_demo_and_upload_routes_use_shared_reader_and_aggregate_validator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    reader_calls = 0
    validator_calls = 0
    original_reader = dataset_module.read_csv
    original_validator = dataset_module.validate_input_tables

    def read_spy(source):
        nonlocal reader_calls
        reader_calls += 1
        return original_reader(source)

    def validation_spy(expression, metadata, de_results):
        nonlocal validator_calls
        validator_calls += 1
        return original_validator(expression, metadata, de_results)

    monkeypatch.setattr(dataset_module, "read_csv", read_spy)
    monkeypatch.setattr(dataset_module, "validate_input_tables", validation_spy)

    assert load_demo_candidate().status == "valid"
    assert load_uploaded_candidate(*_valid_uploaded_sources()).status == "valid"
    assert reader_calls == 6
    assert validator_calls == 2


def test_default_demo_path_is_independent_of_current_working_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.chdir(tmp_path)

    candidate = load_demo_candidate()

    assert candidate.status == "valid"
    assert candidate.report.issues == ()


@pytest.mark.parametrize(
    ("filename", "table_name"),
    [
        ("expression_matrix.csv", "Expression matrix"),
        ("metadata.csv", "Sample metadata"),
        ("deg_results.csv", "Differential-expression results"),
    ],
)
def test_each_missing_demo_file_produces_a_controlled_error(
    tmp_path: Path,
    filename: str,
    table_name: str,
) -> None:
    demo_dir = _copy_demo_directory(tmp_path)
    (demo_dir / filename).unlink()

    candidate = load_demo_candidate(demo_dir)

    assert candidate.status == "invalid"
    assert candidate.tables is None
    assert candidate.report.has_errors
    assert len(candidate.report.errors) == 1
    issue = candidate.report.errors[0]
    assert issue.code is IssueCode.DEMO_FILE_MISSING
    assert issue.severity is Severity.ERROR
    assert issue.table == table_name
    assert filename in issue.message


def test_malformed_demo_csv_produces_a_controlled_read_error(
    tmp_path: Path,
) -> None:
    demo_dir = _copy_demo_directory(tmp_path)
    (demo_dir / "expression_matrix.csv").write_text(
        'gene_id,s1\n"unterminated,1\n',
        encoding="utf-8",
    )

    candidate = load_demo_candidate(demo_dir)

    assert candidate.status == "invalid"
    assert candidate.tables is None
    assert candidate.report.has_errors
    issue = candidate.report.errors[0]
    assert issue.code is IssueCode.CSV_READ_ERROR
    assert issue.severity is Severity.ERROR
    assert issue.table == "Expression matrix"
    assert "Traceback" not in issue.message


def test_valid_uploaded_sources_return_a_valid_candidate() -> None:
    candidate = load_uploaded_candidate(*_valid_uploaded_sources())

    assert candidate.status == "valid"
    assert candidate.tables is not None
    assert not candidate.report.has_errors
    assert candidate.report.issues == ()


@pytest.mark.parametrize("missing_index", [0, 1])
def test_partial_uploaded_sources_return_incomplete(
    missing_index: int,
) -> None:
    sources: list[object | None] = list(_valid_uploaded_sources())
    sources[missing_index] = None

    candidate = load_uploaded_candidate(*sources)

    assert candidate.status == "incomplete"
    assert candidate.tables is None
    assert candidate.report.issues == ()


def test_expression_and_metadata_without_de_results_are_a_valid_candidate() -> None:
    expression, metadata, _ = _valid_uploaded_sources()

    candidate = load_uploaded_candidate(expression, metadata, None)

    assert candidate.status == "valid"
    assert candidate.tables is not None
    assert candidate.tables[2] is None
    assert not candidate.report.has_errors
    assert [issue.code for issue in candidate.report.information] == [
        IssueCode.DE_RESULTS_NOT_SUPPLIED
    ]


def test_optional_de_candidate_reads_only_supplied_sources(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    expression, metadata, _ = _valid_uploaded_sources()
    reader_calls = 0
    original_reader = dataset_module.read_csv

    def read_spy(source):
        nonlocal reader_calls
        reader_calls += 1
        return original_reader(source)

    monkeypatch.setattr(dataset_module, "read_csv", read_spy)

    candidate = load_uploaded_candidate(expression, metadata, None)

    assert candidate.status == "valid"
    assert reader_calls == 2


def test_partial_upload_does_not_call_reader_or_validator(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def unexpected_call(*args: object, **kwargs: object) -> None:
        raise AssertionError("Partial upload must not be parsed or validated.")

    monkeypatch.setattr(dataset_module, "read_csv", unexpected_call)
    monkeypatch.setattr(dataset_module, "validate_input_tables", unexpected_call)
    _, metadata, de_results = _valid_uploaded_sources()

    candidate = load_uploaded_candidate(None, metadata, de_results)

    assert candidate.status == "incomplete"


def test_parser_error_is_a_structured_report_without_traceback() -> None:
    _, metadata, de_results = _valid_uploaded_sources()
    malformed = NamedStringIO(
        'gene_id,s1\n"unterminated,1\n',
        "broken-expression.csv",
    )

    candidate = load_uploaded_candidate(malformed, metadata, de_results)

    assert candidate.status == "invalid"
    assert candidate.tables is None
    issue = candidate.report.errors[0]
    assert issue.code is IssueCode.CSV_READ_ERROR
    assert issue.severity is Severity.ERROR
    assert issue.table == "Expression matrix"
    assert "broken-expression.csv" in issue.message
    assert "Traceback" not in issue.message


def test_duplicate_header_preserves_existing_structured_issue_code() -> None:
    _, metadata, de_results = _valid_uploaded_sources()
    duplicate = NamedStringIO(
        "gene_id,s1,s1\ng1,1.0,2.0\n",
        "duplicate-header.csv",
    )

    candidate = load_uploaded_candidate(duplicate, metadata, de_results)

    assert candidate.status == "invalid"
    assert candidate.report.errors[0].code is IssueCode.DUPLICATE_COLUMN_NAME
    assert candidate.report.errors[0].table == "Expression matrix"


def test_candidate_result_invariants_are_enforced() -> None:
    error_report = ValidationReport(
        (
            ValidationIssue(
                code=IssueCode.CSV_READ_ERROR,
                severity=Severity.ERROR,
                table="Expression matrix",
                message="Could not read CSV.",
            ),
        )
    )

    with pytest.raises(ValueError, match="incomplete"):
        CandidateResult("incomplete", None, error_report)
    with pytest.raises(ValueError, match="invalid"):
        CandidateResult("invalid", None, ValidationReport())
    with pytest.raises(ValueError, match="valid"):
        CandidateResult("valid", None, ValidationReport())


def test_error_candidate_cannot_build_or_activate_a_bundle() -> None:
    candidate = load_uploaded_candidate(*_invalid_uploaded_sources())
    assert candidate.status == "invalid"
    assert candidate.tables is None
    state: dict[str, object] = {}
    demo_candidate = load_demo_candidate()
    assert demo_candidate.tables is not None

    with pytest.raises(ValueError, match="Errors"):
        build_dataset_bundle(
            demo_candidate.tables,
            source="uploaded",
            source_label="Invalid upload",
            report=candidate.report,
        )

    record_candidate_failure(
        state,
        source="uploaded",
        source_label="Invalid upload",
        report=candidate.report,
    )
    assert get_current_dataset(state) is None


def test_candidate_failure_requires_an_error_report() -> None:
    state: dict[str, object] = {}

    with pytest.raises(ValueError, match="at least one Error"):
        record_candidate_failure(
            state,
            source="uploaded",
            source_label="Non-failing candidate",
            report=ValidationReport(),
        )

    assert state == {}


@pytest.mark.parametrize(
    ("sources_factory", "expected_severity"),
    [
        (_warning_uploaded_sources, Severity.WARNING),
        (_information_uploaded_sources, Severity.INFORMATION),
    ],
)
def test_warning_and_information_candidates_can_activate(
    sources_factory,
    expected_severity: Severity,
) -> None:
    candidate = load_uploaded_candidate(*sources_factory())
    bundle = _bundle_from_candidate(
        candidate,
        source="uploaded",
        source_label="User-uploaded CSV tables",
    )
    state: dict[str, object] = {}

    set_current_dataset(state, bundle)

    assert get_current_dataset(state) is bundle
    assert any(
        issue.severity is expected_severity
        for issue in bundle.validation_report.issues
    )


def test_demo_and_uploaded_source_labels_are_explicit() -> None:
    demo_bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    sources = _valid_uploaded_sources()
    upload_label = uploaded_source_label(*sources)
    uploaded_bundle = _bundle_from_candidate(
        load_uploaded_candidate(*sources),
        source="uploaded",
        source_label=upload_label,
    )

    assert demo_bundle.source == "demo"
    assert "synthetic" in demo_bundle.source_label.lower()
    assert uploaded_bundle.source == "uploaded"
    assert "expression.csv" in uploaded_bundle.source_label
    assert "metadata.csv" in uploaded_bundle.source_label
    assert "deg.csv" in uploaded_bundle.source_label


def test_uploaded_source_label_omits_optional_de_when_not_supplied() -> None:
    expression, metadata, _ = _valid_uploaded_sources()

    label = uploaded_source_label(expression, metadata, None)

    assert "expression.csv" in label
    assert "metadata.csv" in label
    assert "deg.csv" not in label


def test_demo_to_uploaded_replacement_is_a_single_bundle_assignment() -> None:
    demo_bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    sources = _valid_uploaded_sources()
    uploaded_bundle = _bundle_from_candidate(
        load_uploaded_candidate(*sources),
        source="uploaded",
        source_label=uploaded_source_label(*sources),
    )
    state: dict[str, object] = {}

    set_current_dataset(state, demo_bundle)
    set_current_dataset(state, uploaded_bundle)

    assert get_current_dataset(state) is uploaded_bundle
    assert "expression" not in state
    assert "metadata" not in state
    assert "de_results" not in state


def test_uploaded_to_demo_replacement_and_uploader_clear_are_atomic_in_state() -> None:
    sources = _valid_uploaded_sources()
    uploaded_bundle = _bundle_from_candidate(
        load_uploaded_candidate(*sources),
        source="uploaded",
        source_label=uploaded_source_label(*sources),
    )
    demo_bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    state: dict[str, object] = {
        EXPRESSION_UPLOAD_KEY: sources[0],
        METADATA_UPLOAD_KEY: sources[1],
        DE_RESULTS_UPLOAD_KEY: sources[2],
    }

    set_current_dataset(state, uploaded_bundle)
    set_current_dataset(state, demo_bundle)
    clear_uploader_state(state)

    assert get_current_dataset(state) is demo_bundle
    assert all(key not in state for key in (
        EXPRESSION_UPLOAD_KEY,
        METADATA_UPLOAD_KEY,
        DE_RESULTS_UPLOAD_KEY,
    ))


def test_valid_demo_remains_active_after_invalid_upload() -> None:
    demo_bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    invalid = load_uploaded_candidate(*_invalid_uploaded_sources())
    state: dict[str, object] = {}

    set_current_dataset(state, demo_bundle)
    record_candidate_failure(
        state,
        source="uploaded",
        source_label="Invalid uploaded candidate",
        report=invalid.report,
    )

    assert get_current_dataset(state) is demo_bundle
    assert state[CANDIDATE_SOURCE_KEY] == "uploaded"
    assert state[CANDIDATE_REPORT_KEY] is invalid.report


def test_valid_upload_remains_active_after_invalid_demo(
    tmp_path: Path,
) -> None:
    sources = _valid_uploaded_sources()
    uploaded_bundle = _bundle_from_candidate(
        load_uploaded_candidate(*sources),
        source="uploaded",
        source_label=uploaded_source_label(*sources),
    )
    invalid_demo = load_demo_candidate(tmp_path / "missing-demo")
    state: dict[str, object] = {}

    set_current_dataset(state, uploaded_bundle)
    record_candidate_failure(
        state,
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
        report=invalid_demo.report,
    )

    assert get_current_dataset(state) is uploaded_bundle
    assert state[CANDIDATE_SOURCE_KEY] == "demo"
    assert state[CANDIDATE_REPORT_KEY] is invalid_demo.report


def test_candidate_failure_replaces_instead_of_accumulating(
    tmp_path: Path,
) -> None:
    first = load_uploaded_candidate(*_invalid_uploaded_sources()).report
    second = load_demo_candidate(tmp_path / "missing-demo").report
    state: dict[str, object] = {}

    record_candidate_failure(
        state,
        source="uploaded",
        source_label="First failure",
        report=first,
    )
    record_candidate_failure(
        state,
        source="demo",
        source_label="Second failure",
        report=second,
    )

    assert state[CANDIDATE_REPORT_KEY] is second
    assert state[CANDIDATE_SOURCE_KEY] == "demo"
    assert state[CANDIDATE_LABEL_KEY] == "Second failure"


def test_successful_activation_clears_candidate_feedback() -> None:
    invalid = load_uploaded_candidate(*_invalid_uploaded_sources())
    demo_bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    state: dict[str, object] = {}
    record_candidate_failure(
        state,
        source="uploaded",
        source_label="Failed upload",
        report=invalid.report,
    )

    set_current_dataset(state, demo_bundle)

    assert get_current_dataset(state) is demo_bundle
    assert all(
        key not in state
        for key in (
            CANDIDATE_REPORT_KEY,
            CANDIDATE_SOURCE_KEY,
            CANDIDATE_LABEL_KEY,
        )
    )


def test_set_current_dataset_clears_the_shared_group_by_choice() -> None:
    demo_bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    state: dict[str, object] = {ACTIVE_GROUP_COLUMN_KEY: "genotype"}

    set_current_dataset(state, demo_bundle)

    assert ACTIVE_GROUP_COLUMN_KEY not in state


def test_ensure_valid_group_column_state_keeps_a_still_applicable_choice() -> None:
    state: dict[str, object] = {ACTIVE_GROUP_COLUMN_KEY: "genotype"}

    ensure_valid_group_column_state(state, ("genotype", "batch"))

    assert state[ACTIVE_GROUP_COLUMN_KEY] == "genotype"


def test_ensure_valid_group_column_state_defaults_an_unset_choice_to_condition() -> (
    None
):
    state: dict[str, object] = {}

    ensure_valid_group_column_state(state, ("genotype",))

    assert state[ACTIVE_GROUP_COLUMN_KEY] == "condition"


def test_ensure_valid_group_column_state_resets_a_stale_choice() -> None:
    state: dict[str, object] = {ACTIVE_GROUP_COLUMN_KEY: "genotype"}

    ensure_valid_group_column_state(state, ("batch",))

    assert state[ACTIVE_GROUP_COLUMN_KEY] == "condition"


def test_ensure_valid_group_column_state_resets_when_no_additional_columns() -> None:
    state: dict[str, object] = {ACTIVE_GROUP_COLUMN_KEY: "genotype"}

    ensure_valid_group_column_state(state, ())

    assert state[ACTIVE_GROUP_COLUMN_KEY] == "condition"


def test_reset_deletes_only_explicit_application_keys_and_preserves_unrelated() -> None:
    state = {key: object() for key in APPLICATION_DATA_KEYS}
    state["unrelated_state"] = "keep"
    state["pee_unrelated_state"] = "also keep"

    reset_data_state(state)

    assert all(key not in state for key in APPLICATION_DATA_KEYS)
    assert state == {
        "unrelated_state": "keep",
        "pee_unrelated_state": "also keep",
    }


def test_dataset_context_clear_removes_only_documented_context_keys() -> None:
    state = {key: object() for key in DATASET_CONTEXT_KEYS}
    state[CURRENT_DATASET_KEY] = "keep current"
    state[EXPRESSION_UPLOAD_KEY] = "keep upload"
    state["other"] = "keep"

    clear_dataset_context_state(state)

    assert all(key not in state for key in DATASET_CONTEXT_KEYS)
    assert state == {
        CURRENT_DATASET_KEY: "keep current",
        EXPRESSION_UPLOAD_KEY: "keep upload",
        "other": "keep",
    }


def test_reset_is_idempotent_and_leaves_no_application_data_state() -> None:
    state = {key: object() for key in APPLICATION_DATA_KEYS}

    reset_data_state(state)
    reset_data_state(state)

    assert state == {}


def test_legacy_cleanup_removes_only_documented_legacy_keys() -> None:
    state = {
        "expression": object(),
        "metadata": object(),
        "de_results": object(),
        CURRENT_DATASET_KEY: "keep current",
        "other": "keep",
    }

    clear_legacy_data_state(state)

    assert all(key not in state for key in LEGACY_DATA_KEYS)
    assert state[CURRENT_DATASET_KEY] == "keep current"
    assert state["other"] == "keep"


def test_same_uploaded_sources_can_be_read_repeatedly() -> None:
    sources = _valid_uploaded_sources()

    first = load_uploaded_candidate(*sources)
    second = load_uploaded_candidate(*sources)

    assert first.status == second.status == "valid"
    assert first.tables is not None
    assert second.tables is not None
    for first_table, second_table in zip(first.tables, second.tables, strict=True):
        assert not second_table.empty
        assert_frame_equal(first_table, second_table)


def test_uploaded_source_cursor_positions_are_restored() -> None:
    sources = _valid_uploaded_sources()
    original_positions = (4, 7, 2)
    for source, position in zip(sources, original_positions, strict=True):
        source.seek(position)

    candidate = load_uploaded_candidate(*sources)

    assert candidate.status == "valid"
    assert tuple(source.tell() for source in sources) == original_positions


def test_bundle_state_and_preview_helpers_do_not_mutate_dataframes() -> None:
    candidate = load_uploaded_candidate(*_valid_uploaded_sources())
    assert candidate.tables is not None
    originals = tuple(table.copy(deep=True) for table in candidate.tables)
    bundle = build_dataset_bundle(
        candidate.tables,
        source="uploaded",
        source_label="User-uploaded CSV tables",
        report=candidate.report,
    )
    state: dict[str, object] = {}

    set_current_dataset(state, bundle)
    preview = table_preview(bundle.expression, 1)
    preview.iloc[0, 0] = "changed-in-preview"
    get_current_dataset(state)

    for actual, expected in zip(candidate.tables, originals, strict=True):
        assert_frame_equal(actual, expected)


def test_bundle_accepts_optional_de_and_preserves_provenance_identity() -> None:
    expression, metadata, _ = _valid_uploaded_sources()
    candidate = load_uploaded_candidate(expression, metadata, None)
    assert candidate.tables is not None
    provenance = DatasetProvenance(
        dataset_title="  exact title  ",
        expression_scale_description="VST",
        notes="备注",
    )

    bundle = build_dataset_bundle(
        candidate.tables,
        source="uploaded",
        source_label="Expression and metadata only",
        report=candidate.report,
        provenance=provenance,
    )

    assert bundle.de_results is None
    assert bundle.has_de_results is False
    assert bundle.de_row_count == 0
    assert bundle.provenance is provenance
    assert bundle.provenance.dataset_title == "  exact title  "


def test_bundle_rejects_invalid_provenance_without_coercion() -> None:
    candidate = load_demo_candidate()
    assert candidate.tables is not None

    with pytest.raises(TypeError, match="DatasetProvenance"):
        build_dataset_bundle(
            candidate.tables,
            source="demo",
            source_label=DEMO_SOURCE_LABEL,
            report=candidate.report,
            provenance={"organism": "tomato"},  # type: ignore[arg-type]
        )


def test_preview_preserves_row_and_column_order_and_returns_a_copy() -> None:
    table = pd.DataFrame(
        {"second": [2, 1, 3], "first": ["b", "a", "c"]},
        index=[20, 10, 30],
    )

    preview = table_preview(table, 2)

    assert preview.columns.tolist() == ["second", "first"]
    assert preview.index.tolist() == [20, 10]
    assert preview["second"].tolist() == [2, 1]
    assert preview is not table


@pytest.mark.parametrize("rows", [0, -1, -10])
def test_preview_zero_or_negative_rows_returns_empty_copy(rows: int) -> None:
    table = pd.DataFrame({"a": [1, 2], "b": [3, 4]})

    preview = table_preview(table, rows)

    assert preview.empty
    assert preview.columns.tolist() == table.columns.tolist()
    assert preview is not table


def test_preview_more_rows_than_available_returns_all_rows() -> None:
    table = pd.DataFrame({"a": [1, 2]})

    preview = table_preview(table, 10)

    assert_frame_equal(preview, table)
    assert preview is not table


def test_rerun_style_repeated_state_calls_are_idempotent() -> None:
    bundle = _bundle_from_candidate(
        load_demo_candidate(),
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
    )
    state: dict[str, object] = {}

    set_current_dataset(state, bundle)
    set_current_dataset(state, bundle)

    assert state == {CURRENT_DATASET_KEY: bundle}


def test_loading_issue_codes_are_errors_and_not_emitted_by_phase_2_validation(
    tmp_path: Path,
) -> None:
    missing = load_demo_candidate(tmp_path / "missing-demo")
    malformed, metadata, de_results = _valid_uploaded_sources()
    malformed = NamedStringIO(
        'gene_id,s1\n"unterminated,1\n',
        "malformed.csv",
    )
    unreadable = load_uploaded_candidate(malformed, metadata, de_results)

    assert {
        issue.code: issue.severity
        for issue in (*missing.report.errors, *unreadable.report.errors)
    } == {
        IssueCode.DEMO_FILE_MISSING: Severity.ERROR,
        IssueCode.CSV_READ_ERROR: Severity.ERROR,
    }

    expression = pd.DataFrame({"gene_id": ["g1"], "s1": [1.0]})
    metadata_table = pd.DataFrame(
        {"sample_id": ["other"], "condition": ["control"]}
    )
    de_table = pd.DataFrame(
        {
            "gene_id": ["g2"],
            "log2FoldChange": ["bad"],
            "pvalue": [2.0],
            "padj": [None],
        }
    )
    report = validate_input_tables(expression, metadata_table, de_table)
    assert report.has_errors
    assert IssueCode.CSV_READ_ERROR not in {issue.code for issue in report.issues}
    assert IssueCode.DEMO_FILE_MISSING not in {
        issue.code for issue in report.issues
    }


def test_existing_phase_2_issue_code_and_severity_values_are_unchanged() -> None:
    expected_phase_2_codes = (
        "NOT_A_TABLE",
        "EMPTY_TABLE",
        "MISSING_REQUIRED_COLUMN",
        "EMPTY_COLUMN_NAME",
        "DUPLICATE_COLUMN_NAME",
        "NO_SAMPLE_COLUMNS",
        "LOW_SAMPLE_COUNT",
        "MISSING_IDENTIFIER",
        "DUPLICATE_IDENTIFIER",
        "WHITESPACE_IN_IDENTIFIER",
        "WHITESPACE_IN_SAMPLE_NAME",
        "WHITESPACE_IN_CATEGORY_VALUE",
        "MISSING_REQUIRED_VALUE",
        "SINGLE_CONDITION",
        "CONDITION_WITH_SINGLE_SAMPLE",
        "NON_NUMERIC_VALUE",
        "COERCIBLE_NUMERIC_STRING",
        "MISSING_EXPRESSION_VALUE",
        "NON_FINITE_VALUE",
        "NEGATIVE_EXPRESSION_VALUE",
        "MISSING_DE_STATISTIC",
        "NO_USABLE_VALUES",
        "VALUE_OUT_OF_RANGE",
        "ZERO_ADJUSTED_P_VALUE",
        "MISSING_METADATA_SAMPLE",
        "UNEXPECTED_METADATA_SAMPLE",
        "NO_SAMPLE_OVERLAP",
        "SAMPLE_ORDER_DIFFERS",
        "DE_GENE_NOT_IN_EXPRESSION",
        "EXPRESSION_GENE_NOT_IN_DE",
        "NO_GENE_OVERLAP",
    )
    non_phase_2_codes = {
        IssueCode.CSV_READ_ERROR,
        IssueCode.POSSIBLE_DELIMITER_MISMATCH,
        IssueCode.DEMO_FILE_MISSING,
        IssueCode.DE_RESULTS_NOT_SUPPLIED,
    }

    assert tuple(
        code.value for code in IssueCode if code not in non_phase_2_codes
    ) == expected_phase_2_codes
    assert IssueCode.DE_RESULTS_NOT_SUPPLIED.value == "DE_RESULTS_NOT_SUPPLIED"
    assert tuple(severity.value for severity in Severity) == (
        "error",
        "warning",
        "information",
    )
