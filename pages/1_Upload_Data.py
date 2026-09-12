"""Load and validate the three Phase 4 input tables."""

import streamlit as st

from plant_expression_explorer.dataset import (
    CANDIDATE_LABEL_KEY,
    CANDIDATE_REPORT_KEY,
    CANDIDATE_SOURCE_KEY,
    DE_RESULTS_UPLOAD_KEY,
    DEMO_SOURCE_LABEL,
    EXPRESSION_UPLOAD_KEY,
    METADATA_UPLOAD_KEY,
    CandidateResult,
    build_dataset_bundle,
    clear_candidate_feedback,
    clear_legacy_data_state,
    clear_uploader_state,
    get_current_dataset,
    load_demo_candidate,
    load_uploaded_candidate,
    record_candidate_failure,
    reset_data_state,
    set_current_dataset,
    table_preview,
    uploaded_source_label,
)
from plant_expression_explorer.validation import (
    Severity,
    ValidationIssue,
    ValidationReport,
)


def _activate_demo() -> None:
    candidate = load_demo_candidate()
    if candidate.status == "valid":
        bundle = build_dataset_bundle(
            candidate.tables,
            source="demo",
            source_label=DEMO_SOURCE_LABEL,
            report=candidate.report,
        )
        set_current_dataset(st.session_state, bundle)
        clear_uploader_state(st.session_state)
        return

    record_candidate_failure(
        st.session_state,
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
        report=candidate.report,
    )


def _activate_uploaded_candidate(
    candidate: CandidateResult,
    source_label: str,
) -> None:
    if candidate.status == "valid":
        bundle = build_dataset_bundle(
            candidate.tables,
            source="uploaded",
            source_label=source_label,
            report=candidate.report,
        )
        set_current_dataset(st.session_state, bundle)
    elif candidate.status == "invalid":
        record_candidate_failure(
            st.session_state,
            source="uploaded",
            source_label=source_label,
            report=candidate.report,
        )


def _issue_body(issue: ValidationIssue) -> str:
    locations = [issue.table]
    if issue.related_table is not None:
        locations.append(issue.related_table)

    evidence: list[str] = []
    if issue.column is not None:
        evidence.append(f"column: `{issue.column}`")
    if issue.count != 1:
        evidence.append(f"count: `{issue.count}`")
    if issue.row_positions:
        positions = ", ".join(str(value) for value in issue.row_positions)
        evidence.append(f"row positions: `{positions}`")
    if issue.example_values:
        examples = ", ".join(issue.example_values)
        evidence.append(f"examples: `{examples}`")

    body = (
        f"**{issue.severity.value.upper()} · `{issue.code.value}`**  \n"
        f"Location: {' ↔ '.join(locations)}  \n"
        f"{issue.message}"
    )
    if evidence:
        body += "  \nEvidence — " + "; ".join(evidence)
    return body


def render_validation_report(report: ValidationReport) -> None:
    """Render existing structured issues without adding validation rules."""

    sections = (
        ("Errors", Severity.ERROR, st.error),
        ("Warnings", Severity.WARNING, st.warning),
        ("Information", Severity.INFORMATION, st.info),
    )
    for title, severity, render_issue in sections:
        issues = tuple(
            issue for issue in report.issues if issue.severity is severity
        )
        st.markdown(f"#### {title} ({len(issues)})")
        if not issues:
            st.caption("None.")
            continue
        for issue in issues:
            render_issue(_issue_body(issue))


def _render_current_status() -> None:
    current = get_current_dataset(st.session_state)
    failed_report = st.session_state.get(CANDIDATE_REPORT_KEY)

    st.header("Current dataset status")
    if current is None:
        st.warning("No validated dataset is currently active.")
        st.write(
            "A validated dataset is not yet available for the downstream "
            "analysis pages."
        )
    else:
        st.success("A validated dataset is active.")
        st.write(f"**Source:** {current.source}")
        st.write(f"**Source label:** {current.source_label}")
        st.write(
            f"**Counts:** {current.gene_count:,} genes, "
            f"{current.sample_count:,} samples, and "
            f"{current.de_row_count:,} differential-expression rows."
        )
        st.write(
            "A validated dataset is available for the downstream analysis pages."
        )

    if isinstance(failed_report, ValidationReport) and failed_report.has_errors:
        if current is None:
            st.error(
                "The latest candidate failed validation, so no dataset is active."
            )
        else:
            st.warning(
                "The latest candidate failed validation. The active dataset remains "
                "the previous successfully validated dataset shown above; it does "
                "not represent the failed candidate."
            )


def _render_table_summary(
    title: str,
    table,
    *,
    source_label: str,
    count_label: str,
) -> None:
    st.subheader(title)
    st.caption(f"Source: {source_label}")
    st.write(
        f"Shape: {table.shape[0]:,} rows × {table.shape[1]:,} columns · "
        f"{count_label}"
    )
    st.table(table_preview(table))


clear_legacy_data_state(st.session_state)

st.title("📤 Upload Data")
st.write(
    "Load the three preprocessed transcriptomics tables required by the explorer. "
    "This workflow validates table structure and cross-file identifiers; it does "
    "not assess raw sequencing data or perform statistical inference."
)

status_container = st.container()

st.header("Synthetic demo data")
st.warning(
    "The bundled demo values are synthetic, the gene IDs are fictional, and the "
    "p-values are constructed. They do not support real tomato nitrate-response "
    "inference. DESeq2 and other fitted RNA-seq models were not run."
)
demo_requested = st.button(
    "Load synthetic demo data",
    key="pee_load_demo",
    on_click=_activate_demo,
)

st.header("Upload your own data")
st.write(
    "Uploaded tables are retained in the active Streamlit session for application "
    "use. The application does not intentionally write uploaded tables to project "
    "files."
)

expression_file = st.file_uploader(
    "Normalized expression matrix",
    type="csv",
    key=EXPRESSION_UPLOAD_KEY,
    help="Required: gene_id and at least two numeric sample columns.",
)
metadata_file = st.file_uploader(
    "Sample metadata",
    type="csv",
    key=METADATA_UPLOAD_KEY,
    help="Required: sample_id and condition.",
)
de_file = st.file_uploader(
    "Precomputed differential-expression results",
    type="csv",
    key=DE_RESULTS_UPLOAD_KEY,
    help="Required: gene_id, log2FoldChange, pvalue, and padj.",
)

uploaded_sources = (expression_file, metadata_file, de_file)
missing_uploads = [
    label
    for label, source in zip(
        (
            "expression matrix",
            "sample metadata",
            "differential-expression results",
        ),
        uploaded_sources,
        strict=True,
    )
    if source is None
]

if missing_uploads:
    st.info("Still required: " + ", ".join(missing_uploads) + ".")
    if not demo_requested:
        clear_candidate_feedback(st.session_state)
elif not demo_requested:
    upload_label = uploaded_source_label(*uploaded_sources)
    uploaded_candidate = load_uploaded_candidate(*uploaded_sources)
    _activate_uploaded_candidate(uploaded_candidate, upload_label)

with status_container:
    _render_current_status()

st.header("Candidate validation results")
candidate_report = st.session_state.get(CANDIDATE_REPORT_KEY)
candidate_source = st.session_state.get(CANDIDATE_SOURCE_KEY)
candidate_label = st.session_state.get(CANDIDATE_LABEL_KEY)
current = get_current_dataset(st.session_state)

if isinstance(candidate_report, ValidationReport) and candidate_report.has_errors:
    st.subheader("Failed candidate report")
    st.write(f"**Candidate source:** {candidate_source}")
    st.write(f"**Candidate label:** {candidate_label}")
    render_validation_report(candidate_report)
else:
    st.info("No failed candidate is currently recorded.")

if current is not None:
    st.subheader("Current active dataset report")
    st.write(f"**Current source:** {current.source_label}")
    render_validation_report(current.validation_report)
else:
    st.caption("There is no current active dataset report.")

st.header("Current valid table summaries and previews")
if current is None:
    st.info("Load a valid complete dataset to view table summaries and previews.")
else:
    _render_table_summary(
        "Expression matrix",
        current.expression,
        source_label=current.source_label,
        count_label=f"{current.gene_count:,} genes",
    )
    _render_table_summary(
        "Sample metadata",
        current.metadata,
        source_label=current.source_label,
        count_label=f"{current.sample_count:,} samples",
    )
    _render_table_summary(
        "Differential-expression results",
        current.de_results,
        source_label=current.source_label,
        count_label=f"{current.de_row_count:,} rows",
    )

st.header("Reset data")
st.button(
    "Reset Data",
    key="pee_reset_data",
    on_click=reset_data_state,
    args=(st.session_state,),
)
st.caption(
    "Reset removes application-owned session data and uploader values. It does "
    "not delete user files or bundled demo files from disk."
)

st.header("Upload Data page scope")
st.info(
    "This page loads and validates data only. Sample Quality Control, PCA, and "
    "Sample Correlation are available on their own pages for the active "
    "dataset. Descriptive exploratory threshold classification of supplied, "
    "precomputed differential-expression results is also available, together "
    "with exact, descriptive single-gene expression lookup. Dedicated CSV "
    "downloads of current descriptive result tables are available on the "
    "analysis pages. Differential-expression modelling and volcano plots are "
    "not implemented."
)
