"""Load, validate, and activate preprocessed transcriptomics tables."""

import streamlit as st

from plant_expression_explorer.dataset import (
    CANDIDATE_LABEL_KEY,
    CANDIDATE_REPORT_KEY,
    CANDIDATE_SOURCE_KEY,
    CONTEXT_NOTES_KEY,
    DATASET_TITLE_KEY,
    DE_CONTRAST_KEY,
    DE_RESULTS_UPLOAD_KEY,
    DEMO_SOURCE_LABEL,
    EXPRESSION_SCALE_KEY,
    EXPRESSION_UPLOAD_KEY,
    FEATURE_LEVEL_KEY,
    METADATA_UPLOAD_KEY,
    ORGANISM_KEY,
    REFERENCE_ANNOTATION_KEY,
    UPSTREAM_NORMALIZATION_KEY,
    CandidateResult,
    DatasetBundle,
    build_dataset_bundle,
    clear_candidate_feedback,
    clear_dataset_context_state,
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
from plant_expression_explorer.exports import CSV_MEDIA_TYPE
from plant_expression_explorer.provenance import (
    DEMO_PROVENANCE,
    DatasetProvenance,
    provenance_display_rows,
)
from plant_expression_explorer.templates import (
    build_example_de_results_template,
    build_example_expression_template,
    build_example_metadata_template,
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
            provenance=DEMO_PROVENANCE,
        )
        set_current_dataset(st.session_state, bundle)
        clear_uploader_state(st.session_state)
        clear_dataset_context_state(st.session_state)
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
    provenance: DatasetProvenance | None,
) -> None:
    if candidate.status == "valid":
        bundle = build_dataset_bundle(
            candidate.tables,
            source="uploaded",
            source_label=source_label,
            report=candidate.report,
            provenance=provenance,
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


def _render_dataset_context(provenance: DatasetProvenance | None) -> None:
    with st.expander("Dataset context (descriptive only)"):
        st.caption(
            "Context is displayed verbatim and is not scientifically verified, "
            "parsed, or used in any calculation."
        )
        for label, value in provenance_display_rows(provenance):
            st.caption(label)
            st.code(value, language=None)


def _render_overview_card(current: DatasetBundle) -> None:
    """Show an at-a-glance summary of the active dataset's shape."""

    overview_columns = st.columns(3)
    overview_columns[0].metric("Genes", f"{current.gene_count:,}")
    overview_columns[1].metric("Samples", f"{current.sample_count:,}")
    overview_columns[2].metric(
        "DE rows",
        f"{current.de_row_count:,}" if current.has_de_results else "Not supplied",
        help="Differential-expression results rows; DE is optional.",
    )
    condition_counts = (
        current.metadata["condition"]
        .value_counts(sort=False)
        .rename_axis("condition")
        .reset_index(name="sample_count")
    )
    st.caption("Samples per condition")
    st.bar_chart(
        condition_counts,
        x="condition",
        y="sample_count",
        x_label="Condition",
        y_label="Samples",
        sort=False,
        height=180,
    )


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
        _render_overview_card(current)
        _render_dataset_context(current.provenance)
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
    "Load a preprocessed expression matrix and matching sample metadata. A "
    "precomputed differential-expression table is optional and enables its own "
    "exploration page. This workflow validates table structure and cross-file "
    "identifiers; it does not assess raw sequencing data or perform statistical "
    "inference."
)

status_container = st.container()

demo_tab, upload_tab = st.tabs(["Use demo data", "Upload my own data"])

with demo_tab:
    st.warning(
        "The bundled demo values are synthetic, the gene IDs are fictional, and the "
        "p-values are constructed. They do not support real tomato nitrate-response "
        "inference. DESeq2 and other fitted RNA-seq models were not run. In this "
        "demo, larger preset effect sizes are constructed to produce smaller p-values "
        "and replicate noise is balanced within each condition; neither pattern is "
        "guaranteed in real experiments."
    )
    demo_requested = st.button(
        "Load synthetic demo data",
        key="pee_load_demo",
        on_click=_activate_demo,
    )

with upload_tab:
    st.write(
        "Uploaded tables are retained in the active Streamlit session for "
        "application use. The application does not intentionally write "
        "uploaded tables to project files."
    )

    with st.expander("Need the exact column format? Download example CSV templates"):
        st.caption(
            "Each template has fabricated placeholder values only, to illustrate "
            "the required column names and layout. Replace every value before "
            "uploading; these templates are not real biological data and are "
            "never loaded as a dataset."
        )
        template_columns = st.columns(3)
        template_columns[0].download_button(
            "Expression matrix template",
            data=build_example_expression_template().to_csv(index=False),
            file_name="expression_matrix_template.csv",
            mime=CSV_MEDIA_TYPE,
            key="pee_download_expression_template",
        )
        template_columns[1].download_button(
            "Sample metadata template",
            data=build_example_metadata_template().to_csv(index=False),
            file_name="sample_metadata_template.csv",
            mime=CSV_MEDIA_TYPE,
            key="pee_download_metadata_template",
        )
        template_columns[2].download_button(
            "DE results template",
            data=build_example_de_results_template().to_csv(index=False),
            file_name="de_results_template.csv",
            mime=CSV_MEDIA_TYPE,
            key="pee_download_de_template",
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
        help=(
            "Required: sample_id and condition. Optional design columns such as "
            "genotype, tissue, developmental_stage, timepoint, dose, batch, block, "
            "and biological_replicate are preserved but not interpreted."
        ),
    )
    de_file = st.file_uploader(
        "Precomputed differential-expression results",
        type="csv",
        key=DE_RESULTS_UPLOAD_KEY,
        help=(
            "Optional: gene_id, log2FoldChange, pvalue, and padj. Supplying this "
            "file enables the Differential Expression page."
        ),
    )
    st.caption(
        "Each file may be up to 500 MB. A typical plant genome (around "
        "30,000-70,000 genes) with up to a few hundred samples is well within "
        "this limit; whole-genome matrices with very many samples may approach "
        "it."
    )

    with st.expander("Dataset context (optional)"):
        st.caption(
            "Non-empty entries are carried and displayed exactly as typed. "
            "Untouched empty fields are recorded as not supplied. Context is "
            "not scientifically verified, parsed, or used in any calculation."
        )
        dataset_title = st.text_input("Dataset title", key=DATASET_TITLE_KEY)
        organism = st.text_input("Organism / taxon", key=ORGANISM_KEY)
        expression_scale = st.text_input(
            "Expression scale or preprocessing description",
            key=EXPRESSION_SCALE_KEY,
            help=(
                "For example: VST, rlog, log2(TPM + 1), logCPM, or another exact "
                "description from the upstream workflow. The application does "
                "not infer or validate this text."
            ),
        )
        upstream_normalization = st.text_input(
            "Upstream normalization method",
            key=UPSTREAM_NORMALIZATION_KEY,
        )
        reference_annotation = st.text_input(
            "Reference genome / annotation release",
            key=REFERENCE_ANNOTATION_KEY,
        )
        feature_level = st.text_input(
            "Feature level",
            key=FEATURE_LEVEL_KEY,
            help="For example: gene, transcript, or another supplied feature unit.",
        )
        de_contrast = st.text_input(
            "Differential-expression contrast description",
            key=DE_CONTRAST_KEY,
            help=(
                "Describe the supplied coefficient or contrast, including "
                "direction and reference level if known. The application does "
                "not infer them."
            ),
        )
        context_notes = st.text_area("Notes", key=CONTEXT_NOTES_KEY)

context_values = (
    dataset_title,
    organism,
    expression_scale,
    upstream_normalization,
    reference_annotation,
    feature_level,
    de_contrast,
    context_notes,
)
uploaded_provenance = (
    None
    if all(value == "" for value in context_values)
    else DatasetProvenance(
        dataset_title=dataset_title,
        organism=organism,
        expression_scale_description=expression_scale,
        upstream_normalization_method=upstream_normalization,
        reference_genome_annotation=reference_annotation,
        feature_level=feature_level,
        de_contrast_description=de_contrast,
        notes=context_notes,
    )
)

uploaded_sources = (expression_file, metadata_file, de_file)
missing_uploads = [
    label
    for label, source in zip(
        (
            "expression matrix",
            "sample metadata",
        ),
        uploaded_sources[:2],
        strict=True,
    )
    if source is None
]

if missing_uploads:
    st.info(
        "To activate your own upload, still required: "
        + ", ".join(missing_uploads)
        + ". The differential-expression file is optional."
    )
    if not demo_requested:
        clear_candidate_feedback(st.session_state)
elif not demo_requested:
    upload_label = uploaded_source_label(*uploaded_sources)
    with st.spinner("Validating uploaded files…"):
        uploaded_candidate = load_uploaded_candidate(*uploaded_sources)
    _activate_uploaded_candidate(
        uploaded_candidate,
        upload_label,
        uploaded_provenance,
    )

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
    if current.de_results is None:
        st.subheader("Differential-expression results")
        st.info(
            "No precomputed differential-expression results were supplied for "
            "this dataset. The other descriptive analysis pages remain available."
        )
    else:
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
    "precomputed differential-expression results is available when that optional "
    "table is supplied, together "
    "with exact, descriptive single-gene expression lookup. Dedicated CSV "
    "downloads of current descriptive result tables are available on the "
    "analysis pages. Differential-expression modelling and volcano plots are "
    "not implemented."
)
