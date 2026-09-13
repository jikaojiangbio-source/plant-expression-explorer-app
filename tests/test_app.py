"""Smoke test for the Streamlit home page."""

import csv
import io
import json
from pathlib import Path

import pandas as pd
from streamlit.testing.v1 import AppTest

import plant_expression_explorer.exports as exports_module
from plant_expression_explorer.consistency import validate_input_tables
from plant_expression_explorer.dataset import (
    CURRENT_DATASET_KEY,
    DE_RESULTS_UPLOAD_KEY,
    DEMO_SOURCE_LABEL,
    EXPRESSION_UPLOAD_KEY,
    METADATA_UPLOAD_KEY,
    DatasetBundle,
    build_dataset_bundle,
    load_demo_candidate,
)
from plant_expression_explorer.differential_expression import STATUS_COLUMN
from plant_expression_explorer.exports import (
    CsvExportError,
    CsvExportErrorReason,
    build_csv_export,
)
from plant_expression_explorer.provenance import (
    DEMO_PROVENANCE,
    NOT_SUPPLIED,
    DatasetProvenance,
)
from plant_expression_explorer.validation import (
    IssueCode,
    Severity,
    ValidationIssue,
    ValidationReport,
)


def _visible_text(app: AppTest) -> str:
    element_groups = (
        app.title,
        app.header,
        app.subheader,
        app.markdown,
        app.info,
        app.warning,
        app.success,
        app.caption,
    )
    return " ".join(
        str(element.value)
        for group in element_groups
        for element in group
    )


def _download_labels(app: AppTest) -> list[str]:
    return [button.label for button in app.get("download_button")]


def _demo_bundle() -> DatasetBundle:
    candidate = load_demo_candidate()
    assert candidate.tables is not None
    return build_dataset_bundle(
        candidate.tables,
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
        report=candidate.report,
    )


def _without_de_bundle(
    provenance: DatasetProvenance | None = None,
) -> DatasetBundle:
    bundle = _uploaded_bundle()
    report = validate_input_tables(bundle.expression, bundle.metadata, None)
    assert not report.has_errors
    return build_dataset_bundle(
        (bundle.expression, bundle.metadata, None),
        source="uploaded",
        source_label="Expression and metadata only",
        report=report,
        provenance=provenance,
    )


def _provenance_bundle() -> DatasetBundle:
    bundle = _uploaded_bundle()
    provenance = DatasetProvenance(
        dataset_title="  盐胁迫 RNA-seq  ",
        organism="Solanum lycopersicum 🍅",
        expression_scale_description="log2(TPM + 1)",
        upstream_normalization_method="TMM + exact supplied wording",
        reference_genome_annotation="SL4.0 / ITAG4.1",
        feature_level="gene",
        de_contrast_description="treated - control; control reference",
        notes="\t**display as text, not Markdown**  ",
    )
    return build_dataset_bundle(
        (bundle.expression, bundle.metadata, bundle.de_results),
        source=bundle.source,
        source_label=bundle.source_label,
        report=bundle.validation_report,
        provenance=provenance,
    )


def _uploaded_bundle(*, one_gene: bool = False) -> DatasetBundle:
    gene_ids = ["g1"] if one_gene else ["g1", "g2", "g3", "g4"]
    expression = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "sample_b": ["2.0"] if one_gene else ["2", "3", "4", "5"],
            "sample_a": ["1.0"] if one_gene else ["1", "2", "3", "4"],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b"],
            "condition": ["control", "treated"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "log2FoldChange": [0.0] * len(gene_ids),
            "pvalue": [0.5] * len(gene_ids),
            "padj": [0.5] * len(gene_ids),
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="User-uploaded CSV tables (test inputs)",
        report=ValidationReport(),
    )


def _uploaded_bundle_with_genotype() -> DatasetBundle:
    bundle = _uploaded_bundle()
    metadata = bundle.metadata.assign(genotype=["WT", "mutant"])
    report = validate_input_tables(bundle.expression, metadata, bundle.de_results)
    assert not report.has_errors
    return build_dataset_bundle(
        (bundle.expression, metadata, bundle.de_results),
        source=bundle.source,
        source_label=bundle.source_label,
        report=report,
    )


def _invalid_qc_bundle() -> DatasetBundle:
    bundle = _uploaded_bundle()
    expression = bundle.expression.copy(deep=True)
    expression.loc[0, "sample_a"] = "not-numeric"
    return build_dataset_bundle(
        (expression, bundle.metadata, bundle.de_results),
        source="uploaded",
        source_label="Defensive QC test bundle",
        report=ValidationReport(),
    )


def _constant_correlation_bundle(
    *,
    all_constant: bool = False,
) -> DatasetBundle:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "sample_b": [2, 3, 4, 5],
            "sample_a": (
                [1, 1, 1, 1] if all_constant else [1, 2, 3, 4]
            ),
            "sample_constant": [7, 7, 7, 7],
        }
    )
    if all_constant:
        expression["sample_b"] = [2, 2, 2, 2]
    metadata = pd.DataFrame(
        {
            "sample_id": [
                "sample_a",
                "sample_constant",
                "sample_b",
            ],
            "condition": ["control", "treated", "control"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": expression["gene_id"],
            "log2FoldChange": [0.0] * len(expression.index),
            "pvalue": [0.5] * len(expression.index),
            "padj": [0.5] * len(expression.index),
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="Constant-sample test bundle",
        report=ValidationReport(),
    )


def _one_sample_pca_bundle() -> DatasetBundle:
    expression = pd.DataFrame(
        {"gene_id": ["g1", "g2", "g3"], "sample_a": ["1", "2", "3"]}
    )
    metadata = pd.DataFrame({"sample_id": ["sample_a"], "condition": ["control"]})
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "log2FoldChange": [0.0] * 3,
            "pvalue": [0.5] * 3,
            "padj": [0.5] * 3,
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="One-sample PCA test bundle",
        report=ValidationReport(),
    )


def _zero_gene_bundle() -> DatasetBundle:
    expression = pd.DataFrame(columns=["gene_id", "sample_a", "sample_b"])
    metadata = pd.DataFrame(
        {"sample_id": ["sample_a", "sample_b"], "condition": ["control", "treated"]}
    )
    de_results = pd.DataFrame(
        columns=["gene_id", "log2FoldChange", "pvalue", "padj"]
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="Zero-gene test bundle",
        report=ValidationReport(),
    )


def _all_genes_constant_bundle() -> DatasetBundle:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3", "g4"],
            "sample_a": [10, 20, 30, 40],
            "sample_b": [10, 20, 30, 40],
            "sample_c": [10, 20, 30, 40],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": expression["gene_id"],
            "log2FoldChange": [0.0] * 4,
            "pvalue": [0.5] * 4,
            "padj": [0.5] * 4,
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="All-genes-constant test bundle",
        report=ValidationReport(),
    )


def _colinear_pca_bundle() -> DatasetBundle:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "sample_a": [1, 2],
            "sample_b": [2, 4],
            "sample_c": [5, 10],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2"],
            "log2FoldChange": [0.0, 0.0],
            "pvalue": [0.5, 0.5],
            "padj": [0.5, 0.5],
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="Colinear PCA test bundle",
        report=ValidationReport(),
    )


def _large_finite_value_pca_bundle() -> DatasetBundle:
    expression = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "sample_a": [1e300, 2.0, 3.0],
            "sample_b": [-1e300, 4.0, 6.0],
            "sample_c": [0.0, 5.0, 9.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b", "sample_c"],
            "condition": ["control", "control", "treated"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1", "g2", "g3"],
            "log2FoldChange": [0.0] * 3,
            "pvalue": [0.5] * 3,
            "padj": [0.5] * 3,
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="Large finite value PCA test bundle",
        report=ValidationReport(),
    )


def _run_qc_page(bundle: DatasetBundle | None = None) -> AppTest:
    app = AppTest.from_file("pages/2_Sample_Quality_Control.py")
    if bundle is not None:
        app.session_state[CURRENT_DATASET_KEY] = bundle
    return app.run()


def _run_correlation_page(bundle: DatasetBundle | None = None) -> AppTest:
    app = AppTest.from_file("pages/4_Sample_Correlation.py")
    if bundle is not None:
        app.session_state[CURRENT_DATASET_KEY] = bundle
    return app.run()


def _correlation_session_keys(app: AppTest) -> list[str]:
    return [
        str(key)
        for key in app.session_state.filtered_state
        if "correlation" in str(key).lower()
    ]


def _run_pca_page(bundle: DatasetBundle | None = None) -> AppTest:
    app = AppTest.from_file("pages/3_PCA.py")
    if bundle is not None:
        app.session_state[CURRENT_DATASET_KEY] = bundle
    return app.run()


def _de_uploaded_bundle(*, gene_mismatch: bool = False) -> DatasetBundle:
    expression_gene_ids = ["g1", "g2", "g3", "g4", "g5"]
    de_gene_ids = (
        ["g1", "g2", "g3", "g4", "g6"]
        if gene_mismatch
        else expression_gene_ids
    )
    expression = pd.DataFrame(
        {
            "gene_id": expression_gene_ids,
            "sample_a": [1, 2, 3, 4, 5],
            "sample_b": [2, 3, 4, 5, 6],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_b"],
            "condition": ["control", "treated"],
        }
    )
    de_results = pd.DataFrame(
        {
            "gene_id": de_gene_ids,
            "annotation": ["a", "b", "c", "d", "e"],
            "log2FoldChange": [1.0, -1.0, 0.5, None, 1.5],
            "pvalue": [0.01, 0.02, 0.4, 0.8, None],
            "padj": [0.05, 0.05, 0.5, 0.2, 0.03],
        },
        index=pd.Index([9, 4, 9, 2, 7], name="source_index"),
    )
    report = validate_input_tables(expression, metadata, de_results)
    assert not report.has_errors
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="User-uploaded CSV tables (Phase 8 test inputs)",
        report=report,
    )


def _blocking_de_bundle() -> DatasetBundle:
    bundle = _uploaded_bundle()
    report = ValidationReport(
        (
            ValidationIssue(
                code=IssueCode.VALUE_OUT_OF_RANGE,
                severity=Severity.ERROR,
                table="Differential-expression results",
                column="padj",
                message="Defensive blocking test Error.",
            ),
        )
    )
    return DatasetBundle(
        expression=bundle.expression,
        metadata=bundle.metadata,
        de_results=bundle.de_results,
        source=bundle.source,
        source_label="Defensive differential-expression test bundle",
        validation_report=report,
    )


def _high_precision_de_bundle() -> DatasetBundle:
    bundle = _de_uploaded_bundle()
    de_results = bundle.de_results.copy(deep=True)
    de_results["log2FoldChange"] = [
        1.0000001,
        -1.0000001,
        1.00000009,
        -1.0000001,
        None,
    ]
    de_results["padj"] = [
        0.050000001,
        0.050000001,
        0.050000001,
        0.0500000011,
        0.03,
    ]
    report = validate_input_tables(bundle.expression, bundle.metadata, de_results)
    assert not report.has_errors
    return build_dataset_bundle(
        (bundle.expression, bundle.metadata, de_results),
        source="uploaded",
        source_label="High-precision Phase 8 test inputs",
        report=report,
    )


def _run_de_page(bundle: DatasetBundle | None = None) -> AppTest:
    app = AppTest.from_file("pages/5_Differential_Expression.py")
    if bundle is not None:
        app.session_state[CURRENT_DATASET_KEY] = bundle
    return app.run()


def _gene_uploaded_bundle() -> DatasetBundle:
    expression = pd.DataFrame(
        {
            "gene_id": ["GeneA", "genea"],
            "sample_b": ["2.0", "8.0"],
            "sample_a": ["1.0", "4.0"],
            "sample_c": ["3.0", "6.0"],
        },
        index=pd.Index([7, 3], name="source_index"),
    )
    metadata = pd.DataFrame(
        {
            "sample_id": ["sample_a", "sample_c", "sample_b"],
            "condition": ["Control", "Treated", "Control"],
        },
        index=[20, 30, 10],
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["GeneA", "genea"],
            "log2FoldChange": [0.0, 0.0],
            "pvalue": [0.5, 0.5],
            "padj": [0.5, 0.5],
        }
    )
    report = validate_input_tables(expression, metadata, de_results)
    assert not report.has_errors
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="User-uploaded CSV tables (Phase 9 test inputs)",
        report=report,
    )


def _invalid_gene_options_bundle() -> DatasetBundle:
    bundle = _gene_uploaded_bundle()
    expression = bundle.expression.copy(deep=True)
    expression.loc[expression.index[0], "sample_a"] = "not-numeric"
    return build_dataset_bundle(
        (expression, bundle.metadata, bundle.de_results),
        source="uploaded",
        source_label="Defensive invalid gene-options bundle",
        report=ValidationReport(),
    )


def _invalid_gene_lookup_bundle() -> DatasetBundle:
    bundle = _gene_uploaded_bundle()
    metadata = bundle.metadata.drop(columns="condition")
    return build_dataset_bundle(
        (bundle.expression, metadata, bundle.de_results),
        source="uploaded",
        source_label="Defensive invalid gene-lookup bundle",
        report=ValidationReport(),
    )


def _one_sample_gene_bundle() -> DatasetBundle:
    expression = pd.DataFrame({"gene_id": ["g1"], "only": [5.0]})
    metadata = pd.DataFrame({"sample_id": ["only"], "condition": ["Control"]})
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [0.0],
            "pvalue": [0.5],
            "padj": [0.5],
        }
    )
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="Defensive one-sample gene bundle",
        report=ValidationReport(),
    )


def _numerical_range_gene_bundle() -> DatasetBundle:
    maximum = float.fromhex("0x1.fffffffffffffp+1023")
    expression = pd.DataFrame(
        {"gene_id": ["g1"], "s1": [maximum], "s2": [-maximum]}
    )
    metadata = pd.DataFrame(
        {"sample_id": ["s1", "s2"], "condition": ["A", "A"]}
    )
    de_results = pd.DataFrame(
        {
            "gene_id": ["g1"],
            "log2FoldChange": [0.0],
            "pvalue": [0.5],
            "padj": [0.5],
        }
    )
    report = validate_input_tables(expression, metadata, de_results)
    assert not report.has_errors
    return build_dataset_bundle(
        (expression, metadata, de_results),
        source="uploaded",
        source_label="Numerical-range gene test bundle",
        report=report,
    )


def _run_gene_page(bundle: DatasetBundle | None = None) -> AppTest:
    app = AppTest.from_file("pages/6_Gene_Expression.py")
    if bundle is not None:
        app.session_state[CURRENT_DATASET_KEY] = bundle
    return app.run()


def _run_selected_gene_page(bundle: DatasetBundle, gene_id: str) -> AppTest:
    app = _run_gene_page(bundle)
    return app.selectbox[0].select(gene_id).run()


def _metric_values(app: AppTest) -> dict[str, str]:
    return {metric.label: metric.value for metric in app.metric}


def _pca_session_keys(app: AppTest) -> list[str]:
    return [
        str(key)
        for key in app.session_state.filtered_state
        if "pca" in str(key).lower()
    ]


def _pca_scatter_chart_count(app: AppTest) -> int:
    # st.bar_chart also renders as a vega_lite_chart element (mark "bar");
    # the PC1-versus-PC2 scatter is the only one using a "circle" mark.
    return sum(
        1
        for chart in app.get("vega_lite_chart")
        if '"circle"' in chart.proto.spec
    )


def test_home_page_loads_and_describes_scope() -> None:
    app = AppTest.from_file("app.py").run()

    assert not app.exception
    assert app.title[0].value == "🌱 Plant Expression Explorer"
    visible_text = " ".join(
        element.value for element in [*app.subheader, *app.markdown, *app.info]
    )
    assert "preprocessed plant transcriptomics data" in visible_text
    assert "does not process FASTQ files" in visible_text


def test_home_page_workflow_order_matches_sidebar_page_sequence() -> None:
    app = AppTest.from_file("app.py").run()

    workflow_markdown = next(
        element.value for element in app.markdown if "Upload data" in element.value
    )

    def position(term: str) -> int:
        index = workflow_markdown.lower().find(term.lower())
        assert index != -1, f"{term!r} not found in workflow markdown"
        return index

    upload_position = position("Upload data")
    quality_control_position = position("quality-control")
    pca_position = position("PCA")
    correlation_position = position("sample correlation")
    differential_expression_position = position("differential-expression")
    gene_expression_position = position("gene expression")

    assert (
        upload_position
        < quality_control_position
        < pca_position
        < correlation_position
        < differential_expression_position
        < gene_expression_position
    )


def test_home_page_no_longer_describes_pca_as_planned() -> None:
    app = AppTest.from_file("app.py").run()

    workflow_markdown = next(
        element.value for element in app.markdown if "Upload data" in element.value
    )
    pca_line = next(
        line for line in workflow_markdown.splitlines() if "PCA" in line
    )
    assert "planned" not in pca_line.lower()
    assert "descriptive" in pca_line.lower()

    info_text = " ".join(element.value for element in app.info)
    assert "PCA" not in info_text


def test_home_page_has_pca_page_link() -> None:
    app = AppTest.from_file("app.py").run()

    page_links = app.get("page_link")
    matching = [link for link in page_links if link.proto.page == "PCA"]
    assert len(matching) == 1
    assert matching[0].proto.label == "Review PCA"


def test_home_page_analysis_link_order_includes_phase_9_gene_expression() -> None:
    app = AppTest.from_file("app.py").run()

    page_links = app.get("page_link")
    pages_in_order = [link.proto.page for link in page_links]

    assert pages_in_order == [
        "Upload_Data",
        "Sample_Quality_Control",
        "PCA",
        "Sample_Correlation",
        "Differential_Expression",
        "Gene_Expression",
    ]


def test_home_page_lists_phase_12_scientific_context_as_available() -> None:
    app = AppTest.from_file("app.py").run()

    workflow_markdown = next(
        element.value for element in app.markdown if "Upload data" in element.value
    )
    deg_line = next(
        line
        for line in workflow_markdown.splitlines()
        if "differential-expression" in line.lower()
    )
    gene_expression_line = next(
        line
        for line in workflow_markdown.splitlines()
        if "gene expression" in line.lower()
    )
    assert "planned" not in deg_line.lower()
    assert "descriptive thresholds" in deg_line.lower()
    assert "planned" not in gene_expression_line.lower()
    assert "exact supplied gene" in gene_expression_line.lower()
    assert "dedicated utf-8 csv" in workflow_markdown.lower()

    info_text = " ".join(element.value for element in app.info).lower()
    assert "differential-expression exploration" not in info_text
    assert "volcano plots" in info_text
    assert "dedicated exports" not in info_text
    assert "gene lookup" not in info_text
    visible_text = _visible_text(app).lower()
    assert "phase 12" in visible_text
    assert "differential-expression results are now optional" in visible_text
    assert "descriptive dataset context" in visible_text


def test_upload_page_initial_state_and_synthetic_disclaimer() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py").run()

    assert not app.exception
    assert app.title[0].value == "📤 Upload Data"
    visible_text = _visible_text(app)
    assert "No validated dataset is currently active" in visible_text
    assert "values are synthetic" in visible_text
    assert "gene IDs are fictional" in visible_text
    assert "p-values are constructed" in visible_text
    assert "DESeq2" in visible_text
    assert "replicate noise is balanced within each condition" in visible_text
    assert "neither pattern is guaranteed in real experiments" in visible_text
    assert "differential-expression file is optional" in visible_text
    assert len(app.text_input) == 7
    assert len(app.text_area) == 1
    assert "Sample Quality Control, PCA, and Sample Correlation are available" in visible_text
    assert "single-gene expression lookup" in visible_text


def test_upload_page_offers_example_csv_template_downloads() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py").run()

    assert not app.exception
    template_buttons = app.get("download_button")
    assert [button.label for button in template_buttons] == [
        "Expression matrix template",
        "Sample metadata template",
        "DE results template",
    ]
    visible_text = _visible_text(app)
    assert "fabricated placeholder values" in visible_text
    assert "never loaded as a dataset" in visible_text


def test_upload_page_lists_phase_8_exploration_as_available() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py").run()

    assert not app.exception
    visible_text = _visible_text(app)

    for unimplemented_phrase in (
        "Sample QC, PCA, sample correlation",
        "Analysis features are not implemented yet",
    ):
        assert unimplemented_phrase not in visible_text

    assert "Sample Quality Control, PCA, and Sample Correlation are available" in visible_text
    assert "exploratory threshold classification" in visible_text
    assert "supplied, precomputed differential-expression results" in visible_text
    assert "exact, descriptive single-gene expression lookup" in visible_text
    assert "Dedicated CSV downloads of current descriptive result tables" in visible_text
    assert "Differential-expression modelling and volcano plots are not implemented" in visible_text


def test_upload_page_demo_load_and_reset_are_repeatable() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py").run()

    app.button[0].click().run()

    assert not app.exception
    bundle = app.session_state[CURRENT_DATASET_KEY]
    assert isinstance(bundle, DatasetBundle)
    assert bundle.source == "demo"
    assert bundle.gene_count == 120
    assert bundle.sample_count == 6
    assert bundle.de_row_count == 120
    assert bundle.has_de_results is True
    assert bundle.provenance == DEMO_PROVENANCE
    visible_text = _visible_text(app)
    assert "A validated dataset is active" in visible_text
    assert "Bundled synthetic demonstration data" in visible_text
    overview_metrics = {metric.label: metric.value for metric in app.metric}
    assert overview_metrics["Genes"] == "120"
    assert overview_metrics["Samples"] == "6"
    assert overview_metrics["DE rows"] == "120"
    assert DEMO_PROVENANCE.expression_scale_description in [
        element.value for element in app.get("code")
    ]

    app.button[1].click().run()

    assert not app.exception
    assert CURRENT_DATASET_KEY not in app.session_state
    assert app.session_state[EXPRESSION_UPLOAD_KEY] is None
    assert app.session_state[METADATA_UPLOAD_KEY] is None
    assert app.session_state[DE_RESULTS_UPLOAD_KEY] is None
    assert "No validated dataset is currently active" in _visible_text(app)

    app.button[0].click().run()
    app.button[1].click().run()

    assert not app.exception
    assert CURRENT_DATASET_KEY not in app.session_state


def test_upload_page_overview_card_shows_shape_and_condition_distribution() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py")
    app.session_state[CURRENT_DATASET_KEY] = _uploaded_bundle()
    app.run()

    assert not app.exception
    overview_metrics = {metric.label: metric.value for metric in app.metric}
    assert overview_metrics["Genes"] == "4"
    assert overview_metrics["Samples"] == "2"
    condition_chart = next(
        chart for chart in app.get("vega_lite_chart") if '"bar"' in chart.proto.spec
    )
    assert '"field": "condition"' in condition_chart.proto.spec
    assert '"field": "sample_count"' in condition_chart.proto.spec


def test_upload_page_overview_card_reports_de_not_supplied_when_omitted() -> None:
    bundle = _uploaded_bundle()
    metadata = bundle.metadata
    expression = bundle.expression
    report = validate_input_tables(expression, metadata, None)
    assert not report.has_errors
    without_de = build_dataset_bundle(
        (expression, metadata, None),
        source=bundle.source,
        source_label=bundle.source_label,
        report=report,
    )
    app = AppTest.from_file("pages/1_Upload_Data.py")
    app.session_state[CURRENT_DATASET_KEY] = without_de
    app.run()

    assert not app.exception
    overview_metrics = {metric.label: metric.value for metric in app.metric}
    assert overview_metrics["DE rows"] == "Not supplied"


def test_qc_page_no_data_state_is_clear_and_does_not_create_dataset_state() -> None:
    app = _run_qc_page()

    assert not app.exception
    assert app.title[0].value == "🧪 Sample Quality Control"
    visible_text = _visible_text(app)
    assert "No validated dataset is currently loaded" in visible_text
    assert "Upload Data" in visible_text
    assert "Matrix overview" not in visible_text
    assert CURRENT_DATASET_KEY not in app.session_state


def test_qc_page_demo_state_renders_summaries_charts_and_disclaimer() -> None:
    app = _run_qc_page(_demo_bundle())

    assert not app.exception
    visible_text = _visible_text(app)
    assert "Bundled synthetic demonstration data" in visible_text
    assert "120 genes, 6 samples" in visible_text
    assert "Matrix overview" in visible_text
    assert "Condition and replicate summary" in visible_text
    assert "Per-sample statistics" in visible_text
    assert "Descriptive charts" in visible_text
    assert "Median expression by sample" in visible_text
    assert "Interquartile range by sample" in visible_text
    assert "Sample count by condition" in visible_text
    assert len(app.dataframe) == 2
    assert app.dataframe[0].value["condition"].tolist() == [
        "Control",
        "High_nitrate",
    ]
    assert app.dataframe[1].value["sample_id"].tolist() == [
        "Control_1",
        "Control_2",
        "Control_3",
        "High_nitrate_1",
        "High_nitrate_2",
        "High_nitrate_3",
    ]
    assert "Every sample has a zero-value count of 0" in visible_text
    assert "Every sample has a negative-value count of 0" in visible_text
    assert "values are synthetic" in visible_text
    assert "gene IDs are fictional" in visible_text
    assert "p-values are constructed" in visible_text
    assert "does not establish real QC thresholds" in visible_text
    prohibited = visible_text.lower()
    assert "outlier" not in prohibited
    assert "failed sample" not in prohibited
    assert "poor-quality sample" not in prohibited
    assert "analysis passed" not in prohibited
    assert "analysis failed" not in prohibited


def test_qc_page_uploaded_state_maps_conditions_without_mutating_bundle() -> None:
    bundle = _uploaded_bundle()
    original_expression = bundle.expression.copy(deep=True)
    original_metadata = bundle.metadata.copy(deep=True)

    app = _run_qc_page(bundle)

    assert not app.exception
    visible_text = _visible_text(app)
    assert "User-uploaded CSV tables (test inputs)" in visible_text
    assert "Source type:** uploaded" in visible_text
    assert "4 genes, 2 samples" in visible_text
    assert app.dataframe[1].value["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
    ]
    assert app.dataframe[1].value["condition"].tolist() == [
        "treated",
        "control",
    ]
    assert "differs; conditions were mapped by exact sample ID" in visible_text
    assert "values are synthetic" not in visible_text
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    pd.testing.assert_frame_equal(bundle.metadata, original_metadata)


def test_qc_page_one_gene_state_displays_na_without_mutating_result_source() -> None:
    bundle = _uploaded_bundle(one_gene=True)
    original_expression = bundle.expression.copy(deep=True)

    app = _run_qc_page(bundle)

    assert not app.exception
    visible_text = _visible_text(app)
    sample_display = app.dataframe[1].value
    assert sample_display["standard_deviation"].tolist() == ["N/A", "N/A"]
    assert (
        "The matrix contains only one gene, so sample-spread and "
        "constant-column diagnostics are not informative."
    ) in visible_text
    assert "this diagnostic is not informative" in visible_text
    pd.testing.assert_frame_equal(bundle.expression, original_expression)


def test_qc_page_collapses_limitations_into_an_expander() -> None:
    app = _run_qc_page(_demo_bundle())

    assert not app.exception
    expander_labels = [element.label for element in app.get("expander")]
    assert "Scientific and statistical limitations" in expander_labels
    visible_text = _visible_text(app)
    assert "scale-dependent meanings" in visible_text


def test_qc_page_offers_no_grouping_selector_without_additional_metadata() -> None:
    app = _run_qc_page(_uploaded_bundle())

    assert not app.exception
    assert len(app.selectbox) == 0
    assert app.dataframe[1].value["condition"].tolist() == ["treated", "control"]


def test_qc_page_grouping_selector_relabels_summaries_without_recomputation() -> None:
    bundle = _uploaded_bundle_with_genotype()

    app = _run_qc_page(bundle)
    assert len(app.selectbox) == 1
    selector = app.selectbox[0]
    assert selector.options == ["condition", "genotype"]
    assert selector.value == "condition"

    app.selectbox[0].select("genotype").run()

    assert not app.exception
    per_sample = app.dataframe[1].value
    assert "genotype" in per_sample.columns
    assert "condition" not in per_sample.columns
    assert dict(zip(per_sample["sample_id"], per_sample["genotype"])) == {
        "sample_a": "WT",
        "sample_b": "mutant",
    }
    visible_text = _visible_text(app)
    assert "Grouped by metadata column 'genotype'" in visible_text


def test_qc_page_shows_only_controlled_qc_error_without_traceback() -> None:
    app = _run_qc_page(_invalid_qc_bundle())

    assert not app.exception
    assert len(app.error) == 1
    error_text = app.error[0].value
    assert "NON_COERCIBLE_EXPRESSION_VALUE" in error_text
    assert "1 non-coercible sample value" in error_text
    assert "Traceback" not in error_text


def test_correlation_page_no_data_state_is_clear_and_has_no_side_effects() -> None:
    app = _run_correlation_page()

    assert not app.exception
    assert app.title[0].value == "🔥 Sample Correlation"
    visible_text = _visible_text(app)
    assert "No validated dataset is currently loaded" in visible_text
    assert "Upload Data" in visible_text
    assert "Correlation matrix" not in visible_text
    assert CURRENT_DATASET_KEY not in app.session_state
    assert _correlation_session_keys(app) == []


def test_correlation_page_demo_renders_matrix_heatmap_and_summaries() -> None:
    app = _run_correlation_page(_demo_bundle())

    assert not app.exception
    visible_text = _visible_text(app)
    assert "Bundled synthetic demonstration data" in visible_text
    assert "Source type:** demo" in visible_text
    assert "120 genes and 6 samples" in visible_text
    assert "Pearson correlation" in visible_text
    assert "Correlation matrix" in visible_text
    assert "Correlation heatmap" in visible_text
    assert "Per-sample correlation summary" in visible_text
    assert "Unique sample-pair summary" in visible_text
    assert "Condition-pair descriptive summary" in visible_text
    assert len(app.dataframe) == 4
    assert app.dataframe[0].value.shape == (6, 6)
    assert app.dataframe[1].value["sample_id"].tolist() == [
        "Control_1",
        "Control_2",
        "Control_3",
        "High_nitrate_1",
        "High_nitrate_2",
        "High_nitrate_3",
    ]
    assert len(app.dataframe[2].value.index) == 15
    assert app.dataframe[3].value[
        ["condition_a", "condition_b"]
    ].values.tolist() == [
        ["Control", "Control"],
        ["Control", "High_nitrate"],
        ["High_nitrate", "High_nitrate"],
    ]
    assert len(app.get("vega_lite_chart")) == 1
    assert "values are synthetic" in visible_text
    assert "gene IDs are fictional" in visible_text
    assert "p-values are constructed" in visible_text
    assert "does not define real correlation thresholds" in visible_text
    assert _correlation_session_keys(app) == []


def test_correlation_page_uploaded_maps_conditions_and_preserves_bundle() -> None:
    bundle = _uploaded_bundle()
    original_expression = bundle.expression.copy(deep=True)
    original_metadata = bundle.metadata.copy(deep=True)

    app = _run_correlation_page(bundle)

    assert not app.exception
    visible_text = _visible_text(app)
    assert "User-uploaded CSV tables (test inputs)" in visible_text
    assert "Source type:** uploaded" in visible_text
    assert "4 genes and 2 samples" in visible_text
    assert app.dataframe[0].value.index.tolist() == [
        "sample_b",
        "sample_a",
    ]
    assert app.dataframe[0].value.columns.tolist() == [
        "sample_b",
        "sample_a",
    ]
    assert app.dataframe[1].value["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
    ]
    assert app.dataframe[1].value["condition"].tolist() == [
        "treated",
        "control",
    ]
    assert "mapped by exact sample ID" in visible_text
    assert len(app.dataframe[2].value.index) == 1
    assert app.dataframe[3].value[
        ["condition_a", "condition_b"]
    ].values.tolist() == [
        ["treated", "treated"],
        ["treated", "control"],
        ["control", "control"],
    ]
    assert len(app.get("vega_lite_chart")) == 1
    assert "values are synthetic" not in visible_text
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    pd.testing.assert_frame_equal(bundle.metadata, original_metadata)
    assert _correlation_session_keys(app) == []


def test_correlation_page_collapses_method_and_limitations_into_expanders() -> None:
    app = _run_correlation_page(_demo_bundle())

    assert not app.exception
    expander_labels = [element.label for element in app.get("expander")]
    assert "Method" in expander_labels
    assert "Scientific and statistical limitations" in expander_labels
    visible_text = _visible_text(app)
    assert "Pearson correlation is calculated between every sample pair" in visible_text
    assert "does not prove replicate validity" in visible_text


def test_correlation_page_offers_no_grouping_selector_without_additional_metadata() -> (
    None
):
    app = _run_correlation_page(_uploaded_bundle())

    assert not app.exception
    assert len(app.selectbox) == 0


def test_correlation_page_grouping_selector_relabels_summaries_without_recomputation() -> (
    None
):
    bundle = _uploaded_bundle_with_genotype()

    app = _run_correlation_page(bundle)
    assert len(app.selectbox) == 1
    selector = app.selectbox[0]
    assert selector.options == ["condition", "genotype"]
    assert selector.value == "condition"

    app.selectbox[0].select("genotype").run()

    assert not app.exception
    pair_summary = app.dataframe[2].value
    assert "genotype_a" in pair_summary.columns
    assert "genotype_b" in pair_summary.columns
    assert "condition_a" not in pair_summary.columns
    assert pair_summary[["genotype_a", "genotype_b"]].values.tolist() == [
        ["mutant", "WT"]
    ]
    condition_pair_summary = app.dataframe[3].value
    assert "genotype_a" in condition_pair_summary.columns
    visible_text = _visible_text(app)
    assert "Grouped by metadata column 'genotype'" in visible_text


def test_correlation_page_constant_sample_retains_undefined_values() -> None:
    bundle = _constant_correlation_bundle()
    original_expression = bundle.expression.copy(deep=True)

    app = _run_correlation_page(bundle)

    assert not app.exception
    visible_text = _visible_text(app)
    assert "sample_constant" in visible_text
    assert "undefined Pearson correlations" in visible_text
    assert "Undefined unique non-self sample pairs:** 2" in visible_text
    matrix = app.dataframe[0].value
    assert matrix.loc["sample_constant"].eq("N/A").all()
    assert matrix.loc[:, "sample_constant"].eq("N/A").all()
    sample_summary = app.dataframe[1].value
    constant_row = sample_summary.loc[
        sample_summary["sample_id"].eq("sample_constant")
    ].iloc[0]
    assert constant_row["defined_pair_count"] == 0
    assert constant_row["undefined_pair_count"] == 2
    assert constant_row[
        [
            "minimum_correlation",
            "median_correlation",
            "mean_correlation",
            "maximum_correlation",
        ]
    ].eq("N/A").all()
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    assert _correlation_session_keys(app) == []


def test_correlation_page_all_constant_state_renders_without_exception() -> None:
    app = _run_correlation_page(
        _constant_correlation_bundle(all_constant=True)
    )

    assert not app.exception
    visible_text = _visible_text(app)
    assert "3 constant sample column(s)" in visible_text
    assert "Undefined unique non-self sample pairs:** 3" in visible_text
    assert app.dataframe[0].value.eq("N/A").all().all()
    assert app.dataframe[1].value["defined_pair_count"].eq(0).all()
    assert app.dataframe[1].value["undefined_pair_count"].eq(2).all()
    assert len(app.get("vega_lite_chart")) == 1
    assert _correlation_session_keys(app) == []


def test_correlation_page_one_gene_shows_only_controlled_error() -> None:
    app = _run_correlation_page(_uploaded_bundle(one_gene=True))

    assert not app.exception
    assert len(app.error) == 1
    error_text = app.error[0].value
    assert "INSUFFICIENT_GENE_ROWS" in error_text
    assert "1 gene row" in error_text
    assert "Traceback" not in error_text
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    assert _correlation_session_keys(app) == []


def test_correlation_page_has_no_prohibited_classification_wording() -> None:
    for bundle in (
        _demo_bundle(),
        _uploaded_bundle(),
        _constant_correlation_bundle(),
    ):
        app = _run_correlation_page(bundle)
        assert not app.exception
        visible_text = _visible_text(app).lower()
        for prohibited in (
            "confirmed outlier",
            "failed sample",
            "poor-quality sample",
            "correlation passed",
            "correlation failed",
            "remove this sample",
        ):
            assert prohibited not in visible_text


def test_pca_page_no_data_state_is_clear_and_has_no_side_effects() -> None:
    app = _run_pca_page()

    assert not app.exception
    assert app.title[0].value == "📊 PCA"
    visible_text = _visible_text(app)
    assert "No validated dataset is currently loaded" in visible_text
    assert "Upload Data" in visible_text
    assert "Explained variance" not in visible_text
    assert CURRENT_DATASET_KEY not in app.session_state
    assert _pca_session_keys(app) == []


def test_pca_page_demo_renders_scores_and_variance() -> None:
    app = _run_pca_page(_demo_bundle())

    assert not app.exception
    visible_text = _visible_text(app)
    assert "Bundled synthetic demonstration data" in visible_text
    assert "Source type:** demo" in visible_text
    assert "120 genes and 6 samples" in visible_text
    assert "Explained variance" in visible_text
    assert "Sample scores" in visible_text
    assert "PC1-versus-PC2 sample plot" in visible_text
    assert len(app.dataframe) == 2
    variance_table = app.dataframe[0].value
    assert variance_table["component"].tolist() == ["PC1", "PC2", "PC3", "PC4", "PC5"]
    score_table = app.dataframe[1].value
    assert score_table["sample_id"].tolist() == [
        "Control_1",
        "Control_2",
        "Control_3",
        "High_nitrate_1",
        "High_nitrate_2",
        "High_nitrate_3",
    ]
    assert _pca_scatter_chart_count(app) == 1
    assert "values are synthetic" in visible_text
    assert "gene IDs are fictional" in visible_text
    assert "p-values are constructed" in visible_text
    assert "does not define real PCA structure" in visible_text
    assert _pca_session_keys(app) == []


def test_pca_page_uploaded_maps_conditions_preserves_bundle_and_omits_scatter() -> None:
    bundle = _uploaded_bundle()
    original_expression = bundle.expression.copy(deep=True)
    original_metadata = bundle.metadata.copy(deep=True)

    app = _run_pca_page(bundle)

    assert not app.exception
    visible_text = _visible_text(app)
    assert "User-uploaded CSV tables (test inputs)" in visible_text
    assert "Source type:** uploaded" in visible_text
    assert "4 genes and 2 samples" in visible_text
    score_table = app.dataframe[1].value
    assert score_table["sample_id"].tolist() == ["sample_b", "sample_a"]
    assert score_table["condition"].tolist() == ["treated", "control"]
    assert list(score_table.columns) == ["sample_id", "condition", "pc1"]
    assert "Exactly two samples are present" in visible_text
    assert "PC1-versus-PC2 plot is not shown" in visible_text
    assert _pca_scatter_chart_count(app) == 0
    assert "values are synthetic" not in visible_text
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    pd.testing.assert_frame_equal(bundle.metadata, original_metadata)
    assert _pca_session_keys(app) == []


def test_pca_page_all_genes_zero_variance_shows_only_controlled_error() -> None:
    bundle = _all_genes_constant_bundle()
    original_expression = bundle.expression.copy(deep=True)

    app = _run_pca_page(bundle)

    assert not app.exception
    assert len(app.error) == 1
    error_text = app.error[0].value
    assert "ZERO_TOTAL_VARIANCE" in error_text
    assert "no variation" in error_text.lower()
    assert "invalid" not in error_text.lower()
    assert "Traceback" not in error_text
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    assert app.session_state[CURRENT_DATASET_KEY] is bundle
    assert _pca_session_keys(app) == []


def test_pca_page_large_finite_values_shows_only_controlled_numerical_range_error() -> None:
    bundle = _large_finite_value_pca_bundle()
    original_expression = bundle.expression.copy(deep=True)
    original_metadata = bundle.metadata.copy(deep=True)

    app = _run_pca_page(bundle)

    assert not app.exception
    assert len(app.error) == 1
    error_text = app.error[0].value
    assert "NUMERICAL_RANGE_ERROR" in error_text
    assert "float64" in error_text
    assert "not modified" in error_text
    assert "Traceback" not in error_text
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    pd.testing.assert_frame_equal(bundle.metadata, original_metadata)
    assert app.session_state[CURRENT_DATASET_KEY] is bundle
    assert _pca_session_keys(app) == []


def test_pca_page_one_gene_renders_with_observation_not_error() -> None:
    app = _run_pca_page(_uploaded_bundle(one_gene=True))

    assert not app.exception
    assert len(app.error) == 0
    visible_text = _visible_text(app)
    assert "Only one gene is present" in visible_text
    assert len(app.dataframe) == 2
    assert _pca_session_keys(app) == []


def test_pca_page_one_sample_shows_only_controlled_error_without_traceback() -> None:
    bundle = _one_sample_pca_bundle()
    original_expression = bundle.expression.copy(deep=True)
    original_metadata = bundle.metadata.copy(deep=True)

    app = _run_pca_page(bundle)

    assert not app.exception
    assert len(app.error) == 1
    error_text = app.error[0].value
    assert "INSUFFICIENT_SAMPLE_COUNT" in error_text
    assert "at least two sample columns" in error_text.lower()
    assert "Traceback" not in error_text
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    pd.testing.assert_frame_equal(bundle.expression, original_expression)
    pd.testing.assert_frame_equal(bundle.metadata, original_metadata)
    assert app.session_state[CURRENT_DATASET_KEY] is bundle
    assert _pca_session_keys(app) == []


def test_pca_page_zero_genes_shows_only_controlled_error() -> None:
    bundle = _zero_gene_bundle()

    app = _run_pca_page(bundle)

    assert not app.exception
    assert len(app.error) == 1
    error_text = app.error[0].value
    assert "EMPTY_EXPRESSION_TABLE" in error_text
    assert "Traceback" not in error_text
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    assert app.session_state[CURRENT_DATASET_KEY] is bundle
    assert _pca_session_keys(app) == []


def test_pca_page_component2_zero_variance_shows_plot_with_observation() -> None:
    app = _run_pca_page(_colinear_pca_bundle())

    assert not app.exception
    visible_text = _visible_text(app)
    assert _pca_scatter_chart_count(app) == 1
    assert "Component 2 accounts for 0% of total variance" in visible_text
    variance_table = app.dataframe[0].value
    second_row = variance_table.loc[variance_table["component"] == "PC2"].iloc[0]
    assert second_row["explained_variance_ratio"] == 0.0
    assert _pca_session_keys(app) == []


def _colinear_pca_bundle_with_genotype() -> DatasetBundle:
    bundle = _colinear_pca_bundle()
    metadata = bundle.metadata.assign(genotype=["WT", "mutant", "mutant"])
    report = validate_input_tables(bundle.expression, metadata, bundle.de_results)
    assert not report.has_errors
    return build_dataset_bundle(
        (bundle.expression, metadata, bundle.de_results),
        source=bundle.source,
        source_label=bundle.source_label,
        report=report,
    )


def test_pca_page_collapses_method_and_limitations_into_expanders() -> None:
    app = _run_pca_page(_demo_bundle())

    assert not app.exception
    expander_labels = [element.label for element in app.get("expander")]
    assert "Method" in expander_labels
    assert "Scientific and statistical limitations" in expander_labels
    visible_text = _visible_text(app)
    assert "mean-centred across samples" in visible_text
    assert "does not calculate a hypothesis" in visible_text


def test_pca_page_offers_no_grouping_selector_without_additional_metadata() -> None:
    app = _run_pca_page(_colinear_pca_bundle())

    assert not app.exception
    assert len(app.selectbox) == 0
    assert _pca_scatter_chart_count(app) == 1


def test_pca_page_grouping_selector_defaults_to_condition() -> None:
    app = _run_pca_page(_colinear_pca_bundle_with_genotype())

    assert not app.exception
    assert len(app.selectbox) == 1
    selector = app.selectbox[0]
    assert selector.options == ["condition", "genotype"]
    assert selector.value == "condition"
    chart = next(
        chart for chart in app.get("vega_lite_chart") if '"circle"' in chart.proto.spec
    )
    assert '"field": "condition"' in chart.proto.spec


def test_pca_page_grouping_selector_switches_chart_colour_field() -> None:
    app = _run_pca_page(_colinear_pca_bundle_with_genotype())

    app.selectbox[0].select("genotype").run()

    assert not app.exception
    chart = next(
        chart for chart in app.get("vega_lite_chart") if '"circle"' in chart.proto.spec
    )
    assert '"field": "genotype"' in chart.proto.spec
    assert "Genotype" in chart.proto.spec


def test_pca_page_has_no_prohibited_classification_wording() -> None:
    for bundle in (
        _demo_bundle(),
        _uploaded_bundle(),
        _colinear_pca_bundle(),
    ):
        app = _run_pca_page(bundle)
        assert not app.exception
        visible_text = _visible_text(app).lower()
        for prohibited in (
            "confirmed outlier",
            "failed sample",
            "poor-quality sample",
            "pca passed",
            "pca failed",
            "remove this sample",
            "samples are clustered",
            "cluster assignment",
            "hypothesis test performed",
            "confidence interval calculated",
        ):
            assert prohibited not in visible_text


def test_de_page_no_data_state_is_clear_and_creates_no_analysis_state() -> None:
    app = _run_de_page()

    assert not app.exception
    assert app.title[0].value == "🧬 Differential Expression"
    assert "No validated dataset is currently loaded" in _visible_text(app)
    assert CURRENT_DATASET_KEY not in app.session_state
    assert len(app.number_input) == 0
    assert len(app.dataframe) == 0


def test_de_page_without_optional_results_has_truthful_empty_state() -> None:
    bundle = _without_de_bundle()

    app = _run_de_page(bundle)

    assert not app.exception
    visible_text = _visible_text(app)
    assert "No differential-expression results were supplied" in visible_text
    assert "Sample Quality Control, PCA, Sample Correlation, and Gene Expression" in visible_text
    assert len(app.error) == 0
    assert len(app.number_input) == 0
    assert len(app.dataframe) == 0
    assert _download_labels(app) == []
    assert app.session_state[CURRENT_DATASET_KEY] is bundle


def test_de_page_collapses_limitations_into_an_expander() -> None:
    app = _run_de_page(_demo_bundle())

    assert not app.exception
    expander_labels = [element.label for element in app.get("expander")]
    assert "Scientific and statistical limitations" in expander_labels
    visible_text = _visible_text(app)
    assert "does not calculate, adjust, replace, or modify p-values" in visible_text


def test_expression_pages_compute_the_same_results_without_optional_de() -> None:
    with_de = _uploaded_bundle()
    without_de = _without_de_bundle()
    originals = (
        without_de.expression.copy(deep=True),
        without_de.metadata.copy(deep=True),
    )

    page_pairs = (
        (_run_qc_page(with_de), _run_qc_page(without_de)),
        (_run_correlation_page(with_de), _run_correlation_page(without_de)),
        (_run_pca_page(with_de), _run_pca_page(without_de)),
        (
            _run_selected_gene_page(_gene_uploaded_bundle(), "GeneA"),
            _run_selected_gene_page(
                build_dataset_bundle(
                    (
                        _gene_uploaded_bundle().expression,
                        _gene_uploaded_bundle().metadata,
                        None,
                    ),
                    source="uploaded",
                    source_label="Gene expression without DE",
                    report=validate_input_tables(
                        _gene_uploaded_bundle().expression,
                        _gene_uploaded_bundle().metadata,
                        None,
                    ),
                ),
                "GeneA",
            ),
        ),
    )

    for with_de_app, without_de_app in page_pairs:
        assert not with_de_app.exception
        assert not without_de_app.exception
        assert len(with_de_app.dataframe) == len(without_de_app.dataframe)
        for with_table, without_table in zip(
            with_de_app.dataframe,
            without_de_app.dataframe,
            strict=True,
        ):
            pd.testing.assert_frame_equal(
                with_table.value,
                without_table.value,
                check_exact=True,
            )
        assert _download_labels(with_de_app) == _download_labels(without_de_app)

    pd.testing.assert_frame_equal(without_de.expression, originals[0], check_exact=True)
    pd.testing.assert_frame_equal(without_de.metadata, originals[1], check_exact=True)


def test_dataset_context_is_rendered_verbatim_on_every_analysis_page() -> None:
    bundle = _provenance_bundle()
    apps = (
        _run_qc_page(bundle),
        _run_pca_page(bundle),
        _run_correlation_page(bundle),
        _run_de_page(bundle),
        _run_gene_page(bundle),
    )
    expected_values = {
        value
        for value in (
            bundle.provenance.dataset_title,
            bundle.provenance.organism,
            bundle.provenance.expression_scale_description,
            bundle.provenance.upstream_normalization_method,
            bundle.provenance.reference_genome_annotation,
            bundle.provenance.feature_level,
            bundle.provenance.de_contrast_description,
            bundle.provenance.notes,
        )
        if value is not None
    }

    for app in apps:
        assert not app.exception
        rendered_values = {element.value for element in app.get("code")}
        assert expected_values <= rendered_values
        assert "not scientifically verified" in _visible_text(app)


def test_missing_dataset_context_is_explicit_on_every_analysis_page() -> None:
    bundle = _uploaded_bundle()

    for app in (
        _run_qc_page(bundle),
        _run_pca_page(bundle),
        _run_correlation_page(bundle),
        _run_de_page(bundle),
        _run_gene_page(bundle),
    ):
        assert not app.exception
        values = [element.value for element in app.get("code")]
        assert values.count(NOT_SUPPLIED) == 8


def test_de_page_blocks_aggregate_validation_errors_before_controls() -> None:
    bundle = _blocking_de_bundle()

    app = _run_de_page(bundle)

    assert not app.exception
    assert len(app.error) == 2
    assert "blocking validation Errors" in app.error[0].value
    assert "VALUE_OUT_OF_RANGE" in app.error[1].value
    assert len(app.number_input) == 0
    assert len(app.dataframe) == 0
    assert app.session_state[CURRENT_DATASET_KEY] is bundle


def test_de_page_demo_defaults_have_exact_counts_and_views() -> None:
    app = _run_de_page(_demo_bundle())

    assert not app.exception
    metrics = _metric_values(app)
    assert metrics["Positive threshold matches"] == "20"
    assert metrics["Negative threshold matches"] == "20"
    assert metrics["Other evaluable rows"] == "80"
    assert metrics["Not-evaluable rows"] == "0"
    assert len(app.number_input) == 2
    assert app.number_input[0].value == 0.05
    assert app.number_input[1].value == 1.0
    assert len(app.dataframe) == 6
    annotated = app.dataframe[1].value
    assert len(annotated.index) == 120
    assert annotated[STATUS_COLUMN].value_counts().to_dict() == {
        "DOES_NOT_MEET_COMBINED_THRESHOLDS": 80,
        "POSITIVE_THRESHOLD_MATCH": 20,
        "NEGATIVE_THRESHOLD_MATCH": 20,
    }
    assert len(app.dataframe[2].value.index) == 20
    assert len(app.dataframe[3].value.index) == 20
    assert len(app.dataframe[4].value.index) == 80
    assert len(app.dataframe[5].value.index) == 0
    assert "p-values are constructed" in _visible_text(app)


def test_de_page_uploaded_data_preserves_extra_columns_and_shows_every_status() -> None:
    bundle = _de_uploaded_bundle()

    app = _run_de_page(bundle)

    assert not app.exception
    annotated = app.dataframe[1].value
    assert annotated.columns.tolist() == [*bundle.de_results.columns, STATUS_COLUMN]
    assert annotated.index.tolist() == bundle.de_results.index.tolist()
    assert annotated["annotation"].tolist() == ["a", "b", "c", "d", "e"]
    assert annotated[STATUS_COLUMN].tolist() == [
        "POSITIVE_THRESHOLD_MATCH",
        "NEGATIVE_THRESHOLD_MATCH",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
        "NOT_EVALUABLE",
        "POSITIVE_THRESHOLD_MATCH",
    ]
    assert app.dataframe[2].value["gene_id"].tolist() == ["g1", "g5"]
    assert app.dataframe[3].value["gene_id"].tolist() == ["g2"]
    assert app.dataframe[4].value["gene_id"].tolist() == ["g3"]
    assert app.dataframe[5].value["gene_id"].tolist() == ["g4"]
    assert "MISSING_DE_STATISTIC" in _visible_text(app)


def test_de_page_inclusive_boundaries_and_threshold_widget_changes() -> None:
    app = _run_de_page(_de_uploaded_bundle())

    assert app.dataframe[2].value["gene_id"].tolist() == ["g1", "g5"]
    assert app.dataframe[3].value["gene_id"].tolist() == ["g2"]

    app.number_input[0].set_value(0.01).run()

    assert not app.exception
    metrics = _metric_values(app)
    assert metrics["Positive threshold matches"] == "0"
    assert metrics["Negative threshold matches"] == "0"
    assert metrics["Other evaluable rows"] == "4"
    assert metrics["Not-evaluable rows"] == "1"

    app.number_input[0].set_value(0.05).run()
    app.number_input[1].set_value(1.5).run()

    assert not app.exception
    assert app.dataframe[2].value["gene_id"].tolist() == ["g5"]
    assert len(app.dataframe[3].value.index) == 0


def test_de_page_displays_exact_high_precision_thresholds_used_for_classification() -> None:
    app = _run_de_page(_high_precision_de_bundle())

    app.number_input[0].set_value(0.050000001).run()
    app.number_input[1].set_value(1.0000001).run()

    assert not app.exception
    visible_text = _visible_text(app)
    assert "padj ≤ 0.050000001" in visible_text
    assert "log2FoldChange ≥ 1.0000001" in visible_text
    assert "log2FoldChange ≤ -1.0000001" in visible_text
    assert app.dataframe[2].value["gene_id"].tolist() == ["g1"]
    assert app.dataframe[3].value["gene_id"].tolist() == ["g2"]
    assert app.dataframe[4].value["gene_id"].tolist() == ["g3", "g4"]
    assert app.dataframe[5].value["gene_id"].tolist() == ["g5"]


def test_de_page_zero_fold_change_threshold_is_a_controlled_error() -> None:
    app = _run_de_page(_de_uploaded_bundle())

    app.number_input[1].set_value(0.0).run()

    assert not app.exception
    assert any(
        "INVALID_ABSOLUTE_LOG2_FOLD_CHANGE_THRESHOLD" in error.value
        for error in app.error
    )
    assert len(app.dataframe) == 0


def test_de_page_renders_expression_deg_gene_mismatch_warnings() -> None:
    app = _run_de_page(_de_uploaded_bundle(gene_mismatch=True))

    assert not app.exception
    visible_text = _visible_text(app)
    assert "DE_GENE_NOT_IN_EXPRESSION" in visible_text
    assert "EXPRESSION_GENE_NOT_IN_DE" in visible_text
    assert "The rows have been retained" in visible_text
    assert "No expression rows were removed" in visible_text


def test_de_page_repeated_runs_preserve_active_bundle_and_all_input_tables() -> None:
    bundle = _de_uploaded_bundle()
    originals = [
        table.copy(deep=True)
        for table in (bundle.expression, bundle.metadata, bundle.de_results)
    ]
    app = _run_de_page(bundle)

    app.number_input[0].set_value(0.01).run()
    app.number_input[1].set_value(2.0).run()
    app.number_input[0].set_value(0.05).run()

    assert not app.exception
    assert app.session_state[CURRENT_DATASET_KEY] is bundle
    for actual, expected in zip(
        (bundle.expression, bundle.metadata, bundle.de_results),
        originals,
        strict=True,
    ):
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)


def test_de_page_has_required_scientific_wording_and_no_later_features() -> None:
    for bundle in (_demo_bundle(), _de_uploaded_bundle()):
        app = _run_de_page(bundle)
        assert not app.exception
        visible_text = _visible_text(app)
        lower_text = visible_text.lower()
        assert "user-selected exploratory thresholds" in lower_text
        assert "does not fit a differential-expression model" in lower_text
        assert "calculate or modify p-values" in lower_text
        assert "infer the experimental contrast or reference level" in lower_text
        assert "statistical significance or biological importance" in lower_text
        assert "positive and negative labels refer only to the sign" in lower_text
        assert "missing adjusted p-values remain missing" in lower_text
        assert "not a claim of statistical significance" in lower_text
        assert "no volcano plot or gene lookup" in lower_text
        assert "threshold matches are not claims of statistical significance" in lower_text
        assert "upregulated" not in lower_text
        assert "downregulated" not in lower_text
        assert "deseq2 was run" not in lower_text
        assert "fastq processing is available" not in lower_text
        assert len(app.button) == 0
        assert len(app.get("vega_lite_chart")) == 0


def test_gene_page_no_data_state_is_clear_and_has_no_analysis_side_effects() -> None:
    app = _run_gene_page()

    assert not app.exception
    assert app.title[0].value == "🌿 Gene Expression"
    visible_text = _visible_text(app)
    assert "No validated dataset is currently loaded" in visible_text
    assert "Upload Data" in visible_text
    assert len(app.selectbox) == 0
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    assert CURRENT_DATASET_KEY not in app.session_state


def test_gene_page_blocks_aggregate_validation_errors_before_options() -> None:
    bundle = _blocking_de_bundle()

    app = _run_gene_page(bundle)

    assert not app.exception
    assert len(app.error) == 2
    assert "blocking validation Errors" in app.error[0].value
    assert "VALUE_OUT_OF_RANGE" in app.error[1].value
    assert len(app.selectbox) == 0
    assert len(app.dataframe) == 0
    assert app.session_state[CURRENT_DATASET_KEY] is bundle


def test_gene_page_catches_controlled_error_while_building_options() -> None:
    bundle = _invalid_gene_options_bundle()
    original = bundle.expression.copy(deep=True)

    app = _run_gene_page(bundle)

    assert not app.exception
    assert len(app.error) == 1
    assert "Gene options could not be built" in app.error[0].value
    assert "NON_COERCIBLE_EXPRESSION_VALUE" in app.error[0].value
    assert len(app.selectbox) == 0
    assert len(app.dataframe) == 0
    pd.testing.assert_frame_equal(bundle.expression, original, check_exact=True)


def test_gene_page_catches_controlled_error_during_lookup() -> None:
    bundle = _invalid_gene_lookup_bundle()
    app = _run_gene_page(bundle)

    app.selectbox[0].select("GeneA").run()

    assert not app.exception
    assert len(app.error) == 1
    assert "Gene-expression values could not be displayed" in app.error[0].value
    assert "MISSING_REQUIRED_COLUMN" in app.error[0].value
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0


def test_gene_page_numerical_range_error_shows_no_partial_summary() -> None:
    bundle = _numerical_range_gene_bundle()
    original_expression = bundle.expression.copy(deep=True)
    original_metadata = bundle.metadata.copy(deep=True)
    app = _run_gene_page(bundle)

    app.selectbox[0].select("g1").run()

    assert not app.exception
    assert len(app.error) == 1
    assert "NUMERICAL_RANGE_ERROR" in app.error[0].value
    assert "no partial summary was returned" in app.error[0].value
    assert len(app.dataframe) == 0
    assert len(app.get("vega_lite_chart")) == 0
    pd.testing.assert_frame_equal(
        bundle.expression,
        original_expression,
        check_exact=True,
    )
    pd.testing.assert_frame_equal(bundle.metadata, original_metadata, check_exact=True)


def test_gene_page_excludes_expression_deg_coverage_observations() -> None:
    app = _run_gene_page(_de_uploaded_bundle(gene_mismatch=True))

    assert not app.exception
    initial_text = _visible_text(app)
    assert "DE_GENE_NOT_IN_EXPRESSION" not in initial_text
    assert "EXPRESSION_GENE_NOT_IN_DE" not in initial_text
    assert "differential-expression gene identifier(s)" not in initial_text

    app.selectbox[0].select("g1").run()

    assert not app.exception
    visible_text = _visible_text(app)
    assert "DE_GENE_NOT_IN_EXPRESSION" not in visible_text
    assert "EXPRESSION_GENE_NOT_IN_DE" not in visible_text
    assert "differential-expression gene identifier(s)" not in visible_text
    assert len(app.dataframe) == 2


def test_gene_page_uploaded_selection_preserves_values_and_visual_order() -> None:
    bundle = _gene_uploaded_bundle()
    originals = [
        table.copy(deep=True)
        for table in (bundle.expression, bundle.metadata, bundle.de_results)
    ]
    app = _run_gene_page(bundle)

    assert not app.exception
    assert len(app.selectbox) == 1
    assert app.selectbox[0].value is None
    assert len(app.dataframe) == 0

    app.selectbox[0].select("GeneA").run()

    assert not app.exception
    visible_text = _visible_text(app)
    assert "Expression values for GeneA" in visible_text
    assert "differs; conditions were mapped by exact sample ID" in visible_text
    assert len(app.dataframe) == 2
    sample_table = app.dataframe[0].value
    assert sample_table["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
        "sample_c",
    ]
    assert sample_table["condition"].tolist() == [
        "Control",
        "Control",
        "Treated",
    ]
    assert sample_table["expression_value"].tolist() == ["2.0", "1.0", "3.0"]
    condition_summary = app.dataframe[1].value
    assert condition_summary["condition"].tolist() == ["Control", "Treated"]
    assert condition_summary["sample_ids"].tolist() == [
        "sample_b, sample_a",
        "sample_c",
    ]
    assert condition_summary["standard_deviation"].tolist()[1] == "N/A"

    charts = app.get("vega_lite_chart")
    assert len(charts) == 1
    chart_spec = json.loads(charts[0].proto.spec)
    assert chart_spec["mark"]["type"] == "point"
    assert chart_spec["encoding"]["x"]["sort"] == [
        "sample_b",
        "sample_a",
        "sample_c",
    ]
    assert chart_spec["encoding"]["order"] == {
        "field": "sample_position",
        "type": "quantitative",
    }
    assert "line" not in chart_spec["mark"]

    app.selectbox[0].select("genea").run()
    app.selectbox[0].select("GeneA").run()
    assert not app.exception
    assert app.session_state[CURRENT_DATASET_KEY] is bundle
    for actual, expected in zip(
        (bundle.expression, bundle.metadata, bundle.de_results),
        originals,
        strict=True,
    ):
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)


def test_gene_page_dataset_replacement_clears_a_stale_gene_selection() -> None:
    app = _run_gene_page(_gene_uploaded_bundle())
    app.selectbox[0].select("GeneA").run()
    assert len(app.dataframe) == 2

    replacement = _demo_bundle()
    app.session_state[CURRENT_DATASET_KEY] = replacement
    app.run()

    assert not app.exception
    assert app.session_state[CURRENT_DATASET_KEY] is replacement
    assert app.selectbox[0].value is None
    assert app.selectbox[0].options[0] == "SYN_Solyc_0001"
    assert "GeneA" not in app.selectbox[0].options
    assert len(app.dataframe) == 0
    assert "Select one exact supplied gene ID" in _visible_text(app)
    assert _download_labels(app) == []


def test_gene_page_demo_selection_has_exact_options_and_synthetic_disclaimer() -> None:
    app = _run_gene_page(_demo_bundle())

    assert not app.exception
    assert app.selectbox[0].options[0] == "SYN_Solyc_0001"
    assert app.selectbox[0].options[-1] == "SYN_Solyc_0120"
    app.selectbox[0].select("SYN_Solyc_0001").run()

    assert not app.exception
    visible_text = _visible_text(app)
    assert "expression values are synthetic" in visible_text
    assert "gene IDs are fictional" in visible_text
    assert "does not support conclusions about tomato biology" in visible_text
    assert len(app.dataframe[0].value.index) == 6
    assert len(app.get("vega_lite_chart")) == 1


def test_gene_page_collapses_limitations_into_an_expander() -> None:
    app = _run_gene_page(_demo_bundle())
    app.selectbox[0].select("SYN_Solyc_0001").run()

    assert not app.exception
    expander_labels = [element.label for element in app.get("expander")]
    assert "Scientific and statistical limitations" in expander_labels
    visible_text = _visible_text(app)
    assert "Condition labels provide display context only" in visible_text


def test_gene_page_one_sample_has_table_summary_and_truthful_no_chart_state() -> None:
    app = _run_gene_page(_one_sample_gene_bundle())
    app.selectbox[0].select("g1").run()

    assert not app.exception
    visible_text = _visible_text(app)
    assert len(app.dataframe) == 2
    assert app.dataframe[0].value["expression_value"].tolist() == [5.0]
    assert app.dataframe[1].value["standard_deviation"].tolist() == ["N/A"]
    assert "across-sample expression plot is not shown" in visible_text
    assert "Only one supplied condition label" in visible_text
    assert len(app.get("vega_lite_chart")) == 0


def test_gene_page_scientific_disclaimers_are_negative_not_affirmative_claims() -> None:
    app = _run_gene_page(_gene_uploaded_bundle())
    app.selectbox[0].select("GeneA").run()

    assert not app.exception
    visible_text = _visible_text(app).lower()
    for required in (
        "does not test condition effects",
        "does not infer a reference level",
        "no group estimate or statistical comparison is shown",
        "do not establish biological replication, a condition effect, or statistical significance",
        "does not infer a reference level, contrast direction, regulation, a condition effect, statistical significance, or biological importance",
    ):
        assert required in visible_text

    for affirmative_claim in (
        "is statistically significant",
        "shows a condition effect",
        "is upregulated",
        "is downregulated",
        "control is the reference level",
        "is biologically important",
    ):
        assert affirmative_claim not in visible_text


def test_phase_10_downloads_are_absent_without_successful_results() -> None:
    no_data_apps = (
        _run_qc_page(),
        _run_correlation_page(),
        _run_pca_page(),
        _run_de_page(),
        _run_gene_page(),
    )
    controlled_error_apps = (
        _run_qc_page(_invalid_qc_bundle()),
        _run_correlation_page(_uploaded_bundle(one_gene=True)),
        _run_pca_page(_all_genes_constant_bundle()),
        _run_de_page(_blocking_de_bundle()),
        _run_gene_page(_invalid_gene_options_bundle()),
    )

    for app in (*no_data_apps, *controlled_error_apps):
        assert not app.exception
        assert _download_labels(app) == []

    gene_range_error = _run_gene_page(_numerical_range_gene_bundle())
    gene_range_error.selectbox[0].select("g1").run()
    assert not gene_range_error.exception
    assert _download_labels(gene_range_error) == []

    invalid_de_threshold = _run_de_page(_de_uploaded_bundle())
    invalid_de_threshold.number_input[1].set_value(0.0).run()
    assert not invalid_de_threshold.exception
    assert _download_labels(invalid_de_threshold) == []


def test_phase_10_successful_pages_offer_only_contextual_csv_downloads() -> None:
    expected_by_app = (
        (
            _run_qc_page(_demo_bundle()),
            [
                "Download condition and sample membership (CSV)",
                "Download per-sample statistics (CSV)",
            ],
        ),
        (
            _run_correlation_page(_demo_bundle()),
            [
                "Download correlation matrix in long form (CSV)",
                "Download per-sample correlation summary (CSV)",
                "Download unique sample pairs (CSV)",
                "Download condition-pair summary (CSV)",
            ],
        ),
        (
            _run_pca_page(_demo_bundle()),
            [
                "Download explained variance (CSV)",
                "Download sample scores (CSV)",
            ],
        ),
        (
            _run_de_page(_demo_bundle()),
            [
                "Download category counts (CSV)",
                "Download complete annotated results (CSV)",
            ],
        ),
    )

    for app, expected_labels in expected_by_app:
        assert not app.exception
        assert _download_labels(app) == expected_labels
        assert all(
            button.proto.ignore_rerun for button in app.get("download_button")
        )

    gene_page = _run_gene_page(_gene_uploaded_bundle())
    assert _download_labels(gene_page) == []
    gene_page.selectbox[0].select("GeneA").run()
    assert not gene_page.exception
    assert _download_labels(gene_page) == [
        "Download per-sample expression values (CSV)",
        "Download condition-grouped summary (CSV)",
    ]
    assert all(
        button.proto.ignore_rerun for button in gene_page.get("download_button")
    )


def test_phase_10_pages_export_underlying_ordered_tables(monkeypatch) -> None:
    captured: list[tuple[pd.DataFrame, str, tuple[str, ...]]] = []

    def record_export(
        table: pd.DataFrame,
        *,
        filename: str,
        json_sequence_columns: tuple[str, ...] = (),
    ):
        captured.append(
            (table.copy(deep=True), filename, tuple(json_sequence_columns))
        )
        return build_csv_export(
            table,
            filename=filename,
            json_sequence_columns=json_sequence_columns,
        )

    monkeypatch.setattr(exports_module, "build_csv_export", record_export)

    qc_bundle = _uploaded_bundle(one_gene=True)
    qc_app = _run_qc_page(qc_bundle)
    assert not qc_app.exception
    assert [filename for _, filename, _ in captured] == [
        "sample-qc-condition-membership.csv",
        "sample-qc-sample-statistics.csv",
    ]
    assert captured[0][0]["sample_ids"].tolist() == [
        ("sample_a",),
        ("sample_b",),
    ]
    assert captured[0][2] == ("sample_ids",)
    assert captured[1][0]["standard_deviation"].isna().all()

    captured.clear()
    correlation_app = _run_correlation_page(_uploaded_bundle())
    assert not correlation_app.exception
    assert [filename for _, filename, _ in captured] == [
        "sample-correlation-matrix-long.csv",
        "sample-correlation-sample-summary.csv",
        "sample-correlation-unique-pairs.csv",
        "sample-correlation-condition-pairs.csv",
    ]
    assert captured[0][0].columns.tolist() == [
        "row_sample",
        "column_sample",
        "correlation",
        "defined",
    ]
    assert captured[0][0][
        ["row_sample", "column_sample"]
    ].values.tolist() == [
        ["sample_b", "sample_b"],
        ["sample_b", "sample_a"],
        ["sample_a", "sample_b"],
        ["sample_a", "sample_a"],
    ]

    captured.clear()
    pca_app = _run_pca_page(_demo_bundle())
    assert not pca_app.exception
    assert [filename for _, filename, _ in captured] == [
        "pca-explained-variance.csv",
        "pca-sample-scores.csv",
    ]
    assert captured[0][0]["component"].tolist() == [1, 2, 3, 4, 5]
    assert captured[1][0]["sample_id"].tolist() == [
        "Control_1",
        "Control_2",
        "Control_3",
        "High_nitrate_1",
        "High_nitrate_2",
        "High_nitrate_3",
    ]

    captured.clear()
    de_bundle = _de_uploaded_bundle()
    de_app = _run_de_page(de_bundle)
    assert not de_app.exception
    assert [filename for _, filename, _ in captured] == [
        "differential-expression-category-counts.csv",
        "differential-expression-annotated-results.csv",
    ]
    assert captured[1][0].columns.tolist() == [
        *de_bundle.de_results.columns,
        STATUS_COLUMN,
        "applied_adjusted_p_value_threshold",
        "applied_absolute_log2_fold_change_threshold",
    ]
    assert captured[1][0].index.tolist() == de_bundle.de_results.index.tolist()
    for table, _, _ in captured:
        assert table["applied_adjusted_p_value_threshold"].tolist() == [
            0.05
        ] * len(table.index)
        assert table["applied_absolute_log2_fold_change_threshold"].tolist() == [
            1.0
        ] * len(table.index)

    captured.clear()
    de_app.number_input[0].set_value(0.01).run()
    assert not de_app.exception
    assert captured[0][0]["row_count"].tolist() == [1, 0, 0, 4]
    assert captured[1][0][STATUS_COLUMN].tolist() == [
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
        "NOT_EVALUABLE",
        "DOES_NOT_MEET_COMBINED_THRESHOLDS",
    ]
    for table, _, _ in captured:
        assert table["applied_adjusted_p_value_threshold"].tolist() == [
            0.01
        ] * len(table.index)
        assert table["applied_absolute_log2_fold_change_threshold"].tolist() == [
            1.0
        ] * len(table.index)

    captured.clear()
    gene_app = _run_gene_page(_gene_uploaded_bundle())
    assert captured == []
    gene_app.selectbox[0].select("GeneA").run()
    assert not gene_app.exception
    assert [filename for _, filename, _ in captured] == [
        "gene-expression-sample-values.csv",
        "gene-expression-condition-summary.csv",
    ]
    assert captured[0][0]["gene_id"].tolist() == ["GeneA", "GeneA", "GeneA"]
    assert captured[0][0]["sample_id"].tolist() == [
        "sample_b",
        "sample_a",
        "sample_c",
    ]
    assert captured[1][2] == ("sample_ids",)

    captured.clear()
    gene_app.selectbox[0].select("genea").run()
    assert not gene_app.exception
    assert captured[0][0]["gene_id"].tolist() == ["genea", "genea", "genea"]


def test_phase_10_pages_catch_controlled_export_errors(monkeypatch) -> None:
    def fail_export(*args, **kwargs):
        raise CsvExportError(
            CsvExportErrorReason.SERIALIZATION_ERROR,
            "Forced export failure for page regression testing.",
        )

    monkeypatch.setattr(exports_module, "build_csv_export", fail_export)

    apps = (
        _run_qc_page(_uploaded_bundle()),
        _run_correlation_page(_uploaded_bundle()),
        _run_pca_page(_uploaded_bundle()),
        _run_de_page(_de_uploaded_bundle()),
    )
    gene_page = _run_gene_page(_gene_uploaded_bundle())
    gene_page.selectbox[0].select("GeneA").run()

    for app in (*apps, gene_page):
        assert not app.exception
        assert len(app.error) == 1
        assert "SERIALIZATION_ERROR" in app.error[0].value
        assert "No download buttons were shown" in app.error[0].value
        assert _download_labels(app) == []


def test_phase_10_page_download_sets_are_atomic_for_nth_export_failure(
    monkeypatch,
) -> None:
    original_export = build_csv_export
    page_cases = (
        ("QC", 2, lambda: _run_qc_page(_uploaded_bundle())),
        (
            "correlation",
            4,
            lambda: _run_correlation_page(_uploaded_bundle()),
        ),
        ("PCA", 2, lambda: _run_pca_page(_uploaded_bundle())),
        ("DE", 2, lambda: _run_de_page(_de_uploaded_bundle())),
        (
            "gene",
            2,
            lambda: _run_selected_gene_page(_gene_uploaded_bundle(), "GeneA"),
        ),
    )

    for page_name, export_count, run_page in page_cases:
        positions = {1, export_count}
        if export_count > 2:
            positions.add(2)
        for failure_position in sorted(positions):
            call_count = 0

            def fail_nth_export(*args, **kwargs):
                nonlocal call_count
                call_count += 1
                if call_count == failure_position:
                    raise CsvExportError(
                        CsvExportErrorReason.SERIALIZATION_ERROR,
                        f"Forced {page_name} export failure at {failure_position}.",
                    )
                return original_export(*args, **kwargs)

            monkeypatch.setattr(
                exports_module,
                "build_csv_export",
                fail_nth_export,
            )
            app = run_page()

            assert not app.exception
            assert call_count == failure_position
            assert len(app.error) == 1
            assert "SERIALIZATION_ERROR" in app.error[0].value
            assert "No download buttons were shown" in app.error[0].value
            assert _download_labels(app) == []


def test_de_download_csvs_embed_exact_current_thresholds(monkeypatch) -> None:
    captured: dict[str, bytes] = {}
    bundle = _high_precision_de_bundle()
    original = bundle.de_results.copy(deep=True)

    def record_export(
        table: pd.DataFrame,
        *,
        filename: str,
        json_sequence_columns: tuple[str, ...] = (),
    ):
        artifact = build_csv_export(
            table,
            filename=filename,
            json_sequence_columns=json_sequence_columns,
        )
        captured[filename] = artifact.data
        return artifact

    monkeypatch.setattr(exports_module, "build_csv_export", record_export)
    app = _run_de_page(bundle)
    captured.clear()
    adjusted_threshold = 0.010000000000000002
    fold_change_threshold = 1.0000000000000002
    app.number_input[0].set_value(adjusted_threshold).run()
    captured.clear()
    app.number_input[1].set_value(fold_change_threshold).run()

    assert not app.exception
    assert set(captured) == {
        "differential-expression-category-counts.csv",
        "differential-expression-annotated-results.csv",
    }
    for data in captured.values():
        rows = list(csv.DictReader(io.StringIO(data.decode("utf-8"))))
        assert rows
        assert {
            row["applied_adjusted_p_value_threshold"] for row in rows
        } == {repr(adjusted_threshold)}
        assert {
            row["applied_absolute_log2_fold_change_threshold"] for row in rows
        } == {repr(fold_change_threshold)}
    pd.testing.assert_frame_equal(bundle.de_results, original, check_exact=True)


def test_de_export_threshold_column_conflict_is_controlled_and_atomic() -> None:
    bundle = _de_uploaded_bundle()
    de_results = bundle.de_results.copy(deep=True)
    de_results["applied_adjusted_p_value_threshold"] = "supplied"
    original = de_results.copy(deep=True)
    conflicting_bundle = build_dataset_bundle(
        (bundle.expression, bundle.metadata, de_results),
        source="uploaded",
        source_label="DE threshold export column conflict",
        report=bundle.validation_report,
    )

    app = _run_de_page(conflicting_bundle)

    assert not app.exception
    assert len(app.error) == 1
    assert "DUPLICATE_COLUMNS" in app.error[0].value
    assert "No download buttons were shown" in app.error[0].value
    assert _download_labels(app) == []
    pd.testing.assert_frame_equal(de_results, original, check_exact=True)


def test_phase_10_download_wording_discloses_csv_boundaries() -> None:
    apps = (
        _run_qc_page(_uploaded_bundle()),
        _run_correlation_page(_uploaded_bundle()),
        _run_pca_page(_uploaded_bundle()),
        _run_de_page(_de_uploaded_bundle()),
    )
    gene_page = _run_gene_page(_gene_uploaded_bundle())
    gene_page.selectbox[0].select("GeneA").run()

    for app in (*apps, gene_page):
        assert not app.exception
        text = _visible_text(app).lower()
        assert "utf-8 csv" in text
        assert "pandas index" in text
        assert "dtype metadata" in text
        assert "empty fields" in text
        assert "spreadsheet software may interpret" in text
        assert "csv text is not prefixed or rewritten" in text

    de_text = _visible_text(apps[3]).lower()
    assert "threshold matches are not claims of statistical significance" in de_text
    assert "both files also add export-only fields" in de_text
    assert "exact applied adjusted-p-value" in de_text
    gene_text = _visible_text(gene_page).lower()
    assert "gene id is added as an explicit context column" in gene_text


def test_phase_11_navigation_and_documentation_keep_later_work_planned() -> None:
    home = AppTest.from_file("app.py").run()
    page_links = home.get("page_link")
    assert [link.proto.page for link in page_links] == [
        "Upload_Data",
        "Sample_Quality_Control",
        "PCA",
        "Sample_Correlation",
        "Differential_Expression",
        "Gene_Expression",
    ]
    matching = [link for link in page_links if link.proto.page == "Gene_Expression"]
    assert len(matching) == 1
    assert matching[0].proto.label == "Explore Gene Expression"

    repository = Path(__file__).parents[1]
    readme = (repository / "README.md").read_text(encoding="utf-8")
    guide = (repository / "AGENTS.md").read_text(encoding="utf-8")
    assert "Phases 1–11 provide" in readme
    assert "Phase 8 descriptive differential-expression exploration" in readme
    assert "Phase 9 descriptive gene-expression lookup" in readme
    assert "Phase 10 descriptive result exports" in readme
    assert "exact adjusted-p-value and absolute log2-fold-change thresholds" in readme
    assert "Phase 11 GitHub release readiness" in readme
    assert "Phase 12 real-world usability" in readme
    assert "As of Phase 12" in guide
    for text in (readme, guide):
        normalized = " ".join(text.split())
        assert "volcano plots" in normalized
        assert "UTF-8 CSV" in normalized
    assert "no partial set of download buttons" in " ".join(readme.split())
