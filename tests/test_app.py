"""Smoke test for the Streamlit home page."""

import pandas as pd
from streamlit.testing.v1 import AppTest

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
from plant_expression_explorer.validation import ValidationReport


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


def _demo_bundle() -> DatasetBundle:
    candidate = load_demo_candidate()
    assert candidate.tables is not None
    return build_dataset_bundle(
        candidate.tables,
        source="demo",
        source_label=DEMO_SOURCE_LABEL,
        report=candidate.report,
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


def test_home_page_analysis_link_order_is_upload_qc_pca_correlation() -> None:
    app = AppTest.from_file("app.py").run()

    page_links = app.get("page_link")
    pages_in_order = [link.proto.page for link in page_links]

    assert pages_in_order == [
        "Upload_Data",
        "Sample_Quality_Control",
        "PCA",
        "Sample_Correlation",
    ]


def test_home_page_deg_and_gene_expression_remain_planned() -> None:
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
    assert "planned" in deg_line.lower()
    assert "planned" in gene_expression_line.lower()

    info_text = " ".join(element.value for element in app.info).lower()
    assert "differential-expression" in info_text
    assert "gene lookup" in info_text
    assert "dedicated exports" in info_text


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
    assert "Sample Quality Control, PCA, and Sample Correlation are available" in visible_text
    assert "not implemented yet" in visible_text


def test_upload_page_no_longer_lists_implemented_features_as_unimplemented() -> None:
    app = AppTest.from_file("pages/1_Upload_Data.py").run()

    assert not app.exception
    visible_text = _visible_text(app)

    for unimplemented_phrase in (
        "Sample QC, PCA, sample correlation",
        "Analysis features are not implemented yet",
    ):
        assert unimplemented_phrase not in visible_text

    assert "Sample Quality Control, PCA, and Sample Correlation are available" in visible_text
    assert "DEG filtering, significance classification, volcano plots, gene " in visible_text
    assert "lookup, and exports are not implemented yet" in visible_text


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
    visible_text = _visible_text(app)
    assert "A validated dataset is active" in visible_text
    assert "Bundled synthetic demonstration data" in visible_text
    assert "120 genes, 6 samples" in visible_text

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
