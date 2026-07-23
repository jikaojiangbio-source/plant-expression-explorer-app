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


def _run_qc_page(bundle: DatasetBundle | None = None) -> AppTest:
    app = AppTest.from_file("pages/2_Sample_Quality_Control.py")
    if bundle is not None:
        app.session_state[CURRENT_DATASET_KEY] = bundle
    return app.run()


def test_home_page_loads_and_describes_scope() -> None:
    app = AppTest.from_file("app.py").run()

    assert not app.exception
    assert app.title[0].value == "🌱 Plant Expression Explorer"
    visible_text = " ".join(
        element.value for element in [*app.subheader, *app.markdown, *app.info]
    )
    assert "preprocessed plant transcriptomics data" in visible_text
    assert "does not process FASTQ files" in visible_text


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
    assert "Analysis features are not implemented yet" in visible_text


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
