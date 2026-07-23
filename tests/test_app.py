"""Smoke test for the Streamlit home page."""

from streamlit.testing.v1 import AppTest

from plant_expression_explorer.dataset import (
    CURRENT_DATASET_KEY,
    DE_RESULTS_UPLOAD_KEY,
    EXPRESSION_UPLOAD_KEY,
    METADATA_UPLOAD_KEY,
    DatasetBundle,
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
