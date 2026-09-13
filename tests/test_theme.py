"""Tests for the shared, presentation-only theme helpers."""

from streamlit.testing.v1 import AppTest


def _html_bodies(app: AppTest) -> list[str]:
    return [element.proto.body for element in app.get("html")]


def test_inject_global_styles_renders_a_style_block() -> None:
    app = AppTest.from_string(
        """
import streamlit as st
from plant_expression_explorer.theme import inject_global_styles
inject_global_styles()
"""
    ).run()

    assert not app.exception
    bodies = _html_bodies(app)
    assert len(bodies) == 1
    assert "<style>" in bodies[0]
    assert "--pee-brand" in bodies[0]


def test_render_hero_illustration_renders_an_svg() -> None:
    app = AppTest.from_string(
        """
import streamlit as st
from plant_expression_explorer.theme import render_hero_illustration
render_hero_illustration()
"""
    ).run()

    assert not app.exception
    assert len(app.markdown) == 1
    markdown_value = app.markdown[0].value
    assert markdown_value.strip().startswith("<svg")
    assert "not real data" in markdown_value


def test_render_chip_row_renders_one_chip_per_label() -> None:
    app = AppTest.from_string(
        """
import streamlit as st
from plant_expression_explorer.theme import render_chip_row
render_chip_row(["Alpha", "Beta", "Gamma"])
"""
    ).run()

    assert not app.exception
    bodies = _html_bodies(app)
    assert len(bodies) == 1
    assert bodies[0].count('class="pee-chip"') == 3
    assert "Alpha" in bodies[0]
    assert "Beta" in bodies[0]
    assert "Gamma" in bodies[0]


def test_render_chip_row_handles_an_empty_list() -> None:
    app = AppTest.from_string(
        """
import streamlit as st
from plant_expression_explorer.theme import render_chip_row
render_chip_row([])
"""
    ).run()

    assert not app.exception
    bodies = _html_bodies(app)
    assert len(bodies) == 1
    assert 'class="pee-chip"' not in bodies[0]
