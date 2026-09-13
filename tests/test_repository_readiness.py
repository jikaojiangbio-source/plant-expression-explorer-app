"""Regression checks for the committed GitHub release configuration."""

from pathlib import Path


REPOSITORY = Path(__file__).parents[1]
WORKFLOW = REPOSITORY / ".github" / "workflows" / "tests.yml"


def test_python_and_dependency_versions_are_pinned() -> None:
    assert (REPOSITORY / ".python-version").read_text(encoding="utf-8") == "3.13.9\n"
    requirements = (REPOSITORY / "requirements.txt").read_text(encoding="utf-8")
    assert requirements.splitlines() == [
        "streamlit==1.60.0",
        "pandas==2.3.3",
        "numpy==2.5.3",
        "pytest==9.1.1",
    ]


def test_ci_uses_pinned_actions_and_read_only_permissions() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")

    assert "permissions:\n  contents: read\n" in workflow
    assert "pull_request_target:" not in workflow
    assert "uses: actions/checkout@3d3c42e5aac5ba805825da76410c181273ba90b1" in workflow
    assert "uses: actions/setup-python@5fda3b95a4ea91299a34e894583c3862153e4b97" in workflow
    assert "persist-credentials: false" in workflow
    assert "python-version-file: .python-version" in workflow
    assert "runs-on: ubuntu-24.04" in workflow


def test_ci_runs_dependency_and_complete_test_checks() -> None:
    workflow = WORKFLOW.read_text(encoding="utf-8")

    expected_commands = (
        "python -m pip install --upgrade pip==26.2.1",
        "python -m pip install -r requirements.txt",
        "python -m pip check",
        "python -m pytest",
    )
    positions = [workflow.index(command) for command in expected_commands]
    assert positions == sorted(positions)


def test_repository_guidance_preserves_scientific_boundaries() -> None:
    contributing = (REPOSITORY / "CONTRIBUTING.md").read_text(encoding="utf-8")
    normalized = " ".join(contributing.split()).lower()

    required_boundaries = (
        "processes fastq files",
        "normalizes or transforms raw counts",
        "corrects batch effects",
        "fits a differential-expression model",
        "source dataframes remain unchanged",
        "never commit private, unpublished, or identifiable biological data",
    )
    for boundary in required_boundaries:
        assert boundary in normalized
