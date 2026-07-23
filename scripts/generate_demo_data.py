"""Generate the bundled synthetic demonstration CSV files."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from plant_expression_explorer.demo_data import DEMO_SEED, write_demo_data

DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "demo"


def main() -> None:
    """Parse command-line options and write the synthetic demo tables."""

    parser = argparse.ArgumentParser(
        description="Generate deterministic synthetic transcriptomics demo data."
    )
    parser.add_argument("--seed", type=int, default=DEMO_SEED)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the three generated CSV files.",
    )
    arguments = parser.parse_args()

    paths = write_demo_data(arguments.output_dir, seed=arguments.seed)
    print("Wrote synthetic demo data:")
    for path in paths:
        print(f"- {path}")


if __name__ == "__main__":
    main()
