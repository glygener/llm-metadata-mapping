"""Run the existing disease pipeline on Luis's assigned review rows."""

from pathlib import Path

import disease_processing


DATA_DIR = Path(__file__).resolve().parent.parent / "data"

disease_processing.DATASETS = [
    (
        DATA_DIR / "mapping_Luis_Disease-mapped.xlsx",
        DATA_DIR / "mapping_Luis_Disease_LLM_OLS_results.xlsx",
    )
]

disease_processing.main()
