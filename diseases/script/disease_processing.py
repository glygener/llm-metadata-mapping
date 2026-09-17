from difflib import SequenceMatcher
from pathlib import Path
from time import sleep
import warnings

import pandas as pd
import requests


# Hide the harmless formatting warning from the source spreadsheet.
warnings.filterwarnings(
    "ignore",
    message="Workbook contains no default style.*",
)

# File and service locations.
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR.parent / "data" / "mapping_Disease-mapped.xlsx"
OUTPUT_FILE = (
    SCRIPT_DIR.parent
    / "data"
    / "mapping_Disease_OLS_results.xlsx"
)

OLS_SEARCH_URL = "https://www.ebi.ac.uk/ols4/api/search"
MAX_RETRIES = 3
REQUEST_TIMEOUT = 60


def normalize_text(text):
    """Prepare text for comparisons."""
    return str(text).strip().casefold()


def get_synonyms(document):
    """Return synonyms as a list."""
    synonyms = document.get("synonym") or []

    if isinstance(synonyms, str):
        return [synonyms]

    return synonyms


def similarity_score(term, document):
    """Calculate how similar an OLS result is to the input term."""
    candidate_texts = [
        document.get("label", ""),
        *get_synonyms(document),
    ]

    scores = [
        SequenceMatcher(
            None,
            normalize_text(term),
            normalize_text(candidate),
        ).ratio()
        for candidate in candidate_texts
        if candidate
    ]

    return max(scores, default=0)


def request_ols_documents(term, exact=False):
    """Request Disease Ontology results, retrying temporary failures."""
    parameters = {
        "q": term,
        "ontology": "doid",
        "type": "class",
        "queryFields": "label,synonym",
        "rows": 100,
    }

    if exact:
        parameters["exact"] = "true"

    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.get(
                OLS_SEARCH_URL,
                params=parameters,
                timeout=REQUEST_TIMEOUT,
            )
            response.raise_for_status()

            documents = (
                response.json()
                .get("response", {})
                .get("docs", [])
            )

            return [
                document
                for document in documents
                if str(document.get("obo_id", "")).startswith("DOID:")
            ]

        except requests.RequestException as error:
            last_error = error

            if attempt < MAX_RETRIES:
                print(
                    f"  OLS attempt {attempt}/{MAX_RETRIES} failed. "
                    "Retrying..."
                )
                sleep(2)

    raise last_error


def select_best_document(term, documents):
    """Select the most appropriate document returned by OLS."""
    if not documents:
        return None

    normalized_term = normalize_text(term)

    # Prefer an exact official label.
    for document in documents:
        if normalize_text(document.get("label", "")) == normalized_term:
            return document

    # Then prefer an exact synonym.
    for document in documents:
        if any(
            normalize_text(synonym) == normalized_term
            for synonym in get_synonyms(document)
        ):
            return document

    # Then prefer a label beginning with the complete input term.
    prefix_matches = [
        document
        for document in documents
        if normalize_text(document.get("label", "")).startswith(
            normalized_term + " "
        )
    ]

    if prefix_matches:
        return min(
            prefix_matches,
            key=lambda document: len(
                normalize_text(document.get("label", ""))
            ),
        )

    # Otherwise, select the result with the most similar wording.
    return max(
        documents,
        key=lambda document: similarity_score(term, document),
    )


def search_disease_ontology(term):
    """Search for one term in Disease Ontology."""
    # Try an exact search first.
    exact_documents = request_ols_documents(term, exact=True)

    if exact_documents:
        selected_document = select_best_document(
            term,
            exact_documents,
        )

        return (
            selected_document.get("label"),
            selected_document.get("obo_id"),
        )

    # Use a broader search when no exact result exists.
    documents = request_ols_documents(term, exact=False)
    selected_document = select_best_document(term, documents)

    if selected_document is None:
        return None, None

    return (
        selected_document.get("label"),
        selected_document.get("obo_id"),
    )


# Read the source spreadsheet.
data = pd.read_excel(INPUT_FILE, sheet_name="Mappings")

required_columns = ["name", "namespacename", "namespaceid"]
missing_columns = [
    column for column in required_columns if column not in data.columns
]

if missing_columns:
    raise ValueError(f"Missing required columns: {missing_columns}")

print(f"Loaded {len(data)} disease entries.")

# Process every disease entry.
results = []
total_entries = len(data)

for position, (_, row) in enumerate(data.iterrows(), start=1):
    original_name = row["name"]
    expected_name = row["namespacename"]

    expected_id = (
        ""
        if pd.isna(row["namespaceid"])
        else str(row["namespaceid"]).strip()
    )

    print(f"[{position}/{total_entries}] Searching: {original_name}")

    if pd.isna(original_name) or not str(original_name).strip():
        ols_name = None
        ols_id = None
    else:
        try:
            ols_name, ols_id = search_disease_ontology(
                str(original_name).strip()
            )
        except requests.RequestException as error:
            print(f"  OLS request failed after retries: {error}")
            ols_name = None
            ols_id = None

    id_match = (
        "yes"
        if ols_id
        and normalize_text(ols_id) == normalize_text(expected_id)
        else "no"
    )

    results.append(
        {
            "name": original_name,
            "namespacename": expected_name,
            "namespaceid": expected_id,
            "ols_name": ols_name or "",
            "ols_id": ols_id or "",
            "namespaceid_matches_ols_id": id_match,
        }
    )

# Write the requested output spreadsheet.
results_data = pd.DataFrame(results)

results_data.to_excel(
    OUTPUT_FILE,
    sheet_name="OLS Results",
    index=False,
)

yes_count = (
    results_data["namespaceid_matches_ols_id"] == "yes"
).sum()

print()
print(f"Finished processing {total_entries} disease entries.")
print(f"Matching IDs: {yes_count}/{total_entries}")
print(f"Output file: {OUTPUT_FILE}")