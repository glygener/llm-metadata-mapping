from difflib import SequenceMatcher
from pathlib import Path
from time import sleep
import json
import os
import warnings

import pandas as pd
import requests
from openai import OpenAI


warnings.filterwarnings(
    "ignore",
    message="Workbook contains no default style.*",
)

SCRIPT_DIR = Path(__file__).resolve().parent
DISEASE_DIR = SCRIPT_DIR.parent
DATA_DIR = DISEASE_DIR / "data"
PROMPT_FILE = (
    DISEASE_DIR
    / "LLM prompts"
    / "diseases_LLM_prompt.txt"
)

DATASETS = [
    (
        DATA_DIR / "mapping_Disease-mapped.xlsx",
        DATA_DIR / "mapping_Disease_LLM_OLS_results.xlsx",
    ),
    (
        DATA_DIR / "mapping_BS_Disease-mapped.xlsx",
        DATA_DIR / "mapping_BS_Disease_LLM_OLS_results.xlsx",
    ),
]

OLS_SEARCH_URL = "https://www.ebi.ac.uk/ols4/api/search"
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o")
BATCH_SIZE = 10
MAX_RETRIES = 3
REQUEST_TIMEOUT = 60
REQUIRED_COLUMNS = ["name", "namespacename", "namespaceid"]


def normalize_text(text):
    """Prepare text for comparisons."""
    return str(text).strip().casefold()


def get_synonyms(document):
    """Return synonyms as a list."""
    synonyms = document.get("synonym") or []

    if isinstance(synonyms, str):
        return [synonyms]

    return synonyms


def create_openai_client():
    """Create the OpenAI client without storing a key in the repository."""
    api_key = os.getenv("OPENAI_API_KEY")

    if not api_key:
        raise RuntimeError(
            "OPENAI_API_KEY is not set. Add it to the current PowerShell "
            "session before running this script."
        )

    return OpenAI(api_key=api_key)


def remove_json_fence(raw_output):
    """Remove optional Markdown fences from an LLM JSON response."""
    cleaned_output = raw_output.strip()

    if cleaned_output.startswith("```"):
        cleaned_output = cleaned_output.removeprefix("```json")
        cleaned_output = cleaned_output.removeprefix("```")
        cleaned_output = cleaned_output.removesuffix("```")
        cleaned_output = cleaned_output.strip()

    return cleaned_output


def translate_diseases(disease_names, client, prompt_template):
    """Translate a batch of input names to standard disease terms."""
    prompt = prompt_template.replace(
        "<<DISEASES>>",
        json.dumps(disease_names, indent=2),
    )

    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = client.chat.completions.create(
                model=OPENAI_MODEL,
                temperature=0,
                messages=[{"role": "user", "content": prompt}],
            )

            raw_output = response.choices[0].message.content

            if not raw_output:
                raise ValueError(
                    "The ChatGPT API returned an empty response."
                )

            records = json.loads(remove_json_fence(raw_output))

            if not isinstance(records, list):
                raise ValueError(
                    "Expected the LLM response to be a JSON array."
                )

            if len(records) != len(disease_names):
                raise ValueError(
                    f"Expected {len(disease_names)} results, "
                    f"got {len(records)}."
                )

            translated_names = {}
            required_keys = {"input", "scientific_name"}

            for position, record in enumerate(records):
                if not isinstance(record, dict):
                    raise ValueError(
                        f"Result {position} is not a JSON object."
                    )

                if set(record.keys()) != required_keys:
                    raise ValueError(
                        "Every LLM result must contain exactly "
                        "'input' and 'scientific_name'."
                    )

                expected_input = disease_names[position]
                returned_input = str(record["input"]).strip()

                if returned_input != expected_input:
                    raise ValueError(
                        f"Input mismatch. Expected '{expected_input}', "
                        f"got '{returned_input}'."
                    )

                translated_names[expected_input] = str(
                    record["scientific_name"]
                ).strip()

            return translated_names

        except Exception as error:
            last_error = error

            if attempt < MAX_RETRIES:
                print(
                    f"  ChatGPT attempt {attempt}/{MAX_RETRIES} failed. "
                    "Retrying..."
                )
                sleep(2)

    raise last_error


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

    for document in documents:
        if normalize_text(document.get("label", "")) == normalized_term:
            return document

    for document in documents:
        if any(
            normalize_text(synonym) == normalized_term
            for synonym in get_synonyms(document)
        ):
            return document

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

    return max(
        documents,
        key=lambda document: similarity_score(term, document),
    )


def search_disease_ontology(term):
    """Search for one standard term in Disease Ontology."""
    exact_documents = request_ols_documents(term, exact=True)

    if exact_documents:
        selected_document = select_best_document(term, exact_documents)
        return (
            selected_document.get("label"),
            selected_document.get("obo_id"),
        )

    documents = request_ols_documents(term, exact=False)
    selected_document = select_best_document(term, documents)

    if selected_document is None:
        return None, None

    return (
        selected_document.get("label"),
        selected_document.get("obo_id"),
    )


def load_dataset(input_file):
    """Load and validate one mapped disease spreadsheet."""
    data = pd.read_excel(input_file, sheet_name="Mappings")
    missing_columns = [
        column
        for column in REQUIRED_COLUMNS
        if column not in data.columns
    ]

    if missing_columns:
        raise ValueError(
            f"{input_file.name} is missing columns: {missing_columns}"
        )

    print(f"Loaded {input_file.name}: {len(data)} entries.")
    return data


def collect_unique_terms(datasets):
    """Collect unique non-empty names from all input spreadsheets."""
    unique_terms = []
    seen_terms = set()

    for _, _, data in datasets:
        for value in data["name"]:
            if pd.isna(value) or not str(value).strip():
                continue

            term = str(value).strip()

            if term not in seen_terms:
                seen_terms.add(term)
                unique_terms.append(term)

    return unique_terms


def build_lookup_cache(unique_terms, client, prompt_template):
    """Translate and validate every unique term only once."""
    lookup_cache = {}

    for start in range(0, len(unique_terms), BATCH_SIZE):
        batch_names = unique_terms[start:start + BATCH_SIZE]
        batch_end = start + len(batch_names)

        print(
            f"Processing LLM batch {start + 1}-{batch_end} "
            f"of {len(unique_terms)}"
        )

        llm_results = translate_diseases(
            batch_names,
            client,
            prompt_template,
        )

        for original_name in batch_names:
            scientific_name = llm_results[original_name]
            print(f"  Input: {original_name}")
            print(f"  LLM term: {scientific_name}")

            if normalize_text(scientific_name) == "no match found":
                ols_name = None
                ols_id = None
            else:
                try:
                    ols_name, ols_id = search_disease_ontology(
                        scientific_name
                    )
                except requests.RequestException as error:
                    print(f"  OLS request failed: {error}")
                    ols_name = None
                    ols_id = None

            lookup_cache[original_name] = (
                scientific_name,
                ols_name,
                ols_id,
            )
            print(f"  OLS result: {ols_name} | {ols_id}")

        sleep(1)

    return lookup_cache


def create_output(data, lookup_cache):
    """Create the output rows for one input spreadsheet."""
    results = []

    for _, row in data.iterrows():
        original_value = row["name"]
        original_name = (
            ""
            if pd.isna(original_value)
            else str(original_value).strip()
        )

        scientific_name, ols_name, ols_id = lookup_cache.get(
            original_name,
            (None, None, None),
        )

        expected_id = (
            ""
            if pd.isna(row["namespaceid"])
            else str(row["namespaceid"]).strip()
        )

        if not expected_id:
            id_match = ""
        elif (
            ols_id
            and normalize_text(ols_id) == normalize_text(expected_id)
        ):
            id_match = "yes"
        else:
            id_match = "no"

        results.append(
            {
                "name": original_value,
                "namespacename": row["namespacename"],
                "namespaceid": expected_id,
                "scientific_name": scientific_name or "",
                "ols_name": ols_name or "",
                "ols_id": ols_id or "",
                "namespaceid_matches_ols_id": id_match,
            }
        )

    return pd.DataFrame(results)


def main():
    """Run the LLM and OLS workflow on both mapped disease files."""
    prompt_template = PROMPT_FILE.read_text(encoding="utf-8")
    client = create_openai_client()
    datasets = []

    for input_file, output_file in DATASETS:
        data = load_dataset(input_file)
        datasets.append((input_file, output_file, data))

    unique_terms = collect_unique_terms(datasets)
    print(f"Unique disease terms across both files: {len(unique_terms)}")

    lookup_cache = build_lookup_cache(
        unique_terms,
        client,
        prompt_template,
    )

    for input_file, output_file, data in datasets:
        output_data = create_output(data, lookup_cache)
        output_data.to_excel(
            output_file,
            sheet_name="LLM OLS Results",
            index=False,
        )

        comparable_rows = (
            output_data["namespaceid_matches_ols_id"] != ""
        ).sum()
        yes_count = (
            output_data["namespaceid_matches_ols_id"] == "yes"
        ).sum()

        print()
        print(f"Finished {input_file.name}.")
        print(f"Matching IDs: {yes_count}/{comparable_rows}")
        print(f"Output file: {output_file}")


if __name__ == "__main__":
    main()
