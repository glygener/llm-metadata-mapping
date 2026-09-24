import json
import os
import time
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from openai import OpenAI
from ols_client import Client

# File locations

BASE_DIR = Path(__file__).resolve().parents[2]

XLSX_IN = BASE_DIR / "tissues/data/mapping_BS_Tissue-mapped.xlsx"

XLSX_OUT = BASE_DIR / "tissues/data/mapping_BS_Tissue-mapped_processed.xlsx"

# Input columns

INPUT_NAME_COL = "name"
NAMESPACENAME_COL = "namespacename"
NAMESPACEID_COL = "namespaceid"

# OLS client

client = Client("https://www.ebi.ac.uk/ols4/api")

# LLM prompt

PROMPT_FILE = BASE_DIR / "tissues/llm prompts/tissues_llm_prompt.txt"
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")

# OpenAI client

load_dotenv(BASE_DIR / ".env")

api_key = os.getenv("OPENAI_API_KEY")

if not api_key:
    raise ValueError("OPENAI_API_KEY was not found in .env")

llm_client = OpenAI(api_key=api_key)

def translate_tissues(tissue_names):
    prompt = prompt_template.replace(
        "<<TISSUES>>",
        json.dumps(tissue_names, indent=2)
    )

    response = llm_client.chat.completions.create(
        model="gpt-4o",
        messages=[
            {
                "role": "user",
                "content": prompt
            }
        ]
    )

    raw_output = response.choices[0].message.content.strip()

    if raw_output.startswith("```"):
        raw_output = raw_output.strip()
        raw_output = raw_output.removeprefix("```json")
        raw_output = raw_output.removeprefix("```")
        raw_output = raw_output.removesuffix("```")
        raw_output = raw_output.strip()

    records = json.loads(raw_output)

    if not isinstance(records, list):
        raise ValueError("Expected the LLM response to be a JSON array.")

    if len(records) != len(tissue_names):
        raise ValueError(
            f"Expected {len(tissue_names)} results, got {len(records)}"
        )

    required_keys = {"input", "scientific_name"}

    results = {}

    for i, record in enumerate(records):

        if not isinstance(record, dict):
            raise ValueError(
                f"Result {i} is not a JSON object."
            )

        if set(record.keys()) != required_keys:
            raise ValueError(
                "Each LLM result must contain exactly "
                "'input' and 'scientific_name'."
            )

        expected_input = tissue_names[i]
        returned_input = str(record["input"]).strip()

        if expected_input != returned_input:
            raise ValueError(
                f"Input mismatch. Expected '{expected_input}', "
                f"got '{returned_input}'."
            )

        results[returned_input] = record["scientific_name"]

    return results

df = pd.read_excel(XLSX_IN)

print("Loaded spreadsheet.")
print(f"Rows: {len(df)}")
print(f"Columns: {list(df.columns)}")
print()


def lookup_tissue(tissue_name):
    try:
        results = client.search(
            tissue_name,
            query_fields=["label", "synonym"],
            params={
                "ontology": "uberon",
                "exact": "true"
            }
        )

        # Check for exact label match
        for result in results:
            if result.get("ontology_name") != "uberon":
                continue

            label = result.get("label")

            if (
                label is not None
                and label.strip().lower() == tissue_name.strip().lower()
            ):
                return label, result.get("obo_id")

        # Check for exact synonym match
        for result in results:
            if result.get("ontology_name") != "uberon":
                continue

            synonym_fields = [
                "exact_synonyms",
                "related_synonyms",
                "narrow_synonyms",
                "broad_synonyms"
            ]

            for field in synonym_fields:
                synonyms = result.get(field) or []

                for synonym in synonyms:
                    if synonym.strip().lower() == tissue_name.strip().lower():
                        return result.get("label"), result.get("obo_id")

        return None, None

    except Exception as e:
        print(f"OLS lookup failed for {tissue_name}: {e}")
        return None, None


lookup_cache = {}

unique_tissues = (
    df[INPUT_NAME_COL]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)

print(f"Unique tissue terms: {len(unique_tissues)}")
print()

BATCH_SIZE = 10

for start in range(0, len(unique_tissues), BATCH_SIZE):

    batch_names = unique_tissues[start:start + BATCH_SIZE]

    print(
        f"Processing LLM batch {start + 1}-"
        f"{start + len(batch_names)} of {len(unique_tissues)}",
        flush=True
    )

    llm_results = translate_tissues(
        batch_names.tolist()
        if hasattr(batch_names, "tolist")
        else list(batch_names)
    )

    for tissue_name in batch_names:

        scientific_name = llm_results[tissue_name]

        print(
            f"Processing: {tissue_name}",
            flush=True
        )

        print(
            f"LLM scientific name: {scientific_name}",
            flush=True
        )

        if scientific_name == "No Match Found":
            ols_name, ontology_id = None, None

        else:
            ols_name, ontology_id = lookup_tissue(scientific_name)

        lookup_cache[tissue_name] = (
            scientific_name,
            ols_name,
            ontology_id
        )

        print(
            f"OLS result: {ols_name} | {ontology_id}",
            flush=True
        )
        print()

    time.sleep(1)

results = []

for _, row in df.iterrows():

    tissue_name = str(row[INPUT_NAME_COL]).strip()

    scientific_name, ols_name, ontology_id = lookup_cache.get(
        tissue_name,
        (None, None, None)
    )

    original_namespace_id = row[NAMESPACEID_COL]

    if (
        pd.isna(original_namespace_id)
        or str(original_namespace_id).strip() == ""
    ):
        match = ""

    elif (
        ontology_id is not None
        and str(original_namespace_id).strip().upper()
        == str(ontology_id).strip().upper()
    ):
        match = "Yes"

    else:
        match = "No"

    results.append({
        "name": tissue_name,
        "namespacename": row[NAMESPACENAME_COL],
        "namespaceid": row[NAMESPACEID_COL],
        "scientific_name": scientific_name,
        "ols_name": ols_name,
        "ols_ontology_id": ontology_id,
        "namespaceid_match": match
    })


output_df = pd.DataFrame(results)

print("\nFinal results:\n")

print(
    output_df[
        [
            "name",
            "namespacename",
            "namespaceid",
            "ols_name",
            "ols_ontology_id",
            "namespaceid_match"
        ]
    ].to_string(index=False)
)

output_df.to_excel(
    XLSX_OUT,
    index=False
)

print()
print(f"Saved output to: {XLSX_OUT}")