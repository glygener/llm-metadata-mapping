import pandas as pd
from pathlib import Path
from ols_client import Client

# File locations

BASE_DIR = Path(__file__).resolve().parents[2]

XLSX_IN = BASE_DIR / "tissues/data/mapping_BS_Tissue-mapped.xlsx"

XLSX_OUT = BASE_DIR / "tissues/data/tissue_processing_output.xlsx"

# Input columns

INPUT_NAME_COL = "name"
NAMESPACENAME_COL = "namespacename"
NAMESPACEID_COL = "namespaceid"

# OLS client

client = Client("https://www.ebi.ac.uk/ols4/api")

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

for tissue_name in unique_tissues:
    print(f"Looking up: {tissue_name}", flush=True)

    lookup_cache[tissue_name] = lookup_tissue(tissue_name)

    ols_name, ontology_id = lookup_cache[tissue_name]

    print(
        f"Result: {ols_name} | {ontology_id}",
        flush=True
    )
    print()


results = []

for _, row in df.iterrows():

    tissue_name = str(row[INPUT_NAME_COL]).strip()

    ols_name, ontology_id = lookup_cache.get(
        tissue_name,
        (None, None)
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