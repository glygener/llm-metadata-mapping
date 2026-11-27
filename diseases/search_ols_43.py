import time
import requests
import pandas as pd
from pathlib import Path

INPUT = Path("diseases/data/mapping_BS_Disease.xlsx")
OUTPUT = Path("diseases/data/ols_direct_search_43.xlsx")

IDS = [
    2487,2494,2497,2501,2502,2503,2505,2506,2507,2508,
    2509,2510,2512,2513,2514,2515,2518,2519,2521,2522,
    2523,2526,2527,2529,2534,2535,2538,2539,2540,2543,
    2544,2547,2548,2549,2550,2554,2555,2557,2558,2559,
    2560,2562,2563
]

API = "https://www.ebi.ac.uk/ols4/api/search"

def normalize(value):
    return " ".join(str(value).strip().casefold().split())

def search_ols(term, exact=False):
    params = {
        "q": term,
        "ontology": "doid",
        "queryFields": "label,synonym",
        "fieldList": "iri,label,short_form,obo_id,ontology_name,ontology_prefix,description,synonym",
        "rows": 10,
        "exact": str(exact).lower()
    }

    response = requests.get(API, params=params, timeout=30)
    response.raise_for_status()
    return response.json().get("response", {}).get("docs", [])

df = pd.read_excel(INPUT)
df["ID"] = pd.to_numeric(df["ID"], errors="coerce")
mine = df[df["ID"].isin(IDS)].copy()

if len(mine) != 43:
    raise ValueError(f"Expected 43 rows, found {len(mine)}")

summary = []
all_candidates = []

for _, row in mine.iterrows():
    disease_id = int(row["ID"])
    name = str(row["name"]).strip()

    print(f"Searching {disease_id}: {name}")

    try:
        exact_results = search_ols(name, exact=True)
        time.sleep(0.15)

        broad_results = search_ols(name, exact=False)
        time.sleep(0.15)

        # Combine results and remove duplicate ontology terms.
        combined = {}
        for result in exact_results + broad_results:
            key = result.get("iri") or result.get("obo_id")
            if key:
                combined[key] = result

        candidates = list(combined.values())

        for rank, result in enumerate(candidates, start=1):
            label = result.get("label", "")
            synonyms = result.get("synonym", []) or []

            if isinstance(synonyms, str):
                synonyms = [synonyms]

            if normalize(label) == normalize(name):
                match_type = "Exact label"
            elif any(normalize(s) == normalize(name) for s in synonyms):
                match_type = "Exact synonym"
            else:
                match_type = "Search candidate"

            candidate = {
                "ID": disease_id,
                "original_name": name,
                "ols_name": label,
                "ols_id": result.get("obo_id", ""),
                "ols_iri": result.get("iri", ""),
                "description": result.get("description", ""),
                "synonyms": "; ".join(synonyms),
                "match_type": match_type,
                "search_rank": rank
            }
            all_candidates.append(candidate)

        exact_match = next(
            (c for c in all_candidates
             if c["ID"] == disease_id
             and c["match_type"] in ["Exact label", "Exact synonym"]),
            None
        )

        top = candidates[0] if candidates else {}

        summary.append({
            "ID": disease_id,
            "name": name,
            "ols_name": exact_match["ols_name"] if exact_match else top.get("label", ""),
            "ols_id": exact_match["ols_id"] if exact_match else top.get("obo_id", ""),
            "match_type": exact_match["match_type"] if exact_match else (
                "Candidate only" if candidates else "No results"
            ),
            "candidate_count": len(candidates),
            "review_required": "Yes"
        })

    except Exception as error:
        print(f"  ERROR: {error}")
        summary.append({
            "ID": disease_id,
            "name": name,
            "ols_name": "",
            "ols_id": "",
            "match_type": "API error",
            "candidate_count": 0,
            "review_required": "Yes"
        })

summary_df = pd.DataFrame(summary)
candidates_df = pd.DataFrame(all_candidates)

with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    candidates_df.to_excel(writer, sheet_name="All candidates", index=False)

print("\nFinished!")
print("Rows searched:", len(summary_df))
print("Output:", OUTPUT)
print("\nSummary:")
print(summary_df.to_string(index=False))
