import sys
import json
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from disease_processing import (
    request_ols_documents,
    similarity_score,
)

SOURCE = Path("diseases/data/assigned_43_spreadsheet_ready.csv")
OUTPUT = Path("diseases/data/assigned_43_broad_candidates.csv")

df = pd.read_csv(SOURCE)
candidates = {}
errors = []

for _, row in df.iterrows():
    disease_id = int(row["ID"])
    original = str(row["name"]).strip()
    gpt_name = str(row["chatgpt_name"]).strip()

    queries = [("original", original)]

    if gpt_name.casefold() != "no match found":
        if gpt_name.casefold() != original.casefold():
            queries.append(("chatgpt", gpt_name))

    for source, term in queries:
        print(f"Searching ID {disease_id}: {term}")

        try:
            documents = request_ols_documents(term, exact=False)

            ranked = sorted(
                documents,
                key=lambda doc: similarity_score(term, doc),
                reverse=True,
            )

            for document in ranked[:5]:
                ols_id = document.get("obo_id")

                if not ols_id:
                    continue

                key = (disease_id, ols_id)
                score = similarity_score(term, document)

                synonyms = document.get("synonym") or []
                if isinstance(synonyms, str):
                    synonyms = [synonyms]

                candidate = {
                    "ID": disease_id,
                    "original_name": original,
                    "chatgpt_name": gpt_name,
                    "ols_name": document.get("label"),
                    "ols_id": ols_id,
                    "similarity_score": score,
                    "synonyms": json.dumps(synonyms),
                    "search_queries": set(),
                }

                candidate["search_queries"].add(
                    f"{source}: {term}"
                )

                if key not in candidates:
                    candidates[key] = candidate
                else:
                    existing = candidates[key]
                    existing["search_queries"].update(
                        candidate["search_queries"]
                    )

                    if score > existing["similarity_score"]:
                        candidate["search_queries"] = (
                            existing["search_queries"]
                        )
                        candidates[key] = candidate

        except Exception as error:
            errors.append({
                "ID": disease_id,
                "query": term,
                "error": str(error),
            })
            print(f"  Search failed: {error}")

        time.sleep(0.2)

results = []

for candidate in candidates.values():
    candidate["search_queries"] = "; ".join(
        sorted(candidate["search_queries"])
    )
    results.append(candidate)

result_df = pd.DataFrame(results)

if not result_df.empty:
    result_df = result_df.sort_values(
        ["ID", "similarity_score"],
        ascending=[True, False],
    )

result_df.to_csv(OUTPUT, index=False, encoding="utf-8-sig")

print("\nSEARCH SUMMARY")
print("Diseases searched:", len(df))
print("Unique candidates:", len(result_df))
print("Search errors:", len(errors))
print("Diseases with candidates:", result_df["ID"].nunique() if not result_df.empty else 0)
print("Saved:", OUTPUT)

if errors:
    pd.DataFrame(errors).to_csv(
        "diseases/data/assigned_43_search_errors.csv",
        index=False,
    )
    print("Errors saved to diseases/data/assigned_43_search_errors.csv")

print("\nTOP CANDIDATES")
if not result_df.empty:
    print(
        result_df[
            ["ID", "original_name", "ols_name", "ols_id", "similarity_score"]
        ].to_string(index=False)
    )
