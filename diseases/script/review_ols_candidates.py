from pathlib import Path
import json

import pandas as pd
from dotenv import load_dotenv
from openai import OpenAI

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from diseases.script import disease_processing as dp


ROOT = Path(__file__).resolve().parents[2]
INPUT_FILE = (
    ROOT / "diseases/data/assigned_43_ols_candidates_review.csv"
)
OUTPUT_FILE = (
    ROOT / "diseases/data/assigned_43_gpt_validation_review.csv"
)

BATCH_SIZE = 10

SYSTEM_PROMPT = """
You are assisting with disease ontology mapping validation.

For each record, determine whether the proposed Disease Ontology
candidate represents the same disease concept as the ORIGINAL DISEASE.

The GPT-GENERATED NAME is a suggested normalization, not ground truth.
The OLS CANDIDATE is also not automatically correct.

Choose exactly one assessment:
- equivalent: same disease concept
- broader: candidate represents a broader disease concept
- narrower: candidate represents a more specific disease concept
- related_not_equivalent: related concepts, but not interchangeable
- unrelated: different concepts
- uncertain: insufficient information to decide

Rules:
1. Do not assume similar names mean equivalent diseases.
2. Do not accept a more specific subtype as equivalent to a broad disease.
3. Do not accept a broad disease category as equivalent to a specific subtype.
4. Consider whether the original term may refer to a historical name,
   synonym, disease family, animal model, or phenotype.
5. If the term is ambiguous, choose uncertain.
6. Do not invent supporting facts.
7. Keep explanations concise.

Return ONLY a JSON array. Each object must contain:
candidate_id, assessment, rationale.

Preserve the order of the input records.
"""


def evaluate_batch(client, records):
    prompt = (
        SYSTEM_PROMPT
        + "\n\nEvaluate these records:\n"
        + json.dumps(records, indent=2)
    )

    response = client.chat.completions.create(
        model=dp.OPENAI_MODEL,
        temperature=0,
        messages=[{"role": "user", "content": prompt}],
    )

    raw = response.choices[0].message.content

    if not raw:
        raise ValueError("GPT returned an empty response.")

    results = json.loads(dp.remove_json_fence(raw))

    if not isinstance(results, list) or len(results) != len(records):
        raise ValueError("GPT returned an unexpected number of results.")

    allowed = {
        "equivalent",
        "broader",
        "narrower",
        "related_not_equivalent",
        "unrelated",
        "uncertain",
    }

    for original, result in zip(records, results):
        if result.get("candidate_id") != original["candidate_id"]:
            raise ValueError("Candidate ID mismatch in GPT response.")

        if result.get("assessment") not in allowed:
            raise ValueError("GPT returned an invalid assessment.")

        if not isinstance(result.get("rationale"), str):
            raise ValueError("GPT returned an invalid rationale.")

    return results


def main():
    load_dotenv(ROOT / ".env")

    if not INPUT_FILE.exists():
        raise FileNotFoundError(INPUT_FILE)

    if OUTPUT_FILE.exists():
        raise FileExistsError(
            f"{OUTPUT_FILE} already exists. "
            "Rename or move it before running again."
        )

    df = pd.read_csv(INPUT_FILE)

    records = []

    for _, row in df.iterrows():
        records.append({
            "ID": int(row["ID"]),
            "original_name": row["original_name"],
            "scientific_name": row["scientific_name"],
            "candidate_name": row["candidate_name"],
            "candidate_id": row["candidate_id"],
            "similarity_score": row["similarity_score"],
        })

    client = OpenAI()
    reviewed = []

    for start in range(0, len(records), BATCH_SIZE):
        batch = records[start:start + BATCH_SIZE]
        assessments = evaluate_batch(client, batch)

        for original, assessment in zip(batch, assessments):
            reviewed.append({
                **original,
                **assessment,
                "validation_status": "Requires human verification",
            })

        print(
            f"Reviewed {min(start + len(batch), len(records))}"
            f"/{len(records)} candidates."
        )

    pd.DataFrame(reviewed).to_csv(
        OUTPUT_FILE,
        index=False,
        encoding="utf-8-sig",
    )

    print("\nReview file created:", OUTPUT_FILE)
    print("Total candidates reviewed:", len(reviewed))


if __name__ == "__main__":
    main()
