import pandas as pd
from pathlib import Path
import json

from openai import OpenAI



CSV_IN = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/mapping_BS_GS-mapped.csv")
CSV_OUT = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/updated_BS_GS-mapped_with_reasoning2.csv")
PROMPT_FILE = Path("modified_prompt.txt")
API_FILE = Path("api_key.txt")

SPECIES_COL = "name"
CHATGPT_COL = "chatgpt_name"
TAXON_COL = "taxon_id"
REASON_COL = "reasoning"

MODEL = "gpt-4o-mini"



api_key = API_FILE.read_text(encoding="utf-8").strip()
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")

df = pd.read_csv(CSV_IN)
df = df.iloc[1:100]


for col in [CHATGPT_COL, TAXON_COL, REASON_COL]:
    if col not in df.columns:
        df[col] = pd.NA

client = OpenAI(api_key=api_key)



for idx, row in df.iterrows():
    species = str(row.get(SPECIES_COL, "")).strip()


    if not species or species.lower() == "nan":
        continue


    if pd.notna(row.get(CHATGPT_COL)) and str(row.get(CHATGPT_COL)).strip():
        print(f"Row {idx}: already processed, skipping.")
        continue

    prompt = prompt_template.replace("<<SPECIES>>", species)
    try:

        response = client.chat.completions.create(
            model=MODEL,
            messages=[
                {"role": "user", "content": prompt}
            ],
            temperature=0
        )

        raw_output = response.choices[0].message.content.strip()
        print(f"\nRow {idx} | Input species: {species}")
        print("Raw model output:")
        print(raw_output)



        data = json.loads(raw_output)


        if isinstance(data, list) and len(data) > 0:
            record = data[0]
        elif isinstance(data, dict):
            record = data
        else:
            print(f"Row {idx}: Unexpected JSON structure, skipping.")
            continue

        taxon_id = record.get("taxon_id", "")
        species_name = record.get("species_name", "")
        reasoning = record.get("reasoning", "")


        df.at[idx, CHATGPT_COL] = species_name
        df.at[idx, TAXON_COL] = taxon_id
        df.at[idx, REASON_COL] = reasoning

        print(f"Parsed -> taxon_id: {taxon_id}, species_name: {species_name}, reasoning: {reasoning}")


        df.to_csv(CSV_OUT, index=False)
        print(f"Row {idx}: saved to {CSV_OUT}")

    except json.JSONDecodeError as e:
        print(f"Row {idx}: JSON decode error: {e}")
        print("Raw output that failed to parse:")
        print(raw_output)

        continue

    except Exception as e:
        print(f"Row {idx}: Error calling API or processing result: {e}")
        continue

print("Done.")
