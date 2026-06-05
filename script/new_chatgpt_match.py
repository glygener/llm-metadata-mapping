import pandas as pd
from pathlib import Path
import json

from openai import OpenAI

#file locations
CSV_IN = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/mapping_BS_GS-mapped.csv")
CSV_OUT = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/updated_BS_GS-mapped_with_reasoning2.csv")
PROMPT_FILE = Path("modified_prompt.txt")
API_FILE = Path("api_key.txt")

#input column name
SPECIES_COL = "name"
#output column name
CHATGPT_COL = "chatgpt_name"
#output column name
TAXON_COL = "taxon_id"
#output column name
REASON_COL = "reasoning"

#API model name
MODEL = "gpt-4o-mini"

#reading files
api_key = API_FILE.read_text(encoding="utf-8").strip()
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")

#loads csv into pandas
df = pd.read_csv(CSV_IN)

#skips row 0 and only keeps 1-99
df = df.iloc[1:100]

#makes sure the output columns exist
for col in [CHATGPT_COL, TAXON_COL, REASON_COL]:
 #if the column doesn't exist yet, it creates it
    if col not in df.columns:
        df[col] = pd.NA

client = OpenAI(api_key=api_key)

#processes each row
for idx, row in df.iterrows():
 #gets the species name from "name" (SPECIES_COL)
    species = str(row.get(SPECIES_COL, "")).strip()

#skips empty rows
    if not species or species.lower() == "nan":
        continue

#skips already processed rows
    if pd.notna(row.get(CHATGPT_COL)) and str(row.get(CHATGPT_COL)).strip():
        print(f"Row {idx}: already processed, skipping.")
        continue

#inserts "name" into the prompt
    prompt = prompt_template.replace("<<SPECIES>>", species)
    try:

#puts into API
        response = client.chat.completions.create(
            model=MODEL,
            messages=[
                {"role": "user", "content": prompt}
            ],
            temperature=0
        )

#gets API response and prints it
        raw_output = response.choices[0].message.content.strip()
        print(f"\nRow {idx} | Input species: {species}")
        print("Raw model output:")
        print(raw_output)

#converts JSON text into python
        data = json.loads(raw_output)

#not sure
        if isinstance(data, list) and len(data) > 0:
            record = data[0]
        elif isinstance(data, dict):
            record = data
        else:
            print(f"Row {idx}: Unexpected JSON structure, skipping.")
            continue

#gets chatgpt_taxon_id
        taxon_id = record.get("taxon_id", "")
#gets chatgpt_name
        species_name = record.get("species_name", "")
#gets chatgpt_reasoning
        reasoning = record.get("reasoning", "")

#stores results
        df.at[idx, CHATGPT_COL] = species_name
        df.at[idx, TAXON_COL] = taxon_id
        df.at[idx, REASON_COL] = reasoning

        print(f"Parsed -> taxon_id: {taxon_id}, species_name: {species_name}, reasoning: {reasoning}")

#saves progress after every row
        df.to_csv(CSV_OUT, index=False)
        print(f"Row {idx}: saved to {CSV_OUT}")

#catches JSON errors
    except json.JSONDecodeError as e:
        print(f"Row {idx}: JSON decode error: {e}")
        print("Raw output that failed to parse:")
#shows what went wrong
        print(raw_output)
#moves to next row
        continue

#catches every other error
    except Exception as e:
#shows what went wrong
        print(f"Row {idx}: Error calling API or processing result: {e}")
#moves to next row
        continue

#prints "Done." when everything is finished
print("Done.")