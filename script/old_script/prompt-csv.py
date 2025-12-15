import pandas as pd
from pathlib import Path

from openai import OpenAI

CSV_IN = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/mapping_BS_GS-mapped.csv")
CSV_OUT = Path("/Users/harivinaygujjula/Documents/GitHub/llm-metadata-mapping/data/updated__BS_GS-mapped10.csv")
PROMPT_FILE = Path("prompt2.txt")
API_FILE = Path("api_key.txt")
SPECIES_COL = "name"
TARGET_COL = "chatgpt_name"
MODEL = "gpt-4o-mini"

api_key = API_FILE.read_text(encoding="utf-8").strip()
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")
df = pd.read_csv(CSV_IN)
df = df.iloc [90:100]

if TARGET_COL not in df.columns:
    df[TARGET_COL] = pd.NA

client = OpenAI(api_key=api_key)

for idx, row in df.iterrows():
    species = str(row.get(SPECIES_COL, "")).strip()
    if not species or species.lower() == "nan":
        continue

    prompt = prompt_template.replace("{species}", species)
    response = client.chat.completions.create(
        model=MODEL,
        messages=[{"role": "user", "content": prompt}],
        temperature = 0.1,
    )
    result = (response.choices[0].message.content or "").strip()
    df.at[idx, TARGET_COL] = result

df.to_csv(CSV_OUT, index=False)
print(f"\nSaved updated file to: {CSV_OUT}")
