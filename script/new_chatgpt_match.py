import pandas as pd
from pathlib import Path
import json
from openai import OpenAI

#tells the script where the prompt file is located
BASE_DIR = Path(__file__).parent

#file locations
CSV_IN = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\list_of_species.csv")
CSV_OUT = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\mapped_with_ollama.csv")
PROMPT_FILE = BASE_DIR/"modified_prompt.txt"

#API_FILE = Path("api_key.txt")
#will use later when we get api key

#input column name
INPUT_NAME_COL = "name"
#input column name
NAMESPACENAME_COL = "namespacename"
#input column name
NAMESPACEID_COL = "namespaceid"
#output column name
CHATGPT_COL = "chatgpt_name"
#output column name
NAME_MATCH_COL = "name_match_?"
#output column name
TAXON_COL = "chatgpt_taxon_id"
#output column name
TAXON_MATCH_COL = "taxon_id_match_?"
#output column name
REASON_COL = "chatgpt_reasoning"

#API model name, for now we're using ollama
MODEL = "llama3.2"

#reading files
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")

#api_key = API_FILE.read_text(encoding="utf-8").strip()
#will use later when we get api key

#loads csv into pandas
df = pd.read_csv(CSV_IN)

if CHATGPT_COL in df.columns:
    df[CHATGPT_COL] = df[CHATGPT_COL].astype("object")

#skips row 0 and only keeps 1-9 for testing ollama
df = df.iloc[1:10]

#makes sure the output columns exist
for col in [CHATGPT_COL, NAME_MATCH_COL, TAXON_COL, TAXON_MATCH_COL, REASON_COL]:
 #if the column doesn't exist yet, it creates it
    if col not in df.columns:
        df[col] = pd.NA

client = OpenAI(
     api_key="ollama",
     base_url="http://localhost:11434/v1"
)

#client = OpenAI(api_key=api_key)
#will use later after we get the api key

#processes each row
for idx, row in df.iterrows():
    
 #gets the species name from "name" (INPUT_NAME_COL)
    species = str(row.get(INPUT_NAME_COL, "")).strip()

    #skips empty rows
    if not species or species.lower() == "nan":
        continue
    
    #skips already processed rows
    if pd.notna(row.get(CHATGPT_COL)) and str(row.get(CHATGPT_COL)).strip():
        print(f"Row {idx}: already processed, skipping.")
        continue
    
    #inserts "name" into the prompt
    prompt = prompt_template.replace("<<SPECIES>>", species)

#calls the API with the prompt and gets a response
    try:
         response = client.chat.completions.create(
              model=MODEL,
              messages=[
            {
                 #tells the model to only return valid JSON, with no explanations or extra text
                 "role": "user",
                 "content": (
                      "Return ONLY valid JSON."
                      "No explanations."
                      "No text before or after the JSON."
             )
            },
            {
                 #gives the prompt with the species name
                 "role": "user",
                 "content": prompt}
            ],
            #gets the most likely answer (not a random one)
            temperature=0
            )
         
         #gets API response and prints it
         raw_output = response.choices[0].message.content.strip()

         print(f"\nRow {idx} | Input species: {species}")
         print("Raw model output:")
         print(raw_output)

         #converts JSON text into python
         data = json.loads(raw_output)

         if isinstance(data, list) and len(data) > 0:
                record = data[0]
         elif isinstance(data, dict):
                record = data
         else:
                print(f"Row {idx}: Unexpected JSON structure, skipping.")
                continue
         
         #gets chatgpt_taxon_id
         chatgpt_taxon_id = record.get("taxon_id", "")
         #gets chatgpt_name
         chatgpt_species_name = record.get("species_name", "")
         #gets chatgpt_reasoning
         chatgpt_reasoning = record.get("reasoning", "")

         #stores results
         df.at[idx, CHATGPT_COL] = chatgpt_species_name
         df.at[idx, TAXON_COL] = chatgpt_taxon_id
         df.at[idx, REASON_COL] = chatgpt_reasoning
         
         print(f"Parsed -> chatgpt_taxon_id: {chatgpt_taxon_id}," 
               f"species_name: {chatgpt_species_name},"
               f"chatgpt_reasoning: {chatgpt_reasoning}"
               )
         
         #comparing values
         input_species_name = str(row.get(NAMESPACENAME_COL, "")).strip()
         input_taxon_id = str(row.get(NAMESPACEID_COL, "")).strip()
         
         returned_name = str(chatgpt_species_name).strip()
         returned_taxon_id = str(chatgpt_taxon_id).strip()
         
         name_match = input_species_name.lower() == returned_name.lower()
         taxon_id_match = input_taxon_id == returned_taxon_id
         
         df.at[idx, NAME_MATCH_COL] = "Yes" if name_match else "No"
         df.at[idx, TAXON_MATCH_COL] = "Yes" if taxon_id_match else "No"
         
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

#prints "Done testing ollama." when everything is finished
print("Done testing ollama.")