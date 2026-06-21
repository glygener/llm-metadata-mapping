import pandas as pd
from pathlib import Path
import json
from openai import OpenAI
import time
from Bio import Entrez

#file locations
CSV_IN = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\with_mapping_names\data_with_mapping_names.csv")
CSV_OUT = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\with_mapping_names\chatgpt_matched_with_mapping_names.csv")
PROMPT_FILE = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\script\LLM_prompt.txt")
API_FILE = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\script\api_key.txt")

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
NCBI_TAXON_COL = "ncbi_taxon_id"
#output column name
NCBI_TAXON_MATCH_COL = "ncbi_taxon_id_match_?"
#output column name
REASON_COL = "chatgpt_reasoning"

#API model name
MODEL = "gpt-5.5"

#reading files
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")

#api key
api_key = API_FILE.read_text(encoding="utf-8").strip()

#email for biopython's entrez module
Entrez.email = "taylor.dimenna@gwmail.gwu.edu"

#loads csv into pandas
df = pd.read_csv(CSV_IN)

#converts columns to object
output_cols= [
     CHATGPT_COL,
     NAME_MATCH_COL,
     TAXON_COL,
     TAXON_MATCH_COL,
     NCBI_TAXON_COL,
     NCBI_TAXON_MATCH_COL,
     REASON_COL,
]
for col in output_cols:
     if col in df.columns:
          df[col] = df[col].astype("object")

#tests rows 20-39
df = df.iloc[20:40]

#makes sure the output columns exist
for col in [CHATGPT_COL, NAME_MATCH_COL, TAXON_COL, TAXON_MATCH_COL, NCBI_TAXON_COL, NCBI_TAXON_MATCH_COL, REASON_COL]:
 #if the column doesn't exist yet, it creates it
    if col not in df.columns:
        df[col] = pd.NA

#biopython's entrez module to get ncbi taxonomy IDs
def lookup_taxonomy_id(species):
     try:
          #strict search first
          handle = Entrez.esearch(
               db="taxonomy",
               term=f'"{species}"[Scientific Name]'
          )

          record = Entrez.read(handle)
          handle.close()

          ids = record.get("IdList", [])
          if ids:
               return ids [0]
          
          #fallback search
          handle = Entrez.esearch(
               db="taxonomy",
               term=species
          )
          record = Entrez.read(handle)
          handle.close()

          ids = record.get("IdList", [])
          return ids [0] if ids else None
     
     except Exception as e:
          print(f"NCBI lookup failed for {species}: {e}")
          return None

client = OpenAI(api_key=api_key)

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
              #forces the model into JSON mode
              response_format={"type":"json_object"},
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
         #gets ncbi_taxon_id from chatgpt_species_name
         ncbi_taxon_id = lookup_taxonomy_id(chatgpt_species_name)
         
         if chatgpt_species_name == "No Match Found":
              ncbi_taxon_id = None
         else:
              ncbi_taxon_id = lookup_taxonomy_id(chatgpt_species_name)

         #stores results
         df.at[idx, CHATGPT_COL] = chatgpt_species_name
         df.at[idx, TAXON_COL] = chatgpt_taxon_id
         df.at[idx, NCBI_TAXON_COL] = ncbi_taxon_id
         df.at[idx, REASON_COL] = chatgpt_reasoning
         
         print(f"chatgpt_taxon_id: {chatgpt_taxon_id}," 
               f"NCBI_taxon_id: {ncbi_taxon_id},"
               f"species_name: {chatgpt_species_name},"
               f"chatgpt_reasoning: {chatgpt_reasoning}"
               )
         
         #comparing values
         input_species_name = str(row.get(NAMESPACENAME_COL, "")).strip()
         input_taxon_id = str(row.get(NAMESPACEID_COL, "")).strip().replace(".0", "")
         
         returned_name = str(chatgpt_species_name).strip()
         returned_chatgpt_taxon_id = str(chatgpt_taxon_id).strip().replace(".0", "")
         returned_ncbi_taxon_id = (
              str(ncbi_taxon_id).strip().replace(".0", "")
              if ncbi_taxon_id not in [None, "None", ""]
              else ""
         )

         name_match = input_species_name.lower() == returned_name.lower()
         
         taxon_id_match = input_taxon_id == returned_chatgpt_taxon_id
         ncbi_taxon_id_match = input_taxon_id == returned_ncbi_taxon_id
         
         df.at[idx, NAME_MATCH_COL] = "Yes" if name_match else "No"
         df.at[idx, TAXON_MATCH_COL] = "Yes" if taxon_id_match else "No"
         df.at[idx, NCBI_TAXON_MATCH_COL] ="Yes" if ncbi_taxon_id_match else "No"

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
    #delays 1 second before moving to the next row
    finally:
         time.sleep(1)

#prints "Done." when everything is finished
print("Done.")