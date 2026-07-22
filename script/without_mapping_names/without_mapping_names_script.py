import pandas as pd
from pathlib import Path
import json
import re
from openai import OpenAI
import time
from Bio import Entrez

#file locations
CSV_IN = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\without_mapping_names\manually_curated_data\updated_data_without_mapping_names.csv")
CSV_OUT = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\without_mapping_names\chatgpt_matched_without_mapping_names.csv")
PROMPT_FILE = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\script\LLM_prompts\with_JSON_array\LLM_prompt_without_reasoning_with_array.txt")
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
NCBI_SPECIES_COL = "ncbi_species_name"
#output column name
NCBI_SPECIES_MATCH_COL = "ncbi_species_name_match_?"
#output column name
NCBI_TAXON_COL = "ncbi_taxon_id"
#output column name
NCBI_TAXON_MATCH_COL = "ncbi_taxon_id_match_?"

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
     NCBI_SPECIES_COL,
     NCBI_SPECIES_MATCH_COL,
     NCBI_TAXON_COL,
     NCBI_TAXON_MATCH_COL,
]
for col in output_cols:
     if col in df.columns:
          df[col] = df[col].astype("object")

#what rows we're testing
#examples: [1:21] tested rows 3-22; [9:29] tested rows 11-30
df = df.iloc[0:100]

#makes sure the output columns exist
for col in [CHATGPT_COL, NAME_MATCH_COL, NCBI_SPECIES_COL, NCBI_SPECIES_MATCH_COL, NCBI_TAXON_COL, NCBI_TAXON_MATCH_COL]:
 #if the column doesn't exist yet, it creates it
    if col not in df.columns:
        df[col] = pd.NA

#biopython's entrez module to get ncbi taxonomy IDs
def lookup_taxonomy(species):
     try:
          #strict search first
          handle = Entrez.esearch(
               db="taxonomy",
               term=f'"{species}"[Scientific Name]'
          )

          record = Entrez.read(handle)
          handle.close()

          ids = record.get("IdList", [])
          
          #fallback search
          if not ids:
               handle = Entrez.esearch(
                    db="taxonomy",
                    term=species
                    )
               record = Entrez.read(handle)
               handle.close()
               
               ids = record.get("IdList", [])
          
          if not ids:
               return None, None
          
          taxid = ids[0]

          #get the taxonomy record
          handle = Entrez.efetch(
               db="taxonomy",
               id=taxid,
               retmode="xml"
          )

          records = Entrez.read(handle)
          handle.close()

          scientific_name = records[0]["ScientificName"]

          return taxid, scientific_name
     
     except Exception as e:
          print(f"NCBI lookup failed for {species}: {e}")
          return None, None

client = OpenAI(api_key=api_key)

#processes rows with JSON
BATCH_SIZE = 10
for start in range(0, len(df), BATCH_SIZE):
    batch = df.iloc[start:start+BATCH_SIZE]

    species_list = (
         batch[INPUT_NAME_COL]
         .fillna("")
         .astype(str)
         .str.strip()
         .tolist()
    )
    #inserts "name" into the prompt
    prompt = prompt_template.replace(
         "<<SPECIES>>",
         json.dumps(species_list, indent=2)
         )

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
            ]
          )
         
         #will answer with api usage
         if response.usage:
              print("\nToken usage:")
              print(f"Prompt tokens: {response.usage.prompt_tokens}")
              print(f"Completion tokens: {response.usage.completion_tokens}")
              print(f"Total tokens: {response.usage.total_tokens}")

         #gets API response and prints it
         raw_output = response.choices[0].message.content.strip()

         print(f"\nProcessing batch starting at row {start}")
         print("Raw model output:")
         print(raw_output)

         #converts JSON text into a clean python
         match = re.search(r"\[.*\]", raw_output, re.S)

         if not match:
              raise ValueError("No JSON array found in model output")
         
         json_text = match.group(0)
         records = json.loads(json_text)

         if not isinstance(records, list):
              raise ValueError("Expected a JSON array.")

         if len(records) != len(batch):
              raise ValueError(
                   f"Expected {len(batch)} results, got {len(records)}"
              )
         
         #makes sure the outputs are in the correct format
         required = {
              "species_name",
              "input_name"
         }

         for i, record in enumerate(records):
              missing = required - record.keys()
              
              if missing:
                   raise ValueError(
                        f"Record {i} is missing keys: {missing}"
                   )
         
         for i, ((idx, row), record) in enumerate(
              zip(batch.iterrows(), records)
              ):
               
               expected_input = str(species_list[i]).strip()
               returned_input = str(
                    record.get("input_name", "")
                    ).strip()

               if expected_input != returned_input:
                    raise ValueError(
                         f"Input order mismatch."
                         f"Expected '{expected_input}'."
                         f"but got '{returned_input}'."
                    )
               
               #gets chatgpt_taxon_id
               chatgpt_species_name = record.get("species_name", "")
               #gets ncbi_taxon_id from chatgpt_species_name
               if chatgpt_species_name == "No Match Found":
                    ncbi_species_name = None
                    ncbi_taxon_id = None
               else:
                    ncbi_taxon_id, ncbi_species_name = lookup_taxonomy(chatgpt_species_name)
               
               #stores results
               df.at[idx, CHATGPT_COL] = chatgpt_species_name
               df.at[idx, NCBI_SPECIES_COL] = ncbi_species_name
               df.at[idx, NCBI_TAXON_COL] = ncbi_taxon_id
               
               #prints results
               print(f"NCBI_species_name: {ncbi_species_name},"
                     f"NCBI_taxon_id: {ncbi_taxon_id},"
                     f"species_name: {chatgpt_species_name}"
                     )
               
               #comparing values
               input_species_name = (
                    str(row.get(NAMESPACENAME_COL, ""))
                    .strip()
               )
               input_taxon_id = (
                    str(row.get(NAMESPACEID_COL, ""))
                    .strip()
                    .replace(".0", "")
               )

               returned_name = str(chatgpt_species_name).strip()
               returned_ncbi_species_name = (
                    str(ncbi_species_name).strip()
                    if ncbi_species_name not in [None, "None", ""]
                    else ""
               )
               returned_ncbi_taxon_id = (
                    str(ncbi_taxon_id).strip().replace(".0", "")
                    if ncbi_taxon_id not in [None, "None", ""]
                    else ""
                    )
               
               name_match = input_species_name.lower() == returned_name.lower()
               ncbi_species_name_match = input_species_name.lower() == returned_ncbi_species_name.lower()
               ncbi_taxon_id_match = input_taxon_id == returned_ncbi_taxon_id
               
               df.at[idx, NAME_MATCH_COL] = "Yes" if name_match else "No"
               df.at[idx, NCBI_SPECIES_MATCH_COL] = "Yes" if ncbi_species_name_match else "No"
               df.at[idx, NCBI_TAXON_MATCH_COL] ="Yes" if ncbi_taxon_id_match else "No"
               
         #saves progress after every row
         df.to_csv(CSV_OUT, index=False)
         print(f"Batch starting at row {start} saved to {CSV_OUT}")

    #catches JSON errors
    except json.JSONDecodeError as e:
         print(f"Batch starting at row {start}: JSON decode error")
         print(f"Error: {e}")
         print("Raw output that failed to parse:")
         #shows what went wrong
         print(raw_output)
         #moves to next row
         continue
    #catches every other error
    except Exception as e:
         #shows what went wrong
         print(f"Batch starting at row {start}: {e}")
         #moves to next row
         continue
    #delays 1 second before moving to the next row
    finally:
         time.sleep(1)

#prints "Done." when everything is finished
print("Done.")