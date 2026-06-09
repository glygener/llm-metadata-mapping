import pandas as pd
from pathlib import Path

#tells the script where the prompt file is located
BASE_DIR = Path(__file__).parent

#file locations
CSV_IN = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\updated_mapping_BS_GS-mapped.csv")
CSV_OUT = Path(r"C:\Users\taylo\OneDrive\Desktop\Github\Test\llm-metadata-mapping\data\Test_data\test.csv")
PROMPT_FILE = BASE_DIR/"test_prompt.txt"

#input column name
SPECIES_COL = "name"
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

#reading files
prompt_template = PROMPT_FILE.read_text(encoding="utf-8")

#loads csv into pandas
df = pd.read_csv(CSV_IN)

if CHATGPT_COL in df.columns:
    df[CHATGPT_COL] = df[CHATGPT_COL].astype("object")

#skips row 0 and only keeps 1-149
df = df.iloc[1:150]

#makes sure the output columns exist
for col in [CHATGPT_COL, NAME_MATCH_COL, TAXON_COL, TAXON_MATCH_COL, REASON_COL]:
 #if the column doesn't exist yet, it creates it
    if col not in df.columns:
        df[col] = pd.Series(index=df.index, dtype="object")

#curate results (fakes the API portion)
def curate(prompt, speciesname):
     results = {
          "taxon_id": "455045",
          "species_name": "Abelmoschus esculentus",
          "reasoning": "Fixed output for testing",
          "input": speciesname
     }
     return results

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
         record = curate(prompt, species)

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
         input_name = str(row.get(SPECIES_COL, "")).strip()
         input_taxon_id = str(row.get(NAMESPACEID_COL, "")).strip()
         
         returned_name = str(chatgpt_species_name).strip()
         returned_taxon_id = str(chatgpt_taxon_id).strip()
         
         name_match = input_name.lower() == returned_name.lower()
         taxon_id_match = input_taxon_id == returned_taxon_id
         
         df.at[idx, NAME_MATCH_COL] = "Yes" if name_match else "No"
         df.at[idx, TAXON_MATCH_COL] = "Yes" if taxon_id_match else "No"
         
         #saves progress after every row
         df.to_csv(CSV_OUT, index=False)
         print(f"Row {idx}: saved to {CSV_OUT}")

    #catches every other error
    except Exception as e:
         #shows what went wrong
         print(f"Row {idx}: Error calling API or processing result: {e}")
         #moves to next row
         continue

#prints "Done." when everything is finished
print("Done.")
