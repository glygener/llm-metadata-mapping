# llm-metadata-mapping
Mapping database metadata (publication, species, tissue, disease) to established ontologies and dictionaries using an LLM.

## Repository Contents
### data
The data folder contains all of the data/datasets. Within this folder are two separate folders and one CSV file.
- The "Species with mapping names" folder which contains the CSV files with species that have mapping name:
    - "updated_data_with_mapping_names.csv" is the file with the total list of species with mapping names (CSV_IN).
    - "chatgpt_matched_with_mapping_names.csv" is the file including the rows that were ran through the script/LLM (CSV_OUT).
- The "Species without mapping names" folder which contains the CSV files with species that do not have mapping names:
    - "updated_data_without_mapping_names.csv" is the file with the total list of species without mapping names (CSV_IN).
    - "chatgpt_matched_without_mapping_names.csv" is the file including the rows that were ran through the script/LLM (CSV_OUT).
- The "updated_list_of_species.csv" file which is the list of the total species.

### doc
The doc folder contains any useful spreadsheet or powerpoint presentation completed in conjunction to this project.

### LLM prompts
The LLM prompts folder contains all of the LLM prompts.
- "LLM_prompt_without_reasoning_with_array.txt" is the most recent prompt used to generate the CSV files.

### script
The script folder contains the python scripts.
- "mapping_names_script.py" is the most recent script used to generate the CSV files with species including mapping names. 
- "without_mapping_names_script" is the the most recent script used to generate the CSV files with species that do not have mapping names.

## CSV Header Descriptions
Here are descriptions detailing what data is located in each column within the CSV files.

### Existing Data
The following columns had existing data from the updated_list_of_species.csv: 
ID, count, name, namespacename, namespaceid, mappingname, rank, matchCount.
- The "name" column is the species name being input into the LLM.
- The "namespacename" column is the correct scientific species name that should be returned by the LLM.
- The "namespaceid" column is the correct NCBI taxonomy ID that should be returned by the LLM.
- The "mappingname" column is the scientific name that was manually mapped to each species. This was only added to species that were problematic to curate, and was done by previous volunteers. The mapping name tends to be the same as the name within the "name" column.
    - This column only exists in datasets that include these problematic species. For example, all the datasets within the "without_mapping_names" folder located in the data folder do not have mapping names, and therefore do not have this column.

### New data
The following headers/columns are filled by the LLM:
chatgpt_name, chatgpt_reasoning.
- The "chatgpt_name" column is the species name returned by the LLM.
- The "name_match_?" column is comparing the "namespacename" and "chatgpt_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
- The "chatgpt_reasoning" column contains the reasoning why the LLM produced that species name. Only some of the CSV files contain this column.

The following headers/columns are filled by the Bio.Entrez package:
ncbi_species_name and ncbi_taxon_id.
- The species name returned within "chatgpt_species_name" is put through this package and its species name output is returned in "ncbi_species_name". This column is the actual scientific species name from the NCBI database that corresponds to the "chatgpt_species_name"
- The "ncbi_species_name_match_?" column is comparing the "namespacename" and "ncbi_species_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
- The species name returned within "chatgpt_species_name" is put through this package and its NCBI taxonomy ID output is returned in "ncbi_taxon_id". This column is the actual taxonomy ID from the NCBI database that corresponds to the "chatgpt_species_name".
- The "ncbi_taxon_id_match_?" column is comparing the "namespaceid" and "ncbi_taxon_id" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.

Bio.Entrez Package: https://biopython.org/docs/1.76/api/Bio.Entrez.html
