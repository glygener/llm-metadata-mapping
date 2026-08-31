# llm-metadata-mapping
Mapping database metadata (publication, species, tissue, disease) to established ontologies and dictionaries using an LLM.

## Repository Contents
### diseases
#### data
The data folder contains all of the data/datasets.
- "mapping_BS_Disease-mapped.xlsx" is the list of total diseases.

#### LLM prompts
The LLM prompts folder contains all of the LLM prompts.
- "diseases_LLM_prompt.txt" is the most recent prompt used to generate the CSV files.

#### script
The script folder contains the python scripts.
- Contains an empty .gitkeep file as a placeholder so that the folder could be committed.

### doc
The doc folder contains any useful spreadsheet or powerpoint presentation completed in conjunction to this project.

### species
#### data
The data folder contains all of the data/datasets. Within this folder are two separate folders and one CSV file.
- The "with mapping names" folder which contains the CSV files with species that have mapping name:
    - "chatgpt_matched_with_mapping_names.csv" is the file including the rows that were ran through the script/LLM (CSV_OUT).
    - "original_data_with_mapping_names.csv" is the file with the total list of species with mapping names (CSV_IN).
- The "without mapping names" folder which contains the CSV files with species that do not have mapping names:
    - "chatgpt_matched_without_mapping_names.csv" is the file including the rows that were ran through the script/LLM (CSV_OUT).
    - "original_data_without_mapping_names.csv" is the file with the total list of species without mapping names (CSV_IN).
- The "updated_list_of_species.csv" file which is the list of the total species.

#### LLM prompts
The LLM prompts folder contains all of the LLM prompts.
- "species_LLM_prompt_without_reasoning_with_array.txt" is the most recent prompt used to generate the CSV files.

#### script
The script folder contains the python scripts.
- "mapping_names_script.py" is the most recent script used to generate the CSV files with species including mapping names. 
- "without_mapping_names_script" is the the most recent script used to generate the CSV files with species that do not have mapping names.

### tissues
#### data
The data folder contains all of the data/datasets.
- "mapping_BS_Tissue-mapped.xlsx" is the list of total tissues.

#### LLM prompts
The LLM prompts folder contains all of the LLM prompts.
- "tissues_LLM_prompt.txt" is the most recent prompt used to generate the CSV files.

#### script
The script folder contains the python scripts.
- Contains an empty .gitkeep file as a placeholder so that the folder could be committed.

## CSV Header Descriptions
Here are descriptions detailing what data is located in each column within the CSV files.

### Species Data
#### Existing Data
The following columns include existing data from "updated_list_of_species.csv":

ID, count, name, namespacename, namespaceid, mappingname, rank, matchCount
- The "ID" column is the ID number provided to each species by the past volunteers who curated updated_list_of_species.csv. This ID is not relevant to this project and was solely for their own identification.
- The "count" column is information provided to each species by the past volunteers who curated updated_list_of_species.csv. This ID is not relevant to this project and was solely for their own identification.
- The "name" column is the species name being input into the LLM.
- The "namespacename" column is the correct scientific species name that should be returned by the LLM.
- The "namespaceid" column is the correct NCBI taxonomy ID that should be returned by the LLM.
- The "mappingname" column is the scientific name that was manually mapped to each species. This was only added to species that were problematic to curate, and was done by previous volunteers. The mapping name tends to be the same as the name within the "name" column.
    - This column only exists in datasets that include these problematic species. For example, all the datasets within the "without_mapping_names" folder located in the data folder do not have mapping names, and therefore do not have this column.
- The "rank" column is information provided to each species by the past volunteers who curated updated_list_of_species.csv. This ID is not relevant to this project and was solely for their own identification.
- The "matchCount" column is information provided to each species by the past volunteers who curated updated_list_of_species.csv. This ID is not relevant to this project and was solely for their own identification.
 
#### New data
The following columns are filled by the LLM:

chatgpt_name and chatgpt_reasoning
- The "chatgpt_name" column is the species name returned by the LLM.
- The "name_match_?" column is comparing the "namespacename" and "chatgpt_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
- The "chatgpt_reasoning" column contains the reasoning why the LLM produced that species name. Only some of the CSV files contain this column.


The following columns are filled by the Bio.Entrez package:

ncbi_species_name and ncbi_taxon_id
- The species name returned within "chatgpt_species_name" is put through this package and its species name output is returned in "ncbi_species_name". This column is the actual scientific species name from the NCBI database that corresponds to the "chatgpt_species_name"
- The "ncbi_species_name_match_?" column is comparing the "namespacename" and "ncbi_species_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
- The species name returned within "chatgpt_species_name" is put through this package and its NCBI taxonomy ID output is returned in "ncbi_taxon_id". This column is the actual taxonomy ID from the NCBI database that corresponds to the "chatgpt_species_name".
- The "ncbi_taxon_id_match_?" column is comparing the "namespaceid" and "ncbi_taxon_id" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.

Bio.Entrez Package: https://biopython.org/docs/1.76/api/Bio.Entrez.html
