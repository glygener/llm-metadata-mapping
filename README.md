# llm-metadata-mapping
Mapping database metadata (publication, species, tissue, disease) to established ontologies and dictionaries using an LLM.

The data folder contains all of the data/datasets. It is broken up into five separate folders. These are: 
- "ollama_data" which contains the ollama CSV files
- "test_data" which contains the test CSV files
- "total_list_of_species" which contains the CSV file "updated_list_of_species.csv" that includes all of the species
- "with_mapping_names" which contains the CSV files with species that have mapping name:
    - "updated_data_with_mapping_names.csv" is the file with the total list of species with mapping names (CSV_IN)
    - "chatgpt_matched_with_mapping_names.csv" is the file including the rows that were ran through the script/LLM (CSV_OUT)
- "without_mapping_names" which contains the CSV files with species that do not have mapping names:
    - "updated_data_without_mapping_names.csv" is the file with the total list of species without mapping names (CSV_IN)
    - "chatgpt_matched_without_mapping_names.csv" is the file including the rows that were ran through the script/LLM (CSV_OUT)

The old_files folder contains five separate folders containing several old files. These are:
- "csv_rw" which contains scripts used to combine, read, and write CSV files
- "old_csv" which contains old CSV files
- "old_data" which contains old species lists and reference data
- "old_prompt" which contains old LLM prompts
- "old_script" which contains old scripts

The script folder contains all of the python scripts and LLM prompts. It is broken up into five separate folders. These are:
- "general_python_scripts" which contains all of the general python scripts. This means that the CSV_IN & CSV_OUT paths are blank, the API model name is left blank, and the rows being tested is left blank. This was done to keep each relevent script avaliable and able to be easily edited for each task. When edited, it was done to the scripts in the "with_mapping_names" and "without_mapping_names" folders. These general python scripts are used more references rather than to run tests. They are in separate folders depending on what specific code they contain. Other than these small changes, they are the same:
    - The "with_JSON_array" folder has the scripts that contain code for reporting the results as a JSON array. The two scripts located here then differ in whether they get the LLM to return the "chatgpt_reasoning" column.
    - The "without_JSON_array" folder has the scripts that do not contain code reporting the results as a regular JSON output. The two scripts located here then differ in whether they get the LLM to return the "chatgpt_reasoning" column.
- "LLM_prompts" which contains all of the LLM prompts. They are in two separate folders depending on their contents:
    - "with_JSON_array" has the prompts with details on reporting the data as a JSON array. The two prompts located here then differ in whether they return the "chatgpt_reasoning" column.
    - "without_JSON_array" has the prompts with details on reporting the data as a regular JSON output. The two prompts located here then differ in whether they return the "chatgpt_reasoning" column.
- "test_script_and_prompt" which contains the test LLM prompt and script
- "with_mapping_names" which contains the most recent script used to generate the CSV files with species including mapping names. This was the script that was edited and used to run tests, not the general python scripts. 
- "without_mapping_names" which contains the most recent script used to generate the CSV files with species that do not have mapping names. This was the script that was edited and used to run tests, not the general python scripts.

# CSV Headers
The following headers/columns already had existing data: ID, count, name, namespacename, namespaceid, mappingname, rank, matchCount.
- The "name" column is the species name being input into the LLM.
- The "namespacename" column is the correct scientific species name that should be returned by the LLM.
- The "namespaceid" column is the correct NCBI taxonomy ID that should be returned by the LLM.
- The "mappingname" column is the scientific name that was manually mapped to each species. This was only added to species that were problematic to curate, and was done by previous volunteers. The mapping name tends to be the same as the name within the "name" column.
    - This column only exists in datasets that include these problematic species. For example, all the datasets within the "without_mapping_names" folder located in the data folder do not have mapping names, and therefore do not have this column.

The following headers/columns are filled by the LLM: chatgpt_name, chatgpt_reasoning.
- The "chatgpt_name" column is the species name returned by the LLM.
- The "name_match_?" column is comparing the "namespacename" and "chatgpt_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
- The "chatgpt_reasoning" column contains the reasoning why the LLM produced that species name. Only some of the CSV files contain this column.

Bio.Entrez Package: https://biopython.org/docs/1.76/api/Bio.Entrez.html

The following headers/columns are filled by the Bio.Entrez package: ncbi_species_name and ncbi_taxon_id.
- The species name returned within "chatgpt_species_name" is put through this package and its species name output is returned in "ncbi_species_name". This column is the actual scientific species name from the NCBI database that corresponds to the "chatgpt_species_name"
- The "ncbi_species_name_match_?" column is comparing the "namespacename" and "ncbi_species_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
- The species name returned within "chatgpt_species_name" is put through this package and its NCBI taxonomy ID output is returned in "ncbi_taxon_id". This column is the actual taxonomy ID from the NCBI database that corresponds to the "chatgpt_species_name".
- The "ncbi_taxon_id_match_?" column is comparing the "namespaceid" and "ncbi_taxon_id" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
