# llm-metadata-mapping
Mapping database metadata (publication, species, tissue, disease) to established ontologies and dictionaries using an LLM.

The data folder contains the total list of species and four additional folders. These are: test_data which contains test CSV files, ollama_data which contains the ollama CSV files, with_mapping_names which contains the CSV files with species that have mapping names, and without_mapping_names which contains the CSV files with species that do not have mapping names.

The old_files folder present contains five separate folders containing several old files. The csv_rw folder contains scripts used to read and write CSV files, old_csv contains old CSV files, old_data ccontains old and reference data, old_prompt contains old LLM prompts, and old_script contains old scripts.

The script folder present contains five total folders: LLM_prompts, Test_script_and_prompt, general_python_scripts, with_mapping_names, and without_mapping_names. The LLM_prompts folder contains two folders: with_JSON_array that has the prompt with details on reporting the data as a JSON array, and without_JSON_array that has the prompt with details on reporting the data as a JSON output. The test_script_and_prompt folder contains the test LLM prompt and script. The general_python_scripts folder contains two folders: with_JSON_array that has the script that reports the data as a JSON array, and without_JSON_array that has the script that reports the data as a JSON output. The with_mapping_names folder has the most recent script used to generate the CSV files with species including mapping names.  The without_mapping_names folder has the most recent script used to generate the CSV files with species that do not have mapping names. 

# CSV Headers
The following headers/columns already had existing data: ID, count, name, namespacename, namespaceid, mappingname, rank, matchCount.

The "name" column is the species name being input into the LLM.

The "namespacename" column is the correct scientific species name that should be returned by the LLM.

The "namespaceid" column is th correct NCBI taxonomy ID that should be returned by the LLM.

The following headers/columns are filled by the LLM: chatgpt_name and chatgpt_taxon_id.

The "chatgpt_name" column is the species name returned by the LLM.

The "name_match_?" column is comparing the "namespacename" and "chatgpt_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.

The "chatgpt_taxon_id" column is the NCBI taxonomy ID returned by the LLM.

The "taxon_id_match_?" column is comparing the "namespaceid" and "chatgpt_taxon_id" columms to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.

The following header/column is filled in by hand: chatgpt_taxon_id_species_name.

The "chatgpt_taxon_id_species_name" is the actual name of the species that the "chatgpt_taxon_id" returned. If the LLM returned the correct taxon ID, then this column is left blank. This column is only filled out if the LLM returned the incorrect taxonomy ID.

The following headers/columns are filled by the Bio.Entrez package: ncbi_species_name and ncbi_taxon_id.

Bio.Entrez Package: https://biopython.org/docs/1.76/api/Bio.Entrez.html

The species name returned within "chatgpt_species_name" is put through this package and its species name output is returned in "ncbi_species_name". This column is the actual scientific species name from the NCBI database that corresponds to the "chatgpt_species_name"

The "ncbi_species_name_match_?" column is comparing the "namespacename" and "ncbi_species_name" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.

The species name returned within "chatgpt_species_name" is put through this package and its NCBI taxonomy ID output is returned in "ncbi_taxon_id". This column is the actual taxonomy ID from the NCBI database that corresponds to the "chatgpt_species_name".

The "ncbi_taxon_id_match_?" column is comparing the "namespaceid" and "ncbi_taxon_id" columns to see if their contents match. If they are a match, 'Yes' is returned. If not, 'No' is returned.
