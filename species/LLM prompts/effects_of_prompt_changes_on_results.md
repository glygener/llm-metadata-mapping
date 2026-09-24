# Timeline of LLM Prompt Edits and Their Effects on Response Quality
This document tracks changes made to ```species_LLM_prompt_without_reasoning_with_array.txt``` and the results produced by each.

## Change 1 - Initial prompt
**Commit**: [`867f67d`](https://github.com/glygener/llm-metadata-mapping/commit/867f67d29026d735582d0f33bed2c27a4669f32f)

**Commit Name**: "Added to LLM prompt"

**Date**: 6/11/26

**Author**: Taylor DiMenna

**Changes**: No changes yet

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/867f67d29026d735582d0f33bed2c27a4669f32f/script/modified_prompt.txt)

**Results**: No results yet

## Change 2 - JSON Output
**Prompt Commit**: [`3d145fe`](https://github.com/glygener/llm-metadata-mapping/commit/3d145fe825b1dc413923206a71d12fe4e0078e45)

**Commit Name**: "Added JSON return specifics to prompt"

**Date**: 6/12/26

**Author**: Taylor DiMenna

**Changes**:
- Report output as a JSON

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/3d145fe825b1dc413923206a71d12fe4e0078e45/script/modified_prompt.txt)

**Results Commit**: [`f35e271`](https://github.com/glygener/llm-metadata-mapping/commit/f35e2713dda50951968c91110fb8a2331b21dba0)
- Tested rows 1-19
- GPT-5.5
- Name match 50%
- ID match 33%
- Both match 16%

**Commit Name**: "Tested rows 1-19"

**Date**: 6/18/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/f35e2713dda50951968c91110fb8a2331b21dba0/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 3 - Added Distinct Steps
**Prompt Commit**: [`e17d50e`](https://github.com/glygener/llm-metadata-mapping/commit/e17d50e2be9e0d1887e580d0e223dec483d9b1ab)

**Commit Name**: "Altered prompt"

**Date**: 6/18/26

**Author**: Taylor DiMenna

**Changes**:
- Added more defined procedural steps to follow
- Deleted 'naming priority' section
- Added section names to each heading which were directly referenced in procedure.
    - Ex. "Read and note the information in "Section G- DISCLAIMER"..."
- Deleted an example result

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/e17d50e2be9e0d1887e580d0e223dec483d9b1ab/script/LLM_prompt.txt)

**Results Commit**: [`c42a422`](https://github.com/glygener/llm-metadata-mapping/commit/c42a422bcc88f4fdc6a745b91ed87d48e3647d62)
- Tested rows 1-19
- GPT-5.5
- Name match 47%
- ID match 37%
- Both match 21%

**Commit Name**: "Results with altered prompt"

**Date**: 6/18/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/c42a422bcc88f4fdc6a745b91ed87d48e3647d62/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 4 - Removed Repetitiveness
**Prompt Commit**: [`021dc1a`](https://github.com/glygener/llm-metadata-mapping/commit/021dc1a70603e5453e5e06a1a283c8a47cbb7dd2)

**Commit Name**: "Reduced redundancy in prompt"

**Date**: 6/20/26

**Author**: Taylor DiMenna

**Changes**:
- Combined definitions for a "match", "no match", and "alternate" names within one Definitions section
- Simplified procedure
- Added a required rules section
- Added a consistency check section
- Added "No Match Found" concept
- Deleted section names to each heading which were directly referenced in procedure.
- Deleted an example result

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/021dc1a70603e5453e5e06a1a283c8a47cbb7dd2/script/LLM_prompt.txt)

**Results Commit**: [`b442876`](https://github.com/glygener/llm-metadata-mapping/commit/b442876b620b15469f4f19ed4ed6e4e087608a9e)
- Tested rows 1-19
- GPT-5.5
- Name match 42%
- ID match 32%
- Both match 21%

**Commit Name**: "Data from shorter prompt"

**Date**: 6/20/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/b442876b620b15469f4f19ed4ed6e4e087608a9e/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 5 - Added Bio.Entrez
**Prompt Commit**: No prompt change

**Changes**:
- Added Bio.Entrez package to Python script

**Results Commit**: [`c98a45a`](https://github.com/glygener/llm-metadata-mapping/commit/c98a45ad76a8544de855c991bf54905747de8203)
- Tested rows 1-19
- GPT-5.5
- Not referencing the Bio.Entrez results yet
- Name match 47%
- ID match 32%
- Both match 32%

**Commit Name**: "Updated with Bio.Entrez taxon IDs"

**Date**: 6/20/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/c98a45ad76a8544de855c991bf54905747de8203/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 6 - Shortened Prompt & Fixed Bio.Entrez
**Prompt Commit**: [`45a12dc`](https://github.com/glygener/llm-metadata-mapping/commit/45a12dc7ff508b98ae3907268f121caedbbcc96c)

**Commit Name**: "Shortened prompt"

**Date**: 6/21/26

**Author**: Taylor DiMenna

**Changes**:
- Deleted definitions for a "match", "no match", and "alternate" names
- Requested every input name to be searched for alternate names (instead of only the previously defined 'alternate' names)
- Added an order of names to search
    - Ex. First for the exact name, then a synonym, then an equivalent name, etc.
- Added sentence asking for each species name and taxon id to be treated as a single record
- Added "No Match Found"
- Removed 'reasoning'
- Switched from JSON to JSON array

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/45a12dc7ff508b98ae3907268f121caedbbcc96c/script/LLM_prompt.txt)

**Results Commit**: [`f30ae49`](https://github.com/glygener/llm-metadata-mapping/commit/f30ae49a45f4365592d2be9fbb1fd4e57b961463)
- Tested rows 1-19
- GPT-5.5
- Not referencing Bio.Entrez results yet
- Name match 42%
- ID match 42%
- Both match 32%

**Commit Name**: "Retested with shorter prompt"

**Date**: 6/21/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/f30ae49a45f4365592d2be9fbb1fd4e57b961463/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 7 - Shortened Prompt & Fixed Bio.Entrez in Script
**Prompt Commit**: [`96a4f2a`](https://github.com/glygener/llm-metadata-mapping/commit/96a4f2a73470bef351a407c399a1785e998aa534)

**Commit Name**: "Adjusted to work better with gpt-4o"

**Date**: 6/24/26

**Author**: Taylor DiMenna

**Changes**:
- Specified "No Match Found" concept
    - Asking for "No Match Found" to be returned for species_name and for taxon_id to be left blank
- Add back 'alternate names' defintion
- Add order to search for alternate names
    - Ex. Synonyms, equivalent names, basionyms, etc.
- Combine some 'required rules'
- Delete 'consitency check'
- Added explainations on returning IDs

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/96a4f2a73470bef351a407c399a1785e998aa534/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [`be4a338`](https://github.com/glygener/llm-metadata-mapping/commit/be4a3387256443838c564521f11ca2afe35ed139)
- Tested rows 1-20
- GPT-4o
- Not referencing Bio.Entrez results yet
- Name match 40%
- ID match 0%
- Both match 0%

**Commit Name**: "Results with gpt-4o"

**Date**: 6/24/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/be4a3387256443838c564521f11ca2afe35ed139/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 8 - Multiple Changes
**Prompt Commit**: [`6b873a3`](https://github.com/glygener/llm-metadata-mapping/commit/6b873a3dd723c5dac116e213c5aead6f32862b47)

**Commit Name**: "Test new prompt with gpt-4o"

**Date**: 6/29/26

**Author**: Taylor DiMenna

**Changes**:
- Added details to 'role'
- Added details to 'objective'
- Added more details about reporting the JSON array
    - Ex. Not to use the previous processed name when determining another name
- Delete order for searching for alternate names
- Add failure example
- Reinforce "No Match Found" result
- Simplify directions about each species_name and taxon_id to be treated as the same record
- Add more directions on JSON array

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/6b873a3dd723c5dac116e213c5aead6f32862b47/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/6b873a3dd723c5dac116e213c5aead6f32862b47)
- Tested rows 1-20
- GPT-4o
- Not referencing Bio.Entrez results yet
- Name match 25%
- ID match 10%
- Both match 0%

**Date**: 6/29/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/6b873a3dd723c5dac116e213c5aead6f32862b47/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 9 - Strengthened "No Match Found" Directions
**Prompt Commit**: [`a1fc7a6`](https://github.com/glygener/llm-metadata-mapping/commit/a1fc7a6f9a325b10350d9e3bc85564c004583c55)

**Commit Name**: "Added Bio.Entrez to return species names"

**Date**: 6/29/26

**Author**: Taylor DiMenna

**Changes**:
- Added more directions asking not to return the input name because no other name was found. Emphasizing that "No Match Found" is a better result

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/a1fc7a6f9a325b10350d9e3bc85564c004583c55/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/a1fc7a6f9a325b10350d9e3bc85564c004583c55)
- Tested rows 1-20
- GPT-4o
- Not referencing Bio.Entrez results yet, just added the comparison rows (ncbi_species_name_match_? and ncbi_taxon_id_match_? rows) to CSV file 
- Name match 30%
- ID match 10%
- Both match 5%

**Date**: 6/29/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/a1fc7a6f9a325b10350d9e3bc85564c004583c55/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 10 - Force "No Match Found"
**Prompt Commit**: [`0ba9ce2`](https://github.com/glygener/llm-metadata-mapping/commit/0ba9ce2da0d265788da0a45186614600712b169c)

**Commit Name**: "Updated to force "No Match Found" instead of guessing"

**Date**: 7/1/26

**Author**: Taylor DiMenna

**Changes**:
- Added details to 'role' about returning "No Match Found"
- Added directions to check if the name was reclassified to a different taxonomic name
- Add details to the procedure about returning "No Match Found"
- Added an example of an accepted name conversion
- Added another name conversion failure example
- Added a "No Match Found" result example

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/0ba9ce2da0d265788da0a45186614600712b169c/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [`8ffee74`](https://github.com/glygener/llm-metadata-mapping/commit/8ffee746904c80cc5d5c5663660835e819e46edf)
- Tested rows 9-29
- GPT-4o
- Referencing Bio.Entrez results; the species name and taxon ids are fetched from NCBI, so they always refer to the same record. They are either both correct, both incorrect, or "No Match Found" was returned
- Correct 70%
- Incorrect 5%
- No Match Found 15%

**Commit Name**: "Results forcing "No Match Found""

**Date**: 7/1/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/8ffee746904c80cc5d5c5663660835e819e46edf/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 11 - Multiple Changes
**Prompt Commit**: [`6f19132`](https://github.com/glygener/llm-metadata-mapping/commit/6f191323cecce066423e6d1de028f6933ff04863)

**Commit Name**: "Tested rows again with updated prompt"

**Date**: 7/5/26

**Author**: Taylor DiMenna

**Changes**:
- Simplified details in 'role'
- Added orthographic varient and recombination of species names as examples of alternative names
- Added another accepted name conversion example
- Added directions not to stop at a broader accepted taxon if a more specific name exists
- Added directions to check if the input name was transferred to another genus or reclassified before returning "No Match Found"

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/6f191323cecce066423e6d1de028f6933ff04863/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/6f191323cecce066423e6d1de028f6933ff04863)
- Tested rows 30-50
- GPT-4o
- Referencing Bio.Entrez results; the species name and taxon ids are fetched from NCBI, so they always refer to the same record. They are either both correct, both incorrect, or "No Match Found" was returned
- Correct 75%
- Incorrect 20%
- No Match Found 5%

**Date**: 7/5/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/6f191323cecce066423e6d1de028f6933ff04863/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Change 12 - Multiple Changes
**Prompt Commit**: [`4877a77`](https://github.com/glygener/llm-metadata-mapping/commit/4877a770b45d166e1f562f62813e75fd24569eab)

**Commit Name**: "Tested rows 50-70 with more examples within the prompt"

**Date**: 7/5/26

**Author**: Taylor DiMenna

**Changes**:
- Added an example with correct and incorrect names
- Added rules for not substituting names within the same genus
- Added directions to return the entire species name, even it is very long
- Deleted 'taxonomy ids must' section

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/4877a770b45d166e1f562f62813e75fd24569eab/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/4877a770b45d166e1f562f62813e75fd24569eab)
- Tested rows 50-70
- GPT-4o
- Referencing Bio.Entrez results; the species name and taxon ids are fetched from NCBI, so they always refer to the same record. They are either both correct, both incorrect, or "No Match Found" was returned
- Correct 55%
- Incorrect 20%
- No Match Found 25%

**Date**: 7/5/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/4877a770b45d166e1f562f62813e75fd24569eab/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)
