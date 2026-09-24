# Timeline of LLM Prompt Edits and Their Effects on Response Quality
This document tracks changes made to ```species_LLM_prompt_without_reasoning_with_array.txt``` and the results produced by each.

## Summer 2026
### Change 1 - Initial prompt
**Commit**: [`867f67d`](https://github.com/glygener/llm-metadata-mapping/commit/867f67d29026d735582d0f33bed2c27a4669f32f)

**Commit Name**: "Added to LLM prompt"

**Date**: 6/11/26

**Author**: Taylor DiMenna

**Changes**: No changes yet

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/867f67d29026d735582d0f33bed2c27a4669f32f/script/modified_prompt.txt)

**Results**: No results yet

### Change 2 - JSON Output
**Prompt Commit**: [`3d145fe`](https://github.com/glygener/llm-metadata-mapping/commit/3d145fe825b1dc413923206a71d12fe4e0078e45)

**Commit Name**: "Added JSON return specifics to prompt"

**Date**: 6/12/26

**Author**: Taylor DiMenna

**Changes**:
- Report output as a JSON

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/3d145fe825b1dc413923206a71d12fe4e0078e45/script/modified_prompt.txt)

**Results Commit**: [`f35e271`](https://github.com/glygener/llm-metadata-mapping/commit/f35e2713dda50951968c91110fb8a2331b21dba0)
- Tested rows 1-19
- With mapping names dataset
- GPT-5.5
- Name match 50%
- ID match 33%
- Both match 16%

**Commit Name**: "Tested rows 1-19"

**Date**: 6/18/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/f35e2713dda50951968c91110fb8a2331b21dba0/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 3 - Added Distinct Steps
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
- With mapping names dataset
- GPT-5.5
- Name match 47%
- ID match 37%
- Both match 21%

**Commit Name**: "Results with altered prompt"

**Date**: 6/18/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/c42a422bcc88f4fdc6a745b91ed87d48e3647d62/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 4 - Removed Repetitiveness
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
- With mapping names dataset
- GPT-5.5
- Name match 42%
- ID match 32%
- Both match 21%

**Commit Name**: "Data from shorter prompt"

**Date**: 6/20/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/b442876b620b15469f4f19ed4ed6e4e087608a9e/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 5 - Added Bio.Entrez
**Prompt Commit**: No prompt change

**Changes**:
- Added Bio.Entrez package to Python script

**Results Commit**: [`c98a45a`](https://github.com/glygener/llm-metadata-mapping/commit/c98a45ad76a8544de855c991bf54905747de8203)
- Tested rows 1-19
- With mapping names dataset
- GPT-5.5
- Not referencing the Bio.Entrez results yet
- Name match 47%
- ID match 32%
- Both match 32%

**Commit Name**: "Updated with Bio.Entrez taxon IDs"

**Date**: 6/20/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/c98a45ad76a8544de855c991bf54905747de8203/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 6 - Shortened Prompt & Fixed Bio.Entrez
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
- With mapping names dataset
- GPT-5.5
- Not referencing Bio.Entrez results yet
- Name match 42%
- ID match 42%
- Both match 32%

**Commit Name**: "Retested with shorter prompt"

**Date**: 6/21/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/f30ae49a45f4365592d2be9fbb1fd4e57b961463/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 7 - Shortened Prompt & Fixed Bio.Entrez in Script
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
- With mapping names dataset
- GPT-4o
- Not referencing Bio.Entrez results yet
- Name match 40%
- ID match 0%
- Both match 0%

**Commit Name**: "Results with gpt-4o"

**Date**: 6/24/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/be4a3387256443838c564521f11ca2afe35ed139/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 8 - Multiple Changes
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
- With mapping names dataset
- GPT-4o
- Not referencing Bio.Entrez results yet
- Name match 25%
- ID match 10%
- Both match 0%

**Date**: 6/29/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/6b873a3dd723c5dac116e213c5aead6f32862b47/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 9 - Strengthened "No Match Found" Directions
**Prompt Commit**: [`a1fc7a6`](https://github.com/glygener/llm-metadata-mapping/commit/a1fc7a6f9a325b10350d9e3bc85564c004583c55)

**Commit Name**: "Added Bio.Entrez to return species names"

**Date**: 6/29/26

**Author**: Taylor DiMenna

**Changes**:
- Added more directions asking not to return the input name because no other name was found. Emphasizing that "No Match Found" is a better result

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/a1fc7a6f9a325b10350d9e3bc85564c004583c55/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/a1fc7a6f9a325b10350d9e3bc85564c004583c55)
- Tested rows 1-20
- With mapping names dataset
- GPT-4o
- Not referencing Bio.Entrez results yet, just added the comparison rows (ncbi_species_name_match_? and ncbi_taxon_id_match_? rows) to CSV file 
- Name match 30%
- ID match 10%
- Both match 5%

**Date**: 6/29/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/a1fc7a6f9a325b10350d9e3bc85564c004583c55/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 10 - Force "No Match Found"
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
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results; the species name and taxon ids are fetched from NCBI, so they always refer to the same record. They are either both correct, both incorrect, or "No Match Found" was returned
- Correct 70%
- Incorrect 5%
- No Match Found 15%

**Commit Name**: "Results forcing "No Match Found""

**Date**: 7/1/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/8ffee746904c80cc5d5c5663660835e819e46edf/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 11 - Multiple Changes
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
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Correct 75%
- Incorrect 20%
- No Match Found 5%

**Date**: 7/5/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/6f191323cecce066423e6d1de028f6933ff04863/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 12 - Multiple Changes
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
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Correct 55%
- Incorrect 20%
- No Match Found 25%

**Date**: 7/5/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/4877a770b45d166e1f562f62813e75fd24569eab/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 13 - Add Back Taxonomy ID Requirements
**Prompt Commit**: [`93bc451`](https://github.com/glygener/llm-metadata-mapping/commit/93bc4511dee68b3e918af43eb8852d054c529c8a)

**Commit Name**: "Add back taxon ID prompt details"

**Date**: 7/6/26

**Author**: Taylor DiMenna

**Changes**:
- Added 'taxonomy ids must' section back

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/93bc4511dee68b3e918af43eb8852d054c529c8a/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/93bc4511dee68b3e918af43eb8852d054c529c8a)
- Tested rows 50-70
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Correct 70%
- Incorrect 15%
- No Match Found 15%

**Date**: 7/6/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/93bc4511dee68b3e918af43eb8852d054c529c8a/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 14 - Add Another Accepted Name Conversion
**Prompt Commit**: [`34f94f2`](https://github.com/glygener/llm-metadata-mapping/commit/34f94f28b7cf925c6c4d4611ce7ba62bf5a5a8e4)

**Commit Name**: "Added another positive example"

**Date**: 7/10/26

**Author**: Taylor DiMenna

**Changes**:
- Added another positive accepted name conversion example

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/34f94f28b7cf925c6c4d4611ce7ba62bf5a5a8e4/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [`5124293`](https://github.com/glygener/llm-metadata-mapping/commit/512429335c20c19d45daa7ca8299c4f630d62e9d)
- Tested rows 21-41
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Correct 80%
- Incorrect 10%
- No Match Found 10%

**Commit Name**: "Tested rows 21-41"

**Date**: 7/10/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/512429335c20c19d45daa7ca8299c4f630d62e9d/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

### Change 15 - Limit Conservative Results & Shorten Prompt
**Prompt Commit**: [`4f5668c`](https://github.com/glygener/llm-metadata-mapping/commit/4f5668cf5bb297b35b70e1c11b33095c072bc135)

**Commit Name**: "Updated to try and limit conservative results"

**Date**: 7/14/26

**Author**: Taylor DiMenna

**Changes**:
- Changed wording of 'role' to be more specific in when to return "No Match Found
- Changed wording of 'objective' to be more specific in when to return "No Match Found
- Make procedure more straightforward by shortening the wording of each step
- Remove "taxon_id" from example JSON result. At this point, we stopped asking the LLM to try and return a taxon id and only using what Bio.Entrez returned.
    - This was altered in the script
- Combine steps in the 'required rules' section
- Delete repetitive steps in the 'required rules' section
- Added 2 example results

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/4f5668cf5bb297b35b70e1c11b33095c072bc135/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**:[`16709e2`](https://github.com/glygener/llm-metadata-mapping/commit/16709e2d44fc98cfcf3bb2cba7b62cd65b40d417)
- GPT-4o
- Referencing Bio.Entrez results
- With mapping names dataset:
    - Tested rows 1-100
    - Correct 71%
    - Incorrect 12%
    - No Match Found 17%
- Without mapping names dataset:
    - Tested rows 1-100
    - Correct 83%
    - Incorrect 0%
    - No Match Found 17%

**Commit Name**: "Retested rows 1-100 with updated prompt"

**Date**: 7/14/26

**Author**: Taylor DiMenna

[`With Mapping Names Results file`](https://github.com/glygener/llm-metadata-mapping/blob/16709e2d44fc98cfcf3bb2cba7b62cd65b40d417/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

[`Without Mapping Names Results file`](https://github.com/glygener/llm-metadata-mapping/blob/16709e2d44fc98cfcf3bb2cba7b62cd65b40d417/data/without_mapping_names/chatgpt_matched_without_mapping_names.csv)

### Change 16 - Add Back Taxonomy ID Requirements
**Prompt Commit**: [`816f48d`](https://github.com/glygener/llm-metadata-mapping/commit/816f48d0ca9bafa5da1b730701fd884ccb1bcaf0)

**Commit Name**: "Test again with the prompt with more details"

**Date**: 7/14/26

**Author**: Taylor DiMenna

**Changes**:
- Minimally change wording of 'objective' and procedure
- Add details to 'required rules' on not prefering older/more commonly used names over a less common one
- Add details to 'required rules' on returning the name accepted by NCBI when there is more than one taxonomic classification for the same taxon
- Switich a example result to a different species

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/816f48d0ca9bafa5da1b730701fd884ccb1bcaf0/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/816f48d0ca9bafa5da1b730701fd884ccb1bcaf0)
- Tested rows 50-70
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Correct 74%
- Incorrect 12%
- No Match Found 14%

**Date**: 7/14/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/816f48d0ca9bafa5da1b730701fd884ccb1bcaf0/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)

## Fall 2026
These prompt changes were all in an effort to address the following issues:

#9- [NCBI taxonomy focus for the prompt](https://github.com/glygener/llm-metadata-mapping/issues/9)

#10- [Genus transfer issue](https://github.com/glygener/llm-metadata-mapping/issues/10)

#11- [Solely because a name is unfamiliar](https://github.com/glygener/llm-metadata-mapping/issues/11)

The same changes were made for each issue and they were the only parts altered with the prommpt each time. This series of tests was more to determine which combination of edits would be the most effective.

To address [#9](https://github.com/glygener/llm-metadata-mapping/issues/9), the following was added to the prompt:
- A returned species_name must exist in the NCBI Taxonomy database as either:
    1. the scientific name of a taxon, or
    2. a synonym attached to an NCBI taxon.
- If you cannot determine that the returned name exists in NCBI Taxonomy, return "No Match Found".

To address [#10](https://github.com/glygener/llm-metadata-mapping/issues/10), the following was added to the prompt:
- Do not create a new genus/species combination. Do not infer a taxonomic transfer just because similar transfers exist in the literature.

To address [#11](https://github.com/glygener/llm-metadata-mapping/issues/11), the following was added to the prompt:
- If the exact NCBI accepted name cannot be determined with high confidence, return "No Match Found".

Each result produced by the edited prompt was compared to a baseline (results before the changes) to determine effectiveness.

**GPT-4o Baseline**: [`bdabb36`](https://github.com/glygener/llm-metadata-mapping/commit/bdabb360f235ee28c61ccd23c65def15162832b8)
- Tested rows 100-150
- With mapping names dataset
- Referencing Bio.Entrez results
- Correct 76%
- Incorrect 20%
- No Match Found 4%

**GPT-5.5 Baseline**: [`6324077`](https://github.com/glygener/llm-metadata-mapping/commit/632407706ecf0bb614e21c0114f2d1c470434f0b)
- Tested rows 100-150
- With mapping names dataset
- Referencing Bio.Entrez results
- Correct 86%
- Incorrect 14%
- No Match Found 0%

### Change 17 - Testing [#9](https://github.com/glygener/llm-metadata-mapping/issues/9) & [#11](https://github.com/glygener/llm-metadata-mapping/issues/11) With GPT-4o
**Prompt Commit**: [`90d7f6f`](https://github.com/glygener/llm-metadata-mapping/commit/90d7f6f13ea7e9b553aad29d3143a0334ac0a760)

**Commit Name**: "Edited to address issues #9 and #11"

**Date**: 9/1/26

**Author**: Taylor DiMenna

**Changes**:
- Added edits for #9
- Added edits for #11

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/90d7f6f13ea7e9b553aad29d3143a0334ac0a760/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/90d7f6f13ea7e9b553aad29d3143a0334ac0a760)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 68% correct
    - 20% incorrect
    - 12% no match found

**Date**: 9/1/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/90d7f6f13ea7e9b553aad29d3143a0334ac0a760/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 18 - Testing [#9](https://github.com/glygener/llm-metadata-mapping/issues/9), [#10](https://github.com/glygener/llm-metadata-mapping/issues/10), & [#11](https://github.com/glygener/llm-metadata-mapping/issues/11) With GPT-4o
**Prompt Commit**: [`2d0a99c`](https://github.com/glygener/llm-metadata-mapping/commit/2d0a99c07b3cd56e4b9ca0e078084c8502f05038)

**Commit Name**: "Edited to address issues #9, #10, and #11"

**Date**: 9/1/26

**Author**: Taylor DiMenna

**Changes**:
- Added edits for #10

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/2d0a99c07b3cd56e4b9ca0e078084c8502f05038/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/2d0a99c07b3cd56e4b9ca0e078084c8502f05038)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 60% correct
    - 28% incorrect
    - 12% no match found

**Date**: 9/1/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/2d0a99c07b3cd56e4b9ca0e078084c8502f05038/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 19 - Testing [#10](https://github.com/glygener/llm-metadata-mapping/issues/10) & [#11](https://github.com/glygener/llm-metadata-mapping/issues/11) With GPT-4o
**Prompt Commit**: [`c3efaa3`](https://github.com/glygener/llm-metadata-mapping/commit/c3efaa3ce9788819cf747d6b0b0a228da94b7075)

**Commit Name**: "Edited to address issues #10 and #11"

**Date**: 9/1/26

**Author**: Taylor DiMenna

**Changes**:
- Deleted edits for #9

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/c3efaa3ce9788819cf747d6b0b0a228da94b7075/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/c3efaa3ce9788819cf747d6b0b0a228da94b7075)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 58% correct
    - 20% incorrect
    - 22% no match found

**Date**: 9/1/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/c3efaa3ce9788819cf747d6b0b0a228da94b7075/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 20 - Testing [#11](https://github.com/glygener/llm-metadata-mapping/issues/11) With GPT-4o
**Prompt Commit**: [`f126dea`](https://github.com/glygener/llm-metadata-mapping/commit/f126deabb7656020377ec2530ff466ae3dc27593)

**Commit Name**: "Edited to only address issue #11"

**Date**: 9/2/26

**Author**: Taylor DiMenna

**Changes**:
- Deleted edits for #10

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/f126deabb7656020377ec2530ff466ae3dc27593/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/f126deabb7656020377ec2530ff466ae3dc27593)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 64% correct
    - 40% incorrect
    - 6% no match found

**Date**: 9/2/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/f126deabb7656020377ec2530ff466ae3dc27593/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 21 - Testing [#10](https://github.com/glygener/llm-metadata-mapping/issues/10) With GPT-4o
**Prompt Commit**: [`11e1a37`](https://github.com/glygener/llm-metadata-mapping/commit/11e1a37ae342587f9f82a357be2b6a68adb50d22)

**Commit Name**: "Edited to only address issue #10"

**Date**: 9/2/26

**Author**: Taylor DiMenna

**Changes**:
- Deleted edits for 11
- Add edits for 10

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/11e1a37ae342587f9f82a357be2b6a68adb50d22/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/11e1a37ae342587f9f82a357be2b6a68adb50d22)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 64% correct
    - 22% incorrect
    - 14% no match found

**Date**: 9/2/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/11e1a37ae342587f9f82a357be2b6a68adb50d22/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)


### Change 21 - Testing [#9](https://github.com/glygener/llm-metadata-mapping/issues/9) & [#10](https://github.com/glygener/llm-metadata-mapping/issues/10) With GPT-4o
**Prompt Commit**: [`22d70eb`](https://github.com/glygener/llm-metadata-mapping/commit/22d70eb6d0f08bf6aba49da025bcf3e409e5c920)

**Commit Name**: "Edited to address issues #9 and #10"

**Date**: 9/2/26

**Author**: Taylor DiMenna

**Changes**:
- Add edits for 9

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/22d70eb6d0f08bf6aba49da025bcf3e409e5c920/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/22d70eb6d0f08bf6aba49da025bcf3e409e5c920)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 74% correct
    - 16% incorrect
    - 10% no match found

**Date**: 9/2/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/22d70eb6d0f08bf6aba49da025bcf3e409e5c920/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 22 - Testing [#9](https://github.com/glygener/llm-metadata-mapping/issues/9) With GPT-4o
**Prompt Commit**: [`2a2fdc6`](https://github.com/glygener/llm-metadata-mapping/commit/2a2fdc6b94d869cd25325d5452633cc3f787b2bb)

**Commit Name**: "Edited to only address issue #9"

**Date**: 9/2/26

**Author**: Taylor DiMenna

**Changes**:
- Deleted edits for #10

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/2a2fdc6b94d869cd25325d5452633cc3f787b2bb/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/2a2fdc6b94d869cd25325d5452633cc3f787b2bb)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 56% correct
    - 30% incorrect
    - 14% no match found

**Date**: 9/2/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/2a2fdc6b94d869cd25325d5452633cc3f787b2bb/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 23 - Testing [#9](https://github.com/glygener/llm-metadata-mapping/issues/9) & [#10](https://github.com/glygener/llm-metadata-mapping/issues/10) With GPT-5.5
**Prompt Commit**: [`f746312`](https://github.com/glygener/llm-metadata-mapping/commit/f746312d8c15b67f5ee092e135d5879265819f45)

**Commit Name**: "Edited to address issues #9 and #10 and test with GPT-5.5"

**Date**: 9/14/26

**Author**: Taylor DiMenna

**Changes**:
- Added edits for #10

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/f746312d8c15b67f5ee092e135d5879265819f45/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/f746312d8c15b67f5ee092e135d5879265819f45)
- Tested rows 100-150
- With mapping names dataset
- GPT-5.5
- Referencing Bio.Entrez results
- Baseline results:
    - 86% correct
    - 14% incorrect
    - 0% no match found
- New results:
    - 88% correct
    - 10% incorrect
    - 2% no match found

**Date**: 9/14/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/f746312d8c15b67f5ee092e135d5879265819f45/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 24 - Testing [#10](https://github.com/glygener/llm-metadata-mapping/issues/10) & [#11](https://github.com/glygener/llm-metadata-mapping/issues/11) With GPT-5.5
**Prompt Commit**: [`892c8b9`](https://github.com/glygener/llm-metadata-mapping/commit/892c8b9434ee91417388f2eb21cca1dfb49d0262)

**Commit Name**: "Edited to address issues #10 and #11 and test with GPT-5.5"

**Date**: 9/16/26

**Author**: Taylor DiMenna

**Changes**:
- Deleted edits for #9
- Added edits for #11

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/892c8b9434ee91417388f2eb21cca1dfb49d0262/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/892c8b9434ee91417388f2eb21cca1dfb49d0262)
- Tested rows 100-150
- With mapping names dataset
- GPT-5.5
- Referencing Bio.Entrez results
- Baseline results:
    - 86% correct
    - 14% incorrect
    - 0% no match found
- New results:
    - 88% correct
    - 12% incorrect
    - 0% no match found

**Date**: 9/16/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/892c8b9434ee91417388f2eb21cca1dfb49d0262/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

### Change 25 - Testing Suggestions from ChatGPT With GPT-4o & GPT-5.5
**Prompt Commit**: [`4713563`](https://github.com/glygener/llm-metadata-mapping/commit/471356314b87bff835f87f80861f18d6ad30b3ff)

**Commit Name**: "Test suggestions from ChatGPT with GPT-4o"

**Date**: 9/17/26

**Author**: Taylor DiMenna

**Changes**:
- Major rewording of 'role' and 'objective' sections
- Deleted procedure section and incorporated the directions throughout the other added sections
- Restructured the examples of accepted name conversion section
- Added 'core rules', 'candidate generation', 'candidate validation', 'orthographic and spelling varients', accepted NCBI names', 'taxonomic rank', 'NCBI requirement', and 'abstention rule' sections. These are essentially the procedure section from the previous prompts split into multiple different sections with added details
- Reworded 'required rules' section
- Added details in the 'final instructin' section about enforcing "No Match Found"

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/471356314b87bff835f87f80861f18d6ad30b3ff/species/LLM%20prompts/species_LLM_prompt_without_reasoning_with_array.txt)

**GPT-4o Results Commit**: [Same as prompt](https://github.com/glygener/llm-metadata-mapping/commit/471356314b87bff835f87f80861f18d6ad30b3ff)
- Tested rows 100-150
- With mapping names dataset
- GPT-4o
- Referencing Bio.Entrez results
- Baseline results:
    - 76% correct
    - 20% incorrect
    - 4% no match found
- New results:
    - 58% correct
    - 18% incorrect
    - 24% no match found

**Date**: 9/17/26

**Author**: Taylor DiMenna

[`GPT-4o Results file`](https://github.com/glygener/llm-metadata-mapping/blob/471356314b87bff835f87f80861f18d6ad30b3ff/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

**GPT-5.5 Results Commit**: [38071ed](https://github.com/glygener/llm-metadata-mapping/commit/38071edaf1490d8aa59e59410d7d0916ad7b7e2a)
- Tested rows 100-150
- With mapping names dataset
- GPT-5.5
- Referencing Bio.Entrez results
- Baseline results:
    - 86% correct
    - 14% incorrect
    - 0% no match found
- New results:
    - 84% correct
    - 14% incorrect
    - 2% no match found

**Commit Name**: "Test suggestions from ChatGPT with GPT-5.5"

**Date**: 9/17/26

**Author**: Taylor DiMenna

[`GPT-5.5 Results file`](https://github.com/glygener/llm-metadata-mapping/blob/38071edaf1490d8aa59e59410d7d0916ad7b7e2a/species/data/with%20mapping%20names/chatgpt_matched_with_mapping_names.csv)

## Current Best Prompt
**Prompt Commit**: [`816f48d`](https://github.com/glygener/llm-metadata-mapping/commit/816f48d0ca9bafa5da1b730701fd884ccb1bcaf0)

**Commit Name**: "Test again with the prompt with more details"

**Date**: 7/14/26

**Author**: Taylor DiMenna

[`Prompt file`](https://github.com/glygener/llm-metadata-mapping/blob/816f48d0ca9bafa5da1b730701fd884ccb1bcaf0/script/LLM_prompts/with_JSON_array/LLM_prompt_without_reasoning_with_array.txt)

**Without Mapping Names Results Commit**: [`d9956e6`](https://github.com/glygener/llm-metadata-mapping/commit/d9956e6d015641c503f5be25931703e7da4e091e)
- Tested rows 1-100
- Without mapping names dataset
- GPT-5.5
- Referencing Bio.Entrez results
- 99% correct
- 0% incorrect
- 1% no match found

**Date**: 7/22/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/d9956e6d015641c503f5be25931703e7da4e091e/data/without_mapping_names/chatgpt_matched_without_mapping_names.csv)

**With Mapping Names Results Commit**: [`a98330f`](https://github.com/glygener/llm-metadata-mapping/commit/a98220f34352577577ad09e84bf6d48dff31cc4a)
- Tested rows 1-200
- With mapping names dataset
- GPT-5.5
- Referencing Bio.Entrez results
- 92% correct
- 7% incorrect
- 1% no match found

**Date**: 7/22/26

**Author**: Taylor DiMenna

[`Results file`](https://github.com/glygener/llm-metadata-mapping/blob/a98220f34352577577ad09e84bf6d48dff31cc4a/data/with_mapping_names/chatgpt_matched_with_mapping_names.csv)
