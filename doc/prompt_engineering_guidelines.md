# Prompt Engineering Guidelines

## Introduction
The purpose of these guidelines are to help those interested/starting in prompt engineering get a faster understanding of what to avoid/include within prompts. These guidelines were created based on the results from the Summer 2026 GlyGen Biocuration Project, rather than simply being broad prompt engineering advice. Since these were lessons observed when completing a specific project, they may not be universally applicable to all projects.

## Key Lessons Learned
### Give the LLM a role.
- By providing a role, the LLM is able to understand the context of the prompt so it is able to do its job more tailored toward the results you want.
- Ex. “You are responsible for validating scientific names before they are entered into the NCBI taxonomy database.”
    - By telling the LLM that we are entering these scientific names into the NCBI database, it is able to comprehend the types of results we are looking for. We are looking for more advanced scientific terms rather than common names.
### Don’t assume that giving more detailed directions is better.
- The longer prompt you give the LLM, the worse the returned results.
- You need to find a balance between adding details and prompt length.
- Repetitive instructions tend to also make the prompt harder to follow.
    - While you want to be sure to emphasize which directions in the prompt are important, restating them multiple times will not automatically improve the results.
- At the same time, removing too many details can have a negative impact.
- The goal is not to have the shortest possible prompt, but rather the shortest prompt with necessary instructions.
### Provide examples.
- Helps the LLM understand what you’re asking.
- Give both correct and incorrect result examples.
- Don’t give an extended list of examples, this will just lengthen the prompt with no benefit.
- I found that the best results were returned when they included examples for difficult names, like genus switches or taxonomy reclassifications.
    - Examples demonstrating when to return “No Match Found” were especially helpful
### Enforce when to return “No Match Found”.
- This was the hardest aspect for the LLM to understand.
- Oftentimes, the LLM would produce an answer that’s incorrect, solely because it was the best possible answer.
- I found that the best way to limit this was to include details on when to return “No Match Found”, but also that this is a valid result. The LLM is not failing by returning this.
- It’s better for the LLM to be overly cautious and return less correct results and more “No Match Found”, than return more correct results with more incorrect ones.
### Give the LLM instructions on what not to do.
- Even if you give the best possible procedure for the LLM to follow, it may still produce answers you don’t want.
- By telling the LLM explicitly not to do something, the results will be more accurate.
    - Ex. “Do not invent a scientific name.”
    - While it seems like this should be self explanatory, the LLM won’t know not to do this unless you tell it.
- Even the most obvious directions may need to be included.
### LLMs are not good at returning NCBI taxonomy IDs.
- Initially, this project had the LLM return both the scientific name and taxonomy ID. The LLM was decently good at finding the scientific name but very bad at finding the correct taxonomy ID.
    - Most of the time, the LLM would hallucinate an ID, even after explicitly telling it not to invent IDs.
- Because of this, we altered the workflow to include the Bio.Entrez package to retrieve the ID from the NCBI database that corresponds to the scientific name returned by the LLM.
- For different projects, methods like these may need to be implemented depending on what the LLM is good at, however you will not come to these conclusions without first testing with prompt engineering.
### API model type will have an impact on the results.
- Even with the same prompt, using different models will inherently produce different results.
- More expensive models will produce better results
- Because of this, we decided on using GPT-4o throughout prompt development to lower costs, and then switch to GPT-5.5 to get the final results

## Prompt Design Principles
### Use direct language, not vague/unambiguous instructions.
- Ex. Do not tell the LLM to “use your best judgement”. This is very vague and will often result in worse results.
### Be specific and include definitions
- What counts as an accepted scientific name
- What counts as an alternative name
- What constitutes a “No Match Found” result
- Etc.
### Give context
- By including a role and/or an objective, the LLM will give results more tailored to its task.
### Order the instructions within sections
- Role
- Objective
- Directions on multiple inputs
- Procedure
- Examples
- Etc.
### Be specific when requesting an output format
- Include the output format desired
- Give examples of outputs
### Be consistent with terminology
- When referring to the NCBI Taxonomy Database, write the full name. If you write another name, like only writing ‘the database’, even though it is referring to the same thing, the LLM may get confused and return worse results.

## Recommended Prompt Patterns
### Define a role
- “You are responsible for…”
- Gives the LLM important context
### Ordered procedure section
- Write the procedure in step by step format, rather than general definitions
### Examples
- Show correct results, incorrect results, and expected output
- Especially important if you show examples of names that are more difficult to map
### Strict directions on returning “No Match Found”
- Make clear that false positives are worse than not returning a name
### Output
- Explain the type of output desired
- Give exact structure with an example

## Common Anti-Patterns
### Overloading the prompt
- Adding extra directions, repeating the same rules, multiple sections saying the same thing, etc. will not result in better results
- More instructions does not equal better results
### Giving contradictory instructions
- If one section says something than contradicts another, the model has conflicting directions and may return worse results
    - Ex. If one section says “return the best possible name” but then another is enforcing the LLM to return “No Match Found” when the name cannot be established.
    - These conflict each other and may confuse the LLM on what to return
### Being too specific
- Don’t attempt to include every possible case.
- If every hard name is attempted to be defined within the prompt, it will be extremely long and inefficient
- Define the general rules, and then use examples to show the especially hard names.
### Asking the LLM to generate taxon IDs
- No matter what directions were given, the LLM oftentimes hallucinated these
- IDs could be returned with Bio.Entrez, so the workflow was altered
- Don’t be afraid to add a package to the script because it changes the workflow. This drastically helped improve results, so similar changes may be needed in the future.
### Assuming the LLM will return “No Match Found”
- If you don’t explicitly tell the model to do this, it will not.
- It will continue to return the best possible result unless specified otherwise
### Changing too much at once
- This affects how well you are able to evaluate the results
- If you change too much at once, it will be hard to tell which added/removed sections affected the results.
- Incremental changes will make evaluating the results and deciding what changes to make to the prompt, easier
