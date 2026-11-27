# Disease Mapping: 43 Assigned Terms

**Contributor:** Shashank Kesineni

**Date:** November 26, 2026

**Branch:** disease-sk

**Related issue:** #44

## Objective

Investigate 43 assigned unmapped disease terms, identify potential Disease Ontology (DOID) mappings, and update the shared disease-mapping spreadsheet.

## Workflow

### 1. Prepare the assigned terms

Prepared a dataset containing the 43 assigned disease records, including:
- Original disease names
- ChatGPT-generated name suggestions
- Proposed OLS names and identifiers
- Validation status and notes

### 2. Search the Disease Ontology through OLS4

Queried the EBI OLS4 API using the Disease Ontology (`doid`).

Searches included:
- Exact searches
- Broader searches
- Label and synonym fields

Results were saved in `diseases/data/ols_direct_search_43.xlsx`.

### 3. Retrieve additional candidates

Searched using both original disease names and ChatGPT-generated name suggestions.

Retrieved and ranked candidate ontology terms using textual similarity. These scores were used for candidate discovery, not as proof of semantic equivalence.

### 4. GPT-assisted candidate review

Reviewed a subset of OLS candidates using GPT.

The review classified candidates as:
- Equivalent
- Broader
- Narrower
- Related but not equivalent
- Unrelated
- Uncertain

All GPT assessments were marked as requiring human verification.

## Results

### Direct OLS search

| Outcome | Count |
|---|---:|
| Disease terms searched | 43 |
| Candidate records returned | 96 |
| Diseases with candidates | 19 |
| Candidate-only results | 19 |
| No results | 20 |
| API errors | 4 |

All 43 records were marked for manual review.

### GPT-assisted review

| Assessment | Count |
|---|---:|
| Narrower | 14 |
| Equivalent | 3 |
| Unrelated | 3 |
| Related, not equivalent | 2 |
| Total candidates reviewed | 22 |
| Unique diseases represented | 8 |

All 22 assessments were marked as requiring human verification.

## Spreadsheet Completion

Completed the assigned 43 rows in the shared `Ongoing mappings` Google Sheet.

The spreadsheet work included reviewing proposed disease names and ontology mappings. The shared spreadsheet is the human-reviewed deliverable; the exploratory CSV and GPT outputs are not substitutes for it.

## Interpretation and Limitations

- OLS search results are candidate mappings, not confirmed mappings.
- Textual similarity does not establish disease equivalence.
- GPT assessments are AI-assisted suggestions and are not independent ground truth.
- Four direct searches encountered API errors and remain unresolved by that search.
- The direct search did not return candidates for 20 terms.
- This work did not constitute a controlled prompt benchmark.
- No before-and-after accuracy metrics or prompt acceptance decision were produced.

## Files

| File | Purpose |
|---|---|
| `diseases/search_ols_43.py` | Direct OLS search for the 43 assigned terms |
| `diseases/script/search_assigned_43_candidates.py` | Broader candidate discovery |
| `diseases/script/review_ols_candidates.py` | GPT-assisted candidate review |
| `diseases/data/ols_direct_search_43.xlsx` | Direct OLS search results |
| `diseases/data/assigned_43_spreadsheet_ready.csv` | Prepared 43-row working dataset |
| `diseases/data/assigned_43_gpt_validation_review.csv` | GPT-assisted review results |

## Follow-up

- Resolve the four API errors.
- Investigate terms with no OLS candidates using additional synonyms and alternative search terms.
- Verify proposed mappings against Disease Ontology definitions and disease scope.
- Continue human review of uncertain or potentially broader/narrower mappings.
