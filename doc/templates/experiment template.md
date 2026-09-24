# Experiment ID

EXP-023

# Volunteer

John Smith

# Baseline Commit

[f5ed638](https://github.com/glygener/llm-metadata-mapping/commit/f5ed638a3eea2f2144fc79502e6cf6a462eb2046)

# Experiment Commit

[409b4d7](https://github.com/glygener/llm-metadata-mapping/commit/409b4d76e3bb76f05771f363f7bae53816e4ae36)

# Change Category

Constraint Addition

# Intent

Reduce incorrect ontology assignments caused by selection of broader or narrower ontology terms.

# Expected Outcome

Decrease Dangerous Error Rate.

# Benchmark Dataset

Curated benchmark dataset - 20 rows

# Results

## Before

Coverage: 92.1%

Ontology Coverage: 88.7%

Correct Mapping Rate: 81.3%

Dangerous Error Rate: 4.6%

Preferred Label Usage: 71.0%

Additional Correct Resolutions: 0

Resolution Regressions: 0

## After

Coverage: 90.3%

Ontology Coverage: 87.4%

Correct Mapping Rate: 83.2%

Dangerous Error Rate: 2.1%

Preferred Label Usage: 70.2%

Additional Correct Resolutions: 4

Resolution Regressions: 0

# Interpretation

The hierarchy restriction reduced incorrect ontology assignments while increasing the number of correctly resolved terms. No previously correct annotations were lost.

# Decision

ACCEPT

# Reason

Dangerous Error Rate improved substantially while introducing no Resolution Regressions.

# Follow-up Idea

Add examples for ambiguous disease names.