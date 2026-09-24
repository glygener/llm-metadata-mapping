# Prompt Optimization Metrics

All percentage (%) metrics are calculated using the total number of terms in the benchmark test dataset as the denominator unless otherwise specified.


| Metric | Definition |
|---------|------------|
| Coverage | Percentage of benchmark terms for which the LLM proposes a scientific term. |
| No Scientific Term | Percentage of benchmark terms for which the LLM does not propose a scientific term. |
|---------|------------|
| Ontology Coverage | Percentage of benchmark terms for which the proposed scientific term can be resolved to an identifier in the target ontology. |
| Unresolvable Scientific Term | Percentage of benchmark terms for which the LLM proposes a scientific term that cannot be resolved to an identifier in the target ontology. |
|---------|------------|
| Correct Mapping Rate | Percentage of benchmark terms that resolve to the same ontology identifier as the human-curated annotation. |
| Dangerous Error Rate | Percentage of ontology-resolvable terms that resolve to an incorrect ontology identifier. Calculated as: Incorrect ontology mappings ÷ Ontology-resolvable terms. |
|---------|------------|
| Preferred Label Usage | Percentage of ontology-resolvable terms for which the LLM directly proposes the ontology preferred label rather than a synonym. |
| Additional Correct Resolutions | Number of terms that were not correctly resolved by the baseline prompt but are correctly resolved after the prompt modification. |
| Resolution Regressions | Number of terms that were correctly resolved by the baseline prompt but become unresolved or incorrectly resolved after the prompt modification. |


## Optimization Priorities
 
Prompt modifications should be evaluated using the following priority order:
 
1. Minimize Resolution Regressions.
2. Minimize Dangerous Error Rate.
3. Maximize Additional Correct Resolutions.
4. Maximize Correct Mapping Rate.
5. Maximize Ontology Coverage.
6. Maximize Coverage.
7. Maximize Preferred Label Usage.
 
A prompt modification that introduces Resolution Regressions or increases Dangerous Error Rate should generally be rejected unless there is a compelling project-specific justification.