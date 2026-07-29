# Edits applied

- Replaced net change in infectious prevalence with susceptible-to-infectious incidence in the deterministic model.
- Updated all deterministic fitting and plotting calls to use the corrected incidence function.
- Aligned one model transition with each observed daily incidence value.
- Replaced silent vector truncation with explicit length assertions.
- Added validation for model parameters, initial proportions, horizons and population sizes.
- Added a reproducible random seed for the 500-path stochastic simulation.
- Rebuilt the simulation summary without relying on row names.
- Corrected references to infectious prevalence and pointwise simulation envelopes.
- Added an explicit explanation of the illustrative effective population size.
- Configured Quarto to render HTML only while debugging.
- Updated the README with requirements, rendering instructions and limitations.
- Added `tests/smoke-tests.R` for basic functional validation.

## Local checks

Run from the repository root:

```bash
Rscript tests/smoke-tests.R
quarto render epidemic-modelling.qmd --to html
```

R and Quarto were not available in the editing environment, so the project was statically checked but not executed here.
