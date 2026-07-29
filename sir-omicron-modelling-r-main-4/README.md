# Epidemic Modelling in R

This project investigates deterministic and stochastic compartmental models of infectious-disease transmission using R and Quarto.

## Methods

- Deterministic SIR and SIRV models
- Euler numerical approximation and ODE solutions
- Incidence defined as susceptible-to-infectious transitions
- Least-squares fitting of the transmission parameter
- Discrete-time stochastic Markov modelling
- Monte Carlo simulation and pointwise simulation envelopes

## Requirements

Install R, Quarto, and the required R packages:

```r
install.packages(c("ggplot2", "deSolve"))
```

## Render

From the repository root:

```bash
quarto render epidemic-modelling.qmd --to html
```

The report is configured for HTML while debugging, avoiding a separate LaTeX dependency.

## Repository contents

- `epidemic-modelling.qmd` — complete analysis
- `R-code/` — model, fitting and simulation functions
- `data/omicron.csv` — observed daily case data
- `tests/smoke-tests.R` — lightweight validation tests

## Interpretation

The stochastic ribbon is a pointwise 95% simulation envelope conditional on fixed fitted parameters. It represents process variability, not a confidence interval and not parameter-estimation uncertainty. An effective population size of 50,000 is used to make stochastic variation visible; this is illustrative rather than calibrated to UK-wide uncertainty.

## Data provenance

Before treating the analysis as formal empirical evidence, add the original source URL, extraction date, geographic coverage and definition of the population denominator for `data/omicron.csv`.

## Limitations

The models assume homogeneous mixing and constant parameters. They do not model reporting delays, under-ascertainment, demographic structure, behavioural adaptation or time-varying policy interventions.
