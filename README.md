# Epidemic Modelling in R

This project investigates deterministic and stochastic models of infectious
disease transmission using R.

The analysis begins with a deterministic SIR model solved using Euler's method,
then extends the model to include vaccination. Model parameters are estimated
using observed Omicron case data. A stochastic Markov model is subsequently
used to generate simulation paths and uncertainty bands.

## Methods

- SIR and SIRV compartmental models
- Euler numerical approximation
- Ordinary differential equation solvers
- Least-squares parameter estimation
- Discrete-time Markov modelling
- Monte Carlo simulation
- Simulation-based uncertainty intervals

## Technologies

- R
- Quarto
- ggplot2
- deSolve

  
## Main findings

- Coarse Euler steps produced unstable and epidemiologically invalid trajectories.
- Increasing the vaccination-rate parameter reduced the simulated infectious peak.
- Treating the net change in infectious prevalence as incidence produced a poor
  model-data comparison.
- Redefining incidence as new susceptible-to-infectious transitions and refitting
  the model gave an estimated transmission parameter of approximately 0.136.
- Repeated stochastic simulations were used to construct pointwise 95% simulation
  envelopes.

## Repository contents

- `epidemic-modelling.qmd` — complete Quarto analysis
- `R-code` — the underlying r code document
- `data/omicron.csv` — data used in the analysis

## Reproducing the analysis

1. Download or clone the repository.
2. Open the Quarto file in RStudio.
3. Install the required R packages.
4. Render the Quarto document.

## Limitations

This is a simplified compartmental model. It assumes homogeneous mixing and
constant model parameters, and it does not account for demographic structure,
reporting delays or changing public-health behaviour.
