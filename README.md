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

## Main results

The project compares deterministic model incidence with observed Omicron case
incidence. It also investigates how vaccination rates affect epidemic peaks and
uses repeated stochastic simulations to quantify model uncertainty.

## Repository contents

- `project.qmd` — complete Quarto analysis
- `project.pdf` — rendered report
- `omicron.csv` — data used in the analysis

## Reproducing the analysis

1. Download or clone the repository.
2. Open the Quarto file in RStudio.
3. Install the required R packages.
4. Render the Quarto document.

## Limitations

This is a simplified compartmental model. It assumes homogeneous mixing and
constant model parameters, and it does not account for demographic structure,
reporting delays or changing public-health behaviour.
