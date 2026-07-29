source("R-code/single-Euler-step.r")
source("R-code/daily-change-in-cases.r")
source("R-code/model-error-function.r")
source("R-code/one-step-Markov.r")
source("R-code/Markov-simulation.r")
source("R-code/beta-param-fitting-Markov.r")

sir <- simulate_sir_euler(0.15, 0.1, 0.996, 0.004, T = 10, h = 1)
stopifnot(nrow(sir) == 11, all(is.finite(unlist(sir))))

inc <- calculate_incidence(0.15, 0.1, 0.996, 0.004, T = 10, h = 1)
stopifnot(nrow(inc) == 10, all(inc$incidence >= 0))

set.seed(1)
markov <- simulate_markov(0.15, 0.1, 0, 0.996, 0.004, 0, T = 10, N = 50000)
stopifnot(nrow(markov) == 11)
stopifnot(all(markov$S + markov$I + markov$R + markov$V == 50000))
stopifnot(all(markov$incidence >= 0))

expected <- markov_expected_incidence(0.15, 0.1, 0, 0.996, 0.004, 0, T = 10)
stopifnot(nrow(expected) == 11, all(expected$incidence >= 0))

cat("All smoke tests passed.\n")
