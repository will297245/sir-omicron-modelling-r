markov_expected_incidence <- function(beta, gamma, mu, x0, y0, v0, T) {
  if (any(!is.finite(c(beta, gamma, mu))) || any(c(beta, gamma, mu) < 0)) {
    stop("beta, gamma and mu must be finite and non-negative.")
  }
  if (T < 1 || T != as.integer(T)) stop("T must be a positive integer.")
  if (any(c(x0, y0, v0) < 0) || x0 + y0 + v0 > 1) {
    stop("Initial proportions must be non-negative and sum to at most one.")
  }

  x <- numeric(T + 1); y <- numeric(T + 1)
  r <- numeric(T + 1); v <- numeric(T + 1)
  incidence <- numeric(T + 1)
  x[1] <- x0; y[1] <- y0; v[1] <- v0
  r[1] <- 1 - x0 - y0 - v0

  for (i in seq_len(T)) {
    lambda_inf <- beta * y[i]
    lambda_vac <- mu
    lambda_total <- lambda_inf + lambda_vac
    if (lambda_total > 0) {
      p_leave_S <- 1 - exp(-lambda_total)
      p_SI <- (lambda_inf / lambda_total) * p_leave_S
      p_SV <- (lambda_vac / lambda_total) * p_leave_S
    } else {
      p_SI <- 0; p_SV <- 0
    }
    p_IR <- 1 - exp(-gamma)
    new_infections <- x[i] * p_SI
    new_vaccinations <- x[i] * p_SV
    new_recoveries <- y[i] * p_IR
    x[i + 1] <- x[i] - new_infections - new_vaccinations
    y[i + 1] <- y[i] + new_infections - new_recoveries
    r[i + 1] <- r[i] + new_recoveries
    v[i + 1] <- v[i] + new_vaccinations
    incidence[i + 1] <- new_infections
  }

  data.frame(t = 0:T, x = x, y = y, r = r, v = v, incidence = incidence)
}

markov_error_beta <- function(beta, data, gamma, mu, x0, y0, v0, N_UK) {
  observed_incidence <- data$cases / N_UK
  n_intervals <- length(observed_incidence)
  model_incidence <- markov_expected_incidence(
    beta = beta, gamma = gamma, mu = mu,
    x0 = x0, y0 = y0, v0 = v0, T = n_intervals
  )$incidence[-1]

  if (length(model_incidence) != length(observed_incidence)) {
    stop("Model and observed incidence lengths differ.")
  }
  sum((model_incidence - observed_incidence)^2)
}
