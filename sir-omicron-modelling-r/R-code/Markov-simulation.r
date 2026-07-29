initial_populations <- function(x0, y0, v0, T, N = 56490048) {
  if (length(T) != 1 || !is.finite(T) || T < 1 || T != as.integer(T)) {
    stop("T must be a positive integer.")
  }
  if (length(N) != 1 || !is.finite(N) || N < 1 || N != as.integer(N)) {
    stop("N must be a positive integer.")
  }
  if (any(!is.finite(c(x0, y0, v0))) ||
      any(c(x0, y0, v0) < 0) || x0 + y0 + v0 > 1) {
    stop("Initial proportions must be finite, non-negative and sum to at most one.")
  }

  T <- as.integer(T)
  N <- as.integer(N)
  starting_susceptible <- as.integer(round(N * x0))
  starting_infectious <- as.integer(round(N * y0))
  starting_vaccinated <- as.integer(round(N * v0))
  starting_removed <- N - starting_susceptible - starting_infectious - starting_vaccinated

  populations <- data.frame(
    t = 0:T,
    S = integer(T + 1), I = integer(T + 1),
    R = integer(T + 1), V = integer(T + 1),
    x = numeric(T + 1), y = numeric(T + 1),
    r = numeric(T + 1), v = numeric(T + 1),
    incidence = integer(T + 1)
  )

  populations[1, c("S", "I", "R", "V")] <- c(
    starting_susceptible, starting_infectious,
    starting_removed, starting_vaccinated
  )
  populations[1, c("x", "y", "r", "v")] <-
    populations[1, c("S", "I", "R", "V")] / N
  populations
}

simulate_markov <- function(beta, gamma, mu, x0, y0, v0, T, N = 56490048) {
  if (any(!is.finite(c(beta, gamma, mu))) || any(c(beta, gamma, mu) < 0)) {
    stop("beta, gamma and mu must be finite and non-negative.")
  }

  df <- initial_populations(x0 = x0, y0 = y0, v0 = v0, T = T, N = N)
  for (i in seq_len(T)) {
    step <- one_step_markov(
      S = df$S[i], I = df$I[i], R = df$R[i], V = df$V[i],
      beta = beta, gamma = gamma, mu = mu, N = N
    )
    df[i + 1, c("S", "I", "R", "V")] <- c(step$S, step$I, step$R, step$V)
    df$incidence[i + 1] <- step$incidence
    df[i + 1, c("x", "y", "r", "v")] <-
      df[i + 1, c("S", "I", "R", "V")] / N
  }
  df
}

simulate_many_markov <- function(M, beta, gamma, mu, x0, y0, v0, T, N_eff) {
  if (length(M) != 1 || !is.finite(M) || M < 1 || M != as.integer(M)) {
    stop("M must be a positive integer.")
  }
  simulations <- vector("list", M)
  for (m in seq_len(M)) {
    sim <- simulate_markov(
      beta = beta, gamma = gamma, mu = mu,
      x0 = x0, y0 = y0, v0 = v0, T = T, N = N_eff
    )[-1, ]
    simulations[[m]] <- data.frame(
      simulation = m,
      t = sim$t,
      incidence_prop = sim$incidence / N_eff
    )
  }
  do.call(rbind, simulations)
}
