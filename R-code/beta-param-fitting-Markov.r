markov_expected_incidence <- function(beta, gamma, mu, x0, y0, v0, T) {
  
  # Store proportions, not integer counts
  x <- numeric(T + 1)
  y <- numeric(T + 1)
  r <- numeric(T + 1)
  v <- numeric(T + 1)
  incidence <- numeric(T + 1)
  
  # Initial conditions
  x[1] <- x0
  y[1] <- y0
  v[1] <- v0
  r[1] <- 1 - x0 - y0 - v0
  incidence[1] <- 0
  
  for (i in 1:T) {
    
    # Markov hazards for one susceptible individual
    lambda_inf <- beta * y[i]
    lambda_vac <- mu
    lambda_total <- lambda_inf + lambda_vac
    
    # One-day probabilities for S -> I and S -> V
    if (lambda_total > 0) {
      p_leave_S <- 1 - exp(-lambda_total)
      p_SI <- (lambda_inf / lambda_total) * p_leave_S
      p_SV <- (lambda_vac / lambda_total) * p_leave_S
    } else {
      p_SI <- 0
      p_SV <- 0
    }
    
    # One-day probability for I -> R
    p_IR <- 1 - exp(-gamma)
    
    # Expected proportions moving between compartments
    new_infections <- x[i] * p_SI
    new_vaccinations <- x[i] * p_SV
    new_recoveries <- y[i] * p_IR
    
    # Update proportions
    x[i + 1] <- x[i] - new_infections - new_vaccinations
    y[i + 1] <- y[i] + new_infections - new_recoveries
    r[i + 1] <- r[i] + new_recoveries
    v[i + 1] <- v[i] + new_vaccinations
    
    # Correct incidence: new infections, not y[i+1] - y[i]
    incidence[i + 1] <- new_infections
  }
  
  data.frame(
    t = 0:T,
    x = x,
    y = y,
    r = r,
    v = v,
    incidence = incidence
  )
}


markov_error_beta <- function(beta, data, gamma, mu, x0, y0, v0, N_UK) {
  
  # Observed daily case proportion
  observed_incidence <- data$cases / N_UK
  
  # We simulate one interval per observed daily case
  T_fit <- length(observed_incidence)
  
  model <- markov_expected_incidence(
    beta = beta,
    gamma = gamma,
    mu = mu,
    x0 = x0,
    y0 = y0,
    v0 = v0,
    T = T_fit
  )
  
  # Remove incidence[1], because it is artificial
  model_incidence <- model$incidence[-1]
  
  # Safety check
  n <- min(length(observed_incidence), length(model_incidence))
  
  sum((model_incidence[1:n] - observed_incidence[1:n])^2)
}

