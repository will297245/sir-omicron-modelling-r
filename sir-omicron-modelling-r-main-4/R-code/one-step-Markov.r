
one_step_markov <- function(S, I, R, V, beta, gamma, mu, N) {
  
  # Infection and vaccination hazards
  lambda_inf <- beta * I / N
  lambda_vac <- mu
  
  # Total hazard out of S
  lambda_total <- lambda_inf + lambda_vac
  
  # Convert hazards into one-day transition probabilities
  if (lambda_total > 0) {
    p_leave_S <- 1 - exp(-lambda_total)
    
    p_SI <- (lambda_inf / lambda_total) * p_leave_S
    p_SV <- (lambda_vac / lambda_total) * p_leave_S
  } else {
    p_SI <- 0
    p_SV <- 0
  }
  
  p_stay <- 1 - p_SI - p_SV
  
  # Removal probability
  p_IR <- 1 - exp(-gamma)
  
  # Random transitions from S
  S_transitions <- as.vector(
    rmultinom(
      n = 1,
      size = S,
      prob = c(p_SI, p_SV, p_stay)
    )
  )
  
  N_SI <- S_transitions[1]  # new infections
  N_SV <- S_transitions[2]  # new vaccinations
  
  # Random recoveries from I
  N_IR <- rbinom(
    n = 1,
    size = I,
    prob = p_IR
  )
  
  # Update compartments
  S_new <- S - N_SI - N_SV
  I_new <- I + N_SI - N_IR
  R_new <- R + N_IR
  V_new <- V + N_SV
  
  # Return new state and incidence
  return(
    list(
      S = S_new,
      I = I_new,
      R = R_new,
      V = V_new,
      incidence = N_SI,
      new_vaccinations = N_SV,
      new_recoveries = N_IR
    )
  )

  
}
