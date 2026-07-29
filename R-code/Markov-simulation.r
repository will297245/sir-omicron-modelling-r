initial_populations <- function(
    x0,
    y0,
    v0,
    T,
    N = 56490048
) {
  starting_susceptible <- as.integer(round(N * x0))
  starting_infectious <- as.integer(round(N * y0))
  starting_vaccinated <- as.integer(round(N * v0))

  starting_removed <- N -
    starting_susceptible -
    starting_infectious -
    starting_vaccinated

  populations <- data.frame(
    t = 0:T,
    S = integer(T + 1),
    I = integer(T + 1),
    R = integer(T + 1),
    V = integer(T + 1),
    x = numeric(T + 1),
    y = numeric(T + 1),
    r = numeric(T + 1),
    v = numeric(T + 1),
    incidence = integer(T + 1)
  )

  populations$S[1] <- starting_susceptible
  populations$I[1] <- starting_infectious
  populations$R[1] <- starting_removed
  populations$V[1] <- starting_vaccinated

  populations$x[1] <- populations$S[1] / N
  populations$y[1] <- populations$I[1] / N
  populations$r[1] <- populations$R[1] / N
  populations$v[1] <- populations$V[1] / N

  populations
}



simulate_markov <- function(beta, gamma, mu, x0, y0, v0, T, N = 56490048) {
  
  df <- initial_populations(
    x0 = x0,
    y0 = y0,
    v0 = v0,
    T = T,
    N = N
  )
  
  for (i in 1:T) {
    
    step <- one_step_markov(
      S = df$S[i],
      I = df$I[i],
      R = df$R[i],
      V = df$V[i],
      beta = beta,
      gamma = gamma,
      mu = mu,
      N = N
    )
    
    df$S[i + 1] <- step$S
    df$I[i + 1] <- step$I
    df$R[i + 1] <- step$R
    df$V[i + 1] <- step$V
    
    df$incidence[i + 1] <- step$incidence
    
    df$x[i + 1] <- step$S / N
    df$y[i + 1] <- step$I / N
    df$r[i + 1] <- step$R / N
    df$v[i + 1] <- step$V / N
  }
  
  return(df)
}







simulate_many_markov <- function(M, beta, gamma, mu, x0, y0, v0, T, N_eff) {
  
  simulations <- vector("list", M)
  
  for (m in 1:M) {
    
    sim <- simulate_markov(
      beta = beta,
      gamma = gamma,
      mu = mu,
      x0 = x0,
      y0 = y0,
      v0 = v0,
      T = T,
      N = N_eff
    )
    
    # Remove first artificial incidence value
    sim <- sim[-1, ]
    
    simulations[[m]] <- data.frame(
      simulation = m,
      t = sim$t,
      incidence_prop = sim$incidence / N_eff
    )
  }
  
  do.call(rbind, simulations)
}





