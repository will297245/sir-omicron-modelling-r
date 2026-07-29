

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


set.seed(123)

M <- 5000
T_fit <- nrow(omicron_data)




many_sims <- simulate_many_markov(
  M = M,
  beta = beta_best_markov,
  gamma = 0.1,
  mu = 0,
  x0 = 1 - 0.004,
  y0 = 0.004,
  v0 = 0,
  T = T_fit,
  N_eff = N_eff
)




simulation_summary <- do.call(
  rbind,
  lapply(split(many_sims$incidence_prop, many_sims$t), function(z) {
    data.frame(
      mean = mean(z),
      lower = as.numeric(quantile(z, 0.025)),
      upper = as.numeric(quantile(z, 0.975))
    )
  })
)
