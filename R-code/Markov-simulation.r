

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
