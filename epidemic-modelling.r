omicron_data <- read.csv("/Users/willlodder/Desktop/Math245/omicron.csv")

#print(head(omicron_data))
T <- 1000


#making a function which sets up a data frame given our starting params
initial_populations <- function(x0, y0, v0, T, N = 56490048) {
  
  starting_susceptible <- as.integer(round(N * x0))
  starting_infectious  <- as.integer(round(N * y0))
  starting_vaccinated  <- as.integer(round(N * v0))
  starting_removed     <- N - starting_susceptible - starting_infectious - starting_vaccinated
  
  populations_df <- data.frame(
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
  
  populations_df$S[1] <- starting_susceptible
  populations_df$I[1] <- starting_infectious
  populations_df$R[1] <- starting_removed
  populations_df$V[1] <- starting_vaccinated
  
  populations_df$x[1] <- populations_df$S[1] / N
  populations_df$y[1] <- populations_df$I[1] / N
  populations_df$r[1] <- populations_df$R[1] / N
  populations_df$v[1] <- populations_df$V[1] / N
  
  populations_df$incidence[1] <- 0
  
  return(populations_df)
}







#next step will be to implement the markov process and iterate over


one_step_markov <- function(S, I, R, V, beta, gamma, mu, N) {
  
  # 1. Infection and vaccination hazards
  lambda_inf <- beta * I / N
  lambda_vac <- mu
  
  # 2. Total hazard out of S
  lambda_total <- lambda_inf + lambda_vac
  
  # 3. Convert hazards into one-day transition probabilities
  if (lambda_total > 0) {
    p_leave_S <- 1 - exp(-lambda_total)
    
    p_SI <- (lambda_inf / lambda_total) * p_leave_S
    p_SV <- (lambda_vac / lambda_total) * p_leave_S
  } else {
    p_SI <- 0
    p_SV <- 0
  }
  
  p_stay <- 1 - p_SI - p_SV
  
  # 4. Removal probability
  p_IR <- 1 - exp(-gamma)
  
  # 5. Random transitions from S
  S_transitions <- as.vector(
    rmultinom(
      n = 1,
      size = S,
      prob = c(p_SI, p_SV, p_stay)
    )
  )
  
  N_SI <- S_transitions[1]  # new infections
  N_SV <- S_transitions[2]  # new vaccinations
  
  # 6. Random recoveries from I
  N_IR <- rbinom(
    n = 1,
    size = I,
    prob = p_IR
  )
  
  # 7. Update compartments
  S_new <- S - N_SI - N_SV
  I_new <- I + N_SI - N_IR
  R_new <- R + N_IR
  V_new <- V + N_SV
  
  # 8. Return new state and incidence
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

  #i would do df$columnname[i] and replace it 
}


#defining intital vectors as was in the project
df1 <- initial_populations(1-0.004, 0.004, 0, T = T)
print(head(df1))


for (i in 1:T){
    
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








T <- 1000

#df1 <- simulate_markov(
  #beta = 0.15,
  #gamma = 0.1,
  #mu = 0.001,
  #x0 = 1 - 0.004,
  #y0 = 0.004,
  #v0 = 0,
  #T = T,
  #N = 56490048
#)

#print(head(df1))


#plot(df1$t, df1$incidence / 56490048, type = "l")




library(ggplot2)

# Real UK population used for observed case proportions
N_UK <- 56490048

# Effective population size used inside stochastic Markov model
# Smaller N_eff = more visible stochastic variation
N_eff <- 50000

T <- 1000

set.seed(4)

df1 <- simulate_markov(
    
  beta = 0.13,
  gamma = 0.1,
  mu = 0.000,
  x0 = 1 - 0.004,
  y0 = 0.004,
  v0 = 0,
  T = T,
  N = N_eff
)

omicron_data <- read.csv("/Users/willlodder/Desktop/Math245/omicron.csv")
omicron_data$date <- as.Date(omicron_data$date)

# Remove initial model row because incidence[1] is artificial
model_incidence <- df1[-1, c("t", "incidence")]

n <- min(nrow(omicron_data), nrow(model_incidence))

comparison_df <- data.frame(
  date = omicron_data$date[1:n],
  observed_cases = omicron_data$cases[1:n],
  observed_prop = omicron_data$cases[1:n] / N_UK,
  model_incidence = model_incidence$incidence[1:n],
  model_prop = model_incidence$incidence[1:n] / N_eff
)

p1 <- ggplot(comparison_df, aes(x = date)) +
  geom_line(aes(y = observed_prop, linetype = "Observed daily cases / N")) +
  geom_line(aes(y = model_prop, linetype = "Markov model incidence / N_eff")) +
  labs(
    x = "Date",
    y = "Daily incidence proportion",
    title = "Observed Omicron Daily Cases vs Markov Model Incidence",
    linetype = "Series"
  ) +
  theme_minimal()




print(p1)

