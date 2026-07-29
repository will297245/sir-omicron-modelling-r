calculate_incidence <- function(beta, gamma, x0, y0, T, h) {
  model_output <- simulate_sir_euler(beta, gamma, x0, y0, T, h)

  # Incidence is the flow of new infections over each interval [t, t + h].
  interval_incidence <- h * beta *
    model_output$x[-nrow(model_output)] *
    model_output$y[-nrow(model_output)]

  data.frame(
    t = model_output$t[-1],
    incidence = interval_incidence
  )
}
