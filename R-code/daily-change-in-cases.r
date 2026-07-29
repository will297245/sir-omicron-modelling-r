calculate_infectious_change <- function(beta, gamma, x0, y0, T, h) {
    model_output <- simulate_sir_euler(beta, gamma, x0, y0, T, h)

    #adding 0 to the start of it to make the vector lengths natch up to t vector
    incidence_values <- c(0, diff(model_output$y))

    data.frame(
        t = model_output$t,
        infectious_change = incidence_values
    )
}
