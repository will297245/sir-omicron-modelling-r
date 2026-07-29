calculate_incidence_sse <- function(beta, gamma, x0, y0, data) {
  if (!all(c("time", "incidence") %in% names(data))) {
    stop("data must contain time and incidence columns.")
  }
  if (nrow(data) < 2) {
    stop("At least two observations are required.")
  }

  h <- data$time[2] - data$time[1]
  if (!is.finite(h) || h <= 0 || any(abs(diff(data$time) - h) > 1e-10)) {
    stop("Observation times must be equally spaced and increasing.")
  }

  observed_incidence <- data$incidence
  # One model transition is required for each observed incidence value.
  T <- length(observed_incidence) * h
  model_incidence <- calculate_incidence(
    beta = beta,
    gamma = gamma,
    x0 = x0,
    y0 = y0,
    T = T,
    h = h
  )$incidence

  if (length(model_incidence) != length(observed_incidence)) {
    stop("Model and observed incidence lengths differ.")
  }

  sum((model_incidence - observed_incidence)^2)
}

factory_error <- function(gamma, x0, y0, data) {
  function(beta) {
    calculate_incidence_sse(beta, gamma, x0, y0, data)
  }
}
