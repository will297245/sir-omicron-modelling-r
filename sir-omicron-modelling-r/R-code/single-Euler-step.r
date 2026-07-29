simulate_sir_euler <- function(beta, gamma, x0, y0, T, h) {
  if (any(!is.finite(c(beta, gamma, x0, y0, T, h)))) {
    stop("All arguments must be finite.")
  }
  if (beta < 0 || gamma < 0) {
    stop("beta and gamma must be non-negative.")
  }
  if (h <= 0 || T <= 0) {
    stop("T and h must be positive.")
  }
  if (x0 < 0 || y0 < 0 || x0 + y0 > 1) {
    stop("x0 and y0 must be non-negative and sum to at most one.")
  }

  n <- as.integer(round(T / h))
  if (!isTRUE(all.equal(n * h, T, tolerance = 1e-10))) {
    stop("h must divide T exactly.")
  }

  t <- seq(0, T, by = h)
  x <- numeric(n + 1)
  y <- numeric(n + 1)
  exact_solution <- numeric(n + 1)

  x[1] <- x0
  y[1] <- y0
  exact_solution[1] <- y0

  for (i in seq_len(n)) {
    x[i + 1] <- x[i] - h * beta * x[i] * y[i]
    y[i + 1] <- y[i] + h * (beta * x[i] * y[i] - gamma * y[i])
    exact_solution[i + 1] <- y0 * exp(-gamma * t[i + 1])
  }

  data.frame(t = t, x = x, y = y, exact_y_part_b = exact_solution)
}
