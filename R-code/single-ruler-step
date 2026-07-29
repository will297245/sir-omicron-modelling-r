
euler <- function(beta, gamma, x0, y0, T, h) {
    #defining our number of steps
    n <- as.integer(T / h)

    #checking that if all n as integers are equal, that they multiply to T
    #basically just cehcking h divides T nicely

    if (!isTRUE(all.equal(n * h, T))) {
        stop("h doesnt divide T exactly")
    }


    #defining some initial vectors whcih we will replace later
    t <- seq(0, T, by = h)
    x <- rep(0, n + 1)
    y <- rep(0, n + 1)
    exact_solution <- rep(0, n + 1)


    
    x[1] <- x0
    y[1] <- y0
    exact_solution[1] <- y0


    
    x_prime <- function(x_value, y_value) {
        -beta * x_value * y_value
    }

    y_prime <- function(x_value, y_value) {
        beta * x_value * y_value - gamma * y_value
    }


    #updating our vectors for each i
    for (i in 1:n) {
        x[i + 1] <- x[i] + h * x_prime(x[i], y[i])
        y[i + 1] <- y[i] + h * y_prime(x[i], y[i])
        exact_solution[i + 1] <- y0 * exp(-gamma * t[i + 1])
    }

    data.frame(t = t, x = x, y = y, exact_y_part_b = exact_solution)
}
