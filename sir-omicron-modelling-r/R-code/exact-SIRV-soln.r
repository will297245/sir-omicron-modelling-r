
solve_sirv_ode <- function(beta, gamma, mu, x0, y0, v0, T, h){

    state <- c(x = x0, y = y0, v = v0)
    parameters <- c(beta = beta , gamma = gamma , mu = mu )

    SIRV <- function(t, state, parameters) {

        with(as.list(c(state, parameters)), {
        dx <-  -beta * x * y - mu * x
        dy <-  beta * x * y - gamma * y
        dv <- mu * x
        list(c(dx, dy, dv))
    })
}

tim <- seq(0, T, by = h)



sol <- ode(y = state, times = tim, func = SIRV, parms = parameters)


return(sol)

}
