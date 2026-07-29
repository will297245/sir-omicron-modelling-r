function_factory_beta <- function(gamma, x0, y0, T, h) {
    function(beta) {
        incidence(beta, gamma, x0, y0, T, h)
    }
}

beta_variable_function <- function_factory_beta(0.2, 0.9, 0.1, 40, 0.1)
