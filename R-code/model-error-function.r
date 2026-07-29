calculate_incidence_sse <- function(beta, gamma, x0, y0, data) {
    #largest time value in data
    T <- data$time[nrow(data)]
    #assuming equal step size
    h <- data$time[2] - data$time[1]

    #incidence from data set
    dt <- data$incidence
    #estimated incidence from previous  function
    mt <- calculate_infectious_change(beta, gamma, x0, y0, T, h)$incidence

    #checking lengths of mt and dt are same
    if (length(mt) != length(dt)) {
        stop("Model output and data incidence have different lengths.")
    }

    sum((mt - dt)^2)
}



#will later use the below function to fit beta 

factory_error <- function(gamma, x0, y0, data) {
    output <- function(beta) {
        calculate_incidence_sse(beta, gamma, x0, y0, data)
    }

    return(output)
}
