
findmin_factory <- function(gamma, x0, y0, data){
    #we want to find minimum beta between beta_min, nbeta_max with N steps
    #we wish to return a function of beta_min, beta_max, N


    error_function <- factory_error(gamma, x0, y0, data)


    #will end up returning this funciton
    findmin <- function(beta_min, beta_max,N){

        errors_variable <- rep(0, N)


        #will find errors for each of these
        #then will return beta with lowest error
        beta_sequence <- seq(from = beta_min, to = beta_max, length.out = N)


        for (i in seq_along(beta_sequence)) {
            errors_variable[i] <- error_function(beta_sequence[i])
        }


        #finding where minimum beta value is, then returning it
        min_beta_position <- which.min(errors_variable)
        return(beta_sequence[min_beta_position])


    }
    
return(findmin)

}
