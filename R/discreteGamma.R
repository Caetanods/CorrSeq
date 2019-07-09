## This is the previous version of the function.
## Values are NOT the same as in Yang, 1994.
## discreteGamma <- function(shape, ncats){
##     quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
##     return(diff(c(0, pgamma(quantiles * shape, shape + 1, rate=shape), 1)) * ncats)
## }

## Corrected version of the function. This directly compute the integration of the Gamma density.
discreteGamma <- function(shape, ncats){
    ## Get the quantiles.    
    quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
    quantiles <- c(0, quantiles, Inf)

    ## Function to compute the average of the density function.
    meanGamma <- function(x){
        ## Using 'sapply' to make sure it works for vectors.
        res <- sapply(x, function(y) y * dgamma(y, shape = shape, rate = shape) )
        return( res )
    }

    rates <- vector( mode = "numeric", length = ncats )
    for( i in 1:ncats ){
        ## This is equation (10) in Yang (1994).
        int1 <- integrate(f = meanGamma, lower = quantiles[i], upper = quantiles[i+1])
        int2 <- integrate(f = dgamma, lower = quantiles[i], upper = quantiles[i+1], shape = shape, rate = shape)
        rates[i] <- int1$value / int2$value
    }

    return( rates )
}
