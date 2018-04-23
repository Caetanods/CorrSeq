discreteGamma <- function(shape, ncats){
    quantiles <- qgamma((1:(ncats - 1))/ncats, shape = shape, rate = shape)
    return(diff(c(0, pgamma(quantiles * shape, shape + 1, rate=shape), 1)) * ncats)
}
