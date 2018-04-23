uniformPrior <- function(x, min=0.00001, max=100){
    ## Returns 1 if within bounds or 0 otherwise.
    ## Works with vectors.
    return( log(as.numeric( all(min <= x & x <= max) ) ) )
}
