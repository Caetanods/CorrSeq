##' Defines a very important trick to compute the log( sum( exp( a ) ) ) to prevent underflow!
##'
##' Note that here we are not concerned with overflow.
##' @title Trick to compute the log(sum(exp(x))).
##' @param x a numeric vector in log
##' @return the quantity back in log scale
##' @author daniel
##' @noRd
#logSumExp <- function(x){
#    ## Sets the minimum quantity to be equal to 1e-15
#    cc <- log(1e-30) - min(x)
#    return( log( sum( exp( x + cc ) ) ) - cc )
#}

# May 28 2026
# Found a blog post that uses a different configuration. Trying this one.
# https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/
# This is the more robust implementation.
logSumExp <- function(x){
    cc <- max(x)
    return( cc + log( sum( exp( x - cc ) ) ) )
}

