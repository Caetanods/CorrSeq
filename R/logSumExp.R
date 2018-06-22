## 
## 
##' Defines a very important trick to compute the log( sum( exp( a ) ) ) to prevent underflow!
##'
##' Note that here we are not concerned with overflow.
##' @title Trick to compute the log(sum(exp(x))).
##' @param x a numeric vector in log
##' @return the quantity back in log scale
##' @author daniel
##' @noRd
logSumExp <- function(x){
    ## Sets the minimum quantity to be equal to 0.5
    cc <- log(0.5) - min(x)
    return( log( sum( exp( x + cc ) ) ) - cc )
}

