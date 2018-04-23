findMaxBeta <- function(ncats){
    ## Function to find the max realized value for beta.
    ## Dependent on the number of categories beta has no effect if larger than a threshould value.
    beta <- 0
    gamma.vec <- rep(0, times=ncats)
    while( gamma.vec[1] < ncats ){
        beta <- beta+1
        gamma.vec <- discreteGamma(shape=beta, ncats=ncats)
    }
    return( beta )
}
