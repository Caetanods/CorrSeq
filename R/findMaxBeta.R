findMaxBeta <- function(ncats, precision = 0.01){
    ## Finds the lower and upper bounds for the beta parameter for the Gamma distribution.
    ## Prevents the generation of the 0 valued rates.
    
    lb <- 0.01 ## This will result in a very low min(rate). This is hard coded.
    
    ## The upper bound will depend on the ncat.
    ub <- lb
    new <- discreteGamma(shape=ub, ncats=ncats)[ncats]
    while( new > 0 ){
        ub <- ub + precision
        new <- discreteGamma(shape=ub, ncats=ncats)[ncats]
    }
    ub <- ub - precision ## One step back.
    bound <- setNames(c(lb, ub), c("lower","upper"))
    return( bound )
}

findMinBeta <- function(ncats, precision = 0.01){
    ## Finds the lower and upper bounds for the beta parameter for the Gamma distribution.
    ## Prevents the generation of the 0 valued rates.
    ## This function searches for the min value for the Beta parameter.
    ## The max value is set to 1000. Seems to be a good guess.

    init.lb <- 0.01 ## VERY low value!
    ## Search for the minimum value of Beta that has all categories larger than 0.
    min.rate <- .Machine$double.eps * 10 # One order of magnitude larger than the minimum double.
    new <- discreteGamma(shape=init.lb, ncats=ncats)
    while( any(new < min.rate) ){
        init.lb <- init.lb + precision
        new <- discreteGamma(shape=init.lb, ncats=ncats)
    }

    lb <- init.lb
    ub <- 1000 ## Hard coded, but could be an option.
    
    bound <- setNames(c(lb, ub), c("lower","upper"))
    return( bound )
}
