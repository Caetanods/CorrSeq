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
