slideWindowReflective <- function(x, w, max=1){
    ## Sliding window proposal for positive real trait.
    ## This is using the reflection at the boundary.
    ## In this case is as easy as getting the absolute value of the quantity.
    ## x = the current value.
    ## w = the width parameter of the proposal.
    y <- abs( runif(1, min = x - (w/2), max = x + (w/2) ) )
    if( y > 1 ){
        y <- (2*max) - y
    }    
    return(y)
}

multiplierProposal <- function(x, a){
    ## This proposal scheme will perform steps in the log space of the parameter.
    ## This is much more efficient for the cases in which small changes in large
    ##    values effect less than small changes in small values.
    ## Note this proposal will never change the sign of the parameter.
    ## x = the current value.
    ## a = the scale of the multiplier.
    lambda <- 2 * log(a)
    m <- exp( lambda * runif(1, min=-0.5, max=0.5) )
    ## Note that 'm' here is the proposal ratio. So need to spit this out.
    return( setNames(c(m * x, m), c("prop","prop.ratio") ) )
}

makePropQ.symm <- function(Q, scale = 2){
    ## Make proposal steps for Q. This assumes Q is symmetric.
    ## This is using the multiplier proposal.
    curr <- Q[1,2] ## All elements should be equal.
    new <- multiplierProposal(x=curr, a=scale)
    Q[] <- new[1]
    diag(Q) <- -(new[1]*(ncol(Q)-1))
    return( list(Q=Q, prop.ratio=new[2]) )
}

makePropPi <- function(pi, width = 0.5){
    curr <- pi[1,] ## Each column is a trait.
    new.0 <- sapply(curr, function(x) slideWindowReflective(x, width, max=1) )
    new.1 <- 1 - new.0
    pi.new <- rbind(new.0, new.1)
    return( pi.new )
}
