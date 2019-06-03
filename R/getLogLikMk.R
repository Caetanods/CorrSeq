##' Internal function to get the likelihood of the model.
##'
##' No details.
##' @title Compute the loglik for the model.
##' @param phy the phylogeny
##' @param data the data matrix
##' @param Q the transition matrix
##' @param root.type the root type
##' @return Return the loglik for the model.
##' @author daniel
getLogLikMk <- function(phy, data, Q, root.type){
    ## Will get the data as a named vector and transform it to the correct format prior to the computation of the likelihood.
    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)
    dt <- make.data.tips(data)
    ## Check if Q is a numeric. If correct, then build the Q matrix.
    if( is.matrix(Q) ){
        Q.mat <- Q
    } else{
        if( !is.numeric(Q) ) stop( " Q needs to be a matrix or a single number. " )
        Q.mat <- matrix(Q, nrow=ncol(dt), ncol=ncol(dt))
        diag(Q.mat) <- -(colSums(Q.mat) - Q)
    }
    
    lik <- logLikMk(phy = phy, X = dt, Q = Q.mat, root.type = root.type)
    return( lik )
}
