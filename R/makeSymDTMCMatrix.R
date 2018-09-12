## Function to do reject sampling of Doubly Stochastic Matrices.
## This will be updated with a smarter sampler soon.

makeSymDTMCMatrix <- function(size, pars){
    ## Function to generate a Discrete Time Markov Matrix from a vector of parameters.
    ## Will return an empty matrix if the matrix is ill-conditioned.
    
    mat <- matrix(0, nrow = size, ncol = size)
    pars.number <- size + sum( upper.tri(mat) )
    lower.id <- lower.tri(mat)
    
    start <- 1
    end <- size
    mat[1,1:size] <- pars[start:end]
    mat[1,1:size] <- mat[1,1:size] / sum(mat[1,1:size])
    start <- end+1
    end <- start+(size-2)
    for( i in 2:(size) ){
        mat[lower.id] <- t(mat)[lower.id] ## Update lower.tri
        mat[i,i:size] <- pars[start:end]
        mat[i,i:size] <- ( mat[i,i:size] * (1 - sum(mat[i,1:(i-1)])) ) / sum( mat[i,i:size] )
        start <- end+1
        end <- start+(size-i-1)
    }
    ## Return mat or return an empty matrix.
    ## If empty matrix then I can set logLik to -Inf.
    ## This will equivalent to rejecting a step into a parameter that is invalid.
    ## A improvement to this would be to find a way to sample this matrix without error.
    if( any( is.na( mat ) ) | any( is.nan( mat ) ) ){
        return( mat <- matrix(0, nrow = size, ncol = size) )
    }
    if( any( mat < 0 ) | any( apply(mat, 1, sum) + 1 != 2 ) ){
        return( mat <- matrix(0, nrow = size, ncol = size) )
    } else{    
        return(mat)
    }    
}
