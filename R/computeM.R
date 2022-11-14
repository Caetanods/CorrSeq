computeM <- function(rate_cat, alpha, rho, k){
    ## This function will generate the transition probabilities matrix for the likelihood of the model.
    ## Need to have a matrix that makes sense in all possible cases.
    ## It is possible that the length of some of the categories is effectivelly 0. In this case computing the density under the curve will result in NaN. Need to find these cases and work with them.
    ## One easy way to do this would be to shrink the matrix and keep only rate categories with non-zero.
    ## Then we can get the size of the matrix when computing the likelihood. Note that we are skipping categories with zero width anyways...

    ## The first bound is always 0.0:
    rate_cat[1] <- 0.0
    
    ## Check the width of the categories:
    width_rate <- sapply(2:length(rate_cat), function(i) round(rate_cat[i] - rate_cat[i-1], digits = 10) )
    if( any(width_rate == 0.0) ){
        keep <- which(!width_rate == 0.0) ## Within the tolerance of the rounding.
        bounds_to_keep <- sapply(keep, function(x) c(x, x+1))
        ## bounds is a matrix, each column with a rate category to keep.
        if( ncol(bounds_to_keep) == 1 ){
            return( 1.0 ) ## A single category model. This will not be a matrix.
        }
        ## At least 2 categories. We can do a matrix with this.
        new_k <- ncol(bounds_to_keep)
        M <- matrix(nrow = new_k, ncol = new_k)
        for( i in 1:(new_k-1) ){
            for( j in 1:(new_k-1) ){
                ## Define the vector of observations we will need:
                rate_vec <- rate_cat[c(bounds_to_keep[,i], bounds_to_keep[,j])]
                ## Compute the square region:
                A <- pbivGamma(x=rate_vec[2], y=rate_vec[4], alpha=alpha, corr=rho, log.p=FALSE)
                B <- pbivGamma(x=rate_vec[1], y=rate_vec[3], alpha=alpha, corr=rho, log.p=FALSE)
                C <- pbivGamma(x=rate_vec[1], y=rate_vec[4], alpha=alpha, corr=rho, log.p=FALSE)
                D <- pbivGamma(x=rate_vec[2], y=rate_vec[3], alpha=alpha, corr=rho, log.p=FALSE)
                ## Some of the densities will be 0.0 (but the function will return NA.)
                if( is.na(A) ) A <- 0.0
                if( is.na(B) ) B <- 0.0
                if( is.na(C) ) C <- 0.0
                if( is.na(D) ) D <- 0.0
                marg_vec <- pgamma(q = rate_vec[1:2], shape = alpha, rate = alpha)
                M[i,j] <- (A + B - C - D) / (marg_vec[2] - marg_vec[1])
            }
        }
        ## We will always have an issue to compute the last category.
        for( i in 1:(new_k-1) ){
            M[new_k,i] <- 1 - sum(M[,i], na.rm = TRUE)
            M[i,new_k] <- M[new_k,i]
        }
        ## Need to finish the M matrix here:    
        M[new_k,new_k] <- 1 - sum(M[new_k,], na.rm = TRUE)
        
    } else{
        ## Compute the M matrix normally.    
        M <- matrix(nrow = k, ncol = k)
        for( i in 1:(k-1) ){
            for( j in 1:(k-1) ){
                ## Define the vector of observations we will need:
                rate_vec <- rate_cat[c(i, i+1, j, j+1)]
                ## Compute the square region:
                A <- pbivGamma(x=rate_vec[2], y=rate_vec[4], alpha=alpha, corr=rho, log.p=FALSE)
                B <- pbivGamma(x=rate_vec[1], y=rate_vec[3], alpha=alpha, corr=rho, log.p=FALSE)
                C <- pbivGamma(x=rate_vec[1], y=rate_vec[4], alpha=alpha, corr=rho, log.p=FALSE)
                D <- pbivGamma(x=rate_vec[2], y=rate_vec[3], alpha=alpha, corr=rho, log.p=FALSE)
                ## Some of the densities will be 0.0 (but the function will return NA.)
                if( is.na(A) ) A <- 0.0
                if( is.na(B) ) B <- 0.0
                if( is.na(C) ) C <- 0.0
                if( is.na(D) ) D <- 0.0
                marg_vec <- pgamma(q = rate_vec[1:2], shape = alpha, rate = alpha)
                M[i,j] <- (A + B - C - D) / (marg_vec[2] - marg_vec[1])
            }
        }
        ## We will always have an issue to compute the last category. This is because the last category is bounded at Inf.
        ## But we can use the attributes of the M matrix, which is doubly stochastic, to derive these probabilities:
        ## By the law of total probability, the last category will simply be 1 - the sum of prob of all others.
        for( i in 1:(k-1) ){
            M[k,i] <- 1 - sum(M[,i], na.rm = TRUE)
            M[i,k] <- M[k,i]
        }
        ## Need to finish the M matrix here:    
        M[k,k] <- 1 - sum(M[k,], na.rm = TRUE)
    }
    
    ## Hopefully this matrix is bullet-proof.
    return( M )
}
