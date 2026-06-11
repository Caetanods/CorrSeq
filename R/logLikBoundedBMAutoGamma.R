logLikBoundedBMSimple <- function(lik_fn, bm, n.cores){
    ## Likelihood for the simple Bound BM across the sites of the sequence trait.
    ## Because this is not a Gamma model, we don't need to compute the weighted average at the end.
    bm_rate <- to_bound_BM_lik(bm)
    lik <- parallel::mclapply(1:length(lik_fn), function(site) lik_fn[[site]](x = bm_rate)
                            , mc.cores = n.cores )
    final.lik <- do.call(sum, lik)
    if( is.na( final.lik ) ){
        ## Return bad lik.
        return( Inf )
    }
    return(final.lik)
}

getMarginalBoundedBMSimpleGamma <- function(lik_fn, solution, k, n.cores){
    ## Function to estimate the marginal rate of the bounded BM model for each sequence position.
    ## This is to be used with the maximum likelihood value for the parameters and NOT part of the MLE estimation.
    gamma.rates <- discreteGamma(shape = solution[2], ncats = k)
    gamma.rates <- gamma.rates[!round(gamma.rates, digits=20) == 0] ## The strict test with 0.0 will never trigger.
    gamma_bm_rates <- gamma.rates * solution[1] # A vector with the scaled rates with length k
    scaled_lik_gamma_bm <- to_bound_BM_lik(gamma_bm_rates)
    gamma.lik <- lapply(1:length(lik_fn), function(site) sapply(scaled_lik_gamma_bm, function(r) lik_fn[[site]](x=r)))
    rel.lik <- lapply(1:length(lik_fn), function(x) exp(gamma.lik[[x]]) / sum( exp(gamma.lik[[x]]) ) )
    marginal.bm <- vector(mode = "numeric", length = length(lik_fn)) # Making sure the output is a vector
    for( i in 1:length(lik_fn) ){
        marginal.bm[i] <- sum(rel.lik[[i]] * gamma.rates) * solution[1]
    }
    ## Returning only the marginal estimate of the rates across the sites.
    return( marginal.bm )
}

logLikBoundedBMSimpleGamma <- function(lik_fn, bm, beta, k, n.cores){
    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    ## Need to protect if any of the rates has 0 value. Showing warning message here.
    ## Rate of 0 will set the likelihood to 0. So we can just skip it.
    gamma.rates <- gamma.rates[!round(gamma.rates, digits=20) == 0] ## The strict test with 0.0 will never trigger.
    gamma_bm_rates <- gamma.rates * bm # A vector with the scaled rates with length k
    scaled_lik_gamma_bm <- to_bound_BM_lik(gamma_bm_rates)
    ## The code here is parallel on the number of sites.
    gamma.lik <- parallel::mclapply(1:length(lik_fn), function(site) sapply(scaled_lik_gamma_bm, function(r) lik_fn[[site]](x=r)), mc.cores = n.cores)
    ## We can use 'gamma.lik' to get the averaged transition matrix for the site.
    rel.lik <- lapply(1:length(lik_fn), function(x) exp(gamma.lik[[x]]) / sum( exp(gamma.lik[[x]]) ) )
    final.lik <- sum( sapply(1:length(lik_fn), function(x) log( sum( exp( gamma.lik[[x]] ) ) / k ) ) )

    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final.lik ) ){
        ## print("Bad proposal! Rejecting")
        return( log(0) )
    } else{
        return( final.lik )
    }

    ## Returns the loglik, the estimated Q per site and the relative likelihood for each of the rate categories on the site.
    ## return( list(log.lik=final.lik, Q=real.Q, rel.lik=rel.lik) )
    ## UPDATE: We have a function for the marginal estimate: getMarginalBMSimpleGamma
}

logLikBoundedBMAutoGamma <- function(lik_fn, bm, M, beta, k, n.cores){
    ## lik_fn is the likelihood function for the model.
    ## M is the M matrix
    ## beta is the parameter for the Gamma
    ## k is the number of categories for the Gamma
    ## n.cores is the number of cores to run the likelihood computation
    
    ## Here the rates are autocorrelated among the sites using the method proposed by Yang (1995).
    gamma.rates <- discreteGamma(shape = beta, ncats = k)

    ## The likelihood for the model can be problematic to compute if transitions on the M matrix are absolute 0.0 . So here we can change the 0.0 values to a very small value.
    M[ round(M, digits = 10) == 0.0 ] <- .Machine$double.eps
    
    ## Need to match the rates with the dimension of the M matrix:
    ## We know that, because this is a Gamma distribution, that if rates are collapsed, they are the ones closer to 0.
    ## We are filtering before, so 'ncol(M)' is at least 2.
    effective_rates <- ncol(M)
    gamma.rates <- gamma.rates[ (k+1-effective_rates):k ]
    gamma_bm_rates <- gamma.rates * bm # A vector with the scaled rates with length k
    scaled_lik_gamma_bm <- to_bound_BM_lik(gamma_bm_rates)

    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(lik_fn), function(site) sapply(scaled_lik_gamma_bm, function(r) lik_fn[[site]](x = r) ), mc.cores = n.cores)

    ## Store number of sites
    nsites <- length(lik_fn)
    N_unit <- sapply(1:effective_rates, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites)) ## Last site.

    lik_unit <- N_unit ## Start the loop.
    for( site in (nsites-1):2){ ## Next to last up to the second site.
        ## Loop will compute the cumulative probabilites across the sites using the recursive algorithm.
        lik_unit <- sapply(1:effective_rates, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                     , n=site, left_unit=lik_unit))
    }

    ## The final sum is the likelihood for the model.    
    final_lik <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final_lik ) ){
        ## print( paste0("Bad proposal ID: ", id) )
        return( Inf )
    }

    ## This model will only return the logLikelihood.
    ## Need another funtion to compute the marginal for the Q matrices per site.
    return( final_lik )
}
