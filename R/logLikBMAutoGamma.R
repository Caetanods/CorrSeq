logLikBMSimple <- function(lik_fn, sigma, mu, n.cores){
    ## Likelihood for the simple BM across the sites of the sequence trait.
    ## Because this is not a Gamma model, we don't need to compute the weighted average at the end.
    lik <- parallel::mclapply(1:length(lik_fn), function(site) lik_fn[[site]](sigma = sigma
                                                                              , mu = mu[site])
                            , mc.cores = n.cores )
    final.lik <- do.call(sum, lik)
    if( is.na( final.lik ) ){
        ## Return bad lik.
        return( log(0) )
    }
    return(final.lik)
}

logLikBMSimpleGamma <- function(lik_fn, sigma, mu, beta, k, n.cores){
    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    ## Need to protect if any of the rates has 0 value. Showing warning message here.
    ## Rate of 0 will set the likelihood to 0. So we can just skip it.
    ## Test below is safe for floating point numbers. Note that this will also exclude very small numbers.
    keep_gamma <- sapply(gamma.rates, function(x) !isTRUE(all.equal(x, 0)))
    gamma.rates <- gamma.rates[keep_gamma]
    gamma_bm_rates <- gamma.rates * sigma # A vector with the scaled rates with length k
    
    ## The code here is parallel on the number of sites.
    gamma.lik <- parallel::mclapply(1:length(lik_fn), function(site) sapply(gamma_bm_rates, function(r) lik_fn[[site]](sigma = r, mu = mu[site])), mc.cores = n.cores)
    ## We can use 'gamma.lik' to get the averaged transition matrix for the site.
    
    ## NOTE: The rel.lik is computed below but it is not used for anything. Why is this value computed?
    rel.lik <- lapply(1:length(lik_fn), function(x) exp(gamma.lik[[x]]) / sum( exp(gamma.lik[[x]]) ) )
    
    #final.lik <- sum( sapply(1:length(lik_fn), function(x) log( sum( exp( gamma.lik[[x]] ) ) / k ) ) )
    ## Re-writing to use the log-sum-exponential trick for underflow:
    final.lik <- sum( sapply(1:length(lik_fn), function(x) logSumExp( gamma.lik[[x]] ) - log(k) ) )
    
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final.lik ) ){
        ## print("Bad proposal! Rejecting")
        return( log(0) )
    } else{
        return( final.lik )
    }
}

logLikBMAutoGamma <- function(lik_fn, sigma, mu, M, beta, k, n.cores){
    ## lik_fn is the likelihood function for the model.
    ## M is the M matrix
    ## beta is the parameter for the Gamma
    ## k is the number of categories for the Gamma
    ## n.cores is the number of cores to run the likelihood computation
    
    ## Here the rates are autocorrelated among the sites using the method proposed by Yang (1995).
    gamma.rates <- discreteGamma(shape = beta, ncats = k)

    ## The likelihood for the model can be problematic to compute if transitions on the M matrix are absolute 0.0 . Set a minimum to .Machine$double.eps
    M[ round(M, digits = 10) < .Machine$double.eps ] <- .Machine$double.eps
    
    ## Need to match the rates with the dimension of the M matrix:
    ## We know that, because this is a Gamma distribution, that if rates are collapsed, they are the ones closer to 0.
    ## We are filtering before, so 'ncol(M)' is at least 2.
    effective_rates <- ncol(M)
    gamma.rates <- gamma.rates[ (k+1-effective_rates):k ]
    gamma_bm_rates <- gamma.rates * sigma # A vector with the scaled rates with length k

    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(lik_fn), function(site) sapply(gamma_bm_rates, function(r) lik_fn[[site]](sigma = r, mu = mu[site]) ), mc.cores = n.cores)

    ## Store number of sites
    nsites <- length(lik_fn)
    N_unit <- sapply(1:effective_rates, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i
                                                                , n=nsites)) ## Last site.

    lik_unit <- N_unit ## Start the loop.
    for( site in (nsites-1):2){ ## Next to last up to the second site.
        ## Loop will compute the cumulative probability across the sites using the recursive algorithm.
        lik_unit <- sapply(1:effective_rates, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                     , n=site, left_unit=lik_unit))
    }

    ## The final sum is the likelihood for the model.    
    final_lik <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final_lik ) ){
        ## print( paste0("Bad proposal ID: ", id) )
        return( log(0) )
    }

    ## This model will only return the logLikelihood.
    ## Need another function to compute the marginal the rate for each site.
    return( final_lik )
}
