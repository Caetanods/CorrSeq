##' Returns the relative likelihood and the estimated transition rates for the sites.
##'
##' Uses the MLE under the auto-discrete-Gamma distribution of rates among sites as described by Yang 1995 to estimate the rates for each site by computing the marginal likelihood across the sites.
##' @title Returns the relative Log-likelihood and estimated rates per site for the auto-discrete-Gamma.
##' @param phy phylogeny
##' @param X the formated data
##' @param Q a list of Q matrices
##' @param M a matrix of transition probabilities
##' @param root.type root
##' @param beta the rate for the Gamma model
##' @param k the number of Gamma discrete categories
##' @param n.cores the number of cores
##' @return The log-likelihood for the model.
##' @author daniel
##' @noRd
getSiteRatesAutoDiscGamma <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, X, Q, M, root_type, beta, k, n.cores){
    ## Function will compute the marginal estimates for the rates per site.
    ## This is important for the analyses of the results.

    ## These rates CANNOT contain zero. If they are 0, then break and return error.
    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    if( any(gamma.rates == 0) ) stop("Rates scalers cannot contain 0 when the model is auto-correlated.")

    ## Keep track of 0 probabilities in the M matrix:
    keep <- is.finite( log( M ) )
    
    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) seqtraits:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states[site], edge_len = edge_len, edge_mat = edge_mat, parents = parents, X = X[[site]], Q = r*Q[[site]], root_node = root_node, root_type = root_type) ), mc.cores = n.cores )
    ## gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) logLikMk(phy, X=X[[site]], Q=r*Q[[site]], root.type=root.type ) ), mc.cores = n.cores )

    ## Store number of sites
    nsites <- length(X)

    ## Make matrix to store the marginals for each rate at each site:
    marg <- matrix(NA, nrow = k, ncol = nsites)
    
    ## Marginal for the last site.
    for( cat in 1:k ){
        N_unit <- sapply(1:k, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
        lik_unit <- N_unit ## Start the loop.
        for( site in (nsites-1):2){ ## Next to last up to the second site.
            lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit, keep = keep))
        }
        marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }

    ## Marginal for the second to last site. [Will be the same for all middle sites.]
    for( cat in 1:k ){
        lik_unit <- sapply(1:k, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=nsites, keep = keep))
        lik_unit <- sapply(1:k, function(i) log(M[i,cat]) + gamma.lik[[nsites-1]][cat] + lik_unit[cat])
        for( site in (nsites-2):2){ ## Need to loop over all the rest.
            lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit, keep = keep))
        }
        marg[cat, nsites-1] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }

    ## A loop from site nsites-2 to site 3.
    for( marg.site in (nsites-2):3 ){ ## The big loop for the marginal.
        for( cat in 1:k ){
            lik_unit <- sapply(1:k, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=nsites, keep = keep))
            for( site in (nsites-1):marg.site+1){ ## The sites to the right of the focus site.
                lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit, keep = keep))
            }
            ## Compute the marginal for the focus site:
            lik_unit <- sapply(1:k, function(i) log(M[i,cat]) + gamma.lik[[marg.site]][cat] + lik_unit[cat])
            for( site in (marg.site-1):2){ ## The sites to the left of the focus site.
                lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit, keep = keep))
            }
            marg[cat, marg.site] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }
    }

    ## Marginal for site 2.
    for( cat in 1:k ){
        lik_unit <- sapply(1:k, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=nsites, keep = keep))
        for( site in (nsites-1):3){ ## The sites to the right of the focus site.
            lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit, keep = keep))
        }
        ## Compute the marginal for the focus site:
        lik_unit <- sapply(1:k, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
        marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }

    ## Marginal for site 1.
    for( cat in 1:k ){
        lik_unit <- sapply(1:k, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=nsites, keep = keep))
        for( site in (nsites-1):2){ ## The sites to the right of the focus site.
            lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit, keep = keep))
        }
        ## Compute the marginal for the focus site:
        marg[cat, 1] <- log(1/k) + gamma.lik[[1]][cat] + lik_unit[cat]
    }

    ## To return the rates we need to compute a weighted average of the Gamma rates following the proportion from the marginal likelihood for each of the hidden rates.
    ## This is equal to exp( log( exp(x) / sum( exp(x) ) ) ), but we cannot compute the exponential of the logLikelihood because of underflow issues.
    rel.marg.lik <- apply( marg, 2, function(x) exp(x - logSumExp(x)) )
    real.Q <- lapply(1:length(X), function(x) sum(rel.marg.lik[,x] * gamma.rates) * Q[[x]] )
    
    ## Return the matrix of marginal probabilities by now.
    return( real.Q )
}
