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

    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    eff_cat <- ncol(M)
    ## Transition probabilities of zero on this matrix can cause problems on log space.
    ## Here trying to just set them to a very small number.
    M[ round(M, digits = 10) == 0.0 ] <- .Machine$double.eps
    
    gamma.rates <- gamma.rates[ (k+1-eff_cat):k ]
    
    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) seqtraits:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states[site], edge_len = edge_len, edge_mat = edge_mat, parents = parents, X = X[[site]], Q = r*Q[[site]], root_node = root_node, root_type = root_type) ), mc.cores = n.cores )

    ## Store number of sites
    nsites <- length(X)

    ## Need special code for nsites < 5.
    if( nsites > 4 ){
        ## This will work for 5 or more sites.

        ## Make matrix to store the marginals for each rate at each site:
        marg <- matrix(NA, nrow = eff_cat, ncol = nsites)
        
        ## Marginal for the last site.
        for( cat in 1:eff_cat ){
            N_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
            lik_unit <- N_unit ## Start the loop.
            for( site in (nsites-1):2){ ## Next to last up to the second site.
                ## The 'getIntUnit" function already deals with the probs of 0.0
                lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=site
                                                                   , left_unit=lik_unit))
            }
            marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for the second to last site. [Will be the same for all middle sites.]
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites-1]][cat] + lik_unit[cat])
            for( site in (nsites-2):2){ ## Need to loop over all the rest.
                lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=site
                                                                   , left_unit=lik_unit))
            }
            marg[cat, nsites-1] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## A loop from site nsites-2 to site 3.
        for( marg.site in (nsites-2):3 ){ ## The big loop for the marginal.
            for( cat in 1:eff_cat ){
                lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
                for( site in (nsites-1):marg.site+1){ ## The sites to the right of the focus site.
                    lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                                       , n=site, left_unit=lik_unit))
                }
                ## Compute the marginal for the focus site:
                lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[marg.site]][cat] + lik_unit[cat])
                for( site in (marg.site-1):2){ ## The sites to the left of the focus site.
                    lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                                       , n=site, left_unit=lik_unit))
                }
                marg[cat, marg.site] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
            }
        }

        ## Marginal for site 2.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            for( site in (nsites-1):3){ ## The sites to the right of the focus site.
                lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                                   , n=site, left_unit=lik_unit))
            }
            ## Compute the marginal for the focus site:
            lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
            marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 1.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            for( site in (nsites-1):2){ ## The sites to the right of the focus site.
                lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                                   , n=site, left_unit=lik_unit))
            }
            ## Compute the marginal for the focus site:
            marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
        }
    }
    if( nsites == 4){
        ## Special case with 4 sites in the data.
        
        ## Make matrix to store the marginals for each rate at each site:
        marg <- matrix(NA, nrow = eff_cat, ncol = nsites)
        
        ## Marginal for the last site.
        for( cat in 1:eff_cat ){
            N_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
            lik_unit <- N_unit ## Start the loop.
            ## For site 2:
            lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                               , left_unit=lik_unit))
            ## For site 1:
            marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 3
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            ## Compute the marginal for the focus site (site 3):
            lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[3]][cat] + lik_unit[cat])
            ## The sites to the left of the focus site (1 and 2):
            lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                               , left_unit=lik_unit))
            marg[cat, 3] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 2 (i.e., second to last site).
        ## In this case we only have a single middle site.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
            marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 1.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            ## The sites to the right of the focus site (i.e., site 2).
            lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                               , left_unit=lik_unit))
            ## Compute the marginal for the focus site:
            marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
        }
    }
    if( nsites == 3){
        ## Special case with 3 sites in the data.
        
        ## Make matrix to store the marginals for each rate at each site:
        marg <- matrix(NA, nrow = eff_cat, ncol = nsites)
        
        ## Marginal for the last site.
        for( cat in 1:eff_cat ){
            N_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
            lik_unit <- N_unit ## Start the loop.
            ## For site 2:
            lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                               , left_unit=lik_unit))
            ## For site 1:
            marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 2 (i.e., second to last site).
        ## In this case we only have a single middle site.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
            marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 1.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            ## The sites to the right of the focus site (i.e., site 2).
            lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                               , left_unit=lik_unit))
            ## Compute the marginal for the focus site:
            marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
        }
    }
    if( nsites == 2){
        ## Special case with 2 sites in the data.
        
        ## Make matrix to store the marginals for each rate at each site:
        marg <- matrix(NA, nrow = k, ncol = nsites)
        
        ## Marginal for the last site.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
            ## For site 1:
            marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        }

        ## Marginal for site 1.
        for( cat in 1:eff_cat ){
            lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
            ## Compute the marginal for the focus site:
            marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
        }
    }

    ## To return the rates we need to compute a weighted average of the Gamma rates following the proportion from the marginal likelihood for each of the hidden rates.
    ## This is equal to exp( log( exp(x) / sum( exp(x) ) ) ), but we cannot compute the exponential of the logLikelihood because of underflow issues.
    rel.marg.lik <- apply( marg, 2, function(x) exp(x - logSumExp(x)) )
    real.Q <- lapply(1:length(X), function(x) sum(rel.marg.lik[,x] * gamma.rates) * Q[[x]] )
    
    ## Return the matrix of marginal probabilities by now.
    return( real.Q )
}
