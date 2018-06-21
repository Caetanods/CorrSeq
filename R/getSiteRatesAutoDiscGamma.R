##' Returns the relative likelihood and the estimated transition rates for the sites.
##'
##' Uses the MLE under the auto-discrete-Gamma distribution of rates among sites as described by Yang 1995 to estimate the rates for each site by computing the marginal likelihood across the sites.
##' @title Returns the relative Log-likelihood and estimated rates per site for the auto-discrete-Gamma.
##' @param phy phylogeny
##' @param X data
##' @param obj the vector of MLE for the parameters
##' @param model the type of model: "ER" or "DEL", used to build the Q matrices.
##' @param gap.key object needed in the case of a "DEL" model.
##' @param root.type root
##' @param k the number of Gamma discrete categories
##' @param n.cores the number of cores
##' @return The log-likelihood for the model.
##' @author daniel
##' @noRd
getSiteRatesAutoDiscGamma <- function(phy, X, obj, model, gap.key=NULL, root.type, k, n.cores){
    ## Here the rates are autocorrelated among the sites using the method proposed by Yang (1995).

    ## Reconstruction of the Q matrices depends on the model.
    if( model == "ER" ){
        make.Q.list <- function(rate, size){
            Q <- matrix(rate, nrow=size, ncol=size)
            diag(Q) <- -(colSums(Q) - rate)
            return(Q)
        }
        x <- exp(obj[1])
        beta <- obj[2]
        Q <- lapply(nstates, function(i) make.Q.list(rate=x[1], size=i) )
    }
    if( model == "DEL" ){
        if( is.null( gap.key ) ) stop("Object 'gap.key' was not found.")
        make.Q.list.DEL <- function(global.rate, del.rate, size, is.gap){
            ## The '-' state is always the first state in the matrix.
            Q <- matrix(global.rate, nrow=size, ncol=size)
            ## The separate transition is the Q[2+,1] and Q[1,2+]
            if( is.gap ){
                Q[1,] <- del.rate
                Q[,1] <- del.rate
            }
            diag(Q) <- sapply(1:size, function(x) -(sum(Q[x,]) - Q[x,x]))
            return(Q)
        }
        x.obs <- exp(obj[1])
        x.del <- exp(obj[2])
        beta <- obj[3]
        Q <- lapply(nstates, function(i) make.Q.list.DEL(global.rate=x.obs, del.rate=x.del, size=i
                                                       , is.gap=gap.key[i])
                    )
    }
    
    ## Get the rates from the beta parameter.
    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    if( any( gamma.rates <= 0 ) ) warning("MLE contains Gamma rate category <= 0.")

    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) logLikMk(phy, X=X[[site]], Q=r*Q[[site]], root.type=root.type ) )
                                  , mc.cores = n.cores )

    ## Define function to compute the site likelihood conditioned on the probability of that site.
    getSumLik <- function(gamma.lik, M, i, n, curr_sum){
        ## gamma.lik: the list of the site likelihoods.
        ## M: the probability matrix.
        ## i: the row id for M (rate category for site n-1)
        ## n: the current site. (note that this CANNOT be the first or the last site.
        ## curr_sum: the value for the cumulative recursive sum.
        ## This is a sum of probabilities, so we need to exponentiate.
        return( gamma.lik[[n]][j] + log( sum( sapply(1:ncol(M), function(j) M[i,j] * exp(curr_sum) ) ) ) )
    }

    ## Need to get the marginal estimate across hidden rates for each of the sites.
    ## For each site need to hold each rate class constant and compute the likelihood.
    ## Then the marginal estimate for the rate at that site will be the weighted average of rates across hidden classes.

    ## Computing for the last site:
    lik.marg <- matrix(nrow=k, ncol=k)
    for( hold.marg in 1:k){
        lik.sum <- rep(NA, times=k)
        for( i in 1:k ){
            lik.sum[i] <- gamma.lik[[length(X)]][hold.marg] ## The last site.
            ## Below is a recursion code. Going from the 'last-1' site to the first site.
            ## See likelihood formula in Yang (1995).
            for(n in (length(X)-1):1){
                lik.sum[i] <- getSumLik(gamma.lik=gamma.lik, M=M, i=i, n=n, curr_sum=lik.sum[i])
            }
        }
        lik.marg[hold.marg,] <- lik.sum
    }
    lik.marg.vec <- log( colSums( exp(lik.marg) ) )
    siteN_marg <- sum( sapply(1:k, function(i) lik.marg.vec[i] * gamma.rates[i] ) ) / sum( lik.marg.vec )

    ## Computing for the next to last site:
    ## Here it is more complicated because fixing a rate has an effect on the right and left element when computing the likelihood. The transition matrix for the rates need to be adapted.
    lik.marg <- matrix(nrow=k, ncol=k)
    for( hold.marg in 1:k){
        lik.sum <- rep(NA, times=k)
        for( i in 1:k ){
            lik.sum[i] <- gamma.lik[[length(X)]][i] ## The last site.
            ## Below is a recursion code. Going from the 'last-1' site to the first site.
            ## See likelihood formula in Yang (1995).
            for(n in (length(X)-1):1){
                lik.sum[i] <- getSumLik(gamma.lik=gamma.lik, M=M, i=i, n=n, curr_sum=lik.sum[i])
            }
        }
        lik.marg[hold.marg,] <- lik.sum
    }
    lik.marg.vec <- log( colSums( exp(lik.marg) ) )
    siteN_marg <- sum( sapply(1:k, function(i) lik.marg.vec[i] * gamma.rates[i] ) ) / sum( lik.marg.vec )
    
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final.lik ) ){
        print( "Bad proposal! Rejecting." )
        print( paste(c(c(Q[1,2]), beta), collapse="; ") )
        final.lik <- log(0)
    }

    ## This model will only return the logLikelihood.
    ## Need another funtion to compute the marginal for the Q matrices per site.
    return( final.lik )
}
