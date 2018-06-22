getLastUnit <- function(gamma.lik, M, i, k, n){
    ## This function computes a minimum independent unit in the recursive algorithm.
    ## gamma.lik: the list of the site likelihoods.
    ## M: the probability matrix.
    ## i: the rate category of the n-1 neighbor site.
    ## n: the current site.
    return( logSumExp( sapply(1:k, function(j) log(M[i,j]) + gamma.lik[[n]][j] ) ) )
}

getIntUnit <- function(gamma.lik, M, i, k, n, left_unit){
    ## This function computes a minimum independent unit given the unit for the site on the left.
    ## gamma.lik: the list of the site likelihoods.
    ## M: the probability matrix.
    ## i: the rate category of the n-1 neighbor site.
    ## n: the current site.
    ## left_unit: the unit computed for the site on the left.
    return( logSumExp( sapply(1:k, function(j) log(M[i,j]) + gamma.lik[[n]][j] + left_unit[j] ) ) )
}

## And so on, until the first site:
getFirstUnit <- function(gamma.lik, k, second_unit){
    ## This function computes the final lik sum, given the accumulated lik across the sites.
    ## It is used only once, but it helps to read the code.
    ## gamma.lik: the list of the site likelihoods.
    ## second_unit: the unit computed for the site on the left.
    ## The probability of each rate category is the same here.
    ## Following Yang 1995 parameterization of the Gamma distribution.
    return( logSumExp( sapply(1:k, function(j) log(1/k) + gamma.lik[[1]][j] + second_unit[j] ) ) )
}

##' The log-likelihood function for using the auto-discrete-Gamma distribution.
##'
##' Computes the log-likelihood for the model using the auto-discrete-Gamma distribution of rates among sites as described by Yang 1995.
##' @title Log-likelihood for the markov model using a auto-discrete-Gamma among site rate heterogeneity.
##' @param phy phylogeny
##' @param X data
##' @param Q a list of transition matrices
##' @param M the probability matrix for the auto-correlation among sites.
##' @param root.type root
##' @param beta the parameter for the Gamma distribution
##' @param k the number of Gamma discrete categories
##' @param n.cores the number of cores
##' @return The log-likelihood for the model.
##' @author daniel
##' @noRd
logLikAutoDiscGamma <- function(phy, X, Q, M, root.type, beta, k, n.cores){
    ## phy = phylo
    ## X = data frame with number of columns equal to the number of traits and rownames equal to the species names in the phylogeny.
    ## Q = transition matrix among the states. Here we use a global transition matrix, meaning that all the traits have the same number of states.
    ## pi = data frame with the probabilities at the root for each of the states for each of the traits.
    ## beta = the shape parameter for the Gamma distribution used to model among site rate variation.
    ## k = the parameter for the number of rate categories.

    ## Here the rates are autocorrelated among the sites using the method proposed by Yang (1995).

    ## These rates CANNOT contain zero. If they are 0, then break and return error.
    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    if( any(gamma.rates == 0) ) stop("Rates scalers cannot contain 0 when the model is auto-correlated.")

    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) logLikMk(phy, X=X[[site]], Q=r*Q[[site]], root.type=root.type ) ), mc.cores = n.cores )

    ## Store number of sites
    nsites <- length(X)
    N_unit <- sapply(1:k, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=nsites)) ## Last site.

    lik_unit <- N_unit ## Start the loop.
    for( site in (nsites-1):2){ ## Next to last up to the second site.
        ## Loop will compute the cumulative probabilites across the sites using the recursive algorithm.
        lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k, n=site, left_unit=lik_unit))
    }

    ## The final sum is the likelihood for the model.    
    final_lik <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final_lik ) ){
        print( "Bad proposal! Rejecting." )
        final_lik <- log(0)
    }

    ## This model will only return the logLikelihood.
    ## Need another funtion to compute the marginal for the Q matrices per site.
    return( final_lik )
}
