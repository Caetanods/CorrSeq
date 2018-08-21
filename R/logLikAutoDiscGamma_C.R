getLastUnit <- function(gamma.lik, M, i, k, n, keep){
    ## This function computes a minimum independent unit in the recursive algorithm.
    ## gamma.lik: the list of the site likelihoods.
    ## M: the probability matrix.
    ## i: the rate category of the n-1 neighbor site.
    ## n: the current site.
    ## keep: a matrix with elements to keep

    ## The keep matrix protects the likelihood when elements of M are 0.
    M.vec <- M[i,keep[i,]]
    lik.vec <- gamma.lik[[n]][keep[i,]]
    return( logSumExp( log(M.vec) + lik.vec ) )
}

getIntUnit <- function(gamma.lik, M, i, k, n, left_unit, keep){
    ## This function computes a minimum independent unit given the unit for the site on the left.
    ## gamma.lik: the list of the site likelihoods.
    ## M: the probability matrix.
    ## i: the rate category of the n-1 neighbor site.
    ## n: the current site.
    ## left_unit: the unit computed for the site on the left.
    M.vec <- M[i,keep[i,]]
    lik.vec <- gamma.lik[[n]][keep[i,]]
    return( logSumExp( log(M.vec) + lik.vec + left_unit[keep[i,]] ) )
}

## And so on, until the first site:
getFirstUnit <- function(gamma.lik, k, second_unit){
    ## This function computes the final lik sum, given the accumulated lik across the sites.
    ## It is used only once, but it helps to read the code.
    ## gamma.lik: the list of the site likelihoods.
    ## second_unit: the unit computed for the site on the left.
    ## The probability of each rate category is the same here.
    ## Following Yang 1995 parameterization of the Gamma distribution.
    return( logSumExp( log(1/k) + gamma.lik[[1]] + second_unit ) )
}

##' The log-likelihood function for using the auto-discrete-Gamma distribution.
##'
##' Computes the log-likelihood for the model using the auto-discrete-Gamma distribution of rates among sites as described by Yang 1995.
##' @title Log-likelihood for the markov model using a auto-discrete-Gamma among site rate heterogeneity.
##' @param n_nodes number of nodes, numeric
##' @param n_tips number of tips, numeric
##' @param n_states number of states in the site, a vector
##' @param edge_len vector with the length of the branches
##' @param edge_mat the matrix with ancestral and descendant nodes
##' @param parents a vector with the order for the ancestral nodes in edge_mat
##' @param root_node the root node in the matrix
##' @param X data, this is a list of matrices
##' @param Q a list of transition matrices
##' @param M the probability matrix for the auto-correlation among sites.
##' @param root_type a 0 for equal root probabilities and a 1 for madfitz.
##' @param beta the parameter for the Gamma distribution
##' @param k the number of Gamma discrete categories
##' @param n.cores the number of cores
##' @return The log-likelihood for the model.
##' @author daniel
##' @noRd
logLikAutoDiscGamma_C <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, X, Q, M, root_type, beta, k, n.cores){

    ## n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node: these are parameters derived from the phylogy.
    ## root_type: 0 = equal probabilities and 1: madfitz // Derived from the original argument that is a character type.
    
    ## Here the rates are autocorrelated among the sites using the method proposed by Yang (1995).

    ## These rates CANNOT contain zero. If they are 0, then break and return error.
    gamma.rates <- discreteGamma(shape = beta, ncats = k)

    ## Keep track of 0 probabilities in the M matrix:
    keep <- is.finite( log( M ) )
    
    ## This computes the likelihood for the sites given all the rate categories.
    gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) seqtraits:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states[site], edge_len = edge_len, edge_mat = edge_mat, parents = parents, X = X[[site]], Q = r*Q[[site]], root_node = root_node, root_type = root_type) ), mc.cores = n.cores )
    ## gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) logLikMk(phy, X=X[[site]], Q=r*Q[[site]], root.type=root.type ) ), mc.cores = n.cores )

    ## Store number of sites
    nsites <- length(X)
    N_unit <- sapply(1:k, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, k=k
                                                , n=nsites, keep=keep)) ## Last site.

    lik_unit <- N_unit ## Start the loop.
    for( site in (nsites-1):2){ ## Next to last up to the second site.
        ## Loop will compute the cumulative probabilites across the sites using the recursive algorithm.
        lik_unit <- sapply(1:k, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, k=k
                                                     , n=site, left_unit=lik_unit, keep=keep))
    }

    ## The final sum is the likelihood for the model.    
    final_lik <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
        
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final_lik ) ){
        ## id <- paste(sample(x=0:9, size = 5, replace = TRUE), collapse="")
        ## save(Q, M, beta, k, phy, X, file = paste0("error_",id,".RData"))
        ## print( paste0("Bad proposal ID: ", id) )
        final_lik <- log(0)
    }

    ## This model will only return the logLikelihood.
    ## Need another funtion to compute the marginal for the Q matrices per site.
    return( final_lik )
}
