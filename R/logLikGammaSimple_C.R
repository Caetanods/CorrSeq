##' The log-likelihood function for using the Gamma distribution.
##'
##' Computes the log-likelihood for the model.
##' @title Log-likelihood for the markov model using a Gamma-distributed rate heterogeneity.
##' @param n_nodes number of nodes, numeric
##' @param n_tips number of tips, numeric
##' @param n_states number of states in the site, a vector
##' @param edge_len vector with the length of the branches
##' @param edge_mat the matrix with ancestral and descendant nodes
##' @param parents a vector with the order for the ancestral nodes in edge_mat
##' @param root_node the root node in the matrix
##' @param X data, this is a list of matrices
##' @param Q a list of transition matrices 
##' @param beta 
##' @param k 
##' @param n.cores 
##' @return The log-likelihood for the model.
##' @author daniel
##' @noRd
loglikGammaSimple_C <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, X, Q, root_type, beta, k, n.cores){
    ## phy = phylo
    ## X = data frame with number of columns equal to the number of traits and rownames equal to the species names in the phylogeny.
    ## Q = transition matrix among the states. Here we use a global transition matrix, meaning that all the traits have the same number of states.
    ## pi = data frame with the probabilities at the root for each of the states for each of the traits.
    ## beta = the shape parameter for the Gamma distribution used to model among site rate variation.
    ## k = the parameter for the number of rate categories.
    ## This function assumes that all sites have the same number of states. This is our simplest case.
    
    ## For each site we need to sum the probabilities with the Q scale by the r factor from the gamma distribution.
    ## scale = beta and shape = alpha -> Comparing between the function and Yang papers.
    ## Here we follow Yang, 1993 and set beta = alpha, such that the mean of the distribution is equal to 1.
    gamma.rates <- discreteGamma(shape = beta, ncats = k)
    ## Need to protect if any of the rates has 0 value. Showing warning message here.
    ## Rate of 0 will set the likelihood to 0. So we can just skip it.
    gamma.rates <- gamma.rates[!round(gamma.rates, digits=20) == 0] ## The strict test with 0.0 will never trigger.
    ## The code here is parallel on the number of sites.
    gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) seqtraits:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states[site], edge_len = edge_len, edge_mat = edge_mat, parents = parents, X = X[[site]], Q = r*Q[[site]], root_node = root_node, root_type = root_type) )
                                  , mc.cores = n.cores )
    ## gamma.lik <- parallel::mclapply(1:length(X), function(site) sapply(gamma.rates, function(r) logLikMk(phy, X=X[[site]], Q=r*Q[[site]], root.type=root.type ) )
    ##                               , mc.cores = n.cores )
    ## We can use 'gamma.lik' to get the averaged transition matrix for the site.
    rel.lik <- lapply(1:length(X), function(x) exp(gamma.lik[[x]]) / sum( exp(gamma.lik[[x]]) ) )
    ## Q is a list of matrices with different dimensions.
    ## Note that this step is only needed if you need to get the matrices back.
    ## So we could just separate this into its own function.
    real.Q <- lapply(1:length(X), function(x) sum(rel.lik[[x]] * gamma.rates) * Q[[x]] )
    final.lik <- sum( sapply(1:length(X), function(x) log( sum( exp( gamma.lik[[x]] ) ) / k ) ) )
    
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final.lik ) ){
        ## print("Bad proposal! Rejecting")
        final.lik <- log(0)
    }

    ## Returns the loglik, the estimated Q per site and the relative likelihood for each of the rate categories on the site.
    return( list(log.lik=final.lik, Q=real.Q, rel.lik=rel.lik) )
}
