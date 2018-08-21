##' The log-likelihood function with equal rates.
##'
##' Computes the log-likelihood for the model.
##' @title Log-likelihood for the markov model using a Gamma-distributed rate heterogeneity.
##' @param phy the phylogeny
##' @param X data frame with number of columns equal to the number of traits and rownames equal to the species names in the phylogeny.
##' @param Q transition matrix among the states. Here we use a global transition matrix, meaning that all the traits have the same number of states.
##' @param root.type the type of root
##' @param n.cores the cores to compute the lik for each site
##' @return The log-likelihood for the model.
##' @author daniel
##' @noRd
loglikSingleRate_C <- function(n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, X, Q, root_type, n.cores){
    ## Simple function computes the lik for a MK model for each site. Returns the joint probability.
    lik <- parallel::mclapply(1:length(X), function(site) seqtraits:::logLikMk_C(n_nodes = n_nodes, n_tips = n_tips, n_states = n_states[site], edge_len = edge_len, edge_mat = edge_mat, parents = parents, X = X[[site]], Q = Q[[site]], root_node = root_node, root_type = root_type)
                            , mc.cores = n.cores )
    ## lik <- parallel::mclapply(1:length(X), function(site) logLikMk(phy, X=X[[site]], Q=Q[[site]], root.type=root.type)
    ##                         , mc.cores = n.cores )
    ## Return the sum of the log likelihoods across the sites.
    final.lik <- do.call(sum, lik)
    
    ## A crappy wrap to catch cases in which the likelihood is just bad.
    ## This is because of bad proposals. Need to improve this.
    if( is.na( final.lik ) ){
        ## print("Bad proposal! Rejecting")
        final.lik <- log(0)
    }

    ## Return the log lik for the whole model
    return( final.lik )
}
