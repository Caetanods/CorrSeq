##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Find the maximum likelihood for the model with the gamma function.
##' @param data matrix with species names as rownames.
##' @param phy phylogeny
##' @param root.type one of "madfitz", or "equal". Need to extend to accept a observed vector of probabilities.
##' @param ncat categories for the gamma function
##' @param bounds 
##' @param opts the list of options for nloptr. If null it will use the default parameters.
##' @param n.cores 
##' @param bonds a numeric vector with the lower and upper bonds for the rates
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @export
##' @author daniel
fitMLGammaModel <- function(data, phy, root.type = "madfitz", ncat = 4, bounds = c(0,100), opts = NULL, n.cores = 1){
    ## At the moment the model assumes that all transitions have the same rate.
    ## So just a single rate is estimated for each of the transition matrices. Of course, this rate is changed by the gamma distribution.

    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)

    ## Make a function to call the nloptr function.
    wrapLogLikGammaSimple <- function(x, phy, data, nstates, k, root.type){
        ## x is a vector of parameters.
        ## x[1] = Q[1,2], x[2] = beta.
        ## phy = phylogeny
        ## data = data matrix.
        ## nstates = number of states in each of the sites in the data.
        Q <- matrix(x[1], nrow=nstates, ncol=nstates)
        diag(Q) <- -(colSums(Q) - x[1])
        beta <- x[2]
        ## Loglik function for the model.
        lik <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=k, it=1, n.cores=n.cores)[[1]]
        return( -lik ) ## Remember that NLOPT is minimizying the function!
    }

    ## We can find the maximum 'practical' value for beta conditioned in the number of categories. Beta values larger than this will not change the rate multiplier.
    maxBeta <- findMaxBeta(ncat)
    ## Create the vectors for the upper and lower bound for nloptr.
    if( bounds[1] < 0 ) stop( "The lower bound cannot be negative." )
    lb <- c(bounds[1], 0)
    ub <- c(bounds[2], maxBeta)
    ## Sample the initial parameters for the search.
    init.pars <- c( runif(1, min=lb[1], max=ub[1]), runif(1, min=lb[2], max=ub[2]) )

    ## Get the number of states. At the moment all the sites need to have the same number of states.
    ## No state can be missing in the tips.
    st.count <- apply(data, 2, function(x) length( unique(x) ) )
    if( any( !st.count == st.count[1] ) ) stop( "Some of the sites has a different number of elements." )
    ## Create the list of options or nloptr:
    if( is.null(opts) ){
        nlopt.opts <- list(algorithm="NLOPT_LN_SBPLX", "ftol_rel"=1e-08, "maxtime"=170000000, "maxeval"=10000)
    } else{
        if( !is.list(opts) ) stop( "The argument 'opts' needs to be a list format" )
        nlopt.opts <- opts
    }

    print( "Starting MLE search. " )
    fit <- nloptr(x0=init.pars, eval_f=function(x) wrapLogLikGammaSimple(x, phy=phy, data=data, nstates=st.count[1], k=ncat
                                                                       , root.type=root.type)
                , lb=lb, ub=ub, opts = nlopt.opts )

    print( "Reconstructing site-wise Q matrices." )
    ## After finishing the search. Need to reconstruct the Q matrices for each of the sites.
    Q <- matrix(fit$solution[1], nrow=st.count[1], ncol=st.count[1])
    diag(Q) <- -(colSums(Q) - fit$solution[1])
    beta <- fit$solution[2]
    res <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=ncat, it=1, n.cores=n.cores)[[2]]
    ## Here 'res' is a list with each of the Q matrices. The length of the list is equal to the number of sites in the data.
    out <- list( log.lik=-fit$objective, Q=res, start.par=init.pars, nlopt.opts=fit$options, nlopt.message=fit$message )
    return( out )
}

