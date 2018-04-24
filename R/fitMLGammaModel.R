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
##' @param init 
##' @param verbose 
##' @param n.cores 
##' @param bonds a numeric vector with the lower and upper bonds for the rates
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @export
##' @author daniel
fitMLGammaModel <- function(data, phy, root.type = "madfitz", ncat = 4, bounds = c(0,100), opts = NULL, init = NULL, verbose = TRUE, n.cores = 1){
    ## At the moment the model assumes that all transitions have the same rate.
    ## So just a single rate is estimated for each of the transition matrices. Of course, this rate is changed by the gamma distribution.

    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)

    ## We will search for the rate in log-space and for the beta value in normal space.
    ## Since beta has a upper and lower limit it becomes easy to traverse the whole space, especially using a global search together with the local search.
    wrapLogLikGammaSimple <- function(obj, phy, data, nstates, k, root.type){
        ## obj is a vector of 2 parameters.
        ## obj[1] = Q[1,2] (shared for all sites), obj[2] = beta.
        ## phy = phylogeny
        ## data = data matrix.
        ## nstates = number of states in each of the sites in the data.
        x <- exp(obj[1])
        beta <- obj[2]
        Q <- matrix(x[1], nrow=nstates, ncol=nstates)
        diag(Q) <- -(colSums(Q) - x[1])
        ## Loglik function for the model.
        lik <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=k, it=1, n.cores=n.cores)[[1]]
        return( -lik ) ## Remember that NLOPT is minimizying the function!
    }

    ## We can find the maximum 'practical' value for beta conditioned in the number of categories. Beta values larger than this will not change the rate multiplier.
    maxBeta <- findMaxBeta(ncat)
    ## Create the vectors for the upper and lower bound for nloptr.
    ## The bound is the same for each of the sites.
    if( bounds[1] < 0 ) stop( "The lower bound cannot be negative." )
    ## Log the bounds in order to search in log space.
    if( bounds[1] == 0 ){
        ## Cannot be log(0)
        bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
    }
    ## Make nlopt bounds vectors.
    log_lb <- c( log( bounds[1] ), 0 )
    log_ub <- c( log( bounds[2] ), maxBeta )

    ## Sample the initial parameters for the search.
    ## Here the user can provide a custom start.
    if( is.null(init) ){
        init.pars <- c(  log(runif(1, min=bounds[1], max=bounds[2]))
                       , runif(1, min=0, max=maxBeta  ) )
    } else{
        if( !length( init ) == 2 ) stop("Length of init need to be 2. init[1] is for the rate and init[2] is for the Gamma function parameter (beta).")
        if( any(init[1] < bounds[1]) | any(init[1] > bounds[2]) ) stop("Value for init[1] is out of bounds (defined by 'bounds').")
        if( any(init[2] < 0) | any(init[2] > maxBeta) ) stop( paste0("Value for beta (init[2]) is outside bounds. min = 0 and max = ", maxBeta,".") )
        init.pars <- c(log(init[1]), init[2])
    }

    ## Get the number of states. At the moment all the sites need to have the same number of states.
    ## No state can be missing in the tips.
    st.count <- apply(data, 2, function(x) length( unique(x) ) )
    if( any( !st.count == st.count[1] ) ) stop( "Some of the sites has a different number of elements." )
    ## Create the list of options for local search of nloptr:
    if( is.null(opts) ){
        ## nlopt.opts <- list(algorithm="NLOPT_LN_SBPLX", "ftol_rel"=1e-08, "maxtime"=170000000, "maxeval"=10000)
        global.opts <- list("algorithm"="NLOPT_GN_DIRECT_L_RAND_NOSCAL", "maxeval"=1000000, "ftol_rel"=.Machine$double.eps^0.5)
        local.opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"=1000000, "ftol_rel"=.Machine$double.eps^0.5)
    } else{
        if( !is.list(opts) ) stop( "The argument 'opts' needs to be a list format" )
        local.opts <- opts
        global.opts <- opts
        ## Set the algorithm for the global search.
        global.opts$algorithm <- "NLOPT_GN_DIRECT_L_RAND_NOSCAL"
    }

    if( verbose ){
        print( "Starting global MLE search. (First pass)" )
        global <- nloptr(x0=init.pars, eval_f=wrapLogLikGammaSimple, lb=log_lb, ub=log_ub, opts=global.opts
                       , phy=phy, data=data, nstates=st.count[1], k=ncat, root.type=root.type)
        print( "Global search solution:" )
        print( paste("Log-lik ", -global$objective
                   , "; rate ", exp(global$solution[1])
                   , "; alpha ", global$solution[2], collapse="") )
        print( "Starting local MLE search. (Second pass)" )
        fit <- nloptr(x0=global$solution, eval_f=wrapLogLikGammaSimple, lb=log_lb, ub=log_ub, opts=local.opts
                    , phy=phy, data=data, nstates=st.count[1], k=ncat, root.type=root.type)
        print( "Local search solution:" )
        print( paste("Log-lik ", -fit$objective
                   , "; rate ", exp(fit$solution[1])
                   , "; alpha ", fit$solution[2], collapse="") )
    }
    
    print( "Reconstructing site-wise Q matrices." )
    ## After finishing the search. Need to reconstruct the Q matrices for each of the sites.
    solution <- c(exp(fit$solution[1]), fit$solution[2])
    Q <- matrix(solution[1], nrow=st.count[1], ncol=st.count[1])
    diag(Q) <- -(colSums(Q) - solution[1])
    beta <- solution[2]
    res <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=ncat, it=1, n.cores=n.cores)[[2]]
    ## Here 'res' is a list with each of the Q matrices. The length of the list is equal to the number of sites in the data.
    out <- list( log.lik=-fit$objective, Q=res, global.rate=solution[1], alpha=solution[2]
              , start.par=c(exp(init.pars[1]),init.pars[2])
              , nlopt.global.search=global.opts, nlopt.local.search=local.opts
              , nlopt.message=fit$message )
    return( out )
}

