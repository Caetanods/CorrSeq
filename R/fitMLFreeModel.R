##' .. content for \description{} (no empty lines) ..
##'
##' Here each of the sites in the data is fit a different MK model. So this is the most parameter rich model.
##' @title Fit the MLE for the free model.
##' @param data 
##' @param phy 
##' @param init 
##' @param root.type 
##' @param bounds 
##' @param opts 
##' @param n.cores 
##' @return A list with the fit, solution and other information.
##' @author daniel
##' @export
##' @importFrom nloptr nloptr
fitMLFreeModel <- function(data, phy, init=NULL, root.type = "madfitz", bounds = c(0,100), opts = NULL, n.cores = 1){
    ## This model does not use a function to model rate variation. This is a joint estimate for every site.
    ## This is the most complex model we can use here.
    
    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)

    ## Make a function to call the nloptr function.
    wrapLogLik <- function(x, phy, data.list, nstates, nsites, root.type){
        ## x is a vector of parameters.
        ## length of x is equal to the number of sites.
        ## phy = phylogeny
        ## data = data matrix.
        ## nstates = number of states in each of the sites in the data.
        Q <- list()
        for( i in 1:nsites){
            Q[[i]] <- matrix(x[i], nrow=nstates, ncol=nstates)
            diag(Q[[i]]) <- -(colSums(Q[[i]]) - x[i])
        }
        lik.site <- sapply(1:nsites, function(x) logLikMk(phy=phy, X=data.list[[x]], Q=Q[[x]], root.type))
        lik <- sum(lik.site)
        return( -lik ) ## Remember that NLOPT is minimizying the function!
    }

    ## Create a list of data with length equal to the number of sites.
    nsites <- ncol(data)
    data.list <- lapply(1:nsites, function(x) make.data.tips(setNames(as.numeric(data[,x]), rownames(data))) )

    ## Create the vectors for the upper and lower bound for nloptr.
    ## The bound is the same for each of the sites.
    if( bounds[1] < 0 ) stop( "The lower bound cannot be negative." )
    lb <- rep(bounds[1], times = nsites)
    ub <- rep(bounds[2], times = nsites)
    
    ## Sample the initial parameters for the search.
    ## Here the user can provide a custom start.
    if( is.null(init) ){
        init.pars <- runif(nsites, min=bounds[1], max=bounds[2])
    } else{
        if( any(init < bounds[1]) | any(init > bounds[2]) ) stop("Provided 'init' value is outside search bounds.")
        if( !length(init) == nsites ) stop(" Length of 'init' needs to be the same as the length of the sequence." )
        init.pars <- init
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

    print( "Starting global MLE search. (First pass)" )
    global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=lb, ub=ub, opts=global.opts
                   , phy=phy, data.list=data.list, nstates=st.count[1], root.type=root.type
                   , nsites=nsites)
    print( "Global search solution:" )
    print( paste("Log-lik:", -global$objective, "Pars:", global$solution, sep=" ") )
    print( "Starting local MLE search. (Second pass)" )
    fit <- nloptr(x0=global$solution, eval_f=wrapLogLik, lb=lb, ub=ub, opts=local.opts
                , phy=phy, data.list=data.list, nstates=st.count[1], root.type=root.type
                , nsites=nsites)
    print( "Finished." )

    ## Construct the matrices again to return.
    Q.res <- list()
    for( i in 1:nsites) {
        Q.res[[i]] <- matrix(fit$solution[i], nrow=st.count[1], ncol=st.count[1])
        diag(Q.res[[i]]) <- -(colSums(Q.res[[i]]) - fit$solution[i])
    }
    ## Here 'res' is a list with each of the Q matrices. The length of the list is equal to the number of sites in the data.
    out <- list( log.lik=-fit$objective, Q=Q.res, start.par=init.pars, nlopt.opts=fit$options, nlopt.message=fit$message )
    return( out )
}
