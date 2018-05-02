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
    ## Update the function to take an alignment. This will size the Q matrix for each site depending on the
    ##     number of different elements in the site.
    ## ( Need to be extended to accept sites with multiple states - to account for uncertainty. )

    ## Check argument choice.
    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)
    ## Check phylogeny and data:
    if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
    match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
    if( !match.names ) stop("Some species do not match between data and phylogeny!")

    ## Make data checks and get information from the matrix.
    nsites <- ncol(data)
    nstates <- apply(data, 2, function(x) length( unique(x) ) )
    names.data <- rownames(data)
    ## Translate states to numeric, but keep the labels safe.
    states.key <- lapply(1:nsites, function(x)
        cbind( unique(data[,x]), 1:length(unique(data[,x])) )
        )
    for(i in 1:nsites){
        for(j in 1:length(states.key[[i]][,1])){
            id <- data[,i] == states.key[[i]][,1][j]
            data[id,i] <- as.numeric( states.key[[i]][,2] )[j]
        }
    }
    ## Because R is mega dumb:
    data <- apply(as.matrix(data), 2, as.numeric)
    rownames(data) <- names.data
    ## Now we have a data matrix with numbers.

    ## Elements in this data.list can have a different number of columns.
    data.list <- lapply(1:nsites, function(x) make.data.tips.numeric(setNames(as.numeric(data[,x]), rownames(data))) )

    ## Make a function to call the nloptr function.
    wrapLogLik <- function(logx, phy, data.list, nsites, nstates, root.type){
        ## x is a vector of parameters.
        ## length of x is equal to the number of sites.
        ## Search will be in log-space.
        x <- exp(logx) ## Transform back to evaluate the likelihood.
        Q <- list()
        for( i in 1:nsites){
            Q[[i]] <- matrix(x[i], nrow=nstates[i], ncol=nstates[i])
            diag(Q[[i]]) <- -(colSums(Q[[i]]) - x[i])
        }
        lik.site <- sapply(1:nsites, function(x) logLikMk(phy=phy, X=data.list[[x]], Q=Q[[x]], root.type))
        lik <- sum(lik.site)
        return( -lik ) ## Remember that NLOPT is minimizying the function!
    }

    ## Create the vectors for the upper and lower bound for nloptr.
    ## The bound is the same for each of the sites.
    if( bounds[1] < 0 ) stop( "The lower bound cannot be negative." )
    ## Log the bounds in order to search in log space.
    if( bounds[1] == 0 ){
        ## Cannot be log(0)
        bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
    }
    log_lb <- log( rep(bounds[1], times = nsites) )
    log_ub <- log( rep(bounds[2], times = nsites) )
    
    ## Sample the initial parameters for the search.
    ## Here the user can provide a custom start.
    if( is.null(init) ){
        init.pars <- runif(nsites, min=bounds[1], max=bounds[2])
    } else{
        if( any(init < bounds[1]) | any(init > bounds[2]) ) stop("Provided 'init' value is outside search bounds.")
        if( !length(init) == nsites ) stop(" Length of 'init' needs to be the same as the length of the sequence." )
        init.pars <- init
    }
    ## Log the initial value to search in log-space.
    init.pars <- log( init.pars )
        
    ## Get the number of states. At the moment all the sites need to have the same number of states.
    ## No state can be missing in the tips.
    ## Cannot check this anymore. This should work even if this is true.
    ## st.count <- apply(data, 2, function(x) length( unique(x) ) )
    ## if( any( !st.count == st.count[1] ) ) stop( "Some of the sites has a different number of elements." )
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
    global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                   , phy=phy, data.list=data.list, nstates=nstates, root.type=root.type
                   , nsites=nsites)
    print( "Global search solution:" )
    print( paste("Log-lik:", -global$objective, "Pars:", exp(global$solution), sep=" ") )
    print( "Starting local MLE search. (Second pass)" )
    fit <- nloptr(x0=global$solution, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=local.opts
                , phy=phy, data.list=data.list, nstates=nstates, root.type=root.type
                , nsites=nsites)
    print( "Finished." )

    ## Note that the parameters for the model returned will be in log.
    normal.solution <- exp( fit$solution )

    ## Construct the matrices again to return.
    Q.res <- list()
    for( i in 1:nsites) {
        Q.res[[i]] <- matrix(normal.solution[i], nrow=st.count[1], ncol=st.count[1])
        diag(Q.res[[i]]) <- -(colSums(Q.res[[i]]) - normal.solution[i])
    }
    ## Here 'res' is a list with each of the Q matrices. The length of the list is equal to the number of sites in the data.
    out <- list( log.lik=-fit$objective, Q=Q.res, start.par=exp(init.pars), nlopt.opts=fit$options, nlopt.message=fit$message )
    return( out )
}
