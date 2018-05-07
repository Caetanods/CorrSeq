##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Find the maximum likelihood for the model with the gamma function.
##' @param data matrix with species names as rownames.
##' @param phy phylogeny
##' @param model Possible models are "ER" (single global rate) and "DEL" (a global rate for transitions between observed states and another rate for gains and loses of states).
##' @param root.type one of "madfitz", or "equal". Need to extend to accept a observed vector of probabilities.
##' @param ncat categories for the gamma function
##' @param bounds 
##' @param state.space can be "max.div", "site.div" or a numeric (length of 1 or equal to the number of columns in the data). The size of the state space is equal to the largest observed or equal to the per site state diversity plus the number.
##' @param gap.char the character identified as the gap in the data. Default: "-"
##' @param opts the list of options for nloptr. If null it will use the default parameters.
##' @param init 
##' @param verbose 
##' @param n.cores 
##' @param bonds a numeric vector with the lower and upper bonds for the rates
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @export
##' @author daniel
fitMLGammaSpanSpace <- function(data, phy, model = "ER", root.type = "madfitz", ncat = 4, bounds = c(0,100), state.space = "max.div", gap.char = "-", opts = NULL, init = NULL, verbose = TRUE, n.cores = 1){
    ## At the moment the model assumes that all transitions have the same rate.
    ## So just a single rate is estimated for each of the transition matrices. Of course, this rate is changed by the gamma distribution.
    ## This is a fork of the original function.
    ## Here we expand the state space to allow for transitions to states that are not observed among the tips.
    ## Need to think on some rules as how to deal with these. Assuming that the transition rate is the same among the states, we might have different possible options to deal with the sizes of the states.
    
    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)
    model <- match.arg(model, choices=c("ER","DEL"), several.ok=FALSE)
    
    ## Check phylogeny and data:
    ## Consider adding step to organize the species names in the data following the phy names.
    if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
    match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
    if( !match.names ) stop("Secies names do not match between data and phylogeny!")

    ## Re-order the species in the data matrix to match the tree:
    data.order <- match(x=phy$tip.label, table=rownames(data))
    data <- data[data.order,]

    ## Check for invariant sites. Mark these, avoid estimation.
    has.invar <- any( apply(data, 2, function(x) length(unique(x)) == 1) )
    if( has.invar ){
        print("Data contain invariant positions. Returning 'invariant' as the transition matrix for these positions.")
        which.invar <- apply(data, 2, function(x) length(unique(x)) == 1)
        data <- data[,!which.invar]
    }
    
    ## Make data checks and get information from the matrix.
    nsites <- ncol(data)
    names.data <- rownames(data)

    ## nstates is a vector with the sizes of the observed states.
    nstates <- apply(data, 2, function(x) length( unique(x) ) )
    
    ## Translate states to numeric, but keep the labels safe.
    ## If the model is "DEL" then I need to code states such that '-' is always the first state.
    states.key <- list()
    for(i in 1:nsites){
        site.div <- unique(data[,i])
        is.gap <- site.div == gap.char ## The symbol for a gap.
        order.site.div <- c(site.div[is.gap], site.div[!is.gap])
        states.key[[i]] <- cbind(order.site.div, 1:length(order.site.div))
        ## Translate to numeric following the key:
        for(j in 1:length(states.key[[i]][,1])){
            id <- data[,i] == states.key[[i]][,1][j]
            data[id,i] <- as.numeric( states.key[[i]][,2] )[j]
        }
    }

    ## Because R is mega dumb. Need to make sure that the table is really numeric.
    data <- apply(as.matrix(data), 2, as.numeric)
    rownames(data) <- names.data

    ## Expand the state space for each of the sites following the argument choice.
    if( is.character(state.space) ){
        state.space <- match.arg(state.space, choices = c("max.div", "site.div"))
        if(state.space == "max.div"){
            max.div <- max(nstates)
            nstates <- rep(max.div, times=length(nstates))
        }
        if(state.space == "site.div"){
            nstates <- nstates
        }
    } else if( is.numeric( state.space ) ){
        if( !length( state.space ) == 1 & !length( state.space ) == length(nstates) ){
            ## If not a single number. Then need to be a vector with the same length as 'nstates'.
            stop("state.space needs to be 'max.div', 'obs.div', a single number, or a numeric vector with length equal to ncol of data.")
        }
        ## Add the extra state space for all the sites.
        nstates <- nstates + state.space
    } else {
        ## Not any of above, then return error.
        stop("state.space needs to be 'max.div', 'obs.div' or a single number.")
    }

    ## Likelihood function will depend on the model now:
    if( model == "ER" ){        
        make.Q.list <- function(rate, size){
            Q <- matrix(rate, nrow=size, ncol=size)
            diag(Q) <- -(colSums(Q) - rate)
            return(Q)
        }    
        wrapLogLikGammaSimple <- function(obj, phy, data, nstates, k, root.type){
            ## obj is a vector of 2 parameters.
            ## obj[1] = Q[1,2] (shared for all sites), obj[2] = beta.
            ## phy = phylogeny
            ## data = data matrix.
            ## nstates = number of states in each of the sites in the data (now this is a vector).
            x <- exp(obj[1])
            beta <- obj[2]
            ## Will need to be a list.
            Q <- lapply(nstates, function(i) make.Q.list(rate=x[1], size=i) )
            ## Loglik function for the model.
            lik <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=k, it=1, n.cores=n.cores)[[1]]
            return( -lik ) ## Remember that NLOPT is minimizying the function!
        }
    }
    if( model == "DEL" ){
        make.Q.list.DEL <- function(global.rate, del.rate, size){
            ## The '-' state is always the first state in the matrix.
            Q <- matrix(global.rate, nrow=size, ncol=size)
            ## The separate transition is the Q[2+,1] and Q[1,2+]
            Q[1,] <- del.rate
            Q[,1] <- del.rate
            diag(Q) <- sapply(1:size, function(x) -(sum(Q[x,]) - Q[x,x]))
            return(Q)
        }    
        wrapLogLikGammaSimple <- function(obj, phy, data, nstates, k, root.type){
            ## obj is a vector of 3 parameters.
            ## obj[1:2] = has the rate for the observed states obj[1] and for the deletion and insertions obj[2].
            ## obj[3] = the beta rate (for the Gamma function).
            ## phy = phylogeny
            ## data = data matrix.
            ## nstates = number of states in each of the sites in the data (now this is a vector).
            x.obs <- exp(obj[1])
            x.del <- exp(obj[2])
            beta <- obj[3]
            ## Will need to be a list.
            Q <- lapply(nstates, function(i) make.Q.list.DEL(global.rate=x.obs, del.rate=x.del, size=i) )
            ## Loglik function for the model.
            lik <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=k
                                   , it=1, n.cores=n.cores)[[1]]
            return( -lik ) ## Remember that NLOPT is minimizying the function!
        }
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

    ## The "DEL" model has one more parameter.
    ## Need to adjust the length of init and bounds parameters.
    if( model == "DEL" ){
        init.pars <- c(init.pars[1], init.pars[1], init.pars[2])
        log_lb <- c( log( bounds[1] ), log( bounds[1] ), 0 )
        log_ub <- c( log( bounds[2] ), log( bounds[2] ), maxBeta )
    }    

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

    ## Register search time.
    start.time <- Sys.time()
    if( verbose ){
        print( "Starting global MLE search. (First pass)" )
        global <- nloptr(x0=init.pars, eval_f=wrapLogLikGammaSimple, lb=log_lb, ub=log_ub, opts=global.opts
                       , phy=phy, data=data, nstates=nstates, k=ncat, root.type=root.type)
        print( "Global search finished." )
        ## print( paste("Log-lik ", -global$objective
        ##            , "; rate ", exp(global$solution[1])
        ##            , "; alpha ", global$solution[2], collapse="") )
        print( "Starting local MLE search. (Second pass)" )
        fit <- nloptr(x0=global$solution, eval_f=wrapLogLikGammaSimple, lb=log_lb, ub=log_ub, opts=local.opts
                    , phy=phy, data=data, nstates=nstates, k=ncat, root.type=root.type)
        print( "Local search solution:" )
        ## print( paste("Log-lik ", -fit$objective
        ##            , "; rate ", exp(fit$solution[1])
        ##            , "; alpha ", fit$solution[2], collapse="") )
        print( "Reconstructing site-wise Q matrices." )
    } else{
        global <- nloptr(x0=init.pars, eval_f=wrapLogLikGammaSimple, lb=log_lb, ub=log_ub, opts=global.opts
                       , phy=phy, data=data, nstates=nstates, k=ncat, root.type=root.type)
        fit <- nloptr(x0=global$solution, eval_f=wrapLogLikGammaSimple, lb=log_lb, ub=log_ub, opts=local.opts
                    , phy=phy, data=data, nstates=nstates, k=ncat, root.type=root.type)
    }
    ## Register finish search time.
    finish.time <- Sys.time()
    total.time <- format( difftime(finish.time, start.time) )
    
    ## Need to reconstruct the Q matrices for each of the sites.
    ## For this I need to compute the likelihood again.
    ## Apply the original state names to the Q matrices.
    if( model == "ER"){
        solution <- c(exp(fit$solution[1]), fit$solution[2])
        Q <- lapply(nstates, function(i) make.Q.list(rate=solution[1], size=i) )
        beta <- solution[2]
        res <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=ncat, it=1, n.cores=n.cores)[[2]]
    }
    if( model == "DEL"){
        solution <- c(exp(fit$solution[1]), exp(fit$solution[2]), fit$solution[3])
        Q <- lapply(nstates, function(i) make.Q.list.DEL(global.rate=solution[1], del.rate=solution[2], size=i) )
        beta <- solution[3]
        res <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=ncat, it=1, n.cores=n.cores)[[2]]
    }
    ## Paste the state names on the Q matrices.
    for( i in 1:nsites ){
        if(nstates[i] == length(states.key[[i]][,1])){ ## Do not need to append to names.
            rownames(res[[i]]) <- states.key[[i]][,1]
            colnames(res[[i]]) <- states.key[[i]][,1]
        } else{ ## Need to append to names.
            dif.len <- nstates[i] - length(states.key[[i]][,1])
            append.name <- paste("state.", LETTERS[1:dif.len], sep="")
            rownames(res[[i]]) <- c(states.key[[i]][,1], append.name)
            colnames(res[[i]]) <- c(states.key[[i]][,1], append.name)
        }
    }
    ## Before returning the list we need to include a flag for the invariant site positions.
    if( has.invar ){
        res.complete <- list()
        complete.length <- length(which.invar) ## The original size of the data matrix.
        count <- 1 ## A counter for the loop.
        for(i in 1:complete.length){
            if(which.invar[i]){ ## If invariant site.
                res.complete[[i]] <- "invariant"
            }else{
                res.complete[[i]] <- res[[count]]
                count <- count + 1
            }
        }
        res <- res.complete
    }
    
    ## Here 'res' is a list with each of the Q matrices. The length of the list is equal to the number of sites in the data.
    ## Some parameters of the model depend on the choice of model estimate.
    del.rate <- NULL
    start.par <- c(exp(init.pars[1]),init.pars[2])
    alpha <- solution[2]
    if( model == "DEL" ){
        del.rate <- solution[2]
        start.par <- c(exp(init.pars[1]), exp(init.pars[2]), init.pars[3])
        alpha <- solution[3]
    }
    
    out <- list( log.lik=-fit$objective, Q=res, global.rate=solution[1], del.rate=del.rate, alpha=solution[3]
              , start.par=start.par, nlopt.global.search=global.opts, nlopt.local.search=local.opts
              , nlopt.message=fit$message, search.time=total.time, model=model)
    return( out )
}

