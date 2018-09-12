##' Fits the state diversity model.
##'
##' This model is different from the sequence trait model in the sense that it does not focus on rates of transitions between states on each of the positions of a sequence trait. This models make inferences about the rate in which the diversity of a given sequence position increase or decrease compared to the other rates across the sequence trait.
##' @title Find the maximum likelihood for the state diversity model.
##' @param data matrix with species names as rownames.
##' @param phy phylogeny
##' @param rate.model options are "correlated", "gamma", and "single.rate"
##' @param symmetric if the model have equal rates for birth and death events.
##' @param root.type one of "madfitz", or "equal". Need to extend to accept a observed vector of probabilities.
##' @param add.invar TRUE or FALSE. Whether to compute transition rates for the invariant sites or not.
##' @param ncat categories for the gamma function
##' @param bounds a 2x2 numeric matrix. First column are the limits for the birth rate and second column are limits for the death rate. First line is the lower limit and second line is the upper limit.
##' @param div.max can be "observed" for the limits given the observed diversity on the position or a numeric value with the number of states to be added above and below the observed max and min diversity in each position.
##' @param opts the list of options for nloptr. If null it will use default parameters.
##' @param init parameters to initialize the search
##' @param verbose whether to print messages to the screen.
##' @param n.cores number of cores to split the computation of the likelihood
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @importFrom expm expm
##' @importFrom ape reorder.phylo
##' @export
##' @author daniel
fitMLStateDiversity <- function(data, phy, rate.model = "gamma", symmetric = TRUE, root.type = "madfitz", add.invar = FALSE, ncat = 4, bounds = NULL, div.max = "observed", opts = NULL, init = NULL, verbose = TRUE, n.cores = 1){
    
    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)
    rate.model <- match.arg(rate.model, choices=c("correlated", "gamma", "single.rate"), several.ok=FALSE)
    if( !is.logical(add.invar) ) stop("Argument 'add.invar' need to TRUE or FALSE.")

    ## Check phylogeny and data:
    if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
    match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
    if( !match.names ) stop("Secies names do not match between data and phylogeny!")
    ## Two positions in the trait is not enough for this model:
    if( ncol(data) < 3 ) warning("Less than 3 positions in the sequence. Maybe not enough information.")
    if( ncol(data) == 1 ) stop("Cannot work with a sequence of a single position!")

    ## Check if the 'bounds' argument has the correct format:
    if( !is.null( bounds ) ){
        if( !is.matrix(bounds) ) stop( "Argument 'bounds' needs to be a matrix." )
        if( !ncol( bounds ) == 2 ) stop( "Argument 'bounds' is a 2x2 matrix." )
        if( !nrow( bounds ) == 2 ) stop( "Argument 'bounds' is a 2x2 matrix." )        
        if( any(bounds[1,] < 0) ) stop( "The lower bound cannot be negative." )
    } else{
        ## Set the bounds of the search to defaults.
        bounds <- cbind( c(0,10), c(0,10) )
    }
    
    ## Re-order the species in the data matrix to match the tree:
    data.order <- match(x=phy$tip.label, table=rownames(data))
    data <- data[data.order,]

    if( !add.invar ){
        ## Check for invariant sites. Mark these, avoid estimation.
        has.invar <- any( apply(data, 2, function(x) length(unique(x)) == 1) )
        if( has.invar ){
            print("Data contain invariant positions. Returning 'invariant' as the transition matrix for these positions.")
            which.invar <- apply(data, 2, function(x) length(unique(x)) == 1)
            data <- data[,!which.invar]
        }
    } else{
        ## Create 'has.invar' as FALSE.
        has.invar <- FALSE
    }
    
    ## Make data checks and get information from the matrix.
    nsites <- ncol(data)
    names.data <- rownames(data)

    ## Compute the max diversity observed on the sites.
    if( div.max == "observed" ){
        add.div <- 0
    } else{
        if( is.numeric( div.max ) ){
            ## The added diversity need to be integer.
            add.div <- floor(div.max)
        } else{
            stop( "Argument 'div.max' need to be a character or a numeric value." )
        }
    }
    
    ## Produce the data matrices for the likelihood of the model.
    ## The format is the same of the sequence state model, but follows a different data type.
    if( !all( apply(data, 2, is.numeric) ) ){
        stop("data need to be numeric integer values.")
    }
    
    ## Make the list of data matrices:
    Xlist <- list()
    for( i in 1:nsites ){
        sites.range <- ceiling( range( data[,i] ) )
        div.seq <- sites.range[1]:sites.range[2]
        site.matrix <- matrix(0, nrow = nrow(data), ncol = length(div.seq))
        for( j in 1:length(div.seq) ){
            ## Populate the matrix:
            site.matrix[,j] <- as.numeric( data[,i] == div.seq[j] )
        }
        ## Drop unobserved gaps. Leave one state between.
        obs.gaps <- as.logical( colSums( site.matrix ) )
        drop <- vector()
        for(j in 1:(length(obs.gaps)-1)){
            if(!obs.gaps[j] && !obs.gaps[j+1]){
                drop <- c(drop, j+1)
            }
        }
        if( length(drop) > 0 ){
            ## Drop has length 0 if nothing to drop.
            site.matrix <- site.matrix[,-drop]
        }        
        ## Add the padded state to the left and to the right.
        if( add.div > 0 ){
            pad.matrix <- matrix(0, ncol = add.div, nrow = nrow(data))
            site.matrix <- cbind(pad.matrix, site.matrix, pad.matrix)
        }
        Xlist[[i]] <- site.matrix
    }
    
    ## Compute some quantities related to the phylogeny that will be used by the likelihood function:
    prun.phy <- reorder.phylo(x = phy, order = "postorder")
    edge_mat <- prun.phy$edge
    ## The root type need to be a numeric.
    root_type <- as.numeric( switch(root.type, "madfitz" = 1, "equal" = 0) )
    n_nodes <- prun.phy$Nnode
    n_tips <- length( prun.phy$tip.label )
    root_node <- n_tips + 1
    n_states <- sapply(Xlist, function(i) ncol(i) )
    edge_len <- prun.phy$edge.length
    parents <- unique( prun.phy$edge[,1] )

    ## The type of model will depend if symmetric or assymetric.
    if( symmetric ){        
        make.Q.list <- function(rate, size){
            Q <- matrix(0, nrow=size, ncol=size)
            for( i in 2:size ){
                Q[i-1,i] <- rate
                Q[i,i-1] <- rate
            }            
            diag(Q) <- -colSums(Q)
            return(Q)
        }
        
        ## Check if the model is the autocorrelated:
        if( rate.model == "correlated" ){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type){
                ## obj is a vector of variable length.
                ## obj[1] = Q[1,2] (shared for all sites), obj[2] = beta, obj[3:n] = the elements of the M matrix (for the autocorrelation).
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                ## gap.key is ignored in this case.
                x <- exp(obj[1])
                beta <- obj[2]
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list(rate=x[1], size=i) )
                ## The M matrix for the autocorrelation:
                ## This is a probability matrix and is different from a Q matrix:
                M <- makeSymDTMCMatrix(size = k, obj[3:length(obj)])
                if( sum(M) < 0.001 ){ ## Protect from bad behavior in sum function.
                    return( Inf ) ## NLOPT is minimizying the function!
                }
                ## Loglik function for the model.
                lik <- logLikAutoDiscGamma_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                            , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                            , root_node=root_node, X = data, Q = Q, M = M
                                            , root_type=root_type, beta=beta, k=k, n.cores=n.cores)
                ## lik <- logLikAutoDiscGamma(phy=phy, X=data, Q=Q, M=M, root.type=root.type
                ##                          , beta=beta, k=k, n.cores=n.cores)
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
        if(rate.model == "gamma"){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type){
                ## obj is a vector of 2 parameters.
                ## obj[1] = Q[1,2] (shared for all sites), obj[2] = beta.
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                ## gap.key is ignored in this case.
                x <- exp(obj[1])
                beta <- obj[2]
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list(rate=x[1], size=i) )
                ## Loglik function for the model.
                lik <- loglikGammaSimple_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                            , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                            , root_node=root_node, X = data, Q = Q
                                            , root_type=root_type, beta=beta, k=k, n.cores=n.cores)[[1]]
                ## lik <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=k, n.cores=n.cores)[[1]]
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
        if( rate.model == "single.rate" ){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type){
                ## obj is a vector of 1 parameter.
                ## obj[1] = Q[1,2] (shared for all sites).
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                ## k is ignored in this case.
                ## gap.key is ignored in this case.
                x <- exp(obj[1])
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list(rate=x[1], size=i) )
                ## Loglik function for the model.
                lik <- loglikSingleRate_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                        , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                        , root_node=root_node, X = data, Q = Q
                                        , root_type=root_type, n.cores=n.cores)
                ## lik <- loglikSingleRate(phy=phy, X=data, Q=Q, root.type=root.type, n.cores=n.cores)
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
    } else{
        ## Separated rates for birth and death of elements.
        make.Q.list.asym <- function(birth.rate, death.rate, size){
            Q <- matrix(0, nrow=size, ncol=size)
            for( i in 2:size ){
                Q[i-1,i] <- birth.rate
                Q[i,i-1] <- death.rate
            }            
            diag(Q) <- - colSums(Q)
            return(Q)
        }
        ## Check if the model is the autocorrelated:
        if( rate.model == "correlated" ){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type){
                ## obj is a vector of 3 parameters.
                ## obj[1:2] = has the rate for the observed states obj[1] and for the deletion and insertions obj[2]
                ## obj[3] = the beta rate (for the Gamma function).
                ## obj[4:n] = the elements of the M matrix (for the autocorrelation).
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                x.birth <- exp(obj[1])
                x.death <- exp(obj[2])
                beta <- obj[3]
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list.asym(birth.rate=x.birth, death.rate=x.death
                                                                 , size=i)
                            )
                ## The M matrix for the autocorrelation:
                ## This is a probability matrix and is different from a Q matrix:
                M <- makeSymDTMCMatrix(size = k, obj[3:length(obj)])
                if( sum(M) < 0.001 ){ ## Protect from bad behavior in sum function.
                    return( Inf )
                }
                
                ## Loglik function for the model.
                lik <- logLikAutoDiscGamma_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                           , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                           , root_node=root_node, X = data, Q = Q, M = M
                                           , root_type=root_type, beta=beta, k=k, n.cores=n.cores)
                ## lik <- logLikAutoDiscGamma(phy=phy, X=data, Q=Q, M=M, root.type=root.type, beta=beta, k=k, n.cores=n.cores)
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
        if(rate.model == "gamma"){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type){
                ## obj is a vector of 3 parameters.
                ## obj[1:2] = has the rate for the observed states obj[1] and for the deletion and insertions obj[2].
                ## obj[3] = the beta rate (for the Gamma function).
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                x.birth <- exp(obj[1])
                x.death <- exp(obj[2])
                beta <- obj[3]
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list.asym(birth.rate=x.birth, death.rate=x.death
                                                                 , size=i)
                            )
                ## Loglik function for the model.
                lik <- loglikGammaSimple_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                         , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                         , root_node=root_node, X = data, Q = Q
                                         , root_type=root_type, beta=beta, k=k, n.cores=n.cores)[[1]]
                ## lik <- loglikGammaSimple(phy=phy, X=data, Q=Q, root.type=root.type, beta=beta, k=k
                ##                        , n.cores=n.cores)[[1]]
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
        if(rate.model == "single.rate"){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type){
                ## obj is a vector of 2 parameters.
                ## obj[1:2] = has the rate for the observed states obj[1] and for the deletion and insertions obj[2].
                x.birth <- exp(obj[1])
                x.death <- exp(obj[2])
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list.asym(birth.rate=x.birth, death.rate=x.death
                                                                 , size=i)
                            )
                ## Loglik function for the model.
                lik <- loglikSingleRate_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                        , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                        , root_node=root_node, X = data, Q = Q
                                        , root_type=root_type, n.cores=n.cores)
                ## lik <- loglikSingleRate(phy=phy, X=data, Q=Q, root.type=root.type, n.cores=n.cores)
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
    }
    
    ## Search for upper and lower bounds for the beta parameter for the Gamma rates that do not produce 0 rate values.
    if( rate.model %in% c("correlated","gamma") ){
        ## The single rate model does not have a beta parameter.
        beta.bounds <- findMaxBeta(ncat)
    }

    ## Bounce back the bounds if they are too low:
    if( any(bounds[1,] == 0) ){
        ## Cannot be log(0)
        if( bounds[1,1] < .Machine$double.eps ){
            bounds[1,1] <- .Machine$double.eps
        }
        if( bounds[1,2] < .Machine$double.eps ){
            bounds[1,2] <- .Machine$double.eps
        }            
    }
    
    ## Create the vectors for the upper and lower bound for nloptr.
    ## The bound is the same for each of the sites.
    if( rate.model == "correlated" ){
        ## Set bounds of 0 to a very small number.
        
        ## Finde the size of the M matrix
        M.vec.size <- sum( upper.tri( diag(ncat), diag = TRUE ) )

        if( symmetric ){
            ## Make nlopt bounds vectors. The M matrix has always bound between 0 and 1.
            log_lb <- c( log( bounds[1,1] ), beta.bounds[1], rep(0, times=M.vec.size) )
            log_ub <- c( log( bounds[2,1] ), beta.bounds[2], rep(1, times=M.vec.size) )
        } else{
            log_lb <- c( log( bounds[1,1] ), log( bounds[1,2] ), beta.bounds[1], rep(0, times=M.vec.size) )
            log_ub <- c( log( bounds[2,1] ), log( bounds[2,2] ), beta.bounds[2], rep(1, times=M.vec.size) )
        }
    }
    if( rate.model == "gamma"){
        ## Log the bounds in order to search in log space.
        if( symmetric ){
            log_lb <- c( log( bounds[1,1] ), beta.bounds[1] )
            log_ub <- c( log( bounds[2,1] ), beta.bounds[2] )
        } else{
            log_lb <- c( log( bounds[1,1] ), log( bounds[1,2] ), beta.bounds[1] )
            log_ub <- c( log( bounds[2,1] ), log( bounds[2,2] ), beta.bounds[2] )
        }
    }
    if( rate.model == "single.rate" ){
        if( symmetric ){
            log_lb <- log( bounds[1,1] )
            log_ub <- log( bounds[2,1] )
        } else{
            log_lb <- c( log( bounds[1,1] ), log( bounds[1,2] ) )
            log_ub <- c( log( bounds[2,1] ), log( bounds[2,2] ) )
        }
    }

    ## Sample the initial parameters for the search.
    ## Here the user can provide a custom start.
    if( is.null(init) ){
        if( rate.model == "correlated" ){
            ## Reject sample for a good starting value:
            M <- 0
            while( sum(M) < 0.001 ){
                init.pars.M <- runif(M.vec.size, min=0, max=1)
                M <- makeSymDTMCMatrix(size = ncat, pars = init.pars.M)
            }
            if( symmetric ){
                init.pars <- c( log(runif(1, min=bounds[1,1], max=bounds[2,1]))
                             , runif(1, min=beta.bounds[1], max=beta.bounds[2])
                             , init.pars.M )
            }else{
                init.pars <- c( log(runif(1, min=bounds[1,1], max=bounds[2,1]))
                             , log(runif(1, min=bounds[1,2], max=bounds[2,2]))
                             , runif(1, min=beta.bounds[1], max=beta.bounds[2])
                             , init.pars.M )
            }
        }
        if(rate.model == "gamma"){
            if( symmetric ){
                init.pars <- c(  log(runif(1, min=bounds[1,1], max=bounds[2,1]))
                             , runif(1, min=beta.bounds[1], max=beta.bounds[2]) )
            }else{
                init.pars <- c(  log( runif(1, min=bounds[1,1], max=bounds[2,1]) )
                               , log( runif(1, min=bounds[1,2], max=bounds[2,2]) )
                             , runif(1, min=beta.bounds[1], max=beta.bounds[2]) )
            }
        }
        if(rate.model == "single.rate"){
            if( symmetric ){
                init.pars <- log(runif(1, min=bounds[1,1], max=bounds[2,1]))
            }else{
                init.pars <- c(  log( runif(1, min=bounds[1,1], max=bounds[2,1]) )
                             , log( runif(1, min=bounds[1,2], max=bounds[2,2]) ) )
            }
        }
    } else{
        stop("Custom init not implemented yet.")
    }

    ## Create the list of options for local search of nloptr:
    if( is.null(opts) ){
        ## nlopt.opts <- list(algorithm="NLOPT_LN_SBPLX", "ftol_rel"=1e-08, "maxtime"=170000000, "maxeval"=10000)
        ## Increasing the tolerance of the Global search here because it is going to be followed by a local search.
        global.opts <- list("algorithm"="NLOPT_GN_DIRECT", "maxeval"=10000, "ftol_rel"=0.0001)
        local.opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"=1000000, "ftol_rel"=.Machine$double.eps^0.5)
    } else{
        if( !is.list(opts) ) stop( "The argument 'opts' needs to be a list format" )
        local.opts <- opts
        global.opts <- opts
        ## Set the algorithm for the global search.
        global.opts$algorithm <- "NLOPT_GN_DIRECT"
    }

    ## Register search time.
    start.time <- Sys.time()
    if( verbose ){
        print( "Starting global MLE search. (First pass)" )
        global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                       , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                       , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                       , k=ncat, root_type=root_type)
        print( "Global search finished." )
        ## print( paste("Log-lik ", -global$objective
        ##            , "; rate ", exp(global$solution[1])
        ##            , "; alpha ", global$solution[2], collapse="") )
        print( "Starting local MLE search. (Second pass)" )
        fit <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                       , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                       , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                       , k=ncat, root_type=root_type)
        print( "Local search solution:" )
        ## print( paste("Log-lik ", -fit$objective
        ##            , "; rate ", exp(fit$solution[1])
        ##            , "; alpha ", fit$solution[2], collapse="") )
        print( "Reconstructing site-wise Q matrices." )
    } else{
        global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                       , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                       , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                       , k=ncat, root_type=root_type)
        fit <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                       , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                       , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                       , k=ncat, root_type=root_type)
    }
    ## Register finish search time.
    finish.time <- Sys.time()
    total.time <- format( difftime(finish.time, start.time) )

    ## Need to reconstruct the Q matrices for each of the sites.
    ## For this I need to compute the likelihood again.
    ## Apply the original state names to the Q matrices.
    if( symmetric ){
        solution <- c(exp(fit$solution[1]), fit$solution[2])
        Q <- lapply(n_states, function(i) make.Q.list(rate=solution[1], size=i) )
        beta <- solution[2]
        if( rate.model == "correlated" ){
            ## In this case need to append to the solution object.
            solution <- c(solution, fit$solution[3:length(fit$solution)])
            M <- makeSymDTMCMatrix(size = ncat, pars = solution[3:length(solution)] )
            ## Need to make sure that the format for the output is the same between models.
            res <- getSiteRatesAutoDiscGamma(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                           , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                           , root_node=root_node, root_type=root_type, k=ncat
                                           , X=Xlist, Q=Q, M=M, beta=beta, n.cores=n.cores)
        }
        if(rate.model == "gamma"){
            res <- loglikGammaSimple_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                     , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                     , root_node=root_node, X = Xlist, Q = Q
                                     , root_type=root_type, beta=beta, k=ncat, n.cores=n.cores)[[2]]
        }
        if(rate.model == "single.rate"){
            res <- Q
        }
    } else{
        ## The model is assymmetric:
        solution <- c(exp(fit$solution[1]), exp(fit$solution[2]), fit$solution[3])
        Q <- lapply(n_states, function(i) make.Q.list.asym(birth.rate=solution[1], death.rate=solution[2]
                                                       , size=i) )
        beta <- solution[3]
        if( rate.model == "correlated" ){
            ## Need to append the rest of the parameters to the solution object.
            solution <- c(solution, fit$solution[3:length(fit$solution)])
            M <- makeSymDTMCMatrix(size = ncat, pars = solution[3:length(solution)] )
            ## Need to make sure that the format for the output is the same between models.
            res <- getSiteRatesAutoDiscGamma(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                           , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                           , root_node=root_node, root_type=root_type, k=ncat
                                           , X=Xlist, Q=Q, M=M, beta=beta, n.cores=n.cores)
        }
        if(rate.model == "gamma"){
            res <- loglikGammaSimple_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                     , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                     , root_node=root_node, X = Xlist, Q = Q
                                     , root_type=root_type, beta=beta, k=ncat, n.cores=n.cores)[[2]]
        }
        if(rate.model == "single.rate"){
            res <- Q
        }            
    }

    ## Add the state names for the Q matrices.
    for( i in 1:nsites ){
        rownames(res[[i]]) <- colnames(Xlist[[i]])
        colnames(res[[i]]) <- colnames(Xlist[[i]])
    }
    
    ## Before returning the list we need to include a flag for the invariant site positions.
    if( !add.invar & has.invar ){
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

    ## Some parameters of the model depend on the choice of model estimate.
    ## Placeholders for the parameters:
    death.rate <- NULL
    auto.matrix <- NULL
    alpha <- NULL
    if( rate.model == "correlated" ){
        auto.matrix <- M
    }
    if( rate.model == "single.rate" ){
        start.par <- exp(init.pars[1])
        if( !symmetric ){
            ## assymetric model
            death.rate <- solution[2]
            start.par <- c(exp(init.pars[1]), exp(init.pars[2]))
        }
    } else{        
        start.par <- c( exp(init.pars[1]), init.pars[2] )
        alpha <- solution[2]
        if( !symmetric ){
            ## assymetric model            
            death.rate <- solution[2]
            start.par <- c(exp(init.pars[1]), exp(init.pars[2]), init.pars[3])
            alpha <- solution[3]
        }
    }
    
    out <- list( log.lik=-fit$objective, Q=res, birth.rate=solution[1], auto.matrix=auto.matrix
              , death.rate=death.rate, alpha=alpha, start.par=start.par, nlopt.global.search=global.opts
              , nlopt.local.search=local.opts, nlopt.message=fit$message, search.time=total.time
              , symmetric=symmetric, rate.model=rate.model)
    return( out )
}
