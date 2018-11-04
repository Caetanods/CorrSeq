##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Find the maximum likelihood for the model with the gamma function.
##' @param data matrix with species names as rownames.
##' @param phy phylogeny
##' @param Q.model Possible models are "ER" (single global rate) and "DEL" (a global rate for transitions between observed states and another rate for gains and loses of states).
##' @param rate.model options are "correlated", "gamma", and "single.rate"
##' @param root.type one of "madfitz", or "equal". Need to extend to accept a observed vector of probabilities.
##' @param add.invar TRUE or FALSE. Whether to compute transition rates for the invariant sites or not.
##' @param poly.key a named list. Each element is a vector with the states represented by the polymorphism symbol. Names of the list need to match the polymorphism symbols in the data.
##' @param ncat categories for the gamma function
##' @param init.M whether the initial state for the M matrix (for the correlated model) should have starting state equal to the gamma model (i.e., equal probabilities for all the transitions).
##' @param bounds a numeric vector of length 2 with the lower and upper bonds for the rates
##' @param state.space can be "max.div", "site.div" or a numeric (length of 1 or equal to the number of columns in the data). The size of the state space is equal to the largest observed or equal to the per site state diversity plus the number.
##' @param gap.char the character identified as the gap in the data. Default: "-"
##' @param opts the list of options for nloptr. If null it will use the default parameters.
##' @param search.global whether to perform a global MLE search before the local MLE search. Default is TRUE.
##' @param init 
##' @param verbose 
##' @param n.cores 
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @importFrom expm expm
##' @importFrom ape reorder.phylo
##' @export
##' @author daniel
fitMLGammaSpanSpace <- function(data, phy, Q.model = "ER", rate.model = "gamma", root.type = "madfitz", add.invar = FALSE, poly.key = NULL, ncat = 4, init.M = FALSE, bounds = NULL, state.space = "max.div", gap.char = "-", opts = NULL, search.global = TRUE, init = NULL, verbose = TRUE, n.cores = 1){
    ## At the moment the model assumes that all transitions have the same rate.
    ## So just a single rate is estimated for each of the transition matrices. Of course, this rate is changed by the gamma distribution.
    ## This is a fork of the original function.
    ## Here we expand the state space to allow for transitions to states that are not observed among the tips.
    ## Need to think on some rules as how to deal with these. Assuming that the transition rate is the same among the states, we might have different possible options to deal with the sizes of the states.
    
    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)
    rate.model <- match.arg(rate.model, choices=c("correlated", "gamma", "single.rate"), several.ok=FALSE)
    Q.model <- match.arg(Q.model, choices=c("ER","DEL"), several.ok=FALSE)
    if( !is.logical(add.invar) ) stop("Argument 'add.invar' need to TRUE or FALSE.")

    ## Check the poly.key argument.
    if( !is.null(poly.key) ){
        if( is.list(poly.key) ){
            poly.states <- names(poly.key)
            if( !all( poly.states %in% c(data) ) ) stop( "Name(s) of 'poly.key' not present in data." )
        } else{
            stop( "Argument 'poly.key' needs to be NULL or a named list. See 'Details'." )
        }
    }
    
    ## Check phylogeny and data:
    ## Consider adding step to organize the species names in the data following the phy names.
    if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
    match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
    if( !match.names ) stop("Species names do not match between data and phylogeny!")
    ## Two positions in the trait is not enough for this model:
    if( ncol(data) < 3 ) warning("Less than 3 positions in the sequence. Maybe not enough information.")
    if( ncol(data) == 1 ) stop("Cannot work with a sequence of a single position!")

    ## If the Q.model = "DEL" then the data need to show at least 1 instance of the gap.char:
    if( Q.model == "DEL" & !any(c(data) == gap.char) ) stop( "Chosen model assumes 'gaps' but symbol in 'gap.char' not found in the data. Change 'gap.char' or chose other model option." )

    ## Check if the 'bounds' argument has the correct format:
    if( !is.null( bounds ) ){
        if( !length( bounds ) == 2) stop( "Wrong format for the 'bounds' argument." )
        if( bounds[1] < 0 ) stop( "The lower bound cannot be negative." )
    } else{
        ## Set the bounds of the search to defaults.
        bounds <- c(0,100)
    }
    
    ## Re-order the species in the data matrix to match the tree:
    data.order <- match(x=phy$tip.label, table=rownames(data))
    data <- data[data.order,]

    if( !add.invar ){
        ## Check for invariant sites. Mark these, avoid estimation.
        ## This check needs to be expanded for the case of polymorphisms.
        ## The idea is that if a polymorphism would cause a site to be invariant, then there is strong evidence to do so. We will only "collapse" the site to invariant if the number of polymorphisms is low, lower than 10% of the sites.
        if( !is.null(poly.key) ){
            ## Check if any is invariant:
            which.invar <- apply(data, 2, function(x) get.invar(x = x, poly.key = poly.key) )
        } else{
            which.invar <- apply(data, 2, function(x) length(unique(x)) == 1)
        }
        
        if( any(which.invar) ){
            print("Data contain invariant positions. Returning 'invariant' as the transition matrix for these positions.")
            data <- data[,!which.invar]
            has.invar <- TRUE ## Flag used downstream.
        } else{
            has.invar <- FALSE ## Flag used downstream.
        }        
    } else{
        ## Create 'has.invar' as FALSE.
        has.invar <- FALSE
    }
    
    ## Make data checks and get information from the matrix.
    nsites <- ncol(data)
    names.data <- rownames(data)

    ## If the Q.model is "DEL" then I need to code states such that '-' is always the first state.
    nstates <- rep(0, times=nsites) ## Number of observed states in each site.
    gap.key <- rep(FALSE, times=nsites) ## The key for sites with gaps
    Xlist <- list() ## The list of matrices to use in the likelihood.

    ## When making the matrix of the data below we need to take care with cases of polymorphic states.
    ## Sometimes the polymorphic state will be "alone" in the position.
    for(i in 1:nsites){
        site.obs <- unique(data[,i])
        check.poly <- FALSE ## This is FALSE by default.
        if( !is.null(poly.key) ){
            site.div <- c()
            site.poly <- c()
            for( z in 1:length( site.obs ) ){
                is.poly <- site.obs[z] %in% names(poly.key) ## The symbol(s) for polymorphism.
                if( is.poly ){
                    ## Check if can be translated.
                    trans.poly <- poly.key[[ names(poly.key) == site.obs[z] ]]
                    keep.poly <- any( site.obs %in% trans.poly )
                    if( keep.poly ){
                        ## Add to the polymorphic vector of states.
                        site.poly <- c(site.poly, site.obs[z])
                    } else{
                        ## Add to the vector of fixed states.
                        ## This is because we cannot translate it.
                        site.div <- c(site.div, site.obs[z])
                    }
                } else{
                    ## Not polymorphic.
                    site.div <- c(site.div, site.obs[z])
                }
            }
            ## If after checking for translation the vector of polymorphisms is empty, then set flag to ignore polymorphic states for this site.
            if( length( site.poly ) > 0 ){
                check.poly <- TRUE
            }
        }
        
        is.gap <- site.div == gap.char ## The symbol for a gap.
        gap.key[i] <- any(is.gap)
        order.site.div <- c(site.div[is.gap], site.div[!is.gap])
        nstates[i] <- length( order.site.div )
        ## Don't need to translate states to numbers. Can work with the symbols.
        site.mat <- matrix(0, ncol = length(order.site.div), nrow = nrow(data))
        if( !is.null(poly.key) & check.poly ){ ## Need to check for polymorphic states.
            for(j in 1:length(data[,i])){
                if( data[j,i] %in% names(poly.key) ){
                    poly.set <- poly.key[[ names(poly.key) == data[j,i] ]]
                    site.mat[j,] <- as.numeric(order.site.div %in% poly.set)
                    ## Set 1 for all states present in this site.
                } else{
                    site.mat[j,] <- as.numeric( order.site.div == data[j,i] )
                }
            }
        } else{ ## Don't have polymorphic states.
            for(j in 1:length(data[,i])){
                site.mat[j,] <- as.numeric( order.site.div == data[j,i] )
            }
        }
        colnames( site.mat ) <- order.site.div
        rownames( site.mat ) <- rownames( data )
        Xlist[[i]] <- site.mat
    }

    ## Compute the expanded state space for each of the sites given function arguments.
    if( is.character(state.space) ){
        state.space <- match.arg(state.space, choices = c("max.div", "site.div"))
        if(state.space == "max.div"){
            max.div <- max(nstates)
            nstates <- rep(max.div, times=length(nstates))
        }
        if(state.space == "site.div"){
            nstates <- nstates
            if( add.invar ) stop("Invariant sites cannot be included if state.space is 'site.div'.")
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

    ## Fix the dimension of the data matrices given the state space.
    for( i in 1:nsites ){
        obs.state.space <- ncol( Xlist[[i]] )
        add.state.space <- nstates[i] - obs.state.space ## Can be 0 or positive.
        if( add.state.space == 0 ) next
        mat.add.state <- matrix(0, nrow = nrow(Xlist[[i]]), ncol=add.state.space)
        rownames(mat.add.state) <- rownames( Xlist[[i]] )
        tmp.Xlist <- cbind(Xlist[[i]], mat.add.state) ## Temporary copy of the matrix.
        append.col.name <- paste("extra.", LETTERS[1:add.state.space], sep="")
        colnames( tmp.Xlist ) <- c( colnames( Xlist[[i]] ), append.col.name )
        Xlist[[i]] <- tmp.Xlist
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

    ## Likelihood function will depend on the Q.model now:
    if( Q.model == "ER" ){        
        make.Q.list <- function(rate, size){
            Q <- matrix(rate, nrow=size, ncol=size)
            diag(Q) <- -(colSums(Q) - rate)
            return(Q)
        }
        
        ## Check if the model is the autocorrelated:
        if( rate.model == "correlated" ){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type, gap.key){
                ## obj is a vector of variable length.
                ## obj[1] = Q[1,2] (shared for all sites), obj[2] = beta, obj[3] = correlation for the bivariate Gamma.
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                ## gap.key is ignored in this case.
                x <- exp(obj[1])
                beta <- obj[2]
                rho <- obj[3] ## This is the correlation of the bivariate Gamma.
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list(rate=x, size=i) )
                
                ## The M matrix for the autocorrelation:
                cat <- qgamma((1:(k - 1))/k, shape = beta, rate = beta)
                cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
                ## alpha and beta are the same thing.
                M <- computeM(rate_cat = cat, alpha = beta, rho = rho, k = k)
                gamma.rates <- discreteGamma(shape = beta, ncats = k) ## These are the k rates.
                if( !is.matrix(M) ){ ## A single rate category!
                    ## The likelihood of the model considering only the last category.
                    ## Because alpha is such that categories 1 to k-1 are empty (width of 0).
                    scaledQ <- lapply(Q, function(x) x * gamma.rates[k])
                    lik <- loglikSingleRate_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                            , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                            , root_node=root_node, X = data, Q = scaledQ
                                            , root_type=root_type, n.cores=n.cores)
                } else{
                    ## Check if M is a doubly stochastic matrix. Otherwise, reject.
                    if( any(round(rowSums(M), digits=10) != 1.0) | any(round(colSums(M), digits=10) != 1.0) ){
                        ## Bad M matrix, this will happen sometimes. Return bad lik.
                        return( Inf )
                    }
                    ## Loglik function for the model.
                    lik <- logLikAutoDiscGamma_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                               , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                               , root_node=root_node, X = data, Q = Q, M = M
                                               , root_type=root_type, beta=beta, k=k, n.cores=n.cores)
                }
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
        if(rate.model == "gamma"){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type, gap.key){
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
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type, gap.key){
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
    }
    if( Q.model == "DEL" ){
        make.Q.list.DEL <- function(global.rate, del.rate, size, is.gap){
            ## The '-' state is always the first state in the matrix.
            Q <- matrix(global.rate, nrow=size, ncol=size)
            ## The separate transition is the Q[2+,1] and Q[1,2+]
            if( is.gap ){
                Q[1,] <- del.rate
                Q[,1] <- del.rate
            }
            diag(Q) <- sapply(1:size, function(x) -(sum(Q[x,]) - Q[x,x]))
            return(Q)
        }
        ## Check if the model is the autocorrelated:
        if( rate.model == "correlated" ){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type, gap.key){
                ## obj is a vector of 3 parameters.
                ## obj[1] = global rate, obj[2] = rate for deletions, obj[3] = the shape for the Gamma, obj[4] the correlation for the bivariate Gamma.
                ## phy = phylogeny
                ## data = data matrix.
                ## nstates = number of states in each of the sites in the data (now this is a vector).
                x.obs <- exp(obj[1])
                x.del <- exp(obj[2])
                beta <- obj[3]
                rho <- obj[4]
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list.DEL(global.rate=x.obs, del.rate=x.del, size=i
                                                                , is.gap=gap.key[i])
                            )
                ## The M matrix for the autocorrelation:
                cat <- qgamma((1:(k - 1))/k, shape = beta, rate = beta)
                cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
                M <- computeM(rate_cat = cat, alpha = beta, rho = rho, k = k)
                gamma.rates <- discreteGamma(shape = beta, ncats = k)
                if( !is.matrix(M) ){ ## A single rate category!
                    ## The likelihood of the model considering only the last category.
                    scaledQ <- lapply(Q, function(x) x * gamma.rates[k])
                    lik <- loglikSingleRate_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                            , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                            , root_node=root_node, X = data, Q = scaledQ
                                            , root_type=root_type, n.cores=n.cores)
                } else{
                    ## Check if M is a doubly stochastic matrix. Otherwise, reject.
                    if( any(round(rowSums(M), digits=10) != 1.0) | any(round(colSums(M), digits=10) != 1.0) ){
                        ## Bad M matrix, this will happen sometimes. Return bad lik.
                        return( Inf )
                    }
                    ## Loglik function for the model.
                    lik <- logLikAutoDiscGamma_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                               , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                               , root_node=root_node, X = data, Q = Q, M = M
                                               , root_type=root_type, beta=beta, k=k, n.cores=n.cores)
                }
                return( -lik ) ## Remember that NLOPT is minimizying the function!
            }
        }
        if(rate.model == "gamma"){
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type, gap.key){
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
                Q <- lapply(n_states, function(i) make.Q.list.DEL(global.rate=x.obs, del.rate=x.del, size=i
                                                               , is.gap=gap.key[i])
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
            wrapLogLik <- function(obj, n_nodes, n_tips, n_states, edge_len, edge_mat, parents, root_node, data, k, root_type, gap.key){
                ## obj is a vector of 2 parameters.
                ## obj[1:2] = has the rate for the observed states obj[1] and for the deletion and insertions obj[2].
                x.obs <- exp(obj[1])
                x.del <- exp(obj[2])
                ## Will need to be a list.
                Q <- lapply(n_states, function(i) make.Q.list.DEL(global.rate=x.obs, del.rate=x.del, size=i
                                                               , is.gap=gap.key[i])
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
    
    ## Create the vectors for the upper and lower bound for nloptr.
    ## The bound is the same for each of the sites.
    if( rate.model == "correlated" ){
        ## Set bounds of 0 to a very small number.
        if( round(bounds[1], digits = 10) == 0.0 ){
            ## Cannot be log(0)
            bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
        }
        ## The third parameter is a correlation
        log_lb <- c( log( bounds[1] ), beta.bounds[1], 0 )
        log_ub <- c( log( bounds[2] ), beta.bounds[2], 1 )
    }
    if( rate.model == "gamma"){
        ## Log the bounds in order to search in log space.
        if( round(bounds[1], digits = 10) == 0.0 ){
            ## Cannot be log(0)
            bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
        }
        log_lb <- c( log( bounds[1] ), beta.bounds[1] )
        log_ub <- c( log( bounds[2] ), beta.bounds[2] )
    }
    if( rate.model == "single.rate" ){
        ## Log the bounds in order to search in log space.
        if( round(bounds[1], digits = 10) == 0.0 ){
            ## Cannot be log(0)
            bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
        }
        log_lb <- log( bounds[1] )
        log_ub <- log( bounds[2] )
    }

    ## Sample the initial parameters for the search.
    ## Here the user can provide a custom start.
    if( is.null(init) ){
        while( TRUE ){
            print( "Sampling starting state..." )
            ## Keep sampling starting states until the sampled state returns a viable likelihood.
            if(rate.model == "gamma"){
                init.pars <- c(  log(runif(1, min=bounds[1], max=bounds[2]))
                             , runif(1, min=beta.bounds[1], max=beta.bounds[2]) )
            }
            if(rate.model == "correlated"){
                while( TRUE ){
                    init.rho <- ifelse(test=init.M, yes=0.5, no=runif(1, min=0, max=1) )
                    init.pars <- c(  log(runif(1, min=bounds[1], max=bounds[2]))
                                 , runif(1, min=beta.bounds[1], max=beta.bounds[2])
                                 , init.rho ) ## The correlation 'rho'
                    ## When the model is correlated, then the M matrix need to be a good matrix on the starting point.
                    ## Need to keep sampling until the starting point is a valid matrix.
                    cat <- qgamma((1:(ncat-1))/ncat, shape = init.pars[2], rate = init.pars[2])
                    cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
                    M_init <- computeM(rate_cat = cat, alpha = init.pars[2], rho = init.rho, k = ncat)
                    ## If the M matrix is acceptable then keep the initial value, otherwise resample.
                    if( any(round(rowSums(M_init), digits = 10) == 1) & any(round(colSums(M_init), digits = 10) == 1) ){
                        break()
                    }
                }
            }
            if(rate.model == "single.rate"){
                init.pars <- log(runif(1, min=bounds[1], max=bounds[2]))
            }
            if( Q.model == "ER" ){
                start.lik <- wrapLogLik(init.pars, n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                      , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                      , root_node=root_node, data=Xlist, k=ncat, root_type=root_type
                                      , gap.key=gap.key)
            } else{
                init.pars.check.start <- c(init.pars[1], init.pars)
                start.lik <- wrapLogLik(init.pars.check.start, n_nodes=n_nodes, n_tips=n_tips
                                      , n_states=n_states
                                      , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                      , root_node=root_node, data=Xlist, k=ncat, root_type=root_type
                                      , gap.key=gap.key)
            }
            if( is.finite( start.lik ) ){
                break()
            }
        }
    } else{
        if( rate.model == "correlated" ){
            if( !length( init ) == 3 ) stop("Wrong number of init parameters.")
            if( init[1] < bounds[1] | init[1] > bounds[2] ) stop("Value for init[1] is out of bounds (defined by 'bounds').")
            if( init[2] < beta.bounds[1] | init[2] > beta.bounds[2] ) stop( paste0("Value for beta (init[2]) is outside bounds. min = ", beta.bounds[1], " and max = ", beta.bounds[2],".") )
            if( init[3] > 1.0 | init[3] < 0.0 ) stop( paste0("Correlation (init[3]) need to be between 0 and 1.") )
            init[1] <- log( init[1] ) ## Tranform the first element. Keep the rest.
            init.pars <- init
        }
        if(rate.model == "gamma"){
            if( !length( init ) == 2 ) stop("Length of init need to be 2. init[1] is for the rate and init[2] is for the Gamma function parameter (beta).")
            if( any(init[1] < bounds[1]) | any(init[1] > bounds[2]) ) stop("Value for init[1] is out of bounds (defined by 'bounds').")
            if( any(init[2] < beta.bounds[1]) | any(init[2] > beta.bounds[2]) ) stop( paste0("Value for beta (init[2]) is outside bounds. min = ", beta.bounds[1], " and max = ", beta.bounds[2],".") )
            init.pars <- c(log(init[1]), init[2])
        }
        if(rate.model == "single.rate"){
            if( !length( init ) == 1 ) stop("Length of init need to be 1. init is for the rate.")
            if( any(init < bounds[1]) | any(init > bounds[2]) ) stop("Value for init is out of bounds (defined by 'bounds').")
            init.pars <- log(init)
        }
    }

    ## The "DEL" model has one more parameter.
    ## Need to adjust the length of init and bounds parameters.
    ## Will assume the bounds for the indels are the same as for the transitions.
    if( Q.model == "DEL" ){
        init.pars <- c(init.pars[1], init.pars)
        log_lb <- c(log_lb[1], log_lb)
        log_ub <- c(log_ub[1], log_ub)
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
        if( search.global ){
            print( "Starting global MLE search. (First pass)" )
            global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                           , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                           , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                           , k=ncat, root_type=root_type, gap.key=gap.key)
            print( "Global search finished." )
            init.pars <- global$solution
        }
        print( "Starting local MLE search. (Second pass)" )
        fit <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=local.opts
                       , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                       , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                    , k=ncat, root_type=root_type, gap.key=gap.key)
        ## Check if the MLE estimation returned good results.
        print( "Reconstructing site-wise Q matrices." )
    } else{
        if( search.global ){
            global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts
                           , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                           , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                           , k=ncat, root_type=root_type, gap.key=gap.key)
            init.pars <- global$solution
        }
        fit <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=local.opts
                       , n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len
                       , edge_mat=edge_mat, parents=parents, root_node=root_node, data=Xlist
                       , k=ncat, root_type=root_type, gap.key=gap.key)
    }
    ## Register finish search time.
    finish.time <- Sys.time()
    total.time <- format( difftime(finish.time, start.time) )

    ## Need to reconstruct the Q matrices for each of the sites.
    ## For this I need to compute the likelihood again.
    ## Apply the original state names to the Q matrices.
    if( Q.model == "ER"){
        solution <- c(exp(fit$solution[1]), fit$solution[2])
        Q <- lapply(n_states, function(i) make.Q.list(rate=solution[1], size=i) )
        beta <- solution[2]
        if( rate.model == "correlated" ){
            ## In this case need to append to the solution object.
            gamma.rates <- discreteGamma(shape = solution[2], ncats = ncat) ## The gamma rates.
            solution <- fit$solution ## Size of 3
            cat <- qgamma((1:(ncat-1))/ncat, shape = solution[2], rate = solution[2])
            cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
            M <- computeM(rate_cat=cat, alpha = solution[2], rho=solution[3], k=ncat)
            ## Need to make sure that the format for the output is the same between models.
            if( !is.matrix(M) ){ ## A single rate category!
                ## In this case only the last category matters.
                res <- lapply(Q, function(x) x * gamma.rates[ncat])
            } else{
                res <- getSiteRatesAutoDiscGamma(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                               , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                               , root_node=root_node, root_type=root_type, k=ncat, X=Xlist
                                               , Q=Q, M=M, beta=beta, n.cores=n.cores)
            }
        }
        if(rate.model == "gamma"){
            res <- loglikGammaSimple_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                     , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                     , root_node=root_node, X = Xlist, Q = Q
                                     , root_type=root_type, beta=beta, k=ncat, n.cores=n.cores)[[2]]
            ## res <- loglikGammaSimple(phy=phy, X=Xlist, Q=Q, root.type=root.type, beta=beta, k=ncat, n.cores=n.cores)[[2]]
        }
        if(rate.model == "single.rate"){
            res <- Q
        }            
    }
    if( Q.model == "DEL"){
        solution <- c(exp(fit$solution[1]), exp(fit$solution[2]), fit$solution[3])
        Q <- lapply(n_states, function(i) make.Q.list.DEL(global.rate=solution[1], del.rate=solution[2]
                                                       , size=i, is.gap=gap.key[i]) )
        beta <- solution[3]
        if( rate.model == "correlated" ){
            ## Need to append the rest of the parameters to the solution object.
            gamma.rates <- discreteGamma(shape = solution[3], ncats = ncat)
            solution <- c(solution, fit$solution[4]) ## Adds the correlation.
            cat <- qgamma((1:(ncat-1))/ncat, shape = solution[3], rate = solution[3])
            cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
            ## This final matrix can have NaNs. But this is not an issue.
            M <- computeM(rate_cat=cat, alpha = solution[3], rho=solution[4], k=ncat)
            ## Need to make sure that the format for the output is the same between models.
            if( !is.matrix(M) ){ ## A single rate category!
                res <- Q * gamma.rates[ncat] ## In this case only the last category matters.
            } else{
                res <- getSiteRatesAutoDiscGamma(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states, edge_len=edge_len, edge_mat=edge_mat, parents=parents, root_node=root_node, root_type=root_type, k=ncat, X=Xlist, Q=Q, M=M, beta=beta, n.cores=n.cores)
            }
        }
        if(rate.model == "gamma"){
            res <- loglikGammaSimple_C(n_nodes=n_nodes, n_tips=n_tips, n_states=n_states
                                     , edge_len=edge_len, edge_mat=edge_mat, parents=parents
                                     , root_node=root_node, X = Xlist, Q = Q
                                     , root_type=root_type, beta=beta, k=ncat, n.cores=n.cores)[[2]]
            ## res <- loglikGammaSimple(phy=phy, X=Xlist, Q=Q, root.type=root.type, beta=beta, k=ncat, n.cores=n.cores)[[2]]
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
    del.rate <- NULL
    corr <- NULL
    auto.matrix <- NULL
    if( rate.model == "correlated" ){
        auto.matrix <- M
        corr <- solution[length(solution)] ## Last in the vector.
    }
    if( rate.model == "single.rate" ){
        start.par <- exp(init.pars[1])
        alpha <- NULL
        if( Q.model == "DEL" ){
            del.rate <- solution[2]
            start.par <- c(exp(init.pars[1]), exp(init.pars[2]))
            alpha <- NULL
        }
    } else{        
        start.par <- c(exp(init.pars[1]),init.pars[2])
        alpha <- solution[2]
        if( Q.model == "DEL" ){
            del.rate <- solution[2]
            start.par <- c(exp(init.pars[1]), exp(init.pars[2]), init.pars[3])
            alpha <- solution[3]
        }
    }
    
    out <- list( log.lik=-fit$objective, Q=res, global.rate=solution[1], corr=corr, auto.matrix=auto.matrix
              , del.rate=del.rate, alpha=alpha, start.par=start.par, nlopt.global.search=global.opts
              , nlopt.local.search=local.opts, nlopt.message=fit$message, search.time=total.time
              , Q.model=Q.model, rate.model=rate.model)
    return( out )
}

## Define a long convoluted function to test for invariant positions in the presence of polymorphisms.
get.invar <- function(x, poly.key){
    ## Function to check if a site is invariant taking into account the existence of polymorphisms.
    ## This is a long series of conditions.
    hard.invar <- length( unique(x) ) == 1
    if( hard.invar ){
        return( TRUE )
    }
    ## Check if 90% of the site is of a single kind.
    max.prop <- max( table(x) )
    if( max.prop / length( x ) > 0.9 ){
        ## Need to check the role of the polymorphic states.
        has.poly <- any( unique(x) %in% names( poly.key ) )
        if( !has.poly ){
            return( FALSE )
        }
        ## Check if all non-polymorphic are the same:
        not.poly <- x[ !x %in% names( poly.key ) ]
        if( length( unique(not.poly) ) > 1 ){
            return( FALSE )
        }                    
        ## Check if the polymorphisms are less than 10% of the data.
        poly.pos.prop <- sum(x %in% names( poly.key )) / length( x )
        if( poly.pos.prop > 0.1 ){
            ## Too many of these. Better to estimate.
            return( FALSE )
        }
        ## Finally, are all polymorphic sites potentially the same state?
        not.poly <- unique( x[ !x %in% names( poly.key ) ] )
        poly.sites <- unique( x[ x %in% names( poly.key ) ] )
        poly.id <- sapply(poly.sites, function(x) which( x == names( poly.key ) ) )
        can.collapse <- sapply(poly.id, function(x) not.poly %in% poly.key[[x]] )
        if( all( can.collapse ) ){
            return( TRUE )
        } else{
            return( FALSE)
        }        
    } else{
        return( FALSE )
    }
}
