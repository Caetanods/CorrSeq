reconBM <- function(lik_fn, solution, ncat, rate.model, seq_length, n.cores){
  ## M is the matrix for the autocorrelated gamma model. If not provided, then it will be NA. This helps when reporting the results at the end.
  ## seq_length is the length of the sequence.
  ## solution is a vector with length 1, 2, or 3, depending on the model.
  ## Function to make the marginal reconstruction of the rates at the positions of the sequence.
  ## Later needs to fix the computation of the marginal likelihood.
  
  if( rate.model == "correlated" ){
    ## In this case need to append to the solution object.
    gamma.rates <- discreteGamma(shape = solution[2], ncats = ncat) ## The gamma rates.
    cat <- qgamma((1:(ncat-1))/ncat, shape = solution[2], rate = solution[2])
    cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
    M <- computeM(rate_cat=cat, alpha = solution[2], rho=solution[3], k=ncat)
    ## Need to make sure that the format for the output is the same between models.
    if( !is.matrix(M) ){ ## A single rate category!
      ## In this case only the last category matters.
      res_tmp <- solution[1] * gamma.rates[ncat]
      res <- rep(res_tmp, times = seq_length)
    } else{
      res <- getMarginalBMAutoDiscGamma(lik_fn=lik_fn, solution=solution, k=ncat
                                        , M=M, nsites=seq_length, n.cores=n.cores)
    }
  }
  if(rate.model == "gamma"){
    res <- getMarginalBMSimpleGamma(lik_fn=lik_fn, solution=solution, k=ncat, n.cores=n.cores)
  }
  if(rate.model == "single.rate"){
    res <- rep(solution[1], times = seq_length)
  }
  
  ## Return the vector of marginal estimates of BM rates across the sequence positions.
  return( res )
}

##' Function uses the discrete-gamma model of Yang (1994) to estimate variation in rates of evolution independent of sequence position and the auto-discrete-gamma model of Yang (1995) when rates are correlated among sites. This function takes a sequence trait with continuous states and discrete positions. The approach fits a bounded Brownian motion model to the data, with bounds defined by the user. Bounds can be estimated, but that is not recommended. Each sequence position can have their own bound.
##' 
##' @title Maximum likelihood estimate of evolutionary rates for sequence phenotype with continuous states and discrete positions
##' @param data matrix with species names as rownames and sequence positions as the columns. Each sequence position should be a numeric value, within the bounds of the "BM_bounds" matrix.
##' @param phy ultrametric phylogeny in ape 'phylo' format
##' @param se the standard error of the tip data. This should be a named vector with length equal to tip data with the measurement error for each species. Alternatively, this can be a single value that will be used for all species. Set to NULL if unknown.
##' @param rate.model options are "correlated" , "gamma", and "single.rate". See Details below.
##' @param ncat categories for the gamma function.
##' @param init.M whether the initial state for the M matrix (for the correlated model) should have starting state equal to the gamma model (i.e., equal probabilities for all the transitions) or a random starting point.
##' @param bounds a numeric vector of length 2 with the lower and upper bonds for the BM rates.
##' @param opts the list of options for nloptr. If NULL it will use the default parameters of this function (not the same as the default for 'nloptr'). See more information in the help page for 'nloptr'.
##' @param search.global whether to perform a global MLE search before the local MLE search. Default is FALSE.
##' @param init set the initial parameters for the MLE search. The length varies depending on the rate.model .
##' @param verbose whether to print information to the screen.
##' @param n.cores number of cores to perform the likelihood evaluation.
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @importFrom ape reorder.phylo
##' @importFrom stats runif 
##' @export
##' @author Daniel Caetano
fitCorrSeqBM <- function(data, phy, se = NULL, rate.model = "gamma", ncat = 4, init.M = FALSE, bounds = NULL, opts = NULL, search.global = FALSE, init = NULL, verbose = TRUE, n.cores = 1){
  
  rate.model <- match.arg(rate.model, choices=c("correlated", "gamma", "single.rate"), several.ok=FALSE)
  
  ## Check phylogeny and data:
  if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
  match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
  if( !match.names ) stop("Species names do not match between data and phylogeny!")
  ## Two positions in the trait is not enough for this model:
  if( ncol(data) < 3 ) warning("Less than 3 positions in the sequence. Maybe not enough information.")
  if( ncol(data) == 1 ) stop("Cannot work with a sequence of a single position!")
  
  ## Check se:
  if( !is.null(se) ){
    if( length(se) == 1 ){
      se <- setNames(object = rep(se, times = Ntip(phy)), nm = phy$tip.label)
    } else if( length(se) == nrow(data) ){
      if( is.null( names(se) ) ) stop("se need to have rownames as the species names.")
      match.names <- all( names(se) %in% phy$tip.label ) & all( phy$tip.label %in% names(se) )
      if( !match.names ) stop("Species names do not match between data and se vector!")
      se <- se[phy$tip.label]
    } else{
      stop("se needs to be either a named vector with length equal to the number of tips in the data or a single value to be repeated to all species.")
    }
  }
  
  ## Check if the 'bounds' argument has the correct format:
  if( !is.null( bounds ) ){
    if( !length( bounds ) == 2) stop( "Wrong format for the 'bounds' argument." )
    if( bounds[1] < .Machine$double.xmin ) stop( "The lower bound cannot be negative." )
  } else{
    ## Set the bounds of the search to defaults.
    bounds <- c(.Machine$double.xmin, 100)
  }
  
  ## Re-order the species in the data matrix to match the tree:
  data.order <- match(x=phy$tip.label, table=rownames(data))
  data <- data[data.order,]
  
  ## Make data checks and get information from the matrix.
  nsites <- ncol(data)
  names.data <- rownames(data)
  
  ## Search for upper and lower bounds for the beta parameter for the Gamma rates that do not produce 0 rate values.
  if( rate.model %in% c("correlated","gamma") ){
    ## The single rate model does not have a beta parameter.
    beta.bounds <- findMinBeta(ncats = ncat) # Upper bound is set to 1000
  }
  
  ## Make sure the bound for the rates is not 0.0
  ## Avoid "==" tests with numeric values.
  if( bounds[1] < .Machine$double.eps ){
    ## Cannot be log(0)
    bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
  }
  
  ## Create the vectors for the upper and lower bound for nloptr global search.
  ## The bound is the same for each of the sites.
  ## NOTE: bounds are in normal space.
  if( rate.model == "correlated" ){
    ## The third parameter is a correlation
    log_lb <- c(bounds[1], beta.bounds[1], 0)
    log_ub <- c(bounds[2], beta.bounds[2], 1)
  }
  if( rate.model == "gamma"){
    log_lb <- c(bounds[1], beta.bounds[1])
    log_ub <- c(bounds[2], beta.bounds[2])
  }
  if( rate.model == "single.rate" ){
    log_lb <- c(bounds[1])
    log_ub <- c(bounds[2])
  }
  
  ## Prepare the BM likelihood function:
  ## This is a list of functions, one for each sequence position.
  ## NOTE: If se is absent, then it should be NULL
  lik_BM_list <- lapply(1:ncol(data), function(j) prepareBM(dat_vec=data[,j], phy=phy, se=se))
  ## Make a matrix of starting states across all positions.
  ## data vectors on the columns. Two rows.
  start_par <- sapply(1:ncol(data), function(j) prepareBMstart(dat_vec=data[,j], phy=phy))
  sigma.mean.init <- mean(start_par) ## Average MLE sigma across positions.
  
  ## Define the log likelihood functions. One of each model.
  if( rate.model == "correlated" ){
    wrapLogLik <- function(obj){
      ## obj is a vector of variable length.
      ## obj[1] = sigma, obj[2] = correlation for the bivariate Gamma, obj[3] = Gamma rate 
      sigma <- obj[1]
      beta <- obj[2]
      rho <- obj[3]
      
      ## The M matrix for the autocorrelation:
      cat <- qgamma((1:(ncat - 1))/ncat, shape = beta, rate = beta)
      cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
      ## alpha and beta are the same thing.
      
      M <- computeM(rate_cat = cat, alpha = beta, rho = rho, k = ncat)
      if( !is.matrix(M) ){ ## A single rate category!
        ## The likelihood of the model considering only the last category.
        ## Because beta is such that categories 1 to k-1 are empty (width of 0).
        gamma.rates <- discreteGamma(shape = beta, ncats = ncat) ## These are the k rates.
        scaled_sigma <- gamma.rates[ncat] * sigma
        lik <- logLikBMSimple(lik_fn = lik_BM_list, sigma = scaled_sigma, n.cores = n.cores)
        return( -1 * lik ) ## Inverted for NLOPT minimization.
      } else{
        ## Check if M is a doubly stochastic matrix. Otherwise, reject.
        rowCheck <- any( sapply(rowSums(M), function(x) !isTRUE(all.equal(x, 1.0))) )
        colCheck <- any( sapply(colSums(M), function(x) !isTRUE(all.equal(x, 1.0))) )
        if( rowCheck | colCheck ){
          ## Bad M matrix, this will happen sometimes. Return bad lik.
          return( Inf ) 
        }
        ## Loglik function for the model.
        lik <- logLikBMAutoGamma(lik_fn = lik_BM_list, sigma = sigma, M = M, beta = beta
                                 , k = ncat, n.cores = n.cores)
        return( -1 * lik ) ## Inverted for NLOPT minimization.
      }
    }
  }
  if( rate.model == "gamma"){
    wrapLogLik <- function(obj){
      sigma <- obj[1]
      beta <- obj[2]
      lik <- logLikBMSimpleGamma(lik_fn = lik_BM_list, sigma = sigma, beta = beta, k = ncat
                                 , n.cores = n.cores)
      return( -1 * lik )
    }
  }
  if( rate.model == "single.rate" ){
    wrapLogLik <- function(obj){
      sigma <- obj[1]
      lik <- logLikBMSimple(lik_fn = lik_BM_list, sigma = sigma, n.cores = n.cores)
      return( -1 * lik )
    }
  }
  
  ## Sample the initial parameters for the search.
  ## Here the user can provide a custom start.
  if( is.null(init) ){
    if( verbose ) print( "Sampling starting state..." )
    while( TRUE ){
      ## Keep sampling starting states until the sampled state returns a viable likelihood.
      if(rate.model == "gamma"){
        init.pars <- c(sigma.mean.init
                       , runif(1, min=beta.bounds[1], max=beta.bounds[2]))
      }
      if(rate.model == "correlated"){
        while( TRUE ){
          init.rho <- ifelse(test=init.M, yes=0.5, no=runif(1, min=0, max=1) )
          init.pars <- c(sigma.mean.init
                         , runif(1, min=beta.bounds[1], max=beta.bounds[2])
                         , init.rho)
          ## When the model is correlated, then the M matrix needs to be a good matrix on the starting point.
          ## Need to keep sampling until the starting point is a valid matrix.
          cat <- qgamma((1:(ncat-1))/ncat, shape = init.pars[2], rate = init.pars[2])
          cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
          M_init <- computeM(rate_cat = cat, alpha = init.pars[2], rho = init.rho, k = ncat)
          ## If the M matrix is acceptable then keep the initial value, otherwise resample.
          rowCheck <- all( sapply(rowSums(M_init), function(x) isTRUE(all.equal(x, 1.0))) )
          colCheck <- all( sapply(colSums(M_init), function(x) isTRUE(all.equal(x, 1.0))) )
          if( rowCheck & colCheck ){
            break
          }
        }
      }
      if(rate.model == "single.rate"){
        init.pars <- sigma.mean.init
      }
      ## Evaluate the log lik for the model.
      start.lik <- wrapLogLik(init.pars)
      ## Check if it is valid, if yes, then break.
      if( is.finite( start.lik ) ){
        break()
      }
    }
  } else{
    if( rate.model == "correlated" ){
      if( !length( init ) == 3 ) stop("Wrong number of init parameters. Length of init needs to be 3. init[1] is for the rate, init[2] is for the Gamma function parameter (beta), and init[3] is for the correlation of the bivariate Gamma.")
      if( init[1] < bounds[1] | init[1] > bounds[2] ) stop("Value for init[1] is out of bounds (defined by 'bounds').")
      if( init[2] < beta.bounds[1] | init[2] > beta.bounds[2] ) stop( paste0("Value for beta (init[2]) is outside bounds. min = ", beta.bounds[1], " and max = ", beta.bounds[2],".") )
      if( init[3] > 1.0 | init[3] < 0.0 ) stop( paste0("Correlation (init[3]) needs to be between 0 and 1.") )
    }
    if(rate.model == "gamma"){
      if( !length( init ) == 2 ) stop("Wrong number of init parameters. Length of init need to be 2. init[1] is for the rate and init[2] is for the Gamma function parameter (beta).")
      if( any(init[1] < bounds[1]) | any(init[1] > bounds[2]) ) stop("Value for init[1] is out of bounds (defined by 'bounds').")
      if( any(init[2] < beta.bounds[1]) | any(init[2] > beta.bounds[2]) ) stop( paste0("Value for beta (init[2]) is outside bounds. min = ", beta.bounds[1], " and max = ", beta.bounds[2],".") )
    }
    if(rate.model == "single.rate"){
      if( !length( init ) == 1 ) stop("Wrong number of init parameters. Length of init need to be 1. init is for the rate.")
      if( any(init < bounds[1]) | any(init > bounds[2]) ) stop("Value for init is out of bounds (defined by 'bounds').")
    }
    init.pars <- init
    ## Evaluate the log-lik for the model under the start parameters:
    start.lik <- wrapLogLik(init.pars)
    ## Check if it is valid.
    if( !is.finite( start.lik ) ) stop("Likelihood is not finite under the chosen starting values.")
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
  if( search.global ){
    if( verbose ) print( "Starting global MLE search. (First pass)" )
    global <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=global.opts)
    if( verbose ) print( "Global search finished." )
    init.pars <- global$solution
    if( verbose ) print( "Starting local MLE search. (Second pass)" )
  } else{
    if( verbose ) print( "Starting MLE search." )  
  }
  ## Now do the local search.
  fit <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=local.opts)

  if( verbose ) print( "Reconstructing site-wise BM rates." )
  
  solution <- fit$solution
  
  recon_par <- reconBM(lik_fn=lik_BM_list, solution=solution, ncat=ncat, rate.model=rate.model
                       , seq_length=nsites, n.cores=n.cores)
  
  ## Register finish search time.
  finish.time <- Sys.time()
  total.time <- format( difftime(finish.time, start.time) )
  
  ## Complete output for the function:
  if( rate.model == "correlated" ){
    M <- computeM(rate_cat = cat, alpha = solution[2], rho = solution[3], k = ncat)
  } else{
    M <- NA
  }
  
  ## Get the root values to output as the solution.
  get_root <- function(x, n, Cinv){
    one_vec <- rep(1, times = n)
    mu <- (one_vec %*% Cinv %*% x) / (one_vec %*% Cinv %*% one_vec)
    return( mu )
  }
  Cinv <- solve( vcv.phylo(phy) )
  n <- Ntip(phy)
  mu_solution <- vector(mode = "numeric", length = ncol(data))
  for( i in 1:ncol(data) ){
    mu_solution[i] <- get_root(x = data[,i], n = n, Cinv)
  }
  
  ## Output from the model also depends on the type of model.
  if( rate.model == "correlated" ){
    sigma <- solution[1]
    beta <- solution[2]
    rho <- solution[3]
    mu <- mu_solution
  } else if( rate.model == "gamma"){
    sigma <- solution[1]
    beta <- solution[2]
    rho <- NA
    mu <- mu_solution
  } else{ # single-rate
    sigma <- solution[1]
    beta <- NA
    rho <- NA
    mu <- mu_solution
  }
  
  ## Compute AIC and AICc:
  npar <- length( solution ) + length( mu_solution )
  ntips <- length( phy$tip.label )
  loglik <- -1 * fit$objective
  ## AIC = -2 ( ln ( likelihood )) + 2 K
  AIC <- (-2 * loglik) + (2 * npar)
  ## AICc = -2 ( ln ( likelihood )) + 2 K * (n / ( n - K - 1))
  AICc <- (-2 * loglik) + (2 * npar * (ntips / ( ntips - npar - 1)) )
  
  out <- list( log.lik= -1 * fit$objective, AIC = AIC, AICc = AICc
               , global.rate=sigma, root.values=mu
               , corr=rho, auto.matrix=M, alpha=beta, start.par=init.pars
               , recon.rates = recon_par, nlopt.global.search=global.opts
               , nlopt.local.search=local.opts, nlopt.message=fit$message
               , search.time=total.time, rate.model=rate.model )
  ## Give a custom class so we can check for it later.
  class( out ) <- append( class( out ), "bm_seqrates")
  return( out )
}

##' @export
print.bm_seqrates <- function(x){
  if( !inherits(x, "bm_seqrates") ) stop("x is not of class bm_seqrates")
  cat("Type of model:", format(x$rate.model), "\n")
  cat("Loglik:", format(x$log.lik), "\n")
  cat("AIC:", format(x$AIC), "\n")
  cat("AICc:", format(x$AICc), "\n")
  cat("BM rate:", format(x$global.rate), "\n")
  cat("Gamma rate:", format(x$alpha), "\n")
  cat("Autocorrelation parameter:", format(x$corr), "\n")
  cat("This is a list type object. Check 'names()' for more information. \n")
}
