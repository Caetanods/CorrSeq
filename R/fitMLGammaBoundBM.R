## Bound BM likelihood uses transformed parameters. Defining two simple functions:
to_bound_BM_lik <- function(x){
  ## Transforms the BM rate parameter to be used in the lik.
  log(x/2)
}

from_bound_BM_lik <- function(x){
  ## Transforms the BM rate parameter back from the lik.
  2 * exp(x)
}

reconBound <- function(lik_fn, solution, ncat, rate.model, seq_length, n.cores){
  ## M is the matrix for the autocorrelated gamma model. If not provided, then it will be NA. This helps when reporting the results at the end.
  ## seq_length is the length of the sequence.
  ## solution is a vector with length 1, 2, or 3, depending on the model.
  ## Function to make the marginal reconstruction of the rates at the positions of the sequence.
  ## Later needs to fix the computation of the marginal likelihood.
  
  ## Need to transform the first parameter.
  solution[1] <- from_bound_BM_lik(solution[1])
  
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
##' @param BM_bounds numeric matrix with 2 rows and number of columns equal to ncol(data). First row is the lower bound and second row the high bound for each sequence position, matching the order in data. Each value of data should be within the bounds.
##' @param rate.model options are "correlated" , "gamma", and "single.rate". See Details below.
##' @param ncat categories for the gamma function.
##' @param init.M whether the initial state for the M matrix (for the correlated model) should have starting state equal to the gamma model (i.e., equal probabilities for all the transitions) or a random starting point.
##' @param bounds a numeric vector of length 2 with the lower and upper bonds for the BM rates. Note that this is different from the "BM_bounds" matrix, which fix bounds for the state value.
##' @param opts the list of options for nloptr. If NULL it will use the default parameters of this function (not the same as the default for 'nloptr'). See more information in the help page for 'nloptr'.
##' @param search.global whether to perform a global MLE search before the local MLE search. Default is TRUE.
##' @param init set the initial parameters for the MLE search. The length varies depending on the rate.model .
##' @param verbose whether to print information to the screen.
##' @param n.cores number of cores to perform the likelihood evaluation.
##' @return A list with the log-likelihood, initial parameters and the parameter values.
##' @importFrom nloptr nloptr
##' @importFrom ape reorder.phylo
##' @importFrom stats runif 
##' @export
##' @author Daniel Caetano
fitCorrSeqBoundBM <- function(data, phy, BM_bounds, rate.model = "gamma", ncat = 4, init.M = FALSE, bounds = NULL, opts = NULL, search.global = TRUE, init = NULL, verbose = TRUE, n.cores = 1){
  
  rate.model <- match.arg(rate.model, choices=c("correlated", "gamma", "single.rate"), several.ok=FALSE)
  
  ## Check phylogeny and data:
  ## Consider adding step to organize the species names in the data following the phy names.
  if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
  match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
  if( !match.names ) stop("Species names do not match between data and phylogeny!")
  ## Two positions in the trait is not enough for this model:
  if( ncol(data) < 3 ) warning("Less than 3 positions in the sequence. Maybe not enough information.")
  if( ncol(data) == 1 ) stop("Cannot work with a sequence of a single position!")
  
  ## Check if the 'bounds' argument has the correct format:
  if( !is.null( bounds ) ){
    if( !length( bounds ) == 2) stop( "Wrong format for the 'bounds' argument." )
    if( bounds[1] < 0 ) stop( "The lower bound cannot be negative." )
  } else{
    ## Set the bounds of the search to defaults.
    bounds <- c(0,100)
  }
  
  ## Check if the 'BM_bounds' argument has the correct format:
  if( !is.null( BM_bounds ) ){
    if( ncol(BM_bounds) != ncol(data) ) stop( "Wrong format for the 'BM_bounds' argument." )
    ## Check if all upper bounds are higher than lower bounds:
    if( any( BM_bounds[1,] > BM_bounds[2,] ) ) stop( "Error on BM_bounds: Lower bounds need to be lower than higher bounds!" )
  } else{
    stop( "Argument BM_bounds missing and has no default value!" )
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
    beta.bounds <- findMinBeta(ncats = ncat)
  }
  
  ## Make sure the bound for the rates is not 0.0
  if( round(bounds[1], digits = 10) == 0.0 ){
    ## Cannot be log(0)
    bounds[1] <- .Machine$double.eps ## Very small number (smallest possible).
  }
  
  ## Create the vectors for the upper and lower bound for nloptr.
  ## The bound is the same for each of the sites.
  if( rate.model == "correlated" ){
    ## The third parameter is a correlation
    log_lb <- c( to_bound_BM_lik(bounds[1]), beta.bounds[1], 0 )
    log_ub <- c( to_bound_BM_lik(bounds[2]), beta.bounds[2], 1 )
  }
  if( rate.model == "gamma"){
    ## Log the bounds in order to search in log space.
    log_lb <- c( to_bound_BM_lik(bounds[1]), beta.bounds[1] )
    log_ub <- c( to_bound_BM_lik(bounds[2]), beta.bounds[2] )
  }
  if( rate.model == "single.rate" ){
    ## Log the bounds in order to search in log space.
    log_lb <- to_bound_BM_lik( bounds[1] )
    log_ub <- to_bound_BM_lik( bounds[2] )
  }
  
  ## Prepare the bounded BM likelihood function:
  ## This is a list of functions, one for each sequence position.
  lik_BBM_list <- lapply(1:ncol(data), function(j) prepareBoundedBM(dat_vec=data[,j], phy=phy, trait_bounds=BM_bounds[,j]) )
  
  ## Define the log likelihood functions. One of each model.
  if( rate.model == "correlated" ){
    wrapLogLik <- function(obj){
      ## obj is a vector of variable length.
      ## obj[1] = Q[1,2] (shared for all sites), obj[2] = beta, obj[3] = correlation for the bivariate Gamma.
      ## Note that the BM rate will be transformed, so need to a) untransform, b) scale by the gamma, c) transform back.
      bm <- from_bound_BM_lik(obj[1]) ## Searching for the rate of transition in log space.
      beta <- obj[2]
      rho <- obj[3] ## This is the correlation of the bivariate Gamma.
      
      ## The M matrix for the autocorrelation:
      cat <- qgamma((1:(ncat - 1))/ncat, shape = beta, rate = beta)
      cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
      ## alpha and beta are the same thing.
      M <- computeM(rate_cat = cat, alpha = beta, rho = rho, k = ncat)
      if( !is.matrix(M) ){ ## A single rate category!
        ## The likelihood of the model considering only the last category.
        ## Because beta is such that categories 1 to k-1 are empty (width of 0).
        gamma.rates <- discreteGamma(shape = beta, ncats = ncat) ## These are the k rates.
        scaled_bm <- gamma.rates[ncat] * bm
        lik <- logLikBoundedBMSimple(lik_fn = lik_BBM_list, bm = scaled_bm, n.cores = n.cores)
        return( -1 * lik )
      } else{
        ## Check if M is a doubly stochastic matrix. Otherwise, reject.
        if( any(round(rowSums(M), digits=10) != 1.0) | any(round(colSums(M), digits=10) != 1.0) ){
          ## Bad M matrix, this will happen sometimes. Return bad lik.
          return( Inf )
        }
        ## Loglik function for the model.
        lik <- logLikBoundedBMAutoGamma(lik_fn = lik_BBM_list, bm = bm, M = M, beta = beta
                                        , k = ncat, n.cores = n.cores)
      }
      return( -1 * lik ) ## Remember that NLOPT is minimizying the function!
    }
  }
  if( rate.model == "gamma"){
    wrapLogLik <- function(obj){
      bm <- from_bound_BM_lik(obj[1]) ## Searching for the rate of transition in log space.
      beta <- obj[2]
      lik <- logLikBoundedBMSimpleGamma(lik_fn = lik_BBM_list, bm = bm, beta = beta, k = ncat
                                        , n.cores = n.cores)
      return( -1 * lik )
    }
  }
  if( rate.model == "single.rate" ){
    wrapLogLik <- function(obj){
      bm <- from_bound_BM_lik(obj[1]) ## Searching for the rate of transition in log space.
      lik <- logLikBoundedBMSimple(lik_fn = lik_BBM_list, bm = bm, n.cores = n.cores)
      return( -1 * lik )
    }
  }
  
  ## Sample the initial parameters for the search.
  ## Here the user can provide a custom start.
  if( is.null(init) ){
    while( TRUE ){
      print( "Sampling starting state..." )
      ## Keep sampling starting states until the sampled state returns a viable likelihood.
      if(rate.model == "gamma"){
        init.pars <- c(  to_bound_BM_lik( runif(1, min=bounds[1], max=bounds[2]) )
                         , runif(1, min=beta.bounds[1], max=beta.bounds[2]) )
      }
      if(rate.model == "correlated"){
        while( TRUE ){
          init.rho <- ifelse(test=init.M, yes=0.5, no=runif(1, min=0, max=1) )
          init.pars <- c(  to_bound_BM_lik( runif(1, min=bounds[1], max=bounds[2]) )
                           , runif(1, min=beta.bounds[1], max=beta.bounds[2])
                           , init.rho ) ## The correlation 'rho'
          ## When the model is correlated, then the M matrix needs to be a good matrix on the starting point.
          ## Need to keep sampling until the starting point is a valid matrix.
          cat <- qgamma((1:(ncat-1))/ncat, shape = init.pars[2], rate = init.pars[2])
          cat <- c(.Machine$double.eps, cat, Inf) ## These are the bounds of the categories.
          M_init <- computeM(rate_cat = cat, alpha = init.pars[2], rho = init.rho, k = ncat)
          ## If the M matrix is acceptable then keep the initial value, otherwise resample.
          if( all(round(rowSums(M_init), digits = 10) == 1) & all(round(colSums(M_init), digits = 10) == 1) ){
            break()
          }
        }
      }
      if(rate.model == "single.rate"){
        init.pars <- to_bound_BM_lik(runif(1, min=bounds[1], max=bounds[2]))
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
      init[1] <- to_bound_BM_lik( init[1] ) ## Transform the first element. Keep the rest.
      init.pars <- init
    }
    if(rate.model == "gamma"){
      if( !length( init ) == 2 ) stop("Wrong number of init parameters. Length of init need to be 2. init[1] is for the rate and init[2] is for the Gamma function parameter (beta).")
      if( any(init[1] < bounds[1]) | any(init[1] > bounds[2]) ) stop("Value for init[1] is out of bounds (defined by 'bounds').")
      if( any(init[2] < beta.bounds[1]) | any(init[2] > beta.bounds[2]) ) stop( paste0("Value for beta (init[2]) is outside bounds. min = ", beta.bounds[1], " and max = ", beta.bounds[2],".") )
      init[1] <- to_bound_BM_lik( init[1] ) ## Transform the first element. Keep the rest.
      init.pars <- init
    }
    if(rate.model == "single.rate"){
      if( !length( init ) == 1 ) stop("Wrong number of init parameters. Length of init need to be 1. init is for the rate.")
      if( any(init < bounds[1]) | any(init > bounds[2]) ) stop("Value for init is out of bounds (defined by 'bounds').")
      init.pars <- to_bound_BM_lik(init)
    }
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
  }
  if( verbose ) print( "Starting local MLE search. (Second pass)" )
  fit <- nloptr(x0=init.pars, eval_f=wrapLogLik, lb=log_lb, ub=log_ub, opts=local.opts)
  ## The BM parameter needs to be transformed at the end.
  ## fit$x0[1] <- from_bound_BM_lik(fit$x0[1])
  
  if( verbose ) print( "Reconstructing site-wise Q matrices." )
  
  solution <- fit$x0
  solution[1] <- from_bound_BM_lik(solution[1])

  recon_par <- reconBound(lik_fn=lik_BBM_list, solution=solution, ncat=ncat, rate.model=rate.model
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
  out <- list( log.lik=-fit$objective, global.rate=solution[1], corr=solution[3], auto.matrix=M
               , alpha=solution[2], start.par=init.pars, recon.rates = recon_par
               , nlopt.global.search=global.opts, nlopt.local.search=local.opts
               , nlopt.message=fit$message, search.time=total.time, rate.model=rate.model )
  ## Give a custom class so we can check for it later.
  class( out ) <- append( class( out ), "bm_seqrates")
  return( out )
}

##' @export
print.bm_seqrates <- function(x){
  if( !inherits(x, "bm_seqrates") ) stop("x is not of class bm_seqrates")
  cat("Type of model:", format(x$rate.model), "\n")
  cat("Loglik:", format(x$log.lik), "\n")
  cat("Bounded BM rate:", format(x$global.rate), "\n")
  cat("Gamma rate:", format(x$alpha), "\n")
  cat("Autocorrelation parameter:", format(x$corr), "\n")
  cat("This is a list type object. Check 'names()' for more information. \n")
}
