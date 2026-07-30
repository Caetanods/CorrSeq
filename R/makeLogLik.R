## General wrapping function that takes data and outputs the log-likelihood for the choosen model.
## This has the same log-likelihood used in the estimation functions.
## Output can be used for ML and Bayesian estimation.

#' Likelihood function for the models
#'
#' @param data matrix with species names as rownames and sequence positions as the columns.
#' @param phy ultrametric phylogeny in ape 'phylo' format
#' @param se the standard error of the tip data. This should be a named vector with length equal to tip data with the measurement error for each species. Alternatively, this can be a single value that will be used for all species. Set to NULL if unknown.
#' @param model.type the type of model for the likelihood function. Options are "MK_discrete" for data with discrete states and discrete positions or "BM_discrete" for data with continuous states and discrete positions.
#' @param rate.model options are "correlated" , "gamma", and "single.rate". See Details.
#' @param ncat number of categories for the discrete Gamma distribution
#' @param verbose if informative messages should be printed.
#' @param n.cores number of cores used for parallel computation of the likelihood.
#'
#' @returns a likelihood function.
#' @details
#' This is a general likelihood maker that will return a likelihood function for different models implemented in the package. The rate.model controls the type of rate models. The options are "correlated" for models using the bi-variate Discrete Gamma model (auto-discrete-gamma, Yang 1995), "gamma" for models using the independent Discrete Gamma distribution, and "single-rate" for models in which there is a single rate for all sequence positions.
#' 
#' @export
makeLogLik <- function(data, phy, se = NULL, model.type, rate.model, ncat = NULL, n.cores = 1, verbose = TRUE){

  #### Check block ####
  model.type <- match.arg(model.type, choices=c("MK_discrete", "BM_discrete"), several.ok=FALSE)
  rate.model <- match.arg(rate.model, choices=c("correlated", "gamma", "single.rate"), several.ok=FALSE)
  if(rate.model == "correlated" | rate.model == "gamma"){
    if( is.null(ncat) ) stop("ncat needs to be an integer value for rate.model 'correlated' and 'gamma'.")
  }
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
  
  ## Re-order the species in the data matrix to match the tree:
  data.order <- match(x=phy$tip.label, table=rownames(data))
  data <- data[data.order,]
  
  ## Make data checks and get information from the matrix.
  nsites <- ncol(data)
  names.data <- rownames(data)
  
  if( model.type == "BM_discrete" ){
    #### Prepare the BM likelihood function ####
    
    ## This is a list of functions, one for each sequence position.
    ## NOTE: If se is absent, then it should be NULL
    lik_BM_list <- lapply(1:ncol(data), function(j) prepareBM(dat_vec=data[,j], phy=phy, se=se))
    ## Make a matrix of starting states across all positions.
    ## data vectors on the columns. Two rows.
    start_par <- sapply(1:ncol(data), function(j) prepareBMstart(dat_vec=data[,j], phy=phy))
    sigma.mean.init <- mean(start_par) ## Average MLE sigma across positions.
    
    ## Define the log likelihood functions. One of each model.
    if( rate.model == "correlated" ){
      if( verbose ) print( "Function parameter is a vector with length 3 containing the values for sigma (BM rate), alpha (Gamma distribution parameter), and rho (autocorrelation parameter)." )
      
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
          return( lik ) ## Inverted for NLOPT minimization.
        }
      }
    }
    if( rate.model == "gamma"){
      if( verbose ) print( "Function parameter is a vector with length 2 containing the values for sigma (BM rate) and alpha (Gamma distribution parameter)." )
      
      wrapLogLik <- function(obj){
        sigma <- obj[1]
        beta <- obj[2]
        lik <- logLikBMSimpleGamma(lik_fn = lik_BM_list, sigma = sigma, beta = beta, k = ncat
                                   , n.cores = n.cores)
        return( lik )
      }
    }
    if( rate.model == "single.rate" ){
      if( verbose ) print( "Function parameter is sigma (BM rate)." )
      
      wrapLogLik <- function(obj){
        sigma <- obj[1]
        lik <- logLikBMSimple(lik_fn = lik_BM_list, sigma = sigma, n.cores = n.cores)
        return( lik )
      }
    }
    
  }
  
  ## Return the likelihood function.
  return( wrapLogLik )
}