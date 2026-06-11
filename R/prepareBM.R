##' @importFrom ape drop.tip vcv.phylo Ntip
prepareBM <- function(dat_vec, phy, se){
  ## Function will get the data for a single position and prepare the log likelihood.
  ## Will run this only once per position of the sequence.
  ## This will deal with NAs, such that the function adapt for the data at each position.
  ## dat_vec = named vector with the species traits.
  ## phy = phylogeny
  ## se = standard error
  
  if( any( is.na(dat_vec) ) ){
    ## Deal with NAs.
    to_drop <- names( dat_vec[is.na(dat_vec)] )
    dat_vec <- dat_vec[!is.na(dat_vec)]
    phy <- drop.tip(phy = phy, tip = to_drop)
  }
  
  # se is a vector of the length of the tip data and with names.
  if(is.null(se)){
    ## Make vector of zeros.
    se <- setNames(object = rep(0, length(dat_vec)), nm = names(dat_vec))
  }
  
  C <- vcv.phylo(phy)
  n <- Ntip(phy)
  E <- diag(se^2) # Make the error matrix from the standard error.
  
  logLikBMtmp <- function(sigma, mu){
    logL <- as.numeric(-t(dat_vec-mu) %*% solve(sigma * C + E) %*% (dat_vec-mu)/2 -
                         n * log(2*pi) / 2 - determinant(sigma * C + E)$modulus[1]/2)
    return(logL)
  }

  return( logLikBMtmp )
}

prepareBMstart <- function(dat_vec, phy){
  ## Use the vector of the data and the phylogeny to get informative starting states.
  
  if( any( is.na(dat_vec) ) ){
    ## Deal with NAs.
    to_drop <- names( dat_vec[is.na(dat_vec)] )
    dat_vec <- dat_vec[!is.na(dat_vec)]
    phy <- drop.tip(phy = phy, tip = to_drop)
  }
  
  C <- vcv.phylo(phy)
  n <- Ntip(phy)
  
  # Use the analytical solution as the starting state for the search:
  mu_start <- as.numeric( colSums(solve(C)) %*% ( dat_vec/sum(solve(C)) ) )
  sigma_start <- as.numeric( t(dat_vec-mu_start) %*% solve(C) %*% (dat_vec-mu_start)/n )
  
  return( c(sigma_start, mu_start) )
}
