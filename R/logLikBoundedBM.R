##' @importFrom ape drop.tip
prepareBoundedBM <- function(dat_vec, phy, trait_bounds){
  ## Function will get the data for a single position and prepare the log likelihood.
  ## Will run this only once per position of the sequence.
  ## This will deal with NAs, such that the function adapt for the data at each position.
  ## dat_vec = named vector with the species traits.
  ## trait_bounds = vector with upper and lower bounds.
  
  if( any( is.na(dat_vec) ) ){
    ## Deal with NAs.
    to_drop <- names( dat_vec[is.na(dat_vec)] )
    dat_vec <- dat_vec[!is.na(dat_vec)]
    phy <- drop.tip(phy = phy, tip = to_drop)
  }
  
  bound_model <- lnL_BBMV(tree = phy, trait = dat_vec
                          , Npts=50, bounds = trait_bounds
                          , a=0, b=0, c=0)
  ## The output is a function. With a single argument "X".
  to_nloptr <- function(x){
    ## Function is ready for nloptr, which is a minimization algorithm.
    return( -1 * bound_model$fun(X = x) )
  }
  return( to_nloptr )
}
