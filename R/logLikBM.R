logLikBM <- function(sigma, mu, C, x, e){
  # sigma the rate of the BM
  # mu the root value
  # C the phylogenetic variance-covariance matrix
  # x a vector of the tip values
  # e the standard error, vector of zeros if absent.

  n <- length(x)
  E <- diag(e^2) # Make the error matrix from the standard error.
  logL <- as.numeric(-t(x-mu) %*% solve(sigma * C + E) %*% (x-mu)/2 -
                       n * log(2*pi) / 2 - determinant(sigma * C + E)$modulus[1]/2)
  return(logL) # Returning the actual likelihood. No inversion.
}
