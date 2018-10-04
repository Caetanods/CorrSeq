getAIC <- function(loglik, k){
    return( -2 * loglik + (2 * k) )
}

getAICc <- function(loglik, k, n){
    return( -2 * loglik + (2 * k * ( n / (n - k -1) ) ) )
}
