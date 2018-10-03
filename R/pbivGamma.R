pbivGamma <- function(x, y, alpha, corr, log.p=TRUE){
    ## Compute the CDF for a bivariate Gamma distribution.
    ## This follows Yang, 1995 and the approach described in Song, 2000 which relies on the copula of a Gaussian distribution to the Gamma in order to compute the CDF.

    cdf <- vector(mode = "numeric", length = 2)
    cdf[1] <- pgamma(q=x, shape=alpha, rate=alpha)
    cdf[2] <- pgamma(q=y, shape=alpha, rate=alpha)
    normcdf <- qnorm(cdf) ## Original code do some bound opperation to protect here.
    
    ## calculate copula contribution to log-PDF    
    res.copula.add <- log( pbivnorm::pbivnorm(x = normcdf[1], y = normcdf[2], rho = corr) ) ## Not in log.
    res.copula.sub <- sum( pnorm(normcdf, log.p=TRUE) )
    res.copula <- res.copula.add - res.copula.sub ## Division of the probabilities.
    
    ## calculate marginal contribution to log-PDF    
    res.data <- sum( log(cdf) )
    
    ## final calculations and return    
    retn <- res.copula + res.data
    
    if(log.p){
        return(retn)
    }else {
        return( exp(retn) )
    }
    
}
