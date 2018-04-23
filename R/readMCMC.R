##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Read the MCMC from files.
##' @param file.name 
##' @param burn 
##' @param thin 
##' @return A data matrix with the samples.q
##' @author daniel
##' @importFrom readr read_lines
##' @export
readMCMC <- function( file.name, burn = 0.5, thin = 10){
    mcmc <- read_lines( file=file.name )
    header <- mcmc[1]
    header <- as.character( strsplit(x=header, split=";", fixed=TRUE)[[1]] )

    obs.gen <- length( mcmc )-1 ## Compute observed number of samples (Fix for unfinished chains.) Note that header is first line.
    if( burn <= 0 ){
        post <- seq(2, obs.gen+1, by=thin) ## First line is the header.
    } else{
        post <- seq(round(obs.gen * burn)+1, obs.gen+1, by=thin) ## First line is the header.
    }
    mcmc <- mcmc[post]

    ## Parse the posterior samples.
    mcmc <- sapply(mcmc, function(x) as.numeric( strsplit(x=x, split=";", fixed=TRUE)[[1]] )
                 , USE.NAMES=FALSE)
    mcmc <- t( mcmc )
    colnames( mcmc ) <- header

    ## Return the table. Still need to make the part to parse the thing.
    return( mcmc )    
}
