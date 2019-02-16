##' Extracts and summarize the estimates of rates of evolution from the fitted models.
##'
##' Returns a vector if the model has a global rate (i.e., model "ER"). If the estimated model was of type "DEL", this function will return a matrix with the distinct rates for the observed states and between the observed states and gaps.
##'
##' This functions will substitute invariant positions for NA when necessary.
##'
##' Please check the help for 'fitMLGammaSpanSpace' for more information.
##' 
##' @title Extract rates of evolution from fitted models
##' @param x a model fit result.
##' @param Q.model one of "ER" or "DEL". Default is "ER".
##' @param gap.char the charactaed state used to identify gaps. Default is "-".
##' @return a vector or a matrix with the evolutionary rates estimated for each sequence position.
##' @author Daniel Caetano
##' @export
getRates <- function(x, Q.model = "ER", gap.char = "-"){
    ## Gets rates of evolution for the model estimates. Returns a single vector if ER and a matrix with two columns if DEL
    ## Model controls the output if the fitted model is "ER" or "DEL".
    ## In the case of "DEL" we are estimating two rates.
    Q.model <- match.arg(Q.model, choices=c("ER","DEL"), several.ok=FALSE)
    
    if( Q.model == "ER" ){
        rates <- vector(mode = "numeric", length=length(x$Q))
        for( i in 1:length(x$Q) ){
            var.site <- is.matrix(x$Q[[i]])
            if( var.site ){
                ## This is a "ER" model.
                if( ncol(x$Q[[i]]) > 1 ){
                    rates[i] <- x$Q[[i]][1,2]
                }            
            } else{
                rates[i] <- NA
            }
        }
        ## Will return a vector.
        return( rates )
    } else if( Q.model == "DEL"){
        ## Rates will be a matrix with 2 columns.
        rates <- matrix(ncol = 2, nrow = length(x$Q) )
        for( i in 1:length(x$Q) ){
            var.site <- is.matrix(x$Q[[i]])
            if( var.site ){ ## Check if not invariant.

                ## If there is a single col then it is the same as invariant.
                if( ncol(x$Q[[i]]) < 2 ){
                    rates[i,1] <- NA
                    rates[i,2] <- NA
                    next
                }
                   
                ## Check if one of the states is a gap:
                if( "-" %in% rownames(x$Q[[i]]) ){
                    ## Now, we know that the gap state will always be the first row of the Q matrix.
                    ## If the matrix is 2x2, then we only have the gap.rate.
                    ## If the matrix is 3x3 or more, then we have a norm.rate
                    if( ncol( x$Q[[i]] ) < 3 ){
                        rates[i,1] <- x$Q[[i]][1,2]
                        rates[i,2] <- NA
                    } else{
                        rates[i,1] <- x$Q[[i]][1,2] ## Always the first row.
                        rates[i,2] <- x$Q[[i]][2,3] ## At least 3x3.
                    }
                } else{
                    ## We don't have a gap in this position.
                    ## Adding NA to the rate for the gaps.
                    rates[i,1] <- NA
                    rates[i,2] <- x$Q[[i]][1,2] ## This should be fine.
                }
            } else{
                rates[i,1] <- NA
                rates[i,2] <- NA
            }
        }
        ## Will return a vector.
        colnames( rates ) <- c("gap_rate","state_rate")
        return( rates )
    } else{
        stop("Q.Model needs to be one of ER or DEL")
    }
}
