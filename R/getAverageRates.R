##' Function performs model averaging given the AIC weights of the models. This follows a similar procedure to the one implemented in Caetano et al., 2018.
##'
##' Model averaging is performed by multiplying the rate estimated under each of the models by the respective model weight. Model weight is a measure of the relative likelihood of a model. This can be used to associate a relative likelihood to parameter estimates. See more information on Caetano et al., 2018.
##' 
##' @title Model average rates along the positions of the sequence trait
##' @param ... a list of model fit or a series of model fits separated by commas.
##' @param verbose whether to print information to the screen.
##' @param legacy.mode whether to check structure of objects instead of class.
##' @return A list object with a vector of rates, a vector of flags from the model averaging algorithm and a named vector with the AICw for the models.
##' @export
##' @author Daniel Caetano
getAverageRates <- function(..., verbose = TRUE, legacy.mode = FALSE){
    ## Argument "legacy.mode" is used here to run the function in "legacy mode".
    models <- list(...)

    ## Check class of models:
    if( legacy.mode ){
        ## Used for the results of previous version of the model.
        if( verbose ){
            print( "Skipping some format tests. Running on legacy mode." )
        }
        check.model.type <- sapply(models, function(x) "Q.model" %in% names(x) )
        if( !all( check.model.type ) ){
            check.model.type.nested <- sapply(models[[1]], function(x) "Q.model" %in% names(x) )
            if( !all( check.model.type.nested ) ){
                ## Then it is not the correct format.
                stop("Some or all of the models are in the incorrect format. Double check.")
            } else{
                ## Should be just fine.
                models <- models[[1]]
            }
        } else{
            ## All good!
        }        
    } else{
        ## Used for normal results.
        check.class <- sapply(models, function(x) "seqrates" %in% class(x) )
        if( !all( check.class ) ){
            ## Check if we get the correct object by unlisting:
            check.class <- sapply(models[[1]], function(x) "seqrates" %in% class(x) )
            if( !all( check.class ) ){
                stop( "Objects provided to the function need to have class 'seqrates'." )
            } else{
                ## Then we might be just fine.
                models <- models[[1]]
            }
        }
    }

    ## Get info from each of the models:
    Q.model <- sapply(models, function(x) x$Q.model)

    ## Just as an extra defense line:
    if( !all( sapply(Q.model, function(x) x %in% c("ER","DEL") ) ) ){
        stop("Some or all of the models are in the incorrect format. Double check.")
    }    
    
    rate.model <- sapply(models, function(x) x$rate.model)
    log.lik <- sapply(models, function(x) x$log.lik)
    ## If 'gap.char' not present, then use the default for the estimate function.
    gap.char <- sapply(models, function(x) ifelse(test = is.null(x$gap.char), yes = "-", no = x$gap.char) )

    ## Get delta AIC across the models:
    rate.pars <- sapply(rate.model, function(x) switch(x, "correlated"=3, "gamma"=2, "single.rate"=1))
    extra.par <- sapply(Q.model, function(x) switch(x, "ER"=0, "DEL"=1))
    par.number <- rate.pars + extra.par
    aic.vector <- (-2 * log.lik) + (2 * par.number)
    delta.aic.vector <- aic.vector - min( aic.vector )
    aic.w.vector <- exp( -0.5 * delta.aic.vector ) / sum( exp( -0.5 * delta.aic.vector ) )
    
    ## Get the rates across the sequence positions for each of the models:
    seqrate <- list()
    for( i in 1:length(models) ){
        seqrate[[i]] <- getRates(x = models[[i]], Q.model = Q.model[i], gap.char = gap.char[i])
    }

    ## Produce the average rate per sequence position.
    ## This will need to be smart and deal with possible NA.
    ## Also, ER models are vectors and DEL models are matrices.
    ## First, duplicate ER models into matrices:
    for( i in 1:length(seqrate) ){
        if( Q.model[i] == "ER"){
            seqrate[[i]] <- cbind(seqrate[[i]], seqrate[[i]])
        }
    }

    ## Late test, but just to make sure:
    npos <- sapply(seqrate, nrow)
    if( length( unique(npos) ) != 1 ){
        stop( "Length of the sequences are not equal. Cannot model average! Different data!" )
    }
    ## Function here will return two elements, first is a rate and second is a state code.
    ave.res <- sapply(1:npos[1], function(x) make.average.rate(pos=x, seqrate=seqrate, aicw=aic.w.vector) )
    ave.rate.vec <- ave.res[1,]

    ## Put results together:
    ave.state.human <- sapply(ave.res[2,], function(x) switch(x, "good", "invar", "warning") )

    if( verbose ){
    ## Print a summary to the console.
        print( ftable( .=ave.state.human) )
    }

    ## Need to name the vector of AICw and return it:
    if( !is.null( names(models) ) ){
        if( length( names( aic.w.vector ) ) == length( names( models ) ) ){
            ## Checking before to not break in a dumb way.
            names( aic.w.vector ) <- names( models )
        }
    }
    
    out <- list( rates = ave.rate.vec, info = ave.state.human, AICw = aic.w.vector)
    class(out) <- append(class(out), "average.rates")
    return( out )
}

## Helping function:
make.average.rate <- function(pos, seqrate, aicw){
    ## Function gets the model average rate at sequence position 'pos' across models.
    pos.rates <- sapply(seqrate, function(x) x[pos,])
    aicw.mat <- rbind( aicw, aicw )        
    ## First get within model average of state and gap rates.
    pos.rates.row <- apply(pos.rates, 2, function(x) mean(x, na.rm = TRUE) )
    ## Now, if a particular model has NA or NaN for the position. Then we better drop it.
    ## Check if the model that we dropped has hight aicw. If not, then fine. If yes, then need to return
    ##       a warning message. Point is that it would be ok if the model failed to get a good estimate and
    ##       it is a poor model. It would not be cool if the model that failed at some of the positions is
    ##       among the best models.
    bad.rates <- is.na(pos.rates.row) | is.nan(pos.rates.row) | is.infinite(pos.rates.row)
    
    if( any( bad.rates ) ){
        ## If all rates are NA, then the positions is invariant.
        if( all( is.na( pos.rates.row ) ) ){
            return( c(NA, 2) )
        }        
        if( any(aicw[bad.rates] > 0.3) ){ ## using 30% as the cuoff for a warning.
            warning(paste0("Bad estimate for sequence position ", pos, ". Using only valid estimates."))
        }
        final.rate <- mean( pos.rates.row[!bad.rates] * aicw[!bad.rates] )
        return( c(final.rate, 3) )
    } else{
        ## All good.
        final.rate <- mean( pos.rates.row * aicw )
        return( c(final.rate, 1) )
    }
    
}
