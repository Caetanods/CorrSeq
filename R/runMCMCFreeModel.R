##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A simple MCMC sampler for the model.
##' @param phy 
##' @param data 
##' @param root.type 
##' @param gen 
##' @param unif.prior 
##' @return Write posterior samples to the directory.
##' @author daniel
##' @export
runMCMCFreeModel <- function(phy, data, root.type = "equal", gen, unif.prior=c(0.00001, 100), init = NULL){

    ## Assumes that all models are "ER" by now. Need to expand later.

    root.type <- match.arg(root.type, choices=c("madfitz","equal"), several.ok=FALSE)

    ## Create a list of data with length equal to the number of sites.
    nsites <- ncol(data)
    data.list <- lapply(1:nsites, function(x) make.data.tips(setNames(as.numeric(data[,x]), rownames(data))) )
    
    ## Get the number of states in the data. Here we assume every site has the same number of states.
    nstates <- length(unique(data[,1]))

    ## Create a wrap of the likelihood function.
    ## All the parameters are fixed using the data information.
    wrapLikMCMC <- function(x){
        ## x is a vector of parameters.
        ## length of x is equal to the number of sites.
        ## phy = phylogeny
        ## data = data matrix.
        ## nstates = number of states in each of the sites in the data.
        Q <- list()
        for( i in 1:nsites){
            Q[[i]] <- matrix(x[i], nrow=nstates, ncol=nstates)
            diag(Q[[i]]) <- -(colSums(Q[[i]]) - x[i])
        }
        lik.site <- sapply(1:nsites, function(x) logLikMk(phy=phy, X=data.list[[x]], Q=Q[[x]], root.type))
        lik <- sum(lik.site)
        return( lik ) ## Here I need the normal log-lik for the model.
    }
    
    ## Start values:
    if( is.null(init) ){
        Q.curr <- runif(nsites, min=unif.prior[1], max=unif.prior[2])
    } else{
        if( any(init < unif.prior[1]) | any(init > unif.prior[2]) ) stop("Provided 'init' value is outside search bounds.")
        if( !length(init) == nsites ) stop(" Length of 'init' needs to be the same as the length of the sequence." )
        Q.curr <- init
    }

    ## Open file connections. This follows scheme in ratematrix.
    ID <- paste0(sample(1:9, size=5, replace = TRUE), collapse="")
    file.name <- file.path( getwd(), paste0("output_", ID, ".mcmc") )
    mcmc.file <- file(file.name, open="a")
    print( paste0("Writing MCMC to: ", paste0("output_", ID, ".mcmc")) )

    ## Create the header for the files:
    Q.header <- paste(paste("site_", 1:nsites, sep=""), collapse="; ")
    cat(Q.header, "; lik; accept \n", sep="", file=mcmc.file, append=TRUE)

    ## Start the MCMC:
    print( "Starting MCMC..." )

    ## Which number of sites to update each time?
    min.update <- 1
    max.update <- ceiling(nsites/2)

    ## Inital likelihood and parameters value:
    lik.curr <- wrapLikMCMC(x = Q.curr)

    ## Write first generation as the starting state:
    cat(paste(Q.curr, collapse="; "), "; ", lik.curr, "; 1 \n"
      , sep="", file=mcmc.file, append=TRUE)
    
    for( i in 1:(gen-1)){ ## One gen less. First is the starting point.
        ## Here we will update a different number of sites each time.
        ## This can be controlled in the function.
        n.update <- sample(min.update:max.update, size = 1)
        which.update <- sample(1:nsites, size = n.update)

        ## Update step:
        mult.out <- sapply(Q.curr[which.update], function(i) multiplierProposal(x=i, a=2) )
        Q.prop <- Q.curr ## Generate the vector.
        Q.prop[which.update] <- mult.out[1,] ## First line has the parameters.
        prop.ratio <- sum(mult.out[2,]) ## Second line has the prop.ratio.
        ## Compute the likelihood.
        lik.prop <- wrapLikMCMC(x = Q.prop)

        ## Accept and reject
        lik.ratio <- lik.prop - lik.curr
        prior.ratio <- uniformPrior(x=Q.prop, min=unif.prior[1], max=unif.prior[2])
        r <- lik.ratio + prop.ratio + prior.ratio
        if(exp(r) > runif(1)){ ## Accept. Write all props.
            cat(paste(Q.prop, collapse="; "), "; ", lik.prop, "; 1 \n"
              , sep="", file=mcmc.file, append=TRUE)
            ## Update the current values
            Q.curr <- Q.prop
            lik.curr <- lik.prop
        } else{                ## Reject. Write all currents.
            cat(paste(Q.curr, collapse="; "), "; ", lik.curr, "; 0 \n"
              , sep="", file=mcmc.file, append=TRUE)
            ## Do nothing.
        }
    }
    ## Close the files:
    close(mcmc.file)
    
    print( "Finished analyses!" )
    return( file.name )
}
