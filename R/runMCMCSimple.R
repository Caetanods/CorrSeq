##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title A simple MCMC sampler for the model.
##' @param phy 
##' @param data 
##' @param k 
##' @param gen 
##' @param Q.max 
##' @param beta.max 
##' @return Write posterior samples to the directory.
##' @author daniel
##' @export
runMCMCSimple <- function(phy, data, k=4, gen, Q.max=100, beta.max="estimated"){
    ## Get the number of states in the data. Here we assume every site has the same number of states.
    nstates <- length(unique(data[,1]))

    ## Get the maximum value for beta:
    ## IMPORTANT: 'findMaxBeta' returns both the minimum and maximum for the beta parameter.
    stop( "Sorry Daniel, I am affraid you will need to fix this first." )
    if( beta.max == "estimated" ){
        beta.max <- findMaxBeta(ncats=k)
        print( paste("Max shape for Gamma is ", beta.max, sep="") )
    }
    if( !is.numeric(beta.max) ) stop( "beta.max needs to be 'estimated' or a number." )
    
    ## Start values:
    lik.curr <- NA
    while( is.na(lik.curr) ){
        ## Stupid wrap to make sure the random draw to start works.
        pi.0 <- runif(ncol(data)) ## Number of sample here will depend on the size of the data.
        pi.1 <- 1 - pi.0
        pi.curr <- rbind(pi.0, pi.1)
        Q.curr <- matrix(runif(1, min=0, max=Q.max), ncol=nstates, nrow=nstates)
        diag(Q.curr) <- 0
        diag(Q.curr) <- -rowSums(Q.curr)

        beta.curr <- runif(n=1, min=0, max=beta.max)

        lik.curr.list <- loglikGammaSimple(phy=phy, X=data, Q=Q.curr, pi=pi.curr, beta=beta.curr, k=k, it="random starting value. Sampling again.")
        lik.curr <- lik.curr.list[[1]]
        real.Q.curr <- lik.curr.list[[2]]
    }

    ## Open file connections. This follows scheme in ratematrix.
    ID <- paste0(sample(1:9, size=5, replace = TRUE), collapse="")
    file.name <- file.path( getwd(), paste0("output_", ID, ".mcmc") )
    file.name.mat <- file.path( getwd(), paste0("output_mat_", ID, ".mcmc") )
    mcmc.file <- file(file.name, open="a")
    transmat.file <- file(file.name.mat, open="a")
    print( paste0("Writing MCMC to: ", paste0("output_", ID, ".mcmc")) )

    ## Create the header for the files:
    Q.header <- vector(mode="character")
    for( i in 1:ncol(Q.curr) ){
        for( j in 1:nrow(Q.curr) ){
            Q.header <- c( Q.header, paste("Q_", i, j, "; ", sep="") )
        }
    }
    pi.header <- paste0("pi_", 1:ncol(pi.curr), "; ")
    cat(Q.header, pi.header, "beta; lik; accept \n", sep="", file=mcmc.file, append=TRUE)

    ## Header for the matrix file:
    mat.header <- vector(mode="character")
    for( z in 1:ncol(data) ){
        for( i in 1:nrow(Q.curr) ){
            for( j in 1:ncol(Q.curr) ){
                mat.header <- c( mat.header, paste("site", z, ".mat_", i, j, sep="") )
            }
        }
    }
    c.mat.header <- paste(mat.header, collapse="; ")
    cat(c.mat.header, "\n", sep="", file=transmat.file, append=TRUE)

    ## Start the MCMC:
    parnames <- c("Q","pi","beta")
    print( "Starting MCMC..." )
    
    for( i in 1:gen){
        ## Populate the proposal values.
        Q.prop <- Q.curr
        pi.prop <- pi.curr
        beta.prop <- beta.curr
        prop.ratio <- 1
        ## Choose which to update
        update <- sample(parnames, size = sample(1:length(parnames), size = 1))
        ## Make the moves
        if("Q" %in% update){
            Q.prop.list <- makePropQ.symm(Q=Q.curr, scale = 2)
            Q.prop <- Q.prop.list[[1]]
            prop.ratio <- prop.ratio * Q.prop.list[[2]]
        }
        if("pi" %in% update){
            pi.prop <- makePropPi(pi=pi.curr, width = 0.5)
        }
        if("beta" %in% update){
            beta.prop.vec <- multiplierProposal(x=beta.curr, a=2)
            beta.prop <- beta.prop.vec[1]
            prop.ratio <- prop.ratio * beta.prop.vec[2]
        }
        ## Accept and reject
        ## The likelihood for the proposal
        lik.prop.list <- loglikGammaSimple(phy=phy, X=data, Q=Q.prop, pi=pi.prop, beta=beta.prop, k=k, it=i)
        lik.ratio <- lik.prop.list[[1]] - lik.curr
        ## The prior ratio for the proposal. This will be log(1) or log(0), because all priors here are flat.
        prior.ratio <- log(min(uniformPrior(x=Q.prop[1,2], min=0, max=Q.max), uniformPrior(x=beta.prop, min=0, max=beta.max)))
        ## Need to take into account the proposal ratio.
        r <- lik.ratio + prior.ratio + log(prop.ratio)
        if(exp(r) > runif(1)){ ## Accept. Write all props.
            cat(paste(c(Q.prop), collapse="; "), "; ", paste(c(pi.prop[1,]), collapse="; "), "; ", beta.prop, "; ", lik.prop.list[[1]], "; 1 \n"
              , sep="", file=mcmc.file, append=TRUE)
            cat(paste(do.call(c, lik.prop.list[[2]]), collapse="; "), " \n", sep="", file=transmat.file, append=TRUE)            
            ## Update the current values
            Q.curr <- Q.prop
            pi.curr <- pi.prop
            beta.curr <- beta.prop
            lik.curr <- lik.prop.list[[1]]
            real.Q.curr <- lik.prop.list[[2]]
        } else{                ## Reject. Write all currents.
            cat(paste(c(Q.curr), collapse="; "), "; ", paste(c(pi.curr[1,]), collapse="; "), "; ", beta.curr, "; ", lik.curr, "; 0 \n"
              , sep="", file=mcmc.file, append=TRUE)
            cat(paste(do.call(c, real.Q.curr), collapse="; "), " \n", sep="", file=transmat.file, append=TRUE) 
            ## Do nothing.
        }
    }
    ## Close the files:
    close(mcmc.file)
    close(transmat.file)
    
    print( "Finished analyses!" )
    return( file.name )
}
