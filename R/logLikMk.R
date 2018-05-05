## function derived from "phytools::fitMk" by Liam Revell

make.data.tips <- function(X){
    ## Function to correct the format of the data to pass to the likelihood.
    states <- unique( X )
    X.mat <- sapply(states, function(y) as.numeric(X == y) )
    rownames(X.mat) <- names(X)
    colnames(X.mat) <- states
    return( X.mat )
}

make.data.tips.numeric <- function(X, Q){
    ## Function to correct the format of the data to pass to the likelihood.
    ## Same as the other. But assumes the elements are a sequence of numbers.
    ## states <- 1:length(unique(X))
    states <- 1:ncol(Q)
    X.mat <- sapply(states, function(y) as.numeric(X == y) )
    rownames(X.mat) <- names(X)
    colnames(X.mat) <- states
    return( X.mat )
}

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title The log-likelihood for the MK model.
##' @param phy A phylogenetic tree.
##' @param X A named vector.
##' @param Q A matrix with the transition matrix.
##' @param root.type One of "madfitz" or "equal".
##' @return The log-likelihood value.
##' @author daniel
##' @importFrom expm expm
##' @importFrom ape Ntip
##' @noRd
logLikMk <- function(phy, X, Q, root.type){
    ## Log likelihood for the model. Notice that some pre-processing for the tree is possible and likely necessary.
    
    ## phy is a phylo tree
    ## X is a matrix with number of rows equal to the number of tips and columns equal to the number of traits.
    ## Q is a matrix with the Q transition matrix.
    ## root.type is one of "madfitz" or "equal". Need to get a better way to check this later.
    N <- Ntip(phy)
    M <- phy$Nnode
    m <- ncol(X)
    ## Make a matrix to store the probabilities for the nodes. Note that the first rows are the data observed at the tips.
    ## Here X will need to have colnames. I will need to use the colnames of X to be able to do the 'rbind'.
    ## Also, I need to keep these states in the same order. Because the identity of the states will be given by the position.
    xx <- matrix(0, nrow = M, ncol = m )
    colnames(xx) <- colnames(X)
    liks <- rbind(X, xx)
    phy.prun <- reorder(phy, "pruningwise")
    comp <- vector(length=N+M, mode="numeric")
    parents <- unique(phy.prun$edge[,1])
    root <- min(parents)
    
    for(i in 1:length(parents)){
        anc <- parents[i]
        ii <- which(phy.prun$edge[,1] == parents[i])
        desc <- phy.prun$edge[ii,2]
        el <- phy.prun$edge.length[ii]
        v <- vector(length=length(desc), mode="list")
        for(j in 1:length(v)){
            ## Seems better to keep track of this in log space already. Need to check.
            v[[j]] <- expm(Q*el[j], method = "Ward77") %*% liks[desc[j],]
        }
        ## At the root ( when anc == root ) we have a vector of probabilities equal to the number of states.
        ## Right now we take a sum of these. But we can weight then by their probabilities if we use the "madfitz" model.
        if(anc==root){ ## At the root.
                  ## The behavior depends on the chosen model.
                  if(root.type == "equal"){
                      equal.pi <- rep(1, times = ncol(Q)) / ncol(Q)
                      vv <- Reduce('*',v)[,1] * equal.pi
                      comp[anc] <- sum(vv)
                  }
            if(root.type == "madfitz"){
                ## Implement the Maddison and Fitzjohn method.
                      vv.temp <- Reduce('*',v)[,1] ## The probability not yet scaled by the roots.
                      comp.anc <- sum(vv.temp)
                      liks.root <- vv.temp / comp.anc
                      root.p <- liks.root / sum(liks.root)
                      vv <- vv.temp * root.p ## Now we can scale the probability.
                      comp[anc] <- sum(vv)
                  }
              } else{
                  vv <- Reduce('*',v)[,1]
                  comp[anc] <- sum(vv) 
              }
        ## Here comp is the sum of the unnormalized probabilities, scalar, in normal space. This is what we need! We keep track of this quantity for the likelihood at the end.
        liks[anc,] <- vv / comp[anc] ## Vector of probabilities. Normalized. Normal space. This is the equivalent of the "contrast value", that is why it needs to be normalized to sum to 1.
    }

    ## This is the log-lik for the model.
    logL <- sum(log(comp[1:M+N])) ## Sum of the log( unnormalized probabilities ). [negative for convention]
    return( logL )
}
