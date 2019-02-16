##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Function to compute the marginal estimate of root values.
##' @param phy 
##' @param data 
##' @param fit 
##' @return A list with the marginal root probabilities.
##' @author daniel
##' @importFrom expm expm
##' @importFrom ape Ntip
getRootMarginal <- function(phy, data, fit){
    ## Uses the same algorithm as Maddison and Fitzjohn to compute the marginal probability for the root of the tree.
    ## This will require computing the likelihoods for the nodes again.

    ## Get information from the fit object and the data:
    ## This chunk is general. Might migrate to a 'prepare.data' function.
    if( is.null( rownames(data) ) ) stop("data need to have rownames as the species names.")
    match.names <- all( rownames(data) %in% phy$tip.label ) & all( phy$tip.label %in% rownames(data) )
    if( !match.names ) stop("Secies names do not match between data and phylogeny!")

    ## Make data checks and get information from the matrix.
    nsites <- ncol(data)
    nstates <- apply(data, 2, function(x) length( unique(x) ) )
    names.data <- rownames(data)
    ## Translate states to numeric, but keep the labels safe.
    states.key <- lapply(1:nsites, function(x)
        cbind( unique(data[,x]), 1:length(unique(data[,x])) )
        )
    for(i in 1:nsites){
        for(j in 1:length(states.key[[i]][,1])){
            id <- data[,i] == states.key[[i]][,1][j]
            data[id,i] <- as.numeric( states.key[[i]][,2] )[j]
        }
    }
    ## Because R is mega dumb. Need to make sure that the table is really numeric.
    data <- apply(as.matrix(data), 2, as.numeric)
    rownames(data) <- names.data
    Xlist <- lapply(1:ncol(data), function(x) make.data.tips.numeric(setNames(data[,x], rownames(data))) )

    ## Get the Q matrix for each of the sites.
    Q <- fit$Q
    
    ## These are the parameters for the block below.
    ## phy, X=Xlist[[site]], Q=r*Q[[site]], root.type=root.type
    marg.root <- list() ## Return the marginal for the root.
    N <- Ntip(phy)
    M <- phy$Nnode
    phy.prun <- reorder(phy, "pruningwise")
    parents <- unique(phy.prun$edge[,1])
    root <- min(parents)
    for( w in 1:nsites ){
        if( is.character(Q[[w]]) ){ ## Invariant site
            marg.root[[w]] <- "invariant"
            next ## Skip the rest.
        }
        
        m <- ncol(Xlist[[w]])
        xx <- matrix(0, nrow = M, ncol = m )
        colnames(xx) <- colnames(Xlist[[w]])
        liks <- rbind(Xlist[[w]], xx)
        comp <- vector(length=N+M, mode="numeric")
        
        for(i in 1:length(parents)){
            anc <- parents[i]
            ii <- which(phy.prun$edge[,1] == parents[i])
            desc <- phy.prun$edge[ii,2]
            el <- phy.prun$edge.length[ii]
            v <- vector(length=length(desc), mode="list")
            for(j in 1:length(v)){
                ## Seems better to keep track of this in log space already. Need to check.
                v[[j]] <- expm(Q[[w]]*el[j], method = "Ward77") %*% liks[desc[j],]
            }
            if(anc==root){ ## At the root.
                ## This is the same as root.type == "madfitz" for the estimate function.
                vv.temp <- Reduce('*',v)[,1] ## The probability not yet scaled by the roots.
                comp.anc <- sum(vv.temp)
                liks.root <- vv.temp / comp.anc
                root.p <- liks.root / sum(liks.root)
                marg.root[[w]] <- root.p ## The marginal for the root.
                names(marg.root[[w]]) <- rownames(Q[[w]]) ## Set the names.
                vv <- vv.temp * root.p ## Now we can scale the probability.
                comp[anc] <- sum(vv)
            } else{
                vv <- Reduce('*',v)[,1]
                comp[anc] <- sum(vv) 
            }
            liks[anc,] <- vv / comp[anc]
        }
    }
    
    return( marg.root )
}

