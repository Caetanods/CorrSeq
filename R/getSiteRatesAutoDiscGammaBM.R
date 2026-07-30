getMarginalBMSimpleGamma <- function(lik_fn, solution, k, n.cores){
  ## Function to estimate the marginal rate of the bounded BM model for each sequence position.
  ## This is to be used with the maximum likelihood value for the parameters and NOT part of the MLE estimation.
  gamma.rates <- discreteGamma(shape = solution[2], ncats = k)
  gamma.rates <- gamma.rates[!round(gamma.rates, digits=20) == 0] ## The strict test with 0.0 will never trigger.
  gamma_bm_rates <- gamma.rates * solution[1] # A vector with the scaled rates with length k
  gamma.lik <- lapply(1:length(lik_fn), function(site) sapply(gamma_bm_rates, function(r) lik_fn[[site]](sigma = r)))
  # rel.lik <- lapply(1:length(lik_fn), function(x) exp(gamma.lik[[x]]) / sum( exp(gamma.lik[[x]]) ) )
  ## Doing the exponential AFTER taking the proportion below.
  rel.lik <- lapply(1:length(lik_fn), function(x) exp( gamma.lik[[x]] - logSumExp( gamma.lik[[x]] ) ))
  marginal.bm <- vector(mode = "numeric", length = length(lik_fn)) # Making sure the output is a vector
  for( i in 1:length(lik_fn) ){
    marginal.bm[i] <- sum(rel.lik[[i]] * gamma.rates) * solution[1]
  }
  ## Returning only the marginal estimate of the rates across the sites.
  return( marginal.bm )
}

##' @importFrom parallel mclapply
getMarginalBMAutoDiscGamma <- function(lik_fn, solution, M, k, nsites, n.cores){
  ## Computes the marginal estimate for the rates under a BM model.
  gamma.rates <- discreteGamma(shape = solution[2], ncats = k)
  eff_cat <- ncol(M)
  ## Transition probabilities of zero on this matrix can cause problems on log space.
  ## Here trying to just set them to a very small number.
  M[ sapply(c(M), function(x) isTRUE( all.equal(x, 0.0)) ) ] <- .Machine$double.eps
  gamma.rates <- gamma.rates[ (k+1-eff_cat):k ]
  gamma_bm_rates <- gamma.rates * solution[1] # A vector with the scaled rates with length k
  ## This computes the likelihood for each site given all the rate categories.
  gamma.lik <- mclapply(1:nsites, function(site) sapply(gamma_bm_rates, function(r) lik_fn[[site]](sigma = r) ), mc.cores = n.cores)

  ## Need special code for nsites < 5.
  if( nsites > 4 ){
    ## This will work for 5 or more sites.
    
    ## Make matrix to store the marginals for each rate at each site:
    marg <- matrix(NA, nrow = eff_cat, ncol = nsites)
    
    ## Marginal for the last site.
    for( cat in 1:eff_cat ){
      N_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
      lik_unit <- N_unit ## Start the loop.
      for( site in (nsites-1):2){ ## Next to last up to the second site.
        ## The 'getIntUnit" function already deals with the probs of 0.0
        lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=site
                                                             , left_unit=lik_unit))
      }
      marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for the second to last site. [Will be the same for all middle sites.]
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites-1]][cat] + lik_unit[cat])
      for( site in (nsites-2):2){ ## Need to loop over all the rest.
        lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=site
                                                             , left_unit=lik_unit))
      }
      marg[cat, nsites-1] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## A loop from site nsites-2 to site 3.
    for( marg.site in (nsites-2):3 ){ ## The big loop for the marginal.
      for( cat in 1:eff_cat ){
        lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
        for( site in (nsites-1):marg.site+1){ ## The sites to the right of the focus site.
          lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                               , n=site, left_unit=lik_unit))
        }
        ## Compute the marginal for the focus site:
        lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[marg.site]][cat] + lik_unit[cat])
        for( site in (marg.site-1):2){ ## The sites to the left of the focus site.
          lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                               , n=site, left_unit=lik_unit))
        }
        marg[cat, marg.site] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
      }
    }
    
    ## Marginal for site 2.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      for( site in (nsites-1):3){ ## The sites to the right of the focus site.
        lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                             , n=site, left_unit=lik_unit))
      }
      ## Compute the marginal for the focus site:
      lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
      marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 1.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      for( site in (nsites-1):2){ ## The sites to the right of the focus site.
        lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i
                                                             , n=site, left_unit=lik_unit))
      }
      ## Compute the marginal for the focus site:
      marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
    }
  }
  if( nsites == 4){
    ## Special case with 4 sites in the data.
    
    ## Make matrix to store the marginals for each rate at each site:
    marg <- matrix(NA, nrow = eff_cat, ncol = nsites)
    
    ## Marginal for the last site.
    for( cat in 1:eff_cat ){
      N_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
      lik_unit <- N_unit ## Start the loop.
      ## For site 2:
      lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                           , left_unit=lik_unit))
      ## For site 1:
      marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 3
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      ## Compute the marginal for the focus site (site 3):
      lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[3]][cat] + lik_unit[cat])
      ## The sites to the left of the focus site (1 and 2):
      lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                           , left_unit=lik_unit))
      marg[cat, 3] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 2 (i.e., second to last site).
    ## In this case we only have a single middle site.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
      marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 1.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      ## The sites to the right of the focus site (i.e., site 2).
      lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                           , left_unit=lik_unit))
      ## Compute the marginal for the focus site:
      marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
    }
  }
  if( nsites == 3){
    ## Special case with 3 sites in the data.
    
    ## Make matrix to store the marginals for each rate at each site:
    marg <- matrix(NA, nrow = eff_cat, ncol = nsites)
    
    ## Marginal for the last site.
    for( cat in 1:eff_cat ){
      N_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
      lik_unit <- N_unit ## Start the loop.
      ## For site 2:
      lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                           , left_unit=lik_unit))
      ## For site 1:
      marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 2 (i.e., second to last site).
    ## In this case we only have a single middle site.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[2]][cat] + lik_unit[cat])
      marg[cat, 2] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 1.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      ## The sites to the right of the focus site (i.e., site 2).
      lik_unit <- sapply(1:eff_cat, function(i) getIntUnit(gamma.lik=gamma.lik, M=M, i=i, n=2
                                                           , left_unit=lik_unit))
      ## Compute the marginal for the focus site:
      marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
    }
  }
  if( nsites == 2){
    ## Special case with 2 sites in the data.
    
    ## Make matrix to store the marginals for each rate at each site:
    marg <- matrix(NA, nrow = k, ncol = nsites)
    
    ## Marginal for the last site.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) log(M[i,cat]) + gamma.lik[[nsites]][cat] )
      ## For site 1:
      marg[cat, nsites] <- getFirstUnit(gamma.lik=gamma.lik, k=k, second_unit=lik_unit)
    }
    
    ## Marginal for site 1.
    for( cat in 1:eff_cat ){
      lik_unit <- sapply(1:eff_cat, function(i) getLastUnit(gamma.lik=gamma.lik, M=M, i=i, n=nsites))
      ## Compute the marginal for the focus site:
      marg[cat, 1] <- log(1/eff_cat) + gamma.lik[[1]][cat] + lik_unit[cat]
    }
  }
  
  ## To return the rates we need to compute a weighted average of the Gamma rates following the proportion from the marginal likelihood for each of the hidden rates.
  ## This is equal to exp( log( exp(x) / sum( exp(x) ) ) ), but we cannot compute the exponential of the logLikelihood because of underflow issues. So below we compute the exponential only after the proportion is computed.
  rel.marg.lik <- apply( marg, 2, function(x) exp(x - logSumExp(x)) )
  ## The marginal computation below is with the BM rate in the natural scale.
  ## For loop is just to make the output in a vector.
  real.BM <- vector(mode = "numeric", length = ncol(rel.marg.lik))
  for( i in 1:ncol(rel.marg.lik) ){
    real.BM[i] <- sum(rel.marg.lik[,i] * gamma.rates) * solution[1]
  }
  
  ## Return a vector of marginal reconstructions for the bounded BM rate.
  return( real.BM )
}
