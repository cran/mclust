############################################################################### 
## Weights for MCLUST
##
## Written by Thomas Brendan Murphy
## Bugs fix by Luca Scrucca
#############################################################################

me.weighted <- function(data, modelName, z, weights = NULL, prior = NULL, 
                        control = emControl(), Vinv = NULL, warn = NULL, ...)
{
  data <- as.matrix(data)
  nobs <- nrow(data)
  modelName <- switch(EXPR = modelName,
                      "X" = "E",
                      "XII" = "EII",
                      "XXI" = "EEI",
                      "XXX" = "EEE",
                      modelName)
  if(is.null(warn)) warn <- mclust.options("warn")
  if(is.null(weights))
    { weights <- rep(1,nobs) }
  if(any(weights < 0)| any(!is.finite(weights)))
    { stop("Weights must be positive and finite") }
  if(!is.vector(weights))
    { stop("Weights must be a vector") }
  if(max(weights) > 1)
    { if(warn)
        warning("Weights rescaled to have maximum equal to 1")
      weights <- weights/max(weights)
  }
  zw <- z*weights
  llold <- ll <- -Inf
  eps <- .Machine$double.eps
  criterion <- TRUE
  iter <- 0
  while(!is.na(criterion) & criterion)
  {
    iter <- iter+1
    fit.m <- do.call("mstep", list(data = data, 
                                   z = zw,
                                   modelName = modelName, 
                                   prior = prior,
                                   control = control, 
                                   Vinv = Vinv, 
                                   warn = warn))
    fit.m$parameters$pro <- fit.m$parameters$pro/mean(weights)
    fit.e <- do.call("estep", c(list(data = data,
                                     control = control, 
                                     Vinv = Vinv, 
                                     warn = warn),
                                fit.m))
    if(is.na(fit.e$loglik)) 
      { criterion <- FALSE; next() }
    zw <- pmax(fit.e$z*weights, eps)
    ldens <- do.call("dens", c(list(data = data, 
                                    logarithm = TRUE, 
                                    warn = warn), 
                               fit.m))
    ll <- sum(weights*ldens)
    criterion <- criterion & (iter < control$itmax[1])
    criterion <- criterion & ((ll-llold)/(1+abs(ll)) > control$tol[1])
    llold <- ll
  }
  fit <- fit.m
  fit$z <- fit.e$z
  fit$weights <- weights
  fit$loglik <- ll/mean(weights)
  npar <- nMclustParams(modelName = fit$modelName, d = fit$d, G = fit$G)
  fit$bic <- 2*fit$loglik - npar*log(fit$n)  
  return(fit)
}
