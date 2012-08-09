############################################################################### 
## Weights for MCLUST
##
## Written by Thomas Brendan Murphy
##
#############################################################################

me.weighted <- function(modelName, data, z, weights = NULL, prior = NULL, 
control = emControl(), Vinv = NULL, warn = NULL, ...)
{
  N <- nrow(data)
  if (is.null(weights))
  {
    weights <- rep(1,N)
  }
  if (any(weights<0)|any(!is.finite(weights)))
  {
    stop("Weights must be positive and finite")
  }
  if (!is.vector(weights))
  {
    stop("Weights must be a vector")
  }
  if (max(weights)>1)
  {
    warning("Weights rescaled to have maximum equal to 1")
    weights <- weights/max(weights)
  }
  zw <- z*weights
  llold <- -Inf
  criterion <- TRUE
  iter <- 0
  while (criterion)
  {
    iter <- iter+1
    fit.m <- do.call("mstep",list(data=data,z=zw
      ,modelName=modelName,prior=prior
      ,control=control,Vinv=Vinv,warn=warn))
    fit.m$parameters$pro <- fit.m$parameters$pro/mean(weights)
    fit.e <- do.call("estep",c(list(data=data
      ,control=control,Vinv=Vinv,warn=warn),fit.m))
    zw <- fit.e$z*weights
    criterion <- criterion&(iter<control$itmax[1])
    ldens <- do.call("dens",c(list(data=data,logarithm=TRUE
      ,warn=warn),fit.m))
    ll <- sum(weights*ldens)
    criterion <- criterion&(ll-llold>control$tol[1])
    llold <- ll
  }
  fit <- fit.m
  fit$z <- fit.e$z
  fit$weights <- weights
  fit$loglik <- ll
  fit
}
