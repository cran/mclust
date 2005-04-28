##
## Port of mclust (2002 version) to R
## Original program by Chris Fraley and Adrian Raftery
## Port by Ron Wehrens
##
#cat("\nWarning: this is the 2002 version of mclust.",
#    "\n\nIt is not backwards compatible with the previous version,",
#    "\nso old scripts will no longer work. The previous version",
#    "\nof mclust is still available as mclust1998.",
#    "\n\nSince mclust1998 is not actively supported any more,",
#    "\nplease change to the new version.\n\n") 
.First.lib <- function(lib, pkg) {
  library.dynam("mclust", pkg, lib)
}


"[.mclustDAtest" <- function(x, i, j, drop = FALSE)
{
  clx <- oldClass(x)
  oldClass(x) <- NULL
  NextMethod("[")
}


".Mclust" <- 
  list("eps" = .Machine$double.eps, ## 2.2204460492503101e-16,
       "tol" = c(1.0000000000000001e-05, 1.0000000000000001e-05),
       "itmax" = c(Inf, Inf),
       "equalPro" = FALSE,
       "warnSingular" = TRUE,
       "emModelNames" = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE",
         "EEV", "VEV", "VVV"),
       "hcModelName" = c("E", "VVV"),
       "symbols" = c(17, 0, 10, 4, 11, 18, 6, 7, 3, 16, 2, 12, 8, 15,
         1, 9, 14, 13, 5))

"EMclust" <- function(data, G, emModelNames, hcPairs, subset, eps,
                      tol, itmax, equalPro, warnSingular = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
## 
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(oneD) {
    data <- as.vector(data)
    n <- length(data)
    p <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
  }
  if(missing(emModelNames)) {
    if(p == 1) {
      emModelNames <- c("E", "V")
    }
    else {
      emModelNames <- .Mclust$emModelNames
    }
  }
  if(p == 1 && any(nchar(emModelNames) > 1)) {
    Emodel <- any(sapply(emModelNames, function(x)
                         charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(emModelNames, function(x)
                         charmatch("V", x, nomatch = 0)[1]) == 1)
    emModelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  m <- length(emModelNames)
  if(missing(G)) {
    G <- 1:9
  }
  else {
    G <- sort(G)
  }
  if(any(G) <= 0)
    stop("G must be positive")
  l <- length(G)
  Glabels <- as.character(G)
  BIC <- matrix(0, nrow = l, ncol = m, dimnames = list(Glabels, emModelNames))
  if(G[1] == 1) {
    for(mdl in emModelNames) {
      hood <- mvn(modelName = mdl, data = data)$loglik
      BIC[1, mdl] <- bic(modelName = mdl, loglik = hood,
                         n = n, d = p, G = 1, equalPro = equalPro)
    }
    if(l == 1) {
      return(structure(BIC, equalPro = equalPro, args = 
                       as.list(match.call())[-1], class = "EMclust"))
    }
    G <- G[-1]
    Glabels <- Glabels[-1]
  }
  if(missing(subset)) {
    subset <- NULL
    ## #####################################################
    ## all data in initial hierarchical clustering phase
    ## #####################################################
    if(missing(hcPairs)) {
      if(p != 1) {
        hcPairs <- hc(modelName = .Mclust$hcModelName[
                        2], data = data)
      }
      else {
        hcPairs <- hc(modelName = .Mclust$hcModelName[
                        1], data = data)
      }
    }
    clss <- hclass(hcPairs, G)
    for(i in seq(along = G)) {
      z <- unmap(clss[, Glabels[i]])
      for(modelName in emModelNames) {
        hood <- me(modelName = modelName, data = data,
                   z = z, eps = eps, tol = tol, itmax = 
                   itmax, equalPro = equalPro, 
                   warnSingular = warnSingular)$loglik
        BIC[Glabels[i], modelName] <-
          bic(modelName = 
              modelName, loglik = hood, n = n, d = p, G = G[i],
              equalPro = equalPro)
      }
    }
  }
  else {
    ## ####################################################
    ## sample for the initial hierarchical clustering phase
    ## ####################################################
    if(is.logical(subset)) subset <- (1:n)[subset]
    if(missing(hcPairs)) {
      if(p != 1) {
        hcPairs <- hc(modelName = .Mclust$hcModelName[
                        2], data = data[subset,  ])
      }
      else {
        hcPairs <- hc(modelName = .Mclust$hcModelName[
                        1], data = data[subset,  ])
      }
    }
    clss <- hclass(hcPairs, G)
    if(length(tol) > 1)
      tol <- tol[2]
    if(length(itmax) > 1)
      itmax <- itmax[2]
    for(i in seq(along = G)) {
      z <- unmap(clss[, Glabels[i]])
      dimnames(z) <- list(as.character(subset), NULL)
      for(modelName in emModelNames) {
        ms <- mstep(modelName = modelName,
                    data = data[subset,  ], z = z, eps = eps, tol = tol,
                    itmax = itmax, equalPro = equalPro,
                    warnSingular = warnSingular)
        hood <- do.call("em", c(list(data = data, eps = eps,
                                     tol = tol, itmax = itmax, equalPro =
                                     equalPro, warnSingular =
                                     warnSingular), ms))$loglik 
        BIC[Glabels[i], modelName] <-
          bic(modelName = modelName, loglik = hood, n = n, d = p,
              G = G[i], equalPro = equalPro)
      }
    }
  }
  ##
  ## separating hc from its attributes gets around what seems to be a bug
  ##
  attrHC <- attributes(hcPairs)
  attributes(hcPairs) <- NULL
  hcPairs <- matrix(hcPairs, nrow = 2, ncol = length(hcPairs)/2)
  structure(BIC, subset = subset, eps = eps, tol = tol, itmax = itmax,
            equalPro = equalPro, warnSingular = warnSingular, hcPairs = 
            hcPairs, attrHC = attrHC, args = as.list(match.call())[-1],
            class = "EMclust")
}

"EMclustN" <- function(data, G, emModelNames, noise, hcPairs, eps,
                       tol, itmax, equalPro, warnSingular = FALSE, Vinv,
                       ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  ##
  ## noise  is a logical vector in which TRUE indicates noise
  ##
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(missing(tol))
    tol <- .Mclust$tol
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(oneD) {
    data <- as.vector(data)
    n <- length(data)
    p <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
  }
  if(missing(emModelNames)) {
    if(p == 1) {
      emModelNames <- c("E", "V")
    }
    else {
      emModelNames <- .Mclust$emModelNames
    }
  }
  m <- length(emModelNames)
  if(missing(G)) {
    G <- 0:9
  }
  else {
    G <- sort(G)
  }
  if(any(G) < 0)
    stop("G must be non negative")
  l <- length(G)
  Glabels <- as.character(G)
  BIC <- matrix(0, nrow = l, ncol = m, dimnames = list(Glabels, 
                                         emModelNames))
  if(missing(Vinv) || Vinv <= 0)
    Vinv <- hypvol(data, reciprocal = TRUE)
  if(!is.logical(noise))
    noise <- as.logical(match(1:n, noise, nomatch = 0))
  if(!G[1]) {
    hood <- n * logb(Vinv)
    BIC[1,  ] <- 2 * hood - logb(n)
    if(l == 1) {
      return(structure(BIC, equalPro = equalPro, noise = 
                       noise, Vinv = Vinv, args = as.list(match.call(
                                             ))[-1], class = "EMclustN"))
    }
    G <- G[-1]
    Glabels <- Glabels[-1]
  }
  if(missing(hcPairs)) {
    if(p != 1) {
      hcPairs <- hc(modelName = .Mclust$hcModelName[2], data
                    = data[!noise,  ])
    }
    else {
      hcPairs <- hc(modelName = .Mclust$hcModelName[1], data
                    = data[!noise,  ])
    }
  }
  clss <- hclass(hcPairs, G)
  z <- matrix(0, n, max(G) + 1)
  for(i in seq(along = G)) {
    z[!noise, 1:G[i]] <- unmap(clss[, Glabels[i]])
    z[noise, 1:G[i]] <- 0
    G1 <- G[i] + 1
    z[!noise, G1] <- 0
    z[noise, G1] <- 1
    for(modelName in emModelNames) {
      hood <- me(modelName = modelName, data = data, z = z[, 1:G1],
                 eps = eps, tol = tol, itmax = itmax,
                 equalPro = equalPro, noise = TRUE, Vinv = Vinv,
                 warnSingular = warnSingular)$loglik
      BIC[Glabels[i], modelName] <-
        bic(modelName = modelName, loglik = hood, n = n, d = p, G =
            G[i], equalPro = equalPro, noise = TRUE)
    }
  }
  ##
  ## separating hc from its attributes gets around what seems to be a bug
  ##
  attrHC <- attributes(hcPairs)
  attributes(hcPairs) <- NULL
  hcPairs <- matrix(hcPairs, nrow = 2, ncol = length(hcPairs)/2)
  structure(BIC, eps = eps, tol = tol, itmax = itmax, equalPro = equalPro,
            warnSingular = warnSingular, noise = noise, Vinv = Vinv, 
            hcPairs = hcPairs, attrHC = attrHC,
            args = as.list(match.call())[-1], class = "EMclustN")
}

"Mclust" <- function(data, minG = 1, maxG = 9)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  emModelNames <- c("EII", "VII", "EEI", "VVI", "EEE", "VVV")
  G <- minG:maxG
  Bic <- EMclust(data, G = G, emModelNames = emModelNames)
  Sumry <- summary(Bic, data)
  if(!(length(G) == 1)) {
    bestG <- length(unique(Sumry$cl))
    if(bestG == max(G))
      warning("optimal number of clusters occurs at max choice"
              )
    else if(bestG == min(G))
      warning("optimal number of clusters occurs at min choice"
              )
  }
  attr(Bic, "hcPairs") <- attr(Bic, "attrHC") <- NULL
  attr(Bic, "args") <- attr(Bic, "equal") <- attr(Bic, "class") <- NULL
  Sumry$cholsigma <- Sumry$cholSigma <- NULL
  Sumry$Vinv <- Sumry$options <- NULL
  Sumry$bic <- Sumry$bic[1]
  attr(Sumry, "class") <- NULL
  structure(c(list(BIC = Bic), Sumry), class = "Mclust")
}

"[.EMclust" <- function(x, i, j, drop = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
## 
  clx <- class(x)
  attrx <- attributes(x)[c("Vinv", "args", "class", "equal", "fuzzy",
                           "hc", "noise")]
  class(x) <- NULL
  x <- NextMethod("[")
  do.call("structure", c(list(.Data = x), attrx))
}

"[.EMclustN" <- function(x, i, j, drop = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  clx <- class(x)
  attrx <- attributes(x)[c("Vinv", "args", "class", "equal", "fuzzy",
                           "hc", "noise")]
  class(x) <- NULL
  x <- NextMethod("[")
  do.call("structure", c(list(.Data = x), attrx))
}

"bic" <- function(modelName, loglik, n, d, G, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  modelName <- switch(modelName,
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  funcName <- paste("bic", modelName, sep = "")
  do.call(funcName, list(loglik = loglik, n = n, d = d, G = G, ...))
}

"bicE" <- function(loglik, n, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G + 1
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicEEE" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    s <- (d * (d + 1))/2
    nparams <- G * d + s
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicEEI" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G * d + d + 1
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicEEV" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    s <- (d * (d - 1))/2
    nparams <- G * (d + s) + (d - 1) + 1
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicEII" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G * d + 1
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicEMtrain" <- function(data, labels, modelNames)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  z <- unmap(as.numeric(labels))
  G <- ncol(z)
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(oneD || length(dimData) != 2) {
    if(missing(modelNames))
      modelNames <- c("E", "V")
    if(any(!match(modelNames, c("E", "V"), nomatch = 0)))
      stop("modelNames E or V for one-dimensional data")
  }
  else {
    if(missing(modelNames))
      modelNames <- .Mclust$emModelNames
  }
  BIC <- rep(NA, length(modelNames))
  names(BIC) <- modelNames
  for(m in modelNames) {
    mStep <- mstep(modelName = m, data = data, z = z, warnSingular = FALSE)
    eStep <- do.call("estep", c(mStep, list(data = data, 
                                            warnSingular = FALSE)))
    if(is.null(attr(eStep, "warn")))
      BIC[m] <- do.call("bic", eStep)
  }
  BIC
}

"bicEVI" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##  
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
  }
  else {
    nparams <- 2 * (d * G) + 1
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicV" <- function(loglik, n, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G * 2
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicVEI" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G * d + d + G
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicVEV" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    ## d*d - d(d+1)/2 for orientation, 1 for volume
    s <- (d * (d - 1))/2 + 1
    nparams <- G * (d + s) + (d - 1)
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicVII" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G * (d + 1)
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicVVI" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    nparams <- G * d + (d + 1) * G
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"bicVVV" <- function(loglik, n, d, G, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(G == 0) {
    ## one cluster case
    if(!noise) stop("undefined model")
    nparams <- 1
  }
  else {
    s <- (d * (d + 1))/2
    nparams <- G * (d + s)
    if(!equalPro)
      nparams <- nparams + (G - 1)
    if(noise)
      nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}

"cdens" <- function(modelName, data, mu, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## ... sigmasq or sigma, eps
  modName <- switch(EXPR=modelName,
                    X = "E",
                    XII = "EII",
                    XXI = "EEI",
                    XXX = "EEE",
                    modelName)
  funcName <- paste("cdens", modName, sep = "")
  out <- do.call(funcName, list(data = data, mu = mu, ...))
  modName <- switch(EXPR = modelName,
                    X = "X",
                    XII = "XII",
                    XXI = "XXI",
                    XXX = "XXX",
                    modName)
  attr(out, "modelName") <- modName
  out
}

"cdensE" <- function(data, mu, sigmasq, eps, warnSingular,
                     logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  G <- length(mu)
  if(all(is.na(c(mu, sigmasq)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "E", warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(sigmasq <= eps) {
    warn <- "sigma-squared falls below threshold"
    warning("sigma-squared falls below threshold")
    return(structure(matrix(NA, n, G), modelName = "E", warn = warn))
  }
  temp <- .Fortran("den1e",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.integer(n),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[6:7]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    cden[] <- NA
  }
  structure(cden, modelName = "E", warn = warn)
}

"cdensEEE" <- function(data, mu, eps, warnSingular, logarithm=FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) > 2)
    stop("data must be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  cholSigma <- list(...)$cholSigma
  if(is.null(cholSigma)) {
    if(!is.null(decomp <- list(...)$decomp)) {
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      sig <- qr.R(qr(O * sqrt(scale * shape)))
      cholIND <- "U"
    }
    else if(!is.null(Sigma <- list(...)$Sigma)) {
      sig <- Sigma
      cholIND <- "N"
    }
    else if(!missing(sigma)) {
      sig <- sigma
      cholIND <- "N"
    }
    else stop("invalid specification for sigma")
  }
  else {
    sig <- cholSigma
    cholIND <- "U"
  }
  if(any(is.na(c(mu, sig)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "EEE", warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("deneee",
                   as.integer(if (cholIND=="N") 0 else 1),
                   as.double(data),
                   as.double(mu),
                   as.double(sig),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[8:10]
  lapackCholInfo <- temp[[1]][1]
  eps <- temp[[2]]
  cden <- matrix(if(logarithm) temp[[3]] else exp(temp[[3]]), n, G)
  warn <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      warn <- "sigma is not positive definite"
      if(warnSingular)
        warning("sigma is not positive definite")
    }
    else {
      warn <- "input error for LAPACK DPOTRF"
      warning("input error for LAPACK DPOTRF")
    }
    cden[] <- NA
  }
  else if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "EEE", warn = warn)
}

"cdensEEI" <- function(data, mu, decomp, eps, warnSingular,
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "EEI", warn = warn))
  }
  temp <- .Fortran("deneei",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[8:9]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "EEI", warn = warn)
}

"cdensEEV" <- function(data, mu, decomp, eps, warnSingular,
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  G <- ncol(mu)
  mu <- as.matrix(mu)
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "EEV", warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("deneev",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(decomp$orientation),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p),
                   double(p),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[11:12]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "EEV", warn = warn)
}

"cdensEII" <- function(data, mu, sigmasq, eps, warnSingular,
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(sigmasq)) {
    sigmasq <- list(...)$decomp$scale
  }
  if(any(is.na(c(mu, sigmasq)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "EII", warn = 
                     warn))
  }
  if(list(...)$decomp$scale <= eps) {
    warn <- "sigma-squared falls below threshold"
    warning("sigma-squared falls below threshold")
    return(structure(matrix(NA, n, G), modelName = "EII", warn = 
                     warn))
  }
  temp <- .Fortran("deneii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[7:8]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    cden[] - NA
  }
  structure(cden, modelName = "EII", warn = warn)
}

"cdensEVI" <- function(data, mu, decomp, eps, warnSingular,
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "EVI", warn = 
                     warn))
  }
  temp <- .Fortran("denevi",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[8:9]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "EVI", warn = warn)
}

"cdensV" <- function(data, mu, sigmasq, eps, warnSingular, 
                     logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  G <- length(mu)
  if(any(is.na(c(mu, sigmasq)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "V", warn = warn
                     ))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(any(sigmasq <= eps)) {
    warn <- "sigma-squared falls below threshold"
    warning("sigma-squared falls below threshold")
    return(structure(matrix(NA, n, G), modelName = "V", warn = warn
                     ))
  }
  temp <- .Fortran("den1v",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.integer(n),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[6:7]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    cden[] <- NA
  }
  structure(cden, modelName = "V", warn = warn)
}

"cdensVEI" <- function(data, mu, decomp, eps, warnSingular, 
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "VEI", warn = 
                     warn))
  }
  temp <- .Fortran("denvei",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[8:9]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "VEI", warn = warn)
}

"cdensVEV" <- function(data, mu, decomp, eps, warnSingular, 
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "VEV", warn = 
                     warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("denvev",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(decomp$orientation),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p),
                   double(p),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[11:12]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "VEV", warn = warn)
}

"cdensVII" <- function(data, mu, sigmasq, eps, warnSingular,
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(sigmasq)) {
    sigmasq <- list(...)$decomp$scale
  }
  if(any(is.na(c(mu, sigmasq)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "VII", warn = 
                     warn))
  }
  if(any(list(...)$decomp$scale <= eps)) {
    warn <- "sigma-squared falls below threshold"
    warning("sigma-squared falls below threshold")
    return(structure(matrix(NA, n, G), modelName = "VII", warn = 
                     warn))
  }
  temp <- .Fortran("denvii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[7:8]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    cden[] <- NA
  }
  structure(cden, modelName = "VII", warn = warn)
}

"cdensVVI" <- function(data, mu, decomp, eps, warnSingular, 
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "VVI", warn = 
                     warn))
  }
  temp <- .Fortran("denvvi",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[8:9]
  eps <- temp[[1]]
  cden <- matrix(if(logarithm) temp[[2]] else exp(temp[[2]]), n, G)
  warn <- NULL
  if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "VVI", warn = warn)
}

"cdensVVV" <- function(data, mu, eps, warnSingular, 
                       logarithm = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  cholsigma <- list(...)$cholsigma
  if(is.null(cholsigma)) {
    if(!is.null(sigma <- list(...)$sigma)) {
      sig <- sigma
      cholIND <- "N"
    }
    else if(!is.null(decomp <- list(...)$decomp)) {
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      sig <- array(0, c(p, p, G))
      shape <- sqrt(sweep(shape, MARGIN = 2, STATS = scale,
                          FUN = "*"))
      for(k in 1:G)
        sig[,  , k] <- qr.R(qr(O[,  , k] * shape))
      cholIND <- "U"
    }
    else stop("sigma improperly specified")
  }
  else {
    sig <- cholsigma
    cholIND <- "U"
  }
  if(any(is.na(c(mu, sig)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, G), modelName = "VVV", warn = 
                     warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("denvvv",
                   as.integer(if (cholIND == "N") 0 else 1),
                   as.double(data),
                   as.double(mu),
                   as.double(sig),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p),
                   as.double(eps),
                   double(n * G),
                   PACKAGE="mclust")[8:10]
  lapackCholInfo <- temp[[1]][1]
  eps <- temp[[2]]
  cden <- matrix(if(logarithm) temp[[3]] else exp(temp[[3]]), n, G)
  warn <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      warn <- "sigma is not positive definite"
      warning("sigma is not positive definite")
    }
    else {
      warn <- "input error for LAPACK DPOTRF"
      warning("input error for LAPACK DPOTRF")
    }
    cden[] <- NA
  }
  else if(is.infinite(eps) || eps == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    cden[] <- NA
  }
  structure(cden, modelName = "VVV", warn = warn)
}

### Changed clPairs to use the standard R pairs functionality
### Added color argument, because I like that...
### 10/01/09, Ron

"clPairs" <- function(data, classification, symbols,
                      labels = dimnames(data)[[2]], CEX = 1, col, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  data <- as.matrix(data)
  m <- nrow(data)
  n <- ncol(data)
  if(missing(classification))
    classification <- rep(1, m)
  if (!is.factor(classification))
    classification <- as.factor(classification)
  l <- length(levels(classification))
  if(missing(symbols)) {
    if(l == 1) {
      symbols <- "."
    }
    if(l <= length(.Mclust$symbols)) {
      symbols <- .Mclust$symbols
    }
    else if(l <= 9) {
      symbols <- as.character(1:9)
    }
    else if(l <= 26) {
      symbols <- LETTERS[1:l]
    }
    else stop("need more than 26 symbols")
  }
  else if(length(symbols) < l)
    stop("more symbols needed")

  if (missing(col)) col <- 1:l
  if (length(unique(col)) < l & length(unique(col))>1)
    stop("more colors needed")
  if (length(col) == 1) col <- rep(col, l)
  
  pairs(x=data, labels=labels, pch=symbols[classification],
        cex=CEX, col=col[classification], ...)
  invisible()
}

## "clPairs" <- function(data, classification, symbols, 
##                       labels = dimnames(x)[[2]], CEX = 1, PCH = ".", ...)
## {
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##

##   par(pty = "s")
##   x <- as.matrix(data)
##   m <- nrow(x)
##   n <- ncol(x)
##   if(missing(classification))
##     classification <- rep(1, m)
##   u <- sort(unique(classification))
##   l <- length(u)
##   if(missing(symbols)) {
##     if(missing(classification)) {
##       symbols[1] <- PCH
##     }
##     if(l <= length(.Mclust$symbols)) {
##       symbols <- .Mclust$symbols
##     }
##     else if(l <= 9) {
##       symbols <- as.character(1:9)
##     }
##     else if(l <= 26) {
##       symbols <- LETTERS[1:l]
##     }
##     else stop("need more than 26 symbols")
##   }
##   else if(length(symbols) < l)
##     stop("more symbols needed")
##   doaxis <- function(which, dolabel = TRUE)
##     axis(which, outer = TRUE, line = -0.5, labels = dolabel)
##   setup <- function(x, y, ...)
##     .Internal(plot("zplot", range(x[!is.na(x)]), range(y[!is.na(y)]),
##                    type = "n", axes = FALSE, ...), "call_S_Version2")
##   oldpar <- par("oma", "mar", "cex", "tck", "mfg", "mgp", "mex", "mfrow")
##   oldcex <- par("cex")
##   ##
##   ##	CEX <- oldcex * max(7.7/(2. * n + 3.), 0.6)
##   ##
##   par(mfrow = c(n, n), mgp = c(2., 0.8, 0.), oma = rep(3., 4.),
##       mar = rep(0.5, 4.), tck = -0.03/n)  
##   on.exit(par(oldpar))
##   ##
##   ##	par(cex = CEX)
##   ##
##   if(length(labels) < n) labels <- paste(deparse(substitute(data)), "[,",
##                                          1:n, "]", sep = "")
##   if(par("pty") == "s") {
##     dif <- diff(par("fin"))/2
##     if(dif > 0)
##       par(omi = c(dif * n, 0, dif * n, 0) + par("omi"))
##     else par(omi = c(0, ( - dif) * n, 0, ( - dif) * n) + par("omi"))
##   }
##   invert <- list(...)$invert
##   order <- if(is.null(invert) || invert) 1:n else n:1
##   for(i in order) {
##     for(j in 1.:n) {
##       setup(as.vector(x[, j]), as.vector(x[, i]), ...)
##       box()
##       mfg <- par("mfg")
##       if(i == 1)
##         doaxis(3, j %% 2 == 0)
##       if(i == n)
##         doaxis(1, j %% 2 == 1)
##       if(j == 1)
##         doaxis(2, i %% 2 == 0)
##       if(j == n)
##         doaxis(4, i %% 2 == 1)
##       if(i != j) {
##         for(k in 1:l) {
##           r <- (1:m)[classification == u[k]]
##           points(as.vector(x[r, j]), as.vector(x[r, i]), pch = symbols[k],
##                  cex = CEX)
##         }
##       }
##       else {
##         par(usr = c(0., 1., 0., 1.))
##         text(0.5, 0.5, labels[i], cex = CEX)
##       }
##       if(all.equal(par("mfg"), mfg) != TRUE)
##         stop("The panel function made a new plot")
##     }
##   }
##   invisible()
## }


##
## perm is assumed to be a vector of any mode in which no two elements  
## are identical. next.perm produces the next permutation after perm in
## a lexicographic order. If perm is a vector of consecutive positive 
## integers beginning with 1, this order is that of increasing 
## magnitude when each permutation is viewed as an integer in base 
## max(perm)+1 arithmetic.
##
##"nextPerm" <- function(perm)
##{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
##
##n <- length(perm)
##if(n == 1)
##return(perm)
##i <- (n - 1):1
##q <- perm[i + 1] > perm[i]
##if(q[1])
##  return(replace(perm, c(n - 1, n), perm[c(n, n - 1)]))
##if(all(!q))
##  return(rev(perm))
##m <- i[q][1]
##i <- (m + 1):n
##perm[i] <- rev(perm[i])
##l <- i[perm[i] > perm[m]][1]
##replace(perm, c(m, l), perm[c(l, m)])
##}

"coordProj" <- function(data, ..., dimens = c(1, 2),
                        type = c("classification", "uncertainty", "errors"),
                        ask = TRUE, quantiles = c(0.75, 0.95), symbols,
                        scale = FALSE, identify = FALSE, CEX = 1,
                        PCH = ".", xlim, ylim)  
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(scale)
    par(pty = "s")
  aux <- list(...)
  z <- aux$z
  classification <- aux$classification
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  uncertainty <- aux$uncertainty
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  truth <- aux$truth
  mu <- aux$mu
  sigma <- aux$sigma
  decomp <- aux$decomp
  params <- !is.null(mu) && (!is.null(sigma) || !is.null(decomp))
  Data <- data[, dimens]
  if(dim(Data)[2] != 2)
    stop("need two dimensions")
  if(missing(xlim))
    xlim <- range(Data[, 1])
  if(missing(ylim))
    ylim <- range(Data[, 2])
  if(scale) {
    d <- diff(xlim) - diff(ylim)
    if(d > 0) {
      ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
    }
    else {
      xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
    }
  }
  if(is.null(dnames <- dimnames(Data)[[2]]))
    xlab <- ylab <- ""
  else {
    xlab <- dnames[1]
    ylab <- dnames[2]
  }
  if(!is.null(mu)) {
    if(is.null(sigma)) {
      if(is.null(decomp)) {
        params <- FALSE
        warning("covariance not supplied")
      }
      else {
        sigma <- decomp2sigma(decomp)
      }
    }
    G <- ncol(mu)
    mu <- array(mu[dimens, ], c(2,G))
    dimpar <- dim(sigma)
    sigma <- array(sigma[dimens, dimens, ], c(2,2,G))
    if(length(dimpar) != 3) {
      params <- FALSE
      warning("covariance improperly specified")
    }
    if(G != dimpar[3]) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(
                 classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    if(missing(symbols)) {
      if(L <= length(.Mclust$symbols)) {
        symbols <- .Mclust$symbols
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    if(length(symbols) < L) {
      warning("more symbols needed to show classification")
      classification <- NULL
    }
  }
  if(l <- length(type)) {
    choices <- c("classification", "uncertainty", "density", 
                 "errors")
    m <- rep(0, l)
    for(i in 1:l) {
      m[i] <- charmatch(type[i], choices, nomatch = 0)
    }
    choices <- choices[unique(m)]
    if(is.null(classification))
      choices <- choices[choices != "classification"]
    if(is.null(uncertainty))
      choices <- choices[choices != "uncertainty"]
    if(is.null(truth))
      choices <- choices[choices != "errors"]
  }
  else choices <- NULL
  if(length(choices) > 1 && ask)
    choices <- c(choices, "all")
  else {
    if(!length(choices)) {
      plot(Data[, 1], Data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      points(Data[, 1], Data[, 2], pch = PCH, cex = CEX)
      if(identify)
        title(paste("Coordinate Projection: dimens = ",
                    paste(dimens, collapse = ",")), cex = 
              0.5)
      return(invisible())
    }
    if(length(choices) == 1)
      ask <- FALSE
  }
  if(any(choices == "errors")) {
    ERRORS <- classErrors(classification, truth)
  }
  if(!ask)
    pick <- 1:length(choices)
  ALL <- FALSE
  while(TRUE) {
    if(ask) {
      pick <- menu(choices, title = 
                   "\ncoordProj: make a plot selection (0 to exit):\n"
                   )
      if(!pick)
        return(invisible())
      ALL <- any(choices[pick] == "all")
    }
    if(any(choices[pick] == "classification") || (any(choices ==
                    "classification") && ALL)) {
      plot(Data[, 1], Data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      for(k in 1:L) {
        I <- classification == U[k]
        points(Data[I, 1], Data[I, 2], pch = symbols[
                                         k], cex = CEX)
      }
      if(identify)
        title(paste(
                    "Coordinate Projection showing Classification: dimens = ",
                    paste(dimens, collapse = ",")), cex = 
              0.5)
    }
    if(any(choices[pick] == "uncertainty") || (any(choices == 
                    "uncertainty") && ALL)) {
      plot(Data[, 1], Data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      breaks <- quantile(uncertainty, probs = sort(quantiles)
                         )
      I <- uncertainty < breaks[1]
      points(Data[I, 1], Data[I, 2], pch = 16, cex = 0.5 *
             CEX)
      I <- uncertainty < breaks[2] & !I
      points(Data[I, 1], Data[I, 2], pch = 1, cex = 1 * CEX)
      I <- uncertainty >= breaks[2]
      points(Data[I, 1], Data[I, 2], pch = 16, cex = 1.5 *
             CEX)
      if(identify)
        title(paste(
                    "Coordinate Projection showing Uncertainty: dimens = ",
                    paste(dimens, collapse = ",")), cex = 
              0.5)
    }
    if(any(choices[pick] == "errors") || (any(choices == "errors") &&
                    ALL)) {
      plot(Data[, 1], Data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      CLASSES <- unique(as.character(truth))
      symOpen <- c(2, 0, 1, 5)
      symFill <- c(17, 15, 16, 18)
      good <- !ERRORS
      if(L > 4) {
        points(Data[good, 1], Data[good, 2], pch = 1,
               cex = CEX)
        points(Data[!good, 1], Data[!good, 2], pch = 16,
               cex = CEX)
      }
      else {
        for(k in 1:L) {
          K <- truth == CLASSES[k]
          points(Data[K, 1], Data[K, 2], pch = 
                 symOpen[k], cex = CEX)
          if(any(I <- (K & ERRORS))) {
            points(Data[I, 1], Data[I,
                                    2], pch = symFill[
                                          k], cex = CEX)
          }
        }
      }
      if(identify)
        title(paste("Coordinate Projection showing Classification Errors: dimens = ",
                    paste(dimens, collapse = ",")), cex = 0.5)
    }
    if(!ask)
      break
  }
  invisible()
}

"cv1EMtrain" <- function(data, labels, modelNames)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  z <- unmap(as.numeric(labels))
  G <- ncol(z)
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(oneD || length(dimData) != 2) {
    if(missing(modelNames))
      modelNames <- c("E", "V")
    if(any(!match(modelNames, c("E", "V"), nomatch = 0)))
      stop("modelNames E or V for one-dimensional data")
    n <- length(data)
    cv <- matrix(1, nrow = n, ncol = length(modelNames))
    dimnames(cv) <- list(NULL, modelNames)
    for(m in modelNames) {
      for(i in 1:n) {
        mStep <- mstep(modelName = m, data = data[- i], z = z[ - i,],
                       warnSingular = FALSE) 
##        if (m == "V") cat("\n", mStep$sigmasq)
        eStep <- do.call("estep", c(mStep, list(data = data[i],
                                                warnSingular = FALSE)))
        if(is.null(attr(eStep, "warn"))) {
          k <- (1:G)[eStep$z == max(eStep$z)]
          l <- (1:G)[z[i,  ] == max(z[i,  ])]
          cv[i, m] <- as.numeric(!any(k == l))
        } 
      }
    }
  }
  else {
    if(missing(modelNames))
      modelNames <- .Mclust$emModelNames
    n <- nrow(data)
    cv <- matrix(1, nrow = n, ncol = length(modelNames))
    dimnames(cv) <- list(NULL, modelNames)
    for(m in modelNames) {
      for(i in 1:n) {
        mStep <- mstep(modelName = m,
                       data = data[- i,  ], z = z[ - i,  ],
                       warnSingular = FALSE) 
        eStep <- do.call("estep",
                         c(mStep, list(data = data[i,  , drop = FALSE],
                                       warnSingular = FALSE)))
        if(is.null(attr(eStep, "warn"))) {
          k <- (1:G)[eStep$z == max(eStep$z)]
          l <- (1:G)[z[i,  ] == max(z[i,  ])]
          cv[i, m] <- as.numeric(!any(k == l))
        }
      }
    }
  }
  errorRate <- apply(cv, 2, sum)
  ## errorRate[errorRate < 0] <- n
  errorRate/n
}

"decomp2sigma" <- function(d, G, scale, shape, orientation = NULL, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  nod <- missing(d)
  noG <- missing(G)
  lenScale <- length(scale)
  if(lenScale != 1) {
    if(!noG && G != lenScale)
      stop("scale incompatibile with G")
    G <- lenScale
  }
  shape <- as.matrix(shape)
  p <- nrow(shape)
  if(!nod && p != d)
    stop("shape incompatible with d")
  d <- p
  g <- ncol(shape)
  if(g != 1) {
    if(!is.null(G) && g != G)
      stop("shape incompatible with scale")
    if(!noG && g != G)
      stop("shape incompatible with G")
    G <- g
  }
  if(is.null(orientation)) {
    orientName <- "I"
    if(is.null(G)) {
      G <- if(noG) 1 else G
    }
    orientation <- array(diag(d), c(d, d, G))
  }
  else {
    dimO <- dim(orientation)
    l <- length(dimO)
    if(is.null(dimO) || l < 2 || l > 3 || dimO[1] != dimO[2])
      stop("orientation improperly specified")
    if(dimO[1] != d)
      stop("orientation incompatible with shape")
    if(l == 3) {
      orientName <- "V"
      if(is.null(G)) {
        if(!noG && dimO[3] != G)
          stop("orientation incompatible with G")
        G <- dimO[3]
      }
      else if(G != dimO[3])
        stop("orientation incompatible with scale and/or shape"
             )
    }
    else {
      orientName <- "E"
      if(is.null(G)) {
        G <- if(noG) 1 else G
      }
      orientation <- array(orientation, c(d, d, G))
    }
  }
  if(G == 1) {
    scaleName <- shapeName <- "X"
  }
  else {
    scaleName <- if(lenScale == 1) "E" else "V"
    shapeName <- if(g == 1) "E" else "V"
    scale <- rep(scale, G)
    shape <- matrix(shape, nrow = d, ncol = G)
  }
  sigma <- array(0, c(d, d, G))
  for(k in 1:G) {
    sigma[,  , k] <- crossprod(orientation[,  , k] *
                               sqrt(scale[k] * shape[, k]))
  }
  structure(sigma,
            modelName = paste(c(scaleName, shapeName, orientName),
              collapse = ""))
}

"dens" <- function(modelName, data, mu, logarithm=FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  aux <- list(...)
  cden <- do.call("cdens", c(list(modelName = modelName, data = data,
                                  mu = mu),
                             c(list(logarithm=TRUE), aux)))
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  G <- if(oneD) length(mu) else ncol(as.matrix(mu))
  pro <- aux$pro
  if(G > 1) {
    if(is.null(pro))
      stop("mixing proportions must be supplied")
                                        #
    ##	if(is.null(pro <- aux$pro)) pro <- rep(1/G, G)
    ##	apply(sweep(cden, 2, FUN = "*", STATS = pro), 1, sum)
    proz <- !pro
    pro <- pro[!proz]
    cden <- cden[, !proz, drop = FALSE]
    cden <- sweep(cden, 2, FUN = "+", STATS = log(pro))
  }
  maxlog <- apply(cden, 1, max)
  cden <- sweep(cden, 1, FUN = "-", STATS = maxlog)
  den <- log(apply(exp(cden), 1, sum)) + maxlog
  if(is.null(logarithm) || !logarithm)
    den <- exp(den)
  den
}


### Calls R density from base package, unless method=="mclust"
### From R 1.8.1: "density" is part of the "stat" package
### R function density, except for "else" part

"density" <- function(..., method, G)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  aux <- list(...)

  if (missing(G))
    haveG <- FALSE
  else
    haveG <- TRUE

  ## if there is a density function, use it
  if (exists("density", NULL)) {
    densfun <- getFromNamespace("density", ns="base")
  } else {
    if ("stats" %in% .packages(TRUE)) { # it should have a density function
      require(stats, quietly=TRUE)
      densfun <- getFromNamespace("density", ns="stats")
    } else {
      huhn <- getAnywhere("density")
      if (length(huhn$objs) > 0) {
        warning("Using function 'density' from ", huhn$where[1])
        densfun <- huhn$objs[[1]]
      } else {
        stop("Object \"density\" not found")
      }
    }
  }
  val <- do.call("densfun", aux)
    
  if (missing(method)) {
    if (!is.null(class(val))) val$call <- match.call()
    val
  } else {
    if (method != "mclust")
      stop("Unknown method specified")

    x <- aux$x
    if (is.null(x)) {
      if (length(aux)==1 || names(aux)[1] == "") {
        x <- aux[[1]]
      } else {
        stop("argument x is missing, with no default")
      }
    }

    na.rm <- aux$na.rm
    if (is.null(na.rm)) na.rm <- FALSE
    x.na <- is.na(x)
    if (any(x.na)) {
      if (na.rm)
        x <- x[!x.na]
      else stop("x contains missing values")
    }
        
    n <- if (is.null(aux$n)) val$n else aux$n
    bw <- val$bw
    
    name <- deparse(substitute(x))

    if (haveG)
      modl <- summary(EMclust(x, G=G), x)
    else
      modl <- summary(EMclust(x), x)

    x <- seq(min(val$x), max(val$x), length = val$n)
    y <- do.call("dens", c(list(data = x, warnSingular=FALSE), modl))

    structure(list(x=x, y=y, call=match.call(), n=n, data.name=name, bw=bw),
              class="density") 
  } 
}

"em" <- function(modelName, data, mu, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## ... sigmasq or sigma, pro, eps, Vinv
  funcName <- paste("em", modelName, sep = "")
  do.call(funcName, list(data = data, mu = mu, ...))
}

"emE" <- function(data, mu, sigmasq, pro, eps, tol, itmax, equalPro,
                  warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  G <- length(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = 1, G = G, mu = rep(NA, G),
                          sigmasq = NA, pro = rep(NA, K), z = matrix(NA, n, k),
                          loglik = NA, modelName = "E"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(sigmasq <= eps) {
    warn <- "sigma-squared falls below threshold"
    warning("sigma-squared falls below threshold")
    return(structure(list(n = n, d = 1, G = G, mu = rep(NA, G),
                          sigmasq = NA, pro = rep(NA, K), z = matrix(NA, n, k),
                          loglik = NA, modelName = "E"), warn = warn))
  }
  temp <- .Fortran("em1e",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(eps),
                   as.double(tol),
                   double(n * K),
                   PACKAGE="mclust")[6:12]
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  its <- temp[[4]]
  err <- temp[[5]]
  loglik <- temp[[6]]
  z <- matrix(temp[[7]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || sigmasq <= max(eps, 0)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
  }
  else if(its >= itmax) {
    warning("iteration limit reached")
    warn <- "iteration limit reached"
    its <-  - its
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = 1, G = G, mu = mu, sigmasq = sigmasq, pro = 
                 pro, z = z, loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "E"), info = info, warn = warn)
}

"emEEE" <- function(data, mu, Sigma, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  cholSigma <- list(...)$cholSigma
  if(is.null(cholSigma)) {
    if(!is.null(decomp <- list(...)$decomp)) {
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      sig <- qr.R(qr(O * sqrt(scale * shape)))
      cholIND <- "U"
    }
    else if(!missing(Sigma)) {
      sig <- Sigma
      cholIND <- "N"
    }
    else if(!is.null(sigma <- list(...)$sigma)) {
      sig <- sigma
      cholIND <- "N"
    }
    else stop("invalid specification for sigma")
  }
  else {
    sig <- cholSigma
    cholIND <- "U"
  }
  if(any(is.na(c(mu, sig, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = mu,
                          sigma = array(NA, c(p, p, G)),
                          Sigma = matrix(NA, p, p),
                          cholSigma = matrix(NA, p, p),
                          pro = pro, z = matrix(NA, n, K),
                          loglik = NA, modelName = "EEE"), warn = warn))  
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Machine$warnSingular
  storage.mode(mu) <- "double"
  storage.mode(sig) <- "double"
  temp <- .Fortran("emeee",
                   as.integer(if (cholIND == "N") 0 else 1),
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   sig,
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p),
                   double(n * K),
                   PACKAGE="mclust")[8:15]
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholSigma <- structure(temp[[2]], def = 
                         "Sigma = t(cholSigma) %*% cholSigma")
  pro <- temp[[3]]
  its <- temp[[4]]
  err <- temp[[5]]
  loglik <- temp[[6]]
  lapackCholInfo <- temp[[7]][1]
  z <- matrix(temp[[8]], n, K)
  warn <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      if(warnSingular)
        warning("sigma is not positive definite")
      warn <- "sigma is not positive definite"
    }
    else {
      warning("input error for LAPACK DPOTRF")
      warn <- "input error for LAPACK DPOTRF"
    }
    z[] <- loglik <- icond <- NA
    sigma <- array(NA, c(p, p, G))
    cholSigma <- Sigma <- matrix(NA, p, p)
  }
  else {
    Sigma <- unchol(cholSigma, upper = TRUE)
    if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
      if(warnSingular)
        warning("singular covariance")
      warn <- "singular covariance"
      mu[] <- pro[] <- z[] <- loglik <- NA
      sigma <- array(NA, c(p, p, G))
    }
    else {
      sigma <- array(0, c(p, p, G))
      for(k in 1:G)
        sigma[,  , k] <- Sigma
      if(its >= itmax) {
        warning("iteration limit reached")
        warn <- "iteration limit reached"
        its <-  - its
      }
    }
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 Sigma = Sigma, cholSigma = cholSigma, pro = pro,
                 z = z, loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "EEE"), info = info, warn = warn)
}

"emEEI" <- function(data, mu, decomp, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(list(n = n, d = p, G = G, mu = mu,
                sigma = array(NA, c(p, p, G)),
                Sigma = matrix(NA, p, p),
                decomp = list(p = p, G = G, scale = NA, shape = rep(NA, p)),
                pro = pro, z = matrix(NA, n, K),
                loglik = NA, modelName = "EEI"), warn = warn) 
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(mu) <- "double"
  temp <- .Fortran("emeei",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[7:14]
  mu <- temp[[1]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[2]]
  shape <- temp[[3]]
  pro <- temp[[4]]
  its <- temp[[5]]
  err <- temp[[6]]
  loglik <- temp[[7]]
  z <- matrix(temp[[8]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(NA, p, p)
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(rep(scale, p))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, Sigma = 
                 Sigma, decomp = decomp, pro = pro, z = z, loglik = loglik,
                 Vinv = if(noise) Vinv else NULL, modelName = "EEI"), info = 
            info, warn = warn)
}

"emEEV" <- function(data, mu, decomp, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  sigmaIND <- "N"
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = mu,
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(p = p, G = G, scale = NA,
                            shape = rep(NA, p)),
                          pro = pro, z = matrix(NA, n, K), loglik = NA,
                          modelName = "EEV"), warn = warn))  
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p))
  storage.mode(mu) <- "double"
  temp <- .Fortran("emeev",
                   as.integer(if (sigmaIND == "N") 0 else 1),
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(decomp$orientation),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   as.integer(lwork),
                   double(lwork),
                   double(n * K),
                   double(p),
                   PACKAGE="mclust")[8:18]
  mu <- temp[[1]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[2]]
  shape <- temp[[3]]
  O <- array(temp[[4]], c(p, p, G))
  pro <- temp[[5]]
  its <- temp[[6]]
  err <- temp[[7]]
  loglik <- temp[[8]]
  lapackSVDinfo <- temp[[9]]
  z <- matrix(temp[[11]], n, K)
  warn <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      warning("LAPACK DGESVD fails to converge")
      warn <- "LAPACK DGESVD fails to converge"
    }
    else {
      warning("input error for LAPACK DGESVD")
      warn <- "input error for LAPACK DGESVD"
    }
    z[] <- O[] <- shape[] <- NA
    scale <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- scale * shapeO(shape, O, transpose = TRUE)
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  decomp <- structure(list(d = p, G = G, scale = scale, shape = shape,
                           orientation = O), def = 
                      "Sigma = scale * t(O) %*% diag(shape) %*% O")
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 decomp = decomp, pro = pro, z = z, loglik =
                 loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "EEV"), info = info, warn = warn)
}

"emEII" <- function(data, mu, sigmasq, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(sigmasq))
    sigmasq <- list(...)$decomp$scale
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = rep(NA, G),
                          sigma = array(NA, c(p, p, G)), sigmasq = NA,
                          Sigma = matrix(NA, p, p),
                          decomp = list(p = p, G = G, scale = NA),
                          pro = rep(NA, K), z = matrix(NA, n, k),
                          loglik = NA, modelName = "EII"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(sigmasq <= eps) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    return(structure(list(n = n, d = p, G = G, mu = rep(NA, G),
                          sigma = array(NA, c(p, p, G)), sigmasq = NA,
                          Sigma = matrix(NA, p, p), pro = rep(NA, K),
                          z = matrix(NA, n, k), loglik = NA,
                          modelName = "EII"), warn = warn))
  }
  storage.mode(mu) <- "double"
  temp <- .Fortran("emeii",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[7:13]
  mu <- temp[[1]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  Sigma <- diag(rep(sigmasq, p))
  pro <- temp[[3]]
  its <- temp[[4]]
  err <- temp[[5]]
  if(is.infinite(loglik <- temp[[6]]))
    loglik <- NA
  z <- matrix(temp[[7]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax ||
     sigmasq <= max(eps, 0)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = sigmasq)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 sigmasq = sigmasq, Sigma = Sigma, decomp = decomp,
                 pro = pro, z = z, loglik = loglik,
                 Vinv = if(noise) Vinv else NULL,
                 modelName = "EII"), info = info, warn = warn)
}

"emEVI" <- function(data, mu, decomp, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, decomp$scale, decomp$shape, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = mu,
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(p = p, G = G, scale = NA,
                            shape = matrix(NA, p, G)),
                          pro = rep(NA, K), z = matrix(NA, n, K),
                          loglik = NA, modelName = "VEI"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(mu) <- "double"
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  temp <- .Fortran("emevi",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[7:14]
  mu <- temp[[1]]
  scale <- temp[[2]]
  shape <- matrix(temp[[3]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  its <- temp[[5]]
  err <- temp[[6]]
  loglik <- temp[[7]]
  z <- matrix(temp[[8]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, decomp = 
                 decomp, pro = pro, z = z, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "EVI"), info = info, warn
            = warn)
}

"emV" <- function(data, mu, sigmasq, pro, eps, tol, itmax, equalPro,
                  warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensionsal")
  data <- as.vector(data)
  n <- length(data)
  G <- length(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = 1, G = G, mu = rep(NA, G),
                          sigmasq = rep(NA, G), pro = rep(NA, K),
                          z = matrix(NA, n, k), loglik = NA,
                          modelName = "V"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol[1]
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(any(sigmasq <= eps)) {
    warn <- "sigma-squared falls below threshold"
    warning("sigma-squared falls below threshold")
    return(structure(list(n = n, d = 1, G = G, mu = rep(NA, G),
                          sigmasq = rep(NA, G), pro = rep(NA, K),
                          z = matrix(NA, n, k), loglik = NA,
                          modelName = "V"), warn = warn))
  }
  temp <- .Fortran("em1v",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(eps),
                   as.double(tol),
                   double(n * K),
                   PACKAGE="mclust")[6:12]
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  its <- temp[[4]]
  err <- temp[[5]]
  loglik <- temp[[6]]
  z <- matrix(temp[[7]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || any(sigmasq <= max(eps, 0))) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
  }
  else if(its >= itmax) {
    warning("iteration limit reached")
    warn <- "iteration limit reached"
    its <-  - its
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = 1, G = G, mu = mu, sigmasq = sigmasq, pro = 
                 pro, z = z, loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "V"), info = info, warn = warn)
}

"emVEI" <- function(data, mu, decomp, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = mu,
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(p = p, G = G, scale = rep(NA, G),
                            shape = rep(NA, p)), pro = rep(NA, K),
                          z = matrix(NA, n, K), loglik = NA,
                          modelName = "VEI"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  if(length(tol) == 1)
    tol <- c(tol, tol)
  if(missing(itmax))
    itmax <- .Mclust$itmax
  if(length(itmax) == 1)
    itmax <- c(itmax, Inf)
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(mu) <- "double"
  temp <- .Fortran("emvei",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(n * K),
                   double(G),
                   double(p),
                   double(p * G),
                   PACKAGE="mclust")[7:14]
  mu <- temp[[1]]
  scale <- temp[[2]]
  shape <- temp[[3]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  its <- temp[[5]][1]
  inner <- temp[[5]][2]
  err <- temp[[6]][1]
  inerr <- temp[[6]][2]
  loglik <- temp[[7]]
  z <- matrix(temp[[8]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(NA, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= itmax[2]) {
      warning("inner iteration limit reached")
      warn <- "inner iteration limit reached"
      inner <-  - inner
    }
    else if(its >= itmax[1]) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- structure(c(iterations = its, error = err),
                    inner = c(iterations = inner, error = inerr))
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, decomp = 
                 decomp, pro = pro, z = z, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "VEI"), info = info, warn
            = warn)
}

"emVEV" <- function(data, mu, decomp, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  sigmaIND <- "N"
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = mu,
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(p = p, G = G, scale = NA,
                            shape = rep(NA, p)),
                          pro = rep(NA, K), z = matrix(NA, n, K),
                          loglik = NA, modelName = "VEV"),
                     warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  if(length(tol) == 1)
    tol <- c(tol, tol)
  if(missing(itmax))
    itmax <- .Mclust$itmax
  if(length(itmax) == 1)
    itmax <- c(itmax, Inf)
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  storage.mode(mu) <- "double"
  temp <- .Fortran("emvev",
                   as.integer(if (sigmaIND == "N") 0 else 1),
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(decomp$orientation),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   as.integer(lwork),
                   double(lwork),
                   double(n * K),
                   double(p * G),
                   PACKAGE="mclust")[8:18]
  mu <- temp[[1]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[2]]
  shape <- temp[[3]]
  O <- array(temp[[4]], c(p, p, G))
  pro <- temp[[5]]
  its <- temp[[6]][1]
  inner <- temp[[6]][2]
  err <- temp[[7]][1]
  inerr <- temp[[7]][2]
  loglik <- temp[[8]]
  lapackSVDinfo <- temp[[9]]
  z <- matrix(temp[[11]], n, K)
  warn <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      warning("LAPACK DGESVD fails to converge")
      warn <- "LAPACK DGESVD fails to converge"
    }
    else {
      warning("input error for LAPACK DGESVD")
      warn <- "input error for LAPACK DGESVD"
    }
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- shapeO(shape, O, transpose = TRUE)
    sigma <- sweep(sigma, MARGIN = 3, STATS = scale, FUN = "*")
    if(inner >= itmax[2]) {
      warning("inner iteration limit reached")
      warn <- "inner iteration limit reached"
      inner <-  - inner
    }
    else if(its >= itmax[1]) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  decomp <- structure(list(d = p, G = G, scale = scale, shape = shape,
                           orientation = O), def = 
                      "Sigma = scale * t(O) %*% diag(shape) %*% O")
  info <- structure(c(iterations = its, error = err),
                    inner = c(iterations = inner, error = inerr))
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, decomp = 
                 decomp, pro = pro, z = z, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "VEV"), info = info, warn
            = warn)
}

"emVII" <- function(data, mu, sigmasq, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(sigmasq))
    sigmasq <- list(...)$decomp$scale
  if(all(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = rep(NA, G),
                          sigma = array(NA, c(p, p, G)), sigmasq = rep(NA, G),
                          decomp = list(p = p, G = G, scale = rep(NA, G)),
                          pro = rep(NA, K), z = matrix(NA, n, K),
                          loglik = NA, modelName = "VII"), warn = warn))
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    stop("parameters contain missing values")
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(any(sigmasq <= eps)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    return(structure(list(n = n, d = p, G = G, mu = rep(NA, G),
                          sigma = array(NA, c(p, p, G)), sigmasq = rep(NA, G),
                          decomp = list(p = p, G = G, scale = rep(NA, G)),
                          pro = rep(NA, K), z = matrix(NA, n, K), loglik = NA, 
                          modelName = "VII"), warn = warn))
  }
  storage.mode(mu) <- "double"
  temp <- .Fortran("emvii",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[7:13]
  mu <- temp[[1]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  its <- temp[[4]]
  err <- temp[[5]]
  if(is.infinite(loglik <- temp[[6]]))
    loglik <- NA
  z <- matrix(temp[[7]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax ||
     any(sigmasq <= max(eps, 0))) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(rep(sigmasq[k], p))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = sigmasq)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 sigmasq = sigmasq, decomp = decomp, pro = pro, z = z,
                 loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "VII"), info = info, warn = warn)
}

"emVVI" <- function(data, mu, decomp, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G, mu = mu,
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(p = p, G = G, scale = rep(NA, G),
                            shape = matrix(NA, p, G)),
                          pro = rep(NA, K), z = matrix(NA, n, K),
                          loglik = NA, modelName = "VVI"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Machine$warnSingular
  storage.mode(mu) <- "double"
  temp <- .Fortran("emvvi",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[7:14]
  mu <- temp[[1]]
  scale <- temp[[2]]
  shape <- matrix(temp[[3]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  its <- temp[[5]]
  err <- temp[[6]]
  loglik <- temp[[7]]
  z <- matrix(temp[[8]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
                               FUN = "*"), 2, diag), c(p, p, G))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, decomp = 
                 decomp, pro = pro, z = z, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "VVI"), info = info, warn
            = warn)
}

"emVVV" <- function(data, mu, sigma, pro, eps, tol, itmax, equalPro,
                    warnSingular, Vinv, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  cholsigma <- list(...)$cholsigma
  if(is.null(cholsigma)) {
    if(missing(sigma)) {
      decomp$list(...)$decomp
      if(is.null(decomp))
        stop("covariance improperly specified")
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      sig <- array(0, c(p, p, G))
      shape <- sqrt(sweep(shape, MARGIN = 2, STATS = scale,
                          FUN = "*"))
      for(k in 1:G)
        sig[,  , k] <- qr.R(qr(O[,  , k] * shape))
      cholIND <- "U"
    } else {
      sig <- sigma
      cholIND <- "N"
    }
  }
  else {
    sig <- cholsigma
    cholIND <- "U"
  }
  if(all(is.na(c(mu, sig, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(n = n, d = p, G = G,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          cholsigma = array(NA, c(p, p, G)),
                          z = matrix(NA, n, K), loglik = NA,
                          modelName = "VVV"), warn = warn))
  }
  if(any(is.na(c(mu, sig, pro)))) {
    stop("parameters contain missing values")
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(mu) <- "double"
  storage.mode(sig) <- "double"
  temp <- .Fortran("emvvv",
                   as.integer(if (cholIND == "N") 0 else 1),
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   mu,
                   sig,
                   as.double(pro),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p),
                   double(n * K),
                   PACKAGE="mclust")[8:15]
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- structure(array(temp[[2]], c(p, p, G)), def = 
                         "Sigma = t(cholsigma) %*% cholsigma")
  pro <- temp[[3]]
  its <- temp[[4]]
  err <- temp[[5]]
  loglik <- temp[[6]]
  lapackCholInfo <- temp[[7]][1]
  z <- matrix(temp[[8]], n, K)
  warn <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      if(warnSingular)
        warning("sigma is not positive definite")
      warn <- "sigma is not positive definite"
    }
    else {
      warn <- "input error for LAPACK DPOTRF"
      warning("input error for LAPACK DPOTRF")
    }
    z[] <- loglik <- icond <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(cholsigma, 3, unchol, TRUE), c(p, p, G))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 cholsigma = cholsigma, pro = pro, z = z,
                 loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "VVV"), info = info, warn = warn)
}

"estep" <- function(modelName, data, mu, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  funcName <- paste("estep", modelName, sep = "")
  do.call(funcName, list(data = data, mu = mu, ...))
}

"estepE" <- function(data, mu, sigmasq, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  pro <- pro/sum(pro)
  l <- length(pro)
  data <- as.vector(data)
  n <- length(data)
  G <- length(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "E"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(sigmasq <= eps) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    return(structure(matrix(NA, n, K), loglik = NA, modelName = "E",
                     warn = warn))
  }
  temp <- .Fortran("es1e",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = 1, G = G, z = z, loglik = loglik, modelName
		 = "E"), warn = warn)
}

"estepEEE" <- function(data, mu, Sigma, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) > 2)
    stop("data must be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  cholSigma <- list(...)$cholSigma
  if(is.null(cholSigma)) {
    if(!is.null(decomp <- list(...)$decomp)) {
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      sig <- qr.R(qr(O * sqrt(scale * shape)))
      cholIND <- "U"
    }
    else if(!is.null(Sigma <- list(...)$Sigma)) {
      sig <- Sigma
      cholIND <- "N"
    }
    else if(!missing(sigma)) {
      sig <- sigma
      cholIND <- "N"
    }
    else stop("invalid specification for sigma")
  }
  else {
    sig <- cholSigma
    cholIND <- "U"
  }
  if(any(is.na(c(mu, sig, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "EEE"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("eseee",
                   as.integer(if (cholIND == "N") 0 else 1),
                   as.double(data),
                   as.double(mu),
                   as.double(sig),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   double(p),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, K)
  warn <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      warn <- "sigma is not positive definite"
      warning("sigma is not positive definite")
    }
    else {
      warn <- "input error for LAPACK DPOTRF"
      warning("input error for LAPACK DPOTRF")
    }
    z[] <- loglik <- NA
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "EEE"), warn = warn)
}

"estepEEI" <- function(data, mu, decomp, pro, eps, warnSingular, Vinv,
                       ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "EEI"), warn = warn))
  }
  temp <- .Fortran("eseei",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "EEI"), warn = warn)
}

"estepEEV" <- function(data, mu, decomp, pro, eps, warnSingular, Vinv,
                       ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "EEV"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("eseev",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(decomp$orientation),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   double(p),
                   double(p),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "EEV"), warn = warn)
}

"estepEII" <- 
  function(data, mu, sigmasq, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(sigmasq)) {
    sigmasq <- list(...)$decomp$scale
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "EII"), warn = warn))
  }
  ## Changed 09/30/02, Ron
  ##  if(list(...)$decomp$scale <= eps) {
  if(sigmasq <= eps) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    return(structure(matrix(NA, n, K), loglik = NA, modelName = 
                     "EII", warn = warn))
  }
  temp <- .Fortran("eseii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "EII"), warn = warn)
}

"estepEVI" <- function(data, mu, decomp, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "EVI"), warn = warn))
  }
  temp <- .Fortran("esevi",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "EVI"), warn = warn)
}

"estepV" <- function(data, mu, sigmasq, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  pro <- pro/sum(pro)
  l <- length(pro)
  data <- as.vector(data)
  n <- length(data)
  G <- length(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "V"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(any(sigmasq <= eps)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    return(structure(matrix(NA, n, K), loglik = NA, modelName = "V",
                     warn = warn))
  }
  temp <- .Fortran("es1v",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold: second position")
    warn <- "sigma-squared falls below threshold"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = 1, G = G, z = z, loglik = loglik,
                 modelName = "V"), warn = warn)
}

"estepVEI" <- function(data, mu, decomp, pro, eps, warnSingular, Vinv,
                       ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "VEI"), warn = warn))
  }
  temp <- .Fortran("esvei",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "VEI"), warn = warn)
}

"estepVEV" <- function(data, mu, decomp, pro, eps, warnSingular, Vinv,
                       ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "VEV"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("esvev",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(decomp$orientation),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   double(p),
                   double(p),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "VEV"), warn = warn)
}

"estepVII" <- function(data, mu, sigmasq, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(sigmasq)) {
    sigmasq <- list(...)$decomp$scale
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "VII"), warn = warn))
  }
  ##  if(any(list(...)$decomp$scale <= eps)) {
  if(any(sigmasq <= eps)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    return(structure(matrix(NA, n, K), loglik = NA, modelName = 
                     "VII", warn = warn))
  }
  temp <- .Fortran("esvii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "VII"), warn = warn)
}

"estepVVI" <- function(data, mu, decomp, pro, eps, warnSingular, Vinv,
                       ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  if(missing(decomp))
    stop("decomp must be specified")
  if(any(is.na(c(mu, unlist(decomp), pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "VVI"), warn = warn))
  }
  temp <- .Fortran("esvvi",
                   as.double(data),
                   as.double(mu),
                   as.double(decomp$scale),
                   as.double(decomp$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "VVI"), warn = warn)
}

"estepVVV" <- function(data, mu, sigma, pro, eps, warnSingular, Vinv, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(mu)
  G <- ncol(mu)
  pro <- pro/sum(pro)
  l <- length(pro)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- -1
  }
  else {
    K <- G + 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  cholsigma <- list(...)$cholsigma
  if(is.null(cholsigma)) {
    if(missing(sigma)) {
      if(!is.null(sigma <- list(...)$sigma)) {
        sig <- sigma
        cholIND <- "N"
      }
      else if(!is.null(decomp <- list(...)$decomp)) {
        scale <- decomp$scale
        shape <- decomp$shape
        O <- decomp$orientation
        sig <- array(0, c(p, p, G))
        shape <- sqrt(sweep(shape, MARGIN = 2, STATS = 
                            scale, FUN = "*"))
        for(k in 1:G)
          sig[,  , k] <- qr.R(qr(O[,  , k] * 
                                 shape))
        cholIND <- "U"
      }
      else stop("sigma improperly specified")
    }
    else {
      sig <- sigma
      cholIND <- "N"
    }
  }
  else {
    sig <- cholsigma
    cholIND <- "U"
  }
  if(any(is.na(c(mu, sig, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(list(z = matrix(NA, n, K), loglik = NA, 
                          modelName = "VVV"), warn = warn))
  }
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("esvvv",
                   as.integer(if (cholIND == "N") 0 else 1),
                   as.double(data),
                   as.double(mu),
                   as.double(sig),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   double(p),
                   as.double(eps),
                   double(n * K),
                   PACKAGE="mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, K)
  warn <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      warn <- "sigma is not positive definite"
      warning("sigma is not positive definite")
    }
    else {
      warn <- "input error for LAPACK DPOTRF"
      warning("input error for LAPACK DPOTRF")
    }
    z[] <- loglik <- NA
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    z[] <- loglik <- NA
  }
  structure(list(n = n, d = p, G = G, z = z, loglik = loglik, modelName
		 = "VVV"), warn = warn)
}

"grid1" <- function(n, range = c(0., 1.), edge = TRUE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(any(n < 0. | round(n) != n))
    stop("n must be nonpositive and integer")
  G <- rep(0., n)
  if(edge) {
    G <- seq(from = min(range), to = max(range), by = abs(diff(
                                                   range))/(n - 1.))
  }
  else {
    lj <- abs(diff(range))
    incr <- lj/(2. * n)
    G <- seq(from = min(range) + incr, to = max(range) - incr,
             by = 2. * incr)
  }
  G
}

"grid2" <- function(x, y)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  lx <- length(x)
  ly <- length(y)
  xy <- matrix(0, nrow = lx * ly, ncol = 2)
  l <- 0
  for(j in 1:ly) {
    for(i in 1:lx) {
      l <- l + 1
      xy[l,  ] <- c(x[i], y[j])
    }
  }
  xy
}

"hc" <- function(modelName, data, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## ... partition, minclus = 1, 
  do.call(paste("hc", modelName, sep = ""), list(data, ...))
}

"hcE" <- function(data, partition, minclus = 1, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(minclus < 1)
    stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
###====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data mist be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hc1e",
                   data,
                   as.integer(n),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.integer(ld),
                   double(ld),
                   PACKAGE="mclust")[c(1, 3, 7)]
  temp[[1]] <- temp[[1]][1:m]
  temp[[2]] <- temp[[2]][1:m]
  temp[[3]] <- temp[[3]][1:m]
  structure(rbind(temp[[1]], temp[[2]]), change = temp[[3]], 
            initialPartition = partition, dimensions = n, modelName = "E",
            class = "hc")
}

"hcEEE" <- function(data, partition, minclus = 1, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(minclus < 1)
    stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
###=====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p)
    warning("# of observations <= data dimension")
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  temp <- .Fortran("hceee",
                   data,
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   if(p < 3) integer(m) else integer(1),
                   if(p < 4) integer(m) else integer(1),
                   double(p),
                   double(p * p),
                   double(p * p),
                   double(p * p),
                   PACKAGE="mclust")[c(1, 7:10)]
                                        #
                                        # currently temp[[5]] is not output
  temp[[4]] <- temp[[4]][1:2]
  temp[[5]] <- temp[[5]][1:2]
  names(temp[[5]]) <- c("determinant", "trace")
  temp[[1]] <- temp[[1]][1:(m + 1),  ]
  if(p < 3)
    tree <- rbind(temp[[2]], temp[[3]])
  else if(p < 4)
    tree <- rbind(temp[[1]][-1, 3], temp[[3]])
  else tree <- t(temp[[1]][-1, 3:4, drop = FALSE])
  determinant <- temp[[1]][, 1]
  attr(determinant, "breakpoints") <- temp[[4]]
  structure(tree, determinant = determinant, trace = temp[[1]][, 2],
            initialPartition = partition, dimensions = dimdat, modelName = 
            "EEE", class = "hc")
}

"hcEII" <- function(data, partition, minclus = 1, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(minclus < 1)
    stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
###====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  if(n <= p)
    warning("# of observations <= data dimension")
###=============================================================
  storage.mode(data) <- "double"
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hceii",
                   data,
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   double(p),
                   as.integer(ld),
                   double(ld),
                   PACKAGE="mclust")[c(1, 9)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), change = temp[[2]], initialPartition = 
            partition, dimensions = dimdat, modelName = "EII", class = "hc"
            )
}

"hcV" <- function(data, partition, minclus = 1, alpha = 1, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(minclus < 1)
    stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
###=====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  alpha <- alpha * (vecnorm(data - mean(data))^2/n)
  alpha <- min(alpha, .Machine$double.eps)
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hc1v",
                   data,
                   as.integer(n),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.double(alpha),
                   as.integer(ld),
                   double(ld),
                   PACKAGE="mclust")[c(1, 3, 8)]
  temp[[1]] <- temp[[1]][1:m]
  temp[[2]] <- temp[[2]][1:m]
  temp[[3]] <- temp[[3]][1:m]
  structure(rbind(temp[[1]], temp[[2]]), change = temp[[3]], 
            initialPartition = partition, dimensions = n, modelName = "V",
            class = "hc")
}

"hcVII" <- function(data, partition, minclus = 1, alpha = 1, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(minclus < 1)
    stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
###=====================================================================
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p)
    warning("# of observations <= data dimension")
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ld <- max(n, ll, 3 * m)
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvii",
                   data,
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.double(alpha),
                   double(p),
                   as.integer(ld),
                   double(ld),
                   PACKAGE="mclust")[c(1, 10)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), change = temp[[2]], initialPartition = 
            partition, dimensions = dimdat, modelName = "VII", class = "hc"
            )
}

"hcVVV" <- function(data, partition, minclus = 1, alpha = 1, beta = 1,
                    ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(minclus < 1)
    stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p)
    warning("# of observations <= data dimension")
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    stop("initial number of clusters is not greater than minclus")
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ##	dp <- duplicated(partition)
  ##x[c((1:n)[!dp],(1:n)[dp]), ], 
  ##as.integer(c(partition[!dp], partition[dp])), 
  ld <- max(n, ll + 1, 3 * m)
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvvv",
                   cbind(data, 0.),
                   as.integer(n),
                   as.integer(p),
                   as.integer(partition),
                   as.integer(l),
                   as.integer(m),
                   as.double(alpha),
                   as.double(beta),
                   double(p),
                   double(p * p),
                   double(p * p),
                   double(p * p),
                   as.integer(ld),
                   double(ld),
                   PACKAGE="mclust")[c(1, 14)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), change = temp[[2]], initialPartition = 
            partition, dimensions = dimdat, modelName = "VVV", class = "hc"
            )
}

"hclass" <- function(hcPairs, G)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  initial <- attributes(hcPairs)$init
  n <- length(initial)
  k <- length(unique(initial))
  G <- if(missing(G)) k:2 else rev(sort(unique(G)))
  select <- k - G
  if(length(select) == 1 && !select)
    return(matrix(initial, ncol = 1, dimnames = list(NULL, 
                                       as.character(G))))
  bad <- select < 0 | select >= k
  if(all(bad))
    stop("No classification with the specified number of clusters")
  if(any(bad))
    warning("Some selected classifications are inconsistent\n                          with mclust object")
  L <- length(select)
  cl <- matrix(NA, nrow = n, ncol = L, 
               dimnames = list(NULL, as.character(G)))
  if(select[1])
    m <- 1
  else {
    cl[, 1] <- initial
    m <- 2
  }
  for(l in 1:max(select)) {
    ij <- hcPairs[, l]
    i <- min(ij)
    j <- max(ij)
    initial[initial == j] <- i
    if(select[m] == l) {
      cl[, m] <- initial
      m <- m + 1
    }
  }
  apply(cl[, L:1, drop = FALSE], 2, partconv, consec = TRUE)
}

"hypvol" <- function(data, reciprocal = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## finds the minimum hypervolume between principal components and 
  ## variable bounds
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD) {
    ## 1D case
    n <- length(as.vector(data))
    if(reciprocal) {
      ans <- 1/(max(data) - min(data))
    }
    else {
      ans <- max(data) - min(data)
    }
    return(ans)
  }
  if(length(dimdat) != 2)
    stop("data must be a vector or a matrix")
  data <- as.matrix(data)
  dimd <- dim(data)
  n <- dimd[1]
  p <- dimd[2]
  if(FALSE) {
    vol1 <- prod(apply(data, 2, function(z)
                       diff(range(z))))
    V <- matrix(temp[[1]], p, p)
    xbar <- apply(data, 2, mean)
    X <- sweep(data, 2, xbar)
    library(Matrix)
    print(V)
    print(eigen.Hermitian(crossprod(X))$vectors)
    X <- X %*% V
    vol <- prod(apply(X, 2, function(z)
                      diff(range(z))))
  }
  lwgesvd <- max(3 * min(n, p) + max(n, p), 5 * min(n, p) - 4)
                                        # min
  lwsyevd <- p * (3 * p + 2 * ceiling(logb(p, base = 2)) + 5) + 1
                                        # minimum
  lisyevd <- 5 * p + 3
                                        # minimum
  lwsyevx <- 8 * p
  lisyevx <- 5 * p + p
  lwork <- max(lwsyevd, lwsyevx, n)
  liwork <- max(lisyevd,lisyevx)
  temp <- .Fortran("mclvol",
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   double(p),
                   double(p * p),
                   double(p * p),
                   double(lwork),
                   as.integer(lwork),
                   integer(liwork),
                   as.integer(liwork),
                   integer(1),
                   PACKAGE="mclust")[c(4, 11)]
  if(temp[[2]])
    stop("problem in computing principal components")
  if(reciprocal) {
    pcvol <- prod(1/temp[[1]])
    bdvol <- prod(1/(apply(data, 2, max) - apply(data, 2, min)))
    ans <- max(pcvol, bdvol)
  }
  else {
    pcvol <- prod(temp[[1]])
    bdvol <- prod(apply(data, 2, max) - apply(data, 2, min))
    ans <- min(pcvol, bdvol)
  }
  ans
}


"map" <- function(z, warn=TRUE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ##
  ## converts conditional probabilities to a classification
  ##
  nrowz <- nrow(z)
  cl <- numeric(nrowz)
  I <- 1:nrowz
  J <- 1:ncol(z)
  for(i in I) {
    cl[i] <- (J[z[i,  ] == max(z[i,  ])])[1]
  }
  if(warn) {
    K <- as.logical(match(J, sort(unique(cl)), nomatch = 0))
    if(any(!K))
      warning(paste("no assignment to", paste(J[!K], collapse = ",")))
  }
  cl
}

"mclust1Dplot" <- function(data, ..., type = c("classification",
                                        "uncertainty", "density",
                                        "errors"),
                           ask = TRUE, symbols, grid = 100, identify = FALSE,
                           CEX = 1, xlim)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  densNuncer <- function(data, mu, sigmasq, pro)
    {
      cden <- cdensV(data = data, mu = mu, sigmasq = sigmasq, pro = 
                     pro)
      z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = pro)
      den <- apply(z, 1, sum)
      z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
      data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
    }
  p <- ncol(as.matrix(data))
  if(p != 1)
    stop("for one-dimensional data only")
  data <- as.vector(data)
  n <- length(data)
  aux <- list(...)
  z <- aux$z
  classification <- aux$classification
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  uncertainty <- aux$uncertainty
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  truth <- aux$truth
  mu <- as.vector(aux$mu)
  sigma <- as.vector(aux$sigma)
  pro <- aux$pro
  params <- !is.null(mu) && !is.null(sigma) && !is.null(pro)
  if(params) {
    G <- length(mu)
    if((l <- length(sigma)) == 1) {
      sigma <- rep(sigma, G)
    }
    else if(l != G) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(
                 classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- unique(classification)
    L <- length(U)
    if(missing(symbols)) {
      symbols <- rep("|", L)
      if(FALSE) {
        if(L <= length(.Mclust$symbols)) {
          symbols <- .Mclust$symbols
        }
        else if(L <= 9) {
          symbols <- as.character(1:9)
        }
        else if(L <= 26) {
          symbols <- LETTERS
        }
      }
    }
    if(length(symbols) < L) {
      warning("more symbols needed to show classification")
      classification <- NULL
    }
  }
  if(l <- length(type)) {
    choices <- c("classification", "uncertainty", "density", 
                 "errors")
    m <- rep(0, l)
    for(i in 1:l) {
      m[i] <- charmatch(type[i], choices, nomatch = 0)
    }
    choices <- choices[unique(m)]
    if(is.null(classification))
      choices <- choices[choices != "classification"]
    if(is.null(truth))
      choices <- choices[choices != "errors"]
  }
  else choices <- NULL
  if(!params) {
    choices <- choices[choices != "uncertainty"]
    choices <- choices[choices != "density"]
  }
  if(length(choices) > 1 && ask)
    choices <- c(choices, "all")
  else {
    if(!length(choices)) {
      plot(data, rep(0, n), type = "n", xlab = "", ylab = "",
           xlim = xlim, ...)
      points(data, rep(0, n), pch = "|", cex = CEX)
      if(identify)
        title("Point Plot", cex = 0.5)
      return(invisible())
    }
    if(length(choices) == 1)
      ask <- FALSE
  }
  if(any(choices == "errors")) {
    ERRORS <- classErrors(classification, truth)
  }
  if(!ask)
    pick <- 1:length(choices)
  ALL <- FALSE
  if(missing(xlim))
    xlim <- range(data)
  while(TRUE) {
    if(ask) {
      pick <- menu(choices, title = 
                   "\nmclust1Dplot: make a plot selection (0 to exit):\n"
                   )
      if(!pick)
        return(invisible())
      ALL <- any(choices[pick] == "all")
    }
    if(any(choices[pick] == "classification") || (any(choices ==
                    "classification") && ALL)) {
      plot(data, seq(from = 0, to = L, length = n), type = 
           "n", xlab = "", ylab = "", xlim = xlim, yaxt = 
           "n", ...)
      for(k in 1:L) {
        I <- classification == U[k]
        points(data[I], rep(0, length(data[I])), pch = 
               symbols[k], cex = CEX)
        points(data[I], rep(k, length(data[I])), pch = 
               symbols[k], cex = CEX)
      }
      if(identify)
        title("Classification", cex = 0.5)
    }
    if(any(choices[pick] == "errors") || (any(choices == "errors") &&
                    ALL)) {
      plot(data, seq(from = 0, to = L, length = n), type = 
           "n", xlab = "", ylab = "", xlim = xlim, yaxt = 
           "n", ...)
      good <- !ERRORS
      sym <- "|"
      for(k in 1:L) {
        K <- classification == U[k]
        I <- K & good
        if(any(I)) {
          if(FALSE) {
            sym <- if(L > 4) 1 else if(k ==
                                       4)
              5
            else k - 1
          }
          l <- sum(as.numeric(I))
          ##
          ## points(data[I], rep(k, l), pch = sym, cex = CEX)
          ##
          points(data[I], rep(0, l), pch = sym,
                 cex = CEX)
        }
        I <- K & !good
        if(any(I)) {
          if(FALSE)
            sym <- if(L > 5) 16 else k +
              14
          l <- sum(as.numeric(I))
          points(data[I], rep(k, l), pch = sym,
                 cex = CEX)
          points(data[I], rep(0, l), pch = sym,
                 cex = CEX)
        }
      }
      if(identify)
        title("Classification Errors", cex = 0.5)
    }
    if((any(choices == "uncertainty" | choices == "density") && ALL
        ) || any(choices[pick] == "uncertainty" | choices[
                          pick] == "density")) {
      x <- grid1(n = grid, range = xlim, edge = TRUE)
      lx <- length(x)
      Z <- densNuncer(data = x, mu = mu, sigmasq = sigma,
                      pro = pro)
    }
    if(any(choices[pick] == "uncertainty") || (any(choices == 
                    "uncertainty") && ALL)) {
      plot(x, Z$uncertainty, xlab = "", ylab = "uncertainty",
           xlim = xlim, type = "l", ...)
      if(identify)
        title("Uncertainty", cex = 0.5)
    }
    if(any(choices[pick] == "density") || (any(choices == "density"
                    ) && ALL)) {
      plot(x, Z$density, xlab = "", ylab = "density", xlim = 
           xlim, type = "l", ...)
      if(identify)
        title("Density", cex = 0.5)
    }
    if(!ask)
      break
  }
  invisible()
}

"mclust2Dplot" <- function(data, ..., type = c("classification",
                                        "uncertainty", "errors"),
                           ask = TRUE, quantiles = c(0.75, 0.95),
                           symbols, scale = FALSE, identify = FALSE,
                           CEX = 1, PCH = ".", xlim, ylim, swapAxes = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(scale)
    par(pty = "s")
  data <- as.matrix(data)
  p <- ncol(data)
  if(p != 2)
    stop("for two-dimensional data only")
  aux <- list(...)
  z <- aux$z
  classification <- aux$classification
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  uncertainty <- aux$uncertainty
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  truth <- aux$truth
  mu <- aux$mu
  sigma <- aux$sigma
  decomp <- aux$decomp
  params <- !is.null(mu) && (!is.null(sigma) || !is.null(decomp))
  if(!is.null(mu)) {
    if(is.null(sigma)) {
      if(is.null(decomp)) {
        params <- FALSE
        warning("covariance not supplied")
      }
      else {
        sigma <- decomp2sigma(decomp)
      }
    }
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      params <- FALSE
      warning("covariance improperly specified")
    }
    if(G != dimpar[3]) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  }
  if(swapAxes) {
    if(params) {
      mu <- mu[2:1,  ]
      sigma <- sigma[2:1, 2:1,  ]
    }
    data <- data[, 2:1]
  }
  if(!is.null(dnames <- dimnames(data)[[2]])) {
    xlab <- dnames[1]
    ylab <- dnames[2]
  }
  else xlab <- ylab <- ""
  if(missing(xlim))
    xlim <- range(data[, 1])
  if(missing(ylim))
    ylim <- range(data[, 2])
  if(scale) {
    d <- diff(xlim) - diff(ylim)
    if(d > 0) {
      ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
    }
    else {
      xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
    }
  }
  if(!is.null(truth)) {
    truth <- as.character(truth)
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      classification <- as.character(classification)
      if(length(unique(truth)) != length(unique(
                 classification)))
        truth <- NULL
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    if(missing(symbols)) {
      if(L <= length(.Mclust$symbols)) {
        symbols <- .Mclust$symbols
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    if(length(symbols) < L) {
      warning("more symbols needed to show classification")
      classification <- NULL
    }
  }
  if(l <- length(type)) {
    choices <- c("classification", "uncertainty", "errors")
    m <- rep(0, l)
    for(i in 1:l) {
      m[i] <- charmatch(type[i], choices, nomatch = 0)
    }
    choices <- choices[unique(m)]
    if(is.null(classification))
      choices <- choices[choices != "classification"]
    if(is.null(uncertainty))
      choices <- choices[choices != "uncertainty"]
    if(is.null(truth))
      choices <- choices[choices != "errors"]
  }
  else choices <- NULL
  if(length(choices) > 1 && ask)
    choices <- c(choices, "all")
  else {
    if(!length(choices)) {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      points(data[, 1], data[, 2], pch = PCH, cex = CEX)
      if(identify)
        title("Point Plot", cex = 0.5)
      return(invisible())
    }
    if(length(choices) == 1)
      ask <- FALSE
  }
  if(any(choices == "errors")) {
    ERRORS <- classErrors(classification, truth)
  }
  if(!ask)
    pick <- 1:length(choices)
  ALL <- FALSE
  while(TRUE) {
    if(ask) {
      pick <- menu(choices, title = 
                   "\nmclust2Dplot: make a plot selection (0 to exit):\n"
                   )
      if(!pick)
        return(invisible())
      ALL <- any(choices[pick] == "all")
    }
    if(any(choices[pick] == "classification") || (any(choices ==
                    "classification") && ALL)) {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      for(k in 1:L) {
        I <- classification == U[k]
        points(data[I, 1], data[I, 2], pch = symbols[
                                         k], cex = CEX)
      }
      if(identify)
        title("Classification", cex = 0.5)
    }
    if(any(choices[pick] == "uncertainty") || (any(choices == 
                    "uncertainty") && ALL)) {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      breaks <- quantile(uncertainty, probs = sort(quantiles)
                         )
      I <- uncertainty < breaks[1]
      points(data[I, 1], data[I, 2], pch = 16, cex = 0.5 *
             CEX)
      I <- uncertainty < breaks[2] & !I
      points(data[I, 1], data[I, 2], pch = 1, cex = 1 * CEX)
      I <- uncertainty >= breaks[2]
      points(data[I, 1], data[I, 2], pch = 16, cex = 1.5 *
             CEX)
      if(identify)
        title("Classification Uncertainty", cex = 0.5)
    }
    if(any(choices[pick] == "errors") || (any(choices == "errors") &&
                    ALL)) {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab,
           ylab = ylab, xlim = xlim, ylim = ylim, ...)
      if(params) {
        for(k in 1:G) {
          mvn2plot(mu = mu[, k], sigma = sigma[
                                   ,  , k], k = 15)
        }
      }
      CLASSES <- unique(as.character(truth))
      symOpen <- c(2, 0, 1, 5)
      symFill <- c(17, 15, 16, 18)
      good <- ERRORS
      if(L > 4) {
        points(data[good, 1], data[good, 2], pch = 1,
               cex = CEX)
        points(data[!good, 1], data[!good, 2], pch = 16,
               cex = CEX)
      }
      else {
        for(k in 1:L) {
          K <- truth == CLASSES[k]
          points(data[K, 1], data[K, 2], pch = 
                 symOpen[k], cex = CEX)
          if(any(I <- (K & ERRORS))) {
            points(data[I, 1], data[I,
                                    2], pch = symFill[
                                          k], cex = CEX)
          }
        }
      }
      if(identify)
        title("Classification Errors", cex = 0.5)
    }
    if(!ask)
      break
  }
  invisible()
}

"mclustDA" <- function(trainingData, labels, testData, G = 1:6,
                       verbose = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(verbose)
    cat("training ...\n")
  emModelNames <- c("EII", "VII", "EEI", "VVI", "EEE", "VVV")
  trainingModels <- mclustDAtrain(data = trainingData, labels = labels,
                                  G = G, emModelNames = emModelNames,
                                  verbose = verbose)
  if(verbose)
    cat("testing ...\n")
  S <- data.frame(trainClass = as.factor(unique(labels)),
                  mclustModel = as.factor(sapply(trainingModels,
                    function(x) x$modelName)),
                  numGroups = sapply(trainingModels, function(x) x$G))
  test <- mclustDAtest(testData, trainingModels)
  testSumry <- summary(test)
  train <- mclustDAtest(trainingData, trainingModels)
  trainSumry <- summary(train)
  structure(list(testClassification = testSumry$classification, 
                 trainingClassification = trainSumry$classification,
                 summary = S,
                 VofIindex = compareClass(map(trainSumry$z), labels),
                 models = trainingModels, postProb = testSumry$postProb),
            class = "mclustDA")
}

"mclustDAtest" <- function(data, models)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  densfun <- function(model, data)
    {
      do.call("dens", c(list(data = data), model))
    }
  den <- as.matrix(data.frame(lapply(models, densfun, data = data)))
  dimnames(den) <- list(NULL, names(models))
  structure(den, class = "mclustDAtest")
}

"mclustDAtrain" <- function(data, labels, G, emModelNames, eps, tol,
                            itmax, equalPro, warnSingular = FALSE,
                            verbose = TRUE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  if(missing(itmax))
    itmax <- .Mclust$itmax
  
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(oneD) {
    data <- as.vector(data)
    n <- length(data)
    p <- 1
    data <- as.matrix(data)
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
  }
  if(missing(emModelNames)) {
    if(p == 1) {
      emModelNames <- c("E", "V")
    }
    else {
      emModelNames <- .Mclust$emModelNames
    }
  }
  if(p == 1 && any(nchar(emModelNames) > 1)) {
    Emodel <- any(sapply(emModelNames, function(x)
                         charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(emModelNames, function(x)
                         charmatch("V", x, nomatch = 0)[1]) == 1)
    emModelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  ## Old: U <- unique(labels)
  ## New... 10/01/02 Ron
  if (!is.factor(labels)) labels <- as.factor(labels)
  U <- levels(labels)
  
  L <- length(U)
  S <- rep(0, L)
  M <- rep("XXX", L)
  if(missing(G)) {
    G <- 1:9
  }
  else {
    G <- sort(G)
  }
  if(any(G) <= 0)
    stop("G must be positive")
  R <- rep(list(matrix(0, n, max(G) + 1)), L)
  Glabels <- as.character(G)
  for(l in 1:L) {
    I <- labels == U[l]
    BIC <- EMclust(data[I,  ], emModelNames = emModelNames, G = G)
    SUMMARY <- summary(BIC, data[I,  ])
    S[l] <- SUMMARY$G
    M[l] <- SUMMARY$modelName
    R[[l]] <- SUMMARY
  }
  names(S) <- M
  if(verbose) print(S)
  names(R) <- U
  R <- lapply(R, function(x)
              {
		i <- charmatch("Vinv", names(x), nomatch = 0)
		if(i) x[ - i]
		else x
              }
              )
  structure(R, G = G, emModelNames = emModelNames, eps = eps, tol = tol,
            itmax = itmax, equalPro = equalPro, class = "mclustDAtrain")
}

"mclustOptions" <- function(eps = .Machine$double.eps,
                            tol = c(1.0000000000000001e-05, 
                              1.0000000000000001e-05),
                            itmax = c(Inf, Inf), equalPro = FALSE, 
                            warnSingular = TRUE, emModelNames,
                            hcModelName = c("V", "VVV"),
                            symbols) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(emModelNames))
    emModelNames <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI",
                      "EEE", "EEV", "VEV", "VVV")
  if(missing(symbols))
    symbols <- c(17, 0, 10, 4, 11, 18, 6, 7, 3, 16, 2, 12, 8, 15,
                 1, 9, 14, 13, 5)
  list(tol = tol, eps = eps, itmax = itmax, equalPro = equalPro, 
       emModelNames = emModelNames, hcModelName = hcModelName, 
       warnSingular = warnSingular, symbols = symbols)
}

"me" <- function(modelName, data, z, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## ... z, eps, tol, itmax, equal = FALSE, noise = FALSE, Vinv
  funcName <- paste("me", modelName, sep = "")
  do.call(funcName, list(data = data, z = z, ...))
}

"meE" <- function(data, z, eps, tol, itmax, equalPro, warnSingular,
                  noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD)
    stop("data must be 1 dimensional")
  data <- as.vector(data)
  n <- length(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal length of data")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = 1, G = G, z = z, mu = rep(NA, G),
                          sigmasq = NA, pro = rep(NA, K), loglik = NA,
                          modelName = "E"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("me1e",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(G),
                   double(1),
                   double(K),
                   PACKAGE="mclust")[6:12]
  mu <- temp[[5]]
  names(mu) <- as.character(1:G)
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  warn <- NULL
  if(is.infinite(loglik) || sigmasq <= max(eps, 0)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
  }
  else if(its >= itmax) {
    warning("iteration limit reached")
    warn <- "iteration limit reached"
    its <-  - its
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = 1, G = G, z = z, mu = mu, sigmasq = sigmasq,
                 pro = pro, loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "E"), info = info, warn = warn)
}

"meEEE" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          Sigma = matrix(NA, p, p),
                          cholSigma = matrix(NA, p, p), pro = rep(NA, K),
                          loglik = NA, modelName = "EEE"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  storage.mode(z) <- "double"
  temp <- .Fortran("meeee",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(p * p),
                   double(K),
                   double(p),
                   PACKAGE="mclust")[7:13]
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholSigma <- structure(matrix(temp[[6]], p, p), def = 
                         "Sigma = t(cholSigma) %*% cholSigma")
  Sigma <- unchol(cholSigma, upper = TRUE)
  pro <- temp[[7]]
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 Sigma = Sigma, cholSigma = cholSigma, pro = pro, loglik = 
                 loglik, Vinv = if(noise) Vinv else NULL, modelName = "EEE"),
            info = info, warn = warn)
}

"meEEI" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data  should be in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          Sigma = matrix(NA, p, p),
                          decomp = list(d = p, G = G,
                            scale = NA, shape = rep(NA, p)),
                          pro = rep(NA, K), loglik = NA,
                          modelName = "EEI"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("meeei",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(1),
                   double(p),
                   double(K),
                   PACKAGE="mclust")[7:14]
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[6]]
  shape <- temp[[7]]
  pro <- temp[[8]]
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(NA, p, p)
    mu[] <- pro[] <- z[] <- loglik <- NA
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(scale * shape)
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 Sigma = Sigma, decomp = decomp, pro = pro, loglik = loglik,
                 Vinv = if(noise) Vinv else NULL, modelName = "EEI"), info = 
            info, warn = warn)
}

"meEEV" <- function(data, z, eps, tol, itmax, equalPro, warnSingular,
                    noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = NA,
                            shape = rep(NA, p),
                            orientation = array(NA, c(p, p, G))),
                          pro = rep(NA, K), loglik = NA, modelName = "EEV"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p))
  storage.mode(z) <- "double"
  temp <- .Fortran("meeev",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   as.integer(lwork),
                   double(p * G),
                   double(1),
                   double(p),
                   double(p * p * G),
                   double(K),
                   double(lwork),
                   double(p),
                   PACKAGE="mclust")[7:16]
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  lapackSVDinfo <- temp[[5]]
  mu <- matrix(temp[[6]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[7]]
  shape <- temp[[8]]
  O <- array(temp[[9]], c(p, p, G))
  pro <- temp[[10]]
  warn <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      warning("LAPACK DGESVD fails to converge")
      warn <- "LAPACK DGESVD fails to converge"
    }
    else {
      warning("input error for LAPACK DGESVD")
      warn <- "input error for LAPACK DGESVD"
    }
    z[] <- O[] <- shape[] <- NA
    scale <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- scale * shapeO(shape, O, transpose = TRUE)
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  decomp <- structure(list(d = p, G = G, scale = scale, shape = shape,
                           orientation = O), def = 
                      "Sigma = scale * t(O) %*% diag(shape) %*% O")
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 decomp = decomp, pro = pro, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "EEV"), info = info, warn
            = warn)
}

"meEII" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)), sigmasq = NA,
                          Sigma = matrix(NA, p, p),
                          decomp = list(d = p, G = G, scale = NA),
                          pro = rep(NA, K), loglik = NA, modelName = "EII"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("meeii",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(1),
                   double(K),
                   PACKAGE="mclust")[7:13]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  Sigma <- diag(rep(sigmasq, p))
  pro <- temp[[7]]
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax ||
     sigmasq <= max(eps, 0)) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = sigmasq)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 sigmasq = sigmasq, Sigma = Sigma, decomp = decomp, pro = pro,
                 loglik = loglik, Vinv = if(noise) Vinv else NULL, info = info,
                 modelName = "EII"), warn = warn)
}

"meEVI" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = NA,
                            shape = matrix(NA, p, G)), pro = rep(NA, K),
                          loglik = NA, modelName = "EVI"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("meevi",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(1),
                   double(p * G),
                   double(K),
                   PACKAGE="mclust")[7:14]
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 decomp = decomp, pro = pro, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "EVI"), info = info, warn
            = warn)
}

"meV" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal length of data")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warning("z is missing")
    warn <- "z is missing"
    return(structure(list(n = n, d = 1, G = G, z = z, mu = rep(NA, G),
                          sigmasq = rep(NA, G), pro = rep(NA, K), 
                          modelName = "V"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("me1v",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(G),
                   double(G),
                   double(K),
                   PACKAGE="mclust")[6:12]
  mu <- temp[[5]]
  names(mu) <- as.character(1:G)
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  warn <- NULL
  if(is.infinite(loglik) || any(sigmasq <= max(eps, 0))) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
  }
  else if(its >= itmax) {
    warning("iteration limit reached")
    warn <- "iteration limit reached"
    its <-  - its
  }
  info <- c(iterations = its, error = err)
  if(is.infinite(loglik))
    loglik <- NA
  structure(list(n = n, d = 1, G = G, z = z, mu = mu, sigmasq = sigmasq,
                 pro = pro, loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "V"), info = info, warn = warn)
}

"meVEI" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G,
                            scale = rep(NA, G), shape = rep(NA, p)),
                          pro = rep(NA, K), loglik = NA, modelName = "VEI"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  if(length(tol) == 1)
    tol <- c(tol, tol)
  if(missing(itmax))
    itmax <- .Mclust$itmax
  if(length(itmax) == 1)
    itmax <- c(itmax, Inf)
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("mevei",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(G),
                   double(p),
                   double(K),
                   double(G),
                   double(p),
                   double(p * G),
                   PACKAGE="mclust")[7:14]
  z <- temp[[1]]
  its <- temp[[2]][1]
  inner <- temp[[2]][2]
  err <- temp[[3]][1]
  inerr <- temp[[3]][2]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- temp[[7]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- NA
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= itmax[2]) {
      warning("inner iteration limit reached")
      warn <- "inner iteration limit reached"
      inner <-  - inner
    }
    else if(its >= itmax[1]) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  attr(info, "inner") <- c(iterations = inner, error = inerr)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 decomp = decomp, pro = pro, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "VEI"), info = info, warn
            = warn)
}

"meVEV" <- function(data, z, eps, tol, itmax, equalPro, warnSingular,
                    noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z, mu = matrix(
                                                    NA, p, G), sigma = array(NA, c(p, p, G)), decomp = list(
                                                                                                d = p, G = G, scale = rep(NA, G), shape = matrix(NA,
                                                                                                                                    p, G), orientation = array(NA, c(p, p, G))), pro = rep(
                                                                                                                                                                                   NA, K), modelName = "VEV"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  if(length(tol) == 1)
    tol <- c(tol, tol)
  if(missing(itmax))
    itmax <- .Mclust$itmax
  if(length(itmax) == 1)
    itmax <- c(itmax, Inf)
  itmax[is.infinite(itmax)] <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  storage.mode(z) <- "double"
  temp <- .Fortran("mevev",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   as.integer(lwork),
                   double(p * G),
                   double(G),
                   double(p),
                   double(p * p * G),
                   double(K),
                   double(lwork),
                   double(p),
                   PACKAGE="mclust")[7:16]
  z <- temp[[1]]
  its <- temp[[2]][1]
  inner <- temp[[2]][2]
  err <- temp[[3]][1]
  inerr <- temp[[3]][2]
  loglik <- temp[[4]]
  lapackSVDinfo <- temp[[5]]
  mu <- matrix(temp[[6]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[7]]
  shape <- temp[[8]]
  O <- array(temp[[9]], c(p, p, G))
  pro <- temp[[10]]
  warn <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      warning("LAPACK DGESVD fails to converge")
      warn <- "LAPACK DGESVD fails to converge"
    }
    else {
      warning("input error for LAPACK DGESVD")
      warn <- "input error for LAPACK DGESVD"
    }
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    Sigma <- array(NA, c(p, p, G))
  }
  else if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    Sigma <- array(NA, c(p, p, G))
  }
  else {
    Sigma <- shapeO(shape, O, transpose = TRUE)
    Sigma <- sweep(Sigma, MARGIN = 3, STATS = scale, FUN = "*")
    if(inner >= itmax[2]) {
      warning("inner iteration limit reached")
      warn <- "inner iteration limit reached"
      inner <-  - inner
    }
    else if(its >= itmax[1]) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  decomp <- structure(list(d = p, G = G, scale = scale, shape = shape,
                           orientation = O), def = 
                      "Sigma = scale * t(O) %*% diag(shape) %*% O")
  info <- structure(c(iterations = its, error = err), inner = c(
                                                        iterations = inner, error = inerr))
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = Sigma,
                 decomp = decomp, pro = pro, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "VEV"), info = info, warn
            = warn)
}

"meVII" <- function(data, z, eps, tol, itmax, equalPro, warnSingular,
                    noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data must be in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          sigmasq = rep(NA, G),
                          decomp = list(d = p, G = G, scale = rep(NA, G)),
                          pro = rep(NA, K), loglik = NA, modelName = "VII"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("mevii",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(G),
                   double(K),
                   PACKAGE="mclust")[7:13]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax ||
     any(sigmasq <= max(eps, 0))) {
    if(warnSingular)
      warning("sigma-squared falls below threshold")
    warn <- "sigma-squared falls below threshold"
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(rep(sigmasq[k], p))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = sigmasq)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 sigmasq = sigmasq, decomp = decomp, pro = pro,
                 loglik = loglik, Vinv = if(noise) Vinv else NULL,
                 modelName = "VII"), info = info, warn = warn) 
}

"meVVI" <- 
  function(data, z, eps, tol, itmax, equalPro, warnSingular, noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = rep(NA, G),
                            shape = matrix(NA, p, G)), pro = rep(NA, K),
                          loglik = NA, modelName = "VVI"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("mevvi",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(G),
                   double(p * G),
                   double(K),
                   PACKAGE="mclust")[7:14]
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  warn <- NULL
  if(is.infinite(loglik) || abs(loglik) == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    if(loglik < 0)
      shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- NA
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
                               FUN = "*"), 2, diag), c(p, p, G))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 decomp = decomp, pro = pro, loglik = loglik, Vinv = if(noise) 
                 Vinv else NULL, modelName = "VVI"), info = info, warn
            = warn)
}

"meVVV" <- function(data, z, eps, tol, itmax, equalPro, warnSingular,
                    noise = FALSE, Vinv)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
                                        # number of groups
  if(!noise) {
    G <- K
    Vinv <- -1
  }
  else {
    G <- K - 1
    if(missing(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, z = z,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          cholsigma = array(NA, c(p, p, G)),
                          pro = rep(NA, K),
                          loglik = NA, modelName = "VVV"), warn = warn)) 
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- .Mclust$tol
  tol <- tol[1]
  if(missing(itmax))
    itmax <- .Mclust$itmax
  itmax <- itmax[1]
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  storage.mode(z) <- "double"
  temp <- .Fortran("mevvv",
                   as.logical(equalPro),
                   as.double(data),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(Vinv),
                   z,
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(p * p * G),
                   double(K),
                   double(p),
                   PACKAGE="mclust")[7:13]
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- structure(array(temp[[6]], c(p, p, G)), def = 
                         "Sigma = t(cholsigma) %*% cholsigma")
  pro <- temp[[7]]
  warn <- NULL
  if(is.infinite(loglik) || loglik == .Machine$double.xmax) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(cholsigma, 3, unchol, TRUE), c(p, p, G))
    if(its >= itmax) {
      warning("iteration limit reached")
      warn <- "iteration limit reached"
      its <-  - its
    }
  }
  info <- c(iterations = its, error = err)
  structure(list(n = n, d = p, G = G, z = z, mu = mu, sigma = sigma,
                 cholsigma = cholsigma, pro = pro, loglik = loglik,
                 Vinv = if(noise) Vinv else NULL, modelName = "VVV"),
            info = info, warn = warn)
}

"mstep" <- function(modelName, data, z, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## ... eps, tol, itmax, equal = FALSE, noise = FALSE
  funcName <- paste("mstep", modelName, sep = "")
  do.call(funcName, list(data = data, z = z, ...))
}

"mstepE" <- function(data, z, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.vector(data)
  n <- length(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = 1, G = G, mu = rep(NA, G),
                          sigmasq = NA, pro = rep(NA, K), modelName = "E"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  temp <- .Fortran("ms1e",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(G),
                   double(G),
                   double(1),
                   double(G),
                   PACKAGE="mclust")[5:7]
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[3]]
    }
    else {
      pro <- c(temp[[3]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  structure(list(n = n, d = 1, G = G, mu = mu, sigmasq = sigmasq, pro = 
                 pro, modelName = "E"))
}

"mstepEEE" <- function(data, z, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          Sigma = matrix(NA, p, p),
                          cholSigma = matrix(NA, p, p),
                          pro = rep(NA, K), modelName = "EEE"),
                     warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  temp <- .Fortran("mseee",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p),
                   double(p * G),
                   double(p * p),
                   double(G),
                   PACKAGE="mclust")[7:9]
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholSigma <- structure(matrix(temp[[2]], p, p), def = 
                         "Sigma = t(cholSigma) %*% cholSigma")
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[3]]
    }
    else {
      pro <- c(temp[[3]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  sigma <- array(0, c(p, p, G))
  Sigma <- unchol(cholSigma, upper = TRUE)
  for(k in 1:G)
    sigma[,  , k] <- Sigma
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 Sigma = Sigma, cholSigma = cholSigma, pro = pro,
                 modelName = "EEE"))
}

"mstepEEI" <- function(data, z, equalPro, noise = FALSE, eps,
                       warnSingular, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          Sigma = matrix(NA, p, p),
                          decomp = list(d = p, G = G, scale = NA,
                            shape = rep(NA, p)), pro = rep(NA, K),
                          modelName = "EEI"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("mseei",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(p * G),
                   double(1),
                   double(p),
                   double(G),
                   PACKAGE="mclust")[6:10]
  icond <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[3]]
  shape <- temp[[4]]
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[5]]
    }
    else {
      pro <- c(temp[[5]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  warn <- NULL
  if(icond <= eps) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- scale <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(scale * shape)
    for(k in 1:G)
      sigma[,  , k] <- Sigma
  }
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, Sigma = 
                 Sigma, decomp = decomp, pro = pro, modelName = "EEI"), warn = 
            warn)
}

"mstepEEV" <- function(data, z, equalPro, noise = FALSE, eps,
                       warnSingular, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = NA,
                            shape = rep(NA, p),
                            orientation = array(NA, c(p, p, G))),
                          pro = rep(NA, K), modelName = "EEV"),
                     warn = warn))
  }
  ##	shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), G)
  temp <- .Fortran("mseev",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(lwork),
                   as.integer(lwork),
                   as.double(eps),
                   double(p * G),
                   double(1),
                   double(p),
                   double(p * p * G),
                   double(G),
                   PACKAGE="mclust")[7:13]
  lapackSVDinfo <- temp[[1]]
  smin <- temp[[2]]
  mu <- matrix(temp[[3]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[4]]
  shape <- temp[[5]]
  O <- array(temp[[6]], c(p, p, G))
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[7]]
    }
    else {
      pro <- c(temp[[7]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  warn <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      warning("LAPACK DGESVD fails to converge")
      warn <- "LAPACK DGESVD fails to converge"
    }
    else {
      warning("input error for LAPACK DGESVD")
      warn <- "input error for LAPACK DGESVD"
    }
    O[] <- shape[] <- scale <- NA
    Sigma <- array(NA, c(p, p, G))
  }
  else if(smin <= eps) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- scale <- O[] <- shape[] <- NA
    Sigma <- array(NA, c(p, p, G))
  }
  else {
    Sigma <- scale * shapeO(shape, O, transpose = TRUE)
  }
  decomp <- structure(list(d = p, G = G, scale = scale, shape = shape,
                           orientation = O), def = 
                      "Sigma = scale * t(O) %*% diag(shape) %*% O")
  structure(list(n = n, d = p, G = G, mu = mu, sigma = Sigma, decomp = 
                 decomp, pro = pro, modelName = "EEV"), warn = warn)
}

"mstepEII" <- function(data, z, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)), sigmasq = NA,
                          Sigma = matrix(NA, p, p),
                          decomp = list(d = p, G = G, scale = NA),
                          pro = rep(NA, K), modelName = "EII"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  temp <- .Fortran("mseii",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p * G),
                   double(1),
                   double(G),
                   PACKAGE="mclust")[6:8]
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[3]]
    }
    else {
      pro <- c(temp[[3]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  sigma <- array(0, c(p, p, G))
  Sigma <- diag(rep(sigmasq, p))
  for(k in 1:G)
    sigma[,  , k] <- Sigma
  decomp <- list(d = p, G = G, scale = sigmasq)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, sigmasq = 
                 sigmasq, Sigma = Sigma, decomp = decomp, pro = pro, modelName
		 = "EII"))
}

"mstepEVI" <- function(data, z, equalPro, noise = FALSE, eps, warnSingular, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = NA,
                            shape = matrix(NA, p, G)), pro = rep(NA, K),
                          modelName = "EVI"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("msevi",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(p * G),
                   double(1),
                   double(p * G),
                   double(G),
                   PACKAGE="mclust")[6:10]
  icond <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  scale <- temp[[3]]
  shape <- matrix(temp[[4]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[5]]
    }
    else {
      pro <- c(temp[[5]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  warn <- NULL
  if(icond <= eps) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- scale <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
  }
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, decomp = 
                 decomp, pro = pro, modelName = "EVI"), warn = warn)
}

"mstepV" <- function(data, z, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.vector(data)
  n <- length(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = 1, G = G, mu = rep(NA, G),
                          sigmasq = rep(NA, G), pro = rep(NA, K),
                          modelName = "V"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  temp <- .Fortran("ms1v",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(G),
                   double(G),
                   double(G),
                   double(G),
                   PACKAGE="mclust")[5:7]
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[3]]
    }
    else {
      pro <- c(temp[[3]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  structure(list(n = n, d = 1, G = G, mu = mu, sigmasq = sigmasq, pro = 
                 pro, modelName = "V"))
}

"mstepVEI" <- function(data, z, equalPro, noise = FALSE, eps, tol, itmax,
                       warnSingular, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = rep(NA, p),
                            shape = matrix(NA, p, G)), pro = rep(NA, K),
                          modelName = "VEI"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- if(length(.Mclust$tol) > 1) .Mclust$tol[2] else .Mclust$
  tol
  if(missing(itmax))
    itmax <- if(length(.Mclust$itmax) > 1) .Mclust$itmax[2] else 
  Inf
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("msvei",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(G),
                   double(p),
                   double(G),
                   double(G),
                   double(p),
                   double(p * G),
                   PACKAGE="mclust")[6:12]
  inner <- temp[[1]]
  inerr <- temp[[2]]
  icond <- temp[[3]]
  mu <- matrix(temp[[4]], p, G)
  scale <- temp[[5]]
  shape <- temp[[6]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[7]]
    }
    else {
      pro <- c(temp[[7]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  warn <- NULL
  if(icond <= eps) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- shape <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= itmax) {
      warning("inner iteration limit reached")
      warn <- "inner iteration limit reached"
      inner <-  - inner
    }
  }
  info <- c(iterations = inner, error = inerr)
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma,
                 decomp = decomp, pro = pro, modelName = "VEI"),
            info = info, warn = warn)
}

"mstepVEV" <- function(data, z, equalPro, noise = FALSE, eps, tol,
                       itmax, warnSingular, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G,
                          mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = rep(NA, G),
                            shape = rep(NA, p),
                            orientation = array(NA, c(p, p, G))),
                          pro = rep(NA, K), modelName = "VEV"),
                     warn = warn))
  }
  ##	shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(tol))
    tol <- if(length(.Mclust$tol) > 1) .Mclust$tol[2] else .Mclust$tol
  if(missing(itmax))
    itmax <- if(length(.Mclust$itmax) > 1) .Mclust$itmax[2] else 
  .Mclust$itmax
  if(is.infinite(itmax))
    itmax <- .Machine$integer.max
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  temp <- .Fortran("msvev",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(lwork),
                   as.integer(lwork),
                   as.integer(itmax),
                   as.double(tol),
                   as.double(eps),
                   double(p * G),
                   double(G),
                   double(p),
                   double(p * p * G),
                   double(G),
                   PACKAGE="mclust")[7:15]
  lapackSVDinfo <- temp[[1]]
  inner <- temp[[2]]
  inerr <- temp[[3]]
  smin <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[6]]
  shape <- temp[[7]]
  O <- array(temp[[8]], c(p, p, G))
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[9]]
    }
    else {
      pro <- c(temp[[9]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  warn <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      warning("LAPACK DGESVD fails to converge")
      warn <- "LAPACK DGESVD fails to converge"
    }
    else {
      warning("input error for LAPACK DGESVD")
      warn <- "input error for LAPACK DGESVD"
    }
    O[] <- shape[] <- scale[] <- NA
    Sigma <- array(NA, c(p, p, G))
  }
  else if(smin <= eps) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
    Sigma <- array(NA, c(p, p, G))
  }
  else {
    Sigma <- sweep(shapeO(shape, O, transpose = TRUE), MARGIN = 3,
                   STATS = scale, FUN = "*")
    if(inner >= itmax) {
      warning("inner iteration limit reached")
      warn <- "inner iteration limit reached"
      inner <-  - inner
    }
  }
  decomp <- structure(list(d = p, G = G, scale = scale, shape = shape,
                           orientation = O), def = 
                      "Sigma = scale * t(O) %*% diag(shape) %*% O")
  info <- c(iteration = inner, error = inerr)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = Sigma, decomp = 
                 decomp, pro = pro, modelName = "VEV"), info = info,
            warn = warn)
}

"mstepVII" <- function(data, z, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          sigmasq = rep(NA, G),
                          decomp = list(d = p, G = G, scale = rep(NA, G)),
                          pro = rep(NA, K), modelName = "VII"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  temp <- .Fortran("msvii",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p * G),
                   double(G),
                   double(G),
                   PACKAGE="mclust")[6:8]
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[3]]
    }
    else {
      pro <- c(temp[[3]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  sigma <- array(0, c(p, p, G))
  for(k in 1:G)
    sigma[,  , k] <- diag(rep(sigmasq[k], p))
  decomp <- list(d = p, G = G, scale = sigmasq)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, sigmasq = 
                 sigmasq, decomp = decomp, pro = pro, modelName = "VII"))
}

"mstepVVI" <- function(data, z, equalPro, noise = FALSE, eps,
                       warnSingular, ...) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          decomp = list(d = p, G = G, scale = rep(NA, p),
                            shape = matrix(NA, p, G)), pro = rep(NA, K),
                          modelName = "VVI"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Mclust$eps
  if(missing(warnSingular))
    warnSingular <- .Mclust$warnSingular
  temp <- .Fortran("msvvi",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(eps),
                   double(p * G),
                   double(G),
                   double(p * G),
                   double(G),
                   PACKAGE="mclust")[6:10]
  icond <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  scale <- temp[[3]]
  shape <- matrix(temp[[4]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[5]]
    }
    else {
      pro <- c(temp[[5]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  warn <- NULL
  if(icond <= eps) {
    if(warnSingular)
      warning("singular covariance")
    warn <- "singular covariance"
    mu[] <- pro[] <- shape <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
                               FUN = "*"), 2, diag), c(p, p, G))
  }
  decomp <- list(d = p, G = G, scale = scale, shape = shape)
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, decomp = 
                 decomp, pro = pro, modelName = "VVI"), warn = warn)
}

"mstepVVV" <- function(data, z, equalPro, noise = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  if(missing(equalPro))
    equalPro <- .Mclust$equalPro
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
                                        # number of groups
  K <- dimz[2]
  G <- if(noise) K - 1 else K
  ##
  if(all(is.na(z))) {
    warn <- "z is missing"
    warning("z is missing")
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA, p, G),
                          sigma = array(NA, c(p, p, G)),
                          cholsigma = array(NA, c(p, p, G)),
                          pro = rep(NA, K), modelName = "VVV"), warn = warn))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  temp <- .Fortran("msvvv",
                   as.double(data),
                   as.double(z),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   double(p),
                   double(p * G),
                   double(p * p * G),
                   double(G),
                   PACKAGE="mclust")[7:9]
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- structure(array(temp[[2]], c(p, p, G)), def = 
                         "Sigma = t(cholsigma) %*% cholsigma")
  sigma <- array(apply(cholsigma, 3, unchol, TRUE), c(p, p, G))
  if(!equalPro) {
    if(!noise) {
      pro <- temp[[3]]
    }
    else {
      pro <- c(temp[[3]], sum(z[, K])/n)
    }
  }
  else {
    if(!noise) {
      pro <- rep(1/G, G)
    }
    else {
      pron <- sum(z[, K])/n
      pro <- c(rep((1 - pron)/G, G), pron)
    }
  }
  structure(list(n = n, d = p, G = G, mu = mu, sigma = sigma, cholsigma
		 = cholsigma, pro = pro, modelName = "VVV"))
}

"mvn" <- function(modelName, data)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  switch(EXPR=as.character(modelName),
         E = , V = ,
         X = mvnX(data),
	 Spherical = , EII = , VII = ,
         XII = mvnXII(data),
	 Diagonal = , EEI = , VEI = , EVI = , VVI = ,
         XXI = mvnXXI(data),
	 Ellipsoidal = , EEE = , VEE = , EVE = , VVE = , EEV = , VEV = ,
         VVV = , 
         XXX = mvnXXX(data),
	 stop("invalid model name"))
}

"mvn2plot" <- function(mu, sigma, k = 15., alone = FALSE, col=1)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  p <- length(mu)
  if(p != 2.)
    stop("two-dimensional case only")
  if(any(unique(dim(sigma)) != p))
    stop("mu and sigma are incompatible")
  ev <- eigen(sigma, symmetric = TRUE)
  s <- sqrt(rev(sort(ev$values)))
  V <- t(ev$vectors[, rev(order(ev$values))])
  theta <- (0.:k) * (pi/(2. * k))
  x <- s[1.] * cos(theta)
  y <- s[2.] * sin(theta)
  xy <- cbind(c(x,  - x,  - x, x), c(y, y,  - y,  - y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2., STATS = mu, FUN = "+")
  if(alone) {
    xymin <- apply(xy, 2., FUN = "min")
    xymax <- apply(xy, 2., FUN = "max")
    r <- ceiling(max(xymax - xymin)/2.)
    xymid <- (xymin + xymax)/2.
    plot(xy[, 1.], xy[, 2.], xlim = c( - r, r) + xymid[1.], ylim = 
         c( - r, r) + xymid[2.], xlab = "x", ylab = "y", type = 
         "n")
  }
  l <- length(x)
  i <- 1.:l
  for(k in 1.:4.) {
    lines(xy[i,  ], col=col)
    i <- i + l
  }
                                        # semi-major axes
  x <- s[1.]
  y <- s[2.]
  xy <- cbind(c(x,  - x, 0, 0), c(0, 0, y,  - y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2., STATS = mu, FUN = "+")
  lines(xy[1:2,  ], lty = 2., col=col)
  lines(xy[3:4,  ], lty = 2., col=col)
  points(mu[1.], mu[2.], pch = 8, col=col)
  invisible()
}

"mvnX" <- function(data)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one dimensional")
  data <- as.vector(data)
  n <- length(data)
  temp <- .Fortran("mvn1d",
                   as.double(data),
                   as.integer(n),
                   double(1),
                   double(1),
                   double(1),
                   PACKAGE="mclust")[3:5]
  mu <- temp[[1]]
  sigmasq <- temp[[2]]
  loglik <- temp[[3]]
  list(n = n, d = 1, G = 1, mu = mu, sigmasq = sigmasq, loglik = loglik,
       modelName = "X")
}

"mvnXII" <- function(data)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD) {
    data <- as.vector(data)
    n <- length(data)
    temp <- .Fortran("mvn1d",
                     as.double(data),
                     as.integer(n),
                     double(1),
                     double(1),
                     double(1),
                   PACKAGE="mclust")[3:5]
    mu <- temp[[1]]
    sigmasq <- temp[[2]]
    loglik <- temp[[3]]
    list(n = n, d = 1, G = 1, mu = mu, sigmasq = sigmasq, loglik = 
         loglik, modelName = "X")
  }
  else {
    if(length(dimdat) != 2)
      stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    temp <- .Fortran("mvnxii",
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     double(p),
                     double(1),
                     double(1),
                   PACKAGE="mclust")[4:6]
    mu <- temp[[1]]
    sigmasq <- temp[[2]]
    loglik <- temp[[3]]
    Sigma <- sigmasq * diag(p)
    decomp <- list(d = p, G = 1, scale = sigmasq)
    list(n = n, d = p, G = 1, mu = matrix(mu, ncol = 1), sigmasq = 
         sigmasq, sigma = array(Sigma, c(p, p, 1)), Sigma = 
         Sigma, decomp = decomp, loglik = loglik, modelName = 
         "XII")
  }
}

"mvnXXI" <- function(data)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD) {
    data <- as.vector(data)
    n <- length(data)
    temp <- .Fortran("mvn1d",
                     as.double(data),
                     as.integer(n),
                     double(1),
                     double(1),
                     double(1),
                   PACKAGE="mclust")[3:5]
    mu <- temp[[1]]
    sigmasq <- temp[[2]]
    loglik <- temp[[3]]
    list(n = n, d = 1, G = 1, mu = mu, sigmasq = sigmasq, loglik = 
         loglik, modelName = "X")
  }
  else {
    if(length(dimdat) != 2)
      stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    temp <- .Fortran("mvnxxi",
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     double(p),
                     double(1),
                     double(p),
                     double(1),
                   PACKAGE="mclust")[4:7]
    mu <- temp[[1]]
    scale <- temp[[2]]
    shape <- temp[[3]]
    loglik <- temp[[4]]
    Sigma <- diag(scale * shape)
    decomp <- list(p = p, G = 1, scale = scale, shape = shape)
    list(n = n, d = p, G = 1, mu = matrix(mu, ncol = 1), sigma = 
         array(Sigma, c(p, p, 1)), Sigma = Sigma, decomp = 
         decomp, loglik = loglik, modelName = "XXI")
  }
}

"mvnXXX" <- function(data)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD) {
    data <- as.vector(data)
    n <- length(data)
    temp <- .Fortran("mvn1d",
                     as.double(data),
                     as.integer(n),
                     double(1),
                     double(1),
                     double(1),
                   PACKAGE="mclust")[c(3:5)]
    mu <- temp[[1]]
    sigmasq <- temp[[2]]
    loglik <- temp[[3]]
    list(n = n, d = 1, G = 1, mu = mu, sigmasq = sigmasq, loglik = 
         loglik, modelName = "X")
  }
  else {
    if(length(dimdat) != 2)
      stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    temp <- .Fortran("mvnxxx",
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     double(p),
                     double(p * p),
                     double(1),
                   PACKAGE="mclust")[c(4:6)]
    mu <- temp[[1]]
    chol <- matrix(temp[[2]], p, p)
    Sigma <- unchol(chol, upper = TRUE)
    loglik <- temp[[3]]
    list(n = n, d = p, G = 1, mu = matrix(mu, ncol = 1), sigma = 
         array(Sigma, c(p, p,  , 1)), Sigma = Sigma, cholSigma
         = structure(chol, def = "Sigma = t(chol) %*% chol"),
         loglik = loglik, modelName = "XXX")
  }
}

"orth2" <- function(n)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  u <- rnorm(n)
  u <- u/vecnorm(u)
  v <- rnorm(n)
  v <- v/vecnorm(v)
  Q <- cbind(u, v - sum(u * v) * u)
  dimnames(Q) <- NULL
  Q
}

"partconv" <- function(x, consec = TRUE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  n <- length(x)
  y <- numeric(n)
  u <- unique(x)
  if(consec) {
    ## number groups in order of first row appearance
    l <- length(u)
    for(i in 1:l)
      y[x == u[i]] <- i
  }
  else {
    ## represent each group by its lowest-numbered member
    for(i in u) {
      l <- x == i
      y[l] <- (1:n)[l][1]
    }
  }
  y
}

"charconv" <- function(x, sep = "001")
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(!is.data.frame(x))
    x <- data.frame(x)
  do.call("paste", c(as.list(x), sep = sep))
}

"partuniq" <- function(x)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  n <- nrow(x)
  x <- charconv(x)
  k <- duplicated(x)
  partition <- 1.:n
  partition[k] <- match(x[k], x)
  partition
}

### dots are not used but inserted to keep R from complaining
### about generic/method consistency...
### Somewhat simpler function using matplot... inserted 10/01/02,Ron
"plot.EMclust" <- function(x, modelNames, G, symbols,
                           xlab = "number of clusters",
                           ylab = "BIC", ...)
{ 
  n <- ncol(x)
  dnx <- dimnames(x)

  x <- matrix(as.vector(x), ncol = n)
  dimnames(x) <- dnx
  if(missing(modelNames))
    modelNames <- dimnames(x)[[2]]
  if(missing(G))
    G <- dimnames(x)[[1]]
  else G <- as.character(G)
  BIC <- x[G, modelNames, drop = FALSE]
  X <- is.na(BIC)
  nrowBIC <- nrow(BIC)
  ncolBIC <- ncol(BIC)
  if(missing(symbols)) {
     if(n > 9) {
       symbols <- LETTERS[1:n]
     }
     else {
       symbols <- as.character(1:n)
     }
     names(symbols) <- dnx[[2]]
     symbols <- symbols[modelNames]
   }
   xrange <-
     if(!is.null(dn <- dimnames(BIC)[[1]])) as.numeric(dn) else 1:nrowBIC
##   plot(xrange, BIC[, 1], type = "n", ylim = range(as.vector(BIC[!X])),
##        xlim = range(xrange), xlab = "number of clusters", ylab = "BIC"
##        )
##   for(i in 1:ncolBIC) {
##     x <- !X[, i]
##     if(any(x)) {
##       if(all(x)) {
##         points(xrange, BIC[, i], pch = symbols[i])
##         if(nrowBIC > 1)
##           lines(xrange, BIC[, i], lty = i)
##       }
##       else {
##         points(xrange[x], BIC[x, i], pch = symbols[
##                                        i])
##         v <- (1:nrowBIC)[x]
##         n <- length(v)
##         d <- c(diff(c(0, v)), 2)
##         if(any(d) == 1) {
##           f <- 1
##           while(TRUE) {
##             while(f <= n && d[f] != 1) f <-
##               f + 1
##             if(f > n)
##               break
##             l <- f
##             while(d[l + 1] == 1) l <- l +
##               1
##             lines(xrange[v[f:l]], BIC[
##                                       v[f:l], i], lty = i)
##             if(l == n)
##               break
##             f <- l + 1
##           }
##         }
##       }
##     }
##   }
  matplot(xrange, BIC, type="b", xlab=xlab, ylab = ylab, pch=symbols, ...)
            
  symbols <- symbols[1:length(modelNames)]
  names(symbols) <- modelNames
  symbols
}

### dots are not used but inserted to keep R from complaining
### about generic/method consistency...
"plot.EMclustN" <- function(x, modelNames, G, symbols, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  n <- ncol(x)
  dnx <- dimnames(x)
  ##
  x <- matrix(as.vector(x), ncol = n)
  dimnames(x) <- dnx
  if(missing(modelNames))
    modelNames <- dimnames(x)[[2]]
  if(missing(G))
    G <- dimnames(x)[[1]]
  else G <- as.character(G)
  BIC <- x[G, modelNames, drop = FALSE]
  X <- is.na(BIC)
  nrowBIC <- nrow(BIC)
  ncolBIC <- ncol(BIC)
  if(missing(symbols)) {
    if(n > 9) {
      symbols <- LETTERS[1:n]
    }
    else {
      symbols <- as.character(1:n)
    }
    names(symbols) <- dnx[[2]]
    symbols <- symbols[modelNames]
  }
  xrange <- if(!is.null(dn <- dimnames(BIC)[[1]])) as.numeric(dn) else 1:
    nrowBIC
  plot(xrange, BIC[, 1], type = "n", ylim = range(as.vector(BIC[!X])),
       xlim = range(xrange), xlab = "number of clusters", ylab = "BIC"
       )
  for(i in 1:ncolBIC) {
    x <- !X[, i]
    if(any(x)) {
      if(all(x)) {
        points(xrange, BIC[, i], pch = symbols[i])
        if(nrowBIC > 1)
          lines(xrange, BIC[, i], lty = i)
      }
      else {
        points(xrange[x], BIC[x, i], pch = symbols[i])
        v <- (1:nrowBIC)[x]
        n <- length(v)
        d <- c(diff(c(0, v)), 2)
        if(any(d) == 1) {
          f <- 1
          while(TRUE) {
            while(f <= n && d[f] != 1) f <- f + 1
            if(f > n)
              break
            l <- f
            while(d[l + 1] == 1) l <- l + 1
            lines(xrange[v[f:l]], BIC[v[f:l], i], lty = i)
            if(l == n)
              break
            f <- l + 1
          }
        }
      }
    }
  }
  symbols <- symbols[1:length(modelNames)]
  names(symbols) <- modelNames
  symbols
}

"plot.Mclust" <- function(x, data, dimens = c(1, 2), scale = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(data)) {
    warning("data not supplied")
    plot.EMclust(x$BIC)
    return(invisible())
  }
  ### Ask for confirmation before showing the next plot, if users want
  ### ALL plots...
  parSave <- par(no.readonly=TRUE)
  on.exit(par(parSave))
  par(ask=TRUE)
  
  p <- ncol(as.matrix(data))
  if(p > 2) {
    choices <- c("BIC", "Pairs", "Classification (2-D projection)",
                 "Uncertainty (2-D projection)", "All")
    tmenu <- paste("plot:", choices)
    while(TRUE) {
      pick <- menu(tmenu, title = 
                   "\nmake a plot selection (0 to exit):\n")
      if(!pick)
        break
      ALL <- pick == 5
      if(pick == 1 || ALL) {
        plot.EMclust(x$BIC)
      }
      if(pick == 2 || ALL) {
        ##        parSave <- par(no.readonly=TRUE)
        clPairs(data, classification = x$classification, ...)
        ##        par(parSave)
      }
      if(pick == 3 || ALL) {
        do.call("coordProj",
                c(list(data = data, dimens = dimens, scale = scale,
                       identify = FALSE),
                  x[c("classification", "mu", "sigma", "decomp")]))
      }
      if(pick == 4 || ALL) {
        do.call("coordProj",
                c(list(data = data, dimens = dimens, scale = scale,
                       identify = FALSE), 
                  x[c("uncertainty", "mu", "sigma", "decomp")]))
      }
    }
  }
  else if(p == 2) {
    choices <- c("BIC", "Classification", "Uncertainty", "Density",
                 "All")
    tmenu <- paste("plot:", choices)
    while(TRUE) {
      pick <- menu(tmenu, title = 
                   "\nmake a plot selection (0 to exit):\n")
      if(!pick)
        break
      ALL <- pick == 5
      if(pick == 1 || ALL) {
        plot.EMclust(x$BIC)
      }
      if(pick == 2 || ALL) {
        do.call("mclust2Dplot",
                c(list(data = data, scale = scale, identify = FALSE),
                  x[c("classification", "mu", "sigma", "decomp")]))
      }
      if(pick == 3 || ALL) {
        do.call("mclust2Dplot",
                c(list(data = data, scale = scale, identify = FALSE),
                  x[c("uncertainty", "mu", "sigma", "decomp")]))
      }
      if(pick == 4 || ALL) {
        do.call("surfacePlot",
                c(list(data = data, type = "contour", what = "density",
                       transformation = "none", grid = 100,
                       scale = scale, identify = FALSE, verbose = FALSE),
                  x[c("mu", "sigma", "decomp", "pro")]))
      }
    }
  }
  else {
    ##
    ## p == 1
    ##
    choices <- c("BIC", "Classification", "Uncertainty", "Density", "All")
    tmenu <- paste("plot:", choices)
    while(TRUE) {
      pick <- menu(tmenu, title = 
                   "\nmake a plot selection (0 to exit):\n")
      if(!pick)
        break
      ALL <- pick == 5
      if(pick == 1 || ALL) {
        plot.EMclust(x$BIC)
      }
      if(pick == 2 || ALL) {
        do.call("mclust1Dplot",
                c(list(data = data, type = "classification", identify = FALSE),
                  x))
      }
      if(pick == 3 || ALL) {
        do.call("mclust1Dplot",
                c(list(data = data, type = "uncertainty", identify = FALSE),
                  x))
      }
      if(pick == 4 || ALL) {
        do.call("mclust1Dplot",
                c(list(data = data, type = "density", identify = FALSE), x))
      }
    }
  }
  invisible()
}

"plot.mclustDA" <- 
  function(x, trainingData, labels, testData, dimens = c(1, 2), scale = FALSE, 
           identify = FALSE, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(A <- missing(trainingData)) {
    warning("training data not supplied")
  }
  if(B <- missing(labels)) {
    warning("training labels not supplied")
  }
  if(C <- missing(testData)) {
    warning("test data not supplied")
  }
  if(A || B || C)
    return(invisible())
  p <- ncol(as.matrix(trainingData))
  if(p > 2) {
    if(is.null(dim(testData)))
      testData <- matrix(testData, nrow = 1, ncol = p)
    choices <- c("Training and Test Data (Pairs Plot)", 
                 "Training and Test Data (2D projection)", 
                 "Training Data - known classification", 
                 "Test Data - mclustDA classification", 
                 "Training Data - misclassified observations", "All")
    tmenu <- paste("plot:", choices)
    Data <- rbind(as.matrix(testData), as.matrix(trainingData))
    xlim <- range((Data[, dimens])[, 1])
    ylim <- range((Data[, dimens])[, 2])
    cl <- c(rep(1, nrow(testData)), rep(2, nrow(trainingData)))
    while(TRUE) {
      pick <- menu(tmenu, title = 
                   "\nplot.mclustDA : make a plot selection (0 to exit):\n"
                   )
      tmenu <- paste("plot:", choices)
      if(!pick)
        break
      ALL <- pick == 6
      if(pick == 1 || ALL) {
        ##        parSave <- par(no.readonly=TRUE)
        clPairs(Data, classification = cl, symbols = c(1, 3), ...)
        ##        par(parSave)
      }
      if(pick == 2 || ALL) {
        coordProj(data = Data, dimens = dimens, 
                  classification = cl, scale = scale,
                  identify = FALSE, symbols = c(1, 3), ...)
        if(identify)
          title("Training and Test Data", cex = 0.75)
      }
      if(pick == 3 || ALL) {
        coordProj(data = trainingData, dimens = dimens,
                  classification = labels, scale = scale,
                  identify = FALSE, xlim = xlim, ylim = ylim, ...)
        if(identify)
          title("Training Data: known Classification", cex = 0.75)
      }
      if(pick == 4 || ALL) {
        coordProj(data = testData, dimens = dimens,
                  classification = x$testClass, scale = 
                  scale, identify = FALSE, xlim = xlim, ylim = ylim, ...)
        if(identify)
          title("Test Data: mclustDA Classification", cex = 0.75)
      }
      if(pick == 5 || ALL) {
        coordProj(data = trainingData, dimens = dimens,
                  classification = x$trainingClass, truth = labels,
                  type = "errors", scale = scale, identify = FALSE,
                  xlim = xlim, ylim = ylim, ...)
        if(identify)
          title("Training Error", cex = 0.75)
      }
    }
  }
  else if(p == 2) {
    if(is.null(dim(testData)))
      testData <- matrix(testData, nrow = 1, ncol = p)
    choices <- c("Training and Test Data", 
                 "Training Data - known classification", 
                 "Test Data - mclustDA classification", 
                 "Training Data - misclassified observations", "All")
    tmenu <- paste("plot:", choices)
    Data <- rbind(as.matrix(testData), as.matrix(trainingData))
    xlim <- range((Data[, dimens])[, 1])
    ylim <- range((Data[, dimens])[, 2])
    cl <- c(rep(1, nrow(testData)), rep(2, nrow(trainingData)))
    while(TRUE) {
      pick <- menu(tmenu, title = 
                   "\nmake a plot selection (0 to exit):\n")
      if(!pick)
        break
      ALL <- pick == 5
      if(pick == 1 || ALL) {
        mclust2Dplot(data = Data, classification = cl,
                     scale = scale, identify = FALSE, symbols = 
                     c(1, 3), ...)
        if(identify)
          title("Training and Test Data", cex = 
                0.75)
      }
      if(pick == 2 || ALL) {
        mclust2Dplot(data = trainingData, dimens = 
                     dimens, classification = labels, scale
                     = scale, identify = FALSE, xlim = xlim,
                     ylim = ylim, ...)
        if(identify)
          title("Training Data: known Classification",
                cex = 0.75)
      }
      if(pick == 3 || ALL) {
        mclust2Dplot(data = testData, dimens = dimens,
                     classification = x$testClass, scale = 
                     scale, identify = FALSE, xlim = xlim, ylim
                     = ylim, ...)
        if(identify)
          title("Test Data: mclustDA Classification",
                cex = 0.75)
      }
      if(pick == 4 || ALL) {
        mclust2Dplot(data = trainingData, dimens = 
                     dimens, classification = x$
                     trainingClass, truth = labels, type = 
                     "errors", scale = scale, identify = FALSE,
                     xlim = xlim, ylim = ylim, ...)
        if(identify)
          title("Training Error", cex = 0.75)
      }
    }
  }
  else {
    ##
    ## p == 1
    ##
    Data <- c(testData, trainingData)
    xlim <- range(Data)
    cl <- c(rep(1, length(testData)), rep(2, length(trainingData)))
    choices <- c("Training and Test Data", 
                 "Training Data - known classification", 
                 "Test Data - mclustDA classification", 
                 "Training Data - misclassified observations", "All")
    tmenu <- paste("plot:", choices)
    while(TRUE) {
      pick <- menu(tmenu, if(identify) title = 
                   "\nmake a plot selection (0 to exit):\n"
                   )
      if(!pick)
        break
      ALL <- pick == 5
      if(pick == 1 || ALL) {
        mclust1Dplot(data = Data, classification = cl,
                     identify = FALSE, xlim = xlim, ...)
        if(identify)
          title("Training and Test Data", cex = 
                0.75)
      }
      if(pick == 2 || ALL) {
        mclust1Dplot(data = trainingData, 
                     classification = labels, scale = scale,
                     identify = FALSE, xlim = xlim, ...)
        if(identify)
          title("Training Data: known Classification",
                cex = 0.75)
      }
      if(pick == 3 || ALL) {
        mclust1Dplot(data = testData, classification = 
                     x$testClass, scale = scale, identify = 
                     FALSE, xlim = xlim, ...)
        if(identify)
          title("Test Data mclustDA Classification",
                cex = 0.75)
      }
      if(pick == 4 || ALL) {
        mclust1Dplot(data = trainingData, dimens = 
                     dimens, classification = x$
                     trainingClass, truth = labels, type = 
                     "errors", identify = FALSE, xlim = xlim,
                     ...)
        if(identify)
          title("Training Error", cex = 0.75)
      }
    }
  }
  invisible()
}

"print.EMclust" <- function(x, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##

  subset <- !is.null(attr(x, "subset"))
  class(x) <- attr(x, "args") <- NULL
  attr(x, "hcPairs") <- attr(x, "attrHC") <- NULL
  attr(x, "eps") <- attr(x, "tol") <- attr(x, "itmax") <- NULL
  attr(x, "equalPro") <- attr(x, "subset") <- NULL
  attr(x, "warnSingular") <- NULL
  cat("\n BIC:\n")
  NextMethod("print")
  cat("\n")
  invisible()
}

"print.EMclustN" <- function(x, ...)
{
  class(x) <- attr(x, "args") <- NULL
  attr(x, "hcPairs") <- attr(x, "attrHC") <- NULL
  attr(x, "eps") <- attr(x, "tol") <- attr(x, "itmax") <- NULL
  attr(x, "equalPro") <- attr(x, "subset") <- NULL
  attr(x, "warnSingular") <- NULL
  attr(x, "noise") <- attr(x, "Vinv") <- NULL
  cat("\n BIC:\n")
  NextMethod("print")
  cat("\n")
  invisible()
}

"print.Mclust" <- function(x, ndigits = options()$digits, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  M <- switch(EXPR = x$model,
              X = "univariate normal",
              E = "equal variance",
              V = "unequal variance",
              XII = "spherical multivariate normal",
              EII = "spherical, equal volume",
              VII = "spherical, varying volume",
              XXI = "diagonal multivariate normal",
              EEI = "diagonal, equal volume and shape",
              VEI = "diagonal, equal shape",
              EVI = "diagonal, equal volume",
              VVI = "diagonal, varying volume and shape",
              EEE = "elliposidal, equal variance",
              VEE = "elliposidal, equal shape and orientation",
              EVE = "elliposidal, equal volume and orientation",
              VVE = "ellipsoidal, equal orientation",
              EEV = "ellipsoidal, equal volume and shape",
              VEV = "ellipsoidal, equal shape",
              EVV = "ellipsoidal, equal volume",
              XXX = "ellipsoidal multivariate normal",
              VVV = "ellipsoidal, unconstrained",
              stop("invalid model id for EM"))
  G <- length(unique(x$classification))
  cat("\n best model:", M, "with", G, "groups\n")
  aveUncer <- round(mean(x$uncertainty), 3)
  medUncer <- round(median(x$uncertainty), 3)
  cat("\n averge/median classification uncertainty:", aveUncer, "/",
      medUncer, "\n\n")
  invisible()
}

"print.mclustDA" <- function(x, ndigits = options()$digits, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  cat("\n")
  print(x$summary)
  ##temporarily commented out... Ron, May 2003
  ##  cat("\n training error rate:", round(x$errorRate, 3), "\n")

  ##commented out by Chris
  ## uncer <- 1 - apply(x$postProb, 1, max)
  ## aveUncer <- round(mean(uncer), 3)
  ## medUncer <- round(median(uncer), 3)
  ##
  invisible()
}


"print.summary.EMclust" <- function(x, ndigits = options()$digits, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  bic <- x$bic
  l <- length(bic) > 1
  cat("\nclassification table:\n")
  print(table(x$classification), ...)
  cat("\nuncertainty (quartiles):\n")
  print(quantile(x$uncertainty), digits = ndigits, ...)
  if(l)
    cat("\nbest BIC values:\n")
  else cat("\nbest BIC value:\n")
  print(round(bic, ndigits))
  M <- switch(EXPR = x$model,
              X = "univariate normal",
              E = "equal variance",
              V = "unequal variance",
              XII = "spherical multivariate normal",
              EII = "spherical, equal volume",
              VII = "spherical, varying volume",
              XXI = "diagonal multivariate normal",
              EEI = "diagonal, equal volume and shape",
              VEI = "diagonal, equal shape",
              EVI = "diagonal, equal volume",
              VVI = "diagonal, varying volume and shape",
              EEE = "elliposidal, equal variance",
              VEE = "elliposidal, equal shape and orientation",
              EVE = "elliposidal, equal volume and orientation",
              VVE = "ellipsoidal, equal orientation",
              EEV = "ellipsoidal, equal volume and shape",
              VEV = "ellipsoidal, equal shape",
              EVV = "ellipsoidal, equal volume",
              XXX = "ellipsoidal multivariate normal",
              VVV = "ellipsoidal, unconstrained",
              stop("invalid model id for EM"))
  cat("\nbest model:", M, "\n\n")
  ##
  ##	print(x$options)
  invisible()
}

"print.summary.EMclustN" <- 
  function(x, ndigits = options()$digits, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  bic <- x$bic
  l <- length(bic) > 1
  cat("\nclassification table:\n")
  G <- x$G
  y <- rep(0, G)
  names(y) <- as.numeric(1:G)
  tab <- table(x$classification)
  y[names(tab)] <- tab
  print(y, ...)
  cat("\nuncertainty (quartiles):\n")
  print(quantile(x$uncertainty), digits = ndigits, ...)
  if(l)
    cat("\nbest BIC values:\n")
  else cat("\nbest BIC value:\n")
  print(round(bic, ndigits))
  M <- switch(EXPR = x$model,
              X = "univariate normal",
              E = "equal variance",
              V = "unequal variance",
              XII = "spherical multivariate normal",
              EII = "spherical, equal volume",
              VII = "spherical, varying volume",
              XXI = "diagonal multivariate normal",
              EEI = "diagonal, equal volume and shape",
              VEI = "diagonal, equal shape",
              EVI = "diagonal, equal volume",
              VVI = "diagonal, varying volume and shape",
              EEE = "elliposidal, equal variance",
              VEE = "elliposidal, equal shape and orientation",
              EVE = "elliposidal, equal volume and orientation",
              VVE = "ellipsoidal, equal orientation",
              EEV = "ellipsoidal, equal volume and shape",
              VEV = "ellipsoidal, equal shape",
              EVV = "ellipsoidal, equal volume",
              XXX = "ellipsoidal multivariate normal",
              VVV = "ellipsoidal, unconstrained",
              stop("invalid model id for EM"))
  cat("\nbest model:", M, "\n\n")
  ##
  ##	print(x$options)
  G <- x$G + 1
  if(length(tab) != G) {
    I <- (1:G)[!match(1:G, as.numeric(names(tab)), nomatch = 0)]
    cat("***", paste("no assignment to", paste(I, collapse = ",")),
        " ***\n\n")
  }
  invisible()
}

"randProj" <- function(data, seeds = 0, ...,
                       type = c("classification", "uncertainty", "errors"),
                       ask = TRUE, quantiles = c(0.75, 0.95), symbols,
                       scale = FALSE, identify = FALSE, CEX = 1,
                       PCH = ".", xlim, ylim) 
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(scale)
    par(pty = "s")
  aux <- list(...)
  z <- aux$z
  classification <- aux$classification
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  uncertainty <- aux$uncertainty
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  truth <- aux$truth
  mu <- aux$mu
  sigma <- aux$sigma
  decomp <- aux$decomp
  params <- !is.null(mu) && (!is.null(sigma) || !is.null(decomp))
  if(!is.null(mu)) {
    if(is.null(sigma)) {
      if(is.null(decomp)) {
        params <- FALSE
        warning("covariance not supplied")
      }
      else {
        sigma <- decomp2sigma(decomp)
      }
    }
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      params <- FALSE
      warning("covariance improperly specified")
    }
    if(G != dimpar[3]) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
  }
  p <- ncol(data)
  if(params)
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(
                 classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    if(missing(symbols)) {
      if(L <= length(.Mclust$symbols)) {
        symbols <- .Mclust$symbols
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    if(length(symbols) < L) {
      warning("more symbols needed to show classification")
      classification <- NULL
    }
  }
  if(l <- length(type)) {
    choices <- c("classification", "uncertainty", "density", 
                 "errors")
    m <- rep(0, l)
    for(i in 1:l) {
      m[i] <- charmatch(type[i], choices, nomatch = 0)
    }
    choices <- choices[unique(m)]
    if(is.null(classification))
      choices <- choices[choices != "classification"]
    if(is.null(uncertainty))
      choices <- choices[choices != "uncertainty"]
    if(is.null(truth))
      choices <- choices[choices != "errors"]
  }
  else choices <- NULL
  if(length(choices) > 1 && ask)
    choices <- c(choices, "all")
  else if(length(choices) == 1)
    ask <- FALSE
  if(any(choices == "errors")) {
    ERRORS <- classErrors(classification, truth)
  }
  if(!ask)
    pick <- 1:length(choices)
  ALL <- FALSE
  xlimMISS <- missing(xlim)
  ylimMISS <- missing(ylim)
  for(seed in seeds) {
    set.seed(seed)
    Q <- orth2(p)
    Data <- as.matrix(data) %*% Q
    if(xlimMISS)
      xlim <- range(Data[, 1])
    if(ylimMISS)
      ylim <- range(Data[, 2])
    if(scale) {
      d <- diff(xlim) - diff(ylim)
      if(d > 0) {
        ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
      }
      else {
        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
      }
    }
    if(!length(choices)) {
      plot(Data[, 1], Data[, 2], type = "n", xlab = "", ylab
           = "", xlim = xlim, ylim = ylim, ...)
      if(params) {
        Mu <- crossprod(Q, mu)
        Sigma <- array(apply(cho, 3, function(R, Q)
                             crossprod(R %*% Q), Q = Q), c(2, 2, G))
        for(k in 1:G) {
          mvn2plot(mu = Mu[, k], sigma = Sigma[
                                   ,  , k], k = 15)
        }
      }
      points(Data[, 1], Data[, 2], pch = PCH, cex = CEX)
      if(identify)
        title(paste("Random Projection: seed = ", seed,
                    collapse = ""), cex = 0.5)
      next
    }
    while(TRUE) {
      if(ask) {
        pick <-
          menu(choices,
               title = paste("\nrandProj: make a plot selection (0 to exit) seed =",
                 seed, "\n", collapse = ""))
        if(!pick)
          return(invisible())
        ALL <- any(choices[pick] == "all")
      }
      if(any(choices[pick] == "classification") || (any(
                      choices == "classification") && ALL)) {
        plot(Data[, 1], Data[, 2], type = "n", xlab = 
             "", ylab = "", xlim = xlim, ylim = ylim,
             ...)
        if(params) {
          Mu <- crossprod(Q, mu)
          Sigma <- array(apply(cho, 3, function(R,
						Q)
                               crossprod(R %*% Q), Q = Q), c(2, 2,
                                                     G))
          for(k in 1:G) {
            mvn2plot(mu = Mu[, k], sigma = 
                     Sigma[,  , k], k = 15)
          }
        }
        for(k in 1:L) {
          I <- classification == U[k]
          points(Data[I, 1], Data[I, 2], pch = 
                 symbols[k], cex = CEX)
        }
        if(identify)
          title(paste(
                      "Random Projection showing Classification: seed = ",
                      seed, collapse = ""), cex = 0.5
                )
      }
      if(any(choices[pick] == "uncertainty") || (any(choices ==
                      "uncertainty") && ALL)) {
        plot(Data[, 1], Data[, 2], type = "n", xlab = 
             "", ylab = "", xlim = xlim, ylim = ylim,
             ...)
        if(params) {
          Mu <- crossprod(Q, mu)
          Sigma <- array(apply(cho, 3, function(R,
						Q)
                               crossprod(R %*% Q), Q = Q), c(2, 2,
                                                     G))
          for(k in 1:G) {
            mvn2plot(mu = Mu[, k], sigma = 
                     Sigma[,  , k], k = 15)
          }
        }
        breaks <- quantile(uncertainty, probs = sort(
                                          quantiles))
        I <- uncertainty < breaks[1]
        points(Data[I, 1], Data[I, 2], pch = 16, cex = 
               0.5 * CEX)
        I <- uncertainty < breaks[2] & !I
        points(Data[I, 1], Data[I, 2], pch = 1, cex = 1 *
               CEX)
        I <- uncertainty >= breaks[2]
        points(Data[I, 1], Data[I, 2], pch = 16, cex = 
               1.5 * CEX)
        if(identify)
          title(paste(
                      "Random Projection showing Classification Uncertainty: seed = ",
                      seed, collapse = ""), cex = 0.5
                )
      }
      if(any(choices[pick] == "errors") || (any(choices ==
                      "errors") && ALL)) {
        plot(Data[, 1], Data[, 2], type = "n", xlab = 
             "", ylab = "", xlim = xlim, ylim = ylim,
             ...)
        if(params) {
          Mu <- crossprod(Q, mu)
          Sigma <- array(apply(cho, 3, function(R,
						Q)
                               crossprod(R %*% Q), Q = Q), c(2, 2,
                                                     G))
          for(k in 1:G) {
            mvn2plot(mu = Mu[, k], sigma = 
                     Sigma[,  , k], k = 15)
          }
        }
        CLASSES <- unique(as.character(truth))
        symOpen <- c(2, 0, 1, 5)
        symFill <- c(17, 15, 16, 18)
        good <- !ERRORS
        if(L > 4) {
          points(Data[good, 1], Data[good, 2],
                 pch = 1, cex = CEX)
          points(Data[!good, 1], Data[!good,
                                      2], pch = 16, cex = CEX)
        }
        else {
          for(k in 1:L) {
            K <- truth == CLASSES[k]
            points(Data[K, 1], Data[K,
                                    2], pch = symOpen[
                                          k], cex = CEX)
            if(any(I <- (K & ERRORS))) {
              points(Data[I, 1],
                     Data[I, 2],
                     pch = symFill[
                       k], cex = CEX)
            }
          }
        }
        if(identify)
          title(paste(
                      "Random Projection showing Classification Errors: seed = ",
                      seed, collapse = ""), cex = 0.5
                )
      }
      if(!ask)
        break
    }
  }
  invisible()
}

"shapeO" <- function(shape, O, transpose = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimO <- dim(O)
  if(dimO[1] != dimO[2])
    stop("leading dimensions of O are unequal")
  if((ldO <- length(dimO)) != 3) {
    if(ldO == 2) {
      dimO <- c(dimO, 1)
      O <- array(O, dimO)
    }
    else stop("O must be a matrix or an array")
  }
  l <- length(shape)
  if(l != dimO[1])
    stop("dimension of O and length s are unequal")
  storage.mode(O) <- "double"
  .Fortran("shapeo",
           as.integer(if(transpose) 1 else 0),
           as.double(shape),
           O,
           as.integer(l),
           as.integer(dimO[3]),
           double(l * l),
           integer(1),
                   PACKAGE="mclust")[[3]]
}

"sigma2decomp" <- function(sigma, G, tol, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  dimSigma <- dim(sigma)
  if(is.null(dimSigma))
    stop("sigma improperly specified")
  d <- dimSigma[1]
  if(dimSigma[2] != d)
    stop("sigma improperly specified")
  l <- length(dimSigma)
  if(l < 2 || l > 3)
    stop("sigma improperly specified")
  if(missing(G)) {
    if(l == 2) {
      G <- 1
      sigma <- array(sigma, c(dimSigma, 1))
    }
    else {
      G <- dimSigma[3]
    }
  }
  else {
    if(l == 3 && G != dimSigma[3])
      stop("sigma and G are incompatible")
    if(l == 2 && G != 1)
      stop("sigma and G are incompatible")
  }
  decomp <- list(d = d, G = G, scale = rep(0, G),
                 shape = matrix(0, d, G), orientation = array(0, c(d, d, G)))
  for(k in 1:G) {
    ev <- eigen(sigma[,  , k], symmetric = TRUE)
    temp <- logb(ev$values)
    logScale <- sum(temp)/d
    decomp$scale[k] <- exp(logScale)
    decomp$shape[, k] <- exp(temp - logScale)
    decomp$orientation[,  , k] <- t(ev$vectors)
  }
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  scaleName <- "V"
  shapeName <- "V"
  orientName <- "V"
  "uniq" <- 
    function(x, tol = sqrt(.Machine$double.eps))
      {
        abs(max(x) - min(x)) < tol
      }

  if(uniq(decomp$scale)) {
    decomp$scale <- decomp$scale[1]
    scaleName <- "E"
  }
  if(all(apply(decomp$shape, 1, uniq, tol = tol))) {
    decomp$shape <- decomp$shape[, 1]
    if(all(uniq(decomp$shape, tol = tol))) {
      shapeName <- "I"
      decomp$shape <- rep(1, d)
    }
    else {
      shapeName <- "E"
    }
  }
  if(all(apply(matrix(decomp$orientation, nrow = d * d, ncol = G), 1,
               uniq, tol = tol))) {
    decomp$orientation = decomp$orientation[,  , 1]
    if(all(apply(cbind(decomp$orientation, diag(d)), 1, uniq, tol
                 = tol))) {
      orientName <- "I"
      decomp$orientation <- NULL
    }
    else {
      orientName <- "E"
    }
  }
  structure(decomp, modelName = paste(c(scaleName, shapeName, orientName),
                      collapse = ""))
}

"sim" <- function(modelName, mu, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## ... variance parameters, n
  funcName <- paste("sim", modelName, sep = "")
  do.call(funcName, list(mu = mu, ..., seed = seed))
}

"simE" <- function(mu, sigmasq, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  G <- length(mu)
  n <- list(...)$n
  if(all(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(rep(NA, n), modelName = "E"))
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size = n)
  ctabel <- table(clabels)
  x <- rep(0, n)
  sd <- sqrt(sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd)
  }
  structure(x, classification = clabels, modelName = "E")
}

"simEEE" <- function(mu, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  cholSigma <- list(...)$cholSigma
  if(is.null(cholSigma)) {
    if(!is.null(decomp <- list(...)$decomp)) {
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      cholSigma <- qr.R(qr(O * sqrt(scale * shape)))
    }
    else if(!is.null(Sigma <- list(...)$Sigma)) {
      cholSigma <- chol(Sigma)
    }
    else if(!missing(sigma)) {
      cholSigma <- chol(Sigma)
    }
    else stop("invalid specification for sigma")
  }
  if(all(is.na(c(mu, cholSigma)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EEE"))
  }
  if(any(is.na(c(mu, cholSigma)))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
        FUN = "+")
  }
  structure(x, classification = clabels, modelName = "EEE")
}

"simEEI" <- function(mu, decomp, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(decomp))
    stop("decomp must be specified")
  if(all(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EEI"))
  }
  if(any(is.na(c(mu, unlist(decomp))))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size=n)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  shape <- decomp$shape
  if(length(shape != d))
    stop("shape incompatible with mu")
  cholSigma <- diag(sqrt(decomp$scale * shape))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
        FUN = "+")
  }
  structure(x, classification = clabels, modelName = "EEI")
}

"simEEV" <- function(mu, decomp, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(decomp))
    stop("decomp must be specified")
  if(all(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EEV"))
  }
  if(any(is.na(c(mu, unlist(decomp))))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  sss <- sqrt(decomp$scale * decomp$shape)
  for(k in 1:G) {
    m <- ctabel[k]
    cholSigma <- decomp$orientation[,  , k] * sss
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
        FUN = "+")
  }
  structure(x, classification = clabels, modelName = "EEV")
}

"simEII" <- function(mu, sigmasq, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(sigmasq)) {
    sigmasq <- list(...)$decomp$scale
  }
  if(all(is.na(c(mu, sigmasq)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EII"))
  }
  if(any(is.na(c(mu, sigmasq)))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size=n)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  cholSigma <- diag(rep(sqrt(sigmasq), d))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
        FUN = "+")
  }
  structure(x, classification = clabels, modelName = "EII")
}

"simEVI" <- function(mu, decomp, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(decomp))
    stop("decomp must be specified")
  if(all(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EVI"))
  }
  if(any(is.na(c(mu, unlist(decomp))))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size = n)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  shape <- decomp$shape
  if(nrow(shape) != d)
    stop("shape incompatible with mu")
  sss <- sqrt(decomp$scale * shape)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% diag(sss[, k]), MARGIN = 2, STAT = mu[
                                                    , k], FUN = "+")
  }
  structure(x, classification = clabels, modelName = "EVI")
}

"simV" <- 
  function(mu, sigmasq, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  G <- length(mu)
  n <- list(...)$n
  if(all(is.na(c(mu, sigmasq, pro)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(rep(NA, n), modelName = "V"))
  }
  if(any(is.na(c(mu, sigmasq, pro)))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size = n)
  ctabel <- table(clabels)
  x <- rep(0, n)
  sd <- sqrt(sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- rnorm(ctabel[k], sd = sd[k])
  }
  structure(x, classification = clabels, modelName = "V")
}

"simVEI" <- 
  function(mu, decomp, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(decomp))
    stop("decomp must be specified")
  if(all(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VEI"))
  }
  if(any(is.na(c(mu, unlist(decomp))))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size = n)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtscale <- sqrt(decomp$scale)
  rtshape <- sqrt(decomp$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mu")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% diag(rtscale[k] * rtshape), MARGIN = 2,
        STAT = mu[, k], FUN = "+")
  }
  structure(x, classification = clabels, modelName = "VEI")
}

"simVEV" <- function(mu, decomp, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(decomp))
    stop("decomp must be specified")
  if(all(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VEV"))
  }
  if(any(is.na(c(mu, unlist(decomp))))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtscale <- sqrt(decomp$scale)
  rtshape <- sqrt(decomp$shape)
  for(k in 1:G) {
    m <- ctabel[k]
    cholSigma <- decomp$orientation[,  , k] * (rtscale[k] * rtshape
                                               )
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
        FUN = "+")
  }
  structure(x, classification = clabels, modelName = "VEV")
}

"simVII" <- function(mu, sigmasq, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(sigmasq)) {
    sigmasq <- list(...)$decomp$scale
  }
  if(all(is.na(c(mu, sigmasq)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VII"))
  }
  if(any(is.na(c(mu, sigmasq)))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size = n)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% diag(rep(sqrt(sigmasq[k]), d)), MARGIN = 
        2, STAT = mu[, k], FUN = "+")
  }
  structure(x, classification = clabels, modelName = "VII")
}

"simVVI" <- 
  function(mu, decomp, pro, ..., seed = 0)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(decomp))
    stop("decomp must be specified")
  if(all(is.na(c(mu, unlist(decomp))))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VVI"))
  }
  if(any(is.na(c(mu, unlist(decomp))))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, prob = pro, replace=TRUE, size = n)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtscale <- sqrt(decomp$scale)
  rtshape <- sqrt(decomp$shape)
  if(nrow(rtshape) != d)
    stop("shape incompatible with mu")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% diag(rtscale[k] * rtshape[, k]), MARGIN
        = 2, STAT = mu[, k], FUN = "+")
  }
  structure(x, classification = clabels, modelName = "VVI")
}

"simVVV" <- function(mu, pro, ..., seed = 0)
{
  ##
  ## This function is part of the MCLUST software described at
  ##       http://www.stat.washington.edu/mclust
  ## Copyright information and conditions for use of MCLUST are given at
  ##        http://www.stat.washington.edu/mclust/license.txt
  ## Distribution of MCLUST is prohibited except by agreement with the 
  ## University of Washington.
  ##
  mu <- as.matrix(mu)
  d <- nrow(mu)
  G <- ncol(mu)
  n <- list(...)$n
  if(is.null(cholsigma <- list(...)$cholsigma)) {
    if(!is.null(sigma <- list(...)$sigma)) {
      cholsigma <- apply(sigma, 3, chol)
      for(i in 1:ncol(cholsigma))
        sigma[,  , i] <- cholsigma[, i]
      cholsigma <- sigma
    }
    else if(!is.null(decomp <- list(...)$decomp)) {
      scale <- decomp$scale
      shape <- decomp$shape
      O <- decomp$orientation
      cholsigma <- array(0, c(p, p, G))
      shape <- sqrt(sweep(shape, MARGIN = 2, STATS = scale,
                          FUN = "*"))
      for(k in 1:G)
        cholsigma[,  , k] <- qr.R(qr(O[,  , k] * shape)
                                  )
    }
    else stop("sigma improperly specified")
  }
  if(all(is.na(c(mu, cholsigma)))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VVV"))
  }
  if(any(is.na(c(mu, cholsigma)))) {
    stop("parameters contain missing values")
  }
  if(missing(pro))
    pro <- rep(1/G, G)
  set.seed(seed)
  clabels <- sample(1:G, size = n, replace = T, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
        ncol = d) %*% cholsigma[,  , k], MARGIN = 2,
        STAT = mu[, k], FUN = "+")
  }
  structure(x, classification = clabels, modelName = "VVV")
}


"spinProj" <- function(data, ..., angles = c(0, pi/3, (2 * pi)/3, pi),
                       seed = 0, reflection = FALSE,
                       type = c("classification", "uncertainty", "errors"),
                       ask = TRUE,
                       quantiles = c(0.75, 0.95), symbols,
                       scale = FALSE, identify = FALSE, CEX = 1,
                       PCH = ".", xlim, ylim)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(scale)
    par(pty = "s")
  data <- as.matrix(data)
  aux <- list(...)
  z <- aux$z
  classification <- aux$classification
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  uncertainty <- aux$uncertainty
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  truth <- aux$truth
  mu <- aux$mu
  sigma <- aux$sigma
  decomp <- aux$decomp
  params <- !is.null(mu) && (!is.null(sigma) || !is.null(decomp))
  if(!is.null(mu)) {
    if(is.null(sigma)) {
      if(is.null(decomp)) {
        params <- FALSE
        warning("covariance not supplied")
      }
      else {
        sigma <- decomp2sigma(decomp)
      }
    }
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      params <- FALSE
      warning("covariance improperly specified")
    }
    if(G != dimpar[3]) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
  }
  p <- ncol(data)
  if(params)
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(
                 classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    if(missing(symbols)) {
      if(L <= length(.Mclust$symbols)) {
        symbols <- .Mclust$symbols
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    if(length(symbols) < L) {
      warning("more symbols needed to show classification")
      classification <- NULL
    }
  }
  if(l <- length(type)) {
    choices <- c("classification", "uncertainty", "density", 
                 "errors")
    m <- rep(0, l)
    for(i in 1:l) {
      m[i] <- charmatch(type[i], choices, nomatch = 0)
    }
    choices <- choices[unique(m)]
    if(is.null(classification))
      choices <- choices[choices != "classification"]
    if(is.null(uncertainty))
      choices <- choices[choices != "uncertainty"]
    if(is.null(truth))
      choices <- choices[choices != "errors"]
  }
  if(length(choices) > 1 && ask)
    choices <- c(choices, "all")
  else if(length(choices) == 1)
    ask <- FALSE
  if(any(choices == "errors")) {
    ERRORS <- classErrors(classification, truth)
  }
  if(!ask)
    pick <- 1:length(choices)
  all <- FALSE
  set.seed(seed)
  O <- orth2(p)
  xlimMISS <- missing(xlim)
  ylimMISS <- missing(ylim)
  for(angle in angles) {
    cosTheta <- cos(angle)
    sinTheta <- sin(angle)
    if(reflection) {
      Q <- O %*% matrix(c(cosTheta, sinTheta, sinTheta,  -
                          cosTheta), 2, 2)
    }
    else {
      Q <- O %*% matrix(c(cosTheta,  - sinTheta, sinTheta,
                          cosTheta), 2, 2)
    }
    Data <- data %*% Q
    if(xlimMISS)
      xlim <- range(Data[, 1])
    if(ylimMISS)
      ylim <- range(Data[, 2])
    if(scale) {
      d <- diff(xlim) - diff(ylim)
      if(d > 0) {
        ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
      }
      else {
        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
      }
    }
    if(!length(choices)) {
      plot(Data[, 1], Data[, 2], type = "n", xlab = "", ylab
           = "", xlim = xlim, ylim = ylim, ...)
      if(params) {
        Mu <- crossprod(Q, mu)
        Sigma <- array(apply(cho, 3, function(R, Q)
                             crossprod(R %*% Q), Q = Q), c(2, 2, G))
        for(k in 1:G) {
          mvn2plot(mu = Mu[, k], sigma = Sigma[
                                   ,  , k], k = 15)
        }
      }
      points(Data[, 1], Data[, 2], pch = PCH, cex = CEX)
      if(identify)
        title(paste("Spin Projection: seed/angle = ",
                    seed, "/", round(angle, 3), collapse = 
                    ""), cex = 0.5)
      next
    }
    while(TRUE) {
      if(ask) {
        pick <- menu(choices, title = paste(
                                "\nspinProj: make a plot selection (0 to exit), seed/angle = ",
                                seed, "/", round(angle, 3), "\n", 
                                collapse = ""))
        if(!pick)
          break
        ALL <- any(choices[pick] == "all")
      }
      if(!pick)
        break
      if(any(choices[pick] == "classification") || (any(
                      choices == "classification") && ALL)) {
        plot(Data[, 1], Data[, 2], type = "n", xlab = 
             "", ylab = "", xlim = xlim, ylim = ylim,
             ...)
        if(params) {
          Mu <- crossprod(Q, mu)
          Sigma <- array(apply(cho, 3, function(R,
						Q)
                               crossprod(R %*% Q), Q = Q), c(2, 2,
                                                     G))
          for(k in 1:G) {
            mvn2plot(mu = Mu[, k], sigma = 
                     Sigma[,  , k], k = 15)
          }
        }
        for(k in 1:L) {
          I <- classification == U[k]
          points(Data[I, 1], Data[I, 2], pch = 
                 symbols[k], cex = CEX)
        }
        if(identify)
          title(paste(
                      "Spin Projection showing Classification: seed/angle = ",
                      seed, "/", round(angle, 3),
                      collapse = ""), cex = 0.5)
      }
      if(any(choices[pick] == "uncertainty") || (any(choices ==
                      "uncertainty") && ALL)) {
        plot(Data[, 1], Data[, 2], type = "n", xlab = 
             "", ylab = "", xlim = xlim, ylim = ylim,
             ...)
        if(params) {
          Mu <- crossprod(Q, mu)
          Sigma <- array(apply(cho, 3, function(R,
						Q)
                               crossprod(R %*% Q), Q = Q), c(2, 2,
                                                     G))
          for(k in 1:G) {
            mvn2plot(mu = Mu[, k], sigma = 
                     Sigma[,  , k], k = 15)
          }
        }
        breaks <- quantile(uncertainty, probs = sort(
                                          quantiles))
        I <- uncertainty < breaks[1]
        points(Data[I, 1], Data[I, 2], pch = 16, cex = 
               0.5 * CEX)
        I <- uncertainty < breaks[2] & !I
        points(Data[I, 1], Data[I, 2], pch = 1, cex = 1 *
               CEX)
        I <- uncertainty >= breaks[2]
        points(Data[I, 1], Data[I, 2], pch = 16, cex = 
               1.5 * CEX)
        if(identify)
          title(paste(
                      "Spin Projection showing Classification Uncertainty: seed/angle = ",
                      seed, "/", round(angle, 3),
                      collapse = ""), cex = 0.5)
      }
      if(any(choices[pick] == "errors") || (any(choices ==
                      "errors") && ALL)) {
        plot(Data[, 1], Data[, 2], type = "n", xlab = 
             "", ylab = "", xlim = xlim, ylim = ylim,
             ...)
        if(params) {
          Mu <- crossprod(Q, mu)
          Sigma <- array(apply(cho, 3, function(R,
						Q)
                               crossprod(R %*% Q), Q = Q), c(2, 2,
                                                     G))
          for(k in 1:G) {
            mvn2plot(mu = Mu[, k], sigma = 
                     Sigma[,  , k], k = 15)
          }
        }
        CLASSES <- unique(as.character(truth))
        symOpen <- c(2, 0, 1, 5)
        symFill <- c(17, 15, 16, 18)
        good <- !ERRORS
        if(L > 4) {
          points(data[good, 1], data[good, 2],
                 pch = 1, cex = CEX)
          points(data[!good, 1], data[!good,
                                      2], pch = 16, cex = CEX)
        }
        else {
          for(k in 1:L) {
            K <- truth == CLASSES[k]
            points(data[K, 1], data[K,
                                    2], pch = symOpen[
                                          k], cex = CEX)
            if(any(I <- (K & ERRORS))) {
              points(data[I, 1],
                     data[I, 2],
                     pch = symFill[
                       k], cex = CEX)
            }
          }
        }
        if(identify)
          title(paste(
                      "Spin Projection showing Classification Errors: seed/angle = ",
                      seed, "/", round(angle, 3),
                      collapse = ""), cex = 0.5)
      }
      if(!ask)
        break
    }
  }
  invisible()
}

### dots are not used but inserted to keep R from complaining
### about generic/method consistency...
### x is renamed to object for the same reason...
"summary.EMclust" <- function(object, data, G, modelNames, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  x <- object
  n <- if(is.null(dimData <- dim(data))) length(data) else dimData[1]
  hcPairs <- attr(x, "hcPairs")
  attr(hcPairs, "initialPartition") <- attr(x, "attrHC")$initialPartition
  subset <- attr(x, "subset")
  eps <- attr(x, "eps")
  tol <- attr(x, "tol")
  itmax <- attr(x, "itmax")
  equalPro <- attr(x, "equalPro")
  warnSingular <- attr(x, "warnSingular")
  class(x) <- attr(x, "args") <- attr(x, "subset") <- NULL
  attr(x, "hcPairs") <- attr(x, "attrHC") <- attr(x, "equalPro") <- NULL
  attr(x, "warnSingular") <- attr(x, "attrHC") <- NULL
  ##
  options <- c(subset = !is.null(subset), equalPro = equalPro)
  ##
  if(missing(G)) {
    G <- as.numeric(dimnames(x)[[1]])
  }
  else {
    G <- sort(G)
  }
  Glabels <- as.character(G)
  if(missing(modelNames))
    modelNames <- dimnames(x)[[2]]
  x <- x[Glabels, modelNames, drop = FALSE]
  X <- is.na(x)
  if(all(X))
    stop("none of the selected models could be fitted")
  x[X] <-  - .Machine$double.xmax
  ##
  l <- nrow(x)
  m <- ncol(x)
  best <- max(x)
  rowsBest <- (matrix(rep(1:l, m), l, m)[x == best])[1]
  colsBest <- (matrix(rep(1:m, rep(l, m)), l, m)[x == best])[1]
  namesBest <- dimnames(x[rowsBest, colsBest, drop = FALSE])
  bestG <- namesBest[[1]]
  maxG <- max(G)
  minG <- min(G)
  if(minG != maxG) {
    if(bestG == maxG) {
      warning("BIC maximized at upper limit on G")
    }
    else if(minG != 1 && bestG == min(G)) {
      warning("BIC maximized at lower limit on G")
    }
  }
  bestModel <- namesBest[[2]]
  if(min(l, m) > 1) {
    M <- modelNames[modelNames != bestModel]
    y <- x[, M]
    other <- max(y)
    otherG <- (matrix(rep(Glabels, m - 1), l, m - 1)[y == other])[1]
    otherModel <- (matrix(rep(M, rep(l, m - 1)), l, m - 1)[y == 
                                                           other])[1]
    y <- x[, bestModel]
    w <- y[y != best]
    if(length(w) == l - 1) {
      same <- max(w)
      sameG <- (Glabels[y == same])[1]
    }
    else {
      same <- best
      sameG <- (Glabels[y == same])[2]
    }
    nam1 <- paste(bestModel, bestG, sep = ",")
    nam2 <- paste(bestModel, sameG, sep = ",")
    nam3 <- paste(otherModel, otherG, sep = ",")
    bestBICs <- c(nam1 = best, nam2 = same, nam3 = other)
    names(bestBICs) <- c(nam1, nam2, nam3)
  }
  else if(l != 1) {
    ## one model, more than one number of clusters
    w <- x[x != best]
    if(length(w) == l - 1) {
      same <- max(w)
      sameG <- (Glabels[x == same])[1]
    }
    else {
      same <- best
      sameG <- (Glabels[x == same])[2]
    }
    nam1 <- paste(bestModel, bestG, sep = ",")
    nam2 <- paste(bestModel, sameG, sep = ",")
    bestBICs <- c(nam1 = best, nam2 = same)
    names(bestBICs) <- c(nam1, nam2)
  }
  else if(m != 1) {
    ## one number of clusters, more than one model
    M <- (1:m)[modelNames == bestModel]
    y <- x[,  - M]
    other <- max(y)
    otherG <- (matrix(rep(Glabels, m - 1), l, m - 1)[y == other])[
                                                       1]
    otherModel <- (matrix(rep(modelNames[ - M], rep(l, m - 1)),
                          l, m - 1)[y == other])[1]
    nam1 <- paste(bestModel, bestG, sep = ",")
    nam3 <- paste(otherModel, otherG, sep = ",")
    bestBICs <- c(best, other)
    names(bestBICs) <- c(nam1, nam3)
  }
  else {
    nam1 <- paste(bestModel, bestG, sep = ",")
    bestBICs <- best
    names(bestBICs) <- nam1
  }
  G <- as.numeric(bestG)
  if(is.null(subset)) {
    if(G == 1) {
      out <- mvn(modelName = bestModel, data = data)
      return(structure(c(list(bic = bestBICs, options = 
                              options, classification = rep(1, n), 
                              uncertainty = rep(0, n)), out),
                       class = "summary.EMclust"))
    }
    clss <- hclass(hcPairs, G)
    z <- unmap(clss)
    out <- me(modelName = bestModel, data = data, z = z, eps = eps,
              tol = tol, itmax = itmax, equalPro = equalPro, 
              warnSingular = warnSingular)
  }
  else {
    clss <- hclass(hcPairs, G)
    z <- unmap(clss)
    ms <- mstep(modelName = bestModel, data = data[subset,  ],
                z = z, eps = eps, tol = tol, itmax = itmax, equalPro = 
                equalPro, warnSingular)
    out <- do.call("em", c(list(data = data, eps = eps, tol = tol,
                                itmax = itmax, equalPro = equalPro,
                                warnSingular = warnSingular), ms))
  }

  bestBICs[bestBICs == -.Machine$double.xmax] <- NA

  structure(c(list(bic = bestBICs, options = options, classification = 
                   map(out$z), uncertainty = 1 - apply(out$z, 1, max)), out),
            class = "summary.EMclust")
}

### dots are not used but inserted to keep R from complaining
### about generic/method consistency...
### x is renamed to object for the same reason...
"summary.EMclustN" <- function(object, data, G, modelNames, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  x <- object
  n <- if(is.null(dimData <- dim(data))) length(data) else dimData[1]
  hcPairs <- attr(x, "hcPairs")
  attr(hcPairs, "initialPartition") <- attr(x, "attrHC")$initialPartition
  eps <- attr(x, "eps")
  tol <- attr(x, "tol")
  itmax <- attr(x, "itmax")
  warnSingular <- attr(x, "warnSingular")
  equalPro <- attr(x, "equalPro")
  noise <- attr(x, "noise")
  Vinv <- attr(x, "Vinv")
  
  class(x) <- attr(x, "args") <- NULL
  attr(x, "eps") <- attr(x, "tol") <- attr(x, "itmax") <- NULL
  attr(x, "noise") <- attr(x, "Vinv") <- NULL
  attr(x, "hcPairs") <- attr(x, "attrHC") <- attr(x, "equalPro") <- NULL
  attr(x, "warnSingular") <- NULL
  ##
  options <- c(noise = TRUE, equalPro = equalPro)
  ##
  if(missing(G)) {
    G <- as.numeric(dimnames(x)[[1]])
  }
  else {
    G <- sort(G)
  }
  Glabels <- as.character(G)
  if(missing(modelNames))
    modelNames <- dimnames(x)[[2]]
  x <- x[Glabels, modelNames, drop = FALSE]
  X <- is.na(x)
  if(all(X))
    stop("none of the selected models could be fitted")
  x[X] <-  - .Machine$double.xmax
  ##
  l <- nrow(x)
  m <- ncol(x)
  best <- max(x)
  rowsBest <- (matrix(rep(1:l, m), l, m)[x == best])[1]
  colsBest <- (matrix(rep(1:m, rep(l, m)), l, m)[x == best])[1]
  namesBest <- dimnames(x[rowsBest, colsBest, drop = FALSE])
  bestG <- namesBest[[1]]
  maxG <- max(G)
  minG <- min(G)
  if(minG != maxG) {
    if(bestG == maxG) {
      warning("BIC maximized at upper limit on G")
    }
    else if(bestG == min(G)) {
      warning("BIC maximized at lower limit on G")
    }
  }
  bestModel <- namesBest[[2]]
  if(min(l, m) > 1) {
    M <- modelNames[modelNames != bestModel]
    y <- x[, M]
    other <- max(y)
    otherG <- (matrix(rep(Glabels, m - 1), l, m - 1)[y == other])[1]
    otherModel <- (matrix(rep(M, rep(l, m - 1)), l, m - 1)[y == other])[1]
    y <- x[, bestModel]
    w <- y[y != best]
    if(length(w) == l - 1) {
      same <- max(w)
      sameG <- (Glabels[y == same])[1]
    }
    else {
      same <- best
      sameG <- (Glabels[y == same])[2]
    }
    nam1 <- paste(bestModel, bestG, sep = ",")
    nam2 <- paste(bestModel, sameG, sep = ",")
    nam3 <- paste(otherModel, otherG, sep = ",")
    bestBICs <- c(nam1 = best, nam2 = same, nam3 = other)
    names(bestBICs) <- c(nam1, nam2, nam3)
  }
  else if(l != 1) {
    ## one model, more than one number of clusters
    w <- x[x != best]
    if(length(w) == l - 1) {
      same <- max(w)
      sameG <- (Glabels[x == same])[1]
    }
    else {
      same <- best
      sameG <- (Glabels[x == same])[2]
    }
    nam1 <- paste(bestModel, bestG, sep = ",")
    nam2 <- paste(bestModel, sameG, sep = ",")
    bestBICs <- c(nam1 = best, nam2 = same)
    names(bestBICs) <- c(nam1, nam2)
  }
  else if(m != 1) {
    ## one number of clusters, more than one model
    M <- (1:m)[modelNames == bestModel]
    y <- x[,  - M]
    other <- max(y)
    otherG <- (matrix(rep(Glabels, m - 1), l, m - 1)[y == other])[
                                                       1]
    otherModel <- (matrix(rep(modelNames[ - M], rep(l, m - 1)),
                          l, m - 1)[y == other])[1]
    nam1 <- paste(bestModel, bestG, sep = ",")
    nam3 <- paste(otherModel, otherG, sep = ",")
    bestBICs <- c(best, other)
    names(bestBICs) <- c(nam1, nam3)
  }
  else {
    nam1 <- paste(bestModel, bestG, sep = ",")
    bestBICs <- best
    names(bestBICs) <- nam1
  }
  G <- as.numeric(bestG)
  if(G == 0) {
    return(structure(list(Vinv = Vinv, loglik = n * logb(Vinv),
                          options = options)), class = "summary.EMclustN")
  }
  clss <- hclass(hcPairs, G)
  k <- as.numeric(bestG)
  k1 <- k + 1
  z <- matrix(0, n, k1)
  z[!noise, 1:k] <- unmap(clss)
  z[!noise, k1] <- 0
  z[noise, k1] <- 1
  out <- me(modelName = bestModel, data = data, z = z[, 1:k1], eps = eps,
            tol = tol, itmax = itmax, equalPro = equalPro, warnSingular = 
            warnSingular, noise = TRUE, Vinv = Vinv)

  bestBICs[bestBICs == -.Machine$double.xmax] <- NA

  structure(c(list(bic = bestBICs, noise = TRUE, equalPro = equalPro,
                   classification = map(out$z),
                   uncertainty = 1 - apply(out$z, 1, max)), out),
            class = "summary.EMclustN")
}

### x is renamed to object to keep R from complaining
### about generic/method consistency...
"summary.Mclust" <- function(object, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  x <- object
  M <- switch(EXPR=x$model,
              E = "equal variance",
              V = "unequal variance",
              EII = "spherical, equal volume",
              VII = "spherical, varying volume",
              EEI = "diagonal, equal volume and shape",
              VEI = "diagonal, equal shape",
              EVI = "diagonal, equal volume",
              VVI = "diagonal, varying volume and shape",
              EEE = "elliposidal, equal variance",
              VEE = "elliposidal, equal shape and orientation",
              EVE = "elliposidal, equal volume and orientation",
              VVE = "ellipsoidal, equal orientation",
              EEV = "ellipsoidal, equal volume and shape",
              VEV = "ellipsoidal, equal shape",
              EVV = "elliposidal, equal volume",
              VVV = "ellipsoidal, unconstrained",
              stop("invalid model id for EM"))
  G <- length(unique(x$classification))
  cat("\n best model:", M, "with", G, "groups\n")
  if(FALSE) {
    aveUncer <- round(mean(x$uncertainty), 3)
    medUncer <- round(median(x$uncertainty), 3)
    cat("\n averge/median classification uncertainty:", aveUncer,
        "/", medUncer, "\n\n")
  }
  invisible()
}

### dots are not used but inserted to keep R from complaining
### about generic/method consistency...
### test is renamed to object for the same reason...
"summary.mclustDAtest" <- function(object, pro, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  test <- object
  clfun <- function(x)
    {
      cl <- names(x)[x == max(x)]
      if(length(cl) > 1)
        NA
      else cl
    }
  if(!missing(pro)) {
    if(length(pro) != ncol(test))
      stop("wrong number of prior probabilities")
    test <- sweep(test, MARGIN = 2, STATS = pro, FUN = "*")
  }
  cl <- apply(test, 1, clfun)
  z <- sweep(test, MARGIN = 1, STATS = apply(test, 1, sum), FUN = "/")
  attr(z, "class") <- NULL
  list(classification = cl, z = z)
}

### dots are not used but inserted to keep R from complaining
### about generic/method consistency...
### x is renamed to object for the same reason...
"summary.mclustDAtrain" <- function(object, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  x <- object
  L <- length(x)
  M <- max(unlist(lapply(x, function(y)
                         y$n)))
  N <- names(x)
  s <- rep(list(list(model = "XXX,00", classification = rep(0, M))),
           times = L)
  names(s) <- N
  for(l in 1:L) {
    s[[l]]$model <- paste(x[[l]]$modelName, x[[l]]$G, sep = ",")
    cl <- if(!is.null(x[[l]]$z)) map(x[[l]]$z) else rep(1, x[[l]]$
                                                        n)
    s[[l]]$classification <- cl
  }
  s
}

"surfacePlot" <- function(data, mu, pro, ...,
                          type = c("contour", "image", "persp"),
                          what = c("density", "uncertainty", "skip"),
                          transformation = c("none", "log", "sqrt"),
                          grid = 50, nlevels = 20, scale = FALSE,
                          identify = FALSE, verbose = FALSE, xlim,
                          ylim, swapAxes = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  data <- as.matrix(data)
  p <- ncol(data)
  if(p != 2)
    stop("for two-dimensional data only")
  densNuncer <- function(modelName, data, mu, pro, sigma)
    {
      ## ... sigmasq or sigma, pro, eps
      cden <- do.call("cdens", list(modelName = modelName, data = 
                                    data, mu = mu, sigma = sigma))
      z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = pro)
      den <- apply(z, 1, sum)
      z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
      data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
    }
  if(scale)
    par(pty = "s")
  aux <- list(...)
  sigma <- aux$sigma
  decomp <- aux$decomp
  params <- !is.null(mu) && (!is.null(sigma) || !is.null(decomp))
  if(!is.null(mu)) {
    if(is.null(sigma)) {
      if(is.null(decomp)) {
        params <- FALSE
        warning("covariance not supplied")
      }
      else {
        sigma <- decomp2sigma(decomp)
      }
    }
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      params <- FALSE
      warning("covariance improperly specified")
    }
    if(G != dimpar[3]) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  }
  if(length(grid) == 1)
    grid <- c(grid, grid)
  if(swapAxes) {
    if(params) {
      mu <- mu[2:1,  ]
      sigma <- sigma[2:1, 2:1,  ]
    }
    data <- data[, 2:1]
  }
  if(!is.null(dnames <- dimnames(data)[[2]])) {
    xlab <- dnames[1]
    ylab <- dnames[2]
  }
  else xlab <- ylab <- ""
  if(missing(xlim))
    xlim <- range(data[, 1])
  if(missing(ylim))
    ylim <- range(data[, 2])
  if(scale) {
    d <- diff(xlim) - diff(ylim)
    if(d > 0) {
      ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
    }
    else {
      xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
    }
  }
  if(!params)
    stop("need parameters to compute density")
  if(missing(pro))
    stop("need mixing proportions to compute density")
  x <- grid1(n = grid[1], range = xlim, edge = TRUE)
  y <- grid1(n = grid[2], range = ylim, edge = TRUE)
  xy <- grid2(x, y)
  if(verbose)
    cat("\n computing density and uncertainty over grid ...\n")
  Z <- densNuncer(modelName = "VVV", data = xy, mu = mu, pro = pro, sigma
                  = sigma)
  lx <- length(x)
  ly <- length(y)
  CI <- type
  DU <- what
  TRANS <- transformation
  ask <- length(CI) > 1 || length(DU) > 1 || length(TRANS) > 1
  while(TRUE) {
    if(!length(DU))
      break
    if(!length(CI))
      break
    if(!length(TRANS))
      break
    if(length(DU) > 1) {
      du <- menu(DU, title = 
                 "\n density or uncertainty? (0 to exit)\n")
      if(!du)
        return(invisible())
      du <- DU[du]
    }
    else du <- DU
    if(mode(du) != "character")
      break
    if(du == "skip") {
      plot(xy[, 1], xy[, 2], xlab = "", ylab = "", axes = FALSE,
           type = "n")
      next
    }
    if(length(CI) > 1) {
      ci <- menu(CI, title = "\n plot type? (0 to exit):\n")
      if(!ci)
        return(invisible())
      ci <- CI[ci]
    }
    else ci <- CI
    if(length(TRANS) > 1) {
      trans <- menu(TRANS, title = 
                    "\n transformation? (0 to exit):\n")
      if(!trans)
        return(invisible())
      trans <- TRANS[trans]
    }
    else trans <- TRANS
    if(mode(trans) != "character")
      break
    if(du == "density") {
      zz <- matrix(Z$density, lx, ly)
      title2 <- "Density"
    }
    else {
      zz <- matrix(Z$uncertainty, lx, ly)
      title2 <- "Uncertainty"
    }
    if(trans == "log") {
      z <- logb(zz)
      title1 <- "Log"
    }
    else if(trans == "sqrt") {
      z <- sqrt(zz)
      title1 <- "Square Root"
    }
    else {
      z <- zz
      title1 <- ""
    }
    if(ci == "contour") {
      title3 <- "Contour"
      contour(x = x, y = y, z = z, nlevels = nlevels, xlab = 
              xlab, ylab = ylab)#, labcex = 0)
    }
    else if(ci == "image") {
      title3 <- "Image"
      image(x = x, y = y, z = z, xlab = xlab, ylab = ylab)
    }
    else {
      title3 <- "Perspective"
      persp(x = x, y = y, z = z, xlab = xlab, ylab = ylab,
            theta = 60, phi = 30, expand = .6)
    }
    if(identify) {
      TITLE <- paste(c(title1, title2, title3, "Plot"), 
                     collapse = " ")
      title(TITLE, cex = 0.5)
    }
    if(!ask)
      break
  }
  invisible(list(x = x, y = y, z = z))
}

"traceW" <- 
  function(x)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## sum(as.vector(sweep(x, 2, apply(x, 2, mean)))^2)
  dimx <- dim(x)
  n <- dimx[1.]
  p <- dimx[2.]
  .Fortran("mcltrw",
           as.double(x),
           as.integer(n),
           as.integer(p),
           double(p),
           double(1.),
                   PACKAGE="mclust")[[5.]]
}

"uncerPlot" <- function(z, truth, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  parSave <- par(no.readonly=TRUE)
  par(pty = "m")
  uncer <- 1 - apply(z, 1, max)
  ord <- order(uncer)
  ##	plot(uncer[ord], pch = ".", xlab = "", ylab = "uncertainty", 
  ##      ylim = c(- (0.5/32), 0.5))
  M <- max(uncer)
  plot(uncer[ord], ylab = "uncertainty", ylim = c( - (M/32), M), xaxt = 
       "n", type = "n")
  points(uncer[ord], pch = 15, cex = 0.5)
  lines(uncer[ord])
  abline(h = c(0, 0), lty = 4)
  if(!missing(truth)) {
    n <- length(truth)
    result <- map(z)
    ERRORS <- classErrors(result, truth)
    if(any(ERRORS)) {
      I <- (1:n)[ERRORS]
      for(i in I) {
        x <- (1:n)[ord == i]
                                        #
        lines(c(x, x), c( - (0.5/32), uncer[i]), lty = 
              1)
      }
    }
  }
  par(parSave)
  invisible()
}

"unchol" <- function(x, upper)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(missing(upper)) {
    upper <- any(x[row(x) < col(x)])
    lower <- any(x[row(x) > col(x)])
    if(upper && lower)
      stop("not a triangular matrix")
    if(!(upper || lower)) {
      x <- diag(x)
      return(diag(x * x))
    }
  }
  dimx <- dim(x)
  storage.mode(x) <- "double"
  .Fortran("unchol",
           as.integer(if(upper) 1 else 0),
           x,
           as.integer(nrow(x)),
           as.integer(ncol(x)),
           integer(1),
                   PACKAGE="mclust")[[2]]
}

"unmap" <- function(classification, noise, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  ## converts a classification to conditional probabilities
  ## classes are arranged in sorted order
  ## if a noise indicator is specified, that column is placed last
  n <- length(classification)
  u <- sort(unique(classification))
  labs <- as.character(u)
  k <- length(u)
  if(!missing(noise)) {
    l <- u == noise
    if(any(l)) {
      m <- max(u) + 1
      u[l] <- m
      labs <- labs[order(u)]
      u <- sort(u)
      classification[classification == noise] <- m
    }
  }
  z <- matrix(0, n, k)
  for(j in 1:k)
    z[classification == u[j], j] <- 1
  ##
  ## z <- matrix(1, n, k)
  ## for(j in 1:k) z[classification == u[j], j] <- k + 1
  ## z <- sweep(z, 1, apply(z, 1, sum), "/")
  ##
  dimnames(z) <- list(NULL, labs)
  z
}

"vecnorm" <- function(x, p = 2)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  if(is.character(p)) {
    if(charmatch(p, "maximum", nomatch = 0) == 1)
      p <- Inf
    else if(charmatch(p, "euclidean", nomatch = 0) == 1)
      p <- 2
    else stop("improper specification of p")
  }
  if(!is.numeric(x) && !is.complex(x))
    stop("mode of x must be either numeric or complex")
  if(!is.numeric(p))
    stop("improper specification of p")
  if(p < 1)
    stop("p must be greater than or equal to 1")
  if(is.numeric(x))
    x <- abs(x)
  else x <- Mod(x)
  if(p == 2)
    return(.Fortran("d2norm",
                    as.integer(length(x)),
                    as.double(x),
                    as.integer(1),
                    double(1),
                   PACKAGE="mclust")[[4]])
#                    value = double(1))$value)
  if(p == Inf)
    return(max(x))
  if(p == 1)
    return(sum(x))
  xmax <- max(x)
  if(!xmax)
    xmax <- max(x)
  if(!xmax)
    return(xmax)
  x <- x/xmax
  xmax * sum(x^p)^(1/p)
}

"estep2" <- function(logCden, pro, ...)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  aux <- list(...)
  n <- nrow(logCden)
  ncolz <- ncol(logCden)
  if(missing(pro))
    pro <- rep(1/ncolz, ncolz)
  out <- .Fortran("estep2",
                  as.integer(n),
                  as.integer(ncolz),
                  as.double(pro),
                  as.double(logCden),
                  double(1),
                   PACKAGE="mclust")[c(4, 5)]
  list(z = matrix(out[[1]], n, ncolz), loglik = out[[2]])
}

"compareClass" <- function(a, b)
{
  if((l <- length(a)) != length(b))
    stop("unequal lengths")
  ta <- table(a)
  na <- length(ta)
  tb <- table(b)
  nb <- length(tb)
  Pa <- ta/l
  Pb <- tb/l
  Tab <- table(a, b)
  Pab <- Tab/l
  Ha <-  - sum(Pa * log(Pa))
  Hb <-  - sum(Pb * log(Pb))
  Iab <- Pab
  Iab <- sweep(Iab, MARGIN = 1, FUN = "/", STATS = Pa)
  Iab <- sweep(Iab, MARGIN = 2, FUN = "/", STATS = Pb)
  Z <- Pab == 0
  Iab <- sum(log(Iab)[!Z] * Pab[!Z])
  ((Ha + Hb) - 2 * Iab)/log(l)
}

"mapClass" <- function(a, b)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  l <- length(a)
  x <- y <- rep(NA, l)
  if(l != length(b)) {
    warning("unequal lengths")
    return(x)
  }
  aChar <- as.character(a)
  bChar <- as.character(b)
  Tab <- table(a, b)
  Ua <- dimnames(Tab)[[1]]
  Ub <- dimnames(Tab)[[2]]
  aTOb <- rep(list(Ub), length(Ua))
  names(aTOb) <- Ua
  bTOa <- rep(list(Ua), length(Ub))
  names(bTOa) <- Ub
  ## -------------------------------------------------------------
  k <- nrow(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 1, max)
  for(i in 1:k) {
    I <- match(Max[i], Tab[i,  ], nomatch = 0)
    aTOb[[i]] <- Ub[I]
  }
  if(is.numeric(b))
    aTOb <- lapply(aTOb, as.numeric)
  k <- ncol(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 2, max)
  for(j in (1:k)) {
    J <- match(Max[j], Tab[, j])
    bTOa[[j]] <- Ua[J]
  }
  if(is.numeric(a))
    bTOa <- lapply(bTOa, as.numeric)
  list(aTOb = aTOb, bTOa = bTOa)
}

"classError" <- function(classification, truth)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  sum(as.numeric(classErrors(classification, truth)))/length(truth)
}

"classErrors" <- function(classification, truth)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
# Distribution of MCLUST is prohibited except by agreement with the 
# University of Washington.
##
  q <- function(map, len, x)
    {
      x <- as.character(x)
      map <- lapply(map, as.character)
      y <- sapply(map, function(x)
                  x[1])
      best <- y != x
      if(all(len) == 1)
        return(best)
      errmin <- sum(as.numeric(best))
      z <- sapply(map, function(x)
                  x[length(x)])
      mask <- len != 1
      counter <- rep(0, length(len))
      k <- sum(as.numeric(mask))
      j <- 0
      while(y != z) {
        i <- k - j
        m <- mask[i]
        counter[m] <- (counter[m] %% len[m]) + 1
        y[x == name(map)[m]] <- map[[m]][counter[m]]
        temp <- y != x
        err <- sum(as.numeric(temp))
        if(err < errmin) {
          errmin <- err
          best <- temp
        }
        j <- (j + 1) %% k
      }
      best
    }
  MAP <- mapClass(classification, truth)
  len <- sapply(MAP[[1]], length)
  if(all(len) == 1) {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT))
    one <- CtoT[I] != truth
  }
  else {
    one <- q(MAP[[1]], len, truth)
  }
  len <- sapply(MAP[[2]], length)
  if(all(len) == 1) {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(truth), names(TtoC))
    two <- TtoC[I] != classification
  }
  else {
    two <- q(MAP[[2]], len, classification)
  }
  if(sum(as.numeric(one)) > sum(as.numeric(two)))
    as.vector(one)
  else as.vector(two)
}

