.onLoad <- function(libname, pkgname) 
{
  library.dynam("mclust", pkgname, libname)
}

#############################################################################

.mclust <- structure(list(
  emModelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV"), 
  hcModelNames = c("EII", "VII", "EEE", "VVV"), 
  bicPlotSymbols = structure(c(17, 2, 16, 10, 13, 1, 15, 12, 7, 0, 17, 2),
                            .Names = c("EII", "VII", "EEI", "EVI", "VEI", "VVI", "EEE", "EEV", "VEV", "VVV", "E", "V")), 
  bicPlotColors = structure(c("gray", "black", "orange", "brown", "red", "magenta", "forestgreen", "green", "cyan", "blue", "gray", "black"), 
                            .Names = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV", "E", "V")), 
#  classPlotSymbols = c(17, 0, 16, 4, 10, 18, 6, 7, 3, 11, 2, 12, 8, 15, 1, 9, 14, 13, 5), 
  classPlotSymbols = c(16, 0, 17, 3, 15, 4, 1, 8, 2, 7, 5, 9, 10, 11, 12, 13, 14),
  # classPlotColors = c("blue", "red", "green", "cyan", "magenta", "forestgreen", "purple", "orange", "gray", "brown", "black"), 
  classPlotColors = c("dodgerblue2", "red3", "green3", "slateblue", "orange", "skyblue1", "forestgreen", "steelblue4", "gray", "brown", "black"),
  warn = TRUE), 
  .Names = c("emModelNames", "hcModelNames", "bicPlotSymbols", "bicPlotColors", "classPlotSymbols", "classPlotColors", "warn"))

mclust.options <- function(...)
{
  current <- .mclust
  if(nargs() == 0) return(current)
  args <- list(...)
  if(length(args) == 1 && is.null(names(args))) 
    { arg <- args[[1]]
      switch(mode(arg),
             list = args <- arg,
             character = return(.mclust[[arg]]),
             stop("invalid argument: ", dQuote(arg)))
    }
  if(length(args) == 0) return(current)
  n <- names(args)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- args
  if(sys.parent() == 0) env <- asNamespace("mclust") else env <- parent.frame()
  assign(".mclust", current, envir = env)
  invisible(current)
}

#############################################################################

# This old version is inefficient for large vectors
# adjustedRandIndex <- function(x, y)
# {
#   x <- as.vector(x)
#   y <- as.vector(y)
#   xx <- outer(x, x, "==")
#   yy <- outer(y, y, "==")
#   upper <- row(xx) < col(xx)
#   xx <- xx[upper]
#   yy <- yy[upper]
#   a <- sum(as.numeric(xx & yy))
#   b <- sum(as.numeric(xx & !yy))
#   c <- sum(as.numeric(!xx & yy))
#   d <- sum(as.numeric(!xx & !yy))
#   ni <- (b + a)
#   nj <- (c + a)
#   abcd <- a + b + c + d
#   q <- (ni * nj)/abcd
#   (a - q)/((ni + nj)/2 - q)
# }

adjustedRandIndex <- function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y)) 
     stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
         ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}

bicEMtrain <- function(data, labels, modelNames=NULL) 
{
    z <- unmap(as.numeric(labels))
    G <- ncol(z)
    dimData <- dim(data)
    oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 
        1
    if (oneD || length(dimData) != 2) {
        if (is.null(modelNames)) 
            modelNames <- c("E", "V")
        if (any(!match(modelNames, c("E", "V"), nomatch = 0))) 
            stop("modelNames E or V for one-dimensional data")
    }
    else {
        if (is.null(modelNames)) 
            modelNames <- .mclust$emModelNames
    }
    BIC <- rep(NA, length(modelNames))
    names(BIC) <- modelNames
    for (m in modelNames) {
        mStep <- mstep(modelName = m, data = data, z = z, warn = FALSE)
        eStep <- do.call("estep", c(mStep, list(data=data, warn=FALSE)))
        if (is.null(attr(eStep, "warn"))) 
            BIC[m] <- do.call("bic", eStep)
    }
    BIC
}

classError <- function(classification, truth)
{
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
      y[x == names(map)[m]] <- map[[m]][counter[m]]
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
        if (any(isNA <- is.na(classification))) {
          classification <- as.character(classification)
          nachar <- paste(unique(classification[!isNA]),collapse="")
          classification[isNA] <- nachar
        }
  MAP <- mapClass(classification, truth)
  len <- sapply(MAP[[1]], length)
  if(all(len) == 1) {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT), nomatch= 0)               
    one <- CtoT[I] != truth
  }
  else {
    one <- q(MAP[[1]], len, truth)
  }
  len <- sapply(MAP[[2]], length)
  if(all(len) == 1) {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(truth), names(TtoC), nomatch = 0)
    two <- TtoC[I] != classification
  }
  else {
    two <- q(MAP[[2]], len, classification)
  }
  err <- if(sum(as.numeric(one)) > sum(as.numeric(two)))
    as.vector(one)
  else as.vector(two)
        bad <- seq(along = classification)[err]
        list(misclassified = bad, errorRate = length(bad)/length(truth))
}

cv1EMtrain <- function(data, labels, modelNames=NULL) 
{
    z <- unmap(as.numeric(labels))
    G <- ncol(z)
    dimDataset <- dim(data)
    oneD <- is.null(dimDataset) || length(dimDataset[dimDataset > 1]) == 
        1
    if (oneD || length(dimDataset) != 2) {
        if (is.null(modelNames)) 
            modelNames <- c("E", "V")
        if (any(!match(modelNames, c("E", "V"), nomatch = 0))) 
            stop("modelNames E or V for one-dimensional data")
        n <- length(data)
        cv <- matrix(1, nrow = n, ncol = length(modelNames))
        dimnames(cv) <- list(NULL, modelNames)
        for (m in modelNames) {
            for (i in 1:n) {
                mStep <- mstep(modelName = m, data = data[-i], 
                  z = z[-i, ], warn = FALSE)
                eStep <- do.call("estep", c(mStep, list(data = data[i], 
                                            warn = FALSE)))
                if (is.null(attr(eStep, "warn"))) {
                  k <- (1:G)[eStep$z == max(eStep$z)]
                  l <- (1:G)[z[i, ] == max(z[i, ])]
                  cv[i, m] <- as.numeric(!any(k == l))
                }
            }
        }
    }
    else {
        if (is.null(modelNames)) 
            modelNames <- .mclust$emModelNames
        n <- nrow(data)
        cv <- matrix(1, nrow = n, ncol = length(modelNames))
        dimnames(cv) <- list(NULL, modelNames)
        for (m in modelNames) {
            for (i in 1:n) {
                mStep <- mstep(modelName = m, data = data[-i,],
                               z = z[-i, ], warn = FALSE)
                eStep <- do.call("estep", c(mStep, list(data = data[i, 
                                       , drop = FALSE], warn = FALSE)))
                if (is.null(attr(eStep, "warn"))) {
                  k <- (1:G)[eStep$z == max(eStep$z)]
                  l <- (1:G)[z[i, ] == max(z[i, ])]
                  cv[i, m] <- as.numeric(!any(k == l))
                }
            }
        }
    }
    errorRate <- apply(cv, 2, sum)
    errorRate/n
}

# This version is bugged when a quantile is equal to the following
qclass <- function (x, k) 
{
  q <- quantile(x, seq(from = 0, to = 1, by = 1/k))
  cl <- rep(0, length(x))
  q[1] <- q[1] - 1
  for(i in 1:k) 
     cl[x > q[i] & x <= q[i+1]] <- i
  return(cl)
}
# This should correct the above bug
qclass <- function (x, k) 
{
  q <- quantile(x, seq(from = 0, to = 1, by = 1/k))
  q[1] <- q[1] - 1
  q[length(q)] <- q[length(q)] + 1
  cl <- rep(0, length(x))
  for(i in 1:k) 
    cl[x >= q[i] & x < q[i+1]] <- i
  return(cl)
}

hclass <- function(hcPairs, G)
{
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
    warning("Some selected classifications are inconsistent\n                          with mclust object"
      )
  L <- length(select)
  cl <- matrix(NA, nrow = n, ncol = L, dimnames = list(NULL, as.character(
    G)))
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

EMclust <- function(data, G = NULL, modelNames = NULL, prior = NULL, control = 
  emControl(), initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),
  Vinv = NULL, warn = FALSE, x = NULL, ...)
{
  if (!is.null(x)) {
    if (!missing(prior) || !missing(control) || 
        !missing(initialization) || !missing(Vinv)) 
 stop("only G and modelNames may be specified as arguments when x is supplied")
    prior <- attr(x,"prior") 
    control <- attr(x,"control")
    initialization <- attr(x,"initialization")
    Vinv <- attr(x,"Vinv")
    warn <- attr(x,"warn")
  }
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  if (is.null(x)) {
    if (is.null(modelNames)) {
      if (d == 1) {
        modelNames <- c("E", "V")
      }
      else {
        modelNames <- .mclust$emModelNames
        if (n <= d) {
          m <- match(c("EEE","EEV","VEV","VVV"),.mclust$emModelNames,nomatch=0)
          modelNames <- modelNames[-m]
        }
      }
    }
    if (is.null(G)) {
      G <- if (is.null(initialization$noise)) 1:9 else 0:9
    }
    else {
      G <- sort(as.numeric(G))
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if (is.null(G)) G <- Glabels
    if (is.null(modelNames)) modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if (all(Gmatch) && all(Mmatch)) {
      attr( x, "G") <- as.numeric(G)
      attr( x, "modelNames") <- modelNames
      attr( x, "returnCodes") <- 
           attr(x, "returnCodes")[as.character(G),modelNames,drop=FALSE]
      return(x[as.character(G),modelNames,drop=FALSE])
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  if (any(as.logical(as.numeric(G))) < 0) {
    if (is.null(initialization$noise)) {
      stop("G must be positive")
    }
    else {
      stop("G must be nonnegative")
    }
   }
  if (d == 1 && any(nchar(modelNames) > 1)) {
    Emodel <- any(sapply(modelNames, function(x)
    charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(modelNames, function(x)
    charmatch("V", x, nomatch = 0)[1]) == 1)
    modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  l <- length(Gall)
  m <- length(Mall)
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, 
                       dimnames = list(as.character(Gall), as.character(Mall)))
  if (!is.null(x)) {
    BIC[dimnames(x)[[1]],dimnames(x)[[2]]] <- x
    RET[dimnames(x)[[1]],dimnames(x)[[2]]] <- attr(x, "returnCodes")
    BIC <- BIC[as.character(G),modelNames,drop=FALSE]
    RET <- RET[as.character(G),modelNames,drop=FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1", ] == EMPTY]) {
         out <- mvn(modelName = mdl, data = data, prior = prior)
         BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                              n = n, d = d, G = 1, equalPro = FALSE)
         RET["1", mdl] <- attr(out, "returnCode")
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
    #######################################################
    # all data in initial hierarchical clustering phase
    #######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = "VVV", data = data)
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data)
          } 
        }
        else {
          hcPairs <- NULL 
      #   hcPairs <- hc(modelName = "E", data = data)
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs))  clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           z <- unmap(clss[, g]) 
         }
         else {
           z <- unmap(qclass( data, as.numeric(g)))
         }
         for (modelName in modelNames[BIC[g,] == EMPTY]) {
            out <- me(modelName = modelName, data = data, z = z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
        }
       }
    }
    else {
    ######################################################
    # initial hierarchical clustering phase on a subset
    ######################################################
      if (is.logical(initialization$subset)) 
        initialization$subset <- (1:n)[initialization$subset]
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
             hcPairs <- hc(modelName = "VVV", 
                   data = data[initialization$subset,  ])
          }
          else {
             hcPairs <- hc(modelName = "EII", 
                   data = data[initialization$subset,  ])
          }
         }
       else {
          hcPairs <- NULL
     #    hcPairs <- hc(modelName = "E", 
     #                  data = data[initialization$subset])
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           z <- unmap(clss[, g]) 
         }
         else {
           z <- unmap(qclass(data[initialization$subset], as.numeric(g)))
         }
         dimnames(z) <- list(as.character(initialization$subset), NULL)
         for (modelName in modelNames[!is.na(BIC[g,])]) {
            ms <- mstep(modelName = modelName, z = z, 
                        data = as.matrix(data)[initialization$subset,  ],
                        prior = prior, control = control, warn = warn)
#
#  ctrl <- control
#  ctrl$itmax[1] <- 1
#  ms <- me(modelName = modelName, data = as.matrix(data)[
#           initialization$subset,  ], z = z, prior = prior, control = ctrl)
#
            es <- do.call("estep", c(list(data = data, warn = warn), ms))
            out <- me(modelName = modelName, data = data, z = es$z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
         }
       }
     }
   }
 else {
    ######################################################
    # noise case
    ######################################################
    if (!is.null(initialization$subset)) 
      stop("subset option not implemented with noise")
    if (is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
    noise <- initialization$noise
    if (!is.logical(noise))
      noise <- as.logical(match(1:n, noise, nomatch = 0))
    if (!G[1]) {
      hood <- n * logb(Vinv)
      BIC["0",  ] <- 2 * hood - logb(n)
      if (l == 1) {
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, 
     initialization = list(hcPairs = hcPairs, subset = initialization$subset), 
                       warn = warn, n = n, d = d, oneD = oneD,
                       returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$hcPairs)) {
      if (d != 1) {
        if (n > d) {
          hcPairs <- hc(modelName = "VVV", data = data[!noise,  ])
        }
        else {
          hcPairs <- hc(modelName = "EII", data = data[!noise,  ])
        }
      }
      else {
        hcPairs <- NULL 
   #    hcPairs <- hc(modelName = "E", data = data[!noise])
      }
    }
    else hcPairs <- initialization$hcPairs
    if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
    z <- matrix(0, n, max(G) + 1)
    for (g in Glabels) {
       z[] <- 0
       k <- as.numeric(g)
       if (d > 1 || !is.null(hcPairs)) {
         z[!noise, 1:k] <- unmap(clss[, g])
       }
       else {
         z[!noise, 1:k] <- unmap(qclass(data[!noise]))
       }
       z[noise, k+1] <- 1
       K <- 1:(k+1) 
       for (modelName in modelNames[BIC[g,] == EMPTY]) {
         out <- me(modelName = modelName, data = data, z = z[, K], 
                   prior = prior, control = control, Vinv = Vinv, warn = warn)
         BIC[g, modelName] <- bic(modelName = modelName, loglik = out$loglik, 
                                  n = n, d = d, G = k, 
                                  noise = TRUE, 
                                  equalPro = control$equalPro)
         RET[g, modelName] <- attr(out, "returnCode")
       }
     }
  }
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            control = control, 
            initialization = list(hcPairs = hcPairs, 
                                  subset = initialization$subset,
                                  noise = initialization$noise), 
            Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD,
            returnCodes = RET, class = "mclustBIC")
}

mclustBIC <- function(data, G = NULL, modelNames = NULL, 
                      prior = NULL, control = emControl(), 
                      initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),  
                      Vinv = NULL, warn = FALSE, x = NULL, ...)
{
  if(!is.null(x)) 
    { if(!missing(prior) || !missing(control) || 
         !missing(initialization) || !missing(Vinv)) 
      stop("only G and modelNames may be specified as arguments when x is supplied")
      prior <- attr(x,"prior") 
      control <- attr(x,"control")
      initialization <- attr(x,"initialization")
      Vinv <- attr(x,"Vinv")
      warn <- attr(x,"warn")
    }
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  if (is.null(x)) {
    if (is.null(modelNames)) {
      if (d == 1) {
        modelNames <- c("E", "V")
      }
      else {
        modelNames <- .mclust$emModelNames
        if (n <= d) {
          m <- match(c("EEE","EEV","VEV","VVV"),.mclust$emModelNames,nomatch=0)
          modelNames <- modelNames[-m]
        }
      }
    }
    if (is.null(G)) {
      G <- if (is.null(initialization$noise)) 1:9 else 0:9
    }
    else {
      G <- sort(as.integer(unique(G)))
    }
    if (is.null(initialization$noise)) {
      if (any(G > n)) G <- G[G <= n]
    }
    else {
      noise <- initialization$noise
      if (!is.logical(noise)) {
        if (any(match(noise, 1:n, nomatch = 0) == 0))
          stop("numeric noise must correspond to row indexes of data")
        noise <- as.logical(match(1:n, noise, nomatch = 0))
      }
      initialization$noise <- noise
      nnoise <- sum(as.numeric(noise))
      if (any(G > (n-nnoise))) G <- G[G <= n-nnoise]
    }
    if (!is.null(initialization$subset)) {
      subset <- initialization$subset
      if (is.logical(subset)) subset <- which(subset)
      n <- length(subset)
      if (any(G > n)) G <- G[G <= n]
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if (is.null(G)) G <- Glabels
    if (is.null(modelNames)) modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if (all(Gmatch) && all(Mmatch)) {
      attr( x, "G") <- as.numeric(G)
      attr( x, "modelNames") <- modelNames
      attr( x, "returnCodes") <- 
           attr(x, "returnCodes")[as.character(G),modelNames,drop=FALSE]
      return(x[as.character(G),modelNames,drop=FALSE])
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  if (any(as.logical(as.numeric(G))) < 0) {
    if (is.null(initialization$noise)) {
      stop("G must be positive")
    }
    else {
      stop("G must be nonnegative")
    }
   }
  if (d == 1 && any(nchar(modelNames) > 1)) {
    Emodel <- any(sapply(modelNames, function(x)
    charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(modelNames, function(x)
    charmatch("V", x, nomatch = 0)[1]) == 1)
    modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  l <- length(Gall)
  m <- length(Mall)
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, 
                       dimnames = list(as.character(Gall), as.character(Mall)))
  if (!is.null(x)) {
    BIC[dimnames(x)[[1]],dimnames(x)[[2]]] <- x
    RET[dimnames(x)[[1]],dimnames(x)[[2]]] <- attr(x, "returnCodes")
    BIC <- BIC[as.character(G),modelNames,drop=FALSE]
    RET <- RET[as.character(G),modelNames,drop=FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1", ] == EMPTY]) {
         out <- mvn(modelName = mdl, data = data, prior = prior)
         BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                              n = n, d = d, G = 1, equalPro = FALSE)
         RET["1", mdl] <- attr(out, "returnCode")
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
    #######################################################
    # all data in initial hierarchical clustering phase
    #######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = "VVV", data = data)
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data)
          } 
        }
        else {
          hcPairs <- NULL 
      #   hcPairs <- hc(modelName = "E", data = data)
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs))  clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           cl <- clss[,g]
         }
         else {
           cl <- qclass( data, as.numeric(g))
         }
         z <- unmap(cl, groups = 1:max(cl))
         if (any(apply( z, 2, max) == 0)) {
#  missing groups
           warning("there are missing groups")    
           small <- sqrt(.Machine$double.neg.eps)
           z[z < small] <- small
           z <-  t(apply( z, 1, function(x) x/sum(x)))
         }
         for (modelName in modelNames[BIC[g,] == EMPTY]) {
            out <- me(modelName = modelName, data = data, z = z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
        }
       }
    }
    else {
      subset <- initialization$subset
      if (is.logical(subset)) subset <- which(subset)
    ######################################################
    # initial hierarchical clustering phase on a subset
    ######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
             hcPairs <- hc(modelName = "VVV", 
                   data = data[subset,  ])
          }
          else {
             hcPairs <- hc(modelName = "EII", 
                   data = data[subset,  ])
          }
         }
       else {
          hcPairs <- NULL
     #    hcPairs <- hc(modelName = "E", 
     #                  data = data[subset])
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           cl <- clss[, g]
         }
         else {
           cl <- qclass(data[subset], as.numeric(g))
         }
         z <- unmap(cl, groups = 1:max(cl))
         if (any(apply( z, 2, max) == 0)) {
#  missing groups
           warning("there are missing groups")         
           small <- sqrt(.Machine$double.neg.eps)
           z[z < small] <- small
           z <-  t(apply( z, 1, function(x) x/sum(x)))
         }
         for (modelName in modelNames[!is.na(BIC[g,])]) {
            ms <- mstep(modelName = modelName, z = z, 
                        data = as.matrix(data)[initialization$subset,  ],
                        prior = prior, control = control, warn = warn)
#
#  ctrl <- control
#  ctrl$itmax[1] <- 1
#  ms <- me(modelName = modelName, data = as.matrix(data)[
#           initialization$subset,  ], z = z, prior = prior, control = ctrl)
#
            es <- do.call("estep", c(list(data = data, warn = warn), ms))
            out <- me(modelName = modelName, data = data, z = es$z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
         }
       }
     }
   }
 else {
    ######################################################
    # noise case
    ######################################################
    noise <- initialization$noise
    if (!is.null(initialization$subset)) 
      stop("subset option not implemented with noise")
    if (is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
    if (!G[1]) {
      hood <- n * logb(Vinv)
      BIC["0",  ] <- 2 * hood - logb(n)
      if (l == 1) {
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, 
     initialization = list(hcPairs = hcPairs, subset = initialization$subset), 
                       warn = warn, n = n, d = d, oneD = oneD,
                       returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$hcPairs)) {
      if (d != 1) {
        if (n > d) {
          hcPairs <- hc(modelName = "VVV", data = data[!noise,  ])
        }
        else {
          hcPairs <- hc(modelName = "EII", data = data[!noise,  ])
        }
      }
      else {
        hcPairs <- NULL 
   #    hcPairs <- hc(modelName = "E", data = data[!noise])
      }
    }
    else hcPairs <- initialization$hcPairs
    if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
    z <- matrix(0, n, max(G) + 1)
    for (g in Glabels) {
       z[] <- 0
       k <- as.numeric(g)
       if (d > 1 || !is.null(hcPairs)) {
         cl <- clss[, g]
       }
       else {
         cl <- qclass(data[!noise], k = k)
       }
       z[!noise,1:k] <- unmap(cl, groups = 1:max(cl))
       if (any(apply( z[!noise,1:k,drop=FALSE], 2, max) == 0)) {
#           missing groups
          warning("there are missing groups")         
          z[!noise,1:k] <- max( z[!noise,1:k], sqrt(.Machine$double.neg.eps))
          z[!noise,1:k] <- apply( z[!noise,1:k,drop=FALSE], 1, function(z) z/sum(z))
       }
       z[noise, k+1] <- 1
       K <- 1:(k+1) 
       for (modelName in modelNames[BIC[g,] == EMPTY]) {
         out <- me(modelName = modelName, data = data, z = z[, K], 
                   prior = prior, control = control, Vinv = Vinv, warn = warn)
         BIC[g, modelName] <- bic(modelName = modelName, loglik = out$loglik, 
                                  n = n, d = d, G = k, 
                                  noise = TRUE, 
                                  equalPro = control$equalPro)
         RET[g, modelName] <- attr(out, "returnCode")
       }
     }
  }
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            control = control, 
            initialization = list(hcPairs = hcPairs, 
                                  subset = initialization$subset,
                                  noise = initialization$noise), 
            Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD,
            returnCodes = RET, class = "mclustBIC")
}

mclustModel <- function(data, BICvalues, G=NULL, modelNames=NULL, ...)
{
  mc <- match.call(expand.dots = FALSE)
  if (is.null(attr(BICvalues,"initialization")$noise)) {
    mc[[1]] <- as.name("summaryMclustBIC")
  }
  else {
    mc[[1]] <- as.name("summaryMclustBICn")
  }
  nm <- names(mc)
  mc[1:3] <- mc[c(1,3,2)]
  nm[1:3] <- nm[c(1,3,2)]
  nm[nm == "BICvalues"] <- "object" 
  names(mc) <- nm
  ans <- eval(mc, parent.frame())
  ans$classification <- ans$uncertainty <- NULL
  attr( ans, "bestBICvalues") <- NULL
  attr( ans, "prior") <- NULL
  attr( ans, "control") <- NULL
  attr( ans, "initialization") <- NULL
  oldClass(ans) <- "mclustModel"
  ans
}

Mclust <- function (data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = NULL, warn = FALSE, ...) 
{
    call <- match.call()
    mc <- match.call(expand.dots = FALSE)
    mc[[1]] <- as.name("mclustBIC")
    mc[[2]] <- data
    Bic <- eval(mc, parent.frame())
    G <- attr(Bic, "G")
    modelNames <- attr(Bic, "modelNames")
    Sumry <- summary(Bic, data, G = G, modelNames = modelNames)
    if (!(length(G) == 1)) {
        bestG <- length(unique(Sumry$cl))
        if (bestG == max(G)) 
            warning("optimal number of clusters occurs at max choice")
        else if (bestG == min(G)) 
            warning("optimal number of clusters occurs at min choice")
    }
    attr(Bic, "n") <- attr(Bic, "warn") <- NULL
    attr(Bic, "initialization") <- attr(Bic, "control") <- NULL
    attr(Bic, "d") <- attr(Bic, "returnCodes") <- attr(Bic, "class") <- NULL
    oldClass(Sumry) <- NULL
    Sumry$bic <- Sumry$bic[1]
    df <- (2*Sumry$loglik - Sumry$bic)/log(Sumry$n)
    ans <- c(list(call = call, BIC = Bic, df = df), Sumry)
    orderedNames <- c("call", "modelName", "n", "d", "G", 
                       "BIC", "bic", "loglik", "df", "parameters", 
                      "classification", "uncertainty")
    structure(if(Sumry$G > 1) ans[c(orderedNames,"z")] 
              else            ans[orderedNames],
              class = "Mclust")
}

predict.Mclust <- function(object, newdata, ...)
{
  if(!inherits(object, "Mclust")) 
    stop("object not of class \"Mclust\"")
  if(missing(newdata))
    { newdata <- eval.parent(object$call$data) }
  prior <- object$parameters$pro
  z <- do.call("cdens", c(list(data = newdata), object))
  z <- sweep(z, MARGIN = 1, FUN = "/", STATS = apply(z, 1, max))
  z <- sweep(z, MARGIN = 2, FUN = "*", STATS = prior/sum(prior))
  z <- sweep(z, MARGIN = 1, STATS = apply(z, 1, sum), FUN = "/")
  cl <- apply(z, 1, which.max)
  out <- list(classification = cl, z = z)
  return(out) 
}

pickBIC <- function(x, k = 3)
{
  Glabels <- dimnames(x)[[1]]
  modelNames <- dimnames(x)[[2]]
  mis <- is.na(x)
  if (all(mis)) {
    warning("none of the selected models could be fitted")
    return(rep(NA,k))
  }
  x[mis] <-  - .Machine$double.xmax
  x <- data.frame(as.vector(x), Glabels[as.vector(row(x))], 
                                modelNames[as.vector(col(x))])
  x <- x[rev(order(x[,1])),]
  namesx <- apply(x[,-1,drop = FALSE], 1, function(z) 
                  paste(as.character(z[2]), as.character(z[1]), sep = ","))
  k <- min(k, nrow(x))
  x <- x[1:k,1]
  x[x ==  - .Machine$double.xmax] <- NA
  namesx <- namesx[1:k]
  namesx[is.na(x)] <- " "
  names(x) <- namesx
  x
}

print.mclustBIC <- function(x, pick = 3, ...)
{
  subset <- !is.null(attr(x, "subset"))
  oldClass(x) <- attr(x, "args") <- NULL
  attr(x, "control") <- attr(x, "initialization") <- NULL
  attr(x, "oneD") <- attr(x, "warn") <- attr(x, "Vinv") <- NULL
  attr(x, "prior") <- attr(x, "G") <- attr(x, "modelNames") <- NULL
  ret <- attr(x, "returnCodes") == -3
  n <- attr(x, "n")
  d <- attr(x, "d")
  attr(x, "returnCodes") <- attr(x, "n") <- attr(x, "d") <- NULL

  cat("\nBIC:\n")
  NextMethod("print")
  cat("\n")
  cat("Top", pick, "models based on the BIC criterion:\n")
  print(pickBIC(x, pick), ...)
  invisible()
}

mclustModelNames <- function(model)
{
  type <- switch(EXPR = as.character(model),
                 E = "univariate, equal variance",
                 V = "univariate, unequal variance",
                 EII = "spherical, equal volume",
                 VII = "spherical, varying volume",
                 EEI = "diagonal, equal volume and shape",
                 VEI = "diagonal, equal shape",
                 EVI = "diagonal, equal volume, varying shape",
                 VVI = "diagonal, varying volume and shape",
                 EEE = "elliposidal, equal volume, shape and orientation",
                 VEE = "elliposidal, equal shape and orientation",
                 EVE = "elliposidal, equal volume and orientation",
                 VVE = "ellipsoidal, equal orientation",
                 EEV = "ellipsoidal, equal volume and shape",
                 VEV = "ellipsoidal, equal shape",
                 EVV = "elliposidal, equal volume",
                 VVV = "ellipsoidal, varying volume, shape, and orientation",
                 X   = "univariate normal",
                 XII = "spherical multivariate normal",
                 XXI = "diagonal multivariate normal",
                 XXX = "elliposidal multivariate normal",
                 warning("invalid model"))
  return(list(model = model, type = type))
}

print.Mclust <- function(x, digits = getOption("digits"), ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  M <- mclustModelNames(x$model)$type
  G <- length(unique(x$classification))
  cat(" best model: ", M, " (", x$model, ") with ", G, " components\n", sep = "")
  invisible()
}

summary.Mclust <- function(object, parameters = FALSE, classification = FALSE, ...)
{
  # collect info
  G  <- object$G
  pro <- object$parameters$pro
  if(is.null(pro)) pro <- 1
  names(pro) <- 1:G
  mean <- object$parameters$mean
  if(object$d > 1)
    { sigma <- object$parameters$variance$sigma }
  else
    { sigma <- rep(object$parameters$variance$sigmasq, object$G)[1:object$G]
      names(sigma) <- names(mean) }
  if(is.null(object$density))
    title <- paste("Gaussian finite mixture model fitted by EM algorithm")
  else
    title <- paste("Density estimation via Gaussian finite mixture modeling")
  #
  obj <- list(title = title, n = object$n, d = object$d, 
              G = G, modelName = object$modelName, 
              loglik = object$loglik, df = object$df, 
              bic = object$bic, icl = icl(object),
              pro = pro, mean = mean, variance = sigma, 
              prior = attr(object$BIC, "prior"), 
              classification = object$classification, 
              printParameters = parameters, 
              printClassification = classification)
  class(obj) <- "summary.Mclust"
  return(obj)
}

print.summary.Mclust <- function(x, digits = getOption("digits"), ...)
{
  cat(rep("-", nchar(x$title)),"\n",sep="")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)),"\n",sep="")
  #
  if(!is.null(x$prior))
    { cat("\nPrior: ")
      cat(x$prior$functionName, "(", 
            paste(names(x$prior[-1]), x$prior[-1], sep = " = ", 
                  collapse = ", "), ")", sep = "")
     cat("\n")
    }
  #
  cat("\nMclust ", x$modelName, " (", mclustModelNames(x$modelName)$type, 
      ") model with ", x$G, ifelse(x$G > 1, " components:", " component:"), "\n\n",
      sep="")
  tab <- data.frame("log-likelihood" = x$loglik, "n" = x$n, 
                    "df" = x$df, "BIC" = x$bic, "ICL" = x$icl, 
                    row.names = "")
  print(tab, digits = digits)
  #
  cat("\nClustering table:")
  print(table(x$classification), digits = digits)
  #
  if(x$printParameters)
    { cat("\nMixing probabilities:\n")
      print(x$pro, digits = digits)
      cat("\nMeans:\n")
      print(x$mean, digits = digits)
      cat("\nVariances:\n")
      if(x$d > 1) 
        { for(g in 1:x$G)
             { cat("[,,", g, "]\n", sep = "")
               print(x$variance[,,g], digits = digits) }
        }
      else print(x$variance, digits = digits)          
    }
  if(x$printClassification)
    { cat("\nClassification:\n")
      print(x$classification, digits = digits)
    }
  #
  invisible(x)
}


print.summary.mclustBIC <- function(x, digits = getOption("digits"), ...)
{
  cat("\nclassification table:")
  print(table(x$classification), ...)
  #----------------------------------------------------------------------
  #  cat("\nuncertainty (quartiles):\n")
  #  print(quantile(x$uncertainty), digits = digits, ...)
  #----------------------------------------------------------------------
  bic <- attr(x,"bestBICvalues")
  l <- length(bic)
  if(l == 1) {
          cat("\nBIC value:\n")
      print(round(bic, digits))
        }
        else {
          cat("\nbest BIC values:\n")
      print(round(bic, digits))
        }
        if(is.null(x$model)) {
          M <- "noise"
        }
        else {
          M <- mclustModelNames(x$model)$type
        }
##  cat("\nbest model:", M, "\n\n")
  ##
  ##  print(x$options)
  invisible()
}

summaryMclustBICn <- function(object, data, G=NULL, modelNames=NULL, ...)
{
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  initialization <- attr(object, "initialization")
  hcPairs <- initialization$hcPairs
  noise <- initialization$noise
  if(!is.logical(noise))
    noise <- as.logical(match(1:n, noise, nomatch = 0))
  prior <- attr(object, "prior")
  control <- attr(object, "control")
  warn <- attr(object, "warn")
  Vinv <- attr(object, "Vinv")
  oldClass(object) <- NULL
  attr(object, "control") <- attr(object, "initialization") <- NULL
  attr(object, "prior") <- attr(object, "Vinv") <- NULL
  attr(object, "warn") <- NULL
  ##
  if (is.null(G)) G <- dimnames(object)[[1]]
  if (is.null(modelNames)) modelNames <- dimnames(object)[[2]]
  bestBICs <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
  if (all(is.na(bestBICs))) {
    return(structure(NULL, bestBICvalues = bestBICs, prior = prior, control
      = control, initialization = initialization, class = "summary.mclustBIC"))
  }
  temp <- unlist(strsplit(names(bestBICs)[1], ","))
  bestModel <- temp[1] 
  G <- as.numeric(temp[2])
  if(G == 0) {
    ans <- list(bic = bestBICs, classification = rep(0, n), 
                uncertainty = rep(0, n), n = n, d = ncol(data), 
                G = 0, loglik = n * logb(Vinv), Vinv = Vinv)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "Vinv", 
                      "classification", "uncertainty")
    return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                     prior = prior, control = control, 
                     initialization = initialization, 
                     class = "summary.mclustBIC"))
  }
  G1 <- G + 1
  z <- matrix(0, n, G1)
  if(d > 1 || !is.null(hcPairs)) {
    z[!noise, 1:G] <- unmap(hclass(hcPairs, G))
  }
  else {
    z[!noise, 1:G] <- unmap(qclass(data[!noise], G))
  }
  z[noise, G1] <- 1
  out <- me(modelName = bestModel, data = data, z = z, prior = prior, 
    control = control, warn = warn, Vinv = Vinv)
  obsNames <- if (is.null(dim(data))) {
                names(data)
              }  
              else {
                dimnames(data)[[1]]
              }
  classification <- map(out$z)
  # print(G1)
  classification[classification == G1] <- 0
  uncertainty <- 1 - apply(out$z, 1, max)
  names(classification) <- names(uncertainty) <- obsNames
  ans <- c(list(bic = as.vector(bestBICs[1]), classification = classification, 
    uncertainty = uncertainty, Vinv = Vinv), out)
  orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", 
    "z", "Vinv", "classification", "uncertainty")
  structure(ans[orderedNames], 
            bestBICvalues = bestBICs, 
            prior = prior, control = control, 
            initialization = initialization, 
            class = "summary.mclustBIC")
}

summary.mclustBIC <- function(object, dataset, G, modelNames, ...)
{
  mc <- match.call(expand.dots = FALSE)
  if (is.null(attr(object,"initialization")$noise)) {
    mc[[1]] <- as.name("summaryMclustBIC")
  }
  else {
    mc[[1]] <- as.name("summaryMclustBICn")
  }
  ans <- eval(mc, parent.frame())
  Glabels <- dimnames(object)[[1]]
  if (length(Glabels) != 1 && (!missing(G) && length(G) > 1)) {
   Grange <- range(as.numeric(Glabels))
   if (match(ans$G, Grange, nomatch = 0))
  warning("best model occurs at the min or max # of components considered")
  }
  ans
}

summaryMclustBIC <- function (object, data, G = NULL, modelNames = NULL, ...) 
{
    dimData <- dim(data)
    oneD <- (is.null(dimData) || length(dimData[dimData > 1]) == 1)
    if (!oneD && length(dimData) != 2) 
        stop("data must be a vector or a matrix")
    if (oneD) {
        data <- drop(as.matrix(data))
        n <- length(data)
        d <- 1
    }
    else {
        data <- as.matrix(data)
        n <- nrow(data)
        d <- ncol(data)
    }
    initialization <- attr(object, "initialization")
    hcPairs <- initialization$hcPairs
    subset <- initialization$subset
    prior <- attr(object, "prior")
    control <- attr(object, "control")
    warn <- attr(object, "warn")
    oldClass(object) <- NULL
    attr(object, "prior") <- attr(object, "warn") <- NULL
    attr(object, "modelNames") <- attr(object, "oneD") <- NULL
    attr(object, "initialization") <- attr(object, "control") <- NULL
    d <- if (is.null(dim(data))) 1 else ncol(data)
    if (is.null(G)) 
        G <- dimnames(object)[[1]]
    if (is.null(modelNames)) 
        modelNames <- dimnames(object)[[2]]
    bestBICs <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
    if(all(is.na(bestBICs))) 
      {
        return(structure(NULL, bestBICvalues = bestBICs, prior = prior, 
                         control = control, initialization = initialization, 
                         class = "summary.mclustBIC"))
      }
    temp <- unlist(strsplit(names(bestBICs)[1], ","))
    bestModel <- temp[1]
    G <- as.numeric(temp[2])
    if(G == 1) 
      {
        out <- mvn(modelName = bestModel, data = data, prior = prior)
        ans <- c(list(bic = bestBICs, classification = rep(1, n), 
                      uncertainty = rep(0, n)), out)
        orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", 
                          "parameters", "classification", "uncertainty")
        return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                         prior = prior, control = control, initialization = initialization, 
                         class = "summary.mclustBIC"))
    }
    if(is.null(subset)) 
      {
        if(d > 1 || !is.null(hcPairs))
          { z <- unmap(hclass(hcPairs, G)) }
        else 
          { z <- unmap(qclass(data, G), groups = 1:G) }
        out <- me(modelName = bestModel, data = data, z = z, 
                  prior = prior, control = control, warn = warn)
    }
    else 
      {
        if(d > 1 || !is.null(hcPairs)) 
          { z <- unmap(hclass(hcPairs, G)) }
        else 
          { z <- unmap(qclass(data[subset], G)) }
        ms <- mstep(modelName = bestModel, prior = prior, z = z, 
                    data = as.matrix(data)[subset, ], control = control, 
                    warn = warn)
        es <- do.call("estep", c(list(data = data), ms))
        out <- me(modelName = bestModel, data = data, z = es$z, 
                  prior = prior, control = control, warn = warn)
      }
    obsNames <- if (is.null(dim(data))) names(data) else dimnames(data)[[1]]
    classification <- map(out$z)
    uncertainty <- 1 - apply(out$z, 1, max)
    names(classification) <- names(uncertainty) <- obsNames
    ans <- c(list(bic = as.vector(bestBICs[1]), 
                  classification = classification, 
                  uncertainty = uncertainty), 
             out)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", 
                      "parameters", "z", "classification", "uncertainty")
    structure(ans[orderedNames], bestBICvalues = bestBICs, prior = prior, 
              control = control, initialization = initialization, 
              class = "summary.mclustBIC")
}

defaultPrior <- function(data, G, modelName, ...)
{
  aux <- list(...)
  if(is.null(aux$shrinkage)) {
    shrinkage <- 0.01
  }
  else if(is.na(aux$shrinkage) || !aux$shrinkage) {
    shrinkage <- 0
  }
  else if(aux$shrinkage < 0) {
    stop("negative value given for shrinkage")
  }
  else {
    shrinkage <- aux$shrinkage
  }
  if(is.null(aux$mean)) {
    mean <- if (is.null(dim(data))) 
                          mean(data) else colMeans(data)
  }
  else if(any(is.na(aux$mean))) {
    if(shrinkage)
      stop("positive shrinkage with no prior mean specified")
    mean <- if (is.null(dim(data))) 
                          mean(data) else colMeans(data)
  }
  else {
    if(!shrinkage)
      stop("prior mean specified but not shrinkage")
    mean <- aux$mean
  }
  switch(EXPR = modelName,
    E = ,
    V = ,
    X = {
      dof <- 3
      if(is.null(aux$scale)) {
        scale <- var(data)/G^2
      }
      else {
        scale <- aux$scale
      }
      list(shrinkage = shrinkage, mean = mean, dof = dof,
        scale = scale)
    }
    ,
    EII = ,
    VII = ,
    XII = ,
    EEI = ,
    EVI = ,
    VEI = ,
    VVI = ,
    XXI = {
      n <- nrow(data)
      p <- ncol(data)
      dof <- p + 2
      if(is.null(aux$scale)) {
        fac <- (1/G)^(2/p)
        scale <- (fac * sum(apply(data, 2, var)))/
          p
      }
      else {
        scale <- aux$scale
      }
      list(shrinkage = shrinkage, mean = mean, dof = dof,
        scale = scale)
    }
    ,
    EEE = ,
    EEV = ,
    VEV = ,
    VVV = ,
    XXX = {
      n <- nrow(data)
      p <- ncol(data)
      dof <- p + 2
      if(is.null(aux$scale)) {
        fac <- (1/G)^(2/p)
        if(n > p) {
          scale <- fac * var(data)
        }
        else {
          scale <- fac * diag(apply(data,
            2, var))
        }
      }
      else {
        scale <- aux$scale
      }
      list(shrinkage = shrinkage, mean = mean, dof = dof,
        scale = scale)
    }
    ,
    stop("no default prior for this model"))
}

emControl <- function(eps = .Machine$double.eps, tol = c(1.0e-05, sqrt(.Machine$double.eps)), itmax = c(.Machine$integer.max, .Machine$integer.max), equalPro = FALSE)
{
##
# argList <- list(eps=eps, tol=tol, itmax=itmax, equalPro=equalPro)
# nullArgs <- sapply(argList, is.null)
# argList[nullArgs] <- .mclust[names(nullArgs[nullArgs])]
# argList$itmax[argList$itmax == Inf] <- .Machine$integer.max
# argList
##
  if(any(eps < 0)) stop("eps is negative")
  if(any(eps >= 1))
    stop("eps is not less than 1")
  if(any(tol < 0))
    stop("tol is negative")
  if(any(tol >= 1))
    stop("tol is not less than 1")
  if(any(itmax < 0))
    stop("itmax is negative")
  if(length(tol) == 1)
    tol <- rep(tol, 2)
  if(length(itmax) == 1)
    itmax <- c(itmax, .Machine$integer.max)
  i <- is.infinite(itmax)
  if(any(i))
    itmax[i] <- .Machine$integer.max
  list(eps = eps, tol = tol, itmax = itmax, equalPro = equalPro)
}

priorControl <- function(functionName = "defaultPrior", ...)
{
  c(list(functionName = functionName), list(...))
}

cdensEEE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) > 2)
    stop("data must be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
      any(is.null(parameters[c("pro", "mean", "variance")]))) 
    { WARNING <- "parameters are missing"
      if (warn) warning(WARNING)
      z <- matrix(NA,n,G)
      dimnames(z) <- list(dimnames(data)[[1]], NULL)
      return(structure(z, logarithm = logarithm, modelName = "EEE", 
                       WARNING = WARNING, returnCode = 9))
    }
  if(is.null(parameters$variance$cholSigma))
     stop("variance parameters are missing")
  temp <- .Fortran("eseee",
    as.logical(1),
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$cholSigma),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(p),
    double(1),
    double(n * G),
    PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, G)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      if (warn) warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      if (warn) warning(WARNING)
    }
    z[] <- NA
    ret <- -9 
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
        if (!logarithm) z <- exp(z)
        ret <- 0
       }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EEE",
            WARNING = WARNING, retrunCode = ret)
}

emEEE <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepEEE(data, parameters = parameters, warn = warn)$z  
  meEEE(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepEEE <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) > 2)
    stop("data must be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EEE", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholSigma))
           stop("variance parameters are missing")
  temp <- .Fortran("eseee",
    as.logical(1),
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$cholSigma),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(p),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, K)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      warning(WARNING)
                        ret <- -4 
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      warning(WARNING)
                        ret <- -5
    }
    z[] <- loglik <- NA
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EEE", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

hcEEE <- function(data, partition, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
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

## R 2.12.0: 32 bit Windows build fails due to compiler bug
## workaround: removal (hopefully temporary) of hc functionality for EEE

# Luca: commente the next line and uncommented below
#  stop("hc for EEE model is not currently supported")

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
   PACKAGE ="mclust")[c(1, 7:10)]
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
        trace <- temp[[1]][, 2]
  structure(tree,  initialPartition = partition, 
                  dimensions = dimdat, modelName = "EEE", 
                  call = match.call())
}

meEEE <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
        if (!is.null(Vinv)) {
                G <- K - 1
                if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
        }
        else G <- K
  if(all(is.na(z))) {
              WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEE", d = p, G = G, 
                     Sigma = matrix(NA, p, p), cholSigma = matrix(NA, p, p)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EEE", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeee",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(p * p),
      double(K),
      double(p),
                        PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EEE"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("meeeep",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(p * p),
      double(K),
      double(p),
                        PACKAGE = "mclust")[c(11:17, 10)]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
     cholSigma <- matrix(temp[[6]], p, p)
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
                Sigma <- matrix( NA, p, p)
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn)
        warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn)
        warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- logprior <- NA
    sigma <- array(NA, c(p, p, G))
                Sigma <- matrix(NA, p, p)
    ret <- if(control$equalPro) -2 else -3
  }
  else {
             Sigma <- unchol(cholSigma, upper = TRUE)
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- dimnames(cholSigma) <- 
                       list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EEE", d = p, G = G,
                         sigma = sigma, Sigma = Sigma, cholSigma = cholSigma) 
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "EEE", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepEEE <- function(data, z, prior = NULL,  warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEE", d = p, G = G, 
                      sigma <- array(NA, c(p,p, G)), 
                      Sigma = matrix(NA, p, p), cholSigma = matrix(NA, p, p)) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EEE", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
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
                        PACKAGE = "mclust")[7:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EEE"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("mseeep",
                        as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$scale) else priorParams$scale),
      as.double(priorParams$dof),
      double(p),
      double(p * G),
      double(p * p),
      double(G),
                        PACKAGE = "mclust")[11:13]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholSigma <- matrix(temp[[2]], p, p)
  pro <- temp[[3]]
  sigma <- array(0, c(p, p, G))
  Sigma <- unchol(cholSigma, upper = TRUE)
  for(k in 1:G)
    sigma[,  , k] <- Sigma
  WARNING <- NULL
  if(any(mu > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- sigma[] <- Sigma[] <- cholSigma[] <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- dimnames(cholSigma) <- 
               list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EEE", d = p, G = G, 
                         sigma = sigma, Sigma = Sigma, cholSigma= cholSigma)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEE", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

simEEE <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
    parameters[c("mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "EEE"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  if(is.null(cholSigma <- parameters$variance$cholSigma)) {
    if(is.null(Sigma <- parameters$variance$Sigma)) {
      stop("variance parameters must inlcude either Sigma or cholSigma"
        )
    }
    cholSigma <- chol(Sigma)
  }
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
      ncol = d) %*% cholSigma, MARGIN = 2, STATS = mu[, k],
      FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EEE")
}

cdensEEI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EEI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
  temp <- .Fortran("eseei",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[10:11]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1 
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(z, logarithm = logarithm, modelName = "EEI",
                  WARNING = WARNING, returnCode = ret)
}

cdensEII <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EII", 
                                 WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
  if(sigmasq < 0)
    stop("sigma-squared is negative")
  if(!sigmasq) {
    WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EII", 
                                 WARNING = WARNING, returnCode = 9))
  }
  temp <- .Fortran("eseii",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[9:10]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EII",
                  WARNING = WARNING, returnCode = ret)
}

emEEI <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepEEI(data, parameters = parameters, warn = warn)$z  
  meEEI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepEEI <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
                Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  } 
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EEI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
       is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
  temp <- .Fortran("eseei",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EEI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

meEEI <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEI", d = p, G = G, 
                                 scale = NA, shape = rep(NA,p)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeei",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(1),
      double(p),
      double(K),
                        PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EEI"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("meeeip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(1),
      double(p),
      double(K),
                        PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[6]]
  shape <- temp[[7]]
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(NA, p, p)
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn)
        warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn)
        warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(NA, p, p)
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(scale * shape)
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
  variance <- list(modelName = "EEI", d = p, G = G, 
                         sigma = sigma, Sigma = Sigma, 
                         scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
      info = info, WARNING = WARNING, returnCode = ret)
}

mstepEEI <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEI", d = p, G = G, 
                                 scale = NA, shape = rep(NA,p)) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))

  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("mseei",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      double(p * G),
      double(1),
      double(p),
      double(G),
                        PACKAGE = "mclust")[6:9]
  }
  else {
    storage.mode(z) <- "double"
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EEI"), prior[names(
      prior) != "functionName"]))
    temp <- .Fortran("mseeip",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p * G),
      double(1),
      double(p),
      double(G),
                        PACKAGE = "mclust")[10:13]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[2]]
  shape <- temp[[3]]
  pro <- temp[[4]]
  WARNING <- NULL
  if(any(c(shape, scale) > signif(.Machine$double.xmax, 6)) || any(!c(
    scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn)
      warning(WARNING)
    mu[] <- pro[] <- scale <- shape[] <- NA
    sigma <- Sigma <- array(NA, c(p, p, G))
                ret <- -1
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(scale * shape)
    for(k in 1:G)
      sigma[,  , k] <- Sigma
                ret <- 0
  }
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EEI", d = p, G = G, 
                         sigma = sigma, Sigma = Sigma, 
                         scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

simEEI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
    parameters[c("mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "EEI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  shape <- parameters$variance$shape
  if(length(shape) != d)
    stop("shape incompatible with mean")
  cholSigma <- diag(sqrt(parameters$variance$scale * shape))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
      ncol = d) %*% cholSigma, MARGIN = 2, STATS = mu[, k],
      FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EEI")
}

cdensE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
        mu <- drop(parameters$mean)
  G <- length(mu)
        if(any(is.na(unlist(parameters[c("mean", "variance")]))) ||
            any(is.null(parameters[c("mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "E",
                       WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
        if(is.null(sigmasq))
                stop("variance parameters are missing")
        if(length(sigmasq) > 1)
                warning("more than one sigma-squared given")
        if(sigmasq < 0)
                stop("sigma-squared is negative")
        if(!sigmasq) {
                WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "E",
                                 WARNING = WARNING, returnCode = 9))
        }
  temp <- .Fortran("es1e",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(-1),
    as.integer(n),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[8:9]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        } 
        dimnames(z) <- list(names(data),NULL)
  structure(z, logarithm = logarithm, modelName = "E", 
                     WARNING = WARNING, returnCode = ret) 
}

emE <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepE(data, parameters = parameters, warn = warn)$z
  meE(data, z = z, prior = prior, control = control, 
      Vinv = parameters$Vinv, warn = warn)
}

estepE <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)

  mu <- drop(parameters$mean)
  G <- length(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
    return(structure(list(modelName = "E", n=n, d=1, G=G, z=z,
            parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
  }
        sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(length(sigmasq) > 1)
    warning("more than one sigma-squared specified")
  if(sigmasq < 0)
    stop("sigma-squared is negative")
  if(!sigmasq) {
    WARNING <- "sigma-squared vanishes"
    if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
    return(structure(list(modelName = "E", n=n, d=1, G=G, z=z,
            parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("es1e",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(pro),
    as.integer(n),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(names(data),NULL) 
  structure(list(modelName = "E", n = n, d = 1, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

cdensEEV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EEV", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape) ||
                   is.null(parameters$variance$orientation)) 
          stop("variance parameters are missing")
  temp <- .Fortran("eseev",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(aperm(parameters$variance$orientation,c(2,1,3))),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(p),
    double(p),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- NA
                ret <- -1
  }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EEV",
                  WARNING = WARNING, returnCode = ret)
}

emEEV <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepEEV(data, parameters = parameters, warn = warn)$z  
  meEEV(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepEEV <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EEV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape) ||
                   is.null(parameters$variance$orientation)) 
          stop("variance parameters are missing")
  temp <- .Fortran("eseev",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(aperm(parameters$variance$orientation,c(2,1,3))),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(p),
    double(p),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1 
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EEV", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

meEEV <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
              WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEV", d = p, G = G, 
              scale = NA, shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
               parameters <- list(Vinv= Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))

  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p))
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeev",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      as.integer(lwork),
      double(p * G),
      double(1),
      double(p),
      double(p * p * G),
      double(K),
      double(lwork),
      double(p),
                        PACKAGE = "mclust")[7:16]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EEV"), 
                        prior[names(prior) !="functionName"]))
    temp <- .Fortran("meeevp",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      as.integer(lwork),
      double(p * G),
      double(1),
      double(p),
      double(p * p * G),
      double(K),
      double(lwork),
      double(p),
                        PACKAGE = "mclust")[11:20]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  lapackSVDinfo <- temp[[5]]
  mu <- matrix(temp[[6]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[7]]
  shape <- temp[[8]]
  O <- aperm(array(temp[[9]], c(p, p, G)),c(2,1,3))
  pro <- temp[[10]]
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
    }
    z[] <- O[] <- shape[] <- NA
    scale <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "a z column sum fell below threshold"
      if(warn)
        warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn)
        warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- scale * shapeO(shape, O, transpose = FALSE)
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]], 
                            NULL)
## Sigma = scale * O %*% diag(shape) %*% t(O)
  variance <- list(modelName = "EEV", d = p, G = G, sigma = sigma,
                         scale = scale, shape = shape, orientation = O) 
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance) 
  structure(list(modelName = "EEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik),
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepEEV <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEV", d = p, G = G, 
                scale = NA, shape = rep(NA,p), orientation=array(NA,c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))

  }
  #  shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop(
      "improper specification of z")
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), G)
  if(is.null(prior)) {
    temp <- .Fortran("mseev",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      double(lwork),
      as.integer(lwork),
      double(p * G),
      double(1),
      double(p),
      double(p * p * G),
      double(G),
                        PACKAGE = "mclust")[7:12]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EEV"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("mseevp",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      double(lwork),
      as.integer(lwork),
      double(p * G),
      double(1),
      double(p),
      double(p * p * G),
      double(G),
                        PACKAGE = "mclust")[11:16]
  }
  lapackSVDinfo <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[3]]
  shape <- temp[[4]]
  O <- aperm( array(temp[[5]], c(p, p, G)), c(2,1,3))
  pro <- temp[[6]]
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
      warning(WARNING)
                        ret <- -4
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
      warning(WARNING)
                        ret <- -5
    }
    O[] <- shape[] <- scale <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else if(any(c(abs(scale), shape) > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn)
      warning(WARNING)
    mu[] <- pro[] <- scale <- O[] <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
                ret <- -1
  }
  else {
    sigma <- scale * shapeO(shape, O, transpose = FALSE)
                ret <- 0
  }
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(O) <- 
              list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
        variance <- list(modelName = "EEV", d = p, G = G, sigma = sigma,
                       scale = scale, shape = shape, orientation = O)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

simEEV <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
    parameters[c("mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "EEV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  shape <- parameters$variance$shape
  if(length(shape) != d)
    stop("shape incompatible with mean")
  sss <- sqrt(parameters$variance$scale * shape)
  for(k in 1:G) {
    m <- ctabel[k]
    cholSigma <- t(parameters$variance$orientation[,  , k]) * sss
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
      ncol = d) %*% cholSigma, MARGIN = 2, STATS = mu[, k],
      FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EEV")
}

hcE <- function(data, partition, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #====================================================================
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
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hc1e",
    data,
    as.integer(n),
    as.integer(partition),
    as.integer(l),
    as.integer(m),
    as.integer(ld),
    double(ld),
                PACKAGE = "mclust")[c(1, 3, 7)]
  temp[[1]] <- temp[[1]][1:m]
  temp[[2]] <- temp[[2]][1:m]
  temp[[3]] <- temp[[3]][1:m]
        change <- temp[[3]]
  structure(rbind(temp[[1]], temp[[2]]),   initialPartition = partition, 
                  dimensions = n, modelName = "E",
      call = match.call())
}

emEII <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepEII(data, parameters = parameters, warn = warn)$z
  meEII(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepEII <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  p <- ncol(data)
  n <- nrow(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
        mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
                Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EII", n=n, d=p, G=G, z=z,
                            parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
  }
        sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    warning("variance parameters are missing")
  if(sigmasq < 0)
    stop("sigma-squared is negative")
  if(!sigmasq) {
    WARNING <- "sigma-squared vanishes"
    if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EII", n=n, d=p, G=G, z=z,
                           parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("eseii",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared falls below threshold"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EII", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

hcEII <- function(data, partition, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #====================================================================
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
  #=============================================================
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
                PACKAGE = "mclust")[c(1, 9)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
        change <- temp[[2]]
  structure(t(temp[[1]]), initialPartition = partition, 
                  dimensions = dimdat, modelName = "EII", 
                  call =  match.call())
}

meEII <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
         WARNING <- "z is missing"
         if (warn) warning(WARNING)
               variance <- list(modelName = "EII", d = p, G = G, sigmasq = NA)
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeii",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(1),
      double(K),
                        PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
    data, G = G, modelName = "EII"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("meeiip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(1),
      double(K),
                        PACKAGE = "mclust")[c(11:17, 10)]
  }
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  Sigma <- diag(rep(sigmasq, p))
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
           sigmasq <= max(control$eps,0)) {
    WARNING <- "sigma-squared falls below threshold"
    if (warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if (warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if (warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL) 
        dimnames(mu) <- list(dimnames(data)[[2]], NULL) 
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL) 
        variance <- list(modelName = "EII", d = p, G = G, sigma = sigma, 
                         Sigma = Sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance = variance)
  structure(list(modelName = "EII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik),
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepEII <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EII", d = p, G = G, sigmasq = NA)
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mseii",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      double(p * G),
      double(1),
      double(G),
                        PACKAGE = "mclust")[6:8]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EII"), 
                        prior[names(prior) !="functionName"]))
    temp <- .Fortran("mseiip",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p * G),
      double(1),
      double(G),
                        PACKAGE = "mclust")[10:12]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  sigma <- array(0, c(p, p, G))
  Sigma <- diag(rep(sigmasq, p))
  for(k in 1:G)
    sigma[,  , k] <- Sigma
  WARNING <- NULL
  if(sigmasq > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]], NULL) 
        dimnames(mu) <- list(dimnames(data)[[2]], NULL) 
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL) 
        variance <- list(modelName = "EII", d = p, G = G, sigma = sigma, 
                         Sigma = Sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance = variance)
        structure(list(modelName = "EII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)

}

simEII <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EII"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  sigmasq <- parameters$variance$sigmasq
  cholSigma <- diag(rep(sqrt(sigmasq), d))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      cholSigma, MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EII")
}

meE <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be 1 dimensional")
  data <- as.vector(data)
  n <- length(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal length of data")
  K <- dimz[2]
        if (!is.null(Vinv)) {
    G <- K - 1
     if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
         WARNING <- "z is missing"
         if (warn) warning(WARNING)
               variance <- list(modelName = "E", d = 1, G = G, sigmasq = NA)
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="E", prior=prior, n=n, d=1, G=G,
                      z=z, parameters=parameters, control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("me1e",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(G),
      double(1),
      double(K),
                        PACKAGE = "mclust")[6:12]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "E"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("me1ep",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(G),
      double(1),
      double(K),
                        PACKAGE = "mclust")[c(10:16, 9)]
  }
  mu <- temp[[5]]
  names(mu) <- as.character(1:G)
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  ## log post <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
                            sigmasq <= max(control$eps,0)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- logprior <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else if(its >= control$itmax[1]) {
    WARNING <- "iteration limit reached"
    warning(WARNING)
    its <-  - its
    ret <- 1
  }
  else ret <- 0
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(names(data), NULL) 
        variance <- list(modelName = "E", d = 1, G = G, sigmasq = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "E", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)

}

mstepE <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
# number of groups 
        G <- dimz[2]
  ##
  if(all(is.na(z))) {
               WARNING <- "z is missing"
         if (warn) warning(WARNING)
               variance <- list(modelName="E", d=1, G=G, sigmasq=NA)
               parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="E", prior=prior, n=n, d=1, G=G,
                                     z = z, parameters=parameters), 
                                WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("ms1e",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(G),
      double(G),
      double(1),
      double(G),
                        PACKAGE = "mclust")[5:7]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "E"), prior[names(prior) !=
      "functionName"]))
    storage.mode(z) <- "double"
    temp <- .Fortran("ms1ep",
      as.double(data),
      z,
      as.integer(n),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(G),
      double(1),
      double(G),
                        PACKAGE = "mclust")[9:11]
  }
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  WARNING <- NULL
  if(sigmasq > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
                pro[] <- mu[] <- sigmasq <- NA
                ret <- -1
 
  }
        else ret <- 0
        dimnames(z) <- list(names(data), NULL) 
        variance <- list(modelName = "E", d = 1, G = G, sigmasq = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "E", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

simE <- function(parameters, n, seed = NULL, ...)
{
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, 2), modelName = "E"))
  }
  if(!is.null(seed))
    set.seed(seed)
  mu <- parameters$mean
  G <- length(mu)
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- rep(0, n)
  sd <- sqrt(parameters$variance$sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd)
  }
  structure(cbind(group = clabels, "1" = x), modelName = "E")
}

cdensEVI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        mu <- parameters$mean
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EVI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
  temp <- .Fortran("esevi",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- NA
                ret <- -1
  }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
         }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EVI",
                  WARNING = WARNING, returnCode = ret)
}

emEVI <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepEVI(data, parameters = parameters, warn = warn)$z  
  meEVI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepEVI <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
                Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  } 
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EVI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
             is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
  temp <- .Fortran("esevi",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EVI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

meEVI <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EVI", d = p, G = G, 
                                 scale = NA, shape = matrix(NA,p,G)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meevi",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(1),
      double(p * G),
      double(K),
                        PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EVI"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("meevip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(1),
      double(p * G),
      double(K),
                        PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      if(warn) warning("z column sum fell below threshold")
      WARNING <- "z column sum fell below threshold"
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "EVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepEVI <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EVI", d = p, G = G, 
                                 scale = NA, shape = matrix(NA,p,G)) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("msevi",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      double(p * G),
      double(1),
      double(p * G),
      double(G),
                        PACKAGE = "mclust")[6:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "EVI"), prior[names(
      prior) != "functionName"]))
    temp <- .Fortran("msevip",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p * G),
      double(1),
      double(p * G),
      double(G),
                        PACKAGE = "mclust")[10:13]
  }
  mu <- matrix(temp[[1]], p, G)
  scale <- temp[[2]]
  shape <- matrix(temp[[3]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  WARNING <- NULL
  if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!c(
    scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn)
      warning(WARNING)
    mu[] <- pro[] <- scale <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
                ret <- -1 
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
                ret <- 0
  }
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

simEVI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "EVI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  shape <- as.matrix(parameters$variance$shape)
  if(!all(dim(shape) == dim(mean)))
    stop("shape incompatible with mean")
  sss <- sqrt(parameters$variance$scale * shape)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(sss[, k]), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EVI")
}

clPairs <- function (data, classification, symbols=NULL, colors=NULL, 
          labels = dimnames(data)[[2]], CEX = 1, ...) 
{
    data <- as.matrix(data)
    m <- nrow(data)
    n <- ncol(data)
    if (missing(classification)) 
        classification <- rep(1, m)
    if (!is.factor(classification)) 
        classification <- as.factor(classification)
    l <- length(levels(classification))
    if (missing(symbols)) {
        if (l == 1) {
            symbols <- "."
        }
        if (l <= length(.mclust$classPlotSymbols)) {
            symbols <- .mclust$classPlotSymbols
        }
        else if (l <= 9) {
            symbols <- as.character(1:9)
        }
        else if (l <= 26) {
            symbols <- LETTERS[1:l]
        }
        else symbols <- rep( 16,l)
    }
    if (length(symbols) == 1) symbols <- rep(symbols, l)
    if (length(symbols) < l) {
        symbols <- rep( 16, l)
        warning("more symbols needed")
    }
    if (is.null(colors)) {
        if (l <= length(.mclust$classPlotColors)) 
          colors <- .mclust$classPlotColors[1:l]
    }
    if (length(colors) == 1) colors <- rep(colors, l)
    if (length(colors) < l) {
       colors <- rep( "black", l)
       warning("more colors needed")
    }
    pairs(x = data, labels = labels, pch = symbols[classification], 
        cex = CEX, col = colors[classification], ...)
    invisible()
}

coordProj <- function(data, dimens = c(1,2), parameters = NULL, 
   z = NULL, classification = NULL, truth = NULL, uncertainty = NULL, 
   what = c("classification", "errors", "uncertainty"), 
   quantiles = c(0.75, 0.94999999999999996), symbols = NULL, 
   colors = NULL, scale = FALSE, xlim = NULL, ylim = NULL, 
   CEX = 1, PCH = ".", identify = FALSE, ...)
{
  if(is.null(dimens)) dimens <- c(1, 2)
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mean)
    sigma <- parameters$variance$sigma
    haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any(
      is.na(sigma))
  }
  else haveParams <- FALSE
  data <- data[, dimens, drop = FALSE]
  if(dim(data)[2] != 2)
    stop("need two dimensions")
  if(is.null(xlim))
    xlim <- range(data[, 1])
  if(is.null(ylim))
    ylim <- range(data[, 2])
  if(scale) {
    par(pty = "s")
    d <- diff(xlim) - diff(ylim)
    if(d > 0) {
      ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
    }
    else {
      xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
    }
  }
  if(is.null(dnames <- dimnames(data)[[2]]))
    xlab <- ylab <- ""
  else {
    xlab <- dnames[1]
    ylab <- dnames[2]
  }
  if(haveParams) {
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      haveParams <- FALSE
      warning("covariance must be a 3D matrix")
    }
    if(G != dimpar[3]) {
      haveParams <- FALSE
      warning("means and variance parameters are incompatible")
    }
    mu <- array(mu[dimens,  ], c(2, G))
    sigma <- array(sigma[dimens, dimens,  ], c(2, 2, G))
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    noise <- classification[1] == "0"
    if(is.null(symbols)) {
      if(L <= length(.mclust$classPlotSymbols)) {
        symbols <- .mclust$classPlotSymbols
        if(noise) {
          first <- symbols[1]
          symbols[symbols == 16] <- first
          symbols[1] <- 16
        }
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    else if(length(symbols) == 1)
      symbols <- rep(symbols, L)
    if(is.null(colors)) {
      if(L <= length(.mclust$classPlotColors)) {
        colors <- .mclust$classPlotColors[1:L]
        if(noise) {
          first <- colors[1]
          colors[colors == "black"] <- first
          colors[1] <- "black"
        }
      }
    }
    else if(length(colors) == 1)
      colors <- rep(colors, L)
    if(length(symbols) < L) {
      warning("more symbols needed to show classification ")
      symbols <- rep(16,L)
    }
    if(length(colors) < L) {
      warning("more colors needed to show classification ")
      colors <- rep("black",L)
    }
  }
  if(length(what) > 1)
    what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
    badClass <- (what == "errors" && (length(unique(classification)) != length(
      unique(truth))))
    if(badClass && !bad)
      warning("classification and truth differ in number of groups")
    bad <- bad && badClass
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad)
    what <- "bad"
  switch(EXPR = what,
    classification = {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), 
          "Coordinate Projection showing Classification")
        title(main = TITLE)
      }
      for(k in 1:L) {
        I <- classification == U[k]
        points(data[I, 1], data[I, 2], pch = symbols[k], col = colors[k], cex
           = if(U[k] == "0") CEX/4 else CEX)
      }
    }
    ,
    errors = {
      ERRORS <- classError(classification, truth)$misclassified
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), 
          "Coordinate Projection showing Errors")
        title(main = TITLE)
      }
      CLASSES <- unique(as.character(truth))
      symOpen <- c(2, 0, 1, 5)
      symFill <- c(17, 15, 16, 18)
      good <- rep(TRUE, length(classification))
      good[ERRORS] <- FALSE
      if(L > 4) {
        points(data[good, 1], data[good, 2], pch = 1, col = colors, cex = CEX)
        points(data[!good, 1], data[!good, 2], pch = 16, cex = CEX)
      }
      else {
        for(k in 1:L) {
          K <- truth == CLASSES[k]
          if(any(I <- (K & good))) {
            points(data[I, 1], data[I, 2], pch = symOpen[k], col = colors[k], 
              cex = CEX)
          }
          if(any(I <- (K & !good))) {
            points(data[I, 1], data[I, 2], pch = symFill[k], cex = CEX)
          }
        }
      }
    }
    ,
    uncertainty = {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), 
          "Coordinate Projection showing Uncertainty")
        title(main = TITLE)
      }
      breaks <- quantile(uncertainty, probs = sort(quantiles))
      I <- uncertainty <= breaks[1]
      points(data[I, 1], data[I, 2], pch = 16, col = "gray75", cex = 0.5 * CEX)
      I <- uncertainty <= breaks[2] & !I
      points(data[I, 1], data[I, 2], pch = 16, col = "gray50", cex = 1 * CEX)
      I <- uncertainty > breaks[2] & !I
      points(data[I, 1], data[I, 2], pch = 16, col = "black", cex = 1.5 * CEX)
    }
    ,
    {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), "Coordinate Projection")
        title(main = TITLE)
      }
      points(data[, 1], data[, 2], pch = PCH, cex = CEX)
    }
    )
  if(haveParams) {
    ## plot ellipsoids
    for(k in 1:G)
      mvn2plot(mu = mu[, k], sigma = sigma[,  , k], k = 15)
  }
  invisible()
}

imputePairs <- function (x, impx, symbols = c(16,1), colors = c("black", "red"),
          labels, panel = points, ...,  
    lower.panel = panel, 
    upper.panel = panel, diag.panel = NULL, text.panel = textPanel, 
    label.pos = 0.5 + has.diag/3, cex.labels = NULL, font.labels = 1, 
    row1attop = TRUE, gap = 1) 
{
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x, 
        y, txt, cex = cex, font = font)
    localAxis <- function(side, x, y, xpd, bg, col = NULL, main, 
        oma, ...) {
        if (side%%2 == 1) 
            Axis(x, side = side, xpd = NA, ...)
        else Axis(y, side = side, xpd = NA, ...)
    }
    localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
    localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
    localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
    localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)
    dots <- list(...)
    nmdots <- names(dots)
    if (!is.matrix(x)) {
        x <- as.data.frame(x)
        for (i in seq_along(names(x))) {
            if (is.factor(x[[i]]) || is.logical(x[[i]])) 
                x[[i]] <- as.numeric(x[[i]])
            if (!is.numeric(unclass(x[[i]]))) 
                stop("non-numeric argument to 'pairs'")
        }
    }
    else if (!is.numeric(x)) 
        stop("non-numeric argument to 'pairs'")
    panel <- match.fun(panel)
    if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)) 
        lower.panel <- match.fun(lower.panel)
    if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)) 
        upper.panel <- match.fun(upper.panel)
    if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)) 
        diag.panel <- match.fun(diag.panel)
    if (row1attop) {
        tmp <- lower.panel
        lower.panel <- upper.panel
        upper.panel <- tmp
        tmp <- has.lower
        has.lower <- has.upper
        has.upper <- tmp
    }
    nc <- ncol(x)
    if (nc < 2) 
        stop("only one column in the argument to 'pairs'")
    has.labs <- TRUE
    if (missing(labels)) {
        labels <- colnames(x)
        if (is.null(labels)) 
            labels <- paste("var", 1:nc)
    }
    else if (is.null(labels)) 
        has.labs <- FALSE
    oma <- if ("oma" %in% nmdots) 
        dots$oma
    else NULL
    main <- if ("main" %in% nmdots) 
        dots$main
    else NULL
    if (is.null(oma)) {
        oma <- c(4, 4, 4, 4)
        if (!is.null(main)) 
            oma[3] <- 6
    }
    opar <- par(mfrow = c(nc, nc), mar = rep.int(gap/2, 4), oma = oma)
    on.exit(par(opar))
    for (i in if (row1attop) 
        1:nc
    else nc:1) for (j in 1:nc) {
        localPlot(impx[, j], impx[, i], xlab = "", ylab = "", axes = FALSE, 
            type = "n", ...)
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
            box()
            if (i == 1 && (!(j%%2) || !has.upper || !has.lower)) 
                localAxis(1 + 2 * row1attop, impx[, j], impx[, i], 
                  ...)
            if (i == nc && (j%%2 || !has.upper || !has.lower)) 
                localAxis(3 - 2 * row1attop, impx[, j], impx[, i], 
                  ...)
            if (j == 1 && (!(i%%2) || !has.upper || !has.lower)) 
                localAxis(2, impx[, j], impx[, i], ...)
            if (j == nc && (i%%2 || !has.upper || !has.lower)) 
                localAxis(4, impx[, j], impx[, i], ...)
            mfg <- par("mfg")
            if (i == j) {
                if (has.diag) 
                  localDiagPanel(as.vector(impx[, i]), ...)
                if (has.labs) {
                  par(usr = c(0, 1, 0, 1))
                  if (is.null(cex.labels)) {
                    l.wid <- strwidth(labels, "user")
                    cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                  }
                  text.panel(0.5, label.pos, labels[i], cex = cex.labels, 
                    font = font.labels)
                }
            }
            else if (i < j) { 
                classification <- as.numeric(apply(x[,c(i,j)], 1, 
                                             function(x) any(is.na(x)))) + 1
               localLowerPanel(as.vector(impx[, j]), as.vector(impx[, 
                  i]), pch = symbols[classification], 
                       col = colors[classification], ...)
              }
            else {
                classification <- as.numeric(apply(x[,c(i,j)], 1, 
                                             function(x) any(is.na(x)))) + 1
               localUpperPanel(as.vector(impx[, j]), as.vector(impx[, 
                i]), pch = symbols[classification], 
                     col = colors[classification], ...)
            }
            if (any(par("mfg") != mfg)) 
                stop("the 'panel' function made a new plot")
        }
        else par(new = FALSE)
    }
    if (!is.null(main)) {
        font.main <- if ("font.main" %in% nmdots) 
            dots$font.main
        else par("font.main")
        cex.main <- if ("cex.main" %in% nmdots) 
            dots$cex.main
        else par("cex.main")
        mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
    }
    invisible(NULL)
}

mclust1Dplot <- function(data, parameters = NULL, z = NULL,
         classification = NULL, truth = NULL, uncertainty = NULL, 
         what = c("classification", "density", "errors", "uncertainty"), 
         symbols = NULL, colors = NULL, ngrid = length(data),  
         xlab = NULL,  xlim = NULL, CEX  = 1, 
         identify = FALSE, ...) 
{

  grid1 <- function (n, range = c(0, 1), edge = TRUE) 
  {
    if (any(n < 0 | round(n) != n)) 
        stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if(edge) 
      { G <- seq(from = min(range), to = max(range), 
                 by = abs(diff(range))/(n - 1)) }
    else 
      { lj <- abs(diff(range))
        incr <- lj/(2 * n)
        G <- seq(from = min(range) + incr, to = max(range) - incr, 
                 by = 2 * incr) }
    G
  }
  
  densNuncer <- function(data, parameters) 
  {
    cden <- cdensV(data = data, parameters = parameters)
    if(parameters$variance$G != 1) 
      { z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = parameters$pro)
        den <- apply(z, 1, sum)
        z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
        data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
      }
    else 
      { data.frame(density = cden, uncertainty =  rep(NA, length(cden))) }
  }

  if (is.null(xlab)) xlab <- " "
  p <- ncol(as.matrix(data))
  if (p != 1) 
      stop("for one-dimensional data only")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(classification) && !is.null(z))
              classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
              uncertainty <- 1 - apply(z, 1, max)
  if (!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mu)
    sigmasq <- parameters$variance$sigmasq
    haveParams <- !is.null(mu) && !is.null(sigmasq) && 
                  !any(is.na(mu)) && !any(is.na(sigmasq)) 
  }
  else haveParams <- FALSE
  if (is.null(xlim)) xlim <- range(data)
  if (haveParams) {
    G <- length(mu)
    if ((l <- length(sigmasq)) == 1) {
          sigmasq <- rep(sigmasq, G)
    }
    else if (l != G) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
  }
  if (!is.null(truth)) {
    if (is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if (length(unique(truth)) != 
          length(unique(classification))) truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) 
  {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    if(is.null(symbols)) 
      { symbols <- rep("|", L) }
    else if(length(symbols) == 1) 
      { symbols <- rep(symbols, L) }
    else if(length(symbols) < L) 
      { warning("more symbols needed to show classification")
        symbols <- rep("|", L) }
    if(is.null(colors))
      { colors <- .mclust$classPlotColors[1:L] }
    else if(length(colors) == 1) 
      { colors <- rep(colors, L) }
    else if(length(colors) < L)
      { warning("more colors needed to show classification")
        colors <- rep("black", L) }
  }
  if (length(what) > 1) what <- what[1]
  choices <- c("classification", "density", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if (m) { 
    type <- choices[m] 
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && 
              (is.null(classification) || is.null(truth)))
    if (bad) warning("insufficient input for specified plot")
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if (bad) what <- "bad"
  M <- L
  switch(EXPR = what,
  "classification" = 
  { plot(data, seq(from = 0, to = M, length = n), type = "n", 
         xlab = xlab, ylab = "", xlim = xlim, yaxt = "n", main = "", ...)
    axis(side = 2, at = 0:M, labels = c("", sort(unique(classification))))
    if(identify) title("Classification")
    for(k in 1:L) 
       { I <- classification == U[k]
         points(data[I], rep(0, length(data[I])), 
                pch = symbols[k], cex = CEX)
         points(data[I], rep(k, length(data[I])), 
                pch = symbols[k], col = colors[k], cex = CEX)
       }
  },
  "errors" = 
  { ERRORS <- classError(classification, truth)$misclassified
    plot(data, seq(from = 0, to = M, length = n), type = "n", 
         xlab = xlab, ylab = "", xlim = xlim, yaxt = "n", main = "", ...)
    axis(side = 2, at = 0:M, labels = c("", unique(classification)))
    if(identify) title("Classification Errors")
    good <- rep(TRUE, length(classification))
    good[ERRORS] <- FALSE
    sym <- "|"
    for(k in 1:L) 
       { K <- classification == U[k]
         I <- K & good
         if(any(I)) 
           { if(FALSE) 
               { sym <- if (L > 4) 
                        1
                        else if (k == 4) 
                        5
                        else k - 1
               }
             l <- sum(as.numeric(I))
             points(data[I], rep(0, l), pch = sym, 
                    col = colors[k], cex = CEX)
           }
        I <- K & !good
        if(any(I)) 
          { if(FALSE) 
              { sym <- if (L > 5) 
                          16
                       else k + 14 }
            l <- sum(as.numeric(I))
            points(data[I], rep(k, l), pch = sym, 
                   col = colors[k], cex = CEX)
            # points(data[I], rep(0, l), pch = sym, cex = CEX)
            # points(data[I], rep(-0.5, l), pch = sym, cex = CEX)
        }
    }
      },
      "uncertainty" = {
          x <- grid1(n = ngrid, range = xlim, edge = TRUE)
          lx <- length(x)
          Z <- densNuncer(data = x, parameters = parameters)
          plot(x, Z$uncertainty, xlab = xlab, ylab = "uncertainty", 
              xlim = xlim, ylim = c(0,1), type = "l", main = "", ...)
          if (identify)  title("Uncertainty")
      },
      "density" = {
          if (is.null(parameters$pro) && parameters$variance$G != 1) 
              stop("mixing proportions missing")
          x <- grid1(n = ngrid, range = xlim, edge = TRUE)
          lx <- length(x)
          Z <- densNuncer(data = x, parameters = parameters)
          plot(x, Z$density, xlab = xlab, ylab = "density", xlim = xlim, 
              type = "l", main = "", ...)
          if (identify) title("Density")
      },
      {
          plot(data, rep(0, n), type = "n", xlab = "", ylab = "", 
              xlim = xlim, main = "", ...)
          points(data, rep(0, n), pch = "|", cex = CEX)
          if (identify) title("Point Plot")
          return(invisible())
      }
     )
  invisible()
}

mclust2Dplot <- function (data, parameters = NULL, z = NULL,
          classification = NULL, truth = NULL, uncertainty = NULL, 
          what = c("classification", "uncertainty", "errors"), 
          quantiles = c(0.75, 0.95), symbols = NULL, colors = NULL, 
          scale = FALSE, xlim = NULL, ylim = NULL, CEX = 1, PCH = ".", 
          identify = FALSE, swapAxes = FALSE, ...) 
{
        if(dim(data)[2] != 2)
                stop("data must be two dimensional")
        if(is.null(classification) && !is.null(z))
                classification <- map(z)
        if(is.null(uncertainty) && !is.null(z))
                uncertainty <- 1 - apply(z, 1, max)
        if (!is.null(parameters)) {
          mu <- parameters$mean
          L <- ncol(mu)
          sigma <- parameters$variance$sigma
          haveParams <- !is.null(mu) && !is.null(sigma) && 
                         !any(is.na(mu)) && !any(is.na(sigma)) 
        }
        else haveParams <- FALSE
        if (is.null(xlim))
                xlim <- range(data[, 1])
        if (is.null(ylim))
                ylim <- range(data[, 2])
        if (scale) {
                par(pty = "s")
                d <- diff(xlim) - diff(ylim)
                if(d > 0) {
                        ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
                }
                else {
                        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
                }
        }
        if(is.null(dnames <- dimnames(data)[[2]]))
                xlab <- ylab <- ""
        else {
                xlab <- dnames[1]
                ylab <- dnames[2]
        }
       if (haveParams) { 
          G <- ncol(mu)
          dimpar <- dim(sigma)
          if (length(dimpar) != 3) {
            haveParams <- FALSE
            warning("covariance must be a 3D matrix")
          }
          if (G != dimpar[3]) {
            haveParams <- FALSE
            warning("means and variance parameters are incompatible")
          }
          mu <- array(mu, c(2, G))
          sigma <- array(sigma, c(2, 2, G))
        }
      if (swapAxes) {
        if (haveParams) {
            mu <- mu[2:1, ]
            sigma <- sigma[2:1, 2:1, ]
        }
        data <- data[, 2:1]
      }
      if (!is.null(truth)) {
                if (is.null(classification)) {
                    classification <- truth
                    truth <- NULL
                }
                else {
                    if (length(unique(truth)) != 
                        length(unique(classification))) truth <- NULL
                    else truth <- as.character(truth)
                }
        }
        if (charmatch("classification", what, nomatch = 0) && 
            is.null(classification) && !is.null(z))
            classification <- map(z)  
        if (!is.null(classification)) {
                classification <- as.character(classification)
                U <- sort(unique(classification))
                L <- length(U)
                noise <- classification[1] == "0"
                if(is.null(symbols)) {
                        if(L <= length(.mclust$classPlotSymbols)) {
                                symbols <- .mclust$classPlotSymbols
                                if (noise) {
                                   first <- symbols[1]
                                   symbols[symbols == 16] <- first
                                   symbols[1] <- 16
                                }
                        }
                        else if(L <= 9) {
                                symbols <- as.character(1:9)
                        }
                        else if(L <= 26) {
                                symbols <- LETTERS
                        }
                }
                if(is.null(colors)) {
                        if(L <= length(.mclust$classPlotColors)) {
                                colors <- .mclust$classPlotColors[1:L]
                                if (noise) {
                                   first <- colors[1]
                                   colors[colors == "black"] <- first
                                   colors[1] <- "black"
                                }
                        }
                }
                else if (length(colors) == 1) colors <- rep(colors, L)
            if(length(symbols) < L) {
                warning("more symbols needed to show classification ")
                symbols <- rep(16,L)
            }
            if(length(colors) < L) {
                warning("more colors needed to show classification ")
                colors <- rep("black",L)
            }
        }
        if (length(what) > 1) what <- what[1]
        choices <- c("classification", "errors", "uncertainty")
        m <- charmatch(what, choices, nomatch = 0)
        if (m) { 
          what <- choices[m] 
          bad <- what == "classification" && is.null(classification)
          bad <- bad || (what == "uncertainty" && is.null(uncertainty))
          bad <- bad || (what == "errors" && 
                (is.null(classification) || is.null(truth)))
          if (bad) warning("insufficient input for specified plot")
        }
        else {
          bad <- !m
          warning("what improperly specified")
        }
        if (bad) what <- "bad"
        switch(EXPR = what,
            "classification"= {
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Classification")
            for (k in 1:L) {
                I <- classification == U[k]
                points(data[I, 1], data[I, 2], pch = symbols[k], 
                       col = colors[k], cex = if (U[k] == "0") CEX/4 else CEX)
            }
            },
            "errors" = {
            ERRORS <- classError(classification, truth)$misclassified
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Classification Errors")
            CLASSES <- unique(as.character(truth))
            symOpen <- c(2, 0, 1, 5)
            symFill <- c(17, 15, 16, 18)
            good <- rep(TRUE,length(classification))
            good[ERRORS] <- FALSE
            if (L > 4) {
                points(data[good, 1], data[good, 2], pch = 1, 
                  col = colors, cex = CEX)
                points(data[!good, 1], data[!good, 2], pch = 16, 
                  cex = CEX)
            }
            else {
                for (k in 1:L) {
                  K <- truth == CLASSES[k]
                  points(data[K, 1], data[K, 2], pch = symOpen[k], 
                    col = colors[k], cex = CEX)
                  if (any(I <- (K & !good))) {
                    points(data[I, 1], data[I, 2], pch = symFill[k], 
                      cex = CEX)
                  }
                }
            }
        },
        "uncertainty" = {
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Classification Uncertainty")
            breaks <- quantile(uncertainty, probs = sort(quantiles))
            I <- uncertainty < breaks[1]
            points(data[I, 1], data[I, 2], pch = 16, col = "gray75",
                   cex = 0.5 * CEX)
            I <- uncertainty < breaks[2] & !I
            points(data[I, 1], data[I, 2], pch = 16, col = "gray50",
                   cex = 1 * CEX)
            I <- uncertainty >= breaks[2]
            points(data[I, 1], data[I, 2], pch = 16, col = "black",
                   cex = 1.5 * CEX)
        },
        {
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Point Plot")
            points(data[, 1], data[, 2], pch = PCH, cex = CEX)
          }
      )
    if (haveParams) {
## plot ellipsoids
          for (k in 1:G) 
             mvn2plot(mu = mu[, k], sigma = sigma[ ,  , k], k = 15)
    }
    invisible()
}

mvn2plot <- function (mu, sigma, k = 15, alone = FALSE, col = 1) 
{
    p <- length(mu)
    if (p != 2) 
        stop("two-dimensional case only")
    if (any(unique(dim(sigma)) != p)) 
        stop("mu and sigma are incompatible")
    ev <- eigen(sigma, symmetric = TRUE)
    s <- sqrt(rev(sort(ev$values)))
    V <- t(ev$vectors[, rev(order(ev$values))])
    theta <- (0:k) * (pi/(2 * k))
    x <- s[1] * cos(theta)
    y <- s[2] * sin(theta)
    xy <- cbind(c(x, -x, -x, x), c(y, y, -y, -y))
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
    if (alone) {
        xymin <- apply(xy, 2, FUN = "min")
        xymax <- apply(xy, 2, FUN = "max")
        r <- ceiling(max(xymax - xymin)/2)
        xymid <- (xymin + xymax)/2
        plot(xy[, 1], xy[, 2], xlim = c(-r, r) + xymid[1], ylim = c(-r, 
            r) + xymid[2], xlab = "x", ylab = "y", type = "n")
    }
    l <- length(x)
    i <- 1:l
    for (k in 1:4) {
        lines(xy[i, ], col = col)
        i <- i + l
    }
    x <- s[1]
    y <- s[2]
    xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
    lines(xy[1:2, ], lty = 2, col = col)
    lines(xy[3:4, ], lty = 2, col = col)
    points(mu[1], mu[2], pch = 8, col = col)
    invisible()
}

plot.Mclust <- function(x, 
          what = c("BIC", "classification", "uncertainty", "density"), 
          dimens = NULL, xlab = NULL, ylim = NULL,  
          addEllipses = TRUE, identify = TRUE,
          legendArgs = list(x="bottomright", ncol=2, cex=1),
           ...) 
{
  
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "Mclust")) 
    stop("object not of class \"Mclust\"")

  data <- eval.parent(object$call$data)
  data <- as.matrix(data)
  p <- ncol(data)
  if(p == 1) 
    colnames(data) <- deparse(x$call$data)
  if(is.null(dimens)) 
    dimens <- seq(p)
  else
    dimens <- dimens[dimens <= p]
  d <- length(dimens)

  what <- match.arg(what, c("BIC", "classification", 
                            "uncertainty", "density"), 
                    several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  if(length(what) > 1)
    { par(ask = TRUE)
      on.exit(par(oldpar))
    }

  if(any(match("BIC", what, nomatch = FALSE)))
    { plot.mclustBIC(object$BIC, xlab=xlab, ylim=ylim,
                     legendArgs=legendArgs,...)
      if(length(what) == 1)
      return(invisible())
    }

  if(p == 1)
    {
      if(any(match("classification", what, nomatch = FALSE)))
         mclust1Dplot(data = data, 
                      # parameters = object$parameters, 
                      what = "classification",
                      classification = object$classification,
                      z = object$z, 
                      xlab = colnames(data)[dimens], 
                      identify = identify, ...)
      if(any(match("uncertainty", what, nomatch = FALSE)))
         mclust1Dplot(data = data,
                      parameters = object$parameters,
                      z = object$z, what = "uncertainty", 
                      xlab = colnames(data)[dimens], 
                      identify = identify, ...)
      if(any(match("density", what, nomatch = FALSE)))
         mclust1Dplot(data = data,
                      parameters = object$parameters,
                      z = object$z, what = "density", 
                      xlab = colnames(data)[dimens], 
                      identify = identify, ...)
    }
  
  if(p == 2) 
    {
      if(any(match("classification", what, nomatch = FALSE)))
        {
          if(addEllipses)
            { coordProj(data = data, what = "classification", 
                        parameters = object$parameters, z = object$z, 
                        dimens = dimens, identify = FALSE, ...)
              title("Classification") 
            }
          else
          { mclust2Dplot(data = data, what = "classification", 
                         classification = object$classification, 
                         # z = object$z, 
                         identify = FALSE, ...) 
          }
        }
      if(any(match("uncertainty", what, nomatch = FALSE)))
        mclust2Dplot(data = data, parameters = object$parameters, 
                     z = object$z, what = "uncertainty", 
                     identify = identify, ...)
      if(any(match("density", what, nomatch = 0)))
         surfacePlot(data = data, parameters = object$parameters,
                     what = "density", nlevels = 11,
                     transformation = "log",
                     identify = identify, ...)
  }

  if(p > 2) 
    {
      if(any(match("classification", what, nomatch = FALSE)))
        { 
          if(d == 2)
            { if(addEllipses)
                { coordProj(data = data, what = "classification", 
                            parameters = object$parameters, z = object$z, 
                            dimens = dimens, identify = identify, ...) 
                }
              else
                { mclust2Dplot(data = data[,dimens], what = "classification", 
                               classification = object$classification,
                               # z = object$z, 
                               identify = identify, ...)
                }
            }
          else
            { if(addEllipses)
                { on.exit(par(oldpar))
                  par(mfrow = c(d, d), 
                      mar = rep(c(0.3,0.3/2),each=2), 
                      oma = c(4, 4, 4, 4))
                  for(i in seq(d))
                     { for(j in seq(d)) 
                          { if(i == j) 
                              { plot(0,0,type="n",xlab="",ylab="",axes=FALSE)
                                text(0,0, colnames(data[,dimens])[i], 
                                     cex=1.5, adj=0.5)
                                box()
                              } 
                            else 
                              { coordProj(data = data, 
                                          what = "classification", 
                                          parameters = object$parameters,
                                          z = object$z,
                                          dimens = dimens[c(j,i)], 
                                          identify = FALSE, 
                                          xaxt = "n", yaxt = "n", ...)
                              }
                            if(i == 1 && (!(j%%2))) axis(3)
                            if(i == d && (j%%2))   axis(1)
                            if(j == 1 && (!(i%%2))) axis(2)
                            if(j == d && (i%%2))   axis(4)
                          }
                     }
                }
              else
                { clPairs(data[,dimens], gap = 0.3, cex.labels = 1.5,
                          classification = object$classification, ...) 
                }
           }
        }
        
      if(any(match("uncertainty", what, nomatch = FALSE)))
        { 
          if(d == 2)
            { coordProj(data = data, parameters = object$parameters, 
                        z = object$z, what = "uncertainty", 
                        dimens = dimens, identify = identify, ...) }
          else
            { on.exit(par(oldpar))
              par(mfrow = c(d, d), 
                  mar = rep(c(0.3,0.3/2),each=2), 
                  oma = c(4, 4, 4, 4))
              for(i in seq(d))
                 { for(j in seq(d)) 
                      { if(i == j) 
                          { plot(0,0,type="n",xlab="",ylab="",axes=FALSE)
                            text(0,0, colnames(data[,dimens])[i], 
                                 cex=1.5, adj=0.5)
                            box()
                          } 
                        else 
                          { coordProj(data = data, 
                                      what = "uncertainty", 
                                      parameters = object$parameters, 
                                      z = object$z,
                                      dimens = dimens[c(j,i)], 
                                      identify = FALSE, 
                                      xaxt = "n", yaxt = "n", ...)
                          }
                        if(i == 1 && (!(j%%2))) axis(3)
                        if(i == d && (j%%2))   axis(1)
                        if(j == 1 && (!(i%%2))) axis(2)
                        if(j == d && (i%%2))   axis(4)
                      }
                 }
            }
        }
      
      if(any(match("density", what, nomatch = FALSE)))
        { objdens <- object
          objdens$varname <- colnames(data)
          objdens$range <- if(objdens$d > 1) apply(data, 2, range) else range(data)
          plotDensityMclustd(objdens, nlevels = 11, ...)
        }
    }
  
  invisible()
}

plot.mclustBIC <- function(x, G = NULL, modelNames = NULL, 
                           symbols = NULL, colors = NULL, 
                           xlab = NULL, ylab = "BIC", ylim = NULL, 
                           legendArgs = list(x = "bottomright", ncol = 2, cex = 1), 
                           ...)
{

  if(is.null(xlab)) xlab <- "Number of components"
  fill <- FALSE
  subset <- !is.null(attr(x, "initialization")$subset)
  noise <- !is.null(attr(x, "initialization")$noise)
  ret <- attr(x, "returnCodes") == -3
## if(!subset && any(ret) && fill) {
##    x <- bicFill(x, ret, n, d)
##}
  n <- ncol(x)
  dnx <- dimnames(x)
  ##
  x <- matrix(as.vector(x), ncol = n)
  dimnames(x) <- dnx
  if(is.null(modelNames))
    modelNames <- dimnames(x)[[2]]
  if(is.null(G))
    G <- as.numeric(dimnames(x)[[1]])
# BIC <- x[as.character(G), modelNames, drop = FALSE]
# X <- is.na(BIC)
# nrowBIC <- nrow(BIC)
# ncolBIC <- ncol(BIC)
  if(is.null(symbols)) 
    { colNames <- dimnames(x)[[2]]
      m <- length(modelNames)
      if(is.null(colNames)) 
        { symbols <- if(m > 9) LETTERS[1:m] else as.character(1:m)
          names(symbols) <- modelNames
        }
    else 
        { symbols <- .mclust$bicPlotSymbols[modelNames] }
    }
  if(is.null(colors)) 
    { colNames <- dimnames(x)[[2]]
      if(is.null(colNames)) 
        { colors <- 1:m
          names(colors) <- modelNames
        }
      else 
        { colors <- .mclust$bicPlotColors[modelNames] }
    }
  x <- x[,modelNames, drop = FALSE]
  if(is.null(ylim))
    ylim <- range(as.vector(x[!is.na(x)]))
  matplot(as.numeric(dnx[[1]]), x, type = "b", 
          xlim = range(G), ylim = ylim,
          pch = symbols, col = colors, lty = 1,
          xlab = xlab, ylab = ylab, main = "")
 if(!is.null(legendArgs))
   { do.call("legend", c(list(legend = modelNames, col = colors, pch = symbols),
                         legendArgs)) }
 invisible(symbols)
}

randProj <- function(data, seeds = 0, parameters = NULL, z = NULL, classification = NULL, 
  truth = NULL, uncertainty = NULL, what = c("classification", "errors", 
  "uncertainty"), quantiles = c(0.75, 0.94999999999999996), symbols = NULL, 
  colors = NULL, scale = FALSE, xlim = NULL, ylim = NULL, CEX = 1, PCH = ".", 
  identify = FALSE, ...)
{
  if(scale) par(pty = "s")
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mu)
    sigma <- parameters$variance$sigma
    haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any(
      is.na(sigma))
  }
  else haveParams <- FALSE
  xlab <- ylab <- ""
  p <- ncol(data)
  if(haveParams) {
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) {
      haveParams <- FALSE
      warning("covariance must be a 3D matrix")
    }
    if(G != dimpar[3]) {
      haveParams <- FALSE
      warning("means and variance parameters are incompatible")
    }
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    noise <- classification[1] == "0"
    if(is.null(symbols)) {
      if(L <= length(.mclust$classPlotSymbols)) {
        symbols <- .mclust$classPlotSymbols
        if(noise) {
          first <- symbols[1]
          symbols[symbols == 16] <- first
          symbols[1] <- 16
        }
      }
      else if(L <= 9) {
        symbols <- as.character(1:9)
      }
      else if(L <= 26) {
        symbols <- LETTERS
      }
    }
    else if(length(symbols) == 1)
      symbols <- rep(symbols, L)
    if(is.null(colors)) {
      if(L <= length(.mclust$classPlotColors)) {
        colors <- .mclust$classPlotColors[1:L]
        if(noise) {
          first <- colors[1]
          colors[colors == "black"] <- first
          colors[1] <- "black"
        }
      }
    }
    else if(length(colors) == 1)
      colors <- rep(colors, L)
    if(length(symbols) < L) {
      warning("more symbols needed to show classification ")
      symbols <- rep(16,L)
    }
    if (length(colors) < L) {
      warning("more colors needed to show classification ")
      colors <- rep("black",L)
    }
  }
  if(length(what) > 1)
    what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad)
    what <- "bad"
  nullXlim <- is.null(xlim)
  nullYlim <- is.null(ylim)
  if(length(seeds) > 1)
    par(ask = TRUE)
  for(seed in seeds) {
    set.seed(seed)
    Q <- orth2(p)
    Data <- as.matrix(data) %*% Q
    if(dim(Data)[2] != 2)
      stop("need two dimensions")
    if(nullXlim)
      xlim <- range(Data[, 1])
    if(nullYlim)
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
    switch(EXPR = what,
      classification = {
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        for(k in 1:L) {
          I <- classification == U[k]
          points(Data[I, 1], Data[I, 2], pch = symbols[k], col = colors[k], cex
             = CEX)
        }
        if(identify) {
          TITLE <- paste("Random Projection showing Classification: seed = ", 
            seed)
          title(TITLE)
        }
      }
      ,
      errors = {
        ERRORS <- classError(classification, truth)$misclassified
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        if(identify) {
          TITLE <- paste("Random Projection showing Errors: seed = ", seed)
          title(TITLE)
        }
        CLASSES <- unique(as.character(truth))
        symOpen <- c(2, 0, 1, 5)
        symFill <- c(17, 15, 16, 18)
        good <- !ERRORS
        if(L > 4) {
          points(Data[good, 1], Data[good, 2], pch = 1, col = colors, cex = CEX)
          points(Data[!good, 1], Data[!good, 2], pch = 16, cex = CEX)
        }
        else {
          for(k in 1:L) {
            K <- truth == CLASSES[k]
            points(Data[K, 1], Data[K, 2], pch = symOpen[k], col = colors[k], 
              cex = CEX)
            if(any(I <- (K & ERRORS))) {
              points(Data[I, 1], Data[I, 2], pch = symFill[k], cex = CEX)
            }
          }
        }
      }
      ,
      uncertainty = {
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        if(identify) {
          TITLE <- paste("Random Projection showing Uncertainty: seed = ", seed
            )
          title(TITLE)
        }
        breaks <- quantile(uncertainty, probs = sort(quantiles))
        I <- uncertainty <= breaks[1]
        points(Data[I, 1], Data[I, 2], pch = 16, col = "gray75", cex = 0.5 * CEX)
        I <- uncertainty <= breaks[2] & !I
        points(Data[I, 1], Data[I, 2], pch = 16, col = "gray50", cex = 1 * CEX)
        I <- uncertainty > breaks[2] & !I
        points(Data[I, 1], Data[I, 2], pch = 16, col = "black", cex = 1.5 * CEX)
      }
      ,
      {
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        if(identify) {
          TITLE <- paste("Random Projection: seed = ", seed)
          title(TITLE)
        }
        points(Data[, 1], Data[, 2], pch = PCH, cex = CEX)
      }
      )
    if(haveParams) {
      ## plot ellipsoids
      muTrans <- crossprod(Q, mu)
      sigmaTrans <- array(apply(cho, 3, function(R, Q)
      crossprod(R %*% Q), Q = Q), c(2, 2, G))
      for(k in 1:G)
        mvn2plot(mu = muTrans[, k], sigma = sigmaTrans[,  , k], k = 15)
    }
  }
  invisible()
}

summary.mclustBIC <- function(object, data, G, modelNames, ...)
{
  mc <- match.call(expand.dots = FALSE)
  if (is.null(attr(object,"initialization")$noise)) {
    mc[[1]] <- as.name("summaryMclustBIC")
  }
  else {
    mc[[1]] <- as.name("summaryMclustBICn")
  }
  ans <- eval(mc, parent.frame())
  Glabels <- dimnames(object)[[1]]
  if (length(Glabels) != 1 && (!missing(G) && length(G) > 1)) {
   Grange <- range(as.numeric(Glabels))
   if (match(ans$G, Grange, nomatch = 0))
  warning("best model occurs at the min or max # of components considered")
  }
  ans
}

surfacePlot <- function(data, parameters, 
                        type = c("contour", "image", "persp"), 
                        what = c("density", "uncertainty"), 
                        transformation = c("none", "log", "sqrt"), 
                        grid = 50, nlevels = 11, levels = NULL, scale = FALSE, 
                        xlim = NULL, ylim = NULL, identify = FALSE, 
                        verbose = FALSE, swapAxes = FALSE, ...) 
{
  grid1 <- function(n, range = c(0, 1), edge = TRUE) {
    if (any(n < 0 | round(n) != n)) 
        stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if (edge) {
        G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n-1))
    }
    else {
        lj <- abs(diff(range))
        incr <- lj/(2 * n)
        G <- seq(from = min(range) + incr, to = max(range) - incr, by = 2 * incr)
    }
    G
  }
  grid2 <- function(x, y) {
    lx <- length(x)
    ly <- length(y)
    xy <- matrix(0, nrow = lx * ly, ncol = 2)
    l <- 0
    for (j in 1:ly) {
        for (i in 1:lx) {
            l <- l + 1
            xy[l, ] <- c(x[i], y[j])
        }
    }
    xy
  }
  data <- as.matrix(data)
  if (dim(data)[2] != 2) 
    stop("data must be two dimensional")
  densNuncer <- function(modelName, data, parameters) {
    if (is.null(parameters$variance$cholsigma)) {
        parameters$variance$cholsigma <- parameters$variance$sigma
        G <- dim(parameters$variance$sigma)[3]
        for(k in 1:G) 
            parameters$variance$cholsigma[, , k] <- chol(parameters$variance$sigma[, , k])
    }
    cden <- cdensVVV(data = data, parameters = parameters)
    z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = parameters$pro)
    den <- apply(z, 1, sum)
    z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
    data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
  }
  pro <- parameters$pro
  mu <- parameters$mean
  sigma <- parameters$variance$sigma
  haveParams <- !is.null(mu) && !is.null(sigma) && !is.null(pro) && 
    !any(is.na(mu)) && !any(is.na(sigma)) && !(any(is.na(pro)))
  if (haveParams) {
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if (length(dimpar) != 3) {
        haveParams <- FALSE
        warning("covariance must be a 3D matrix")
    }
    if (G != dimpar[3]) {
        haveParams <- FALSE
        warning("means and variance parameters are incompatible")
    }
    mu <- array(mu, c(2, G))
    sigma <- array(sigma, c(2, 2, G))
  }
  if (!haveParams) 
    stop("need parameters to compute density")
  if (swapAxes) {
    if (haveParams) {
        parameters$pro <- pro[2:1]
        parameters$mean <- mu[2:1, ]
        parameters$variance$sigma <- sigma[2:1, 2:1, ]
    }
    data <- data[, 2:1]
  }
  if (is.null(xlim)) 
    xlim <- range(data[, 1])
  if (is.null(ylim)) 
    ylim <- range(data[, 2])
  if (scale) {
    par(pty = "s")
    d <- diff(xlim) - diff(ylim)
    if (d > 0) {
        ylim <- c(ylim[1] - d/2, ylim[2] + d/2)
    }
    else {
        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
    }
  }
  if (is.null(dnames <- dimnames(data)[[2]])) 
    xlab <- ylab <- ""
  else {
    xlab <- dnames[1]
    ylab <- dnames[2]
  }
  if (length(grid) == 1) 
    grid <- c(grid, grid)
  x <- grid1(n = grid[1], range = xlim, edge = TRUE)
  y <- grid1(n = grid[2], range = ylim, edge = TRUE)
  xy <- grid2(x, y)
  if (verbose) 
    cat("\n computing density and uncertainty over grid ...\n")
  Z <- densNuncer(modelName = "VVV", data = xy, parameters = parameters)
  lx <- length(x)
  ly <- length(y)
  CI <- type
  DU <- what
  TRANS <- transformation
  if (length(CI) > 1) 
    CI <- CI[1]
  if (length(DU) > 1) 
    DU <- DU[1]
  if (length(TRANS) > 1) 
    TRANS <- TRANS[1]
  switch(EXPR = DU, density = {
    zz <- matrix(Z$density, lx, ly)
    title2 <- "Density"
  }, uncertainty = {
    zz <- matrix(Z$uncertainty, lx, ly)
    title2 <- "Uncertainty"
  }, stop("what improperly specified"))
  switch(EXPR = TRANS, none = {
    title1 <- ""
  }, log = {
    zz <- logb(zz)
    title1 <- "log"
  }, sqrt = {
    zz <- sqrt(zz)
    title1 <- "sqrt"
  }, stop("transformation improperly specified"))
  
  switch(EXPR = CI, 
  contour = {
    title3 <- "Contour"
    if(is.null(levels)) levels <- pretty(zz, nlevels)
    contour(x = x, y = y, z = zz, levels = levels, xlab = xlab, 
        ylab = ylab, main = "", ...)
  }, 
  image = {
    title3 <- "Image"
    image(x = x, y = y, z = zz, xlab = xlab, ylab = ylab, 
        main = "", ...)
  }, 
  persp = {
    title3 <- "Perspective"
    persp(x = x, y = y, z = zz, 
          xlab = xlab, ylab = ylab, zlab = "Density",
          theta = 60, phi = 30, expand = 0.6, main = "", ...)
  }, stop("type improperly specified"))
  if (identify) {
    TITLE <- paste(c(title1, title2, title3, "Plot"), collapse = " ")
    title(TITLE)
  }
  invisible(list(x = x, y = y, z = zz))
}

uncerPlot <- function (z, truth=NULL, ...) 
{
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(pty = "m")
    uncer <- 1 - apply(z, 1, max)
    ord <- order(uncer)
    M <- max(uncer)
    plot(uncer[ord], ylab = "uncertainty", ylim = c(-(M/32), M), 
         xaxt = "n", xlab = "observations in order of increasing uncertainty", 
         type = "n")
    points(uncer[ord], pch = 15, cex = 0.5)
    lines(uncer[ord])
    abline(h = c(0, 0), lty = 3)
    if (!is.null(truth)) {
        truth <- as.numeric(as.factor(truth))
        n <- length(truth)
        result <- map(z)
        bad <- classError(result, truth)$misclassified
        if(length(bad)) 
          { for(i in bad) 
               { x <- (1:n)[ord == i]
                 lines(c(x, x), c(-(0.5/32), uncer[i]), lty = 1)
               }
          }
    }
    invisible()
}

decomp2sigma <- function(d, G, scale, shape, orientation = NULL, ...)
{
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
    sigma[,  , k] <- crossprod(t(orientation[,  , k]) * sqrt(scale[
      k] * shape[, k]))
  }
  structure(sigma, modelName = paste(c(scaleName, shapeName, orientName),
    collapse = ""))
}

grid1 <- function (n, range = c(0, 1), edge = TRUE) 
{
    if (any(n < 0 | round(n) != n)) 
        stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if (edge) {
        G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n - 
            1))
    }
    else {
        lj <- abs(diff(range))
        incr <- lj/(2 * n)
        G <- seq(from = min(range) + incr, to = max(range) - 
            incr, by = 2 * incr)
    }
    G
}

grid2 <- function (x, y) 
{
    lx <- length(x)
    ly <- length(y)
    xy <- matrix(0, nrow = lx * ly, ncol = 2)
    l <- 0
    for (j in 1:ly) {
        for (i in 1:lx) {
            l <- l + 1
            xy[l, ] <- c(x[i], y[j])
        }
    }
    xy
}
hypvol <- function (data, reciprocal = FALSE) 
{
    dimdat <- dim(data)
    oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 
        1
    if (oneD) {
        n <- length(as.vector(data))
        if (reciprocal) {
            ans <- 1/diff(range(data))
        }
        else {
            ans <- diff(range(data))
        }
        return(ans)
    }
    if (length(dimdat) != 2) 
        stop("data must be a vector or a matrix")
    data <- as.matrix(data)

    sumlogdifcol <- function(x) 
            sum(log(apply(x, 2, function(colm) diff(range(colm)))))

    bdvolog <- sumlogdifcol(data)
    pcvolog <- sumlogdifcol(princomp(data)$scores)

    volog <- min(bdvolog,pcvolog)

    if (reciprocal) {
       minlog <- log(.Machine$double.xmin)
        if (-volog < minlog) {
            warning("hypervolume smaller than smallest machine representable positive number")
            ans <- 0
        }
        else ans <- exp(-volog)
    }
    else {
        maxlog <- log(.Machine$double.xmax)
        if (volog > maxlog) {
            warning("hypervolume greater than largest machine representable number")
            ans <- Inf
        }
        else ans <- exp(volog)
    }

    ans
}

imputeData <- function(x, categorical=NULL, seed=NULL) 
{
  if (!exists("prelim.mix") || ! exists("em.mix") || !exists("da.mix") ||
    !exists("imp.mix") || !exists("rngseed") ) library(mix)

  fac <- apply( x, 2, is.factor)
  if (is.null(categorical)) {
   categorical <- fac
  }
  else {
   if (any(!categorical & fac)) {
     stop("x has a factor that is not designated as categorical")
   }
   if (any(categorical | !fac)) {
     warning("a categorical is not designated as a factor")
     for(i in which(categorical | !fac)) x[[i]] <- as.factor(x[[i]])
   }
  }

  # remove categorical variables and add dummy variable
  if (nocat <- !any(categorical)) {
    x <- cbind(as.factor(1),x)
    categorical <- c(TRUE, categorical)
  }

  ord <- c(which(categorical),which(!categorical))

  # do the imputations
  s <- prelim.mix(x[,ord],p=sum(categorical))
  if (is.null(seed)) seed <- runif(1,min=.Machine$integer.max/1024,max=
                                   .Machine$integer.max)
  rngseed(seed)   # set random number generator seed
  thetahat <- em.mix(s) # find ML estimate
  newtheta <- da.mix(s, thetahat, steps=100, showits=TRUE)
  ximp <- imp.mix(s, newtheta) # impute under newtheta

  if (nocat)  ximp[,-1] else ximp[,order(ord)]
}

map <- function(z, warn = TRUE, ...)
{
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
      warning(paste("no assignment to", paste(J[!K], collapse
         = ",")))
  }
  cl
}

orth2 <- function (n) 
{
  u <- rnorm(n)
  u <- u/vecnorm(u)
  v <- rnorm(n)
  v <- v/vecnorm(v)
  Q <- cbind(u, v - sum(u * v) * u)
  dimnames(Q) <- NULL
  Q
}

partconv <- function(x, consec = TRUE)
{
  n <- length(x)
  y <- numeric(n)
  u <- unique(x)
  if(consec) {
    # number groups in order of first row appearance
    l <- length(u)
    for(i in 1:l)
      y[x == u[i]] <- i
  }
  else {
    # represent each group by its lowest-numbered member
    for(i in u) {
      l <- x == i
      y[l] <- (1:n)[l][1]
    }
  }
  y
}

partuniq <- function(x)
{
  # finds the classification that removes duplicates from x
  charconv <- function(x, sep = "001")
  {
    if(!is.data.frame(x)) x <- data.frame(x)
    do.call("paste", c(as.list(x), sep = sep))
  }

  n <- nrow(x)
  x <- charconv(x)
  k <- duplicated(x)
  partition <- 1.:n
  partition[k] <- match(x[k], x)
  partition
}

shapeO <- function(shape, O, transpose = FALSE)
{
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
           as.logical(transpose),
           as.double(shape),
           O,
           as.integer(l),
           as.integer(dimO[3]),
           double(l * l),
           integer(1),
           PACKAGE="mclust")[[3]]
}

sigma2decomp <- function (sigma, G=NULL, tol=NULL, ...) 
{
    dimSigma <- dim(sigma)
    if (is.null(dimSigma)) 
        stop("sigma improperly specified")
    d <- dimSigma[1]
    if (dimSigma[2] != d) 
        stop("sigma improperly specified")
    l <- length(dimSigma)
    if (l < 2 || l > 3) 
        stop("sigma improperly specified")
    if (is.null(G)) {
        if (l == 2) {
            G <- 1
            sigma <- array(sigma, c(dimSigma, 1))
        }
        else {
            G <- dimSigma[3]
        }
    }
    else {
        if (l == 3 && G != dimSigma[3]) 
            stop("sigma and G are incompatible")
        if (l == 2 && G != 1) 
            sigma <- array(sigma, c(d,d,G))
    }
    decomp <- list(d = d, G = G, scale = rep(0, G), shape = matrix(0, 
        d, G), orientation = array(0, c(d, d, G)))
    for (k in 1:G) {
        ev <- eigen(sigma[, , k], symmetric = TRUE)
        temp <- log(ev$values)
        logScale <- sum(temp)/d
        decomp$scale[k] <- exp(logScale)
        decomp$shape[, k] <- exp(temp - logScale)
        decomp$orientation[, , k] <- ev$vectors
    }
    if (is.null(tol)) 
        tol <- sqrt(.Machine$double.eps)
    scaleName <- "V"
    shapeName <- "V"
    orientName <- "V"
    uniq <- function(x, tol = sqrt(.Machine$double.eps)) {
        abs(max(x) - min(x)) < tol
    }
    if (uniq(decomp$scale)) {
        decomp$scale <- decomp$scale[1]
        scaleName <- "E"
    }
    if (all(apply(decomp$shape, 1, uniq, tol = tol))) {
        decomp$shape <- decomp$shape[, 1]
        if (all(uniq(decomp$shape, tol = tol))) {
            shapeName <- "I"
            decomp$shape <- rep(1, d)
        }
        else {
            shapeName <- "E"
        }
    }
    if (all(apply(matrix(decomp$orientation, nrow = d * d, ncol = G), 
        1, uniq, tol = tol))) {
        decomp$orientation = decomp$orientation[, , 1]
        if (all(apply(cbind(decomp$orientation, diag(d)), 1, 
            uniq, tol = tol))) {
            orientName <- "I"
            decomp$orientation <- NULL
        }
        else {
            orientName <- "E"
        }
    }
    modelName <-  paste(c(scaleName, shapeName, orientName), collapse = "")
    c(list(modelName = modelName, decomp))
}

traceW <- function(x)
{
  # sum(as.vector(sweep(x, 2, apply(x, 2, mean)))^2)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  .Fortran("mcltrw",
          as.double(x),
          as.integer(n),
          as.integer(p),
          double(p),
          double(1),
          PACKAGE = "mclust")[[5]]
}

unchol <- function(x, upper = NULL)
{
  if(is.null(upper)) {
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
           as.logical(upper),
           x,
           as.integer(nrow(x)),
           as.integer(ncol(x)),
           integer(1),
           PACKAGE="mclust")[[2]]
}

unmap <- function(classification, groups=NULL, noise=NULL, ...)
{
  # converts a classification to conditional probabilities
  # classes are arranged in sorted order unless groups is specified
  # if a noise indicator is specified, that column is placed last
  n <- length(classification)
        u <- sort(unique(classification))
        if (is.null(groups)) {
    groups <- u
        }
        else {
          if (any(match( u, groups, nomatch = 0) == 0)) 
            stop("groups incompatible with classification")
           miss <- match( groups, u, nomatch = 0) == 0
        }
        cgroups <- as.character(groups)
        if (!is.null(noise)) {
          noiz <- match( noise, groups, nomatch = 0)
          if (any(noiz == 0)) stop("noise incompatible with classification")
          groups <- c(groups[groups != noise],groups[groups==noise])
          noise <- as.numeric(factor(as.character(noise), levels = unique(groups)))
        }
        groups <- as.numeric(factor(cgroups, levels = unique(cgroups)))
        classification <- as.numeric(factor(as.character(classification), levels = unique(cgroups)))
        k <- length(groups) - length(noise)
        nam <- levels(groups)
  if(!is.null(noise)) {
          k <- k + 1
          nam <- nam[1:k]
          nam[k] <- "noise"
        }
  z <- matrix(0, n, k, dimnames = c(names(classification),nam))
  for(j in 1:k)
    z[classification == groups[j], j] <- 1
  z
}

vecnorm <- function (x, p = 2) 
{
    if (is.character(p)) {
        if (charmatch(p, "maximum", nomatch = 0) == 1) 
            p <- Inf
        else if (charmatch(p, "euclidean", nomatch = 0) == 1) 
            p <- 2
        else stop("improper specification of p")
    }
    if (!is.numeric(x) && !is.complex(x)) 
        stop("mode of x must be either numeric or complex")
    if (!is.numeric(p)) 
        stop("improper specification of p")
    if (p < 1) 
        stop("p must be greater than or equal to 1")
    if (is.numeric(x)) 
        x <- abs(x)
    else x <- Mod(x)
    if (p == 2) 
        return(.Fortran("d2norm", as.integer(length(x)), as.double(x), 
            as.integer(1), double(1), PACKAGE = "mclust")[[4]])
    if (p == Inf) 
        return(max(x))
    if (p == 1) 
        return(sum(x))
    xmax <- max(x)
    if (!xmax) 
        xmax <- max(x)
    if (!xmax) 
        return(xmax)
    x <- x/xmax
    xmax * sum(x^p)^(1/p)
}

"[.mclustBIC" <- function (x, i, j, drop = FALSE) 
{
  ATTR <- attributes(x)[c("G", "modelNames", "prior", "control", 
                          "initialization", "Vinv", "warn", "n", "d", 
                          "oneD", "returnCodes", "class")]
  oldClass(x) <- NULL
  x <- NextMethod("[")
  if (is.null(dim(x))) return(x)
  ATTR$G <- as.numeric(dimnames(x)[[1]])
  ATTR$modelNames <- dimnames(x)[[2]]
  ATTR$returnCodes <- ATTR$returnCodes[dimnames(x)[[1]],dimnames(x)[[2]],
                                        drop=FALSE]
  do.call("structure", c(list(.Data = x), ATTR))
}

"[.mclustDAtest" <- function (x, i, j, drop = FALSE) 
{
  clx <- oldClass(x)
  oldClass(x) <- NULL
  NextMethod("[")
}

bic <- function(modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)
  if (G == 0) {
          ## one cluster case
          if(!noise) stop("undefined model")
          nparams <- 1
  }
  else {
          nparams <- nVarParams(modelName, d, G) + G*d
          if(!equalPro)
                  nparams <- nparams + (G - 1)
          if(noise)
                  nparams <- nparams + 2
  }
  2 * loglik - nparams * logb(n)
}


checkModelName <- function(modelName)
{
  switch(EXPR = modelName,
         E = ,
         V = ,
         EII  = ,
         VII = ,
         EEI  = ,
         VEI = ,
         EVI = ,
         VVI = ,
         EEE = ,
         VEE = ,
         EVE = ,
         VVE = ,
         EEV = ,
         VEV = ,
         VVV = TRUE,
         stop("invalid model name"))
}

em <- function(modelName, data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("em", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

estep <- function(modelName, data, parameters, warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("estep", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

hc <- function(modelName, data, ...)
{
  switch(EXPR = modelName,
         E = ,
         V = ,
         EII  = ,
         VII = ,
         EEE = ,
         VVV = TRUE,
         stop("invalid model name for hierarchical clustering"))
  funcName <- paste("hc", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc[[2]] <- NULL
  eval(mc, parent.frame())
}

mclustVariance <- function(modelName, d=NULL, G=2) 
{
 x <- -1
 if (nchar(modelName) == 1) {
   if (!is.null(d) && d != 1)  stop("modelName and d are incompatible")
   varList <- switch(EXPR = modelName,
                     "X" = list(sigmasq = x),
                     "E" = list(sigmasq = x),
                     "V" = list(sigmasq = rep(x,G)),
                     stop("modelName not recognized"))
 }
 else {
   if (nchar(modelName) != 3) stop("modelName is misspecified")
   if (is.null(d)) d <- 3
   varList <- switch(EXPR = modelName,
                     "XII" = list(sigmasq = x),
                     "EII" = list(sigmasq = x),
                     "VII" = list(sigmasq = rep(x,G)),
                     "XXI" = list(scale = x, shape = rep(x,d)),
                     "EEI" = list(scale = x, shape = rep(x,d)),
                     "EVI" = list(scale = x, shape = matrix(x,d,G)),
                     "VEI" = list(scale = rep(x,G), shape = rep(x,d)),
                     "VVI" = list(scale = rep(x,G), shape = matrix(x,d,G)),
                     "XXX" = {M <- matrix(x,d,d); M[row(M) > col(M)] <- 0;
                              list(cholSigma = M)},
                     "EEE" = {M <- matrix(x,d,d); M[row(M) > col(M)] <- 0;
                              list(cholSigma = M)},
                     "EEV" = list(scale = x, shape = rep(x,d),
                                  orientation = array(x,c(d,d,G))),
                     "VEV" = list(scale = x, shape = matrix(x,d,G),
                                  orientation = array(x,c(d,d,G))),
                     "VVV" = {A <- array(x,c(d,d,G)); 
                              I <- row(A[,,1]) > col(A[,,1])
                              for (k in 1:G) A[,,k][I] <- 0
                              list(cholsigma = A)},
                     stop("modelName not recognized"))
  }
 c(modelName = modelName, d=d, G=G, varList)
}

me <- function(modelName, data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("me", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

mstep <- function(modelName, data, z, prior = NULL, warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("mstep", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

mvn <- function(modelName, data, prior = NULL, warn = NULL, ...)
{
    modelName <- switch(EXPR = modelName,
                        E = "X",
                        V = "X",
                        X  =  "X",
                        Spherical = "XII",
                        EII  = "XII",
                        VII = "XII",
                        XII = "XII",
                        Diagonal =  "XXI",
                        EEI  = "XXI",
                        VEI = "XXI",
                        EVI = "XXI",
                        VVI = "XXI",
                        XXI = "XXI",
                        Ellipsoidal = "XXX",
                        EEE = "XXX",
                        VEE = "XXX",
                        EVE = "XXX",
                        VVE = "XXX",
                        EEV = "XXX",
                        VEV = "XXX",
                        VVV = "XXX",
                        XXX = "XXX",
                        stop("invalid model name"))

  funcName <- paste("mvn", modelName, sep = "")
  mc <- match.call()
  mc[[1]] <- as.name(funcName)
  mc[[2]] <- NULL
  out <- eval(mc, parent.frame())
  varnames <- colnames(as.matrix(data))
  if(!all(is.null(varnames)))
    { rownames(out$parameters$mean) <- varnames
      dimnames(out$parameters$variance$Sigma) <- list(varnames, varnames)
      dimnames(out$parameters$variance$sigma) <- list(varnames, varnames, NULL)
    }
  return(out)
}

nVarParams <- function(modelName, d, G)
{
  modelName <- switch(EXPR = modelName, 
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)

  r <- (d*(d-1))/2
  s <- (d*(d+1))/2 # s = r + d
  switch(EXPR = modelName,
         E = 1,
         V = G,
         EII = 1,
         VII = G,
         EEI = d,
         VEI = G+(d-1),
         EVI = 1 + G*(d-1),
         VVI = G*d,
         EEE = s,
         VEE = G + (d-1) + r,
         EVE = 1 + (d-1) + r,
         VVE = G + G*(d-1) + r,
         EEV = 1 + (d-1) + G*r,
         VEV = G + (d-1) + G*r,
         VVV = G * s,
         stop("invalid model name"))
}

sim <- function(modelName, parameters, n, seed = NULL, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)
  funcName <- paste("sim", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

cdensVEI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VEI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
  temp <- .Fortran("esvei",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- NA
                ret <- -1
  }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],dimnames(mu)[[2]])
  structure(z, logarithm = logarithm, modelName = "VEI",
                  WARNING = WARNING, returnCode = ret)
}

emVEI <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepVEI(data, parameters = parameters, warn = warn)$z  
  meVEI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepVEI <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
                Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  } 
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VEI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
  temp <- .Fortran("esvei",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "VEI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

meVEI <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEI", d = p, G = G, 
                                 scale = rep(NA,G), shape = rep(NA,p)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevei",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax),
      as.double(control$tol),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(p),
      double(K),
      double(G),
      double(p),
      double(p * G),
                        PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VEI"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("meveip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax),
      as.double(control$tol),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(p),
      double(K),
      double(G),
      double(p),
      double(p * G),
                        PACKAGE = "mclust")[11:18]
  }
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
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= control$itmax[2]) {
      WARNING <- "inner iteration limit reached"
      warning(WARNING)
      inner <-  - inner
      ret <- 2
    }
    else if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  attr(info, "inner") <- c(iterations = inner, error = inerr)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
  variance <- list(modelName = "VEI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepVEI <- function(data, z, prior = NULL, warn = NULL, control = NULL,...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
              WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEI", d = p, G = G, 
                                 scale = rep(NA,G), shape = rep(NA,p)) 
               parameters <- list(pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))

  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
        if (is.null(control)) control <- emControl()
  itmax <- if(length(control$itmax) == 1) control$itmax else control$
      itmax[2]
  tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
  if(is.null(prior)) {
    temp <- .Fortran("msvei",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.integer(itmax),
      as.double(tol),
      double(p * G),
      double(G),
      double(p),
      double(G),
      double(G),
      double(p),
      double(p * G),
                        PACKAGE = "mclust")[6:11]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VEI"), prior[names(
      prior) != "functionName"]))
    temp <- .Fortran("msveip",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      as.integer(itmax),
      as.double(tol),
      double(p * G),
      double(G),
      double(p),
      double(G),
      double(G),
      double(p),
      double(p * G),
                        PACKAGE = "mclust")[10:15]
  }
  inner <- temp[[1]]
  inerr <- temp[[2]]
  mu <- matrix(temp[[3]], p, G)
  scale <- temp[[4]]
  shape <- temp[[5]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  pro <- temp[[6]]
  WARNING <- NULL
  if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!
    c(scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn)
      warning(WARNING)
    mu[] <- pro[] <- shape <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
                ret <- -1
  }
  else {
                ret <- 0
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= itmax) {
      WARNING <- "inner iteration limit reached"
      warning(WARNING)
      inner <-  - inner
    }
  }
     info <- c(iterations = inner, error = inerr)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "VEI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control),
                  info = info, WARNING = WARNING, returnCode = ret)

}

simVEI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VEI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(rtscale[k] * rtshape), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VEI")
}

cdensV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
        mu <- drop(parameters$mean)
  G <- length(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")])))
            || any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "V", 
                                 WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
        if(is.null(sigmasq))
                stop("variance parameters are missing")
        if(any(sigmasq < 0))
                stop("sigma-squared is negative")
  if(any(!sigmasq)) {
                WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "V", 
                                 WARNING = WARNING, returnCode = 9))
  }
        if (length(sigmasq) == 1) sigmasq <- rep(sigmasq,G)
  temp <- .Fortran("es1v",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(-1),
    as.integer(n),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[8:9]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(names(data),NULL)
  structure(z, logarithm = logarithm, modelName = "V",
                     WARNING = WARNING, returnCode = ret)
}

emV <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepV(data, parameters = parameters, warn = warn)$z  
  meV(data, z = z, prior = prior, control = control, 
            Vinv = parameters$Vinv, warn = warn)
}

estepV <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
        pro <- parameters$pro 
        pro <- pro/sum(pro)
  l <- length(pro)
        mu <- drop(parameters$mean)
  G <- length(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv 
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")])))
            || any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
                return(structure(list(modelName = "V", n=n, d=1, G=G, z=z,
                       parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
    if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
                return(structure(list(modelName = "V", n=n, d=1, G=G, z=z,
                       parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("es1v",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(pro),
    as.integer(n),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(names(data),NULL)
        structure(list(modelName = "V", n = n, d = 1, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

cdensVEV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VEV", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape) ||
                   is.null(parameters$variance$orientation)) 
          stop("variance parameters are missing")
  temp <- .Fortran("esvev",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(aperm(parameters$variance$orientation,c(2,1,3))),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(p),
    double(p),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- NA
                ret <- -1
  }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL) 
  structure(z, logarithm = logarithm, modelName = "VEV",
                  WARNING = WARNING, returnCode = ret)
}

emVEV <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepVEV(data, parameters = parameters, warn = warn)$z  
  meVEV(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepVEV <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
       if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VEV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape) ||
                    is.null(parameters$variance$orientation)) 
           stop("variance parameters are missing")
  temp <- .Fortran("esvev",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(aperm(parameters$variance$orientation,c(2,1,3))),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(p),
    double(p),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1  
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL) 
        structure(list(modelName = "VEV", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)

}

meVEV <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEV", d = p, G = G, 
           scale=rep(NA,G), shape=rep(NA,p), orientation=array(NA,c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters,
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevev",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax),
      as.double(control$tol),
      as.double(control$eps),
      as.integer(lwork),
      double(p * G),
      double(G),
      double(p),
      double(p * p * G),
      double(K),
      double(lwork),
      double(p),
                        PACKAGE = "mclust")[7:16]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VEV"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("mevevp",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax),
      as.double(control$tol),
      as.double(control$eps),
      as.integer(lwork),
      double(p * G),
      double(G),
      double(p),
      double(p * p * G),
      double(K),
      double(lwork),
      double(p),
                        PACKAGE = "mclust")[11:20]
  }
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
  O <- aperm( array(temp[[9]], c(p, p, G)), c(2,1,3))
  pro <- temp[[10]]
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
    }
    warning(WARNING)
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- shapeO(shape, O, transpose = FALSE)
    sigma <- sweep(sigma, MARGIN = 3, STATS = scale, FUN = "*")
    if(inner >= control$itmax[2]) {
      WARNING <- "inner iteration limit reached"
      warning(WARNING)
      inner <-  - inner
      ret <- 2
    }
    else if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- structure(c(iterations = its, error = err), inner = c(
    iterations = inner, error = inerr))
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(O) <- 
           list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
##  Sigma = scale * O %*% diag(shape) %*% t(O)
  variance <- list(modelName = "VEV", d = p, G = G, sigma = sigma, 
                        scale = scale, shape = shape, orientation = O)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance) 
  structure(list(modelName = "VEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepVEV <- function(data, z, prior = NULL, warn = NULL, control = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEV", d = p, G = G, 
      scale = rep(NA,G), shape = rep(NA,p), orientation = array(NA,c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))

    WARNING <- "z is missing"
    warning(WARNING)
    return(structure(list(n = n, d = p, G = G, mu = matrix(NA,
      p, G), sigma = array(NA, c(p, p, G)), decomp = list(
      d = p, G = G, scale = rep(NA, G), shape = rep(NA, p),
      orientation = array(NA, c(p, p, G))), pro = rep(NA,
      G), modelName = "VEV", prior = prior), WARNING = 
      WARNING))
  }
  #  shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop(
      "improper specification of z")
        if (is.null(control)) control <- emControl()
  itmax <- if(length(control$itmax) == 1) control$itmax else control$
      itmax[2]
  tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  if(is.null(prior)) {
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
      double(p * G),
      double(G),
      double(p),
      double(p * p * G),
      double(G),
                        PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VEV"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("msvevp",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      double(lwork),
      as.integer(lwork),
      as.integer(itmax),
      as.double(tol),
      double(p * G),
      double(G),
      double(p),
      double(p * p * G),
      double(G),
                        PACKAGE = "mclust")[11:18]
  }
  lapackSVDinfo <- temp[[1]]
  inner <- temp[[2]]
  inerr <- temp[[3]]
  mu <- matrix(temp[[4]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[5]]
  shape <- temp[[6]]
  O <- aperm(array(temp[[7]], c(p, p, G)),c(2,1,3))
        pro <- temp[[8]] 
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
      warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
      warning(WARNING)
    }
    O[] <- shape[] <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
                ret <- -9
  }
  else if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(
    !c(scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
                ret <- -1
  }
  else {
    sigma <- sweep(shapeO(shape, O, transpose = FALSE), MARGIN = 3,
      STATS = scale, FUN = "*")
    if(inner >= itmax) {
      WARNING <- "inner iteration limit reached"
      warning(WARNING)
      inner <-  - inner
    }
                ret <- 2
  }
  info <- c(iteration = inner, error = inerr)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(O) <-
           list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
        variance <- list(modelName = "VEV", d = p, G = G, sigma = sigma, 
                         scale = scale, shape = shape, orientation = O)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control),
                  info = info, WARNING = WARNING, returnCode = ret)

}

simVEV <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VEV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    sss <- rtscale[k] * rtshape
    cholSigma <- t(parameters$variance$orientation[,  , k]) * sss
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      cholSigma, MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VEV")
}

hcV <- function(data, partition, minclus = 1, alpha = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
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
                PACKAGE = "mclust")[c(1, 3, 8)]
  temp[[1]] <- temp[[1]][1:m]
  temp[[2]] <- temp[[2]][1:m]
  temp[[3]] <- temp[[3]][1:m]
        change <- temp[[3]]
  structure(rbind(temp[[1]], temp[[2]]),   initialPartition = partition, 
                  dimensions = n, modelName = "V",
      call = match.call())
}

cdensVII <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VII", 
                                 WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq 
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VII", 
                                 WARNING = WARNING, returnCode = 9))
  }
  temp <- .Fortran("esvii",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[9:10]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0 
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "VII",
                  WARNING = WARNING, returnCode = ret)
}

emVII <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepVII(data, parameters = parameters, warn = warn)$z  
  meVII(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepVII <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
        mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VII", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq 
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VII", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("esvii",
    as.double(data),
    as.double(mu),
    as.double(sigmasq),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "VII", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

hcVII <- function(data, partition, minclus = 1, alpha = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
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
                PACKAGE = "mclust")[c(1, 10)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
        change <- temp[[2]]
  structure(t(temp[[1]]), initialPartition = partition, 
                  dimensions = dimdat, modelName = "VII", 
                  call = match.call())
}

meVII <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevii",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(K),
                        PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VII"), prior[names(prior) !=
      "functionName"]))
    storage.mode(z) <- "double"
    temp <- .Fortran("meviip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(K),
                        PACKAGE = "mclust")[c(11:17, 10)]
  }
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
           any(sigmasq <= max(control$eps, 0))) {
    WARNING <- "sigma-squared falls below threshold"
    if (warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if (warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if (warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(rep(sigmasq[k], p))
    if(its >= control$itmax[1]) {
      warning("iteration limit reached")
      WARNING <- "iteration limit reached"
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
  variance <- list(modelName = "VII", d = p, G = G, 
                         sigma = sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

meVVI <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
         WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVI", d = p, G = G, 
                                 scale = rep(NA,G), shape = matrix(NA,p,G)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
        }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevvi",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(p * G),
      double(K),
                        PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VVI"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("mevvip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(p * G),
      double(K),
                        PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if (warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
      FUN = "*"), 2, diag), c(p, p, G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
  variance <- list(modelName = "VVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepVII <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal number of observations")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                                WARNING = WARNING, returnCode = 9))
        }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("msvii",
      as.double(data),
      z,
      as.integer(n),
      as.integer(p),
      as.integer(G),
      double(p * G),
      double(G),
      double(G),
                        PACKAGE = "mclust")[6:8]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VII"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("msviip",
      as.double(data),
      z,
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p * G),
      double(G),
      double(G),
                        PACKAGE = "mclust")[10:12]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  sigma <- array(0, c(p, p, G))
  for(k in 1:G)
    sigma[,  , k] <- diag(rep(sigmasq[k], p))
  WARNING <- NULL
  if(any(sigmasq > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "VII", d = p, G = G, 
                         sigma = sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters),
                  WARNING = WARNING, returnCode = ret)
}

simVII <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VII"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  sigmasq <- parameters$variance$sigmasq
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(rep(sqrt(sigmasq[k]), d)), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VII")
}

meV <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal length of data")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
         WARNING <- "z is missing"
         if (warn) warning(WARNING)
               variance <- list(modelName = "V", d=1, G=G, sigmasq = rep(NA,G))
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="V", prior=prior, n=n, d=1, G=G,
                  z=z, parameters=parameters, control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("me1v",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(G),
      as.double(if(is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(G),
      double(G),
      double(K),
                        PACKAGE = "mclust")[6:12]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "V"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("me1vp",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(G),
      double(G),
      double(K),
                        PACKAGE = "mclust")[c(10:16, 9)]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- temp[[5]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  ## logpost <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || any(sigmasq <= max(
    control$eps, 0))) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else if(its >= control$itmax[1]) {
    WARNING <- "iteration limit reached"
    if (warn) warning(WARNING)
    its <-  - its
    ret <- 1
  }
  else ret <- 0
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(names(data),NULL)
        variance = list(modelName = "V", d = 1, G = G, 
                        sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "V", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik),
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepV <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
# number of groups 
  G <- dimz[2]
  ##
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "V", d=1, G=G, sigmasq=rep(NA,G))
               parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="V", prior=prior, n=n, d=1, G=G,
                                     z=z, parameters=parameters), 
                                WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("ms1v",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(G),
      double(G),
      double(G),
      double(G),
                        PACKAGE = "mclust")[5:7]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "V"), prior[names(prior) !=
      "functionName"]))
    storage.mode(z) <- "double"
    temp <- .Fortran("ms1vp",
      as.double(data),
      z,
      as.integer(n),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(G),
      double(G),
      double(G),
                        PACKAGE = "mclust")[9:11]
  }
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  WARNING <- NULL
  if(any(sigmasq > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
                mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
                print(G)
                print(sigmasq)
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(names(data),NULL)
        variance = list(modelName = "V", d = 1, G = G, 
                        sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "V", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters),
                  WARNING = WARNING, returnCode = ret)
}

simV <- function(parameters, n, seed = NULL, ...)
{
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, 2), modelName = "V"))
  }
  if(!is.null(seed))
    set.seed(seed)
  mu <- parameters$mean
  G <- length(mu)
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- rep(0, n)
  sd <- sqrt(parameters$variance$sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd[k])
  }
  structure(cbind(group = clabels, "1" = x), modelName = "V")
}

cdensVVI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mu", "variance")]))) ||
            any(is.null(parameters[c("pro", "mu", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VVI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
  temp <- .Fortran("esvvi",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- NA
                ret <- -1
  }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)  
  structure(z, logarithm = logarithm, modelName = "VVI",
                  WARNING = WARNING, retrinCode = ret)
}

emVVI <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepVVI(data, parameters = parameters, warn = warn)$z  
  meVVI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepVVI <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
        pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
                Vinv <- parameters$Vinv
    if (is.null(Vinv) || Vinv <= 0) 
                  Vinv <- hypvol(data, reciprocal = TRUE)
  }
        if(any(is.na(unlist(parameters[c("pro", "mu", "variance")]))) ||
            any(is.null(parameters[c("pro", "mu", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VVI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
  temp <- .Fortran("esvvi",
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$scale),
    as.double(parameters$variance$shape),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL) 
        structure(list(modelName = "VVI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

meVVI <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
         WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVI", d = p, G = G, 
                                 scale = rep(NA,G), shape = matrix(NA,p,G)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
        }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevvi",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(p * G),
      double(K),
                        PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VVI"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("mevvip",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(G),
      double(p * G),
      double(K),
                        PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if (warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
      FUN = "*"), 2, diag), c(p, p, G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
  variance <- list(modelName = "VVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepVVI <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))

  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("msvvi",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      double(p * G),
      double(G),
      double(p * G),
      double(G),
                        PACKAGE = "mclust")[6:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VVI"), prior[names(
      prior) != "functionName"]))
    temp <- .Fortran("msvvip",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p * G),
      double(G),
      double(p * G),
      double(G),
                        PACKAGE = "mclust")[10:13]
  }
  mu <- matrix(temp[[1]], p, G)
  scale <- temp[[2]]
  shape <- matrix(temp[[3]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  WARNING <- NULL
  if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!
    c(scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn)
      warning(WARNING)
    mu[] <- pro[] <- shape <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
                ret <- -1
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
      FUN = "*"), 2, diag), c(p, p, G))
                ret <- 0
  }
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "VVI", d = p, G = G, 
                         sigma = sigma, sigmasq = scale, 
                         scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters),
                  WARNING = WARNING, returnCode = ret)
}

simVVI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VVI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(!all(dim(rtshape) == dim(mu)))
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(rtscale[k] * rtshape[, k]), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VVI")
}

cdensVVV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VVV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholsigma))
          stop("variance parameters are missing")
  temp <- .Fortran("esvvv",
    as.logical(1),
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$cholsigma),
    as.double(-1),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(-1),
    double(p),
    double(1),
    double(n * G),
                PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, G)
        WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      if (warn) warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      if (warn) warning(WARNING)
    }
    z[] <- NA
                ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- NA
                ret <- -1
  }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "VVV",
                  WARNING = WARNING, returnCode = ret)
}

emVVV <- function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
  z <- estepVVV(data, parameters = parameters, warn = warn)$z  
  meVVV(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

estepVVV <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VVV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholsigma))
          stop("variance parameters are missing")
  temp <- .Fortran("esvvv",
          as.logical(1),
    as.double(data),
    as.double(mu),
    as.double(parameters$variance$cholsigma),
    as.double(pro),
    as.integer(n),
    as.integer(p),
    as.integer(G),
    as.double(if (is.null(Vinv)) -1 else Vinv),
    double(p),
    double(1),
    double(n * K),
                PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, K)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      warning(WARNING)
    }
    z[] <- loglik <- NA
                ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if (warn) warning(WARNING)
    z[] <- loglik <- NA
                ret <- -1
  }
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "VVV", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

hcVVV <- function(data, partition, minclus = 1, alpha = 1, beta = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
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
  #  dp <- duplicated(partition)
  #x[c((1:n)[!dp],(1:n)[dp]), ], 
  #as.integer(c(partition[!dp], partition[dp])), 
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
                PACKAGE = "mclust")[c(1, 14)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
  temp[[2]] <- temp[[2]][1:m]
        change <- temp[[2]] 
  structure(t(temp[[1]]), initialPartition = partition, 
                  dimensions = dimdat, modelName = "VVV", 
                  call = match.call())
}

meVVV <- function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
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
        if (!is.null(Vinv))  {
    G <- K - 1
    if(Vinv <= 0)Vinv <- hypvol(data, reciprocal = TRUE)
  }
        else G <- K
  if(all(is.na(z))) {
            WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVV", d = p, G = G, 
                sigma = array(NA, c(p,p,G)), cholsigma = array(NA, c(p,p,G))) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VVV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevvv",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(p * p * G),
      double(K),
      double(p),
      double(p*p),
                        PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VVV"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("mevvvp",
      as.logical(control$equalPro),
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(if (is.null(Vinv)) -1 else Vinv),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      z,
      as.integer(control$itmax[1]),
      as.double(control$tol[1]),
      as.double(control$eps),
      double(p * G),
      double(p * p * G),
      double(K),
      double(p),
      double(p*p),
                        PACKAGE = "mclust")[c(11:17, 10)]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- array(temp[[6]], c(p, p, G))
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(cholsigma, 3, unchol, upper = TRUE), 
                               c(p,p,G))
    if(its >= control$itmax[1]) {
      warning("iteration limit reached")
      WARNING <- "iteration limit reached"
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = abs(err))
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(cholsigma) <-
              list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
        variance <- list(modelName = "VVV", d = p, G = G,
                          sigma = sigma, cholsigma = cholsigma)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance) 
  structure(list(modelName = "VVV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

mstepVVV <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if (is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVV", d = p, G = G, 
                      sigma <- array(NA, c(p,p, G)), 
                      cholsigma = array(NA, c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VVV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
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
      double(p * p),
                        PACKAGE = "mclust")[7:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = G, modelName = "VVV"), 
                        prior[names(prior) != "functionName"]))
    temp <- .Fortran("msvvvp",
      as.double(data),
      as.double(z),
      as.integer(n),
      as.integer(p),
      as.integer(G),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      double(p),
      double(p * G),
      double(p * p * G),
      double(G),
      double(p * p),
                        PACKAGE = "mclust")[11:13]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- array(temp[[2]], c(p, p, G))
        pro <- temp[[3]]
  WARNING <- NULL
  if(any(c(mu, cholsigma) > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- sigma[] <- cholsigma[] <- NA
                ret <- -1
  }
  else {
    sigma <- array(apply(cholsigma, 3, unchol, upper = TRUE), 
                                     c(p,p,G))
                ret <- 0
  }
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        dimnames(cholsigma) <- list(dimnames(data)[[2]], 
                                    dimnames(data)[[2]], NULL)
        variance <- list(modelName = "VVV", d = p, G = G, 
                         sigma = sigma, cholsigma= cholsigma)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VVV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)

}

simVVV <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
    parameters[c("mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VVV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  if(is.null(cholsigma <- parameters$variance$cholsigma)) {
    if(is.null(sigma <- parameters$variance$sigma)) {
      stop("variance parameters must inlcude either sigma or cholsigma"
        )
    }
    cholsigma <- apply(sigma, 3, chol)
    for(k in 1:ncol(cholsigma))
      sigma[,  , k] <- cholsigma[, k]
    cholsigma <- sigma
  }
  if(dim(cholsigma)[3] != G)
    stop("variance incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
      ncol = d) %*% cholsigma[,  , k], MARGIN = 2, STATS = mu[
      , k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VVV")
}

mvnXII <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD)
    stop("for multidimensional data only")
  if(length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(prior)) {
    temp <- .Fortran("mvnxii",
      as.double(data),
      as.integer(n),
      as.integer(p),
      double(p),
      double(1),
      double(1),
                        PACKAGE = "mclust")[4:6]
    logpost <- NULL
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = 1, modelName = "XII"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("mnxiip",
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p),
      double(1),
      double(1),
                        PACKAGE = "mclust")[c(8:10, 7)]
    logpost <- temp[[4]]
  }
  mu <- temp[[1]]
  sigmasq <- temp[[2]]
  loglik <- temp[[3]]
  Sigma <- sigmasq * diag(p)
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    loglik <- NA
    ret <- -1
  }
  variance <- list(modelName = "XII", d = p, G = 1,  
                   sigmasq   = sigmasq, Sigma = Sigma, 
                   sigma = array(Sigma, c(p, p, 1)), scale = sigmasq)
  parameters <- list(pro = 1, mean = matrix(mu, ncol = 1), variance = variance) 
  structure(list(modelName = "XII", prior = prior, n = n, d = p, G = 1, 
                 parameters = parameters, loglik = loglik), 
                 WARNING = WARNING, returnCode = ret) 
}

mvnX <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(prior)) {
    temp <- .Fortran("mvn1d",
      as.double(data),
      as.integer(n),
      double(1),
      double(1),
      double(1),
                        PACKAGE = "mclust")[3:5]
    logpost <- NULL
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = 1, modelName = "X"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("mvn1p",
      as.double(data),
      as.integer(n),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(1),
      double(1),
      double(1),
                        PACKAGE = "mclust")[c(7:9, 6)]
    logpost <- temp[[4]]
  }
  mu <- temp[[1]]
  sigmasq <- temp[[2]]
  loglik <- temp[[3]]
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    loglik <- NA
    ret <- -1
  }
  variance = list(modelName= "X", d = 1, G = 1, sigmasq = sigmasq)
  parameters <- list(pro = 1, mean = mu, variance = variance)
  structure(list(modelName = "X", prior = prior, n = n, d = 1, G = 1, 
                 parameters = parameters, loglik = loglik),
                  WARNING = WARNING, returnCode = ret) 
}

mvnXXI <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD)
    stop("for multidimensional data only")
  if(length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(prior)) {
    temp <- .Fortran("mvnxxi",
      as.double(data),
      as.integer(n),
      as.integer(p),
      double(p),
      double(1),
      double(p),
      double(1),
                        PACKAGE = "mclust")[4:7]
    logpost <- NULL
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = 1, modelName = "XXI"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("mnxxip",
      as.double(data),
      as.integer(n),
      as.integer(p),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(priorParams$scale),
      as.double(priorParams$dof),
      double(p),
      double(1),
      double(p),
      double(1),
                        PACKAGE = "mclust")[c(8:11, 7)]
    logpost <- temp[[5]]
  }
  mu <- temp[[1]]
  scale <- temp[[2]]
  shape <- temp[[3]]
  loglik <- temp[[4]]
  Sigma <- diag(scale * shape)
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    loglik <- NA
    ret <- -1
  }
  variance <- list(modelName = "XXI", d = p, G = 1, 
                   Sigma = Sigma, sigma = array(Sigma, c(p, p, 1)),
                   scale = scale, shape = shape)
  parameters <- list(pro = 1, mean = matrix(mu, ncol = 1), 
                     variance = variance)
  structure(list(modelName = "XXI", prior = prior, n = n, d = p, G = 1, 
                 parameters = parameters, loglik = loglik),
                 WARNING = WARNING, returnCode = ret) 
}

mvnXXX <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD)
    stop("for multidimensional data only")
  if(length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(prior)) {
    temp <- .Fortran("mvnxxx",
      as.double(data),
      as.integer(n),
      as.integer(p),
      double(p),
      double(p * p),
      double(1),
      PACKAGE = "mclust")[c(4:6)]
    logpost <- NULL
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
      data, G = 1, modelName = "XXX"), prior[names(prior) !=
      "functionName"]))
    temp <- .Fortran("mnxxxp",
      as.double(data),
      as.integer(n),
      as.integer(p),
      double(p),
      as.double(priorParams$shrinkage),
      as.double(priorParams$mean),
      as.double(if(any(priorParams$scale != 0)) chol(priorParams$
          scale) else priorParams$scale),
      as.double(priorParams$dof),
      double(p),
      double(p * p),
      double(1),
      PACKAGE = "mclust")[c(9:11, 8)]
    logpost <- temp[[4]]
  }
  mu <- temp[[1]]
  cholSigma <- matrix(temp[[2]], p, p)
  Sigma <- unchol(cholSigma, upper = TRUE)
  loglik <- temp[[3]]
## Sigma = t(cholSigma) %*% cholSigma
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn)
      warning(WARNING)
    loglik <- NA
    ret <- -1
  }
  variance <- list(modelName = "XXX", d = p, G = 1,
                   Sigma = Sigma, cholSigma = cholSigma, 
                   sigma = array(Sigma, c(p, p, 1))) 
  parameters <- list(pro = 1, mean = matrix(mu, ncol = 1), 
                     variance = variance)
  structure(list(modelName = "XXX", prior = prior, n = n, d = p, G = 1,
                 parameters = parameters, loglik = loglik), 
                 WARNING = WARNING, returnCode = ret)
}
