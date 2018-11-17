# functions to be removed??

EMclust <- function(data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = list(hcPairs=NULL, subset=NULL, noise=NULL), Vinv = NULL, warn = FALSE, x = NULL, ...)
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
        modelNames <- mclust.options("emModelNames")
        if (n <= d) {          
          # select only spherical and diagonal models
          m <- match(modelNames, c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
                     nomatch = 0)
          modelNames <- modelNames[m]
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
      for (mdl in modelNames[BIC["1",] == EMPTY]) {
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
            hcPairs <- hc(modelName = mclust.options("hcModelName"), data = data)
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
            hcPairs <- hc(data = data[initialization$subset,  ],
                          modelName = mclust.options("hcModelName"))
          }
          else {
            hcPairs <- hc(data = data[initialization$subset,],
                          modelName = "EII")
          }
        }
        else {
          hcPairs <- NULL
          #    hcPairs <- hc(data = data[initialization$subset],
          #                  modelName = "E")
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
      hood <- n * log(Vinv)
      BIC["0",  ] <- 2 * hood - log(n)
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
          hcPairs <- hc(data = data[!noise,],
                        modelName = mclust.options("hcModelName"))
        }
        else {
          hcPairs <- hc(data = data[!noise,], modelName = "EII")
        }
      }
      else {
        hcPairs <- NULL 
        #    hcPairs <- hc(data = data[!noise], modelName = "E")
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

# EMclust <- function(...) .Defunct("mclustBIC", PACKAGE = "mclust")
