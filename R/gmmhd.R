######################################################
##                                                  ##
##  Identifying Connected Components in Gaussian    ##
##     Finite Mixture Models for Clustering         ##
##                                                  ##
## Author: Luca Scrucca                             ##
######################################################


gmmhd <- function(object, 
                  ngrid = min(round((log(nrow(data)))*10), nrow(data)), 
                  dr = list(d = 3, lambda = 1, cumEvalues = NULL, mindir = 2),
                  classify = list(G = 1:5, 
                                  modelNames = mclust.options("emModelNames")[-c(8,10)]),
                  ...)
{
  if(!inherits(object, "Mclust"))
    stop("first argument must be an object of class 'Mclust'")
  
  if(!requireNamespace("geometry", quietly = TRUE))
    stop("Package 'geometry' is required. Please install it.")
  
  data <- object$data 
  n <- nrow(data)
  if(ngrid > n)
    { warning("ngrid too large, set equal to n")
      n.grid <- n }

  mNames <- attr(object$BIC, "modelNames")
  if(is.null(dr$d)) 
    dr$d <- 2
  if(is.null(dr$lambda)) 
    dr$lambda <- 1
  if(is.null(classify$G))
    classify$G <- 1:5
  if(is.null(classify$modelNames)) 
    classify$modelNames <- mNames
  classify$modelNames <- intersect(classify$modelNames, mNames)
  if(is.null(dr$mindir)) 
    dr$mindir <- 2

  if(ncol(data) >= dr$d)
    { 
      # compute GMMDR directions
      DR <- MclustDR(object, lambda = dr$lambda)
      # subset selection of GMMDR directions      
      evalues <- DR$evalues[seq(DR$numdir)]
      if(is.null(dr$cumEvalues))
        { # if dr$cumEvalues not provided 
          # perform suset selection of GMMDR directions   
          DR <- MclustDRsubsel(DR, 
                               G = attr(object$BIC, "G"), 
                               modelNames = mNames, 
                               mindir = dr$mindir, 
                               verbose = FALSE)
          dims <- seq(DR$numdir)
      }
      else
        { # select the smallest subset with cumsum eigenvalues > dr$cumEvalues
          dims <- min(which(cumsum(evalues/sum(evalues)) > dr$cumEvalues))
          dims <- seq(min(dr$mindir, dims))
      }

      # estimate the density from Mclust model on the selected directions
      x <- DR$dir[,dims,drop=FALSE]
      colnames(x) <- paste("GMMDR dir", 1:ncol(x), sep = "")
      mc <- object$call
      mc$data <- x
      mc$modelNames <- mNames
      mc$verbose <- FALSE
      obj <- eval(mc, parent.frame())
      DR$parameters <- obj$parameters
      fdens <- dens(modelName = obj$modelName,
                    data = x, parameters = obj$parameters)
  }
  else
    { x <- data
      DR <- NULL
      fdens <- dens(modelName = object$modelName,
                    data = x, parameters = object$parameters)
    }
  
  p <- ncol(x)
  xscaled <- scale(x, colMeans(x), apply(x,2,sd))

  # if to add vertices of convex envelope
  # xrange <- apply(x, 2, range)
  # xbound <- do.call("expand.grid", matrix2list(xrange))
  # x <- rbind(as.matrix(x), as.matrix(xbound*1.1))
  # fdens <- c(fdens, rep(0,nrow(xbound)))

  # uniform grid of proportions for which quantiles are calculated
  pn <- seq(0, 1, length = ngrid)
  qn <- as.numeric(quantile(fdens[1:n], 1-pn))
  nc <- pc <- rep(0, length(qn))
  con <- vector("list", length = length(qn))

  # Delaunay triangulation matrix of dim (m x p+1), where each row provides a
  # set of indices to the points describing a simplex of dimension p
  mode(xscaled) <- "double" # delaunayn requires a real matrix
  DT <- suppressMessages(geometry::delaunayn(xscaled, options="QJ"))
  # plot(x); for(l in 1:nrow(DT)) polygon(x[DT[l,],], border = grey(.8))
  on.exit(unlink("qhull_out.txt"))

  # Graph of neighborhood for each point
  NB <- vector(mode = "list", length = n)
  for(i in seq(n)) 
     { NB[[i]] <- sort(unique(as.vector(DT[rowSums(DT==i)>0,]))) }
  for(i in seq(length(qn)))
  { 
    c <- qn[i]
    Sc <- which(fdens[1:n] > c); names(Sc) <- NULL
    if(length(Sc) < 1) next()
    pc[i] <- length(Sc)/n    
    
    # select neighborhoods of edges with density > c level 
    nb <- NB[Sc]
    # select within neighborhoods those edges whose density > c level
    nb <- lapply(nb, function(nb) sort(intersect(nb, Sc)))
    nb <- nb[!duplicated(nb)]
    # table(sapply(nb,length))
    # remove neighborhoods which do not share any facet, i.e. having
    # less than p edges/obs
    # nb <- nb[sapply(nb, length) >= p]
    # remove neighborhoods which are not simplices of dim (p+1)
    nb <- nb[sapply(nb, length) > p]
        
    # get connected components
    ConComp <- ConnectComp(nb)
    # sapply(ConComp,length); ConComp
    if(length(ConComp) < 1) next()
    nc[i] <- length(ConComp)
    con[[i]] <- ConComp # lapply(ConComp, sort)
  }
  #
  obj <- list(Mclust = object,
              MclustDA = NULL,
              MclustDR = DR,
              x = x,  # i.e. the input data or GMMDR directions
              density = fdens[1:n],
              con = con, 
              nc = structure(nc, names = format(pn, digit = 3)),
              pc = pc,
              pn = pn,
              qn = structure(qn, names = format(pn, digit = 3)),
              clusterCores = NULL,
              cluster = NULL,
              numClusters = NULL)
  class(obj) <- "gmmhd"
  # cluster cores
  obj$clusterCores <- gmmhdClusterCores(obj)
  # semi-supervised classification
  modClass <- gmmhdClassify(obj, G = classify$G,
                            modelNames = classify$modelNames,
                            verbose = FALSE)
  obj$MclustDA <- modClass$model
  obj$cluster <- modClass$cluster
  obj$numClusters <- length(tabulate(obj$cluster))

  return(obj)
}

print.gmmhd <- function(x, digits = getOption("digits"), ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  cat(paste0(" Mclust initial model = (", 
             x$Mclust$modelName, ",",
             x$Mclust$G, ")\n"))
  if(!is.null(x$MclustDR))
    cat(paste0(" MclustDR projection = (", 
               x$MclustDR$modelName, ",",
               x$MclustDR$G, ")\n"))
  cat(paste0(" GMMHD final number of clusters = ", x$numClusters, "\n"))
  invisible()
}
  
summary.gmmhd <- function(object, ...)
{
  title <- paste("GMM with high-density connected components for clustering")
  out <- with(object, 
              list(title = title,
                   "Mclust" = list("G" = Mclust$G,
                                   "modelName" = Mclust$modelName),
                   "MclustDR" = list("G" = MclustDR$G,
                                     "modelName" = MclustDR$modelName),
                   "clusterCores" = table(clusterCores, 
                                          useNA = "ifany", dnn = NULL),
                   "cluster"      = table(cluster, 
                                          useNA = "ifany", dnn = NULL)))
  if(is.null(object$MclustDR)) out$MclustDR <- NULL
  class(out) <- "summary.gmmhd"
  return(out)
}

print.summary.gmmhd <- function(x, digits = getOption("digits"), ...)
{
  cat(rep("-", nchar(x$title)),"\n",sep="")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)),"\n",sep="")
  #
  cat("\nInitial model:  Mclust (", 
      x$Mclust$modelName, ",", 
      x$Mclust$G, ")",
      "\n", sep = "")
  #
  if(!is.null(x$MclustDR))
    cat("\nModel on projection subspace:  (", 
        x$MclustDR$modelName, ",", 
        x$MclustDR$G, ")",
        "\n", sep = "")
  #
  cat("\nCluster cores:\n")
  print(x$clusterCores)
  #
  cat("\nFinal clustering:\n")
  print(x$cluster)
  #
  invisible()  
}


plot.gmmhd <- function(x, what = c("mode", "cores", "clusters"), ...)
{
  object <- x
  what <- match.arg(what, choices = eval(formals(plot.gmmhd)$what), 
                    several.ok = TRUE)
  if(interactive() & length(what) > 1)
    { title <- "GMM high-density connected components:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "mode")     plot.gmmhd.mode(object, ...)
             if(what[choice] == "cores")    plot.gmmhd.cores(object, ...)
             if(what[choice] == "clusters") plot.gmmhd.clusters(object, ...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  } 
  else 
    { if(any(what == "mode"))      plot.gmmhd.mode(object, ...)
      if(any(what == "cores"))     plot.gmmhd.cores(object, ...) 
      if(any(what == "clusters"))  plot.gmmhd.clusters(object, ...) 
  }

  invisible()
}

plot.gmmhd.mode <- function(object, ...)
{ 
  plot(c(object$pc,1), c(object$nc,0), type = "S",
       xlab = "Proportion of observed data",
       ylab = "Mode function", yaxt = "n") 
  axis(side = 2, at = seq(0, max(object$nc, na.rm = TRUE)))
}

plot.gmmhd.cores <- function(object, 
                             col = c("grey50", mclust.options("classPlotColors")), 
                             pch = c(1, mclust.options("classPlotSymbols")), ...)
{
  x <- object$x
  p <- ncol(x)
  n <- nrow(x)
  clCores <- object$clusterCores
  numClusters <- object$numClusters
  colCores <- col[1]
  col <- col[-1]
  col <- col[clCores]
  col[is.na(col)] <- colCores
  pch <- unique(pch)
  pchCores <- pch[1]
  pch <- pch[-1]
  pch <- pch[clCores]
  pch[is.na(pch)] <- pchCores
  cex <- rep(par("cex"), length(pch))
  cex[is.na(clCores)] <- par("cex")/2
  
  if(p == 1)
    { plot(x, object$density, col = col, pch = pch, cex = cex,
           ylim = range(0,object$density),
           xlab = colnames(x)[1], ylab = "Density", ...) }
  else if(p == 2)
    { plot(x[,1:2,drop=FALSE], col = col, pch = pch, cex = cex, ...) }
  else if(p > 2)
    { pairs(x, col = col, pch = pch, cex = cex, gap = 0, ...) }
  
  invisible()
}

plot.gmmhd.clusters <- function(object,
                                col = mclust.options("classPlotColors"), 
                                pch = mclust.options("classPlotSymbols"),
                                ...)
{
  x <- object$x
  p <- ncol(x)
  n <- nrow(x)
  cluster <- object$cluster
  numClusters <- object$numClusters
  col <- col[cluster]
  pch <- setdiff(pch,22)[cluster]

  if(p == 1)
    { plot(x, object$density, col = col, pch = pch, 
           ylim = range(0,object$density),
           xlab = colnames(x)[1], ylab = "Density", ...) }
  else if(p == 2)
    { plot(x[,1:2,drop=FALSE], col = col, pch = pch, ...) }
  else if(p > 2)
    { pairs(x, col = col, pch = pch, cex = 0.8, gap = 0, ...) }
  
  invisible()
}

gmmhdClusterCores <- function(object, tails = FALSE, ...) 
{
# Identify cluster cores as the first subset of connected components 
# corresponding to the largest local mode

  n <- nrow(object$x)
  nc <- object$nc
  pc <- object$pc
  conComp <- object$con

  # select the subset with largest number of modes ...
  i <- which(diff(c(nc,0)) < 0)
  # i <- i[which(nc[i] == max(nc[i]))] # no to consider only the highest mode
  # remove spurius local modes, i.e. those not identified by at least 
  # two consecutive density level
  # LS:20150107 okmode <- which(nc[i] == nc[i-1])[1]
  # LS:20150107 i <- if(length(okmode) > 0) i[okmode] else length(nc)
  # plot(pc, nc); abline(v = pc[i])
  # ... and consider multiplicity of modes
  # LS: 20150107
  i <- which(nc == max(nc[i]))
  #
  cc <- conComp[i]
  clusterCores <- matrix(as.double(NA), n, length(i))
  for(j in 1:ncol(clusterCores))
     for(cl in 1:length(cc[[j]]))
        { clusterCores[cc[[j]][[cl]],j] <- cl }
  while(ncol(clusterCores) > 1)
  { 
    ncl <- length(unique(na.omit(clusterCores[,2])))
    tmp <- rep(NA, n)
    for(cl in 1:ncl)
       { 
         l <- which(clusterCores[,2] == cl)
         if(all(is.na(clusterCores[l,1])))
           { tmp[l] <- paste(clusterCores[l,2],"*",sep="") }
         else
           {
             if(length(unique(na.omit(clusterCores[l,1]))) > 1)
                tmp[l] <- clusterCores[l,1]
             else
                tmp[l] <- paste(clusterCores[l,2],"*",sep="")
           }
       }
    clusterCores[,2] <- unclass(as.factor(tmp))
    clusterCores <- clusterCores[,-1,drop=FALSE]
  }
  clusterCores <- as.vector(clusterCores)
  return(clusterCores)
  
  # select the last subset with largest number of modes
  # i <- max(which(nc == max(nc))) 
  
  # select the first subset with largest number of modes
  i <- which(diff(c(nc,0)) < 0)
  i <- i[which(nc[i] == max(nc[i]))[1]]
  
  # select the largest subset with the largest number of modes
  # i <- i[max(which(nc[i] == max(nc[i])))]

  conComp <- object$con[[i]]
  clusterCores <- rep(NA, n)
  for(cl in 1:length(conComp))
     { clusterCores[conComp[[cl]]] <- cl }

  return(clusterCores)
}  

gmmhdClassify <- function(object, G = 1:5, 
                          modelNames = mclust.options("emModelNames"), 
                          verbose = TRUE, ...)
{
  if(!inherits(object, "gmmhd"))
    stop("object is not of class 'gmmhd'")
  x <- object$x
  n <- nrow(x)
  p <- ncol(x)
  if(p == 1) 
    modelNames <- unique(substr(modelNames, 1, 1))

  clusterCores <- object$clusterCores
  numClusters <- length(tabulate(clusterCores))
  con <- object$con
  
  # classify unclustered obs based on training cluster cores
  isCore <- (!is.na(clusterCores))

  logRatio <- function(p) 
  { 
    p <- pmax(pmin(p, 1-sqrt(.Machine$double.eps)),sqrt(.Machine$double.eps))
    log(p)-log(1-p)
  }

  # select num. components G to guarantee at least minSize obs per class
  numCompClass <- function(class, G, minSize = 10)
  { 
    classSize <- tabulate(class)
    Gin <- as.vector(G)
    Gmax <- classSize %/% minSize
    Gmax <- pmin(Gmax, max(G))
    G <- vector(length = length(Gmax), mode = "list")
    for(k in 1:length(G))
       { G[[k]] <- intersect(Gin, seq(Gmax[k])) }
    return(G)
  }
  
  inc <- isCore
  cluster <- clusterCores
  while(sum(inc) < n)
  { 
    mod <- MclustDA(data = x[inc,,drop=FALSE], 
                    class = as.character(cluster[inc]), 
                    G = numCompClass(cluster[inc], G), 
                    modelNames = modelNames,
                    verbose = verbose)
    unallocated <- which(!inc)
    # remove those obs with density ~ 0
    dens <- density.MclustDA(mod, newdata=x[unallocated,,drop=FALSE])
    dens <- pmax(dens, .Machine$double.eps)
    i <- (dens/max(dens) > sqrt(.Machine$double.eps))
    if(sum(i) > 0) unallocated <- unallocated[i]
    # 
    pred <- predict(mod, newdata = x[unallocated,,drop=FALSE])

    # questa versione puo' non allocare obs ai clusterCores piccoli
    # zmax <- apply(pred$z,1,max)
    # zclass <- apply(pred$z,1,which.max)
    # log.ratio <- logRatio(zmax) 
    # alloc <- (log.ratio >= quantile(log.ratio, prob = sum(inc)/n))

    # questa versione cerca di ctr per dim clusters e alloca alla classe 
    # predicted iff logRatio is larger than sqrt(sum(inc)/n) quantile
    z <- pred$z
    zclass <- apply(z,1,which.max)
    alloc <- matrix(NA, nrow(z), ncol(z))
    for(k in seq(ncol(z)))
       { log.ratio <- logRatio(z[,k])
         alloc[,k] <- (log.ratio >= quantile(log.ratio, prob = sqrt(sum(inc)/n))) &
                      (zclass == k)
    }
    alloc <- apply(alloc, 1, any)
    
    toclass <- unallocated[alloc]
    cluster[toclass] <- zclass[alloc] 
    inc <- (!is.na(cluster))

  }

  mod <- MclustDA(data = x, class = cluster, 
                  G = numCompClass(cluster[inc], G), 
                  modelNames = modelNames,
                  verbose = verbose)
  cluster <- predict(mod, x)$classification

  out <- list(model = mod, 
              clusterCores = clusterCores,
              cluster = cluster)
  return(out)
}

density.MclustDA <- function(object, newdata, prior, logarithm = FALSE, ...)
{
# Compute the density based on a MclustDA model 
# (later it may be included in the 'mclust' package)
# or it can be obtained from predict.MclustDA
  if(!inherits(object, "MclustDA")) 
    stop("object not of class \"MclustDA\"")
  
  models <- object$models
  nclass <- length(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  if(missing(newdata))
    { newdata <- object$data }
  if(object$d == 1) newdata <- as.vector(newdata)
  if(missing(prior))
    { prior <- n/sum(n) }
  else
    { if(length(prior) != nclass)
        stop("wrong number of prior probabilities")
      if(any(prior < 0))
        stop("prior must be nonnegative")
    }
  
  # compute on log scale for stability
  densfun <- function(mod, data)
  { do.call("dens", c(list(data = data, logarithm = TRUE), mod)) }
  #  
  cden <- as.matrix(data.frame(lapply(models, densfun, data = newdata)))
  cden <- sweep(cden, 2, FUN = "+", STATS = log(prior))
  maxlog <- apply(cden, 1, max)
  cden <- sweep(cden, 1, FUN = "-", STATS = maxlog)
  den <- log(apply(exp(cden), 1, sum)) + maxlog
  if(!logarithm) den <- exp(den)
  return(den) 
}


# old version
ConnectComp_old <- function(nb)
{
# Get connected components
# Example:
# nb <- list(c(1,2,3), c(2,3,4), c(9,10,11), c(9,11,12), c(1,6,5))
# 
  if(length(nb) < 1 | !is.list(nb))
    return(NULL)
  nb <- lapply(nb, function(x) as.integer(x))
  n <- length(nb)
  u <- sort(unique(unlist(nb)))
  nu <- length(u)
  cnb <- cnb.old <- nb
  stable <- FALSE
  # merge the neighbors until the configuration is stable
  while(!stable)
  { i <- 0
    while(i < length(cnb))
         { i <- i + 1
           j <- which(sapply(cnb, function(nbb) any(intersect(cnb[[i]], nbb))))
           cnb[[i]] <- sort(unique(unlist(cnb[j])))
           cnb[setdiff(j, i)] <- NULL
         }
    if(identical(cnb, cnb.old)) stable <- TRUE
    cnb.old <- cnb
  }
  return(cnb)
}

ConnectComp <- function(nb)
{
# Get connected components
# Example:
# nb <- list(c(1,2,3), c(2,3,4), c(9,10,11), c(9,11,12), c(1,6,5))
# ConnectComp(nb)
  
  if(length(nb) < 1 | !is.list(nb))
    return(NULL)
  nb <- lapply(nb, function(x) as.integer(x))
  n <- length(nb)
  u <- sort(unique(unlist(nb)))
  nu <- length(u)
  cnb <- cnb.old <- nb
  stable <- FALSE
  # merge the neighbors until the configuration is stable
  while(!stable)
  { i <- 0
    while(i < length(cnb))
    { i <- i + 1
      j <- which(sapply(cnb, function(nbb) any(is.element(cnb[[i]], nbb))))
      cnb[[i]] <- sort(unique(unlist(cnb[j])))
      cnb[setdiff(j, i)] <- NULL
    }
    if(identical(cnb, cnb.old)) stable <- TRUE
    cnb.old <- cnb
  }
  return(cnb)
}
