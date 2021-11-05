##
## Model-based Agglomerative Hierarchical Clustering (MBAHC)
##

# MBAHC used for EM initialization for d-dim data ----

hc <- function(data, 
               modelName = "VVV", 
               use = "VARS",
               partition = dupPartition(data), 
               minclus = 1, ...)
{
  if(!any(modelName == c("E", "V", "EII", "VII", "EEE", "VVV")))
     stop("invalid 'modelName' argument for model-based hierarchical clustering. See help(mclust.options)")

  if(!any(use == c("VARS", "STD", "SPH", "PCS", "PCR", "SVD", "RND")))
     stop("invalid 'use' argument for model-based hierarchical clustering. See help(mclust.options)")

  funcName <- paste("hc", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc$use <- mc$modelName <- NULL
  data <- data.matrix(data)

  dropCols <- function(x)
  { # select only those columns of matrix x with all finite numerical values
    x[,apply(x, 2, function(x) all(is.finite(x))), drop = FALSE]
  }

  use <- toupper(use[1])
  switch(use,
         "VARS" = { Z <- data },
         "STD" = { Z <- scale(data, center = TRUE, scale = TRUE) 
                   Z <- dropCols(Z) },
         "PCR" = { data <- scale(data, center = TRUE, scale = TRUE)
                   data <- dropCols(data)
                   SVD <- svd(data, nu=0)
                   # evalues <- sqrt(SVD$d^2/(nrow(data)-1))
                   Z <- data %*% SVD$v },
         "PCS" = { data <- scale(data, center = TRUE, scale = FALSE)
                   SVD <- svd(data, nu=0)
                   # evalues <- sqrt(SVD$d^2/(nrow(data)-1))
                   Z <- data %*% SVD$v 
                   Z <- dropCols(Z) },
         "SPH" = { data <- scale(data, center = TRUE, scale = FALSE)
                   n <- nrow(data); p <- ncol(data)
                   Sigma <- var(data) * (n - 1)/n
                   SVD <- svd(Sigma, nu = 0)
                   Z <- data %*% SVD$v %*% diag(1/sqrt(SVD$d), p, p) 
                   Z <- dropCols(Z) },
         "SVD" = { data <- scale(data, center = TRUE, scale = TRUE)
                   data <- dropCols(data)
                   p <- min(dim(data))
                   SVD <- svd(data, nu=0)
                   Z <- data %*% SVD$v %*% diag(1/sqrt(SVD$d), p, p) },
         "RND" = { out <- hcRandomPairs(data, ...) 
                   attr(out, "dimensions") <- dim(data)
                   attr(out, "use") <- use
                   attr(out, "call") <- match.call()
                   class(out) <- "hc"
                   return(out) }
        )
  
  # call the proper hc<funcName> function
  mc$data <- Z 
  mc[[1]] <- as.name(funcName)
  out <- eval(mc, parent.frame())
  attr(out, "use") <- use
  attr(out, "call") <- match.call()
  attr(out, "data") <- mc$data
  class(out) <- "hc"
  return(out)
}

print.hc <- function(x, ...) 
{
  if(!is.null(attr(x, "call"))) 
  { 
    cat("Call:\n")
    catwrap(paste0(deparse(attr(x, "call"))))
    cat("\n")
  }
  catwrap("Model-Based Agglomerative Hierarchical Clustering")
  if(!is.null(attr(x, "modelName")))
    cat(paste("Model name        =", attr(x, "modelName"), "\n"))
  if(!is.null(attr(x, "use")))
    cat(paste("Use               =", attr(x, "use"), "\n"))
  if(!is.null(attr(x, "dimensions")))
    cat(paste("Number of objects =", attr(x, "dimensions")[1], "\n"))
  invisible(x)
}

randomPairs <- function(...)
{
  .Deprecated(old = "randomPairs", 
              new = "hcRandomPairs",
              package = "mclust")
  hcRandomPairs(...)
}

hcRandomPairs <- function(data, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  data <- as.matrix(data)
  n <- nrow(data)
  m <- if(n%%2 == 1) n-1 else n
  tree <- matrix(sample(1:n, m, replace = FALSE), 
                 nrow = 2, ncol = ceiling(m/2))
  tree <- apply(tree, 2, sort)
  ind <- unique(tree[1,])
  while(ncol(tree) < (m-1))
  { 
    addtree <- sort(sample(ind, size = 2, replace = FALSE))
    ind <- setdiff(ind, addtree[2])
    tree <- cbind(tree, addtree)
  }
  dimnames(tree) <- NULL
  structure(tree, initialPartition = 1:n, dimensions = c(n,2))
}

dupPartition <- function(data) 
{
  dup <- duplicated(data)
  if (is.null(dim(data))) 
  {
    data <- as.numeric(data)
    if (!any(dup)) return(1:length(data))
    kmeans(data, centers = data[!dup])$cluster
  } else 
  {
    data <- data.matrix(data)
    if (!any(dup)) return(1:nrow(data))
    kmeans(data, centers = data[!dup,])$cluster
  }
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
  if(any(bad) & mclust.options("warn"))
    { warning("Some selected classifications are inconsistent with mclust object") }
  L <- length(select)
  cl <- matrix(as.double(NA), nrow = n, ncol = L, 
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

hcEII <- function(data, partition = NULL, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || NCOL(data) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = TRUE)
  l <- length(unique(partition))
  attr(partition, "unique") <- l
  m <- l - minclus
  if(m <= 0)
    { stop("initial number of clusters is not greater than minclus") }
  if(n <= p & mclust.options("warn"))
    { warning("# of observations <= data dimension") }
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
  structure(t(temp[[1]]), 
            initialPartition = partition, 
            dimensions = dimdat, 
            modelName = "EII", 
            call =  match.call())
}

hcEEE <- function(data, partition = NULL, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || NCOL(data) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p & mclust.options("warn"))
    warning("# of observations <= data dimension")
  if(is.null(partition))
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
  # Luca: commented the next line and uncommented below
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
                   PACKAGE = "mclust")[c(1, 7:10)]
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
  structure(tree,  
            initialPartition = partition, 
            dimensions = dimdat, 
            modelName = "EEE", 
            call = match.call())
}

hcVII <- function(data, partition = NULL, minclus = 1, alpha = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || NCOL(data) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p & mclust.options("warn"))
    warning("# of observations <= data dimension")
  if(is.null(partition))
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
  structure(t(temp[[1]]), 
            initialPartition = partition, 
            dimensions = dimdat, 
            modelName = "VII", 
            call = match.call())
}

hcVVV <- function(data, partition = NULL, minclus = 1, alpha = 1, beta = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || NCOL(data) == 1)
  #if(oneD || length(dimdat) > 2)
  #  stop("data should in the form of a matrix")
  data <- as.matrix(data)
  dimnames(data) <- NULL
  n <- nrow(data)
  p <- ncol(data)
  if(n <= p & mclust.options("warn"))
    warning("# of observations <= data dimension")
  if(is.null(partition))
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
  #x[c((1:n)[!dp],(1:n)[dp]),], 
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
  structure(t(temp[[1]]), 
            initialPartition = partition, 
            dimensions = dimdat, 
            modelName = "VVV", 
            call = match.call())
}

##
## Plot method (dendrogram) for model-based hierarchical agglomeration ----
##
plot.hc <-
function (x, what=c("loglik","merge"), maxG=NULL, labels=FALSE, hang=0,...) 
{
    stopifnot(inherits(x, "hc"))
    what <- what[1]
    hier <- as.hclust(x, what = what, maxG = maxG, labels = labels)
    switch(what,
              "loglik" = {
	     ylab <- paste("Classification log-likelihood",
	                  paste("(", hier$method, sep = ""), "model)")
	     cloglik <- attr(hier,"cloglik")
	     attr(hier,"cloglik") <- NULL
	     plot( as.dendrogram(hier, hang=hang), axes=F, ylab=ylab)
	     r <- range(cloglik,na.rm=T)
	     par.usr <- par("usr")
             ybot <- max(r)-par.usr[3]
             ytop <- min(r)+par.usr[3]
	    },
	    "merge" = {
	     ylab <- paste("Number of Clusters",
	                paste("(", hier$method, sep = ""), "model)")
	     nclus <- attr(hier,"nclus")
	     attr(hier,"nclus") <- NULL
	     plot( as.dendrogram(hier, hang=hang), axes=F, ylab=ylab)
	     par.usr <- par("usr")
	     ybot<- max(nclus)-par.usr[3]
	     ytop <- 1+par.usr[3]
            },
	    stop("unrecognized what option"))

    par(usr=c(par("usr")[1:2],ybot,ytop))
    at <- pretty(seq(from=ybot,to=ytop,length=100), min = 5, max = 10)
    axis(2, at=at)
    
    invisible(hier)
}

as.hclust.hc <-
function (object, what = c("loglik", "merge"), maxG = NULL, labels = FALSE) 
{
    stopifnot(inherits(object, "hc"))

    if (!is.null(maxG) && maxG < 2) stop("maxG < 2")

    what <- what[1]
    switch( what,
            "loglik" = {
	                obj <- ldend(object,maxG=maxG,labels)
     		        obj <- c(obj, list(dist.method = NULL))
			attr(obj,"cloglik") <- as.vector(obj$cloglik)
			obj$cloglik <- NULL
			class(obj) <- "hclust"
			obj
		       },
            "merge" = {
	               obj <- mdend(object,maxG=maxG,labels)
    		       obj <- c(obj, list(dist.method = NULL))
                       attr(obj,"nclus") <- as.vector(obj$nclus)
		       obj$nclus <- NULL
		       class(obj) <- "hclust"
		       obj
		      },
	    stop("unrecognized what option")
	   )
}

ldend <- function (hcObj, maxG = NULL, labels = FALSE) 
{
 stopifnot(inherits(hcObj,"hc"))

# classification log-likelihood dendrogram setup for MBAHC

 if (!is.null(maxG) && maxG < 2) stop("maxG < 2")

 n <- ncol(hcObj) + 1
 cLoglik <- CLL <- cloglik.hc(hcObj)
 
 maxG <- if (is.null(maxG)) length(CLL) else min(maxG,length(CLL))
 
 na <- is.na(CLL)
 m <- length(CLL)
 d <- diff(CLL)
 if (any(neg <- d[!is.na(d)] < 0)) {
   m <- which(neg)[1]
   CLLmax <- CLL[min(maxG,m)]
   CLL[-(1:min(maxG,m))] <- CLLmax
 }
 else if (any(na)) {
   m <- which(na)[1] - 1
   CLLmax <- CLL[min(maxG,m)]
   CLL[-(1:min(maxG,m))] <- CLLmax
 }
 else {
   CLLmax <- max(CLL[1:maxG])
   CLL[-(1:maxG)] <- CLLmax
 }
 
 height <- CLL
 height <- height[-length(height)]
 height <- rev(-height+max(height))

 mo <- mergeOrder(hcObj)

 nam <-  rownames(as.matrix(attr(hcObj,"data"))) 
 leafLabels <- if (labels) nam  else character(length(nam))
 
 obj <- structure(list(merge = mo$merge, 
                       height = height, 
                       order = mo$order, 
                       labels = leafLabels, 
                       cloglik = cLoglik, 
                       method = attr(hcObj, "model"), 
                       call = attr(hcObj, "call")))
 return(obj)
}

mdend <-
function (hcObj, maxG = NULL, labels = FALSE) 
{
 stopifnot(inherits(hcObj,"hc"))

# uniform height dendrgram setup for MBAHC

 if (!is.null(maxG) && maxG < 2) stop("maxG < 2")

 ni <- length(unique(attr(hcObj,"initialPartition")))

 if (!is.null(maxG)) maxG <- min(maxG, ni) else maxG <- ni
 
 mo <- mergeOrder(hcObj)

 j <- ni - maxG
 n <- ncol(hcObj)
 height <- c(rep(0,j),1:(n-j))
 
 nclus <- maxG:1

 nam <- rownames(as.matrix(attr(hcObj,"data"))) 
 leafLabels <- if (labels) nam else character(length(nam))
 
 obj <- structure(list(merge = mo$merge, order = mo$order, height = height,
           labels = leafLabels, nclus = nclus,
	   method = attr(hcObj, "model"), call = attr(hcObj, "call")))
 obj
}

mergeOrder <-
function(hcObj) 
{
# converts the hc representation of merges to conform with hclust
# and computes the corresponding dendrogram leaf order
# CF: inner code written by Luca Scrucca

    HC <- matrix(as.vector(hcObj), ncol(hcObj), nrow(hcObj), byrow = TRUE)
    HCm <- matrix(NA, nrow(HC), ncol(HC))

    merged <- list(as.vector(HC[1, ]))
    HCm[1, ] <- -HC[1, ]
    for (i in 2:nrow(HC)) {
        lmerged <- lapply(merged, function(m) HC[i, ] %in% m)
        lm <- which(sapply(lmerged, function(lm) any(lm)))
        if (length(lm) == 0) {
            merged <- append(merged, list(HC[i, ]))
            HCm[i, ] <- sort(-HC[i, ])
        }
        else if (length(lm) == 1) {
            merged <- append(merged, list(c(merged[[lm]], HC[i, 
                !lmerged[[lm]]])))
            merged[[lm]] <- list()
            HCm[i, ] <- sort(c(-HC[i, !lmerged[[lm]]], lm))
        }
        else {
            merged <- append(merged, list(unlist(merged[lm])))
            merged[[lm[1]]] <- merged[[lm[2]]] <- list()
            HCm[i, ] <- lm
        }
    }

  list(merge = HCm, order = merged[[length(merged)]])
}

cloglik.hc <-
function(hcObj, maxG = NULL) {

n <- ncol(hcObj) + 1

if (is.null(maxG)) maxG <- n

cl <- hclass(hcObj)
cl <- cbind( "1" = 1, cl)

modelName <- attr(hcObj,"modelName")

LL <- rep(list(NA),maxG)
for (j in 1:maxG) {
   ll <- NULL
   for (k in unique(cl[,j])) {
     i <- which(cl[,j] == k)
     # compute loglik term here
     llnew <- mvn( modelName, attr(hcObj,"data")[i,,drop=FALSE])$loglik
     if (substr(modelName,2,2) != "I") {
         llvii <- mvn( "VII", attr(hcObj,"data")[i,,drop=FALSE])$loglik
         if (substr(modelName,3,3) != "I") {
           llvvi <- mvn( "VVI", attr(hcObj,"data")[i,,drop=FALSE])$loglik
  	   llall <- c("VVV"=llnew,"VVI"=llvvi,"VII"=llvii)
	 }
	 else {
  	   llall <- c("VVI"=llnew,"VII"=llvii)
	 }
	 if (!all(nall <- is.na(llall))) {
           llnew <- llall[!nall][which.max(llall[!nall])]
	 }
     }
     if (is.na(llnew)) break
     ll <- c(ll, llnew)
   }
   if (is.na(llnew)) break
   LL[[j]] <- ll
 }
 CLL <- sapply(LL,sum)
 for (i in seq(along = CLL)) {
    if (is.na(CLL[i])) LL[[i]] <- NA
 }
 attr(CLL,"terms") <- LL
 CLL
}

## Initialization for 1-dim data ----

qclass <- function (x, k) 
{
  x <- as.vector(x)
  # eps <- sqrt(.Machine$double.eps) 
  # numerical accuracy problem if scale of x is large, so make tolerance
  # scale dependent
  eps <- sd(x)*sqrt(.Machine$double.eps)
  q <- NA
  n <- k
  while(length(q) < (k+1))
  { n <- n + 1
    q <- unique(quantile(x, seq(from = 0, to = 1, length = n))) 
  }
  if(length(q) > (k+1))
  { dq <- diff(q)
    nr <- length(q)-k-1
    q <- q[-order(dq)[1:nr]]
  }
  q[1] <- min(x) - eps
  q[length(q)] <- max(x) + eps
  cl <- rep(0, length(x))
  for(i in 1:k) 
     { cl[ x >= q[i] & x < q[i+1] ] <- i }
  return(cl)
}

hcE <- function(data, partition = NULL, minclus = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || NCOL(data) == 1)
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(partition))
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
  structure(rbind(temp[[1]], temp[[2]]),
            initialPartition = partition, 
            dimensions = n, 
            modelName = "E",
            call = match.call())
}

hcV <- function(data, partition = NULL, minclus = 1, alpha = 1, ...)
{
  if(minclus < 1) stop("minclus must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  #=====================================================================
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || NCOL(data) == 1)
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(partition))
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
  structure(rbind(temp[[1]], temp[[2]]),
            initialPartition = partition, 
            dimensions = n, 
            modelName = "V",
            call = match.call())
}
