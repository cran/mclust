##
## Model-based Agglomerative Hierarchical Clustering (MBAHC)
##

# MBAHC used for EM initialization for d-dim data ----

hc <- function(data, 
               modelName = mclust.options("hcModelName"), 
               use = mclust.options("hcUse"), ...)
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
  { # select only those columns of matrix x with all finite numeric values
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
         "RND" = { out <- randomPairs(data, ...) 
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

randomPairs <- function(data, seed, ...)
{
  if(!missing(seed)) set.seed(seed)
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

hcEII <- function(data, partition, minclus = 1, ...)
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
  if(missing(partition))
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
  structure(t(temp[[1]]), initialPartition = partition, 
            dimensions = dimdat, modelName = "EII", 
            call =  match.call())
}

hcEEE <- function(data, partition, minclus = 1, ...)
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
  structure(tree,  initialPartition = partition, 
            dimensions = dimdat, modelName = "EEE", 
            call = match.call())
}

hcVII <- function(data, partition, minclus = 1, alpha = 1, ...)
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
  structure(t(temp[[1]]), initialPartition = partition, 
            dimensions = dimdat, modelName = "VII", 
            call = match.call())
}

hcVVV <- function(data, partition, minclus = 1, alpha = 1, beta = 1, ...)
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
  structure(t(temp[[1]]), initialPartition = partition, 
            dimensions = dimdat, modelName = "VVV", 
            call = match.call())
}

##
## Dendrogram for model-based hierarchical agglomeration ----
##

as.hclust.hc <- function(x, ...)
{
# Convert 'hc' objects to class 'hclust'
  stopifnot(inherits(x, "hc"))
  data <- as.matrix(attr(x, "data"))
  labels <- rownames(data)
  # convert a 'hc' hierarchical clustering structure to 'hclust' structure
  HC <- matrix(as.vector(x), ncol(x), nrow(x), byrow = TRUE)
  HCm <- matrix(NA, nrow(HC), ncol(HC))
  merged <- list(as.vector(HC[1,]))
  HCm[1,] <- -HC[1,]
  for(i in 2:nrow(HC))
     { lmerged <- lapply(merged, function(m) HC[i,] %in% m)
       lm <- which(sapply(lmerged, function(lm) any(lm)))
       if(length(lm) == 0)
         { merged <- append(merged, list(HC[i,]))
           HCm[i,] <- sort(-HC[i,]) }
       else if(length(lm) == 1)
              { merged <- append(merged, list(c(merged[[lm]], HC[i,!lmerged[[lm]]])))
                merged[[lm]] <- list()
                HCm[i,] <- sort(c(-HC[i,!lmerged[[lm]]], lm)) }
       else   { merged <- append(merged, list(unlist(merged[lm])))
                merged[[lm[1]]] <- merged[[lm[2]]] <- list()
                HCm[i,] <- lm }
  }
  # compute heights
  height <- attr(x, "deviance")
  if(is.null(height)) 
    height <- hcCriterion(x, ...)
  # create 'hclust' object
  obj <- structure(list(merge = HCm, 
                        height = rev(height), 
                        order = merged[[length(merged)]],
                        labels = labels, 
                        method = attr(x, "model"), 
                        dist.method = NULL,
                        call = attr(x, "call")),
                   class = "hclust")
  return(obj)
}  
  

as.dendrogram.hc <- function(object, ...)
{
# Convert 'hc' objects to class 'dendrogram'
  stopifnot(inherits(object, "hc"))
  as.dendrogram(as.hclust(object))
}  
  

plot.hc <- function(x, ...)
{
  stopifnot(inherits(x, "hc"))
  # dots <- list(...)
  # if(is.null(dots$hang)) dots$hang <- -1
  # if(is.null(dots$sub))  dots$sub <- NA
  dendro <- as.dendrogram(x)
  # do.call("plot", c(list(hcl), dots))
  plot(dendro)
  invisible(dendro)
}

# Auxiliary functions ----

hcCriterion <- function(hcPairs, Gmax, what = c("deviance", "loglik"), ...)
{
  stopifnot(inherits(hcPairs, "hc"))
  hcPairsName <- deparse(substitute(hcPairs))
  what <- match.arg(what, choices = eval(formals(hcCriterion)$what))
  data <- as.matrix(attr(hcPairs, "data"))
  N <- nrow(data)
  p <- ncol(data)
  model <- attr(hcPairs, "model")
  m <- ifelse(missing(Gmax), ncol(hcPairs), as.integer(Gmax))
  hc <- hclass(hcPairs, seq_len(m))
  Wdata <- var(data)*(N-1)
  trWnp <- tr(Wdata)/(N*p)
  # detS <- det(Wdata/N)
  loglik <- rep(as.double(NA), length = m)
  # loglik[1] <- mvn(model, data)$loglik

  switch(model, 
         "EII" =
         { for(k in 1:m)
              { n <- tabulate(hc[,k], k)
                # mu <- by(data, as.factor(hc[,k]), FUN = colMeans, 
                #          simplify = FALSE)
                W <- WSS(data, hc[,k])
                sigmasq <- sum(apply(W, 3, tr), na.rm=TRUE)/(N*p)
                loglik[k] <- -0.5*p*N*log(2*pi) -0.5*N*p +
                             -0.5*sum(n*log(sigmasq^p + apply(W, 3, tr)/n + trWnp))
           }
          },
         "VII" =
         { for(k in 1:m)
              { n <- tabulate(hc[,k], k)
                W <- WSS(data, hc[,k])
                sigmasq <- apply(W, 3, tr)/(n*p)
                loglik[k] <- -0.5*p*N*log(2*pi) -0.5*N*p +
                             -0.5*sum(n*log(sigmasq^p + apply(W, 3, tr)/n + trWnp))
           }
         },
         "EEE" = 
         { for(k in 1:m)
              { n <- tabulate(hc[,k], k)
                W <- WSS(data, hc[,k])
                Sigma <- apply(W, 1:2, sum)/N
                loglik[k] <- -0.5*p*N*log(2*pi) -0.5*N*p +
                             -0.5*sum(n*log(det(Sigma) + apply(W, 3, tr)/n + trWnp))
           }
         },
         "VVV" = 
         { for(k in 1:m)
              { n <- tabulate(hc[,k], k)
                W <- WSS(data, hc[,k])
                Sigma <- sapply(1:k, function(k) W[,,k]/n[k], simplify = "array")
                loglik[k] <- -0.5*p*N*log(2*pi) -0.5*N*p +
                             -0.5*sum(n*log(apply(Sigma, 3, det) + 
                                            apply(W, 3, tr)/n + trWnp))
           }
          }
  )              
  deviance <- -2*(loglik - max(loglik, na.rm = TRUE))
  # attr(hcPairs, "loglik") <- loglik
  # attr(hcPairs, "deviance") <- deviance
  # assign(hcPairsName, hcPairs, envir = parent.frame())
  out <- switch(what, 
                "deviance" = deviance,
                "loglik" = loglik,
                NULL)
  return(out)
}

WSS <- function(X, group, ...)
{
  X <- as.matrix(X)
  Z <- unmap(as.vector(group))
  n <- nrow(X)
  p <- ncol(X)
  G <- ncol(Z)
  tmp <- .Fortran("covwf",
                  X = as.double(X),
                  Z = as.double(Z),
                  n = as.integer(n),
                  p = as.integer(p),
                  G = as.integer(G),
                  mean = double(p * G), 
                  S = double(p * p * G), 
                  W = double(p * p * G) )
  array(tmp$W, c(p,p,G))
}

tr <- function(x) 
{
  sum(diag(as.matrix(x)))
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

hcE <- function(data, partition, minclus = 1, ...)
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
  structure(rbind(temp[[1]], temp[[2]]),   initialPartition = partition, 
            dimensions = n, modelName = "E",
            call = match.call())
}

hcV <- function(data, partition, minclus = 1, alpha = 1, ...)
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
  structure(rbind(temp[[1]], temp[[2]]),   initialPartition = partition, 
            dimensions = n, modelName = "V",
            call = match.call())
}
