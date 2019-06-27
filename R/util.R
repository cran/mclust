
adjustedRandIndex <- function (x, y) 
{
  x <- as.vector(x)
  y <- as.vector(y)
  if(length(x) != length(y)) 
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}


classError <- function(classification, class)
{
  q <- function(map, len, x)
  {
    x <- as.character(x)
    map <- lapply(map, as.character)
    y <- sapply(map, function(x) x[1])
    best <- y != x
    if(all(len) == 1)
      return(best)
    errmin <- sum(as.numeric(best))
    z <- sapply(map, function(x) x[length(x)])
    mask <- len != 1
    counter <- rep(0, length(len))
    k <- sum(as.numeric(mask))
    j <- 0
    while(y != z) 
    {
      i <- k - j
      m <- mask[i]
      counter[m] <- (counter[m] %% len[m]) + 1
      y[x == names(map)[m]] <- map[[m]][counter[m]]
      temp <- y != x
      err <- sum(as.numeric(temp))
      if(err < errmin) 
      {
        errmin <- err
        best <- temp
      }
      j <- (j + 1) %% k
    }
    best
  }
  if (any(isNA <- is.na(classification))) 
  {
    classification <- as.character(classification)
    nachar <- paste(unique(classification[!isNA]),collapse="")
    classification[isNA] <- nachar
  }
  MAP <- mapClass(classification, class)
  len <- sapply(MAP[[1]], length)
  if(all(len) == 1) 
  {
    CtoT <- unlist(MAP[[1]])
    I <- match(as.character(classification), names(CtoT), nomatch= 0)
    one <- CtoT[I] != class
  } else 
  {
    one <- q(MAP[[1]], len, class)
  }
  len <- sapply(MAP[[2]], length)
  if(all(len) == 1) 
  {
    TtoC <- unlist(MAP[[2]])
    I <- match(as.character(class), names(TtoC), nomatch = 0)
    two <- TtoC[I] != classification
  } else 
  {
    two <- q(MAP[[2]], len, classification)
  }
  err <- if(sum(as.numeric(one)) > sum(as.numeric(two))) 
            as.vector(one) else as.vector(two)
  bad <- seq(along = classification)[err]
  list(misclassified = bad, errorRate = length(bad)/length(class))
}

mapClass <- function(a, b)
{
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
  #
  k <- nrow(Tab)
  Map <- rep(0, k)
  Max <- apply(Tab, 1, max)
  for(i in 1:k) 
  {
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
  #
  return(list(aTOb = aTOb, bTOa = bTOa))
}

map <- function(z, warn = mclust.options("warn"), ...)
{
  nrowz <- nrow(z)
  cl <- numeric(nrowz)
  I <- 1:nrowz
  J <- 1:ncol(z)
  for(i in I) 
     { cl[i] <- (J[z[i,  ] == max(z[i,  ])])[1] }
  if(warn) 
    { K <- as.logical(match(J, sort(unique(cl)), nomatch = 0))
      if(any(!K))
        warning(paste("no assignment to", paste(J[!K], 
                      collapse = ",")))
  }
  return(cl)
}

unmap <- function(classification, groups=NULL, noise=NULL, ...)
{
  # converts a classification to conditional probabilities
  # classes are arranged in sorted order unless groups is specified
  # if a noise indicator is specified, that column is placed last
  n <- length(classification)
  u <- sort(unique(classification))
  if(is.null(groups))
    { groups <- u }
  else 
    { if(any(match( u, groups, nomatch = 0) == 0)) 
      stop("groups incompatible with classification")
      miss <- match( groups, u, nomatch = 0) == 0
    }
  cgroups <- as.character(groups)
  if(!is.null(noise)) 
    { noiz <- match( noise, groups, nomatch = 0)
      if(any(noiz == 0)) stop("noise incompatible with classification")
      groups <- c(groups[groups != noise],groups[groups==noise])
      noise <- as.numeric(factor(as.character(noise), levels = unique(groups)))
  }
  groups <- as.numeric(factor(cgroups, levels = unique(cgroups)))
  classification <- as.numeric(factor(as.character(classification), levels = unique(cgroups)))
  k <- length(groups) - length(noise)
  nam <- levels(groups)
  if(!is.null(noise)) 
    { k <- k + 1
      nam <- nam[1:k]
      nam[k] <- "noise"
  }
  z <- matrix(0, n, k, dimnames = c(names(classification),nam))
  for(j in 1:k)
     { z[classification == groups[j], j] <- 1 }
  return(z)
}

BrierScore <- function(z, class)
{
  z <- as.matrix(z)
  z <- sweep(z, 1, STATS = rowSums(z), FUN = "/")
  cl <- unmap(class, groups = if(is.factor(class)) levels(class) else NULL)
  if(any(dim(cl) != dim(z)))
    stop("input arguments do not match!")
  sum((cl-z)^2)/(2*nrow(cl))
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

randomOrthogonalMatrix <- function(n, d)
{
# Generate a random orthogonal basis matrix of dimension (n x d) using 
# the method in
# Heiberger R. (1978) Generation of random orthogonal matrices. JRSS C, 27,
#   199-206.
  Q <- qr.Q(qr(matrix(rnorm(n*d), nrow = n, ncol = d)))
  return(Q)
}

logsumexp <- function(x)
{ 
# Numerically efficient implementation of log(sum(exp(x)))
  max <- max(x)
  max + log(sum(exp(x-max)))
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

dmvnorm <- function(data, mean, sigma, log = FALSE) 
{
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  if(missing(mean))
    mean <- rep(0, length = d)
  mean <- as.vector(mean)
  if(length(mean) != d)
    stop("data and mean have non-conforming size")
  if(missing(sigma)) 
    sigma <- diag(d)
  sigma <- as.matrix(sigma)
  if(ncol(sigma) != d) 
    stop("data and sigma have non-conforming size")
  if(max(abs(sigma - t(sigma))) > .Machine$double.eps) 
    stop("sigma must be a symmetric matrix")
  
  # - 1st approach
  # cholsigma <- chol(sigma)
  # logdet <- 2 * sum(log(diag(cholsigma)))
  # md <- mahalanobis(data, center = mean, 
  #                   cov = chol2inv(cholsigma), inverted = TRUE)
  # logdens <- -(ncol(data) * log(2 * pi) + logdet + md)/2
  #
  # - 2nd approach
  # cholsigma <- chol(sigma)
  # logdet <- 2 * sum(log(diag(cholsigma)))
  # mean <- outer(rep(1, nrow(data)), as.vector(matrix(mean,d)))
  # data  <- t(data - mean)
  # conc <- chol2inv(cholsigma)
  # Q <- colSums((conc %*% data)* data)
  # logdens <- as.vector(Q + d*log(2*pi) + logdet)/(-2)
  #
  # - 3rd approach (via Fortran code)
  logdens <- .Fortran("dmvnorm", 
                      as.double(data),  # x
                      as.double(mean),  # mu
                      as.double(sigma), # Sigma
                      as.integer(n),    # n
                      as.integer(d),    # p
                      double(d),        # w
                      double(1),        # hood
                      double(n),        # logdens
                      PACKAGE = "mclust")[[8]]
  #  
  if(log) logdens else exp(logdens)
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
           PACKAGE = "mclust")[[3]]
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
  .Fortran("uncholf",
           as.logical(upper),
           x,
           as.integer(nrow(x)),
           as.integer(ncol(x)),
           integer(1),
           PACKAGE = "mclust")[[2]]
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

errorBars <- function(x, upper, lower, width = 0.1, code = 3, angle = 90, horizontal = FALSE, ...) 
{ 
# Draw error bars at x from upper to lower. If horizontal = FALSE (default)
# bars are drawn vertically, otherwise horizontally.
  if(horizontal)
    arrows(upper, x, lower, x, length = width, angle = angle, code = code, ...)
  else  
    arrows(x, upper, x, lower, length = width, angle = angle, code = code, ...)
}

covw <- function(X, Z, normalize = TRUE)
# Given data matrix X(n x p) and weight matrix Z(n x G) computes
# weighted means(p x G), weighted covariance matrices S(p x p x G) and
# weighted scattering matrices W(p x p x G)
{
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  n <- nrow(X)
  p <- ncol(X)
  nZ <- nrow(Z)
  G <- ncol(Z)
  if(n != nZ) 
    stop("X and Z must have same number of rows")
  if(normalize)
    Z <- t( apply(Z, 1, function(z) z/sum(z)) )
  
  tmp <- .Fortran("covwf",
                  X = as.double(X),
                  Z = as.double(Z),
                  n = as.integer(n),
                  p = as.integer(p),
                  G = as.integer(G),
                  mean = double(p*G),
                  S = double(p*p*G),
                  W = double(p*p*G),
                  PACKAGE = "mclust")
  
  out <- list(mean = matrix(tmp$mean, p,G), 
              S = array(tmp$S, c(p,p,G)),
              W = array(tmp$W, c(p,p,G)) )
  return(out)
}

hdrlevels <- function(density, prob)
{
  if(missing(density) | missing(prob))
    stop("Please provide both 'density' and 'prob' arguments to function call!")
  density <- as.vector(density)
  prob <- pmin(pmax(as.numeric(prob), 0), 1)
  alpha <- 1-prob
  lev <- quantile(density, alpha)
  names(lev) <- paste0(round(prob*100),"%")
  return(lev)
}

catwrap <- function(x, width = getOption("width"), ...)
{
# version of cat with wrapping at specified width
  cat(paste(strwrap(x, width = width, ...), collapse = "\n"), "\n")
}


##
## Convert to a from classes 'Mclust' and 'densityMclust' 
##

as.Mclust <- function(x, ...)
{ 
  UseMethod("as.Mclust")
}

as.Mclust.default <- function(x, ...)
{ 
  if(inherits(x, "Mclust")) x
  else stop("argument 'x' cannot be coerced to class 'Mclust'")
}

as.Mclust.densityMclust <- function(x, ...)
{ 
  class(x) <- c("Mclust", class(x)[1])
  return(x)
}

as.densityMclust <- function(x, ...)
{ 
  UseMethod("as.densityMclust")
}

as.densityMclust.default <- function(x, ...)
{ 
  if(inherits(x, "densityMclust")) x
  else stop("argument 'x' cannot be coerced to class 'densityMclust'")
}

as.densityMclust.Mclust <- function(x, ...)
{ 
  class(x) <- c("densityMclust", class(x))
  x$density <- dens(modelName = x$modelName, data = x$data, 
                    parameters = x$parameters, logarithm = FALSE)
  return(x)
}