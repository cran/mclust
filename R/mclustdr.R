######################################################
##                                                  ##
##        Dimension reduction for model-based       ##
##          clustering and classification           ##
##                                                  ##
## Author: Luca Scrucca                             ##
######################################################


# GMMDR dimension reduction -----------------------------------------------

MclustDR <- function(object, lambda = 0.5,
                     normalized = TRUE, Sigma, 
                     tol = sqrt(.Machine$double.eps))
{
#  Dimension reduction for model-based clustering and classification

  call <- match.call()
  if(!any(class(object) %in% c("Mclust", "MclustDA")))
    stop("object must be of class 'Mclust' or 'MclustDA'")
  
  x <- data.matrix(object$data)
  n <- nrow(x)
  p <- ncol(x)
  lambda <- pmax(0, min(lambda, 1))
  #-----------------------------------------------------------------  
  # overall parameters
  mu <- colMeans(x)
  if(missing(Sigma)) Sigma <- var(x)*(n-1)/n
  # within-cluster parameters based on fitted mixture model
  if(inherits(object, "Mclust"))
  {
    type <- "Mclust"
    G <- object$G
    modelName <- object$modelName
    y <- object$classification
    cl2mc <- seq(G)
    class <- as.factor(y)
    par <- object$parameters
    f <- par$pro
    if(is.null(f)) f <- 1
    if(!is.na(object$hypvol)) f <- f[-length(f)]
    # within-group means
    mu.G <- matrix(par$mean,p,G) 
    # within-group covars
    if(p == 1) 
    { 
      Sigma.G <- array(par$variance$sigmasq, c(p,p,G)) 
    } else     
    { 
      Sigma.G <- par$variance$sigma 
    }
  }
  else if(inherits(object, "MclustDA"))
  { 
    type <- object$type
    modelName <- sapply(object$models, function(m) m$modelName)
    class <- object$class
    class <- factor(class, levels = names(object$models))
    y <- rep(NA, length(class))
    for(i in 1:nlevels(class))
    { 
      y[class == levels(class)[i]] <- paste(levels(class)[i], 
                                            object$models[[i]]$classification, sep =":") 
    }
    y <- as.numeric(factor(y))
    cl2mc <- rep(seq(length(object$models)), 
                 sapply(object$models, function(m) m$G))
    
    m <- sapply(object$models, function(mod) mod$n) 
    ncomp <- sapply(object$models, function(mod) mod$G) 
    G <- sum(ncomp)
    f <- vector(length = G)
    mu.G <- matrix(as.double(NA), nrow = p, ncol = G)
    Sigma.G <- array(NA, dim = c(p,p,G))
    for(i in 1:length(object$models))
    {
      ii <- seq(c(0,cumsum(ncomp))[i]+1,c(0,cumsum(ncomp))[i+1])
      par <- object$models[[i]]$parameters
      if(is.null(par$pro)) par$pro <- 1 
      f[ii] <- par$pro * m[i]/sum(m)
      # within-group means
      mu.G[,ii] <- par$mean
      # within-group covars
      if(p == 1) 
      { Sigma.G[,,ii] <- array(par$variance$sigmasq, c(p,p,1)) }
      else     
      { Sigma.G[,,ii] <- par$variance$sigma }
    }
  }
  #-----------------------------------------------------------------
  SVD <- svd(Sigma, nu = 0, nv = min(n,p))
  pos <- which(SVD$d > max(tol*SVD$d[1], 0)) 
  SVD$d <- SVD$d[pos]
  SVD$v <- SVD$v[,pos,drop=FALSE]
  inv.Sigma <- SVD$v %*% (1/SVD$d * t(SVD$v))
  inv.sqrt.Sigma <- SVD$v %*% (1/sqrt(SVD$d) * t(SVD$v)) 
  #-----------------------------------------------------------------
  # pooled within-group covariance
  S <- matrix(0, p, p)
  for(j in seq_len(G))
    S <- S + f[j]*Sigma.G[,,j]
  #-----------------------------------------------------------------
  # kernel matrix
  M.I <- crossprod(t(sweep(mu.G, 1, FUN="-", STATS=mu))*sqrt(f))
  M.II <- matrix(0, p, p)
  if(lambda < 1) 
  {
    for(j in seq_len(G))
      M.II <- M.II + f[j]*crossprod(inv.sqrt.Sigma%*%(Sigma.G[,,j]-S))
  }
  # convex combination of M_I and M_II
  M <- 2*lambda*crossprod(inv.sqrt.Sigma %*% M.I) + 2*(1-lambda)*M.II
  # regularize the M_II
  # M <- M.I + lambda*M.II
  # M <- crossprod(inv.sqrt.Sigma %*% M.I) + 
  #      (1-lambda)*M.II + lambda/p * diag(p)  
  #
  SVD <- eigen.decomp(M, inv.sqrt.Sigma, invsqrt = TRUE)
  l <- SVD$l; l <- (l+abs(l))/2
  numdir <- min(p, sum(l > sqrt(.Machine$double.eps)))
  basis <- as.matrix(SVD$v)[,seq(numdir),drop=FALSE]
  sdx <- diag(Sigma)
  std.basis <- as.matrix(apply(basis, 2, function(x) x*sdx))  
  if(normalized)
  { basis <- as.matrix(apply(basis, 2, normalize)) 
    std.basis <- as.matrix(apply(std.basis, 2, normalize))
  }
  dimnames(basis) <- list(colnames(x), paste("Dir", 1:ncol(basis), sep=""))
  dimnames(std.basis) <- dimnames(basis)
  Z <- scale(x, scale = FALSE) %*% basis
  #
  out = list(call = call, type = type,
             x = x, Sigma = Sigma, 
             classification = class, mixcomp = y, 
             class2mixcomp = cl2mc,
             G = G, modelName = modelName,
             mu = mu.G, sigma = Sigma.G, pro = f,
             M = M, M.I = M.I, M.II = M.II, 
             lambda = lambda, 
             evalues = l, raw.evectors = as.matrix(SVD$v),             
             basis = basis, std.basis = std.basis,
             numdir = numdir, dir = Z)
  class(out) = "MclustDR"
  return(out)
}

print.MclustDR <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste0("\'", class(x)[1], "\' model object: ")
  catwrap(txt)
  cat("\n")
  catwrap("\nAvailable components:\n")
  print(names(x))
  # str(x, max.level = 1, give.attr = FALSE, strict.width = "wrap")
  invisible(x)
}
  
summary.MclustDR <- function(object, numdir, std = FALSE, ...)
{
  if(missing(numdir)) numdir <- object$numdir
  dim <- seq(numdir)
  
  if(object$type == "Mclust")
  { n <- as.vector(table(object$classification))
    G <- object$G }
  else
  { n <- as.vector(table(object$classification))
    G <- as.vector(table(object$class2mixcomp)) }
  
  obj <- list(type = object$type, 
              modelName = object$modelName, 
              classes = levels(object$classification),
              n = n, G = G, 
              basis = object$basis[,seq(dim),drop=FALSE],
              std = std, std.basis = object$std.basis[,seq(dim),drop=FALSE],
              evalues = object$evalues[seq(dim)],
              evalues.cumperc = with(object, 
{ evalues <- evalues[seq(numdir)]
  cumsum(evalues)/sum(evalues)*100 })
  )
class(obj) <- "summary.MclustDR"
return(obj)
}

print.summary.MclustDR <- function(x, digits = max(5, getOption("digits") - 3), ...)
{
  title <- paste("Dimension reduction for model-based clustering and classification")
  txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")
  catwrap(txt)
  catwrap(title)
  catwrap(txt)

  if(x$type == "Mclust")
  { 
    tab <- data.frame(n = x$n)
    rownames(tab) <- x$classes
    tab <- as.matrix(tab)
    names(dimnames(tab)) <- c("Clusters", "")
    cat("\n")
    catwrap(paste0("Mixture model type: ", x$type, 
                   " (", x$modelName, ", ", x$G, ")"))
    print(tab, quote = FALSE, right = TRUE)
  }
  else if(x$type == "MclustDA" | x$type == "EDDA")
  {
    tab <- data.frame(n = x$n, Model = x$modelName, G = x$G)
    rownames(tab) <- x$classes
    tab <- as.matrix(tab)
    names(dimnames(tab)) <- c("Classes", "")
    cat("\n")
    catwrap(paste("Mixture model type:", x$type))
    print(tab, quote = FALSE, right = TRUE)
  }
  else stop("invalid model type")
  
  cat("\n")
  if(x$std) 
  { 
    catwrap("Standardized basis vectors using predictors scaled to have std.dev. equal to one:")
    print(x$std.basis, digits = digits)
  }
  else 
  { 
    catwrap("Estimated basis vectors:")
    print(x$basis, digits = digits)
  }
  cat("\n")
  
  evalues <- rbind("Eigenvalues" = x$evalues, 
                   "Cum. %" = x$evalues.cumperc)
  colnames(evalues) <- colnames(x$basis)
  print(evalues, digits=digits)
  
  invisible()
}


projpar.MclustDR <- function(object, dim, center = TRUE, raw = FALSE)
{
# Transform estimated parameters to projection subspace given by 
# 'dim' directions 
  x <- object$x
  p <- ncol(x)
  n <- nrow(x)
  G <- object$G
  numdir <- object$numdir
  if(missing(dim)) dim <- seq(numdir)
  numdir <- length(dim)
  if(raw) V <- object$raw.evectors[,dim,drop=FALSE]
  else    V <- object$basis[,dim,drop=FALSE]
  #
  mu <- t(object$mu)
  if(center) mu <- scale(mu, center = apply(x,2,mean), scale = FALSE)
  Mu <- mu %*% V
  #
  sigma <- object$sigma 
  cho <- array(apply(sigma, 3, chol), c(p, p, G))  
  Sigma <- array(apply(cho, 3, function(R) 
    crossprod(R %*% V)), c(numdir, numdir, G))
  #
  return(list(mean = Mu, variance = Sigma))
}  

predict.MclustDR <- function(object, dim = 1:object$numdir, newdata, eval.points, ...)
{  
  dim <- dim[dim <= object$numdir]
  if(missing(newdata) & missing(eval.points))
    { dir <- object$dir[,dim,drop=FALSE] 
  } else if(!missing(newdata))
    { newdata <- as.matrix(newdata) 
      newdata <- scale(newdata, center = colMeans(object$x), scale = FALSE)
      dir <- newdata %*% object$basis[,dim,drop=FALSE]
  } else if(!missing(eval.points))
    { dir <- as.matrix(eval.points) }
  
  n <- nrow(dir)
  G <- object$G # num. components
  nclass <- nlevels(object$classification) # num. classes
  par <- projpar.MclustDR(object, dim)
  Mu <- par$mean
  Sigma <- par$variance
  
  # old version
  # cden <- array(NA, c(n, G))
  # for(j in 1:G)
  #    { cden[,j] <- mvdnorm(dir, Mu[j,], Sigma[,,j], log = FALSE) }
  # z <- sweep(cden, 2, FUN = "*", STATS = object$pro)  
  # den <- apply(z, 1, sum)
  # z <- sweep(z, 1, FUN = "/", STATS = den)
  # new version: more efficient and accurate
  z <- array(NA, c(n, G))
  for(j in 1:G)
     { z[,j] <- mvdnorm(dir, Mu[j,], Sigma[,,j], log = TRUE) }
  z <- sweep(z, 2, FUN = "+", STATS = log(object$pro))
  logden <- apply(z, 1, logsumexp)
  z <- sweep(z, 1, FUN = "-", STATS = logden)
  z <- exp(z)
  #
  zz <- matrix(0, n, nclass)
  for(j in seq(nclass))
     { zz[,j] <- rowSums(z[,object$class2mixcomp == j,drop=FALSE]) }
  z <- zz; rm(zz)

  class <- factor(apply(z,1,which.max), 
                  levels = 1:nclass, 
                  labels = levels(object$classification))
  
  out <- list(dir = dir,
              density = exp(logden),
              z = z,
              uncertainty = 1 - apply(z,1,max),
              classification = class)
  return(out)
}

predict2D.MclustDR <- function(object, dim = 1:2, ngrid = 100, xlim, ylim)
{
  dim <- dim[1:2]
  dir <- object$dir[,dim,drop=FALSE]
  G <- object$G
  par <- projpar.MclustDR(object, dim)
  Mu <- par$mean
  Sigma <- par$variance
  if(missing(xlim))
    xlim <- range(dir[,1]) # +c(-1,1)*0.05*diff(range(x)))
  if(missing(ylim))
    ylim <- range(dir[,2]) # +c(-1,1)*0.05*diff(range(x)))
  xygrid <- cbind(seq(xlim[1], xlim[2], length = ngrid),
                  seq(ylim[1], ylim[2], length = ngrid))
  grid <- expand.grid(xygrid[,1], xygrid[,2])
  pred <- predict.MclustDR(object, dim = dim, eval.points = grid)
  out <- list(x = xygrid[,1], y = xygrid[,2],
              density = matrix(pred$density, ngrid, ngrid),
              z = array(pred$z, c(ngrid, ngrid, ncol(pred$z))),
              uncertainty = matrix(pred$uncertainty, ngrid, ngrid),
              classification = matrix(pred$classification, ngrid, ngrid))
  return(out)
}

plot.MclustDR <- function(x, dimens, 
                          what = c("scatterplot", "pairs", "contour", 
                                   "classification", "boundaries", "density", "evalues"), 
                          symbols, colors, 
                          col.contour = gray(0.7), 
                          col.sep = grey(0.4), 
                          ngrid = 200, 
                          nlevels = 5, 
                          asp = NULL, ...)
{ 
  object <- x
  x <- object$x
  p <- ncol(x)
  n <- nrow(x)
  G <- object$G
  y <- object$mixcomp
  class <- as.numeric(object$classification)
  nclass <- length(table(class))
  dir <- object$dir
  numdir <- object$numdir
  dimens <- if(missing(dimens)) seq(numdir)
            else intersect(as.numeric(dimens), seq(numdir))
  
  if(missing(symbols)) 
    { if(G <= length(mclust.options("classPlotSymbols"))) 
        { symbols <- mclust.options("classPlotSymbols") }
      else if(G <= 26) 
             { symbols <- LETTERS }
    }
  if(length(symbols) == 1) symbols <- rep(symbols,nclass)
  if(length(symbols) < nclass)
    { warning("more symbols needed to show classification")
      symbols <- rep(16, nclass) }
  
  if(missing(colors))
    { colors <- mclust.options("classPlotColors") }
  if(length(colors) == 1) colors <- rep(colors,nclass)
  if(length(colors) < nclass) 
    { warning("more colors needed to show mixture components")
      colors <- rep("black", nclass) }
  
  niceRange <- function (x, f = 0.04) 
  { 
    r <- range(x)
    d <- diff(r)
    out <- c(r[1] - d*f, r[2] + d*f)
    return(out)
  }

  ####################################################################
  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  # on.exit(par(oldpar))
  if(any(i <- (what == "pairs")) & (length(dimens) == 2))
    { what[i] <- "scatterplot" }
  if(length(dimens) == 1) 
    { what[!(what == "density" | what == "evalues")] <- "density" }
  what <- unique(what)
  
  plot.MclustDR.scatterplot <- function(...)
  { 
    dir <- dir[,dimens,drop=FALSE]
    plot(dir, col = colors[class], pch = symbols[class],
         xlab = colnames(dir)[1], ylab = colnames(dir)[2], 
         asp = asp, ...)
  }
  
  plot.MclustDR.pairs <- function(...)
  { dir <- dir[,dimens,drop=FALSE]
    pairs(dir, col = colors[class], pch = symbols[class], 
          gap = 0.2, asp = asp, ...)
  }
  
  plot.MclustDR.density <- function(...)
  { dimens <- dimens[1]
    dir <- object$dir[,dimens,drop=FALSE]
    par <- projpar.MclustDR(object, dimens)
    Mu <- par$mean
    Sigma <- par$variance
    q <- seq(min(dir), max(dir), length=2*ngrid)
    dens <- matrix(as.double(NA), length(q), G)
    for(j in 1:G)
      dens[,j] <- dnorm(q, Mu[j,], sqrt(Sigma[,,j]))
    #
    if(object$type == "MclustDA")
    { d <- t(apply(dens, 1, function(x, p = object$pro) p*x))
      dens <- matrix(as.double(NA), length(q), nclass) 
      tab <- table(y, class)
      for(i in 1:ncol(tab))
      { j <- which(tab[,i] > 0)
        dens[,i] <- apply(d[,j,drop=FALSE],1,sum) 
      }
    }
    #
    oldpar <- par(mar = c(0,5.1,1,1), mfrow = par("mfrow"), no.readonly = TRUE)
    on.exit(par(oldpar))
    layout(matrix(1:2,2,1), heights = c(2,1))
    plot(0, 0, type = "n", xlab = colnames(dir), ylab = "Density", 
         xlim = range(q, dir), ylim = range(0, dens*1.1), xaxt = "n")
    for(j in 1:ncol(dens))
      lines(q, dens[,j], col = colors[j])
    dir.class <- split(dir, class) 
    par(mar = c(4.1,5.1,0,1))
    boxplot(dir.class, col = adjustcolor(colors[1:nclass], alpha.f = 0.3),
            border = colors[1:nclass], horizontal = TRUE, 
            pars = list(boxwex = 0.6, staplewex = 0.8, medlwd = 2,
                        whisklty = 3, outlty = 1, outpch = NA),
            ylim = range(q,dir), yaxt = "n", xlab = colnames(dir))
    axis(2, at = 1:nclass, labels = levels(object$classification), 
         tick = FALSE, cex = 0.8, las = 2)
  }
  
  plot.MclustDR.contour <- function(...)
  { 
    dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    par <- projpar.MclustDR(object, dimens)
    Mu <- par$mean
    Sigma <- par$variance
    # draw contours for each class or cluster
    plot(dir, type = "n", asp = asp)
    for(k in seq(nclass))
    {
      i <- which(object$class2mixcomp == k)
      parameters <- list(pro = object$pro[i]/sum(object$pro[i]), 
                         mean = t(par$mean[i,,drop=FALSE]), 
                         variance = list(G = length(i), d = 2, 
                                         sigma = par$variance[,,i,drop=FALSE]))
      surfacePlot(dir, parameters, col = col.contour, 
                  nlevels = nlevels, grid = ngrid,
                  xlim = par("usr")[1:2], ylim = par("usr")[3:4], 
                  asp = asp, add = TRUE)
    }
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  plot.MclustDR.classification.Mclust <- function(...)
  { dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = niceRange(dir[,1]), 
                               ylim = niceRange(dir[,2]))
    pred$classification <- apply(pred$z, 1:2, which.max)
    image(pred$x, pred$y, pred$classification, 
          col = adjustcolor(colors[1:G], alpha.f = 0.1),
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          xaxs = "i", yaxs = "i", asp = asp)
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  plot.MclustDR.classification.MclustDA <- function(...)
  { 
    dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = niceRange(dir[,1]), 
                               ylim = niceRange(dir[,2]))
    pred$classification <- apply(pred$z, 1:2, which.max)
    image(pred$x, pred$y, pred$classification, 
          col = adjustcolor(colors[1:nclass], alpha.f = 0.1),
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          xaxs = "i", yaxs = "i", asp = asp)
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  plot.MclustDR.boundaries.Mclust <- function(...)
  { dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = niceRange(dir[,1]), 
                               ylim = niceRange(dir[,2]))
    image(pred$x, pred$y, pred$uncertainty,
          col = rev(gray.colors(10, start = 0, end = 1)),
          breaks = seq(0, 1-1/nclass, length = 11),
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          xaxs = "i", yaxs = "i", asp = asp)
    points(dir, col = colors[class], pch = symbols[class], ...)                        
  }
  
  plot.MclustDR.boundaries.MclustDA <- function(...)
  { 
    dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = niceRange(dir[,1]), 
                               ylim = niceRange(dir[,2]))
    levels <- seq(0, 1-1/nclass, length = 11)
    col <- rev(gray.colors(10, start = 0, end = 1))
    image(pred$x, pred$y, pred$uncertainty,
          col = col, breaks = levels,
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          xaxs = "i", yaxs = "i", asp = asp)
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  plot.MclustDR.evalues <- function(...)
  { plotEvalues.MclustDR(object, numdir = max(dimens), plot = TRUE) }
  
  if(interactive() & length(what) > 1)
    { title <- "Dimension reduction for model-based clustering and classification plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "scatterplot") plot.MclustDR.scatterplot(...)
             if(what[choice] == "pairs")   plot.MclustDR.pairs(...)
             if(what[choice] == "contour") plot.MclustDR.contour(...)
             if(what[choice] == "classification" & object$type == "Mclust")
                                           plot.MclustDR.classification.Mclust(...)
             if(what[choice] == "classification" &
                (object$type == "EDDA" | object$type == "MclustDA"))
                                           plot.MclustDR.classification.MclustDA(...)
             if(what[choice] == "boundaries" & object$type == "Mclust")
                                           plot.MclustDR.boundaries.Mclust(...)
             if(what[choice] == "boundaries" & 
                (object$type == "EDDA" | object$type == "MclustDA"))
                                           plot.MclustDR.boundaries.MclustDA(...)
             if(what[choice] == "density") plot.MclustDR.density(...)
             if(what[choice] == "evalues") plot.MclustDR.evalues(...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  }
  else 
    { if(any(what == "scatterplot")) plot.MclustDR.scatterplot(...)
      if(any(what == "pairs"))       plot.MclustDR.pairs(...)
      if(any(what == "contour"))     plot.MclustDR.contour(...)
      if(any(what == "classification" & object$type == "Mclust"))
                                     plot.MclustDR.classification.Mclust(...)
      if(any(what == "classification" &
             (object$type == "EDDA" | object$type == "MclustDA")))
                                     plot.MclustDR.classification.MclustDA(...)
      if(any(what == "boundaries" & object$type == "Mclust"))
                                     plot.MclustDR.boundaries.Mclust(...)
      if(any(what == "boundaries" & 
             (object$type == "EDDA" | object$type == "MclustDA")))
                                     plot.MclustDR.boundaries.MclustDA(...)
      if(any(what == "density"))     plot.MclustDR.density(...)
      if(any(what == "evalues"))     plot.MclustDR.evalues(...)
  }
    
  invisible()
}

plotEvalues.MclustDR <- function(x, numdir, plot = FALSE, legend = TRUE, ylim, ...)
{ 
  object <- x
  G <- object$G
  f <- object$pro
  lambda <- object$lambda
  # dim <- if(missing(numdir)) seq(object$numdir) else seq(numdir)
  if(missing(numdir)) numdir <- object$numdir
  dim <- seq(numdir)
  d <- length(dim)
  par <- projpar.MclustDR(object, dim = dim, center = TRUE, raw = TRUE)
  mu <- par$mean
  Sigma.G <- par$variance
  #
  M1 <- t(mu) %*% diag(f) %*% mu
  l1 <- 2*lambda*diag(crossprod(M1))
  #
  S <- matrix(0, d, d)
  for(j in seq(G)) 
    S <- S + f[j]*Sigma.G[,,j]
  M2 <- matrix(0, d, d)
  for(j in 1:G)
     { C <- (Sigma.G[,,j]-S)
       M2 <- M2 + f[j] * tcrossprod(C) }
  l2 <- 2*(1-lambda)*diag(M2)
  # 
  l <- object$evalues[dim]
  #
  if(plot)
  { if(missing(ylim)) ylim <- range(0, max(l)+diff(range(l))*0.05)
    plot(dim, l, type="b", lty = 1, pch = 16, cex = 1.5, 
         xaxt = "n", ylim = ylim, 
         xlab = "MclustDR directions", ylab = "Eigenvalues",
         panel.first = { abline(v = dim, col = "lightgray", lty = "dotted")
                         abline(h = axTicks(2,par("yaxp")), 
                                col = "lightgray", lty = "dotted") } )
    axis(1, at = dim, labels = dim)
    lines(dim, l1, type="b", lty = 2, pch = 22, cex = 1.5)
    lines(dim, l2, type="b", lty = 2, pch = 2, cex = 1.5)
    if(legend)
    { legend("topright", lty = c(1,2,2), pch = c(16,22,2), 
             legend = c("Eigenvalues", 
                        "Means contrib.", 
                        "Vars contrib."),
             bg = ifelse(par("bg")=="transparent", "white", par("bg")),
             inset = 0.01, pt.cex = 1.5) }
  }
  
  out <- list(dim = dim, evalues = l, mean.contrib = l1, var.contrib = l2)
  if(plot) invisible(out)
  else     return(out) 
}


# Auxiliary functions -----------------------------------------------------

mvdnorm <- function(x, mu, sigma, log = FALSE, tol = sqrt(.Machine$double.eps))
{
  if(is.vector(x)) 
  { x <- matrix(x, ncol = length(x)) }
  else
  { x <- as.matrix(x) }
  SVD <- svd(sigma)
  pos <- (SVD$d > max(tol*SVD$d[1], 0)) # in case of not full rank covar matrix
  inv.sigma <- SVD$v[,pos,drop=FALSE] %*% (1/SVD$d[pos] *
                                             t(SVD$u[,pos,drop=FALSE]))
  z <- mahalanobis(x, center = mu, cov = inv.sigma, inverted = TRUE)
  # logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logdet <- sum(log(SVD$d[pos]))
  logdens <- -(ncol(x) * log(2 * pi) + logdet + z)/2
  if(log) return(logdens)
  else    return(exp(logdens))
}

ellipse <- function(c, M, r, npoints = 100)
{
  # Returns the cartesian coordinates of points x on the ellipse 
  #                  (x-c)'M(x-c) = r^2,
  # where x = x(theta) and theta varies from 0 to 2*pi radians in npoints steps.
  
  # local functions
  circle <- function(theta, r) r*c(cos(theta),sin(theta))
  ellip  <- function(theta, r, lambda) lambda*circle(theta, r)
  point  <- function(theta) c+c(gamma %*% ellip(theta, r, lam))
  #
  SVD <- svd(M)
  lam   <- 1/sqrt(SVD$d)
  gamma <- SVD$v
  coord  <- t(sapply(seq(0, 2*pi, length=npoints), function(th) point(th)))
  return(coord)
}

eigen.decomp <- function(A, B, invsqrt = FALSE)
{
#
# Generalized eigenvalue decomposition of A with respect to B.
#
# A generalized eigenvalue problem AV = BLV is said to be symmetric positive 
# definite if A is symmetric and B is positive definite. V is the matrix of
# generalized eigenvectors, and L is the diagonal matrix of generalized 
# eigenvalues (Stewart, 2001, pag. 229-230).
# 
# Properties:
# V'AV = L
# V'BV = I
#
# The algorithm implemented is described in Stewart (2001, pag. 234) and used
# by Li (2000).
#
# References: 
# Li, K.C., 2000. High dimensional data analysis via the SIR-PHD approach,
# Stewart, G.W., 2001. Matrix Algorithms: vol II Eigensystems, SIAM.
  
  if(!invsqrt)
  {
    SVD <- svd(B, nu=0)
    # in case of not full rank covar matrix
    tol <- .Machine$double.eps
    pos <- which(SVD$d > max(tol*SVD$d[1], 0)) 
    SVD$d <- SVD$d[pos]
    SVD$v <- SVD$v[,pos,drop=FALSE]
    # Computes inverse square root matrix such that: 
    # t(inv.sqrt.B) %*% inv.sqrt.B = inv.sqrt.B %*% t(inv.sqrt.B) = solve(B)
    inv.sqrt.B <- SVD$v %*% (1/sqrt(SVD$d) * t(SVD$v))
  } else
  { inv.sqrt.B <-  B }
  # Compute  B^(-1/2)' A B^(-1/2) = UDU'
  # evectors = B^(-1/2) U 
  # evalues  = D
  A <- t(inv.sqrt.B) %*% A %*% inv.sqrt.B
  SVD <- svd(A, nu=0)
  list(l = SVD$d, v = inv.sqrt.B  %*% SVD$v)
}


# Subset selection of GMMDR/GMMDRC directions -----------------------------

MclustDRsubsel <- function(object, G = 1:9,
                           modelNames = mclust.options("emModelNames"), 
                           ..., 
                           bic.stop = 0, bic.cutoff = 0, 
                           mindir = 1, 
                           verbose = interactive())
{
# Subset selection for GMMDR directions based on bayes factors.
#
# object     = a MclustDR object
# G          = a vector of cluster sizes for searching 
# modelNames = a vector of models for searching 
# ...        = further arguments passed through Mclust/MclustDA
# bic.stop   = criterion to stop the search. If maximal BIC difference is
#              less than bic.stop the algorithm stops. 
#              Two tipical values are:
#              0    = stops when BIC difference becomes negative (default)
#              -Inf = stops when all directions have been selected 
# bic.cutoff = select simplest ``best'' model within bic.cutoff from the
#              maximum value achieved. Setting this to 0 (default) simply
#              select the model with the largest BIC difference.
# mindir     = the minimum number of diretions to be estimated
# verbose    = if 0 no trace info is shown; if 1 a trace of each step 
#              of the search is printed; if 2 a detailed trace info is
#              is shown.

  if(class(object) != "MclustDR") 
    stop("Not a 'MclustDR' object")
  
  hcUse <- mclust.options("hcUse")
  mclust.options("hcUse" = "VARS")
  on.exit(mclust.options("hcUse" = hcUse))  
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- switch(object$type, 
                    "Mclust"   = as.name("MclustDRsubsel_cluster"),
                    "EDDA"     = as.name("MclustDRsubsel_classif"),
                    "MclustDA" = as.name("MclustDRsubsel_classif"),
                    stop("Not allowed 'MclustDR' type!"))
  eval(mc, parent.frame())
}

MclustDRsubsel_cluster <- function(object, G = 1:9,
                                   modelNames = mclust.options("emModelNames"), 
                                   ..., 
                                   bic.stop = 0, bic.cutoff = 0, 
                                   mindir = 1, 
                                   verbose = interactive())
{
  drmodel <- object
  mclustType <- drmodel$type
  lambda <- drmodel$lambda
  numdir <- drmodel$numdir
  numdir0 <- numdir+1
  dir <- drmodel$dir[,seq(numdir),drop=FALSE]
  mindir <- max(1, as.numeric(mindir), na.rm = TRUE)
  verbose <- as.numeric(verbose)
  ncycle <- 0
  while(numdir < numdir0)
  { ncycle <- ncycle+1
    if(verbose > 0) cat("\nCycle", ncycle, "...\n")
    out <- MclustDRsubsel1cycle(drmodel, 
                                G, modelNames, 
                                bic.stop = bic.stop, 
                                bic.cutoff = bic.cutoff, 
                                verbose = if(verbose > 1) TRUE else FALSE)
    if(verbose > 0) { cat("\n"); print(out$tab) }
    mod <- do.call("Mclust", 
                   list(data = dir[,out$subset,drop=FALSE], 
                        G = G, 
                        modelNames = if(length(out$subset) > 1)
                                        modelNames else c("E", "V"),
                        verbose = FALSE, ...))
    numdir0 <- numdir
    drmodel0 <- MclustDR(mod, lambda = lambda)
    if(drmodel0$numdir < mindir) break
    drmodel <- drmodel0
    numdir <- drmodel$numdir
    dir <- drmodel$dir[,seq(numdir),drop=FALSE]
  }
  
  # format object using original data
  obj <- drmodel
  obj$basisx <- MclustDRrecoverdir(obj, data = object$x, std = FALSE)
  obj$std.basisx <- MclustDRrecoverdir(obj, data = object$x, std = TRUE)
  class(obj) <- c("MclustDRsubsel", class(obj))
  return(obj)
}
  
MclustDRsubsel1cycle <- function(object, 
                                 G = 1:9, 
                                 modelNames = mclust.options("emModelNames"), 
                                 bic.stop = 0, bic.cutoff = 0, 
                                 verbose = interactive())
{
# Single cycle of subset selection for GMMDR directions based on bayes factors.

  if(class(object) != "MclustDR") 
    stop("Not a 'MclustDR' object")

  d <- object$numdir
  dir <- object$dir[,seq(d),drop=FALSE]
  n <- nrow(dir)
  if(is.null(colnames(dir)))
     colnames(dir) = paste("[,", 1:d, "]", sep="")
  dir.names <- colnames(dir)
  BIC <- Model1 <- Model2 <- tab <- NULL; Model1$bic <- 0
  bic.stop <- bic.stop + sqrt(.Machine$double.eps)
  bic.cutoff <- bic.cutoff + sqrt(.Machine$double.eps)
  inc <- NULL; excl <- seq(1,d)
  model1D <- if(any(grep("V", modelNames))) c("E", "V") else "E"
  #
  hskip <- paste(rep(" ",4),collapse="")
  if(verbose) cat("\n", hskip, "Start greedy search...\n", sep="")
  
  while(length(excl)>0)
       { if(verbose) 
           { cat(hskip, rep("-",15), "\n", sep="")
             cat(paste(hskip, "Step", length(inc)+1, "\n")) }
         for(j in excl)
            { # Select simplest model with smallest num. of components 
              # within bic.cutoff
              bic <- mclustBIC(dir[,sort(c(inc, j)),drop=FALSE], 
                               G = G, 
                               modelNames = if(length(inc)>0) modelNames else model1D,
                               verbose = FALSE)
              bic.tab <- (as.matrix(max(bic, na.rm=TRUE) - bic) < bic.cutoff)*1
              bestG   <- which(rowSums(bic.tab, na.rm=TRUE) > 0)[1]
              bestmod <- which(bic.tab[bestG,,drop=FALSE] > 0)[1]
              out <- data.frame(Variable = dir.names[j],
                                model = colnames(bic.tab)[bestmod],
                                G = G[bestG],
                                bic = c(bic[bestG,bestmod]),
                                bic.diff = c(bic[bestG,bestmod] - 
                                             Model1$bic -
                                             MclustDRBICreg(dir[,j], dir[,inc]))
                               )
              #
              Model2 <- rbind(Model2, out)
            }
         if(verbose) print(cbind(" " = hskip, Model2), row.names = FALSE)
         
         # stop if max BIC difference is < than cut-off bic.stop
         if(max(Model2$bic.diff) < bic.stop & length(inc) > 0) 
           { break }
         # otherwise keep selecting
         i <- which.max(Model2$bic.diff)
         inc <- append(inc, excl[i])
         excl <- setdiff(excl, excl[i])
         tab <- rbind(tab, Model2[i,])
         Model1 <- Model2[i,]
         Model2 <- NULL
       }
  rownames(tab) <- 1:nrow(tab)
  colnames(tab) <- c("Variable", "Model", "G", "BIC", "BIC.dif")
  subsets <- sapply(1:nrow(tab), function(x) list(inc[1:x]))
  #
  return(list(subset = subsets[[length(subsets)]], tab = tab))
}

MclustDRsubsel_classif <- function(object, 
                                   G = 1:9, modelNames = mclust.options("emModelNames"), 
                                   ..., 
                                   bic.stop = 0, bic.cutoff = 0, 
                                   mindir = 1, 
                                   verbose = interactive())
{
  drmodel <- object
  mclustType <- drmodel$type
  lambda <- drmodel$lambda
  numdir <- drmodel$numdir
  numdir0 <- numdir+1
  dir <- drmodel$dir[,seq(numdir),drop=FALSE]
  mindir <- max(1, as.numeric(mindir), na.rm = TRUE)
  verbose <- as.numeric(verbose)
  ncycle <- 0
  while(numdir < numdir0)
  { ncycle <- ncycle+1
    if(verbose > 0) cat("\nCycle", ncycle, "...\n")
    out <- MclustDRCsubsel1cycle(drmodel, 
                                 G, modelNames, 
                                 bic.stop = bic.stop, 
                                 bic.cutoff = bic.cutoff, 
                                 verbose = if(verbose > 1) TRUE else FALSE)
    if(verbose > 0) { cat("\n"); print(out$tab) }
    mod <- do.call("MclustDA", list(data = dir[,out$subset,drop=FALSE], 
                                    class = object$classification,
                                    G = G, 
                                    modelNames = if(length(out$subset) > 1) modelNames
                                                 else if(any(grep("V", modelNames))) 
                                                        c("E", "V") else "E",
                                    modelType = mclustType,
                                    verbose = FALSE, ...))
    numdir0 <- numdir
    drmodel0 <- MclustDR(mod, lambda = lambda)
    if(drmodel0$numdir < mindir) break
    drmodel <- drmodel0
    numdir <- drmodel$numdir
    dir <- drmodel$dir[,seq(numdir),drop=FALSE]
  }
  
  # format object using original data
  obj <- drmodel
  obj$basisx <- MclustDRrecoverdir(obj, data = object$x, std = FALSE)
  obj$std.basisx <- MclustDRrecoverdir(obj, data = object$x, std = TRUE)
  class(obj) <- c("MclustDRsubsel", class(obj))
  return(obj)
}
  
MclustDRCsubsel1cycle <- function(object, 
                                  G = 1:9, 
                                  modelNames = mclust.options("emModelNames"), 
                                  bic.stop = 0, bic.cutoff = 0, 
                                  verbose = TRUE)
{
# Single cycle of subset selection for GMMDRC directions based on bayes factors.

  if(class(object) != "MclustDR") 
    stop("Not a 'MclustDR' object")

  d <- object$numdir
  dir <- object$dir[,seq(d),drop=FALSE]
  n <- nrow(dir)
  if(is.null(colnames(dir)))
     colnames(dir) = paste("[,", seq(d), "]", sep="")
  dir.names <- colnames(dir)
  BIC <- Model1 <- Model2 <- tab <- NULL; Model1$bic <- 0
  bic.stop <- bic.stop + sqrt(.Machine$double.eps)
  bic.cutoff <- bic.cutoff + sqrt(.Machine$double.eps)
  inc <- NULL; excl <- seq(d)
  model1D <- if(any(grep("V", modelNames))) c("E", "V") else "E"
  #
  hskip <- paste(rep(" ",4),collapse="")
  if(verbose) cat("\n", hskip, "Start greedy search...\n", sep="")
  
  while(length(excl)>0)
       { if(verbose) 
           { cat(hskip, rep("-",15), "\n", sep="")
             cat(paste(hskip, "Step", length(inc)+1, "\n")) }
         for(j in excl)
            { # Select simplest model with smallest num. of components 
              # within bic.cutoff
              mod <- MclustDA(dir[,sort(c(inc, j)),drop=FALSE], 
                              class = object$classification,
                              G = G, 
                              modelNames = if(length(inc)>0) modelNames else model1D,
                              modelType = object$type,
                              verbose = FALSE)
              out <- data.frame(Variable = dir.names[j],
                                model = paste(sapply(mod$models, function(m) m$modelName),collapse="|"),
                                G = paste(sapply(mod$models, function(m) m$G),collapse="|"),
                                bic = mod$bic,
                                bic.diff = c(mod$bic - 
                                             # (Model1$bic + bic.reg(z2, z1))
                                             Model1$bic -
                                             MclustDRBICreg(dir[,j], dir[,inc]))
                               )
              #
              Model2 <- rbind(Model2, out)
            }
         if(verbose) print(cbind(" " = hskip, Model2), row.names = FALSE)
         
         # stop if max BIC difference is < than cut-off bic.stop
         if(max(Model2$bic.dif) < bic.stop & length(inc) > 0) 
           { break }
         # otherwise keep selecting
         i <- which.max(Model2$bic.dif)
         inc <- append(inc, excl[i])
         excl <- setdiff(excl, excl[i])
         tab <- rbind(tab, Model2[i,])
         Model1 <- Model2[i,]
         Model2 <- NULL
       }
  rownames(tab) <- 1:nrow(tab)
  colnames(tab) <- c("Variable", "Model", "G", "BIC", "BIC.dif")
  subsets <- sapply(1:nrow(tab), function(x) list(inc[1:x]))
  #
  return(list(subset = subsets[[length(subsets)]], tab = tab))
}

# BICreg <- function(y, x)
# { n <- length(y)
#   mod <- lm.fit(cbind(rep(1,n), x), y)
#   rss <- sum(mod$residuals^2)
#   -n*log(2*pi) - n*log(rss/n) - n - (n - mod$df.residual + 1) * log(n) 
# }

MclustDRBICreg <- function(y, x, stepwise = TRUE)
{ 
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- length(y)
  mod0 <- lm(y ~ 1)
  if(ncol(x) >= 1)
    { mod1 <- lm(y ~ 1+x)
      if(stepwise) 
        { mod <- step(mod0, k = log(n), trace = 0,
                      scope = list(lower = formula(mod0), 
                                   upper = formula(mod1)),
                      direction = "forward") }
       else mod <- mod1
  }
  else mod <- mod0
  rss <- sum(mod$residuals^2)
  p <- (n - mod$df.residual + 1)
  -n*log(2*pi) - n*log(rss/n) - n - p*log(n) 
}

normalize <- function(x)
{
# Normalize the vector x to have unit length
  x <- as.vector(x)
  x <- x/sqrt(as.vector(crossprod(x)))
  return(x)
}

MclustDRrecoverdir <- function(object, data, normalized = TRUE, std = FALSE)
{
# Recover coefficients of the linear combination defining the MclustDR
# directions. This is useful if the directions are obtained from other 
# directions

  if(!any(class(object) == "MclustDR"))
    stop("object must be of class mclustsir")
  if(missing(data)) x <- object$x
  else              x <- as.matrix(data)
  x <- scale(x, center=TRUE, scale=FALSE)
  numdir <- object$numdir
  dir <- object$dir[,seq(numdir),drop=FALSE]
  # B <- as.matrix(coef(lm(dir ~ x)))[-1,,drop=FALSE] # ok but old
  B <- qr.solve(x, dir)  
  if(std) 
    { sdx <- sd(x)
      B <- apply(B, 2, function(x) x*sdx) }
  if(normalized) 
    { B <- as.matrix(apply(B, 2, normalize)) }
  rownames(B) <- colnames(x)
  return(B)
}

## Define print and summary methods for showing basis coefs 
## in the original scale of variables

print.MclustDRsubsel <- function(x, ...)
{
  x$basis <- x$basisx
  class(x) <- class(x)[2]
  NextMethod()
} 
  
summary.MclustDRsubsel <- function(object, ...)
{
  object$basis <- object$basisx
  object$std.basis <- object$std.basisx
  class(object) <- class(object)[2]
  NextMethod()
} 

