######################################################
##                                                  ##
##        Dimension reduction for model-based       ##
##          clustering and classification           ##
##                                                  ##
## Author: Luca Scrucca                             ##
######################################################

MclustDR <- function(object, normalized = TRUE, Sigma, lambda = 0.5, tol = sqrt(.Machine$double.eps))
{
  #  Dimension reduction for model-based clustering and classification
  #
  # object = a object of class '"Mclust' 
  # data = the data used to produce object.
  # normalized = normalize direction coefs to have unit norm
  
  call <- match.call()
  if(!any(class(object) == c("Mclust", "MclustDA")))
    stop("object must be of class 'Mclust' or 'MclustDA'")
  
  normalize <- function(x) 
  { x <- as.vector(x)
    x/sqrt(crossprod(x))
  }
  
  x <- data.matrix(object$data)
  p <- ncol(x)
  n <- nrow(x)
  #-----------------------------------------------------------------  
  # overall parameters
  mu <- colMeans(x)
  if(missing(Sigma)) Sigma <- var(x)*(n-1)/n
  # within-cluster parameters based on fitted mixture model
  if(class(object) == "Mclust")
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
    # within-group means
    mu.G <- matrix(par$mean,p,G) 
    # within-group covars
    if(p == 1) 
    { Sigma.G <- array(par$variance$sigmasq, c(p,p,G)) }
    else     
    { Sigma.G <- par$variance$sigma }
  }
  else if(class(object) == "MclustDA")
  { 
    type <- object$type
    modelName <- sapply(object$models, function(m) m$modelName)
    class <- object$class
    class <- factor(class, levels = names(object$models))
    y <- rep(NA, length(class))
    for(i in 1:nlevels(class))
    { y[class == levels(class)[i]] <- paste(levels(class)[i], 
                                            object$models[[i]]$classification, sep =":") }
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
  SVD <- svd(Sigma)
  pos <- (SVD$d > max(tol*SVD$d[1], 0)) # in case of not full rank covar matrix
  if(all(pos)) 
  { inv.Sigma <- SVD$v %*% (1/SVD$d * t(SVD$u))
    inv.sqrt.Sigma <- SVD$v %*% (1/sqrt(SVD$d) * t(SVD$u)) }
  else
  { inv.Sigma <- SVD$v[,pos,drop=FALSE] %*% 
      (1/SVD$d[pos] * t(SVD$u[,pos,drop=FALSE]))
    inv.sqrt.Sigma <- SVD$v[,pos,drop=FALSE] %*% 
      (1/sqrt(SVD$d[pos]) * t(SVD$u[,pos,drop=FALSE])) }
  #-----------------------------------------------------------------
  # pooled within-group covariance
  S <- matrix(0, p, p)
  for(j in 1:G) 
    S <- S + f[j]*Sigma.G[,,j]
  #-----------------------------------------------------------------
  # kernel matrix
  M.I <- crossprod(t(sweep(mu.G, 1, FUN="-", STATS=mu))*sqrt(f))
  M.II <- matrix(0, p, p)
  for(j in 1:G)
    M.II <- M.II + f[j]*crossprod(inv.sqrt.Sigma%*%(Sigma.G[,,j]-S))
  # convex combiation of M_I and M_II
  M <- 2*lambda*crossprod(inv.sqrt.Sigma %*% M.I) + 2*(1-lambda)*M.II
  # regularize the M_II
  # M <- M.I + lambda*M.II
  # M <- crossprod(inv.sqrt.Sigma %*% M.I) + 
  #      (1-lambda)*M.II + lambda/p * diag(p)  
  #
  SVD <- eigen.decomp(M, Sigma)
  l <- SVD$l; l <- (l+abs(l))/2
  numdir <- min(p, sum(l > sqrt(.Machine$double.eps)))
  basis <- as.matrix(SVD$v)[,1:numdir,drop=FALSE]
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
             class = class, mixcomp = y, class2mixcomp = cl2mc,
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

print.MclustDR <- function(x, ...) 
{
  cat(paste("\'", class(x), "\' object for ", x$type, 
            " mixture model:\n\n", sep = ""))
  tab <- rbind(x$basis[,seq(x$numdir),drop=FALSE],
               "-----------" = rep(NA, x$numdir), 
               Eigenvalues = x$evalues[seq(x$numdir)])
  print(tab, na.print = "", ...)
  invisible()
}

summary.MclustDR <- function(object, numdir, std = FALSE, ...)
{
  if(missing(numdir)) numdir <- object$numdir
  dim <- 1:numdir
  
  if(object$type == "Mclust")
  { n <- as.vector(table(object$class))
    G <- object$G }
  else
  { n <- as.vector(table(object$class))
    G <- as.vector(table(object$class2mixcomp)) }
  
  obj <- list(type = object$type, 
              modelName = object$modelName, 
              classes = levels(object$class),
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
  cat(rep("-", nchar(title)),"\n",sep="")
  cat(title, "\n")
  cat(rep("-", nchar(title)),"\n",sep="")
  
  if(x$type == "Mclust")
  { 
    tab <- data.frame(n = x$n)
    rownames(tab) <- x$classes
    tab <- as.matrix(tab)
    names(dimnames(tab)) <- c("Clusters", "")
    cat(paste("\nMixture model type: ", x$type, 
              " (", x$modelName, ", ", x$G, ")\n", sep = ""))
    print(tab, quote = FALSE, right = TRUE)
  }
  else if(x$type == "MclustDA" | x$type == "EDDA")
  {
    tab <- data.frame(n = x$n, Model = x$modelName, G = x$G)
    rownames(tab) <- x$classes
    tab <- as.matrix(tab)
    names(dimnames(tab)) <- c("Classes", "")
    cat(paste("\nMixture model type:", x$type, "\n"))
    print(tab, quote = FALSE, right = TRUE)
  }
  else stop("invalid model type")
  
  if(x$std) 
  { cat("\nStandardized basis vectors using predictors \nscaled to have std.dev. equal to one:\n")
    print(x$std.basis, digits = digits)
  }
  else 
  { cat("\nEstimated basis vectors:\n")
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
  
  #plot(object$dir[,dim], col = as.numeric(object$class))
  #points(par$mean, col = object$class2mixcomp, pch = 3, cex = 2)
  
  #plot(x[,1:2], col = as.numeric(object$class))
  #points(mu[,1:2], col = object$class2mixcomp, pch = 3, cex = 2)
  
}  

predict.MclustDR <- function(object, dim = 1:object$numdir, newdata, eval.points, ...)
{  
  dim <- dim[dim <= object$numdir]
  if(missing(newdata) & missing(eval.points))
  { dir <- object$dir[,dim,drop=FALSE] }
  else if(!missing(newdata))
  { newdata <- as.matrix(newdata) 
    newdata <- scale(newdata, center = colMeans(object$x), scale = FALSE)
    dir <- newdata %*% object$basis[,dim,drop=FALSE]
  }
  else if(!missing(eval.points))
  { dir <- as.matrix(eval.points) }
  
  n <- nrow(dir)
  G <- object$G # num. components
  nclass <- nlevels(object$class) # num. classes
  par <- projpar.MclustDR(object, dim)
  Mu <- par$mean
  Sigma <- par$variance
  
  cden <- array(NA, c(n, G))
  for(j in 1:G)
  { cden[,j] <- mvdnorm(dir, Mu[j,], Sigma[,,j], log = FALSE) }
  z <- sweep(cden, 2, FUN = "*", STATS = object$pro)  
  den <- apply(z, 1, sum)
  z <- sweep(z, 1, FUN = "/", STATS = den)
  zz <- matrix(0, n, nclass)
  for(j in seq(nclass))
  { zz[,j] <- rowSums(z[,object$class2mixcomp == j,drop=FALSE]) }
  z <- zz; rm(zz)
  
  class <- factor(apply(z,1,which.max), 
                  levels = 1:nclass, 
                  labels = levels(object$class))
  
  out <- list(dir = dir,
              density = den,
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

plot.MclustDR <- function(x, dimens, what = c("scatterplot", "pairs", "contour", "classification", "boundaries", "density", "evalues"), symbols, colors, col.contour = gray(0.7), col.sep = grey(0.4), ngrid = 100, nlevels = 5, asp = NULL, ...)
{ 
  object <- x
  x <- object$x
  p <- ncol(x)
  n <- nrow(x)
  G <- object$G
  y <- object$mixcomp
  class <- as.numeric(object$class)
  nclass <- length(table(class))
  dir <- object$dir
  numdir <- object$numdir
  if(missing(dimens)) dimens <- seq(numdir)
  
  what <- match.arg(what)
  if(what == "pairs")
  { if(length(dimens) == 2) what <- "scatterplot" }  
  if(length(dimens) == 1) 
  { if(!(what == "density" | what == "evalues"))
    what <- "density"
  }
  
  if(missing(symbols)) 
  { if(G <= length(.mclust$classPlotSymbols)) 
  { symbols <- .mclust$classPlotSymbols }
  else if(G <= 26) 
  { symbols <- LETTERS }
  }
  if(length(symbols) == 1) symbols <- rep(symbols,nclass)
  if(length(symbols) < nclass)
  { warning("more symbols needed to show classification")
    symbols <- rep(16, nclass) }
  
  if(missing(colors))
  { colors <- .mclust$classPlotColors }
  if(length(colors) == 1) colors <- rep(colors,nclass)
  if(length(colors) < nclass) 
  { warning("more colors needed to show mixture components")
    colors <- rep("black", nclass) }
  
  if(what == "scatterplot")
  { dir <- dir[,dimens,drop=FALSE]
    plot(dir, col = colors[class], pch = symbols[class],
         xlab = colnames(dir)[1], ylab = colnames(dir)[2], 
         asp = asp, ...)
  }
  
  if(what == "pairs")
  { dir <- dir[,dimens,drop=FALSE]
    pairs(dir, col = colors[class], pch = symbols[class], 
          gap = 0.25, asp = asp, ...)
  }
  
  if(what == "density")
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
    axis(2, at = 1:nclass, labels = levels(object$class), tick = FALSE, cex = 0.8, las = 2)
  }
  
  if(what == "contour") 
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
      surfacePlot(dir, parameters, col = col.contour, nlevels = nlevels,
                  xlim = par("usr")[1:2], ylim = par("usr")[3:4], 
                  asp = asp, add = TRUE)
    }
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  if(what == "classification" & object$type == "Mclust")
  { dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    plot.new()
    plot.window(xlim = range(dir[,1]), ylim = range(dir[,2]), asp = asp)
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = par("usr")[1:2], 
                               ylim = par("usr")[3:4])
    pred$classification <- apply(pred$z, 1:2, which.max)
    #
    image(pred$x, pred$y, pred$classification, 
          col = adjustcolor(colors[1:G], alpha.f = 0.1),
          #xlim = par("usr")[1:2], ylim = par("usr")[3:4],
          xaxs = "i", yaxs = "i",
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          useRaster = TRUE)
    for(j in 1:G)         
    { z <- ifelse(pred$classification == j, 1, -1)
      contour(pred$x, pred$y, z, col = col.sep,
              add = TRUE, levels = 0, drawlabels = FALSE) 
    }
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  if(what == "classification" & 
       (object$type == "EDDA" | object$type == "MclustDA"))
  { 
    dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    plot.new()
    plot.window(xlim = range(dir[,1]), ylim = range(dir[,2]), asp = asp)      
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = par("usr")[1:2], 
                               ylim = par("usr")[3:4])
    pred$classification <- apply(pred$z, 1:2, which.max)
    #
    image(pred$x, pred$y, pred$classification, 
          col = adjustcolor(colors[1:nclass], alpha.f = 0.1),
          xlim = par("usr")[1:2], ylim = par("usr")[3:4],
          xaxs = "i", yaxs = "i",
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          useRaster = TRUE)
    for(j in 1:nclass)
    { z <- ifelse(pred$classification == j, 1, -1)
      contour(pred$x, pred$y, z, col = col.sep,
              add = TRUE, levels = 0, drawlabels = FALSE) 
    }
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  if(what == "boundaries" & object$type == "Mclust")
  { dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    plot.new()
    plot.window(xlim = range(dir[,1]), ylim = range(dir[,2]), asp = asp)      
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = par("usr")[1:2], 
                               ylim = par("usr")[3:4])
    image(pred$x, pred$y, pred$uncertainty, 
          col = rev(gray.colors(10, start = 0, end = 1)),
          breaks = seq(0, 1-1/nclass, length = 11),
          xlim = par("usr")[1:2], ylim = par("usr")[3:4],
          xaxs = "i", yaxs = "i",
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          useRaster = TRUE)
    points(dir, col = colors[class], pch = symbols[class], ...)                        
  }
  
  if(what == "boundaries" & 
       (object$type == "EDDA" | object$type == "MclustDA"))
  { 
    dimens <- dimens[1:2]
    dir <- object$dir[,dimens,drop=FALSE]
    plot.new()
    plot.window(xlim = range(dir[,1]), ylim = range(dir[,2]), asp = asp)
    pred <- predict2D.MclustDR(object, dimens, ngrid,
                               xlim = par("usr")[1:2], 
                               ylim = par("usr")[3:4])
    levels <- seq(0, 1-1/nclass, length = 11)
    col <- rev(gray.colors(10, start = 0, end = 1))
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    #
    # this will add a legend on the right. Useful??
    # w <- (1 + oldpar$mar[2L]) * par("csi") * 2.54
    # layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    # par(las = 1)
    # mar <- oldpar$mar; mar[4L] <- min(2.5,mar[2L]); mar[2L] <- 1; par(mar = mar)
    # plot.new()
    # plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", yaxs = "i")
    # rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    # axis(4); box()
    #
    # mar <- oldpar$mar; mar[4L] <- 1; par(mar = mar)
    image(pred$x, pred$y, pred$uncertainty,
          col = col, breaks = levels,
          xlim = oldpar$usr[1:2], ylim = oldpar$usr[3:4],
          # xlim = range(dir[,1]), ylim = range(dir[,2]), 
          # xaxs = "i", yaxs = "i",
          xlab = colnames(dir)[1], ylab = colnames(dir)[2],
          useRaster = TRUE)
    points(dir, col = colors[class], pch = symbols[class], ...)
  }
  
  if(what=="evalues")
  { plotEvalues.MclustDR(object, numdir = max(dimens), plot = TRUE) }
  
  return(invisible())
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
  for(j in seq(G)) S <- S + f[j]*Sigma.G[,,j]
  M2 <- matrix(0, d, d)
  for(j in 1:G)
  { C <- (Sigma.G[,,j]-S)
    M2 <- M2 + f[j] * C %*% t(C) }
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
                                col = "lightgray", lty = "dotted") 
         })
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


# Auxiliary functions

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

eigen.decomp <- function(A, B)
{
  #
  # Generalized eigenvalue decomposition of A with respect to B
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
  
  svd <- svd(B, nu=0)
  p <- length(svd$d)
  # Computes inverse square root matrix such that: 
  # t(inv.sqrt.B) %*% inv.sqrt.B = inv.sqrt.B %*% t(inv.sqrt.B) = solve(B)
  inv.sqrt.B <- svd$v %*% diag(1/sqrt(svd$d),p,p) %*% t(svd$v)
  # Compute  B^(-1/2)' A B^(-1/2) = UDU'
  # evectors = B^(-1/2) U 
  # evalues  = D
  A <- t(inv.sqrt.B) %*% A %*% inv.sqrt.B
  svd <- svd(A, nu=0)
  list(l = svd$d, v = inv.sqrt.B  %*% svd$v)
}
