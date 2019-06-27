mclust1Dplot <- function(data, parameters = NULL, z = NULL,
                         classification = NULL, truth = NULL, uncertainty = NULL, 
                         what = c("classification", "density", "error", "uncertainty"), 
                         symbols = NULL, colors = NULL, ngrid = length(data),  
                         xlab = NULL, ylab = NULL, 
                         xlim = NULL, ylim = NULL, 
                         CEX  = 1, main = FALSE, ...) 
{
  p <- ncol(as.matrix(data))
  if (p != 1) 
    stop("for one-dimensional data only")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if (!is.null(parameters)) 
  {
    mu <- parameters$mean
    L <- ncol(mu)
    sigmasq <- parameters$variance$sigmasq
    haveParams <- !is.null(mu) && !is.null(sigmasq) && 
      !any(is.na(mu)) && !any(is.na(sigmasq)) 
  }
  else haveParams <- FALSE
  if (is.null(xlim)) xlim <- range(data)
  if (haveParams) 
  {
    G <- length(mu)
    if ((l <- length(sigmasq)) == 1) {
      sigmasq <- rep(sigmasq, G)
    }
    else if (l != G) 
    {
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
      { colors <- mclust.options("classPlotColors")[1:L] }
    else if(length(colors) == 1) 
      { colors <- rep(colors, L) }
    else if(length(colors) < L)
      { warning("more colors needed to show classification")
      colors <- rep("black", L) }
  }
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  
  what <- match.arg(what, choices = eval(formals(mclust1Dplot)$what))
  bad <- what == "classification" && is.null(classification)
  bad <- bad || (what == "uncertainty" && is.null(uncertainty))
  bad <- bad || (what == "error" && 
                 (is.null(classification) || is.null(truth)))
  if(bad) 
    stop("insufficient input for specified plot")

  M <- L
  switch(what,
         "classification" = 
         { 
           plot(data, seq(from = 0, to = M, length = n), type = "n", 
                xlab = if(is.null(xlab)) "" else xlab,
                ylab = if(is.null(ylab)) "Classification" else ylab, 
                xlim = xlim, 
                ylim = if(is.null(ylim)) grDevices::extendrange(r = c(0,M), f = 0.1) else ylim,
                yaxt = "n", main = "", ...)
           axis(side = 2, at = 0:M, labels = c("", sort(unique(classification))))
           if(main) title("Classification")
           for(k in 1:L)
           {
             I <- classification == U[k]
             if(symbols[k] == "|")
             { vpoints(data[I], rep(0, length(data[I])), cex = CEX)
               vpoints(data[I], rep(k, length(data[I])), 
                       col = colors[k], cex = CEX)
             } else
             { points(data[I], rep(0, length(data[I])),
                      pch = symbols[k], cex = CEX)
               points(data[I], rep(k, length(data[I])),
                      pch = symbols[k], col = colors[k], cex = CEX)
             }
           }
         },
         "error" = 
         { 
           ERRORS <- classError(classification, truth)$misclassified
           plot(data, seq(from = 0, to = M, length = n), type = "n", 
                xlab = xlab,
                ylab = if(is.null(ylab)) "Class errors" else ylab, 
                xlim = xlim, 
                ylim = if(is.null(ylim)) grDevices::extendrange(r = c(0,M), f = 0.1) else ylim,
                yaxt = "n", ...)
           axis(side = 2, at = 0:M, labels = c("", unique(classification)))
           if(main) title("Classification error")
           good <- rep(TRUE, length(classification))
           good[ERRORS] <- FALSE
           sym <- "|"
           for(k in 1:L) 
           {
             K <- classification == U[k]
             I <- (K & good)
             if(any(I)) 
             { if(FALSE) 
                 { sym <- if (L > 4) 1 else if (k == 4) 5 else k - 1 }
               l <- sum(as.numeric(I))
               if(sym == "|")
                 vpoints(data[I], rep(0, l), col = colors[k], cex = CEX)
               else
                 points(data[I], rep(0, l), pch = sym, col = colors[k], cex = CEX)
             }
             I <- K & !good
             if(any(I)) 
             { if(FALSE) { sym <- if (L > 5) 16 else k + 14 }
               l <- sum(as.numeric(I))
               if(sym == "|") 
                 vpoints(data[I], rep(k, l), col = colors[k], cex = CEX)
               else
                 points(data[I], rep(k, l), pch = sym, col = colors[k], cex = CEX)
             }
           }
         },
         "uncertainty" = 
         { 
           u <- (uncertainty - min(uncertainty))/
                (max(uncertainty) - min(uncertainty) + sqrt(.Machine$double.eps))
           b <- bubble(u, cex = CEX*c(0.3, 2), alpha = c(0.3, 1))
           if(is.null(classification))
           { 
             classification <- rep(1, length(u))
             U <- 1 
           }
           if(is.null(colors))
             colors <- palette()[1]
           cl <- sapply(classification, function(cl) which(cl == U))
           plot(data, uncertainty, type = "h", 
                xlab = xlab, 
                ylab = if(is.null(ylab)) "Uncertainty" else ylab, 
                xlim = xlim, 
                ylim = if(is.null(ylim)) c(0,1) else ylim,
                col = mapply(adjustcolor, 
                             col = colors[cl], 
                             alpha.f = b$alpha),
                ...)
           rug(data, lwd = 1, col = adjustcolor(par("fg"), alpha.f = 0.8))
           if(main) title("Uncertainty")
         },
         "density" = 
         { 
           if(is.null(parameters$pro) && parameters$variance$G != 1) 
             stop("mixing proportions missing")
           x <- grid1(n = ngrid, range = xlim, edge = TRUE)
           plot(x, dens("V", data = x, parameters = parameters),
                xlab = xlab, 
                ylab = if(is.null(ylab)) "Density" else ylab, 
                xlim = xlim, 
                type = "l", main = "", ...)
           if(main) title("Density")
         },
         { 
           plot(data, rep(0, n), type = "n", xlab = "", ylab = "", 
                xlim = xlim, main = "", ...)
           vpoints(data, rep(0, n), cex = CEX)
           if(main) title("Point Plot")
         }
  )
  invisible()
}

mclust2Dplot <- function(data, parameters = NULL, z = NULL,
                         classification = NULL, truth = NULL, 
                         uncertainty = NULL, 
                         what = c("classification", "uncertainty", "error"), 
                         addEllipses = TRUE,
                         fillEllipses = mclust.options("fillEllipses"),
                         symbols = NULL, colors = NULL, 
                         xlim = NULL, ylim = NULL, 
                         xlab = NULL, ylab = NULL, 
                         scale = FALSE, CEX = 1, PCH = ".", 
                         main = FALSE, swapAxes = FALSE, ...) 
{
  if(dim(data)[2] != 2)
    stop("data must be two dimensional")
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) 
    { mu <- parameters$mean
      L <- ncol(mu)
      sigma <- parameters$variance$sigma
      haveParams <- !is.null(mu) && !is.null(sigma) && 
                    !any(is.na(mu)) && !any(is.na(sigma)) 
  }
  else haveParams <- FALSE
  
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  if(is.null(xlim))
    xlim <- range(data[, 1])
  if(is.null(ylim))
    ylim <- range(data[, 2])
  if(scale) 
    { par(pty = "s")
      d <- diff(xlim) - diff(ylim)
      if(d > 0) { ylim <- c(ylim[1] - d/2, ylim[2] + d/2.) }
      else      { xlim <- c(xlim[1] + d/2, xlim[2] - d/2) }
  }
  
  dnames <- dimnames(data)[[2]]
  if(is.null(xlab)) 
    { xlab <- if(is.null(dnames)) "" else dnames[1] }
  if(is.null(ylab)) 
    { ylab <- if(is.null(dnames)) "" else dnames[2] }
  
  if(haveParams) 
    { G <- ncol(mu)
      dimpar <- dim(sigma)
      if(length(dimpar) != 3) 
        { haveParams <- FALSE
          warning("covariance must be a 3D matrix")
      }
      if(G != dimpar[3])
        { haveParams <- FALSE
          warning("means and variance parameters are incompatible")
      }
      mu <- array(mu, c(2, G))
      sigma <- array(sigma, c(2, 2, G))
    }
  
  if(swapAxes)
    { if(haveParams) 
        { mu <- mu[2:1,]
          sigma <- sigma[2:1, 2:1,]
      }
    data <- data[, 2:1]
  }
  
  if(!is.null(truth)) 
    { if(is.null(classification)) 
        { classification <- truth
          truth <- NULL
      }
      else 
        { if(length(unique(truth)) != 
             length(unique(classification))) 
                truth <- NULL
           else truth <- as.character(truth)
      }
  }
  
  if(charmatch("classification", what, nomatch = 0) && 
       is.null(classification) && !is.null(z))
    { classification <- map(z) }
  
  if(!is.null(classification)) 
    { classification <- as.character(classification)
      U <- sort(unique(classification))
      L <- length(U)
      noise <- (U[1] == "0")
      if(is.null(symbols))
      { if(L <= length(mclust.options("classPlotSymbols"))) 
      { symbols <- mclust.options("classPlotSymbols")[1:L]
        if(noise)
        { symbols <- c(16,symbols)[1:L] }
      }
      else if(L <= 9)
      { symbols <- as.character(1:9) }
      else if(L <= 26) { symbols <- LETTERS }
      }
      if(is.null(colors)) 
      { if(L <= length(mclust.options("classPlotColors"))) 
      { colors <- mclust.options("classPlotColors")[1:L]
        if(noise) 
        { colors <- unique(c("black", colors))[1:L] }
      }
      }
      else if(length(colors) == 1) colors <- rep(colors, L)
      if(length(symbols) < L) 
      { warning("more symbols needed to show classification ")
        symbols <- rep(16,L)
      }
      if(length(colors) < L) 
      { warning("more colors needed to show classification ")
        colors <- rep("black",L)
      }
  }
  
  what <- match.arg(what, choices = eval(formals(mclust2Dplot)$what))
  bad <- what == "classification" && is.null(classification)
  bad <- bad || (what == "uncertainty" && is.null(uncertainty))
  bad <- bad || (what == "error" && 
                   (is.null(classification) || is.null(truth)))
  if(bad) 
    stop("insufficient input for specified plot")

  switch(EXPR = what,
         "classification" = 
         {
           plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
           if(main) title("Classification")
           for(k in 1:L) 
              { I <- classification == U[k]
                points(data[I, 1], data[I, 2], pch = symbols[k], 
                       col = colors[k], cex = if(U[k] == "0") CEX/2 else CEX)
           }
         },
         "error" = 
         {
           ERRORS <- classError(classification, truth)$misclassified
           plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
           if(main) title("Classification Errors")
           CLASSES <- unique(as.character(truth))
           symOpen <- c(2, 0, 1, 5)
           symFill <- c(17, 15, 16, 18)
           good <- rep(TRUE,length(classification))
           good[ERRORS] <- FALSE
           if(L > 4) 
           { points(data[good, 1], data[good, 2], pch = 1, 
                    col = colors, cex = CEX)
             points(data[!good, 1], data[!good, 2], pch = 16, 
                    cex = CEX)
           }
           else 
           { for(k in 1:L) 
           { K <- truth == CLASSES[k]
             points(data[K, 1], data[K, 2], pch = symOpen[k], 
                    col = colors[k], cex = CEX)
             if(any(I <- (K & !good))) 
             { points(data[I, 1], data[I, 2], 
                      pch = symFill[k], cex = CEX)
             }
           }
           }
         },
         "uncertainty" =
         { 
           u <- (uncertainty - min(uncertainty))/
                (max(uncertainty) - min(uncertainty) + sqrt(.Machine$double.eps))
           b <- bubble(u, cex = CEX*c(0.3, 2), alpha = c(0.3, 0.9))
           cl <- sapply(classification, function(cl) which(cl == U))
           plot(data[, 1], data[, 2], pch = 19, 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", 
                cex = b$cex, 
                col = mapply(adjustcolor, 
                             col = colors[cl], 
                             alpha.f = b$alpha),
                ...)
           if(main) title("Uncertainty")
           fillEllipses <- FALSE
         },
         {  plot(data[, 1], data[, 2], type = "n", 
                 xlab = xlab, ylab = ylab, 
                 xlim = xlim, ylim = ylim, main = "", ...)
            if(main) title("Point Plot")
            points(data[, 1], data[, 2], pch = PCH, cex = CEX)
         }
  )
  if(haveParams && addEllipses) 
  { ## plot ellipsoids
    for(g in 1:G) 
      mvn2plot(mu = mu[,g], sigma = sigma[,,g], k = 15,
               fillEllipse = fillEllipses,
               col = if(fillEllipses) colors[g] else rep("grey30",3))
  }

  invisible()
}

# old version
mvn2plot <- function(mu, sigma, k = 15, alone = FALSE, 
                     col = rep("grey30",3), pch = 8, lty = c(1,2), lwd = c(1,1)) 
{
  p <- length(mu)
  if (p != 2) 
    stop("only two-dimensional case is available")
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
  if(alone) 
    { xymin <- apply(xy, 2, FUN = "min")
      xymax <- apply(xy, 2, FUN = "max")
      r <- ceiling(max(xymax - xymin)/2)
      xymid <- (xymin + xymax)/2
      plot(xy[, 1], xy[, 2], type = "n", xlab = "x", ylab = "y", 
           xlim = c(-r, r) + xymid[1], ylim = c(-r, r) + xymid[2])
  }
  l <- length(x)
  i <- 1:l
  for(k in 1:4) 
     { lines(xy[i,], col = col[1], lty = lty[1], lwd = lwd[1])
       i <- i + l
  }
  x <- s[1]
  y <- s[2]
  xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
  lines(xy[1:2,], col = col[2], lty = lty[2], lwd = lwd[2])
  lines(xy[3:4,], col = col[2], lty = lty[2], lwd = lwd[2])
  points(mu[1], mu[2], col = col[3], pch = pch)
  invisible()
}

mvn2plot <- function(mu, sigma, k = 15, alone = FALSE, 
                     fillEllipse = FALSE, alpha = 0.3, 
                     col = rep("grey30", 3), 
                     pch = 8, lty = c(1,2), lwd = c(1,1), ...)
{
  p <- length(mu)
  if(p != 2) 
    stop("only two-dimensional case is available")
  if(any(unique(dim(sigma)) != p)) 
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
  #
  if(alone) 
  { 
    xymin <- apply(xy, 2, FUN = "min")
    xymax <- apply(xy, 2, FUN = "max")
    r <- ceiling(max(xymax - xymin)/2)
    xymid <- (xymin + xymax)/2
    plot(xy[, 1], xy[, 2], type = "n", xlab = "x", ylab = "y", 
         xlim = c(-r, r) + xymid[1], ylim = c(-r, r) + xymid[2])
  }
  # draw ellipses
  if(fillEllipse)
  {
		col <- rep(col, 3)
    polygon(xy[chull(xy),], border = NA, 
            col = adjustcolor(col[1], alpha.f = alpha))
  } else
  {
    l <- length(x)
    i <- 1:l
    for(k in 1:4) 
    { 
      lines(xy[i,], col = col[1], lty = lty[1], lwd = lwd[1])
      i <- i + l
    }
  }
  # draw principal axes and centroid
  x <- s[1]
  y <- s[2]
  xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
  lines(xy[1:2,], col = col[2], lty = lty[2], lwd = lwd[2])
  lines(xy[3:4,], col = col[2], lty = lty[2], lwd = lwd[2])
  points(mu[1], mu[2], col = col[3], pch = pch)
  #
  invisible()
}

clPairs <- function (data, classification, symbols = NULL, colors = NULL, 
                     labels = dimnames(data)[[2]], cex.labels = 1.5, 
                     gap = 0.2, ...) 
{
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  if(missing(classification)) 
    classification <- rep(1, n)
  if(!is.factor(classification)) 
    classification <- as.factor(classification)
  l <- length(levels(classification))
  if(length(classification) != n)
    stop("classification variable must have the same length as nrows of data!")
  if(missing(symbols)) 
    { if(l == 1) 
        { symbols <- "." }
      if(l <= length(mclust.options("classPlotSymbols")))
        { symbols <- mclust.options("classPlotSymbols") }
      else { if(l <= 9) { symbols <- as.character(1:9) }
             else if(l <= 26) { symbols <- LETTERS[1:l] }
                  else symbols <- rep(16,l)
           }
  }
  if(length(symbols) == 1) symbols <- rep(symbols, l)
  if(length(symbols) < l) 
    { symbols <- rep(16, l)
      warning("more symbols needed")
  }
  if(is.null(colors)) 
    { if(l <= length(mclust.options("classPlotColors"))) 
      colors <- mclust.options("classPlotColors")[1:l]
  }
  if(length(colors) == 1) colors <- rep(colors, l)
  if(length(colors) < l) 
    { colors <- rep( "black", l)
      warning("more colors needed")
  }

  if(d > 2)
    { pairs(x = data, labels = labels, 
            pch = symbols[classification], 
            col = colors[classification], 
            gap = gap, 
            cex.labels = cex.labels,
            ...) }
  else if(d == 2)
    { plot(data, 
            pch = symbols[classification], 
            col = colors[classification], 
            ...) }
  
  invisible(list(d = d,
                 class = levels(classification), 
                 col = colors,
                 pch = symbols[seq(l)]))
}

clPairsLegend <- function(x, y, class, col, pch, box = TRUE, ...)
{
  
  usr <- par("usr")
  if(box & all(usr == c(0,1,0,1))) 
  {
    oldpar <- par(mar = rep(0.2, 4), no.readonly = TRUE)
    on.exit(par(oldpar))
    box(which = "plot", lwd = 0.8)
  }
  if(!all(usr == c(0,1,0,1)))
  {
    x <- x*(usr[2]-usr[1])+usr[1]
    y <- y*(usr[4]-usr[3])+usr[3]
  }

  dots <- list(...)
  dots$x <- x
  dots$y <- y
  dots$legend <- class
  dots$text.width <- max(strwidth(dots$title, units = "user"), 
                         strwidth(dots$legend, units = "user"))
  dots$col <- col
  dots$text.col <- col
  dots$pch <- pch
  dots$title.col <- par("fg")
  dots$title.adj <- 0.1
  dots$xpd <- NA
  do.call("legend", dots)
}

coordProj <- function(data, dimens = c(1,2), parameters = NULL, 
                      z = NULL, classification = NULL, 
                      truth = NULL, uncertainty = NULL, 
                      what = c("classification", "error", "uncertainty"), 
                      addEllipses = TRUE, 
                      fillEllipses = mclust.options("fillEllipses"),
                      symbols = NULL, colors = NULL, scale = FALSE, 
                      xlim = NULL, ylim = NULL, 
                      CEX = 1, PCH = ".", main = FALSE, ...)
{
  if(is.null(dimens)) dimens <- c(1, 2)
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) 
    { mu <- parameters$mean
      L <- ncol(mu)
      sigma <- parameters$variance$sigma
      haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any( is.na(sigma))
  }
  else haveParams <- FALSE
  data <- data[, dimens, drop = FALSE]
  if(dim(data)[2] != 2)
    stop("need two dimensions")
  if(is.null(xlim))
    xlim <- range(data[, 1])
  if(is.null(ylim))
    ylim <- range(data[, 2])
  if(scale) 
  {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
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
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
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
  if(!is.null(truth)) 
  {
    truth <- as.factor(truth)
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
  }
  if(!is.null(classification)) 
  {
    classification <- as.factor(classification)
    U <- levels(classification)
    L <- nlevels(classification)
    noise <- (U[1] == "0")
    if(is.null(symbols)) {
      if(L <= length(mclust.options("classPlotSymbols"))) 
      { symbols <- mclust.options("classPlotSymbols")[1:L]
        if(noise) 
          { symbols <- c(16,symbols)[1:L] }
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
      if(L <= length(mclust.options("classPlotColors"))) 
      {
        colors <- mclust.options("classPlotColors")[1:L]
        if(noise) 
          { colors <- unique(c("black", colors))[1:L] }
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
  if(length(what) > 1) what <- what[1]
  choices <- c("classification", "error", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) 
  {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "error" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
    badClass <- (what == "error" && (length(unique(classification)) != length(
      unique(truth))))
    if(badClass && !bad)
      warning("classification and truth differ in number of groups")
    bad <- bad && badClass
  } else 
  {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad) what <- "bad"
  
  switch(EXPR = what,
         "classification" = {
           plot(data[, 1], data[, 2], type = "n", 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", ...)
           if(main) {
             TITLE <- paste(paste(dimens, collapse = ","), 
                            "Coordinate Projection showing Classification")
             title(main = TITLE)
           }
           for(k in 1:L) {
             I <- classification == U[k]
             points(data[I, 1], data[I, 2], 
                    pch = symbols[k], col = colors[k], 
                    cex = if(U[k] == "0") CEX/3 else CEX)
           }
         },
         "error" = { 
           ERRORS <- classError(classification, truth)$misclassified
           plot(data[, 1], data[, 2], type = "n", 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", ...)
           if(main) {
             TITLE <- paste(paste(dimens, collapse = ","), 
                            "Coordinate Projection showing Errors")
             title(main = TITLE)
           }
           CLASSES <- levels(truth)
           symOpen <- symb2open(mclust.options("classPlotSymbols"))
           symFill <- symb2fill(mclust.options("classPlotSymbols"))
           good <- rep(TRUE, length(classification))
           good[ERRORS] <- FALSE
           if(L > length(symOpen)) 
           {
             points(data[good, 1], data[good, 2], pch = 1, col = colors, cex = CEX)
             points(data[!good, 1], data[!good, 2], pch = 16, cex = CEX)
           }
           else {
             for(k in 1:L) {
               K <- truth == CLASSES[k]
               if(any(I <- (K & good))) {
                 points(data[I, 1], data[I, 2], pch = symOpen[k], 
                        col = colors[k], cex = CEX)
               }
               if(any(I <- (K & !good))) {
                 points(data[I, 1], data[I, 2], cex = CEX,
                        pch = symFill[k], col = "black", bg = "black")
               }
             }
           }
         },
         "uncertainty" = { 
           u <- (uncertainty - min(uncertainty)) /
                (max(uncertainty) - min(uncertainty) + sqrt(.Machine$double.eps))
           b <- bubble(u, cex = CEX * c(0.3, 2), alpha = c(0.3, 0.9))
           cl <- sapply(classification, function(cl) which(cl == U))
           plot(data[, 1], data[, 2], pch = 19, main = "", 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim,
                cex = b$cex, 
                col = mapply(adjustcolor, col = colors[cl], alpha.f = b$alpha), 
                ...)
           if(main) 
           { 
             TITLE <- paste(paste(dimens, collapse = ","), 
                            "Coordinate Projection showing Uncertainty")
             title(main = TITLE)
           }
           fillEllipses <- FALSE
         },
         { plot(data[, 1], data[, 2], type = "n", 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", ...)
           if(main) 
             { TITLE <- paste(paste(dimens, collapse = ","), "Coordinate Projection")
               title(main = TITLE) }
           points(data[, 1], data[, 2], pch = PCH, cex = CEX)
         }
  )
  if(haveParams && addEllipses)
  { ## plot ellipsoids
    for(g in 1:G)
      mvn2plot(mu = mu[,g], sigma = sigma[,,g], k = 15,
               fillEllipse = fillEllipses,
               col = if(fillEllipses) colors[g] else rep("grey30",3))
  }

  invisible()
}

symb2open <- function(x)
{
  symb <- 0:18
  open <- c(0:14,0,1,2,5)
  open[sapply(x, function(x) which(symb == x))]
}

symb2fill <- function(x)
{
  symb <- 0:18
  fill <- c(15:17, 3:4, 23, 25, 7:9, 20, 11:18)
  fill[sapply(x, function(x) which(symb == x))]
}

# x <- c(16, 0, 17, 3, 15, 4, 1, 8, 2, 7, 5, 9, 6, 10, 11, 18, 12, 13, 14)
# plot(seq(x), rep(1,length(x)), pch = x, cex = 2, ylim = c(0.8, 3.8), yaxt = "n")
# points(seq(x), rep(2,length(x)), pch = symb2open(x), cex = 2)
# points(seq(x), rep(3,length(x)), pch = symb2fill(x), cex = 2, bg = "black")

randProj <- function(data, seeds = NULL, 
                     parameters = NULL, z = NULL,
                     classification = NULL, truth = NULL, 
                     uncertainty = NULL, 
                     what = c("classification", "error", "uncertainty"), 
                     quantiles = c(0.75, 0.95), 
                     addEllipses = TRUE, 
                     fillEllipses = mclust.options("fillEllipses"),
										 symbols = NULL, colors = NULL, scale = FALSE, 
                     xlim = NULL, ylim = NULL, 
                     xlab = NULL, ylab = NULL,
                     CEX = 1, PCH = ".", 
                     main = FALSE, ...)
{
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) 
  {
    mu <- parameters$mean
    L <- ncol(mu)
    sigma <- parameters$variance$sigma
    haveParams <- !is.null(mu) && !is.null(sigma) && 
                  !any(is.na(mu)) && !any(is.na(sigma))
  } else haveParams <- FALSE
  d <- ncol(data)
  if(haveParams) 
  {
    G <- ncol(mu)
    dimpar <- dim(sigma)
    if(length(dimpar) != 3) 
    {
      haveParams <- FALSE
      warning("covariance must be a 3D matrix")
    }
    if(G != dimpar[3]) 
    {
      haveParams <- FALSE
      warning("means and variance parameters are incompatible")
    }
    cho <- array(apply(sigma, 3, chol), c(d, d, G))
  }
  if(!is.null(truth)) 
  {
    truth <- as.factor(truth)
    if(is.null(classification)) 
    {
      classification <- truth
      truth <- NULL
    } else 
    {
      if(length(unique(truth)) != length(unique(classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) 
  {
    classification <- as.factor(classification)
    U <- levels(classification)
    L <- nlevels(classification)
    noise <- (U[1] == "0")
    if(is.null(symbols)) 
    {
      if(L <= length(mclust.options("classPlotSymbols"))) 
      { symbols <- mclust.options("classPlotSymbols")[1:L]
        if(noise)
        { symbols <- c(16,symbols)[1:L] }
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
    if(is.null(colors)) 
    { if(L <= length(mclust.options("classPlotColors"))) 
      { colors <- mclust.options("classPlotColors")[1:L]
        if(noise) 
          colors <- unique(c("black", colors))[1:L]
      }
    }
    else if(length(colors) == 1)
           colors <- rep(colors, L)
    if(length(symbols) < L) 
    {
      warning("more symbols needed to show classification ")
      symbols <- rep(16,L)
    }
    if (length(colors) < L) 
    {
      warning("more colors needed to show classification ")
      colors <- rep("black",L)
    }
  }
  if(is.null(xlab)) xlab <- "randProj1"
  if(is.null(ylab)) ylab <- "randProj2"
  what <- match.arg(what, choices = eval(formals(randProj)$what))
  bad <- what == "classification" && is.null(classification)
  bad <- bad || (what == "uncertainty" && is.null(uncertainty))
  bad <- bad || (what == "error" && (is.null(classification) || is.null(truth)))
  if(bad)
    stop("insufficient input for specified plot")
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  nullXlim <- is.null(xlim)
  nullYlim <- is.null(ylim)

  if(scale || length(seeds) > 1)
  { 
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if(scale) 
      par(pty = "s")
    if(length(seeds) > 1)
      par(ask = TRUE)
  }
  
  # if not provided get a seed at random
  if(length(seeds) == 0) 
  { 
    seeds <- as.numeric(Sys.time())
    seeds <- (seeds - floor(seeds))*1e8
  }
  for(seed in seeds) 
  {
    set.seed(seed)
    # B <- orth2(d)
    B <- randomOrthogonalMatrix(d, 2)
    dataProj <- as.matrix(data) %*% B
    if(dim(dataProj)[2] != 2)
      stop("need two dimensions")
    if(nullXlim)
      xlim <- range(dataProj[,1])
    if(nullYlim) 
      ylim <- range(dataProj[,2])
    if(scale) 
    {
      d <- diff(xlim) - diff(ylim)
      if(d > 0) 
      {
        ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
      } else 
      {
        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
      }
    }
    switch(what,
           "classification" = 
           {
             plot(dataProj[,1:2], type = "n", 
                  xlab = xlab, ylab = ylab, 
                  xlim = xlim, ylim = ylim, ...)
             for(k in 1:L) 
             {
               I <- classification == U[k]
               points(dataProj[I,1:2], pch = symbols[k], 
                      col = colors[k], cex = CEX)
             }
             if(main) 
             {
               TITLE <- paste("Random Projection showing Classification: seed = ", seed)
               title(TITLE)
             }
           },
           "error" = 
           {
             ERRORS <- classError(classification, truth)$misclassified
             plot(dataProj[, 1:2], type = "n", 
                  xlab = xlab, ylab = ylab, 
                  xlim = xlim, ylim = ylim, ...)
             if(main) 
             {
               TITLE <- paste("Random Projection showing Errors: seed = ", seed)
               title(TITLE)
             }
             CLASSES <- unique(as.character(truth))
             symOpen <- c(2, 0, 1, 5)
             symFill <- c(17, 15, 16, 18)
             good <- !ERRORS
             if(L > 4) 
             {
               points(dataProj[good, 1:2], pch = 1, col = colors, cex = CEX)
               points(dataProj[!good, 1:2], pch = 16, cex = CEX)
             } else 
             {
               for(k in 1:L) 
               {
                 K <- which(truth == CLASSES[k])
                 points(dataProj[K, 1:2], pch = symOpen[k], 
                        col = colors[k], cex = CEX)
                 if(any(I <- intersect(K, ERRORS)))
                   points(dataProj[I,1:2], pch = symFill[k], cex = CEX)
               }
             }
           },
           "uncertainty" = 
           {
             plot(dataProj[, 1:2], type = "n", 
                  xlab = xlab, ylab = ylab, 
                  xlim = xlim, ylim = ylim, main = "", ...)
             if(main) 
             {
               TITLE <- paste("Random Projection showing Uncertainty: seed = ", seed)
               title(TITLE)
             }
             breaks <- quantile(uncertainty, probs = sort(quantiles))
             I <- uncertainty <= breaks[1]
             points(dataProj[I, 1:2], 
                    pch = 16, col = "gray75", cex = 0.5 * CEX)
             I <- uncertainty <= breaks[2] & !I
             points(dataProj[I, 1:2],
                    pch = 16, col = "gray50", cex = 1 * CEX)
             I <- uncertainty > breaks[2] & !I
             points(dataProj[I, 1:2], 
                    pch = 16, col = "black", cex = 1.5 * CEX)
             fillEllipses <- FALSE
           },
           {
             plot(dataProj[, 1:2], type = "n", 
                  xlab = xlab, ylab = ylab, 
                  xlim = xlim, ylim = ylim, ...)
             if(main) 
             {
               TITLE <- paste("Random Projection: seed = ", seed)
               title(TITLE)
             }
             points(dataProj[, 1:2], pch = PCH, cex = CEX)
           }
    )

    muProj    <- crossprod(B, mu)
    sigmaProj <- array(apply(cho, 3, function(R) 
                             crossprod(R %*% B)), 
                       c(2, 2, G))
    
    if(haveParams && addEllipses) 
    { ## plot ellipsoids
      for(g in 1:G)
         mvn2plot(mu = muProj[,g], sigma = sigmaProj[,,g], 
                  k = 15, fillEllipse = fillEllipses,
                  col = if(fillEllipses) colors[g] else rep("grey30",3))
    }
  }

  invisible(list(basis = B,
                 data = dataProj, 
                 mu = muProj,
                 sigma = sigmaProj))
}

surfacePlot <- function(data, parameters, 
                        what = c("density", "uncertainty"), 
                        type = c("contour", "hdr", "image", "persp"), 
                        transformation = c("none", "log", "sqrt"), 
                        grid = 200, nlevels = 11, levels = NULL, 
                        prob = c(0.25, 0.5, 0.75),
                        col = gray(0.7),
                        col.palette = function(...) hcl.colors(..., "blues", rev = TRUE),
                        hdr.palette = blue2grey.colors,
                        xlim = NULL, ylim = NULL, 
                        xlab = NULL, ylab = NULL, 
                        main = FALSE, scale = FALSE, 
                        swapAxes = FALSE,
                        verbose = FALSE,  ...) 
{
  data <- as.matrix(data)
  if(dim(data)[2] != 2) 
    stop("data must be two dimensional")
  if(any(type == "level")) 
    type[type == "level"] <- "hdr" # TODO: to be removed
  type <- match.arg(type, choices = eval(formals(surfacePlot)$type))
  what <- match.arg(what, choices = eval(formals(surfacePlot)$what))
  transformation <- match.arg(transformation, 
                              choices = eval(formals(surfacePlot)$transformation))
  #
  
  densNuncer <- function(modelName, data, parameters) 
  {
    if(is.null(parameters$variance$cholsigma))
    { parameters$variance$cholsigma <- parameters$variance$sigma
      G <- dim(parameters$variance$sigma)[3]
      for(k in 1:G)
        parameters$variance$cholsigma[,,k] <- chol(parameters$variance$sigma[,,k])
    }
    cden <- cdensVVV(data = data, parameters = parameters, logarithm = TRUE)
    pro <- parameters$pro
    if(!is.null(parameters$Vinv))
      pro <- pro[-length(pro)]
    z <- sweep(cden, MARGIN = 2, FUN = "+", STATS = log(pro))
    logden <- apply(z, 1, logsumexp)
    z <- sweep(z, MARGIN = 1, FUN = "-", STATS = logden)
    z <- exp(z)
    data.frame(density = exp(logden),
               uncertainty = 1 - apply(z, 1, max))
  }
  
  pro <- parameters$pro
  mu <- parameters$mean
  sigma <- parameters$variance$sigma
  haveParams <- (!is.null(mu) && !is.null(sigma) && !is.null(pro) && 
                 !any(is.na(mu)) && !any(is.na(sigma)) && !(any(is.na(pro))))
  if(haveParams) 
    { G <- ncol(mu)
      dimpar <- dim(sigma)
      if(length(dimpar) != 3) 
        { haveParams <- FALSE
          warning("covariance must be a 3D matrix")
      }
      if(G != dimpar[3]) 
        { haveParams <- FALSE
          warning("means and variance parameters are incompatible")
      }
      mu <- array(mu, c(2, G))
      sigma <- array(sigma, c(2, 2, G))
  }
  else stop("need parameters to compute density")
  
  if(swapAxes) 
    { if(haveParams) 
        { parameters$pro <- pro[2:1]
          parameters$mean <- mu[2:1,]
          parameters$variance$sigma <- sigma[2:1, 2:1,]
      }
      data <- data[, 2:1]
  }
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  if(is.null(xlim)) xlim <- range(data[, 1])
  if(is.null(ylim)) ylim <- range(data[, 2])
  if(scale)
    { par(pty = "s")
      d <- diff(xlim) - diff(ylim)
      if(d > 0) 
        { ylim <- c(ylim[1] - d/2, ylim[2] + d/2) }
      else 
        { xlim <- c(xlim[1] + d/2, xlim[2] - d/2) }
  }
  
  dnames <- dimnames(data)[[2]]
  if(is.null(xlab)) 
    { xlab <- if(is.null(dnames)) "" else dnames[1] }
  if(is.null(ylab)) 
    { ylab <- if(is.null(dnames)) "" else dnames[2] }
  
  if(length(grid) == 1) 
    grid <- c(grid, grid)
  x <- grid1(n = grid[1], range = xlim, edge = TRUE)
  y <- grid1(n = grid[2], range = ylim, edge = TRUE)
  xy <- grid2(x, y)
  
  if(verbose) 
    message("computing density and uncertainty over grid ...")
  Z <- densNuncer(modelName = "VVV", data = xy, parameters = parameters)
  lx <- length(x)
  ly <- length(y)
  #
  switch(what, 
         "density"     = { zz <- matrix(Z$density, lx, ly)
                           title2 <- "Density" }, 
         "uncertainty" = { zz <- matrix(Z$uncertainty, lx, ly)
                           title2 <- "Uncertainty" }, 
         stop("what improperly specified"))
  #
  switch(transformation, 
         "none" = { title1 <- "" }, 
         "log"  = { zz <- log(zz)
                    title1 <- "log" }, 
         "sqrt" = { zz <- sqrt(zz)
                    title1 <- "sqrt" }, 
         stop("transformation improperly specified"))
  #
  switch(type,
         "contour" = {
           title3 <- "Contour"
           if(is.null(levels)) levels <- pretty(zz, nlevels)
           contour(x = x, y = y, z = zz, levels = levels, 
                   xlab = xlab, ylab = ylab, 
                   col = col, main = "", ...)
         },
         "hdr" = {
           title3 <- "HDR level"
           z <- densNuncer(modelName = "VVV", data = data, 
                           parameters = parameters)$density
           levels <- c(sort(hdrlevels(z, prob)), 1.1*max(z))
           plot(x, y, type = "n",
                xlab = xlab, ylab = ylab, ...)
           fargs <- formals(".filled.contour")
           dargs <- c(list(x = x, y = y, z = zz, 
                           levels = levels,
                           col = hdr.palette(length(levels))), 
                      args)
           dargs <- dargs[names(dargs) %in% names(fargs)]
           fargs[names(dargs)] <- dargs
           do.call(".filled.contour", fargs)
         }, 
         "image" = {
           title3 <- "Image"
           col <- col.palette(nlevels)
           if(length(col) == 1)
             { if(!is.null(levels)) 
                 nlevels <- length(levels)
               col <- mapply(adjustcolor, col = col, 
                             alpha.f = seq(0.1, 1, length = nlevels))
           }
           image(x = x, y = y, z = zz, xlab = xlab, ylab = ylab, 
                 col = col, main = "", ...)
         }, 
         "persp" = {
           title3 <- "Perspective"
           dots <- list(...)
           if(is.null(dots$theta))  dots$theta <- -30
           if(is.null(dots$phi))    dots$phi <- 20
           if(is.null(dots$expand)) dots$expand <- 0.6
           p3d <- do.call("persp", 
                          c(list(x = x, y = y, z = zz, border = NA,
                                 xlab = xlab, ylab = ylab, 
                                 col = col,
                                 zlab = "Density", main = ""), dots))
           ii <- floor(seq(1, length(y), length.out = 2*nlevels))
           for(i in ii[-c(1,length(ii))])
             lines(trans3d(x, y[i], zz[,i], pmat = p3d))
           ii <- floor(seq(1, length(x), length.out = 2*nlevels))
           for(i in ii[-c(1,length(ii))])
             lines(trans3d(x[i], y, zz[i,], pmat = p3d))
         }
  )
  if(main) 
    { TITLE <- paste(c(title1, title2, title3, "Plot"), collapse = " ")
      title(TITLE) }

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
  plot(uncer[ord], type = "n", xaxt = "n", ylim = c(-(M/32), M), 
       ylab = "uncertainty", 
       xlab = "observations in order of increasing uncertainty")
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


blue2grey.colors <- function(n) 
{
  # manually selected
  basecol <- c("#E6E6E6", "#bcc9d1", "#6c7f97", "#3e5264")
  # selected using colorspace::sequential_hcl(5, palette = "blues2")
  # basecol <- c("#023FA5", "#6A76B2", "#A1A6C8", "#CBCDD9", "#E2E2E2")
  palette <- grDevices::colorRampPalette(basecol, space = "Lab")
  palette(n)
}

bubble <- function(x, cex = c(0.2, 3), alpha = c(0.1, 1)) 
{
  x <- as.vector(x)
  cex <- cex[!is.na(cex)]
  alpha <- alpha[!is.na(alpha)]
  x <- (x - min(x))/(max(x) - min(x) + sqrt(.Machine$double.eps))
  n <- length(x)
  r <- sqrt(x/pi)
  r <- (r - min(r, na.rm = TRUE))/
       (max(r, na.rm = TRUE) - min(r, na.rm = TRUE) + sqrt(.Machine$double.eps))
  cex <- r * diff(range(cex)) + min(cex)
  alpha <- x * diff(range(alpha)) + min(alpha)
  return(list(cex = cex, alpha = alpha))
}

grid1 <- function(n, range = c(0, 1), edge = TRUE) 
{
  if(any(n < 0 | round(n) != n)) 
    stop("n must be nonpositive and integer")
  G <- rep(0, n)
  if(edge) 
  {
    G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n-1))
  } else 
  {
    lj <- abs(diff(range))
    incr <- lj/(2 * n)
    G <- seq(from = min(range) + incr, to = max(range) - incr, by = 2 * incr)
  }
  return(G)
}

grid2 <- function(x, y) 
{
  lx <- length(x)
  ly <- length(y)
  xy <- matrix(0, nrow = lx * ly, ncol = 2)
  l <- 0
  for(j in 1:ly) 
  {
    for(i in 1:lx) 
    {
      l <- l + 1
      xy[l,] <- c(x[i], y[j])
    }
  }
  return(xy)
}

vpoints <- function(x, y, col, cex = 1, ...)
{
  xy <- xy.coords(x, y)
  symbols(xy$x, xy$y, add = TRUE, inches = 0.2*cex,
          fg = if(missing(col)) par("col") else col,
          rectangles = matrix(c(0,1), nrow = length(xy$x), ncol = 2, byrow = TRUE), 
          ...)
}
