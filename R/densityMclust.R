densityMclust <- function(data, ...) 
{
  mc <- match.call()
  if(length(dim(data)) > 1) varname <- colnames(data)
  else                      varname <- deparse(substitute(data))
  obj <- Mclust(data, ...)
  obj$call <- mc
  d <- dens(modelName = obj$modelName, 
            data = data, 
            parameters = obj$parameters)
  obj <- c(obj,
           list(varname = varname,
                range = if(obj$d > 1) apply(data, 2, range) else range(data),
                density = d))
  class(obj) <- c("densityMclust", "Mclust")
  return(obj)
}

predict.densityMclust <- function(object, newdata, ...)
{
  if(!inherits(object, "densityMclust")) 
    stop("object not of class \"densityMclust\"")
  if(missing(newdata))
    { newdata <- eval.parent(object$call$data) }
  d <- dens(modelName = object$modelName, 
            data = newdata, 
            parameters = object$parameters)
  return(d)
}

plot.densityMclust <- function(x, data = NULL, what = c("density", "BIC", "diagnostic"), ...) 
{
  object <- x # Argh.  Really want to use object anyway
  what <- match.arg(what)
         
  if(what == "density")
    { 
      if(object$d == 1)      plotDensityMclust1(object, data = data, ...)
      else if(object$d == 2) plotDensityMclust2(object, data = data, ...)
           else              plotDensityMclustd(object, data = data, ...)
    }

  if(what == "BIC")
    { 
      # this add right axis for bic diff
      # oldpar <- par(no.readonly = TRUE)
      # on.exit(par(oldpar))
      # mar <- oldpar$mar
      # mar[4] <- max(mar[4],3)
      # par(mar = mar)
      # plot.mclustBIC(object$BIC, ...)
      # yaxp <- par("yaxp")
      # bicdiff <- seq(0, yaxp[1] - object$bic, length = 100)
      # bicdiff <- pretty(bicdiff, yaxp[3]+1)
      # axis(4, at = object$bic+bicdiff, labels = signif(bicdiff,2))
      plot.mclustBIC(object$BIC, ...)
    }
    
  if(what == "diagnostic")
    { if(missing(data))
         data <- eval.parent(object$call$data)
      densityMclust.diagnostic(object, data, what = c("cdf", "qq"), ...) 
    }
    
  invisible()
}

plotDensityMclust1 <- function(x, data = NULL, hist.col = "lightgrey", hist.border = "grey", breaks = "Sturges", ...) 
{
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$data <- mc$hist.col <- mc$hist.border <- mc$breaks <- NULL
  xlab <- mc$xlab
  if(is.null(xlab)) 
     xlab <- x$varname[1]
  ylab <- mc$ylab
  if(is.null(ylab)) 
     ylab <- "Density"
  #
  xrange <- c(min(object$range) - diff(object$range)/10,
              max(object$range) + diff(object$range)/10)
  xlim <- mc$xlim
  if(!is.null(xlim)) 
     xrange <- range(xlim)
  ylim <- mc$ylim
  #
  eval.points <- seq(from = xrange[1], to = xrange[2], length = 1000)
  d <- predict.densityMclust(object, eval.points)
  #
  if(!is.null(data)) 
    { h <- hist(data, breaks = breaks, plot = FALSE)
      plot(h, freq = FALSE, col = hist.col, border = hist.border, main = "",
           xlim = range(h$breaks, xrange), 
           ylim = range(0, ylim, h$density, max(d)+diff(range(d))*0.1),
           xlab = xlab, ylab = ylab)
      mc[[1]] <- as.name("lines")
      mc$x <- eval.points; mc$y <- d; mc$type <- "l"
      eval(mc, parent.frame())
    }
  else
    { mc[[1]] <- as.name("plot")
      mc$x <- eval.points; mc$y <- d
      mc$type <- "l"; mc$xlim <- xlim
      mc$ylim <- range(0, ylim, max(d)+diff(range(d))*0.1)
      mc$ylab <- ylab; mc$xlab <- xlab
      eval(mc, parent.frame())
    }
  invisible()
}

plotDensityMclust2 <- function(x, data = NULL, col = grey(0.6), nlevels = 11, levels = NULL, points.col = 1, pch = 1, ...) 
{
# This function call surfacePlot() with a suitable modification of arguments
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$points.col <- mc$pch <- NULL
  mc$nlevels <- nlevels; mc$levels <- levels
  if(!is.null(mc$type))
    if(mc$type == "image" & (length(col) < 2)) col <- NULL
  mc$col <- col

  if(is.null(data)) 
    { addPoints <- FALSE
      mc$data <- object$range 
    }
  else
    { addPoints <- TRUE }

  # set mixture parameters
  par <- object$parameters
  # these parameters should be missing 
  par$variance$cholSigma <- par$Sigma <- par$Vinv <- NULL
  if(is.null(par$pro)) par$pro <- 1  # LS: bug?
  #par$variance$d <- 2 # LS: bug?
  par$variance$cholsigma <- par$variance$sigma
  for(k in seq(par$variance$G))
     { par$variance$cholsigma[,,k] <- chol(par$variance$sigma[,,k]) }
  mc$parameters <- par
  # now surfacePlot() is called
  mc[[1]] <- as.name("surfacePlot")
  out <- eval(mc, parent.frame())
  if(addPoints)
    points(data, col = points.col, pch = pch)
  #
  invisible(out)
}

plotDensityMclustd <- function(x, data = NULL, col = grey(0.6), nlevels = 11, levels = NULL, points.col = 1, pch = 1, gap = 0.2, ...) 
{
# This function call surfacePlot() with a suitable modification of arguments
    
  object <- x # Argh.  Really want to use object anyway
  mc <- match.call(expand.dots = TRUE)
  mc$x <- mc$points.col <- mc$pch <- mc$gap <- NULL
  mc$nlevels <- nlevels; mc$levels <- levels
  if(!is.null(mc$type))
    if(mc$type == "image" & (length(col) < 2)) col <- NULL
  mc$col <- col

  if(is.null(data)) 
    { addPoints <- FALSE
      mc$data <- object$range }
  else
    { data <- as.matrix(data)
      addPoints <- TRUE 
      object$range <- apply(data, 2, range) 
      object$varname <- colnames(data) }

  nc <- object$d
  oldpar <- par(mfrow = c(nc, nc), 
                mar = rep(c(gap,gap/2),each=2), 
                oma = c(4, 4, 4, 4),
                no.readonly = TRUE)
  on.exit(par(oldpar))

  for(i in seq(nc))
     { for(j in seq(nc)) 
          { if(i == j) 
              { plot(0,0,type="n",xlab="",ylab="",axes=FALSE)
                text(0,0, object$varname[i], cex=1.5, adj=0.5)
                box()
              } 
            else 
              { # set mixture parameters
                par <- object$parameters
                if(is.null(par$pro)) par$pro <- 1
                par$mean <- par$mean[c(j,i),,drop=FALSE]
                par$Vinv <- NULL
                par$variance$d <- 2
                sigma <- array(dim = c(2, 2, par$variance$G))
                for(g in seq(par$variance$G))
                   sigma[,,g] <- par$variance$sigma[c(j,i),c(j,i),g]
                par$variance$sigma <- sigma
                par$variance$Sigma <- NULL
                par$variance$cholSigma <- NULL
                par$variance$cholsigma <- NULL
                mc$parameters <- par
                mc$data <- object$range[,c(j,i)]
                mc$axes <- FALSE
                mc[[1]] <- as.name("surfacePlot")
                out <- eval(mc, parent.frame())
                box()
                if(addPoints)
                  points(data[,c(j,i)], col = points.col, pch = pch)

              }
            if(i == 1 && (!(j%%2))) axis(3)
            if(i == nc && (j%%2))   axis(1)
            if(j == 1 && (!(i%%2))) axis(2)
            if(j == nc && (i%%2))   axis(4)
          }
     }
  #
  invisible(out)
}

dens <- function(modelName, data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- .mclust$warn
  aux <- list(...)
  cden <- cdens(modelName = modelName, data = data,
                logarithm = TRUE, parameters = parameters, warn = warn)
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  G <- if(oneD) 
         { length(parameters$mean) }
       else 
         { ncol(as.matrix(parameters$mean)) }
  if(G > 1) 
    { pro <- parameters$pro
      if(is.null(pro))
         stop("mixing proportions must be supplied")
      noise <- length(pro) == (G + 1)
      if(is.null(pro)) 
        stop("mixing proportions must be supplied")
      if(noise) 
        { proN <- pro[length(pro)]
          pro <- pro[-length(pro)]
        }
      if(any(proz <- pro == 0)) 
        { pro <- pro[!proz]
          cden <- cden[, !proz, drop = FALSE]
        }
      cden <- sweep(cden, 2, FUN = "+", STATS = log(pro))
  }
  maxlog <- apply(cden, 1, max)
  cden <- sweep(cden, 1, FUN = "-", STATS = maxlog)
  den <- logb(apply(exp(cden), 1, sum)) + maxlog
  if(!logarithm) den <- exp(den)
  den
}

cdens <- function(modelName, data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)
  funcName <- paste("cdens", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}



densityMclust.diagnostic <- function(object, data, what = c("cdf", "qq"), col = c(1,3), lwd = c(2,2), lty = c(1,2), legend = TRUE, ...)
{
# Diagnostic plots for density estimation 
# (only available for the one-dimensional case)
# 
# Arguments:
# object = a 'densityMclust' object
# data = the data vector
# what = type of diagnostic plot:
# cdf = the fitted distribution function vs the empirical distribution function;
# qq = the fitted distribution function evaluated over the observed points vs 
#      the quantile from a uniform distribution.
#
# Reference: 
# Loader C. (1999), Local Regression and Likelihood. New York, Springer, 
#   pp. 87-90)

  if(!any(class(object) == "densityMclust"))
    { stop("first argument must be an object of class 'densityMclust'") }
  if(object$d > 1)
    { warning("only available for one-dimensional data") }  
  if(missing(data))
    { stop("data must be provided") }
  what <- match.arg(what, c("cdf", "qq"), several.ok = TRUE)
                    
  data <- as.numeric(data)
  n <- length(data)
  cdf <- cdfMclust(object, data = data, ...)

  oldpar <- par(no.readonly = TRUE)
  if(length(what) > 1) 
    { par(ask = TRUE)
      on.exit(par(oldpar)) }

  if(any(what == "cdf"))
    { # Fitted CDF vs Emprical CDF    
      empcdf <- ecdf(data)
      plot(empcdf, do.points = FALSE, col = col[2], lwd = lwd[2], lty = lty[2],
           xlab = object$varname, ylab = "Cumulative Distribution Function",
           main = "CDF plot")
      lines(cdf, col = col[1], lwd = lwd[1], lty = lty[1])
      if(legend)
        { legend("bottomright", legend = c("Est.CDF", "Emp.CDF"), 
                 ncol = 2, inset = 0.05, cex = 0.8,
                 col = col, lwd = lwd, lty = lty) }
    }

 if(any(what == "qq"))
   { # Q-Q plot
     inv.cdf <- approx(cdf$y, cdf$x, ppoints(n))$y
     plot(inv.cdf, sort(data),
          xlab = "Quantiles from estimated density", 
          ylab = "Sample Quantiles", 
          main = "Q-Q plot")
     with(list(y = sort(data), x = inv.cdf),
          { i <- (y > quantile(y, 0.25) & y < quantile(y, 0.75))
            abline(lm(y ~ x, subset = i), lty = 2) 
          })
     # P-P plot
     # cdf <- cdfMclust(object, data, ...)
     # plot(seq(1,n)/(n+1), cdf$y, xlab = "Uniform quantiles", 
     #    ylab = "Cumulative Distribution Function",
     #      main = "Diagnostic: P-P plot")
     # abline(0, 1, lty = 2)
   }

  invisible()
} 

cdfMclust <- function(object, data, ngrid = 100, ...)
{
# Cumulative Density Function
# (only available for the one-dimensional case)
#
# Returns the estimated CDF evaluated at points given by the optional
# argument data. If not provided, a regular grid of ngrid points is used. 
#
# Arguments:
# object = a 'densityMclust' object
# data = the data vector
# ngrid = the length of rectangular grid 

  if(!any(class(object) == "densityMclust"))
    { stop("first argument must be an object of class 'densityMclust'") }

  if(missing(data))
    { eval.points <- seq(min(object$range) - diff(object$range)/10,
                         max(object$range) + diff(object$range)/10, 
                         length = ngrid) }
  else
    { eval.points <- sort(as.vector(data))
      ngrid <- length(eval.points) }

  G <- object$parameters$variance$G
  pro <- object$parameters$pro
  mean <- object$parameters$mean
  var <- object$parameters$variance$sigmasq
  if(length(var) < G) var <- rep(var, G)
  
  cdf <- rep(0, ngrid)
  for(k in seq(G))
     { cdf <- cdf + pro[k]*pnorm(eval.points, mean[k], sqrt(var[k])) }
  
  # old: integral of the density function is approximated by adaptive
  # quadrature using the R function integrate().
  # f <- function(x) { predict.densityMclust(object, x) }
  # cdf <- rep(NA, ngrid)
  # for(i in 1:ngrid)
  #    { cdf[i] <- integrate(f, stop.on.error = FALSE,
  #                          rel.tol = sqrt(.Machine$double.eps),
  #                          lower = -Inf,
  #                          upper = eval.points[i])$value }
  
  out <- list(x = eval.points, y = cdf)    
  return(out)
}

