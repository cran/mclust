##
## Resampling methods
##

#
# Bootstrap Likelihood Ratio Test
#

mclustBootstrapLRT <- function(data, modelName = NULL, 
                               nboot = 999, level = 0.05, maxG = NULL, 
                               verbose = interactive(), ...)
{
  if(is.null(modelName))
    stop("A 'modelName' must be provided. Please see help(mclustModelNames) which describes the available models.")
  modelName <- modelName[1]
  if(is.null(maxG)) G <- seq.int(1, 9)
  else { maxG <- as.numeric(maxG); G <- seq.int(1, maxG+1) }
  Bic <- mclustBIC(data, G = G, modelNames = modelName, 
                   warn = FALSE, verbose = FALSE, ...)
  if(!(modelName %in% attr(Bic, "modelNames")))
    stop("'modelName' not compatibile with data. Please see help(mclustModelNames) which describes the available models.")
  if(all(is.na(Bic)))
    stop(paste("no model", modelName, "can be fitted."))
  # select only models that can be fit
  G <- which(!is.na(Bic[, attr(Bic, "modelNames") == modelName]))

  if(verbose) 
    { cat("bootstrapping LRTS ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = (max(G)-1)*nboot, style = 3)
      on.exit(close(pbar))
  }

  obsLRTS <- p.value <- vector("numeric", length = max(G)-1)
  bootLRTS <- matrix(as.double(NA), nrow = nboot, ncol = max(G)-1)
  g <- 0; continue <- TRUE
  while(g < (max(G)-1) & continue)
  { g <- g + 1
    # fit model under H0
    Mod0 <- summary(Bic, data, G = g, modelNames = modelName)
    # fit model under H1
    Mod1 <- summary(Bic, data, G = g+1, modelNames = modelName)
    # observed LRTS
    obsLRTS[g] <- 2*(Mod1$loglik - Mod0$loglik)
    # bootstrap
    b <- 0
    while(b < nboot)
    { 
			b <- b + 1
      # generate 'parametric' bootstrap sample under H0
      bootSample <- sim(Mod0$modelName, Mod0$parameters, n = Mod0$n)
      # fit model under H0
      bootMod0 <- em(data = bootSample[,-1], modelName = Mod0$modelName, 
                     parameters = Mod0$parameters, warn = FALSE, ...)
      # fit model under H1
      bootMod1 <- em(data = bootSample[,-1], modelName = Mod1$modelName, 
                     parameters = Mod1$parameters, warn = FALSE, ...)
      # compute bootstrap LRT
      LRTS <- 2*(bootMod1$loglik - bootMod0$loglik)
      if(is.na(LRTS)) { b <- b - 1; next() }
      bootLRTS[b,g] <- LRTS 
      if(verbose) 
        setTxtProgressBar(pbar, (g-1)*nboot+b)
    }
    p.value[g] <- (1 + sum(bootLRTS[,g] >= obsLRTS[g]))/(nboot+1)
    # check if not-significant when no maxG is provided
    if(is.null(maxG) & p.value[g] > level) 
      { continue <- FALSE
        if(verbose) 
          setTxtProgressBar(pbar, (max(G)-1)*nboot) 
      }
  }

  out <- list(G = 1:g, 
              modelName = modelName,
              obs = obsLRTS[1:g],
              boot = bootLRTS[,1:g,drop=FALSE],
              p.value = p.value[1:g])
  class(out) <- "mclustBootstrapLRT"
  return(out)
}

print.mclustBootstrapLRT <- function(x, ...)
{
  txt <- paste(rep("-", min(61, getOption("width"))), collapse = "")
  catwrap(txt)
  catwrap("Bootstrap sequential LRT for the number of mixture components") 
  catwrap(txt)
  cat(formatC("Model", flag = "-", width = 12), "=", x$modelName, "\n")
  cat(formatC("Replications", flag = "-", width = 12), "=", nrow(x$boot), "\n")
  df <- data.frame(x$obs, x$p.value)
  colnames(df) <- c("LRTS", "bootstrap p-value")
  rownames(df) <- formatC(paste(x$G, "vs", x$G+1), flag = "-", width = 8)
  print(df, ...)
}

plot.mclustBootstrapLRT <- function(x, G = 1, hist.col = "grey", hist.border = "lightgrey", breaks = "Scott", col = "forestgreen", lwd = 2, lty = 3, main = NULL, ...) 
{
  if(!any(G == x$G))
    { warning(paste("bootstrap LRT not available for G =", G)) 
      return() }
  G <- as.numeric(G)[1]
  h <- hist(x$boot[,G], breaks = breaks, plot = FALSE)
  xlim <- range(h$breaks, x$boot[,G], x$obs[G], na.rm = TRUE)
  xlim <- extendrange(xlim, f = 0.05)
  plot(h, xlab = "LRTS", freq = FALSE, xlim = xlim,
       col = hist.col, border = hist.border, main = NULL)
  box()
  abline(v = x$obs[G], lty = lty, lwd = lwd, col = col)
  if(is.null(main) | is.character(main))
    { if(is.null(main)) main <- paste("Bootstrap LRT for model", x$modelName, 
                                      "with", G, "vs", G+1, "components")
      title(main = main, cex.main = 1) }
  invisible()
}

#
# Bootstrap inference (standard errors and percentile confidence intervals) 
#

MclustBootstrap <- function(object, nboot = 999, type = c("bs", "wlbs", "pb", "jk"),
                            max.nonfit = 10*nboot, verbose = interactive(), ...)
{
  
  if(!any(class(object) %in% c("Mclust", "densityMclust")))
    stop("object must be of class 'Mclust' or 'densityMclust'")
  
  if(any(type %in% c("nonpara", "wlb")))
    { type <- gsub("nonpara", "bs", type)
      type <- gsub("wlb", "wlbs", type)
      warning("resampling type converted to \"", type, "\"")
    }
  type <- match.arg(type, choices = eval(formals(MclustBootstrap)$type))
  
  # data <- object$data
  n <- object$n
  d <- object$d
  G <- object$G
  if(type == "jk") nboot <- n
  varnames <- rownames(object$parameters$mean)
  # model parameters
  par <- summary(object)[c("pro", "mean", "variance")]
  if(d == 1)
    { par$mean <- array(par$mean, dim = c(d, G))
      par$variance <- array(par$variance, dim = c(d, d, G)) }
  # bootstrapped parameters 
  pro.boot  <- array(NA, c(nboot,G), 
                     dimnames = list(NULL, seq.int(G)))
  mean.boot <- array(NA, c(nboot,d,G), 
                     dimnames = list(NULL, varnames, seq.int(G)))
  var.boot  <- array(NA, c(nboot,d,d,G),
                     dimnames = list(NULL, varnames, varnames, seq.int(G)))

  if(verbose) 
    { cat("resampling ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = nboot, style = 3) 
      on.exit(close(pbar))
    }
  b <- nonfit <- 0
  while(b < nboot & nonfit < max.nonfit)
  { 
    b <- b + 1
    obj <- object
    switch(type, 
           "bs" = 
           { idx <- sample(seq_len(n), size = n, replace = TRUE)
             obj$data <- object$data[idx,]
             obj$z <- object$z[idx,]
             obj$warn <- FALSE
             mod.boot <- try(do.call("me", obj), silent = TRUE)
           },
           "wlbs" = 
           { w <- rexp(n)
             # w <- w/mean(w)
             w <- w/max(w)
             mod.boot <- try(do.call("me.weighted", 
                                     c(list(weights = w, warn = FALSE), obj)),
                             silent = TRUE)
           },
           "pb" = 
           { obj$data <- do.call("sim", object)[,-1,drop=FALSE]
             obj$z <- predict(obj)$z
             obj$warn <- FALSE
             mod.boot <- try(do.call("me", obj), silent = TRUE)
           },
           "jk" =
           { idx <- seq_len(n)[-b]
             obj$data <- object$data[idx,]
             obj$z <- object$z[idx,]
             obj$warn <- FALSE
             mod.boot <- try(do.call("me", obj), silent = TRUE)
           }
    )

    # check model convergence
    if(inherits(mod.boot, "try-error"))
      { if(type != "jk") b <- b - 1
        nonfit <- nonfit + 1
        next() }
    if(is.na(mod.boot$loglik))
      { if(type != "jk") b <- b - 1
        nonfit <- nonfit + 1
        next() }
    
    if(type == "jk")
      { # pseudovalues ...
        # pro.boot[b,]   <- n*par$pro - (n-1)*mod.boot$parameters$pro
        # mean.boot[b,,] <- n*par$mean - (n-1)*mod.boot$parameters$mean
        # var.boot[b,,,] <- n*par$variance - (n-1)*mod.boot$parameters$variance$sigma
        pro.boot[b,]   <- mod.boot$parameters$pro
        mean.boot[b,,] <- mod.boot$parameters$mean
        var.boot[b,,,] <- mod.boot$parameters$variance$sigma
    } else 
      { # bootstrap values
        pro.boot[b,]   <- mod.boot$parameters$pro
        mean.boot[b,,] <- mod.boot$parameters$mean
        var.boot[b,,,] <- mod.boot$parameters$variance$sigma
    }

    if(verbose) setTxtProgressBar(pbar, b)
  }

  out <- list(G = G, 
              modelName = object$modelName, 
              parameters = par,
              nboot = nboot, 
              type = type,
              nonfit = nonfit,
              pro = pro.boot, 
              mean = mean.boot, 
              variance = var.boot)
  class(out) <- "MclustBootstrap"
  return(out)
}

print.MclustBootstrap <- function(x,  digits = getOption("digits"), ...)
{
  cat("\'", class(x)[1], "\' object:\n", sep = "")
  str(x, max.level = 1, give.attr = FALSE, strict.width = "wrap")
  invisible()
}

summary.MclustBootstrap <- function(object, what = c("se", "ci", "ave"), conf.level = 0.95, ...)
{
  what <- match.arg(what, choices = eval(formals(summary.MclustBootstrap)$what))
  dims <- dim(object$mean)
  # varnames <- dimnames(object$mean)[[2]]
  nboot <- dims[1]
  d <- dims[2]
  G <- dims[3]

  switch(what,
    "se" = { out <- list(pro  = apply(object$pro, 2, sd, na.rm=TRUE),
                         mean = apply(object$mean, c(2,3), sd, na.rm=TRUE),
                         variance = apply(object$variance, c(2,3,4), sd, na.rm=TRUE))
             if(object$type == "jk")
                out <- lapply(out, function(x) 
                              sqrt(x^2*(nboot-object$nonfit-1)^2/(nboot-object$nonfit)))
           },
    "ave" = { out <- list(pro  = apply(object$pro, 2, mean, na.rm=TRUE),
                          mean = apply(object$mean, c(2,3), mean, na.rm=TRUE),
                          variance = apply(object$variance, c(2,3,4), mean, na.rm=TRUE))
           },
    "ci" = { levels <- c((1-conf.level)/2, (1 + conf.level)/2)
             if(object$type == "jk")
             { # bias-corrected ci based on normal-approximation
               ave <- list(pro  = apply(object$pro, 2, mean, na.rm=TRUE),
                           mean = apply(object$mean, c(2,3), mean, na.rm=TRUE),
                           variance = t(sapply(seq.int(d), function(j)
                                        apply(object$variance[,j,j,], 2, mean, na.rm=TRUE),
                                        simplify = "array")))
               se <- list(pro  = apply(object$pro, 2, sd, na.rm=TRUE),
                          mean = apply(object$mean, c(2,3), sd, na.rm=TRUE),
                          variance  = t(sapply(seq.int(d), function(j)
                            apply(object$variance[,j,j,], 2, sd, na.rm=TRUE),
                            simplify = "array")))
               se <- lapply(se, function(x) 
                 sqrt(x^2*(nboot-object$nonfit-1)^2/(nboot-object$nonfit)))
               zq <- qnorm(max(levels))
               lnames <- paste0(formatC(levels * 100, format = "fg", width = 1, 
                                        digits = getOption("digits")), "%")
               # the code above mimic stats:::format_perc(levels) which can't be used
               # because format_perc is not exported from stats
               out <- list(pro = array(as.double(NA), c(2,G),
                                       dimnames = list(lnames, 1:G)),
                           mean = array(as.double(NA), dim = c(2,d,G),
                                        dimnames = list(lnames, 1:d, 1:G)),
                           variance = array(as.double(NA), dim = c(2,d,G),
                                            dimnames = list(lnames, 1:d, 1:G)))
               out$pro[1,] <- ave$pro - zq*se$pro
               out$pro[2,] <- ave$pro + zq*se$pro
               out$mean[1,,] <- ave$mean - zq*se$mean
               out$mean[2,,] <- ave$mean + zq*se$mean
               out$variance[1,,] <- ave$variance - zq*se$variance
               out$variance[2,,] <- ave$variance + zq*se$variance
           } else
             { # percentile-based ci
               out <- list(pro = apply(object$pro, 2, quantile, probs = levels, na.rm=TRUE),
                           mean = apply(object$mean, c(2,3), quantile, probs = levels, na.rm=TRUE))
               v <- array(as.double(NA), dim = c(2,d,G),
                          dimnames = dimnames(out$mean))
               for(j in seq.int(d))
                 v[,j,] <- apply(object$variance[,j,j,], 2, quantile, probs = levels, na.rm=TRUE)
               out$variance <- v
             }
           }
  )

  obj <- append(object[c("modelName", "G", "nboot", "type")],
                list(d = d, what = what))
  if(what == "ci") obj$conf.level <- conf.level
  obj <- append(obj, out)
  class(obj) <- "summary.MclustBootstrap"
  return(obj)
}

print.summary.MclustBootstrap <- function(x, digits = getOption("digits"), ...)
{
  txt <- paste(rep("-", min(58, getOption("width"))), collapse = "")
  catwrap(txt)
  catwrap(paste("Resampling", 
                switch(x$what,
                       "se" = "standard errors",
                       "ave" = "averages",
                       "ci"  = "confidence intervals")))
  catwrap(txt)
  #
  cat(formatC("Model", flag = "-", width = 26), "=", x$modelName, "\n")
  cat(formatC("Num. of mixture components", flag = "-", width = 26), 
      "=", x$G, "\n")
  cat(formatC("Replications", flag = "-", width = 26), "=", x$nboot, "\n")
  cat(formatC("Type", flag = "-", width = 26), "=", 
      switch(x$type, 
             "bs"   = "nonparametric bootstrap",
             "wlbs" = "weighted likelihood bootstrap", 
             "pb"   = "parametric bootstrap",
             "jk" = "jackknife"),
      "\n")
  if(x$what == "ci")
    cat(formatC("Confidence level", flag = "-", width = 26), 
        "=", x$conf.level, "\n")
  #
  cat("\nMixing probabilities:\n")
  print(x$pro, digits = digits)
  #
  cat("\nMeans:\n")
  if(x$d == 1) 
    { if(x$what == "se" | x$what == "ave") 
        print(x$mean[1,], digits = digits)
      else
        print(x$mean[,1,], digits = digits) 
  } else
  if(x$what == "se" | x$what == "ave") 
      print(x$mean, digits = digits)
    else 
    { for(g in seq.int(x$G))
      { cat("[,,", g, "]\n", sep = "")
        print(x$mean[,,g], digits = digits) }
  }
  #
  cat("\nVariances:\n")
  if(x$d == 1)
    { print(x$variance[,1,], digits = digits) }
  else
    { for(g in seq.int(x$G))
         { cat("[,,", g, "]\n", sep = "")
           print(x$variance[,,g], digits = digits) }
  }
  
  invisible(x)
}

plot.MclustBootstrap <- function(x, what = c("pro", "mean", "var"), show.parest = TRUE, show.confint = TRUE, hist.col = "grey", hist.border = "lightgrey", breaks = "Sturges", col = "forestgreen", lwd = 2, lty = 3, xlab = NULL, xlim = NULL, ylim = NULL, ...)
{
  object <- x # Argh.  Really want to use object anyway
  what <- match.arg(what, choices = eval(formals(plot.MclustBootstrap)$what))
  par <- object$parameters
  d <- dim(object$mean)[2]
  varnames <- rownames(par$mean)
  if(show.confint)
    { ci <- summary(object, what = "ci", ...)
      ave <- summary(object, what = "ave", ...) 
    }
  
  histBoot <- function(boot, stat, ci, ave, breaks, xlim, ylim, xlab, ...)
  { 
    hist(boot, breaks = breaks, xlim = xlim, ylim = ylim,
         main = "", xlab = xlab, ylab = "",
         border = hist.border, col = hist.col)
    box()
    if(show.parest)
      abline(v = stat, col = col, lwd = lwd, lty = lty)
    if(show.confint)
      { lines(ci, rep(par("usr")[3]/2,2), lwd = lwd, col = col)
        points(ave, par("usr")[3]/2, pch = 15, col = col)
    }
  }

  switch(what, 
         "pro" = { xlim <- range(if(is.null(xlim)) pretty(object$pro) else xlim)
                   for(k in 1:object$G) 
                       histBoot(object$pro[,k], breaks = breaks, 
                                stat = par$pro[k], 
                                ci = ci$pro[,k],
                                ave = ave$pro[k],
                                xlim = xlim, ylim = ylim,
                                xlab = ifelse(is.null(xlab),
                                              paste("Mix. prop. for comp.",k),
                                              xlab))
         },
         "mean" = { isNull_xlim <- is.null(xlim)
                    for(j in 1:d)
                       { xlim <- range(if(isNull_xlim) pretty(object$mean[,j,]) 
                                       else xlim)
                         for(k in 1:object$G)
                            histBoot(object$mean[,j,k], breaks = breaks, 
                                     stat = par$mean[j,k], 
                                     ci = ci$mean[,j,k],
                                     ave = ave$mean[j,k],
                                     xlim = xlim, ylim = ylim,
                                     xlab = ifelse(is.null(xlab),
                                              paste(varnames[j], "mean for comp.",k),
                                              xlab))
                       }
         },
         "var" = { isNull_xlim <- is.null(xlim)
                   for(j in 1:d)
                      { xlim <- range(if(isNull_xlim) pretty(object$variance[,j,j,]) 
                                       else xlim)
                        for(k in 1:object$G)
                            histBoot(object$variance[,j,j,k], breaks = breaks, 
                                     stat = par$variance[j,j,k], 
                                     ci = ci$variance[,j,k],
                                     ave = ave$variance[j,k],
                                     xlim = xlim, ylim = ylim,
                                     xlab = ifelse(is.null(xlab),
                                              paste(varnames[j], "var. for comp.",k),
                                              xlab))
                       }
         }
        )  
  invisible()
}

