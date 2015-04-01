##
## Bootstrap methods
##

#
# Bootstrap Likelihood Ratio Test
#

mclustBootstrapLRT <- function(data, modelName = NULL, nboot = 999, level = 0.05, maxG = NULL, verbose = TRUE, ...)
{
  if(is.null(modelName))
    stop("A 'modelName' must be provided. Please see help(mclustModelNames) which describes the available models.")
  modelName <- modelName[1]
  if(is.null(maxG)) G <- seq.int(1, 9)
  else { maxG <- as.numeric(maxG); G <- seq.int(1, maxG+1) }
  Bic <- mclustBIC(data, G = G, modelNames = modelName, warn = FALSE, ...)
  if(!(modelName %in% attr(Bic, "modelNames")))
    stop("'modelName' not compatibile with data. Please see help(mclustModelNames) which describes the available models.")
  if(all(is.na(Bic)))
    stop(paste("no model", modelName, "can be fitted."))
  # select only models that can be fit
  G <- which(!is.na(Bic[, attr(Bic, "modelNames") == modelName]))
  # maxG <- max(G)
  # G <- setdiff(G, maxG)
  
  if(verbose & interactive()) 
    { cat("bootstrapping LRTS ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = (max(G)-1)*nboot, style = 3) 
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
    { b <- b + 1
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
      if(verbose & interactive()) 
        setTxtProgressBar(pbar, (g-1)*nboot+b)
    }
    p.value[g] <- (1 + sum(bootLRTS[,g] >= obsLRTS[g]))/(nboot+1)
    # check if not-significant when no maxG is provided
    if(is.null(maxG) & p.value[g] > level) 
      { continue <- FALSE
        if(verbose & interactive()) 
          setTxtProgressBar(pbar, (max(G)-1)*nboot) 
      }
  }
  if(verbose & interactive()) close(pbar)
  out <- list(G = 1:g, 
              modelName = modelName,
              obs = obsLRTS[1:g],
              boot = bootLRTS[,1:g],
              p.value = p.value[1:g])
  class(out) <- "mclustBootstrapLRT"
  return(out)
}

print.mclustBootstrapLRT <- function(x, ...)
{
  cat("Bootstrap sequential LRT for the number of mixture components\n") 
  cat(rep("-", 61), "\n", sep = "")
  cat(formatC("Model", flag = "-", width = 12), "=", x$modelName, "\n")
  cat(formatC("Replications", flag = "-", width = 12), "=", nrow(x$boot), "\n")
  df <- data.frame(x$obs, x$p.value)
  colnames(df) <- c("LRTS", "bootstrap p-value")
  rownames(df) <- formatC(paste(x$G, "vs", x$G+1), flag = "-", width = 8)
  print(df, ...)
}

plot.mclustBootstrapLRT <- function(x, G = 1, col = "lightgrey", border = "white", breaks = "scott", main = NULL, ...) 
{
  if(G > x$G) 
    { warning(paste("bootstrap LRT not available for G =", G)) 
      return() }
  G <- as.numeric(G)[1]
  h <- hist(x$boot[,G], breaks = breaks, plot = FALSE)
  xlim <- range(h$breaks, x$boot[,G], x$obs[G]*1.1, na.rm = TRUE)
  xlim <- c(xlim[1] - diff(xlim) * 0.1, xlim[2] + diff(xlim) * 0.1)
  plot(h, xlab = "LRTS", freq = FALSE, xlim = xlim,
       col = col, border = border, main = NULL)
  abline(v = x$obs[G], lty = 3, lwd = 2, col = "forestgreen")
  if(is.null(main) | is.character(main))
    { if(is.null(main)) main <- paste("Bootstrap LRT for model", x$modelName, 
                                      "with", G, "vs", G+1, "components")
      title(main = main, cex.main = 1) }
  invisible()
}

#
# Bootstrap inference (standard errors and percentile confidence intervals) 
#

MclustBootstrap <- function(object, nboot = 999, verbose = TRUE, ...)
{
  data <- object$data
  n <- object$n
  d <- object$d
  G <- object$G
  varnames <- rownames(object$parameters$mean)
  pro.boot  <- array(NA, c(nboot,G), 
                     dimnames = list(NULL, seq.int(G)))
  mean.boot <- array(NA, c(nboot,d,G), 
                     dimnames = list(NULL, varnames, seq.int(G)))
  var.boot  <- array(NA, c(nboot,d,d,G),
                     dimnames = list(NULL, varnames, varnames, seq.int(G)))

  if(verbose & interactive()) 
    { cat("bootstrapping Mclust model ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0, max = nboot, style = 3) 
    }
  b <- 0
  while(b < nboot)
  { b <- b + 1
    idx <- sample(seq_len(n), size = n, replace = TRUE)
#    w <- tabulate(idx, nbins = n)
#    w <- w/max(w)
#     mod.boot <- try(do.call("me.weighted", 
#                              c(list(weights = w, warn = FALSE), object)),
#                      silent = TRUE)
    obj <- object
    obj$data <- object$data[idx,]
    obj$z <- obj$z[idx,]
    obj$warn <- FALSE
    mod.boot <- try(do.call("me", obj), silent = TRUE)
    # check model convergence
    if(inherits(mod.boot, "try-error"))
      { b <- b - 1; next() }
    if(is.na(mod.boot$loglik))
      { b <- b - 1; next() }
    #
    pro.boot[b,]   <- mod.boot$parameters$pro
    mean.boot[b,,] <- mod.boot$parameters$mean
    var.boot[b,,,] <- mod.boot$parameters$variance$sigma
    if(verbose & interactive()) setTxtProgressBar(pbar, b)
  }
  if(verbose & interactive()) close(pbar)
  
  out <- list(G = G, modelName = object$modelName, 
              nboot = nboot, pro = pro.boot, 
              mean = mean.boot, variance = var.boot)
  class(out) <- "MclustBootstrap"
  return(out)
}

print.MclustBootstrap <- function(x,  digits = getOption("digits"), ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  str(x,1)
  invisible()
}

summary.MclustBootstrap <- function(object, what = c("se", "ci"), conf.level = 0.95, ...)
{
  what <- match.arg(what, several.ok = FALSE)
  dims <- dim(object$mean)
  varnames <- dimnames(object$mean)[[2]]
  nboot <- dims[1]
  d <- dims[2]
  G <- dims[3]
      
  if(what == "se")
    { out <- list(pro = apply(object$pro, 2, sd),
                  mean = apply(object$mean, c(2,3), sd),
                  variance = apply(object$variance, c(2,3,4), sd))
  } else 
  if(what == "ci")
    { levels <- c((1-conf.level)/2, (1 + conf.level)/2)
      out <-  list(pro = apply(object$pro, 2, quantile, probs = levels),
                   mean = apply(object$mean, c(2,3), quantile, probs = levels))
      v <- array(as.double(NA), dim = c(2,d,G), 
                 dimnames = dimnames(out$mean))
      for(j in seq.int(d))
         v[,j,] <- apply(object$variance[,j,j,], 2, quantile, probs = levels)
     out$variance <- v
  }
  
  obj <- append(object[c("modelName", "G", "nboot")], list(d = d, what = what))
  if(what == "ci") obj$conf.level <- conf.level
  obj <- append(obj, out)
  class(obj) <- "summary.MclustBootstrap"
  return(obj)
}

print.summary.MclustBootstrap <- function(x, digits = getOption("digits"), ...)
{
  if(x$what == "se")
    cat("Bootstrap standard errors\n")
  else  
    cat("Bootstrap confidence intervals\n")
  cat(rep("-", 40),"\n",sep="")
  cat(formatC("Model", flag = "-", width = 28), "=", x$modelName, "\n")
  cat(formatC("Num. of mixture components", flag = "-", width = 28), 
      "=", x$G, "\n")
  cat(formatC("Replications", flag = "-", width = 28), "=", x$nboot, "\n")
  if(x$what == "ci")
    cat(formatC("Confidence level", flag = "-", width = 28), 
        "=", x$conf.level, "\n")

  cat("\nMixing probabilities:\n")
  print(x$pro, digits = digits)
  #
  cat("\nMeans:\n")
#   if(x$what == "se") 
#     { if(dim(x$mean)[1] == 1) print(x$mean[1,], digits = digits)
#       else print(x$mean, digits = digits) }
#   else
#     { for(g in seq.int(x$G))
#          { cat("[,,", g, "]\n", sep = "")
#            print(x$mean[,,g], digits = digits) }
#   }
  if(x$d == 1) 
    { if(x$what == "se") print(x$mean[1,], digits = digits)
      else               print(x$mean[,1,], digits = digits) 
  } else
  if(x$what == "se") print(x$mean, digits = digits)
    else { for(g in seq.int(x$G))
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


