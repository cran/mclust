# Semi-Supervised Classification

MclustSSC <- function(data, class, 
	                    G = NULL, modelNames = NULL, 
                      prior = NULL, control = emControl(), 
											warn = mclust.options("warn"),
                      verbose = interactive(), 
											...) 
{
  call <- match.call()
  data <- data.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  oneD <- if(d==1) TRUE else FALSE
  #
  class <- factor(class, exclude = NA)
  nclass <- nlevels(class)
  #
  if(is.null(G))
    G <- nclass
  if(any(G < nclass))
    stop("G cannot be smaller than the number of classes")
  G <- G[G >= nclass][1]
  #
  if(is.null(modelNames)) 
  { 
    modelNames <- if(oneD) c("E", "V") 
                  else     mclust.options("emModelNames")
  }
  #
  if(n <= d) 
  {
    m <- match(c("EEE","EEV","VEV","VVV"),
               mclust.options("emModelNames"), nomatch=0)
    modelNames <- modelNames[-m]
  }
  nModelNames <- length(modelNames)
  
  if(verbose) 
  { 
    cat("fitting ...\n")
    flush.console()
    pbar <- txtProgressBar(min = 0, max = nModelNames, style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  
  args <- list(data = data, class = class, G = G, ...)
  Model <- NULL
  BIC <- rep(as.double(NA), length(modelNames))
  for(m in seq(nModelNames))
  { 
    mod <- try(do.call("MclustSSC.fit",
                       c(args, list(modelName = modelNames[m]))),
               silent = TRUE)
    if(verbose) 
    { ipbar <- ipbar+1; setTxtProgressBar(pbar, ipbar) }
    if(inherits(mod, "try-error")) next()
    BIC[m] <- mod$bic
    if(!is.na(BIC[m]) && BIC[m] >= max(BIC, na.rm = TRUE))
      Model <- mod
  }
  if(all(is.na(BIC)))
  { 
    warning("No model(s) can be estimated!!")
    return() 
  }
  BIC <- matrix(BIC, nrow = 1, dimnames = list(G, modelNames))

  out <- c(list(call = call, data = data, class = class, BIC = BIC,
            		control = control), Model)
  orderedNames <- c("call", "data", "class",
	                  "modelName", "G", "n", "d", 
				            "BIC", "loglik", "df", "bic", 
										"parameters", "z", "classification",
										"prior", "control")
  out <- structure(out[orderedNames], 
		               class = "MclustSSC")  
  return(out)
}

print.MclustSSC <- function(x, ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  cat("\n")
  catwrap("\nAvailable components:\n")
  print(names(x))
  # str(x, max.level = 2, give.attr = FALSE, strict.width = "wrap")
  invisible(x)
}

summary.MclustSSC <- function(object, parameters = FALSE, ...)
{
  # collect info
  nclass <- nlevels(object$class)
  classes <- levels(object$class)
  G <- object$G
  printParameters <- parameters
  class <- object$class
  classif <- object$classification
  classifNames <- levels(object$classification)
  err <- classError(class[!is.na(class)], 
                    classif[!is.na(class)])$errorRate
  # n <- c(table(class, useNA = "always"))
  n <- tabulate(class, nbins = G)
  names(n) <- classifNames
  if(any(is.na(class)))
    n <- c(n, "<NA>" = sum(is.na(class)))
  tab <- table("Class" = class, "Predicted" = classif, useNA = "ifany")
  noise <- FALSE
  # todo:
  # noise <- if(is.na(object$hypvol)) FALSE else object$hypvol
  pro <- object$parameters$pro
  if(is.null(pro)) pro <- 1
  names(pro) <- if(noise) c(classifNames,0) else classifNames
  mean <- object$parameters$mean
  colnames(mean) <- names(pro)
  if(object$d > 1)
  { 
    sigma <- object$parameters$variance$sigma
    dimnames(sigma)[[3]] <-  names(pro)
  } else
  { 
    sigma <- rep(object$parameters$variance$sigmasq, object$G)[1:object$G]
    names(sigma) <- names(mean) 
  }

  obj <- list(n = n, d = object$d,
              loglik = object$loglik, 
              df = object$df, bic = object$bic,
              nclass = nclass, classes = classes,
              G = object$G, modelName = object$modelName,
              pro = pro, mean = mean, variance = sigma,
              noise = noise, prior = object$prior,
              tab = tab, err = err,
              printParameters = printParameters)
  class(obj) <- "summary.MclustSSC"
  return(obj)
}

print.summary.MclustSSC <- function(x, digits = getOption("digits"), ...)
{
  
  title <- paste("Gaussian finite mixture model for semi-supervised classification")
  txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")
  catwrap(txt)
  catwrap(title)
  catwrap(txt)

  cat("\n")
  tab <- data.frame("log-likelihood" = x$loglik,
                     "n" = sum(x$n), "df" = x$df, 
                     "BIC" = x$bic, 
                    row.names = "", check.names = FALSE)
  print(tab, digits = digits)

  tab <- data.frame("n" = x$n, "%" = round(x$n/sum(x$n)*100,2), 
                    "Model" = c(rep(x$modelName, x$G), ""),
                    "G" = c(rep(1, x$G), ""),
                    check.names = FALSE,
                    row.names = ifelse(is.na(names(x$n)), 
                                       "<NA>", names(x$n)))
  tab <- as.matrix(tab)
  names(dimnames(tab)) <- c("Classes", "")
  print(tab, quote = FALSE, right = TRUE)

  if(!is.null(x$prior))
  { cat("\nPrior: ")
    cat(x$prior$functionName, "(", 
          paste(names(x$prior[-1]), x$prior[-1], sep = " = ", collapse = ", "), 
          ")", sep = "")
   cat("\n")
  }

  if(x$printParameters)
  {
    cat("\nMixing probabilities:\n")
    print(x$pro, digits = digits)
    cat("\nMeans:\n")
    print(x$mean, digits = digits)
    cat("\nVariances:\n")
    if(x$d > 1) 
    { 
      for(g in 1:x$G)
         { cat(names(x$pro)[g], "\n")
           print(x$variance[,,g], digits = digits) }
    }
    else print(x$variance, digits = digits)
    if(x$noise)
      { cat("\nHypervolume of noise component:\n")
        cat(signif(x$noise, digits = digits), "\n") }
  }

  cat("\nClassification summary:\n")
  print(x$tab)

  invisible(x)
}

MclustSSC.fit <- function(data, class, 
                          G = NULL, modelName = NULL, 
                          prior = NULL, control = emControl(), 
													warn = NULL, .verbose = FALSE, ...) 
{
  data <- data.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  class <- factor(class, exclude = NA)
  nclass <- nlevels(class)
  known.class <- which(!is.na(class))
  unknown.class <- which(is.na(class))
  if(is.null(G)) G <- nclass
  if(is.null(modelName)) 
    stop("modelName must be specified!")
  #browser()
  
  # initialization of z matrix by filling with 0/1 for observations 
  # with known labels
  z <- matrix(0.0, nrow = n, ncol = G)
  for(k in 1:nclass)
    z[class == levels(class)[k], k] <- 1
  # store the z which should not be updated
  z0 <- z[known.class,,drop=FALSE]
  # initialization of unlabeled data...
  if(G > nclass)
  { 
    # via k-means if unobserved classes 
    km <- kmeans(data[unknown.class,,drop=FALSE],
                 centers = G,
                 nstart = 25, iter.max = 100)
    # z[unknown.class,] <- unmap(km$cluster)
    z[unknown.class,] <- rep(1,length(unknown.class)) %o% km$size/sum(km$size)
  } else
  {
    # by equal proportion otherwise
    z[unknown.class,] <- 1/G
  }
  
  loglik0 <- -Inf
  criterion <- TRUE
  iter <- 0
  if(.verbose)
    cat("\nmodelName =", modelName, "\n")
  #
  while(criterion) 
  {
    iter <- iter + 1
    fit.m <- do.call("mstep", list(data = data, z = z, 
                                   modelName = modelName, 
                                   prior = prior, 
																	 control = control, 
																	 warn = warn))
    fit.e <- do.call("estep", c(list(data = data,
                                     control = control, 
																		 warn = warn), 
                                fit.m))
    z <- fit.e$z
    z[known.class,] <- z0
    ldens <- do.call("dens", c(list(data = data[-known.class,,drop=FALSE],
                                    logarithm = TRUE), fit.m))
    lcdens <- do.call("cdens", c(list(data = data[known.class,,drop=FALSE], 
                                      logarithm = TRUE), fit.m))
    lcdens <- sweep(lcdens, MARGIN = 2, FUN = "+", 
                    STATS = log(fit.m$parameters$pro))
    loglik <- sum(ldens) + sum(lcdens * z0)
    criterion <- ( iter < control$itmax[1] & 
                   (loglik - loglik0) > control$tol[1] )
    # print(loglik - loglik0)
    loglik0 <- loglik
    if(.verbose)
      cat("iter =", iter, "  loglik =", loglik0, "\n")
  }
  fit <- fit.m
  fit$loglik <- loglik
  fitclass <- map(fit$z, warn = FALSE)
  # assign labels of known classes
  fitclass <- factor(fitclass)
  labels <- levels(class)
  if(G > nclass) 
    labels <- c(labels, paste0("class", seq(nclass+1,G)))
  levels(fitclass) <- labels
  fit$classification <- fitclass

  fit$df <- (G-1) + p*nclass + nVarParams(fit$modelName, d = p, G = nclass)
  fit$bic <- 2*fit$loglik - fit$df*log(n)
  #
  return(fit)
}

plot.MclustSSC <- function(x, what = c("BIC", "classification", "uncertainty"), ...)
{
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "MclustSSC")) 
    stop("object not of class 'MclustSSC'")
  class(object) <- c(class(object), "Mclust")
  
  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  
  plot.MclustSSC.bic <- function(...)
  {
    dotchart(rev(object$BIC[1,]), pch = 19, xlab = paste("BIC for G =", object$G), ...)
  }
  
  if(interactive() & length(what) > 1)
    { title <- "Model-based semi-supervised classification plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "BIC")            plot.MclustSSC.bic(...)
             if(what[choice] == "classification") plot.Mclust(object, what = "classification", ...)
             if(what[choice] == "uncertainty")    plot.Mclust(object, what = "uncertainty", ...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  } 
  else 
    { if(any(what == "BIC"))            plot.MclustSSC.bic(...)
      if(any(what == "classification")) plot.Mclust(object, what = "classification", ...)
      if(any(what == "uncertainty"))    plot.Mclust(object, what = "uncertainty", ...) 
  }
    
  invisible()
}
  
predict.MclustSSC <- function(object, newdata, ...)
{
  if(!inherits(object, "MclustSSC")) 
    stop("object not of class 'MclustSSC'")
  if(missing(newdata))
    { newdata <- object$data }
  newdata <- as.matrix(newdata)
  if(ncol(object$data) != ncol(newdata))
    { stop("newdata must match ncol of object data") }
  #
	object$data <- newdata
  z <- do.call("cdens", c(object, list(logarithm = TRUE)))
  pro <- object$parameters$pro
	logpro <- log(pro) - log(sum(pro))
  noise <- FALSE # (!is.na(object$hypvol))
  z <- if(noise) cbind(z, log(object$parameters$Vinv))
       else      cbind(z) # drop redundant attributes
  # TODO: to be removed at a certain point
  # z <- sweep(z, MARGIN = 2, FUN = "+", STATS = logpro)
  # z <- sweep(z, MARGIN = 1, FUN = "-", STATS = apply(z, 1, logsumexp_old))
  # z <- exp(z)
  z <- softmax(z, logpro)
  cl <- c(levels(object$classification), if(noise) 0)
  colnames(z) <- cl
  cl <- factor(cl[apply(z, 1, which.max)], levels = cl)
  out <- list(classification = cl, z = z)
  return(out) 
}
