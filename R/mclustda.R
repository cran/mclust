MclustDA <- function(data, class, G = NULL, modelNames = NULL, 
                     modelType = c("MclustDA", "EDDA"), 
                     prior = NULL, control = emControl(), 
                     initialization = NULL, warn = mclust.options("warn"),
                     verbose = interactive(), ...) 
{
  call <- match.call()
  mc <- match.call(expand.dots = TRUE)
  #
  if(missing(data))
    stop("no training data provided!")
  data <- data.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  oneD <- if(p==1) TRUE else FALSE
  #
  if(missing(class))
    stop("class labels for training data must be provided!")
  class <- as.factor(class)
  classLabel <- levels(class)
  ncl <- nlevels(class)
  if(ncl == 1) G <- 1
  prop <- as.vector(table(class))/n
  names(prop) <- classLabel
  #
  modelType <- match.arg(modelType, 
                         choices = eval(formals(MclustDA)$modelType), 
                         several.ok = FALSE)
  #
  if(is.null(G)) 
    { G <- rep(list(1:5), ncl) }
  else if(is.list(G))
    { G <- lapply(G, sort) }
  else 
    { G <- rep(list(sort(G)), ncl) }
  if(any(unlist(G) <= 0))
    stop("G must be positive")
  #
  if(is.null(modelNames)) 
    { if(oneD) modelNames <- c("E", "V")
      else     modelNames <- mclust.options("emModelNames")
  }
  if(n <= p) 
    { m <- match(c("EEE","EEV","VEV","VVV"), mclust.options("emModelNames"), nomatch=0)
      modelNames <- modelNames[-m]
  }
  if(!is.list(modelNames))
    { modelNames <- rep(list(modelNames), ncl) }
  #
  hcUse <- mclust.options("hcUse")
  mclust.options("hcUse" = "VARS")
  on.exit(mclust.options("hcUse" = hcUse))
  #
  if(modelType == "EDDA")
  { 
    mc[[1]] <- as.name("mstep")
    mc$class <- mc$G <- mc$modelNames <- mc$modelType <- NULL
    mc$warn <- FALSE
    mc$z <- unmap(as.numeric(class))
    G <- 1
    modelNames <- unique(unlist(modelNames))
    BIC <- rep(NA, length(modelNames))
    Model <- NULL
    if(verbose) 
      { cat("fitting ...\n")
        flush.console()
        pbar <- txtProgressBar(min = 0, max = length(modelNames), style = 3) 
        on.exit(close(pbar))
        ipbar <- 0
    }
    for(i in seq(modelNames))
       { mc$modelName <- as.character(modelNames[i])
         mStep <- eval(mc, parent.frame())
         eStep <- do.call("estep", c(mStep, list(data = data, warn = FALSE)))
         BIC[i] <- do.call("bic", c(eStep, list(equalPro = TRUE)))
         if(!is.na(BIC[i]) && BIC[i] >= max(BIC, na.rm = TRUE))
           Model <- eStep
         if(verbose) 
           { ipbar <- ipbar+1; setTxtProgressBar(pbar, ipbar) }
    }
    if(all(is.na(BIC)))
      { warning("No model(s) can be estimated!!")
        return() }
    names(BIC) <- modelNames
    bic <- max(BIC, na.rm = TRUE)
    loglik <- Model$loglik
    df <- (2*loglik - bic)/log(Model$n)
    # there are (nclass-1) more df than real needed
    # equal to logLik(object) but faster
    Model <- c(Model, list("BIC" = BIC))
    Models <- rep(list(Model), ncl)
    names(Models) <- classLabel
    for(l in 1:ncl)
       { I <- (class == classLabel[l]) 
         Models[[l]]$n <- sum(I)
         Models[[l]]$G <- 1
         Models[[l]]$bic <- Models[[l]]$loglik <- NULL
         par <- Models[[l]]$parameters
         par$pro <- 1
         par$mean <- if(oneD) par$mean[l] else par$mean[,l,drop=FALSE]
         par$variance$G <- 1
         if(oneD)
           { # par$variance$sigma <- par$variance$sigma[l]
             if(length(par$variance$sigmasq) > 1)
               par$variance$sigmasq <- par$variance$sigmasq[l]
             else
               par$variance$sigmasq <- par$variance$sigmasq
         }
         else
           { par$variance$sigma <- par$variance$sigma[,,l,drop=FALSE]
             if(length(par$variance$sigmasq) > 1)
               par$variance$sigmasq <- par$variance$sigmasq[l]
             if(length(par$variance$scale) > 1)
               par$variance$scale <- par$variance$scale[l]
             if(length(dim(par$variance$shape)) > 1)
               par$variance$shape <- par$variance$shape[,l]
             if(length(dim(par$variance$orientation)) > 2)  # LS was > 1
               par$variance$orientation <-
                 par$variance$orientation[,,l,drop=FALSE]
             if(length(dim(par$variance$cholSigma)) > 2) 
               par$variance$cholSigma <-
                 par$variance$cholSigma[,,l,drop=FALSE]
             if(length(dim(par$variance$cholsigma)) > 2) 
               par$variance$cholsigma <-
                par$variance$cholsigma[,,l,drop=FALSE]
           }
         Models[[l]]$parameters <- par
         Models[[l]]$z <- NULL # z[I,,drop=FALSE]
         Models[[l]]$classification <- rep(1, sum(I)) # apply(z[I,,drop=FALSE], 1, which.max)
         Models[[l]]$uncertainty <- NULL # 1 - apply(z[I,], 1, max)
         Models[[l]]$observations <- which(I)     
    }
  }
  else
  { # modelType == "MclustDA" i.e. different covariance structures for each class
    Models <- rep(list(NULL), ncl)
    mc[[1]] <- as.name("mclustBIC")
    mc$class <- NULL
    for(l in 1:ncl) 
       { I <- (class == classLabel[l])
         mc[[2]] <- data[I,]
         mc$G <- G[[l]]
         mc$modelNames <- as.character(modelNames[[l]])
         if(verbose) cat(paste0("Class ", classLabel[l], ": "))
         BIC <- eval(mc, parent.frame())
         # slightly adjust parameters if none of the models can be fitted
         while(all(is.na(BIC)))
         { if(length(mc$modelNames) == 1)
             { j <- which(mc$modelNames == mclust.options("emModelNames"))
               if(j == 1) mc$G <- mc$G - 1
               else       mc$modelNames <- mclust.options("emModelNames")[j-1]
           }
           else
             { mc$G <- mc$G - 1 }
           BIC <- eval(mc, parent.frame())
         }
         SUMMARY <- summary(BIC, data[I,])
         SUMMARY$bic <- BIC; 
         names(SUMMARY)[which(names(SUMMARY) == "bic")] <- "BIC"
         Models[[l]] <- c(SUMMARY, list(observations = which(I)))
    }
    bic <- loglik <- df <- NULL
  }
  
  names(Models) <- classLabel
  Models$Vinv <- NULL
  out <- list(call = call, data = data, class = class,
              type = modelType, n = n, d = p, prop = prop, 
              models = Models, bic = bic, loglik = loglik, df = df)
  out <- structure(out, prior = prior, control = control, 
                   class = "MclustDA")
  
  if(modelType == "MclustDA") 
  { 
    l <- logLik.MclustDA(out, data)
    out$loglik <- as.numeric(l)
    out$df <- attr(l, "df")
    out$bic <- 2*out$loglik - log(n)*out$df
  }
  
  return(out)
}

print.MclustDA <- function(x, ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  models <- x$models
  nclass <- length(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  M <- sapply(1:nclass, function(i) models[[i]]$modelName)
  G <- sapply(1:nclass, function(i) models[[i]]$G)
  out <- data.frame(n = n, Model = M, G = G)
  rownames(out) <- names(models)
  out <- as.matrix(out)
  names(dimnames(out)) <- c("Classes", "")
  print(out, quote = FALSE, right = TRUE)
  cat("\n")
  catwrap("\nAvailable components:\n")
  print(names(x))
  # str(x, max.level = 2, give.attr = FALSE, strict.width = "wrap")
  invisible(x)
}

summary.MclustDA <- function(object, parameters = FALSE, newdata, newclass, ...)
{
  # collect info
  models <- object$models
  nclass <- length(models)
  classes <- names(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  G <- sapply(1:nclass, function(i) models[[i]]$G)
  modelName <- sapply(1:nclass, function(i) models[[i]]$modelName)
  prior <- attr(object, "prior")
  printParameters <- parameters
  par <- getParameters.MclustDA(object)
  class <- object$class
  data <- object$data
  pred <- predict(object, newdata = data, ...)
  err <- mean(class != pred$classification)
  brier <- BrierScore(pred$z, class)
  tab <- try(table(class, pred$classification))
  if(class(tab) == "try-error") 
    { err <- tab <- NA }
  else 
    { names(dimnames(tab)) <- c("Class", "Predicted") }
  
  tab.newdata <- err.newdata <- brier.newdata <- NULL
  if(!missing(newdata) & !missing(newclass))
  { 
    pred.newdata <- predict(object, newdata = newdata, ...)
    if(missing(newclass))
    { tab.newdata <- table(pred.newdata$classification)
      names(dimnames(tab.newdata)) <- "Predicted"
    }
    else
    { tab.newdata <- table(newclass, pred.newdata$classification)
      names(dimnames(tab.newdata)) <- c("Class", "Predicted")
      err.newdata <- mean(newclass != pred.newdata$classification)
      brier.newdata <- BrierScore(pred.newdata$z, newclass)
    }
  }
  
  obj <- list(type = object$type, n = n, d = object$d,
              loglik = object$loglik, df = object$df, bic = object$bic,
              nclass = nclass, classes = classes,
              G = G, modelName = modelName,
              prop = object$prop, parameters = par, prior = prior, 
              tab = tab, err = err, brier = brier,
              tab.newdata = tab.newdata, 
              err.newdata = err.newdata, brier.newdata = brier.newdata,
              printParameters = printParameters)
  class(obj) <- "summary.MclustDA"
  return(obj)
}

print.summary.MclustDA <- function(x, digits = getOption("digits"), ...)
{
  title <- paste("Gaussian finite mixture model for classification")
  txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")
  catwrap(txt)
  catwrap(title)
  catwrap(txt)
  
  cat("\n")
  catwrap(paste(x$type, "model summary:"))
  cat("\n")
  #
  tab <- data.frame("log-likelihood" = x$loglik,
                    "n" = sum(x$n), "df" = x$df, 
                    "BIC" = x$bic, 
										row.names = "", check.names = FALSE)
  print(tab, digits = digits)
  
  tab <- data.frame("n" = x$n, "%" = round(x$n/sum(x$n)*100,2), 
                    "Model" = x$modelName, "G" = x$G,
                    check.names = FALSE,
                    row.names = x$classes)
  tab <- as.matrix(tab)
  names(dimnames(tab)) <- c("Classes", "")
  print(tab, digits = digits, quote = FALSE, right = TRUE)
  
  if(!is.null(x$prior))
  { cat("\nPrior: ")
    cat(x$prior$functionName, "(", 
        paste(names(x$prior[-1]), x$prior[-1], sep = " = ", collapse = ", "), 
        ")", sep = "")
    cat("\n")
  }
  
  if(x$printParameters)
  {
    cat("\nClass prior probabilities:\n")
    print(x$prop, digits = digits)
    for(i in seq(x$nclass))
    { cat("\nClass = ", x$class[i], "\n", sep = "")
      par <- x$parameters[[i]]
      if(x$type == "MclustDA")
      {
        cat("\nMixing probabilities: ")
        cat(round(par$pro, digits = digits), "\n")
      }
      cat("\nMeans:\n")
      print(par$mean, digits = digits)
      cat("\nVariances:\n")
      if(x$d > 1)
      { for(g in seq(x$G[i]))
      { cat("[,,", g, "]\n", sep = "")
        print(par$variance[,,g], digits = digits) }
      }
      else print(par$variance, digits = digits)          
    }
  }
  
  cat("\nTraining confusion matrix:\n")
  print(x$tab)
  cat("Classification error =", round(x$err, min(digits,4)), "\n")
  cat("Brier score          =", round(x$brier, min(digits,4)), "\n")
  
  if(!is.null(x$tab.newdata)) 
  {
    cat("\nTest confusion matrix:\n")
    print(x$tab.newdata)
    if(!is.null(x$err.newdata))
    { cat("Classification error =", round(x$err.newdata, min(digits,4)), "\n") 
      cat("Brier score          =", round(x$brier.newdata, min(digits,4)), "\n")
    }
  }
  
  invisible(x)
}

getParameters.MclustDA <- function(object)
{
  # collect info
  models <- object$models
  nclass <- length(models)
  classes <- names(models)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  G <- sapply(1:nclass, function(i) models[[i]]$G)
  modelName <- sapply(1:nclass, function(i) models[[i]]$modelName)
  # prior <- attr(object, "prior")
  par <- vector(mode = "list", length = nclass)
  for(i in seq(nclass))
  { par[[i]] <- models[[i]]$parameters
    if(is.null(par[[i]]$pro)) par$pro <- 1
    if(par[[i]]$variance$d < 2)
    { sigma <- rep(par[[i]]$variance$sigma,
                   models[[i]]$G)[1:models[[i]]$G]
      names(sigma) <- names(par[[i]]$mean)
      par[[i]]$variance$sigma <- sigma
    }
    par[[i]]$variance <- par[[i]]$variance$sigma
  }
  return(par)
}

logLik.MclustDA <- function (object, data, ...) 
{
  if(missing(data)) 
    data <- object$data
  n <- object$n
  d <- object$d
  par <- getParameters.MclustDA(object)
  nclass <- length(par)
  fclass <- sapply(object$models, function(m) m$n)/n
  logfclass <- log(fclass)
  G <- sapply(par, function(x) length(x$pro))
  if(object$type == "EDDA") 
    { df <- d * nclass + nVarParams(object$models[[1]]$modelName, 
                                    d = d, G = nclass)
  }
  else 
    { df <- sum(sapply(object$models, function(mod) with(mod, 
                       (G - 1) + G * d + nVarParams(modelName, d = d, G = G))))
  }
  # ll <- sapply(object$models, function(mod) 
  #       { do.call("dens", c(list(data = data, logarithm = FALSE), mod)) })
  # l <- sum(log(apply(ll, 1, function(l) sum(fclass*l))))
  ll <- sapply(object$models, function(mod) 
        { do.call("dens", c(list(data = data, logarithm = TRUE), mod)) })
  l <- sum(apply(ll, 1, function(l) logsumexp(logfclass+l)))
    
  attr(l, "nobs") <- n
  attr(l, "df") <- df
  class(l) <- "logLik"
  return(l)
}

predict.MclustDA <- function(object, newdata, prop = object$prop, ...)
{
  
  if(!inherits(object, "MclustDA")) 
    stop("object not of class \"MclustDA\"")
  
  models <- object$models
  nclass <- length(models)
  classNames <- levels(object$class)
  n <- sapply(1:nclass, function(i) models[[i]]$n)
  if(missing(newdata))
    { newdata <- object$data }
  if(object$d == 1) newdata <- as.vector(newdata)
  if(is.numeric(prop))
  {
    if(any(prop < 0)) 
      stop("'prop' must be nonnegative")
    if(length(prop) != nclass) 
      stop("'prop' is of incorrect length")
    prop <- prop/sum(prop)
  } else
  {
    prop <- n/sum(n)
  }

  # class density computed on log scale
  densfun <- function(mod, data)
  { 
    do.call("dens", c(list(data = data, logarithm = TRUE), mod)) 
  }
  #
  z <- as.matrix(data.frame(lapply(models, densfun, data = newdata)))
  z <- sweep(z, MARGIN = 2, FUN = "+", STATS = log(prop))
  z <- sweep(z, MARGIN = 1, FUN = "-", STATS = apply(z, 1, logsumexp))
  z <- exp(z)
  colnames(z) <- classNames
  cl <- apply(z, 1, which.max)
  class <- factor(classNames[cl], levels = classNames)
  #
  out <- list(classification = class, z = z)
  return(out) 
}

plot.MclustDA <- function(x, what = c("scatterplot", "classification", "train&test", "error"), 
                          newdata, newclass, 
                          dimens = NULL, 
                          symbols, colors, 
                          main = NULL,
                          ...)
{
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "MclustDA")) 
    stop("object not of class \"MclustDA\"")
  
  data <- object$data
  if(object$d > 1) dataNames <- colnames(data)
  else             dataNames <- deparse(object$call$data)
  n <- nrow(data)
  p <- ncol(data)
  dimens <- if(is.null(dimens)) seq(p) else dimens[dimens <= p]
  d <- length(dimens)
  if(d == 0)
  {
    warning("dimens larger than data dimensionality...")
    return(invisible())
  }
  
  if(missing(newdata))
    { newdata <- matrix(as.double(NA), 0, p) }
  else
    { newdata <- as.matrix(newdata) }
  if(ncol(newdata) != p)
    stop("incompatible newdata dimensionality")
  if(missing(newclass))
    { newclass <- vector(length = 0) }
  else
    { if(nrow(newdata) != length(newclass))
      stop("incompatible newdata and newclass") }
  if(object$d > 1) newdataNames <- colnames(newdata)
  else             newdataNames <- deparse(match.call()$newdata)
  
  what <- match.arg(what, several.ok = TRUE)
  models <- object$models
  M <- length(models)
  if(missing(dimens)) dimens <- 1:p
  trainClass <- object$class
  nclass <- length(unique(trainClass))
  Data <- rbind(data, newdata)
  predClass <- predict(object, Data)$classification
  
  if(missing(symbols)) 
    { if(M <= length(mclust.options("classPlotSymbols"))) 
        { symbols <- mclust.options("classPlotSymbols") }
      else if(M <= 26) 
             { symbols <- LETTERS }
  }
  if(length(symbols) == 1) symbols <- rep(symbols,M)
  if(length(symbols) < M & !any(what == "train&test"))
    { warning("more symbols needed to show classification")
      symbols <- rep(16, M) }
  
  if(missing(colors))
    { colors <- mclust.options("classPlotColors") }
  if(length(colors) == 1) colors <- rep(colors,M)
  if(length(colors) < M & !any(what == "train&test"))
    { warning("more colors needed to show classification")
      colors <- rep("black", M) }
  
  oldpar <- par(no.readonly = TRUE)

  plot.MclustDA.scatterplot <- function(...)
  {
    if(d == 1)
    { 
      mclust1Dplot(data = if(nrow(newdata) == 0) data[,dimens[1],drop=FALSE]
                          else                   newdata[,dimens[1],drop=FALSE],
                   what = "classification",
                   classification = if(nrow(newdata) == 0) trainClass
                                    else                   newclass,
                   xlab = if(nrow(newdata) == 0) dataNames[dimens]
                          else                   newdataNames[dimens], 
                   ylab = "Classes",
                   main = NULL, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) 
                     if(nrow(newdata) == 0) "Training data" else "Test data"
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }

    scatellipses <- function(data, dimens, nclass, symbols, colors, ...)
    {
      m <- lapply(models, function(m) 
      { m$parameters$mean <- array(m$parameters$mean[dimens,], c(2,m$G))
        m$parameters$variance$sigma <- 
          array(m$parameters$variance$sigma[dimens,dimens,], c(2,2,m$G))
        m
      })
      plot(data[,dimens], type = "n", ...)
      for(l in 1:nclass) 
      { 
        I <- m[[l]]$observations
        points(data[I,dimens[1]], data[I,dimens[2]], 
               pch = symbols[l], col = colors[l])
        for(g in 1:(m[[l]]$G))
        { 
          mvn2plot(mu = m[[l]]$parameters$mean[,g], 
                   sigma = m[[l]]$parameters$variance$sigma[,,g],
                   k = 15,
                   fillEllipse = mclust.options("fillEllipses"),
                   col = if(mclust.options("fillEllipses")) 
                           colors[l] else rep("grey30",3))
        }
      }
    }
    
    if(d == 2)
    { 
      scatellipses(if(nrow(newdata) == 0) data else newdata, 
                   dimens = dimens[1:2], 
                   nclass = nclass, 
                   symbols = symbols, colors = colors, 
                   xlab = if(nrow(newdata) == 0) dataNames[dimens[1]]
                          else                   newdataNames[dimens[1]],
                   ylab = if(nrow(newdata) == 0) dataNames[dimens[2]]
                          else                   newdataNames[dimens[2]],
                   ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) 
                     if(nrow(newdata) == 0) "Training data" else "Test data"
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    if(d > 2)
    { 
      on.exit(par(oldpar))
      par(mfrow = c(d, d), 
          mar = rep(0.2/2,4),
          oma = rep(4,4)+c(0,0,1*(!is.null(main)),0))
      for(i in seq(d))
      { 
        for(j in seq(d)) 
        { 
          if(i == j) 
          { 
            plot(if(nrow(newdata) == 0) data[,dimens[c(i,j)]] 
                 else                   newdata[,dimens[c(i,j)]],
                 type="n", xlab = "", ylab = "", axes=FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), 
                 labels = if(nrow(newdata) == 0) dataNames[dimens][i]
                          else                   newdataNames[dimens][i], 
                 cex = 1.5, adj = 0.5)
            box()
          } else 
          { 
            scatellipses(if(nrow(newdata) == 0) data else newdata, 
                         dimens = dimens[c(j,i)], 
                         nclass = nclass, 
                         symbols = symbols, colors = colors, 
                         xaxt = "n", yaxt = "n") 
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) 
                     if(nrow(newdata) == 0) "Training data" else "Test data"
                   else NULL, 
              cex.main = 1.2*oldpar$cex.main, 
              outer = TRUE, line = 3) 
    }
  }      
  
  plot.MclustDA.classification <- function(...)
  { 
    if(nrow(newdata) == 0 && d == 1)
    { 
      mclust1Dplot(data = data[,dimens[1],drop=FALSE], 
                   what = "classification",
                   classification = predClass[1:n], 
                   colors = colors[1:nclass],
                   xlab = dataNames[dimens],
                   main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Training data classification" 
              else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 && d == 2)
    { 
      coordProj(data = data[,dimens], what = "classification",
                classification = predClass[1:n], 
                main = FALSE, 
                colors = colors[1:nclass], 
                symbols = symbols[1:nclass], ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Training data classification" 
              else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 && d > 2)
    { 
      clPairs(data[,dimens], 
              classification = predClass[1:n],
              colors = colors[1:nclass], 
              symbols = symbols[1:nclass],
              cex.labels = 1.5,
              main = if(!is.null(main)) 
                        if(is.character(main)) main 
                        else if(as.logical(main)) "Training data classification" 
                             else NULL,
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 && d == 1)
    { 
      mclust1Dplot(data = newdata[,dimens], 
                   what = "classification",
                   classification = predClass[-(1:n)], 
                   xlab = newdataNames[dimens],
                   main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Test data classification" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 && d == 2)
    { 
      coordProj(data = newdata[,dimens], what ="classification",
                classification = predClass[-(1:n)], 
                colors = colors[1:nclass], 
                symbols = symbols[1:nclass],
                main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Test data classification" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 & length(dimens) > 2)
    { 
      on.exit(par(oldpar))
      # par(oma = c(0,0,10,0))      
      clPairs(data = newdata[,dimens], 
              classification = predClass[-(1:n)], 
              colors = colors[1:nclass], 
              symbols = symbols[1:nclass],
              cex.labels = 1.5, 
              main = if(!is.null(main))
                       if(is.character(main)) main 
                       else if(as.logical(main)) "Test data classification" 
                            else NULL, 
              cex.main = oldpar$cex.lab)
    }
  }

  plot.MclustDA.traintest <- function(...)
  { 
    cl <- factor(rep(c("Train","Test"), 
                     times = c(nrow(data), nrow(newdata))),
                 levels = c("Train", "Test"))

    if(d == 1)
    { 
      mclust1Dplot(data = Data[,dimens], what = "classification",
                   classification = cl,
                   xlab = dataNames[dimens], ylab = "",
                   colors = c("grey20", "grey80"),
                   main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Training and Test data" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(d == 2)
    { 
      coordProj(Data, dimens = dimens[1:2], what = "classification",
                classification = cl, CEX = 0.8,
                symbols = c(19,3), colors = c("grey80", "grey20"),
                main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Training (o) and Test (+) data"
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(d > 2)
    { 
      clPairs(Data[,dimens], classification = cl, 
              symbols = c(19,3), colors = c("grey80", "grey20"),
              main = if(!is.null(main))
                       if(is.character(main)) main 
                       else if(as.logical(main)) "Training (o) and Test (+) data"
                            else NULL,
              cex.main = oldpar$cex.lab)
    }
    
  }

  plot.MclustDA.error <- function(...)
  { 
    if(nrow(newdata) == 0 && d == 1)
    { 
      mclust1Dplot(data = data[,dimens], what = "error", 
                   classification = predClass[1:n], 
                   truth = trainClass, 
                   xlab = dataNames[dimens],
                   main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Training data error" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 && d == 2)
    { 
      coordProj(data = data[,dimens[1:2]], what = "error",
                classification = predClass[1:n], 
                truth = trainClass, main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Training data error" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) == 0 && d > 2)
    { 
      on.exit(par(oldpar))
      par(mfrow = c(d, d), 
          mar = rep(0.2/2,4),
          oma = rep(4,4)+c(0,0,1*(!is.null(main)),0))
      for(i in seq(d))
      { 
        for(j in seq(d)) 
        { 
          if(i == j) 
          { 
            plot(data[,dimens[c(i,j)]], type="n",
                 xlab = "", ylab = "", axes=FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), 
                 dataNames[dimens][i], cex = 1.5, adj = 0.5)
            box()
          } else 
          { 
            coordProj(data = data[,dimens[c(j,i)]], what = "error",
                      classification = predClass[1:n], 
                      truth = trainClass, main = FALSE,
                      xaxt = "n", yaxt = "n")
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }    
      
      if(!is.null(main))
        title(if(is.character(main)) main
              else if(as.logical(main)) "Training data error"
                   else NULL,
              cex.main = 1.2*oldpar$cex.main, 
              outer = TRUE, line = 3)
    }
    
    if(nrow(newdata) > 0 && d == 1)
    { 
      mclust1Dplot(data = newdata[,dimens], what = "error", 
                   classification = predClass[-(1:n)], 
                   truth = newclass, 
                   xlab = newdataNames[dimens],
                   main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Test data error" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }
    
    if(nrow(newdata) > 0 && d == 2)
    { 
      coordProj(data = newdata[,dimens[1:2]], what = "error",
                classification = predClass[-(1:n)], 
                truth = newclass, main = FALSE, ...)
      if(!is.null(main))
        title(if(is.character(main)) main 
              else if(as.logical(main)) "Test data error" 
                   else NULL, 
              cex.main = oldpar$cex.lab)
    }

    if(nrow(newdata) > 0 && d > 2)
    { 
      on.exit(par(oldpar))
      par(mfrow = c(d, d), 
          mar = rep(0.2/2,4),
          oma = rep(4,4)+c(0,0,1*(!is.null(main)),0))
      for(i in seq(d))
      { 
        for(j in seq(d)) 
        { 
          if(i == j) 
          { 
            plot(newdata[,dimens[c(i,j)]], type="n",
                 xlab = "", ylab = "", axes=FALSE)
            text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), 
                 newdataNames[dimens][i], cex = 1.5, adj = 0.5)
            box()
          } else 
          { 
            coordProj(data = newdata[,dimens[c(j,i)]], what = "error",
                      classification = predClass[-(1:n)], 
                      truth = newclass, main = FALSE,
                      xaxt = "n", yaxt = "n")
          }
          if(i == 1 && (!(j%%2))) axis(3)
          if(i == d && (j%%2))    axis(1)
          if(j == 1 && (!(i%%2))) axis(2)
          if(j == d && (i%%2))    axis(4)
        }
      }    
      
      if(!is.null(main))
        title(if(is.character(main)) main
              else if(as.logical(main)) "Test data error"
                   else NULL,
              cex.main = 1.2*oldpar$cex.main, 
              outer = TRUE, line = 3)
    }

        
  }

  if(interactive() & length(what) > 1)
    { title <- "Model-based discriminant analysis plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "scatterplot")    plot.MclustDA.scatterplot(...)
             if(what[choice] == "classification") plot.MclustDA.classification(...)
             if(what[choice] == "train&test")     plot.MclustDA.traintest(...)
             if(what[choice] == "error")          plot.MclustDA.error(...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  }
  else 
    { if(any(what == "scatterplot"))    plot.MclustDA.scatterplot(...)
      if(any(what == "classification")) plot.MclustDA.classification(...)
      if(any(what == "train&test"))     plot.MclustDA.traintest(...) 
      if(any(what == "error"))          plot.MclustDA.error(...)
  }
    
  invisible()
  
}


cvMclustDA <- function(object, nfold = 10, 
                       metric = c("error", "brier"), 
                       prop = object$prop,
                       verbose = interactive(), ...) 
{

  call <- object$call
  nfold <- as.numeric(nfold)
  metric <- match.arg(metric, 
                      choices = eval(formals(cvMclustDA)$metric), 
                      several.ok = FALSE)
  #
  data <- object$data
  class <- as.factor(object$class)
  n <- length(class)
  G <- lapply(object$models, function(mod) mod$G)
  modelName <- lapply(object$models, function(mod) mod$modelName)
  #
  ce <- function(pred, class)
  {
    1 - sum(class == pred, na.rm = TRUE)/length(class)
  }
  #
  folds <- if(nfold == n) lapply(1:n, function(x) x)
           else           balanced.folds(class, nfolds = nfold)
  nfold <- length(folds)
  folds.size <- sapply(folds, length)
  #
  cvmetric <- rep(NA, nfold)
  cvclass <- factor(rep(NA, n), levels = levels(class))
  cvprob  <- matrix(as.double(NA), nrow = n, ncol = nlevels(class),
                    dimnames = list(NULL, levels(class)))
  
  if(verbose)
  { 
    cat("cross-validating ...\n")
    flush.console()
    pbar <- txtProgressBar(min = 0, max = nfold, style = 3)
    on.exit(close(pbar))
  }
  
  for(i in seq(nfold))
  { 
    x <- data[-folds[[i]],,drop=FALSE]
    y <- class[-folds[[i]]]
    call$data <- x
    call$class <- y
    call$G <- G
    call$modelNames <- modelName
    call$verbose <- FALSE
    mod <- eval(call, parent.frame())
    #
    predTest <- predict(mod, data[folds[[i]],,drop=FALSE], prop = prop)
    cvmetric[i] <- if(metric == "error") 
                     ce(predTest$classification, class[folds[[i]]])
                   else 
                     BrierScore(predTest$z, class[folds[[i]]])
    cvclass[folds[[i]]] <- predTest$classification
    cvprob[folds[[i]],] <- predTest$z
    #
    if(verbose) 
      setTxtProgressBar(pbar, i)
  }
  #    
  cv <- sum(cvmetric*folds.size)/sum(folds.size)
  se <- sqrt(var(cvmetric)/nfold)
  #
  out <- list(classification = cvclass, 
              z = cvprob,
              error = if(metric == "error") cv else NA,
              brier = if(metric == "brier") cv else NA,
              se = se)
  return(out)
}


balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) 
{ 
  # Create 'nfolds' balanced folds conditional on grouping variable 'y'.
  # Function useful in evaluating a classifier by balanced cross-validation.
  # Returns a list with 'nfolds' elements containing indexes of each fold.
  # 
  # From package 'pamr' by T. Hastie, R. Tibshirani, Balasubramanian 
  # Narasimhan, Gil Chu.
  
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)     # nice we to get the ids in a list, split by class
  ### create a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(as.double(NA), ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) 
    #   { bigmat[seq(totals[i]), i] <- sample(yids[[i]]) }
    #   Luca: this version has a bug if a class has only 1 obs
  { if (totals[i]==1) 
    bigmat[seq(totals[i]), i] <- yids[[i]]    
    else 
      bigmat[seq(totals[i]), i] <- sample(yids[[i]]) }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat))   
  res <-vector("list", nfolds)
  for(j in 1:nfolds) 
  { jj <- !is.na(smallmat[, j])
    res[[j]] <- smallmat[jj, j] }
  return(res)
}

permute.rows <- function(x)
{
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

# Deprecated functions

cv1EMtrain <- function(data, labels, modelNames=NULL) 
{
  .Deprecated("cvMclustDA", package = "mclust")
  z <- unmap(as.numeric(labels))
  G <- ncol(z)
  dimDataset <- dim(data)
  oneD <- is.null(dimDataset) || length(dimDataset[dimDataset > 1]) == 1
  if (oneD || length(dimDataset) != 2) {
    if (is.null(modelNames)) 
      modelNames <- c("E", "V")
    if (any(!match(modelNames, c("E", "V"), nomatch = 0))) 
      stop("modelNames E or V for one-dimensional data")
    n <- length(data)
    cv <- matrix(1, nrow = n, ncol = length(modelNames))
    dimnames(cv) <- list(NULL, modelNames)
    for (m in modelNames) {
      for (i in 1:n) {
        mStep <- mstep(modelName = m, data = data[-i], 
                       z = z[-i,], warn = FALSE)
        eStep <- do.call("estep", c(mStep, list(data = data[i], 
                                                warn = FALSE)))
        if (is.null(attr(eStep, "warn"))) {
          k <- (1:G)[eStep$z == max(eStep$z)]
          l <- (1:G)[z[i,] == max(z[i,])]
          cv[i, m] <- as.numeric(!any(k == l))
        }
      }
    }
  }
  else {
    if (is.null(modelNames)) 
      modelNames <- mclust.options("emModelNames")
    n <- nrow(data)
    cv <- matrix(1, nrow = n, ncol = length(modelNames))
    dimnames(cv) <- list(NULL, modelNames)
    for (m in modelNames) {
      for (i in 1:n) {
        mStep <- mstep(modelName = m, data = data[-i,],
                       z = z[-i,], warn = FALSE)
        eStep <- do.call("estep", c(mStep, list(data = data[i, 
                                                            , drop = FALSE], warn = FALSE)))
        if (is.null(attr(eStep, "warn"))) {
          k <- (1:G)[eStep$z == max(eStep$z)]
          l <- (1:G)[z[i,] == max(z[i,])]
          cv[i, m] <- as.numeric(!any(k == l))
        }
      }
    }
  }
  errorRate <- apply(cv, 2, sum)
  errorRate/n
}

bicEMtrain <- function(data, labels, modelNames=NULL) 
{
  .Deprecated("MclustDA", package = "mclust")
  
  z <- unmap(as.numeric(labels))
  G <- ncol(z)
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 
    1
  if (oneD || length(dimData) != 2) {
    if (is.null(modelNames)) 
      modelNames <- c("E", "V")
    if (any(!match(modelNames, c("E", "V"), nomatch = 0))) 
      stop("modelNames E or V for one-dimensional data")
  }
  else {
    if (is.null(modelNames)) 
      modelNames <- mclust.options("emModelNames")
  }
  BIC <- rep(NA, length(modelNames))
  names(BIC) <- modelNames
  for (m in modelNames) {
    mStep <- mstep(modelName = m, data = data, z = z, warn = FALSE)
    eStep <- do.call("estep", c(mStep, list(data=data, warn=FALSE)))
    if (is.null(attr(eStep, "warn"))) 
      BIC[m] <- do.call("bic", eStep)
  }
  BIC
}

cv.MclustDA <- function(...) 
{
  .Deprecated("cvMclustDA", package = "mclust")
  cvMclustDA(...)
}

# "[.mclustDAtest" <- function (x, i, j, drop = FALSE) 
# {
#   clx <- oldClass(x)
#   oldClass(x) <- NULL
#   NextMethod("[")
# }

classPriorProbs <- function(object, newdata = object$data, 
                            itmax = 1e3, eps = sqrt(.Machine$double.eps))
{

  if(!inherits(object, "MclustDA")) 
    stop("object not of class \"MclustDA\"")

  z <- predict(object, newdata = newdata)$z
  prop <- object$prop
  p <- colMeans(z)
  p0 <- p+1
  it <- 0
  while(max(abs(p-p0)/abs(p)) > eps & it < itmax)
  {
    it <- it+1
    p0 <- p
    # z_upd <- t(apply(z, 1, function(z) { z <- z*p/prop; z/sum(z) }))
    z_upd <- sweep(z, 2, FUN = "*", STATS = p/prop)
    z_upd <- sweep(z_upd, MARGIN = 1, FUN = "/", STATS = rowSums(z_upd))
    p <- colMeans(z_upd)
  }
  return(p)
}
