# Example
#
# library(mclust)
# data(faithful)
# out = mclustCriteria(faithful)
# print(out)
# plot(out)

icl <- function(object, ...) UseMethod("icl")

icl.Mclust <- function(object, ...)
{
  n <- object$n
  # G <- object$G + ifelse(is.na(object$hypvol),0,1)
  z <- object$z
  if(is.null(z)) z <- matrix(1, nrow = n, ncol = 1)
  C <- matrix(0, n, ncol(z))
  for(i in 1:n) 
    C[i, which.max(z[i,])] <- 1
  object$bic + 2*sum(C * ifelse(z > 0, log(z), 0))
}

icl.MclustDA <- function(object, ...)
{
  n <- object$n
  z <- predict(object)$z
  df <- object$df
  if(is.null(z)) z <- matrix(1, nrow = n, ncol = 1)
  C <- matrix(0, n, ncol(z))
  for(i in 1:n) 
    C[i, which.max(z[i,])] <- 1
  object$bic + 2*sum(C * ifelse(z > 0, log(z), 0))
}

mclustICL <- function(data, G = NULL, modelNames = NULL, 
                      initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),  
                      ...)
{
  data <- data.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  
  if(is.null(modelNames))
    { if(d == 1)
        { modelNames <- c("E", "V") }
      else
        { modelNames <- mclust.options("emModelNames") }
  }
  if(n <= d) 
    { # select only spherical and diagonal models
      m <- match(modelNames, c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
                 nomatch = 0)
      modelNames <- modelNames[m]
  }
  
  if(is.null(G)) 
    { G <- if(is.null(initialization$noise)) 1:9 else 0:9 }
  else 
    { G <- sort(as.integer(G)) }
  
  if(is.null(initialization$hcPairs) & (d > 1))
    { if(n > d)
        { initialization$hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], 
                                       data = data) }
      else 
    { initialization$hcPairs <- hc(modelName = "EII", data = data) }
  }
  
  ICL <- matrix(as.double(NA), nrow = length(G), ncol = length(modelNames), 
                dimnames = list(as.character(G), as.character(modelNames)))
  warn <- getOption("warn")
  on.exit(options("warn" = warn))
  options("warn" = -1)
  for(i in 1:nrow(ICL))
     { for(j in 1:ncol(ICL))
          { mod <- try(Mclust(data, G = G[i], modelNames = modelNames[j], 
                              initialization = initialization, ...),
                       silent = TRUE)
            if(class(mod) == "try-error") next()
            if(all(is.na(mod$BIC))) next()
            ICL[i,j] <- icl(mod)
       }
  }
  ICL <- structure(ICL, criterion = "ICL", 
                   class = "mclustICL")
  return(ICL)
}

print.mclustICL <- function (x, pick = 3, ...) 
{
  oldClass(x) <- attr(x, "args") <- NULL
  cat("Integrated Complete Likelihood (ICL) criterion:\n")
  NextMethod("print")
  cat("\n")
  cat("Top", pick, "models based on the ICL criterion:\n")
  print(pickBIC(x, pick), ...)
  invisible()
}

summary.mclustICL <- function(object, G, modelNames, ...)
{
  if(!missing(G)) 
    object <- object[rownames(object) %in% G,,drop=FALSE]
  if(!missing(modelNames)) 
    object <- object[,colnames(object) %in% modelNames,drop=FALSE]
  structure(pickBIC(object, ...),
            class = "summary.mclustICL")
}

print.summary.mclustICL <- function(x, digits = getOption("digits"), ...)
{
  cat("Best ICL values:\n")
  x <- drop(as.matrix(x))
  x <- rbind(ICL = x, "ICL diff" = x - max(x))
  print(x, digits = digits)
  invisible()
}


plot.mclustICL <- function(x, ylab = "ICL", ...) 
{
  plot.mclustBIC(x, ylab = ylab, ...)  
}

