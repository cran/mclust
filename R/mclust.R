### mclust for R, version 1.1-2: some changes in graphics functions
###     and completion of help files
### Changed dd. Jan. 10th, 2001
### .   made xlim, ylim, xlab and ylab arguments in mixproj
### .   use matplot in plot.emclust, with some extra arguments
### Changed dd. jan 22, 2001
### .   made xlab and ylab arguments in emclust and emclust1
### .   removed function charconv since it was only used in function partuniq
### .   removed function print.mclust since there is no class mclust
###
### mclust for R, version 1.1-1: first working version
### Change dd. July 28th, 2000 (S-incompatibilities?)
### .   changed sigmasq in estep.EI, estep.VI, mstep.EI and mstep.VI to sigma;
###     apparently S can complete names
### .   added if (!is.na(bicval)) to RC statements in emclust and emclust1
###     apparently S does not mind...
### .   added [1:3] after "cols" in summary.emclust in statement 
###           names(best) <- names(rcond) <- paste 
### .   replaced -.Machine$double.xmax with -Inf in summary.emclust

"awe" <- function(tree, data)
{
  data <- as.matrix(data)
  p <- ncol(data)
  n <- nrow(data)
  dof <- switch(attr(tree, "model"),
		EI = p,
		VI = p + 1,
		EEE = p,
		VVV = (p * (p - 1))/2 + 2 * p,
		EEV = (p * (p - 1))/2 + p,
		VEV = (p * (p - 1))/2 + p + 1,
		stop("invalid model id"))
  like <- loglik(tree, data)
  class(tree) <- NULL
  u <- attr(attr(tree, "initial.partition"), "unique")
  s <- ncol(tree)
  nmerge <- attr(like, "nmerge")
  attr(like, "nmerge") <- NULL
  if(!all(good <- !is.na(like))) {
    like <- like[good]
    l <- length(like)
    if(l <= 1)
      return(rep(NA, n - 1))
    nmerge <- nmerge[ - ((1:length(good))[!good])]
  }
  AWE <- -2 * diff(like) - (3 + 2 * log(p * nmerge)) * dof
  c(0, rep(NA, u - s - 1), cumsum(rev(AWE)), 
    rep(NA, (n - 1) - length(AWE) + (u - s)))
}

"bic" <- function(data, modelid, ...)
{
### ... z, eps, equal = F, noise = F, Vinv
  switch(as.character(modelid),
	 EI = bic.EI(data, ...),
	 VI = bic.VI(data, ...),
	 EEE = bic.EEE(data, ...),
	 VVV = bic.VVV(data, ...),
	 EEV = bic.EEV(data, ...),
	 VEV = bic.VEV(data, ...),
	 stop("invalid model id"))
}

"bic.EEE" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(z)) {
### one cluster case
    if(noise) {
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      loglik <- n * log(Vinv)
      bic <- 2 * loglik - log(n)
      attr(bic, "params") <- list(Vinv = Vinv)
      attr(bic, "loglik") <- loglik
    }
    else {
      nparams <- (p * (p + 1))/2
      temp <- one.XXX(data)
      loglik <- attr(temp, "loglik")
      rcond <- attr(temp, "rcond")
      attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
      bic <- 2 * loglik - (p + nparams) * log(n)
      attr(bic, "params") <- temp
      attr(bic, "rcond") <- rcond
      attr(bic, "loglik") <- loglik
    }
  }
  else {
    if(!any(is.na(z))) {
      K <- ncol(z)
      nparams <- (p * (p + 1))/2
      if(noise) {
	G <- K - 1
	if(missing(Vinv))
	  Vinv <- hypvol(data, reciprocal = T)
	temp <- mstep.EEE(data, z, eps = eps, equal = 
			  equal, noise = T, Vinv = Vinv)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- NULL	
	## plus one parameter for 1/V; plus p - 1 parameters for proportions
	if(equal) {
	  bic <- 2 * loglik - (G * p + nparams + 1) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * p + G + nparams + 1) * log(n)
	}
	attr(bic, "params") <- c(temp, Vinv = Vinv)
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
      else {
	G <- K
	temp <- mstep.EEE(data, z, eps = eps, equal = equal)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- NULL	
	## plus p - 1 parameters for proportions
	if(equal) {
	  bic <- 2 * loglik - (G * p + nparams) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * p + (G - 1) + nparams) * log(n)
	}
	attr(bic, "params") <- temp
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
    }
    else {
      bic <- NA
    }
  }
  attr(bic, "model") <- "EEE"
  attr(bic, "class") <- "bic"
  bic
}

"bic.EEV" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- c(.Machine$double.eps, sqrt(.Machine$double.eps))
  else if(length(eps) == 1)
    eps <- c(eps, sqrt(.Machine$double.eps))
  if(missing(z)) {
### one cluster case
    if(noise) {
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      loglik <- n * log(Vinv)
      bic <- 2 * loglik - log(n)
      attr(bic, "params") <- list(Vinv = Vinv)
      attr(bic, "loglik") <- loglik
    }
    else {
      nparams <- (p * (p - 1))/2
      temp <- one.XXX(data)
      loglik <- attr(temp, "loglik")
      rcond <- attr(temp, "rcond")
      attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
      bic <- 2 * loglik - (p + nparams + (p - 1) + 1) * log(n)
      attr(bic, "params") <- temp
      attr(bic, "rcond") <- rcond
      attr(bic, "loglik") <- loglik
    }
  }
  else {
    if(!any(is.na(z))) {
      K <- ncol(z)
      nparams <- (p * (p - 1))/2
      if(noise) {
	G <- K - 1
	if(missing(Vinv))
	  Vinv <- hypvol(data, reciprocal = T)
	temp <- mstep.EEV(data, z, eps = eps, equal = 
			  equal, noise = T, Vinv = Vinv)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + nparams) + (p - 1) + 1 + 1) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * (p+nparams) + G + (p-1) + 1 + 1) * log(n) 
	}
	attr(bic, "params") <- c(temp, list(Vinv = Vinv))
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
      else {
	G <- K
	temp <- mstep.EEV(data, z, eps = eps, equal = equal)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + nparams) + (p - 1) + 1) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * (p+nparams) + (G-1) + (p-1) + 1) * log(n)}
	attr(bic, "params") <- temp
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
    }
    else {
      bic <- NA
    }
  }
  attr(bic, "model") <- "EEV"
  attr(bic, "class") <- "bic"
  bic
}

"bic.EI" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(z)) {
### one cluster case
    if(noise) {
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      loglik <- n * log(Vinv)
      bic <- 2 * loglik - log(n)
      attr(bic, "params") <- list(Vinv = Vinv)
      attr(bic, "loglik") <- loglik
    }
    else {
      temp <- one.XI(data)
      loglik <- attr(temp, "loglik")
      rcond <- attr(temp, "rcond")
      attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
      bic <- 2 * loglik - (p + 1) * log(n)
      attr(bic, "params") <- temp
      attr(bic, "rcond") <- rcond
      attr(bic, "loglik") <- loglik
    }
  }
  else {
    if(!any(is.na(z))) {
      K <- ncol(z)
      if(noise) {
	G <- K - 1
	if(missing(Vinv))
	  Vinv <- hypvol(data, reciprocal = T)
	temp <- mstep.EI(data, z, eps = eps, equal = 
			 equal, noise = T, Vinv = Vinv)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
	if(equal) {
	  bic <- 2 * loglik - (G * p + 1 + 1) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * p + G + 1 + 1) * log(
							n)
	}
	attr(bic, "params") <- c(temp, list(Vinv = Vinv
					    ))
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
      else {
	G <- K
	temp <- mstep.EI(data, z, eps = eps, equal = 
			 equal)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik - (G * p + 1) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * p + (G - 1) + 1) * 
	    log(n)
	}
	attr(bic, "params") <- temp
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
    }
    else {
      bic <- NA
    }
  }
  attr(bic, "model") <- "EI"
  class(bic) <- "bic"
  bic
}

"bic.VEV" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- c(.Machine$double.eps, .Machine$double.eps)
  else if(length(eps) == 1)
    eps <- c(eps, .Machine$double.eps)
  if(missing(z)) {
					## one cluster case
    if(noise) {
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      loglik <- n * log(Vinv)
      bic <- 2 * loglik - log(n)
      attr(bic, "params") <- list(Vinv = Vinv)
      attr(bic, "loglik") <- loglik
    }
    else {
      ## p*p - p(p+1)/2 for orientation, 1 for volume
      nparams <- (p * (p - 1))/2 + 1
      temp <- one.XXX(data)
      loglik <- attr(temp, "loglik")
      rcond <- attr(temp, "rcond")
      attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
      bic <- 2 * loglik - (p + nparams + (p - 1)) * log(n)
      attr(bic, "params") <- temp
      attr(bic, "rcond") <- rcond
      attr(bic, "loglik") <- loglik
    }
  }
  else {
    if(!any(is.na(z))) {
      K <- ncol(z)	
      ## p*p - p(p+1)/2 for orientation, 1 for volume
      nparams <- (p * (p - 1))/2 + 1
      if(noise) {
	G <- K - 1
	if(missing(Vinv))
	  Vinv <- hypvol(data, reciprocal = T)
	temp <- mstep.VEV(data, z, eps = eps, equal = 
			  equal, noise = T, Vinv = Vinv)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + nparams) + (p - 1) + 1) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * (p + nparams) + G + (p - 1) + 1) * log(n)
	}
	attr(bic, "params") <- c(temp, list(Vinv = Vinv))
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
      else {
	G <- K
	temp <- mstep.VEV(data, z, eps = eps, equal = 
			  equal)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + nparams) + (p - 1)) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * (p + nparams) + (G - 1) + (p - 1)) * log(n)
	}
	attr(bic, "params") <- temp
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
    }
    else {
      bic <- NA
    }
  }
  attr(bic, "model") <- "VEV"
  attr(bic, "class") <- "bic"
  bic
}

"bic.VI" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(z)) {
					## one cluster case
    if(noise) {
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      loglik <- n * log(Vinv)
      bic <- 2 * loglik - log(n)
      attr(bic, "params") <- list(Vinv = Vinv)
      attr(bic, "loglik") <- loglik
    }
    else {
      temp <- one.XI(data)
      loglik <- attr(temp, "loglik")
      rcond <- attr(temp, "rcond")
      attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
      bic <- 2 * loglik - (p + 1) * log(n)
      attr(bic, "params") <- temp
      attr(bic, "rcond") <- rcond
      attr(bic, "loglik") <- loglik
    }
  }
  else {
    if(!any(is.na(z))) {
      K <- ncol(z)
      if(noise) {
	G <- K - 1
	if(missing(Vinv))
	  Vinv <- hypvol(data, reciprocal = T)
	temp <- mstep.VI(data, z, eps = eps, equal = 
			 equal, noise = T, Vinv = Vinv)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + 1) + 1) * log(n
						      )
	}
	else {
	  bic <- 2 * loglik - (G * (p + 1) + G + 1) * 
	    log(n)
	}
	attr(bic, "params") <- c(temp, list(Vinv = Vinv
					    ))
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
      else {
	G <- K
	temp <- mstep.VI(data, z, eps = eps, equal = 
			 equal)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + 1)) * log(n)
	}
	else {
	  bic <- 2 * loglik - (G * (p + 1) + (G - 1)) * 
	    log(n)
	}
	attr(bic, "params") <- temp
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
    }
    else {
      bic <- NA
    }
  }
  attr(bic, "model") <- "VI"
  class(bic) <- "bic"
  bic
}
"bic.VVV" <-
function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(z)) {
					## one cluster case
    if(noise) {
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      loglik <- n * log(Vinv)
      bic <- 2 * loglik - log(n)
      attr(bic, "params") <- list(Vinv = Vinv)
      attr(bic, "loglik") <- loglik
    }
    else {
      nparams <- (p * (p + 1))/2
      temp <- one.XXX(data)
      loglik <- attr(temp, "loglik")
      rcond <- attr(temp, "rcond")
      attr(temp, "loglik") <- attr(temp, "rcond") <- NULL
      bic <- 2 * loglik - (p + nparams) * log(n)
      attr(bic, "params") <- temp
      attr(bic, "rcond") <- rcond
      attr(bic, "loglik") <- loglik
    }
  }
  else {
    if(!any(is.na(z))) {
      K <- ncol(z)
      nparams <- (p * (p + 1))/2
      if(noise) {
	G <- K - 1
	if(missing(Vinv))
	  Vinv <- hypvol(data, reciprocal = T)
	temp <- mstep.VVV(data, z, eps = eps, equal = 
			  equal, noise = T, Vinv = Vinv)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik - (G * (p + nparams) + 1) * 
	    log(n)
	}
	else {
	  bic <- 2 * loglik - (G * (p + nparams) + G + 
			       1) * log(n)
	}
	attr(bic, "params") <- c(temp, list(Vinv = Vinv
					    ))
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
      else {
	G <- K
	temp <- mstep.VVV(data, z, eps = eps, equal = 
			  equal)
	loglik <- attr(temp, "loglik")
	rcond <- attr(temp, "rcond")
	attr(temp, "loglik") <- attr(temp, "rcond") <- 
	  NULL
	if(equal) {
	  bic <- 2 * loglik -
            (G * (p + nparams)) * log(n)
	}
	else {
	  bic <- 2 * loglik -
            (G * (p + nparams) + (G - 1)) * log(n)
	}
	attr(bic, "params") <- temp
	attr(bic, "rcond") <- rcond
	attr(bic, "loglik") <- loglik
      }
    }
    else {
      bic <- NA
    }
  }
  attr(bic, "model") <- "VVV"
  attr(bic, "class") <- "bic"
  bic
}

"censcale" <- function(x, tol = 0.0001)
{
  x <- as.matrix(x)
  mu <- apply(x, 2, mean)
  x <- sweep(x, 2, mu)
  sd <- apply(x, 2, function(z)
	      sqrt(var(z)))
  sd[sd <= tol] <- 1
  structure(sweep(x, 2, sd, "/"), mu = mu, sd = sd)
}

#"charconv" <- function(x, sep = "001")
#{
#  if(!is.data.frame(x))
#    x <- data.frame(x)
#  do.call("paste", c(as.list(x), sep = sep))
#}

"clpairs" <- function(x, partition, col, ...) 
{
  x <- as.matrix(x)
  m <- nrow(x)
  n <- ncol(x)
  if(missing(partition))
    partition <- rep(1, m)
  l <- length(unique(partition))
  if(missing(col))
    col <- partition
  else if (length(unique(col)) < l)
    stop("more colors needed")

  pairs(x, col=col, ...)
  invisible()
}

"ctoz" <- function(cl, noise)
{
  ## converts a classification to conditional probabilities
  ## classes are arranged in sorted order
  ## if a noise indicator is specified, that column is placed last
  n <- length(cl)
  u <- sort(unique(cl))
  labs <- as.character(u)
  k <- length(u)
  z <- matrix(0, n, k)
  if(!missing(noise)) {
    l <- u == noise
    if(any(l)) {
      m <- max(u) + 1
      u[l] <- m
      labs <- labs[order(u)]
      u <- sort(u)
      cl[cl == noise] <- m
    }
  }
  for(j in 1:k)
    z[cl == u[j], j] <- 1
  dimnames(z) <- list(NULL, labs)
  if(any((sumz <- as.integer(apply(z, 1, sum))) != 1)) {
    warning("improper z")
    stop("STOP")
  }
  z
}

"emclust" <- function(data, nclus, modelid, k, equal = F, noise, Vinv)
{
  dd <- dim(data)
  d <- dd[dd != 1]
  if(length(d) != 2)
    stop("data must be matrix-conformal")
  data <- if(length(d) != length(dd)) matrix(data, d[1], d[2]) else 
  as.matrix(data)
  n <- d[1]
  if(missing(modelid))
    modelid <- c("EI", "VI", "EEE", "VVV", "EEV", "VEV")
  if(missing(k)) {

    ## use all of the data in the initial hierarchical clustering phase

    if(missing(noise)) {
      if(missing(nclus))
	nclus <- 1:9
      nclus <- sort(nclus)
      l <- length(nclus)	
###---------------------------------------------------------------------------
      tree <- mhtree.VVV(data)
      clss <- mhclass(tree, nclus)	
###---------------------------------------------------------------------------
      gauss <- function(modelid, data, clss, nclus, l, equal
			)
	{
	  bicnam <- paste("bic.", modelid, sep = "")
	  menam <- paste("me.", modelid, sep = "")
	  BIC <- RC <- numeric(l)
	  for(j in 1:l) {
	    i <- nclus[j]
	    if(i == 1) {
	      bicval <- do.call(bicnam, list(data))
	    }
	    else {
	      z <- do.call(menam, list(data, ctoz(clss[, as.character(i)]),
				       equal = equal))
	      bicval <- do.call(bicnam, list(data, z, equal = equal))
	    }
	    BIC[j] <- bicval
	    if (!is.na(bicval)) RC[j] <- attr(bicval, "rcond")
	  }
	  list(bic = BIC, rc = RC)
	}
      all <- lapply(as.list(modelid), gauss, data = data, 
		    clss = clss, nclus = nclus, l = l, equal = 
		    equal)
      all <- list(bic = t(sapply(all, function(z)
		    z$bic)), rc = t(sapply(all, function(z)
			       z$rc)))	
      ##
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      dimnames(all$bic) <- dimnames(all$rc) <- list(modelid, 
						    as.character(nclus))
      structure(all$bic, equal = equal, tree = tree, rcond = 
		all$rc, class = "emclust")
    }
    else {
      ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ## noise case
      ##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      noise <- as.logical(noise)
      n <- nrow(data)	
      ##
      ##----------------------------------------------------------------------
      if(missing(nclus))
	nclus <- 0:9
      nclus <- sort(nclus)
      l <- length(nclus)	
      ##-----------------------------------------------------------------------
      tree <- mhtree.VVV(data[!noise,  ])
      clss <- mhclass(tree, nclus[nclus != 0])
      if(missing(Vinv)) Vinv <- hypvol(data, reciprocal = T)	
      ##-----------------------------------------------------------------------
      gaussn <- function(modelid, data, clss, nclus, l, n, 
			 equal, noise, Vinv)
	{
	  bicnam <- paste("bic.", modelid, sep = "")
	  menam <- paste("me.", modelid, sep = "")
	  BIC <- RC <- numeric(l)
	  if(n != length(noise))
	    stop("STOP")
	  k <- 1
	  if(nclus[1] == 0) {
					# all noise --- same for all models
	    BIC[k] <- bic(data, modelid = modelid, equal
			  = equal, noise = T, Vinv = Vinv)
	    RC[k] <- NA
	    k <- k + 1
	  }
	  if(k < l) {
	    cl <- numeric(n)
	    for(j in k:l) {
	      i <- nclus[j]
	      cl[!noise] <- clss[, as.character(i)]
	      cl[noise] <- i + 1
	      z <- do.call(menam, 
			   list(data, ctoz(cl), 
				equal = equal, noise = T, Vinv = Vinv))
	      bicval <- do.call(bicnam, 
				list(data, z, equal = equal, noise = T,
				     Vinv = Vinv)) 
	      BIC[j] <- bicval
	      if (!is.na(bicval)) RC[j] <- attr(bicval, "rc")
	    }
	  }
	  list(bic = BIC, rc = RC)
	}
      all <- lapply(as.list(modelid), gaussn, data = data, 
		    clss = clss, nclus = nclus, l = l, n = n, equal
		    = equal, noise = noise, Vinv = Vinv)
      all <- list(bic = t(sapply(all, function(z)
		    z$bic)), rc = t(sapply(all, function(z)
			       z$rc)))	
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      dimnames(all$bic) <- dimnames(all$rc) <- list(modelid, 
						    as.character(nclus))
      structure(all$bic, equal = equal, tree = tree, noise = 
		noise, Vinv = Vinv, rcond = all$rc, class = 
		"emclust")
    }
  }
  else {
############################################################################
### use only a sample of the data in the initial hierarchical clustering phase
############################################################################
    if(missing(noise)) {
      if(missing(nclus))
	nclus <- 1:9
      nclus <- sort(nclus)
      l <- length(nclus)	
      ##-----------------------------------------------------------------------
      smpl <- sample(1:n, size = k)
      tree <- mhtree.VVV(data[smpl,  ])
      clss <- mhclass(tree, nclus)	
      ##-----------------------------------------------------------------------
      gaussk <- function(modelid, data, clss, nclus, l, smpl, 
			 equal)
	{
	  bicnam <- paste("bic.", modelid, sep = "")
	  menam <- paste("me.", modelid, sep = "")
	  msnam <- paste("mstep.", modelid, sep = "")
	  esnam <- paste("estep.", switch(modelid,
					  EI = "EI",
					  VI = "VI",
					  EEE = "EEE",
					  VVV = "VVV",
					  EEV = "XEV",
					  VEV = "XEV"), sep = "")
	  BIC <- RC <- numeric(l)
	  for(j in 1:l) {
	    i <- nclus[j]
	    if(i == 1) {
	      bicval <- do.call(bicnam, list(data))
	    }
	    else {
	      pars <- do.call(msnam, list(data[smpl,  ], 
					  ctoz(clss[, as.character(i)]), 
					  equal = equal))
	      z <- do.call(esnam, c(list(data), pars))
	      z <- do.call(menam, list(data, z, equal = 
				       equal))
	      bicval <- do.call(bicnam, list(data, z, 
					     equal = equal))
	    }
	    BIC[j] <- bicval
	    if (!is.na(bicval)) RC[j] <- attr(bicval, "rcond")
	  }
	  list(bic = BIC, rc = RC)
	}
      all <- lapply(as.list(modelid), gaussk, data = data, 
		    clss = clss, nclus = nclus, l = l, smpl = smpl, 
		    equal = equal)
      all <- list(bic = t(sapply(all, function(z)
		    z$bic)), rc = t(sapply(all, function(z)
			       z$rc)))	
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      dimnames(all$bic) <- dimnames(all$rc) <- list(modelid, 
						    as.character(nclus))
      structure(all$bic, equal = equal, tree = tree, subset
		= smpl, rcond = all$rc, class = "emclust")
    }
    else 
      {
	##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	## noise case
	##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      noise <- as.logical(noise)
      m <- sum(as.numeric(!noise))
      n <- nrow(data)	
      ##-----------------------------------------------------------------------
      if(missing(nclus))
	nclus <- 0:9
      nclus <- sort(nclus)
      l <- length(nclus)	
      ##-----------------------------------------------------------------------
      smpl <- sample(1:m, size = k)
      tree <- mhtree.VVV(data[!noise,  ][smpl,  ])
      clss <- mhclass(tree, nclus[nclus != 0])
      if(missing(Vinv)) Vinv <- hypvol(data, reciprocal = T)	
      ##-----------------------------------------------------------------------
      gaussnk <- function(modelid, data, clss, nclus, l, n, 
			  smpl, equal, noise, Vinv)
	{
	  bicnam <- paste("bic.", modelid, sep = "")
	  menam <- paste("me.", modelid, sep = "")
	  msnam <- paste("mstep.", modelid, sep = "")
	  esnam <- paste("estep.", switch(modelid,
					  EI = "EI",
					  VI = "VI",
					  EEE = "EEE",
					  VVV = "VVV",
					  EEV = "XEV",
					  VEV = "XEV"), sep = "")
	  BIC <- RC <- numeric(l)
	  if(n != length(noise))
	    stop("STOP")
	  k <- 1
	  if(nclus[1] == 0) {
					# all noise --- same for all models
	    BIC[k] <- bic(data, modelid = modelid, noise
			  = T, Vinv = Vinv)
	    RC[k] <- NA
	    k <- k + 1
	  }
	  if(k < l) {
	    z <- matrix(0, n, nclus[l] + 1)
	    for(j in k:l) {
	      i <- nclus[j]
	      z[, 1:(i + 1)] <- 0
	      pars <- do.call(msnam, 
			      list(data[!noise,][smpl,], 
				   ctoz(clss[, as.character(i)]), 
				   equal = equal))
	      z[!noise, 1:i] <- do.call(esnam, c(list(data[!noise,]), pars))
	      z[noise, i + 1] <- 1
	      z[, 1:(i + 1)] <- 
		do.call(menam, list(data, z[, 1:(i + 1)], equal = equal,
				    noise = T, Vinv = Vinv)) 
	      bicval <- do.call(bicnam, list(data, z[, 1:(i + 1)],
					     equal = equal, noise = T, 
					     Vinv = Vinv)) 
	      BIC[j] <- bicval
	      if (!is.na(bicval)) RC[j] <- attr(bicval, "rc")
	    }
	  }
	  list(bic = BIC, rc = RC)
	}
      all <- lapply(as.list(modelid), gaussnk, data = data, 
		    clss = clss, nclus = nclus, l = l, n = n, smpl
		    = smpl, equal = equal, noise = noise, Vinv = 
		    Vinv)
      all <- list(bic = t(sapply(all, function(z)
		    z$bic)), rc = t(sapply(all, function(z)
			       z$rc)))	
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      dimnames(all$bic) <- dimnames(all$rc) <- list(modelid, 
						    as.character(nclus))
      structure(all$bic, equal = equal, tree = tree, subset
		= smpl, noise = noise, Vinv = Vinv, rcond = 
		all$rc, class = "emclust")
    }
  }
}

"emclust1" <- function(data, nclus, modelid = c("VVV", "VVV"), k, 
		       equal = F, noise, Vinv) 
{
  data <- as.matrix(data)
  if(length(modelid) == 1)
    modelid <- c(modelid, modelid)
  names(modelid) <- c("HC", "EM")
  if(missing(k)) {
    if(missing(noise)) {
      nclus <- if(missing(nclus)) 1:9 else sort(unique(nclus))
      ## no mhtree.EEV or mhtree.VEV
      tree <- switch(modelid[1],
		     EI = mhtree.EI(data),
		     VI = mhtree.VI(data),
		     EEE = mhtree.EEE(data),
		     VVV = mhtree.VVV(data),
		     stop("invalid model id for HC"))
      l <- length(nclus)
      BIC <- RC <- numeric(l)	
      ##-----------------------------------------------------------------------
      clss <- mhclass(tree, nclus)
      for(j in 1:l) {
	i <- nclus[j]
	if(i == 1) {
	  bicval <- bic(data, modelid = modelid[2])
	}
	else {
	  z <- me(data, modelid = modelid[2], ctoz(clss[
			  , as.character(i)]), equal = equal)
	  bicval <- bic(data, modelid = modelid[2], z, 
			equal = equal)
	}
	BIC[j] <- bicval
	if (!is.na(bicval)) RC[j] <- attr(bicval, "rcond")
      }
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      names(BIC) <- names(RC) <- as.character(nclus)
      structure(BIC, modelid = modelid, equal = equal, tree
		= tree, rcond = RC, class = "emclust1")
    }
    else {
      ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ## noise case
      ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      nclus <- if(missing(nclus)) 0:9 else sort(unique(nclus)
						)
      noise <- as.logical(noise)	## no mhtree.EEV or mhtree.VEV
      tree <- switch(modelid[1],
		     EI = mhtree.EI(data[!noise,  ]),
		     VI = mhtree.VI(data[!noise,  ]),
		     EEE = mhtree.EEE(data[!noise,  ]),
		     VVV = mhtree.VVV(data[!noise,  ]),
		     stop("invalid model id for HC"))
      l <- length(nclus)
      BIC <- RC <- numeric(l)	
      ##-----------------------------------------------------------------------
      clss <- mhclass(tree, nclus[nclus != 0])
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      cl <- numeric(nrow(data))
      k <- 1
      if(nclus[k] == 0) {
					# all noise --- same for all models
	one <- bic.EI(data, noise = T, Vinv = Vinv)
	BIC[k] <- one
	if (!is.na(bicval)) RC[k] <- attr(one, "rc")
	k <- k + 1
      }
      if(k < l) {
	for(j in k:l) {
	  i <- nclus[j]
	  cl[!noise] <- clss[, as.character(i)]
	  cl[noise] <- i + 1
	  z <- me(data, modelid = modelid[2], ctoz(cl), 
		  equal = equal, noise = T, Vinv = Vinv)
	  bicval <- bic(data, modelid = modelid[2], z, 
			equal = equal, noise = T, Vinv = Vinv)
	  BIC[j] <- bicval
	  if (!is.na(bicval)) RC[j] <- attr(bicval, "rc")
	}
      }
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      names(BIC) <- names(RC) <- as.character(nclus)
      structure(BIC, modelid = modelid, equal = equal, tree
		= tree, noise = noise, Vinv = Vinv, rcond = RC,
		class = "emclust1")
    }
  }
  else {
    ##=========================================================================
    ## hierarchical clustering with a sample
    ##=========================================================================
    if(missing(noise)) {
      nclus <- if(missing(nclus)) 1:9 else sort(unique(nclus))	
      ## no mhtree.EEV or mhtree.VEV
      smpl <- sample(1:nrow(data), k)
      tree <- switch(modelid[1],
		     EI = mhtree.EI(data[smpl,  ]),
		     VI = mhtree.VI(data[smpl,  ]),
		     EEE = mhtree.EEE(data[smpl,  ]),
		     VVV = mhtree.VVV(data[smpl,  ]),
		     stop("invalid model id for HC"))
      l <- length(nclus)
      BIC <- RC <- numeric(l)	
      ##-----------------------------------------------------------------------
      clss <- mhclass(tree, nclus)	
      for(j in 1:l) {
	i <- nclus[j]
	if(i == 1) {
	  bicval <- bic(data, modelid = modelid[2])
	}
	else {
	  pars <- mstep(data[smpl,  ], modelid = 
			modelid[2], ctoz(clss[, as.character(i)]), 
			equal = equal)
	  z <- do.call("estep", c(list(data = data, 
				       modelid = modelid[2]), pars))
	  z <- me(data, modelid = modelid[2], z, equal
		  = equal)
	  bicval <- bic(data, modelid = modelid[2], z, 
			equal = equal)
	}
	BIC[j] <- bicval
	if (!is.na(bicval)) RC[j] <- attr(bicval, "rcond")
      }
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      names(BIC) <- names(RC) <- as.character(nclus)
      structure(BIC, modelid = modelid, equal = equal, tree
		= tree, subset = smpl, rcond = RC, class = 
		"emclust1")
    }
    else {
      ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ## noise case
      ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      nclus <- if(missing(nclus)) 0:9 else sort(unique(nclus)
						)
      noise <- as.logical(noise)
      m <- sum(as.numeric(!noise))
      smpl <- sample(1:m, size = k)	## no mhtree.EEV or mhtree.VEV
      tree <- switch(modelid[1],
		     EI = mhtree.EI(data[!noise,  ][smpl,  ]),
		     VI = mhtree.VI(data[!noise,  ][smpl,  ]),
		     EEE = mhtree.EEE(data[!noise,  ][smpl,  ]),
		     VVV = mhtree.VVV(data[!noise,  ][smpl,  ]),
		     stop("invalid model id for HC"))
      l <- length(nclus)
      BIC <- RC <- numeric(l)	
      ##-----------------------------------------------------------------------
      clss <- mhclass(tree, nclus[nclus != 0])
      if(missing(Vinv))
	Vinv <- hypvol(data, reciprocal = T)
      cl <- numeric(nrow(data))
      k <- 1
      if(nclus[k] == 0) {
					# all noise --- same for all models
	one <- bic.EI(data, noise = T, Vinv = Vinv)
	BIC[k] <- one
	if (!is.na(bicval)) RC[k] <- attr(one, "rc")
	k <- k + 1
      }
      z <- matrix(0, nrow(data), nclus[l] + 1)
      if(k < l) {
	for(j in k:l) {
	  i <- nclus[j]
	  pars <- mstep(data[!noise,  ][smpl,  ], 
			modelid = modelid[2], ctoz(clss[, 
			  as.character(i)]), equal = equal)
	  z0 <- do.call("estep.EI", c(list(data[!noise,], 
					   modelid = modelid[2]), pars))
	  z[, 1:(i + 1)] <- 0
	  z[!noise, 1:i] <- z0
	  z[noise, i + 1] <- 1
	  z <- me(data, modelid = modelid[2], z[, 1:(i + 
			  1)], equal = equal, noise = T, Vinv = Vinv)
	  bicval <- bic(data, modelid = modelid[2], z[, 
				1:(i + 1)], equal = equal, noise = T, Vinv
			= Vinv)
	  BIC[j] <- bicval
	  if (!is.na(bicval)) RC[j] <- attr(bicval, "rc")
	}
      }
      ##-----------------------------------------------------------------------
      ## output
      ##-----------------------------------------------------------------------
      names(BIC) <- names(RC) <- as.character(nclus)
      structure(BIC, modelid = modelid, equal = equal, tree
		= tree, subset = smpl, noise = noise, Vinv = 
		Vinv, rcond = RC, class = "emclust1")
    }
  }
}

"estep" <- function(data, modelid, mu, ...)
{
					# ... sigsq or sigma, eps, Vinv
  switch(as.character(modelid),
	 EI = estep.EI(data, mu, ...),
	 VI = estep.VI(data, mu, ...),
	 EEE = estep.EEE(data, mu, ...),
	 VVV = estep.VVV(data, mu, ...),
	 EEV = estep.XEV(data, mu, ...),
	 VEV = estep.XEV(data, mu, ...),
	 stop("invalid model id"))
}

"estep.EEE" <- function(data, mu, sigma, prob, eps, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  G <- ncol(mu)
  if(missing(prob))
    prob <- rep(1/G, G)
  prob <- prob/sum(prob)
  equal <- length(unique(prob)) == 1
  l <- length(prob)
  noise <- l != G
  if(all(is.na(c(sigma, mu, prob)))) {
    z <- matrix(NA, n, if(noise) G + 1 else G)
    attr(z, "loglik") <- NA
    return(z)
  }
  else if(any(is.na(sigma)) || any(is.na(mu)) || any(is.na(prob))) {
    stop("parameters contain missing values")
  }
  if(!noise) {
					# no noise assumed
    K <- G
    temp <- .Fortran("eseee",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(p),
		     double(n * G),
		     as.double(eps))[9:10]
  }
  else {
    K <- G + 1
    if(l != K)
      stop("length(prob) = G+1 for noise")
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    temp <- .Fortran("esneee",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(p),
		     double(n * K),
		     as.double(eps),
		     as.double(Vinv))[9:10]
  }
  z <- matrix(NA, n, K)
  attr(z, "loglik") <- NA
  loglik <- temp[[2]]
  if(loglik == .Machine$double.xmin) {
    warning("sigma is not positive definite")
    return(z)
  }
  if(loglik ==  - .Machine$double.xmin) {
    warning("input error for LAPACK DPOTRF")
    return(z)
  }
  if(loglik == .Machine$double.xmax) {
    warning("sigma is nearly singular")
    return(z)
  }
  z[1:n, 1:K] <- temp[[1]]
  attr(z, "loglik") <- loglik
  z
}
"estep.EI" <-
function(data, mu, sigma, prob, eps, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  G <- ncol(mu)
  if(missing(prob))
    prob <- rep(1/G, G)
  prob <- prob/sum(prob)
  equal <- length(unique(prob)) == 1
  l <- length(prob)
  noise <- l != G
  if(all(is.na(c(sigma, mu, prob)))) {
    z <- matrix(NA, n, if(noise) G + 1 else G)
    attr(z, "loglik") <- NA
    return(z)
  }
  else if(is.na(sigma) || any(is.na(mu)) || any(is.na(prob))) {
    stop("parameters contain missing values")
  }
  if(sigma <= eps) {
    warning("sigma-squared falls below threshold")
    z <- matrix(NA, n, if(noise) G + 1 else G)
    attr(z, "loglik") <- NA
    return(z)
  }
  if(!noise) {
					# no noise assumed
    K <- G
    temp <- .Fortran("esei",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(n * G),
		     as.double(eps))[8:9]
  }
  else {
    K <- G + 1
    if(l != K)
      stop("length(prob) = G+1 for noise")
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    temp <- .Fortran("esnei",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(n * K),
		     as.double(eps),
		     as.double(Vinv))[8:9]
  }
  z <- matrix(NA, n, K)
  attr(z, "loglik") <- NA
  loglik <- temp[[2]]
  if(loglik == .Machine$double.xmax) {
    warning("sigma-squared falls below threshold")
    return(z)
  }
  z[1:n, 1:K] <- temp[[1]]
  attr(z, "loglik") <- loglik
  z
}
"estep.VI" <-
function(data, mu, sigma, prob, eps, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  G <- ncol(mu)
  if(missing(prob))
    prob <- rep(1/G, G)
  prob <- prob/sum(prob)
  equal <- length(unique(prob)) == 1
  l <- length(prob)
  noise <- l != G
  if(all(is.na(c(sigma, mu, prob)))) {
    z <- matrix(NA, n, if(noise) G + 1 else G)
    attr(z, "loglik") <- NA
    return(z)
  }
  else if(any(is.na(sigma)) || any(is.na(mu)) || any(is.na(prob))) {
    stop("parameters contain missing values")
  }
  if(any(sigma) <= eps) {
    warning("sigma-squared falls below threshold")
    z <- matrix(NA, n, length(prob))
    attr(z, "loglik") <- NA
    return(z)
  }
  if(!noise) {
					# no noise assumed
    K <- G
    temp <- .Fortran("esvi",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(n * G),
		     as.double(eps))[8:9]
  }
  else {
    K <- G + 1
    if(l != K)
      stop("length(prob) = G+1 for noise")
    if(missing(Vinv))
      Vinv <- hypvol(x, reciprocal = T)
    temp <- .Fortran("esnvi",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(n * K),
		     as.double(eps),
		     as.double(Vinv))[8:9]
  }
  z <- matrix(NA, n, K)
  attr(z, "loglik") <- NA
  loglik <- temp[[2]]
  if(loglik == .Machine$double.xmax) {
    warning("sigma-squared falls below threshold")
    return(z)
  }
  z[1:n, 1:K] <- temp[[1]]
  attr(z, "loglik") <- loglik
  z
}
"estep.VVV" <-
function(data, mu, sigma, prob, eps, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- .Machine$double.eps
  G <- ncol(mu)
  if(missing(prob))
    prob <- rep(1/G, G)
  prob <- prob/sum(prob)
  equal <- length(unique(prob)) == 1
  l <- length(prob)
  noise <- l != G
  if(all(is.na(c(sigma, mu, prob)))) {
    z <- matrix(NA, n, if(noise) G + 1 else G)
    attr(z, "loglik") <- NA
    return(z)
  }
  else if(any(is.na(sigma)) || any(is.na(mu)) || any(is.na(prob))) {
    stop("parameters contain missing values")
  }
  if(!noise) {
					# no noise assumed
    K <- G
    temp <- .Fortran("esvvv",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(p),
		     double(n * G),
		     as.double(eps))[9:10]
  }
  else {
    K <- G + 1
    if(l != K)
      stop("length(prob) = G+1 for noise")
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    temp <- .Fortran("esnvvv",
		     as.double(data),
		     as.double(mu),
		     as.double(sigma),
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     double(p),
		     double(n * K),
		     as.double(eps),
		     as.double(Vinv))[9:10]
  }
  z <- matrix(NA, n, K)
  attr(z, "loglik") <- NA
  loglik <- temp[[2]]
  if(loglik == .Machine$double.xmin) {
    warning("sigma is not positive definite")
    return(z)
  }
  if(loglik ==  - .Machine$double.xmin) {
    warning("input error for LAPACK DPOTRF")
    return(z)
  }
  if(loglik == .Machine$double.xmax) {
    warning("sigma is nearly singular")
    return(z)
  }
  z[1:n, 1:K] <- temp[[1]]
  attr(z, "loglik") <- loglik
  z
}
"estep.XEV" <-
function(data, mu, sigma, prob, eps, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(missing(eps))
    eps <- c(.Machine$double.eps, .Machine$double.eps)
  else if(length(eps) == 1)
    eps <- c(eps, .Machine$double.eps)
  G <- ncol(mu)
  if(missing(prob))
    prob <- rep(1/G, G)
  prob <- prob/sum(prob)
  equal <- length(unique(prob)) == 1
  l <- length(prob)
  noise <- l != G
  if(all(is.na(c(sigma, mu, prob)))) {
    z <- matrix(NA, n, if(noise) G + 1 else G)
    attr(z, "loglik") <- NA
    return(z)
  }
  else if(any(is.na(sigma)) || any(is.na(mu)) || any(is.na(prob))) {
    stop("parameters contain missing values")
  }
  lwork <- max(4 * p, 5 * p - 4)
  storage.mode(mu) <- "double"
  storage.mode(sigma) <- "double"
  if(!noise) {
					# no noise assumed
    K <- G
    temp <- .Fortran("esxev",
		     as.double(data),
		     mu,
		     sigma,
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(lwork),
		     as.integer(lwork),
		     double(n * G),
		     double(1))[c(2:3, 8, 12:13)]
  }
  else {
    K <- G + 1
    if(l != K)
      stop("length(prob) = G+1 for noise")
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    temp <- .Fortran("esnxev",
		     as.double(data),
		     mu,
		     sigma,
		     as.double(if(equal) 1 else prob),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(lwork),
		     as.integer(lwork),
		     double(n * K),
		     double(1),
		     as.double(Vinv))[c(2:3, 8, 12:13)]
  }
  rcmin <- temp[[3]][2]
  z <- matrix(NA, n, K)
  attr(z, "loglik") <- NA
  attr(z, "var.shape") <- max(apply(temp[[1]], 1, var))
  attr(z, "lambda") <- temp[[2]][1, 1,  ]
  attr(z, "rcmin") <- rcmin
  loglik <- temp[[5]]
  if(rcmin <= eps[2])
    warning("reciprocal condition number falls below threshold")
  if(loglik == .Machine$double.neg.eps) {
    stop("sigma is not positive definite")
  }
  if(loglik ==  - .Machine$double.neg.eps) {
    stop("input error in LAPACK DPOTRF")
  }
  if(loglik == .Machine$double.xmin) {
    stop("LAPACK DGESVD fails to converge")
  }
  if(loglik == .Machine$double.xmin) {
    stop("input error for LAPACK DGESVD")
  }
  if(loglik == .Machine$double.xmax) {
    warning("volume falls below threshold")
    return(z)
  }
  z[1:n, 1:K] <- temp[[4]]
  attr(z, "loglik") <- loglik
  z
}

"hypvol" <- function(data, reciprocal = F)
{
  ## finds the minimum hypervolume between principal components and 
  ## variable bounds
  data <- as.matrix(data)
  dimd <- dim(data)
  n <- dimd[1]
  p <- dimd[2]
#  if(F) {
#    vol1 <- prod(apply(data, 2, function(z)
#		       diff(range(z))))
#    V <- matrix(temp[[1]], p, p)
#    xbar <- apply(data, 2, mean)
#    X <- sweep(data, 2, xbar)
#    library(Matrix)
#    print(V)
#    print(eigen.Hermitian(crossprod(X))$vectors)
#    X <- X %*% V
#    vol <- prod(apply(X, 2, function(z)
#		      diff(range(z))))
#  }
  lwgesvd <- max(3 * min(n, p) + max(n, p), 5 * min(n, p) - 4)	# min
  lwsyevd <- p * (3 * p + 2 * ceiling(log(p, base = 2)) + 5) + 1	# min
  lisyevd <- 5 * p + 2	# minimum
  lwsyevx <- 8 * p	# minimum
  lisyevx <- 5 * p + p
  lwork <- max(lwsyevd, lwsyevx, n)
  liwork <- lisyevx
  temp <- .Fortran("mclvol",
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   double(p),
		   double(p * p),
		   double(p * p),
		   double(lwork),
		   as.integer(lwork),
		   integer(liwork),
		   as.integer(liwork),
		   integer(1))[c(4, 11)]
  if(temp[[2]])
    stop("problem in computing principal components")
  if(reciprocal) {
    pcvol <- prod(1/temp[[1]])
    bdvol <- prod(1/(apply(data, 2, max) - apply(data, 2, min)))
    ans <- max(pcvol, bdvol)
  }
  else {
    pcvol <- temp[[1]]
    bdvol <- prod(apply(data, 2, max) - apply(data, 2, min))
    ans <- min(pcvol, bdvol)
  }
  ans
}

"loglik" <- function(tree, data, ...)
{
  switch(attr(tree, "model"),
	 EI = loglik.EI(tree, data, ...),
	 VI = loglik.VI(tree, data, ...),
	 EEE = loglik.EEE(tree, data, ...),
	 VVV = loglik.VVV(tree, data, ...),
	 EEV = loglik.EEV(tree, data, ...),
	 VEV = loglik.VEV(tree, data, ...),
	 stop("invalid model id"))
}

"loglik.EEE" <- function(tree, data)
{
  n <- nrow(data)
  if(length(attr(tree, "initial.partition")) != n)
    stop("initial partition incompatible with number of observations"
	 )
  data <- as.matrix(data)
  temp <- .Fortran("likeee",
		   as.integer(tree),
		   as.integer(ns <- ncol(tree)),
		   as.double(data),
		   as.integer(n),
		   as.integer(p <- ncol(data)),
		   as.integer(attr(tree, "initial.partition")),
		   as.integer(length(unique(attr(tree, "initial.partition")))),
		   double(p),
		   double(p * p),
		   integer(ns),
		   double(ns + 1))[c(11, 10)]
  nmerge <- temp[[2]]
  temp <- temp[[1]]
  temp[temp ==  - .Machine$double.xmax] <- NA
  structure(temp, nmerge = nmerge)
}

"loglik.EFV" <- function(tree, data, Vinv)
{
  n <- nrow(data)
  if(length(attr(tree, "initial.partition")) != n)
    stop("initial partition incompatible with number of observations"
	 )
  p <- ncol(data)
  data <- as.matrix(data)
  if(missing(Vinv))
    Vinv <- hypvol(data, reciprocal = T)
  temp <- .Fortran("likefv",
		   as.integer(tree),
		   as.integer(ns <- ncol(tree)),
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   as.integer(attr(tree, "initial.partition")),
		   as.integer(length(unique(attr(tree, "initial.partition")))),
		   as.double(sqrt(attr(tree, "shape"))),
		   as.integer(lwork <- max(4 * p, 5 * p - 4)),
		   double(lwork),
		   double(p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   integer(ns),
		   double(ns + 1),
		   integer(1))[c(9, 16:18)]
  if(temp[[4]])
    stop("SVD does not converge")
  lwopt <- temp[[1]]
  temp <- temp[ - c(1, 4)]
  nmerge <- temp[[1]]
  temp <- temp[[2]]
  temp[temp ==  - .Machine$double.xmax] <- NA
  structure(temp, nmerge = nmerge)
}

"loglik.EI" <- function(tree, data)
{
  n <- nrow(data)
  if(length(attr(tree, "initial.partition")) != n)
    stop("initial partition incompatible with number of observations"
	 )
  data <- as.matrix(data)
  temp <- .Fortran("likei",
		   as.integer(tree),
		   as.integer(ns <- ncol(tree)),
		   as.double(data),
		   as.integer(n),
		   as.integer(p <- ncol(data)),
		   as.integer(attr(tree, "initial.partition")),
		   as.integer(length(unique(attr(tree, "initial.partition")))),
		   double(p),
		   integer(ns),
		   double(ns + 1))[c(10, 9)]
  nmerge <- temp[[2]]
  temp <- temp[[1]]
  temp[temp ==  - .Machine$double.xmax] <- NA
  structure(temp, nmerge = nmerge)
}

"loglik.VFV" <- function(tree, data, Vinv)
{
  n <- nrow(data)
  if(length(attr(tree, "initial.partition")) != n)
    stop("initial partition incompatible with number of observations"
	 )
  data <- as.matrix(data)
  p <- ncol(data)
  if(missing(Vinv))
    Vinv <- hypvol(data, reciprocal = T)
  temp <- .Fortran("likvfv",
		   as.integer(tree),
		   as.integer(ns <- ncol(tree)),
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   as.integer(attr(tree, "initial.partition")),
		   as.integer(length(unique(attr(tree, "initial.partition")))),
		   as.double(Vinv),
		   as.double(sqrt(attr(tree, "shape"))),
		   as.integer(lwork <- max(4 * p, 5 * p - 4)),
		   double(lwork),
		   double(p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   integer(ns),
		   double(ns + 1),
		   integer(1))[c(10, 17:19)]
  if(temp[[4]])
    stop("SVD does not converge")
  lwopt <- temp[[1]]
  temp <- temp[ - c(1, 4)]
  nmerge <- temp[[1]]
  temp <- temp[[2]]
  temp[temp ==  - .Machine$double.xmax] <- NA
  structure(temp, nmerge = nmerge)
}

"loglik.VI" <- function(tree, data, Vinv)
{
  n <- nrow(data)
  if(length(attr(tree, "initial.partition")) != n)
    stop("initial partition incompatible with number of observations"
	 )
  data <- as.matrix(data)
  p <- ncol(data)
  if(missing(Vinv))
    Vinv <- hypvol(data, reciprocal = T)
  temp <- .Fortran("likvi",
		   as.integer(tree),
		   as.integer(ns <- ncol(tree)),
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   as.integer(attr(tree, "initial.partition")),
		   as.integer(length(unique(attr(tree, "initial.partition")))),
		   as.double(Vinv),
		   double(p),
		   integer(ns),
		   double(ns + 1))[c(11, 10)]
  nmerge <- temp[[2]]
  temp <- temp[[1]]
  temp[temp ==  - .Machine$double.xmax] <- NA
  structure(temp, nmerge = nmerge)
}

"loglik.VVV" <- function(tree, data, Vinv)
{
  n <- nrow(data)
  if(length(attr(tree, "initial.partition")) != n)
    stop("initial partition incompatible with number of observations"
	 )
  data <- as.matrix(data)
  p <- ncol(data)
  if(missing(Vinv))
    Vinv <- hypvol(data, reciprocal = T)
  temp <- .Fortran("likvvv",
		   as.integer(tree),
		   as.integer(ns <- ncol(tree)),
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   as.integer(attr(tree, "initial.partition")),
		   as.integer(length(unique(attr(tree, "initial.partition")))),
		   as.double(Vinv),
		   double(p),
		   double(p * p),
		   integer(ns),
		   double(ns + 1))[c(12, 11)]
  nmerge <- temp[[2]]
  temp <- temp[[1]]
  temp[temp ==  - .Machine$double.xmax] <- NA
  structure(temp, nmerge = nmerge)
}

"me" <- function(data, modelid, z, ...)
{
  ## ... z, eps, tol, itmax, equal = F, noise = F, Vinv
  switch(as.character(modelid),
	 EI = me.EI(data, z, ...),
	 VI = me.VI(data, z, ...),
	 EEE = me.EEE(data, z, ...),
	 VVV = me.VVV(data, z, ...),
	 EEV = me.EEV(data, z, ...),
	 VEV = me.VEV(data, z, ...),
	 stop("invalid model id"))
}

"me.EEE" <- function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  if(all(is.na(z))) {
    attr(z, "info") <- c(iterations = NA, maxerr = NA, rcond = NA)
    return(z)
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  K <- dimz[2]	# number of groups
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  if(missing(itmax) || is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(!noise) {
    G <- K
    temp <- .Fortran("meeee",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(p * G),
		     double(p * p),
		     double(if(equal) 1 else G),
		     double(p))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, G)
    rc <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("meneee",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(p * G),
		     double(p * p),
		     double(if(equal) 1 else G),
		     double(p),
		     as.double(Vinv))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, K)
    rc <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  if(its >= itmax) {
    warning("iteration limit reached")
    attr(z, "warn") <- list("iteration limit reached")
    its <-  - its
  }
  if(rc <= abs(eps)) {
    warning("reciprocal condition estimate falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"),
	list("reciprocal condition estimate falls below threshold"))
  }
  attr(z, "info") <- c(iterations = its, maxerr = err, rcond = rc)
  z
}

"me.EEV" <- function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  if(all(is.na(z))) {
    attr(z, "info") <- c(iterations = NA, maxerr = NA, rcond = NA)
    return(z)
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  K <- dimz[2]	# number of groups
  if(missing(eps))
    eps <- c(.Machine$double.eps, sqrt(.Machine$double.eps))
  else if(length(eps) == 1)
    eps <- c(eps, sqrt(.Machine$double.eps))
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  if(missing(itmax) || is.infinite(itmax))
    itmax <- .Machine$integer.max
  lwork <- max(4 * p, 5 * p - 4)
  if(!noise) {
    G <- K
    temp <- .Fortran("meeev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(p * G),
		     double(p * p * G),
		     double(if(equal) 1 else G),
		     double(p),
		     double(p),
		     double(lwork),
		     as.integer(lwork))[c(2, 6:8)]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("meneev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(p * G),
		     double(p * p * G),
		     double(if(equal) 1 else G),
		     double(p),
		     double(p),
		     double(lwork),
		     as.integer(lwork),
		     as.double(Vinv))[c(2, 6:8)]
  }
  z <- matrix(NA, n, K)
  temp[2:3] <- lapply(temp[2:3], function(z)
		      {
			z[z == .Machine$double.xmax] <- NA
			z
		      }
		      )
  lamin <- temp[[2]][1]
  rcmin <- temp[[2]][2]
  err <- temp[[3]]
  its <- temp[[4]]
  if(its >= itmax) {
    warning("iteration limit reached")
    attr(z, "warn") <- list("iteration limit reached")
    its <-  - its
  }
  if(lamin == .Machine$double.xmax) {
    warning("LAPACK DGESVD fails to converge")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("LAPACK DGESVD fails to converge"))
  }
  else if(lamin ==  - .Machine$double.xmax) {
    warning("input error for LAPACK DGESVD")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("input error for LAPACK DGESVD"))
  }
  else if(lamin <= eps[1]) {
    warning("volume falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("volume falls below threshold"))
  }
  if(rcmin <= eps[2]) {
    warning("reciprocal condition estimate falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), 
	list("reciprocal condition estimate falls below threshold"))
  }
  else z[1:n, 1:K] <- temp[[1]]
  attr(z, "info") <- c(iterations = its, maxerr = err, rcond = rcmin)
  z
}
"me.EI" <-
function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  if(all(is.na(z))) {
    attr(z, "info") <- c(iterations = NA, maxerr = NA, rcond = NA)
    return(z)
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  K <- dimz[2]	# number of groups
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  if(missing(itmax) || is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(!noise) {
    G <- K
    temp <- .Fortran("meei",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(if(eps > 0) n * G else 1),
		     double(p),
		     double(if(equal) 1 else G))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, G)
    sigsq <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("menei",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(if(eps > 0) n * G else 1),
		     double(p),
		     double(if(equal) 1 else K),
		     as.double(Vinv))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, K)
    sigsq <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  rcond <- sigsq/(1 + sigsq)
  if(its >= itmax) {
    warning("iteration limit reached")
    attr(z, "warn") <- list("iteration limit reached")
    its <-  - its
  }
  if(sigsq <= abs(eps)) {
    warning("sigma-squared falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("sigma-squared falls below threshold"))
  }
  attr(z, "info") <- c(iterations = its, maxerr = err, rcond = rcond)
  z
}
"me.VEV" <- function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  if(all(is.na(z))) {
    attr(z, "info") <- c(iterations = NA, maxerr = NA, rcond = NA)
    return(z)
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  K <- dimz[2]	# number of groups
  if(missing(eps))
    eps <- c(.Machine$double.eps, .Machine$double.eps)
  else if(length(eps) == 1)
    eps <- c(eps, .Machine$double.eps)
  if(missing(tol))
    tol <- rep(sqrt(.Machine$double.eps), 2)
  else if(length(tol) == 1)
    tol <- c(tol, sqrt(.Machine$double.eps))
  if(missing(itmax))
    itmax <- c(Inf, Inf)
  else if(length(itmax) == 1)
    itmax <- c(itmax, Inf)
  itmax[itmax == Inf] <- .Machine$integer.max
  lwork <- max(4 * p, 5 * p - 4)
  if(!noise) {
    G <- K
    lwork <- max(lwork, p + G)
    temp <- .Fortran("mevev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(p * G),
		     double(p * p * G),
		     double(G),
		     double(G),
		     double(p),
		     double(p),
		     double(p * G),
		     double(lwork),
		     as.integer(lwork))[c(2, 6:8)]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    lwork <- max(lwork, p + G)
    temp <- .Fortran("menvev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(p * G),
		     double(p * p * G),
		     double(G),
		     double(G),
		     double(p),
		     double(p),
		     double(p * G),
		     double(lwork),
		     as.integer(lwork),
		     as.double(Vinv))[c(2, 6:8)]
  }
  z <- matrix(NA, n, K)
  temp[2:3] <- lapply(temp[2:3], function(z)
		      {
			z[z == .Machine$double.xmax] <- NA
			z
		      }
		      )
  lamin <- temp[[2]][1]
  rcmin <- temp[[2]][2]
  its <- temp[[4]][1]
  inner <- temp[[4]][2]
  err <- if(its > 1) temp[[3]][1] else NA
  inerr <- temp[[3]][2]
  if(its >= itmax[1]) {
    warning("iteration limit reached")
    attr(z, "warn") <- list("iteration limit reached")
    its <-  - its
  }
  if(itmax[2] > 0 && inner >= itmax[2]) {
    warning("inner iteration limit reached")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("inner iteration limit reached"))
    inner <-  - inner
  }
  if(lamin == .Machine$double.xmax) {
    warning("LAPACK DGESVD fails to converge")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("LAPACK DGESVD fails to converge"))
  }
  else if(lamin ==  - .Machine$double.xmax) {
    warning("input error for LAPACK DGESVD")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("input error for LAPACK DGESVD"))
  }
  else if(lamin <= eps[1]) {
    warning("volume falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("volume falls below threshold"))
  }
  if(rcmin <= eps[2]) {
    warning("reciprocal condition estimate falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), 
	list("reciprocal condition estimate falls below threshold"))
  }
  else z[1:n, 1:K] <- temp[[1]]
  attr(z, "info") <- c(iterations = its, maxerr = err, rcond = rcmin)
#  attr(attr(z, "info"), "inner") <- c(iterations = inner, maxerr = inerr)
  z
}

"me.VI" <- function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  if(all(is.na(z))) {
    attr(z, "info") <- c(iterations = NA, maxerr = NA, rcond = NA)
    return(z)
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  K <- dimz[2]	# number of groups
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  if(missing(itmax) || is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(!noise) {
    G <- K
    temp <- .Fortran("mevi",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(if(eps > 0) n * G else 1),
		     double(p))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, G)
    sigmin <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("menvi",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(if(eps > 0) n * G else 1),
		     double(p),
		     as.double(Vinv))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, K)
    sigmin <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  rcond <- sigmin/(1 + sigmin)
  if(its >= itmax) {
    warning("iteration limit reached")
    attr(z, "warn") <- list("iteration limit reached")
    its <-  - its
  }
  if(sigmin <= abs(eps)) {
    warning("sigma-squared falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), list("sigma-squared falls below threshold"))
  }
  attr(z, "info") <- c(iterations = its, maxerr = err, rcond = rcond)
  z
}
"me.VVV" <-
function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  if(all(is.na(z))) {
    attr(z, "info") <- c(iterations = NA, maxerr = NA, rcond = NA)
    return(z)
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  K <- dimz[2]	# number of groups
  if(missing(eps))
    eps <- .Machine$double.eps
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  if(missing(itmax) || is.infinite(itmax))
    itmax <- .Machine$integer.max
  if(!noise) {
    G <- K
    temp <- .Fortran("mevvv",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(if(eps > 0) n * G else 1),
		     double(p),
		     double(p * p),
		     double(p))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, G)
    rcmin <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("menvvv",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(if(eps > 0) n * G else 1),
		     double(p),
		     double(p * p),
		     double(p),
		     as.double(Vinv))[c(2, 6:8)]
    z <- matrix(temp[[1]], n, K)
    rcmin <- temp[[2]]
    err <- temp[[3]]
    its <- temp[[4]]
  }
  if(its >= itmax) {
    warning("iteration limit reached")
    attr(z, "warn") <- list("iteration limit reached")
    its <-  - its
  }
  if(rcmin <= abs(eps)) {
    warning("reciprocal condition estimate falls below threshold")
    attr(z, "warn") <- 
      c(attr(z, "warn"), 
	list("reciprocal condition estimate falls below threshold"))
  }
  attr(z, "info") <- c(iterations = its, maxerr = err, rcond = rcmin)
  z
}

"mhclass" <- function(tree, nclusters)
{
  initial <- attributes(tree)$init
  n <- length(initial)
  k <- length(unique(initial))
  nclusters <- if(missing(nclusters)) 
    k:2 
  else 
    rev(sort(unique(nclusters)))
  select <- k - nclusters
  if(length(select) == 1 && !select)
    return(matrix(initial, ncol = 1, 
		  dimnames = list(NULL, as.character(nclusters))))
  bad <- select < 0 | select >= k
  if(all(bad))
    stop("No classification with the specified number of clusters")
  if(any(bad))
    warning("Some selected classifications are inconsistent\n                          with mclust object"
	    )	# all stages
  l <- length(select)
  cl <- matrix(NA, nrow = n, ncol = l, 
	       dimnames = list(NULL, as.character(nclusters)))
  if(select[1])
    m <- 1
  else {
    cl[, 1] <- initial
    m <- 2
  }
  for(l in 1:max(select)) {
    ij <- tree[, l]
    i <- min(ij)
    j <- max(ij)
    initial[initial == j] <- i
    if(select[m] == l) {
      cl[, m] <- initial
      m <- m + 1
    }
  }
  apply(cl[, length(select):1, drop = F], 2, partconv, consec = T)
}

"mhtree" <- function(data, modelid, partition, min.clusters = 1, 
		    verbose = F, ...) 
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) 
    warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  if(missing(modelid))
    modelid <- "VVV"
  modelblurb <- switch(modelid,
		       EI = "EI : uniform spherical",
		       VI = "VI : spherical",
		       EEE = "EEE : uniform variance",
		       VVV = "VVV : unconstrained variance",
		       EFV = "EFV : fixed shape, uniform volume",
		       VFV = "VFV : fixed shape",
		       "invalid model id")
  if(verbose) cat("model ", modelblurb, "\n")	
  ##===========================================================================
  switch(modelid,
	 EI = mhtree.EI(data, partition, min.clusters, ...),
	 VI = mhtree.VI(data, partition, min.clusters, ...),
	 EEE = mhtree.EEE(data, partition, min.clusters, ...),
	 VVV = mhtree.VVV(data, partition, min.clusters, ...),
	 EFV = mhtree.EFV(data, partition, min.clusters, ...),
	 VFV = mhtree.VFV(data, partition, min.clusters, ...),
	 stop("invalid model id"))
}

"mhtree.EEE" <- function(data, partition, min.clusters = 1)
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = T)
  l <- length(unique(partition))
  attr(partition, "unique") <- l	
  ##===========================================================================
  m <- l - min.clusters
  if(m <= 0) {
    return(structure(NA, determinant = NA, trace = NA, 
		     initial.partition = partition, dimension = dimx, 
		     call = match.call(), class = "mhtree"))
  }
  storage.mode(data) <- "double"
  temp <- .Fortran("hceee",
		   data,
		   as.integer(n),
		   as.integer(p),
		   as.integer(partition),
		   as.integer(l),
		   as.integer(m),
		   if(p < 3) integer(m) else integer(1),
		   if(p < 4) integer(m) else integer(1),
		   double(p),
		   double(p * p),
		   double(p * p),
		   double(p * p))[c(1, 7:10)]
  temp[[4]] <- temp[[4]][1:2]	## currently temp[[5]] is not output
  temp[[5]] <- temp[[5]][1:2]
  names(temp[[5]]) <- c("determinant", "trace")
  temp[[1]] <- temp[[1]][1:(m + 1),  ]
  if(p < 3)
    tree <- rbind(temp[[2]], temp[[3]])
  else if(p < 4)
    tree <- rbind(temp[[1]][-1, 3], temp[[3]])
  else tree <- t(temp[[1]][-1, 3:4, drop = F])
  determinant <- temp[[1]][, 1]
  attr(determinant, "breakpoints") <- temp[[4]]
  structure(tree, determinant = determinant, trace = temp[[1]][, 2], 
	    initial.partition = partition, dimensions = dimx, 
	    ##call = match.call(),
	    model = "EEE", class = "mhtree")
}

"mhtree.EFV" <- function(data, partition, min.clusters = 1, shape)
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(shape))
    stop("no shape vector specified")
  shape <- sqrt(rev(sort(shape)))
  shape <- shape/shape[1]
  if(length(shape) != p)
    stop("length of shape vector must equal ncol(data)")
  if(any(shape <= 0)) stop("nonpositive shape")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = T)
  l <- length(unique(partition))
  attr(partition, "unique") <- l	
  ##===========================================================================
  m <- l - min.clusters
  if(m <= 0) {
    return(structure(NA, change = NA, initial.partition = partition,
		     shape = shape^2, dimension = dimx, call = match.call(), 
		     class = "mhtree"))
  }
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ld <- max(n, ll + 1, 3 * m)	
  ##	dp <- duplicated(partition)
  ##    x[c((1:n)[!dp],(1:n)[dp]), ], 
  ##    as.integer(c(partition[!dp], partition[dp])), 
  temp <- .Fortran("hcefv",
		   data,
		   as.integer(n),
		   as.integer(p),
		   as.integer(partition),
		   as.integer(l),
		   as.integer(m),
		   as.integer(lwork <- max(4 * p, 5 * p - 4)),
		   as.double(shape),
		   double(p),
		   double(lwork),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   as.integer(ld),
		   double(ld))[c(1, 10, 15, 16)]
  if(temp[[3]])
    stop("SVD does not converge")
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = F]
  temp[[4]] <- temp[[4]][1:m]
  structure(t(temp[[1]]), change = temp[[4]], initial.partition = 
	    partition, shape = shape^2, dimensions = dimx, workspace.svd = 
	    c(use = lwork, opt = temp[[2]]), call = match.call(), model = 
	    "EFV", class = "mhtree")
}

"mhtree.EI" <- function(data, partition, min.clusters = 1)
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = T)
  l <- length(unique(partition))
  attr(partition, "unique") <- l	
  ##===========================================================================
  m <- l - min.clusters
  if(m <= 0) {
    return(structure(NA, change = NA, initial.partition = partition,
		     dimensions = dimx, call = match.call(), class = 
		     "mhtree"))
  }
  storage.mode(data) <- "double"
  ld <- max(c((l * (l - 1))/2, 3 * m))
  temp <- .Fortran("hcei",
		   data,
		   as.integer(n),
		   as.integer(p),
		   as.integer(partition),
		   as.integer(l),
		   as.integer(m),
		   double(p),
		   as.integer(ld),
		   double(ld))[c(1, 9)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = F]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), change = temp[[2]], 
	    initial.partition = partition, dimensions = dimx,
	    ##call = match.call(), 
	    model = "EI", class = "mhtree")
}

"mhtree.VEV" <- "mhtree.VFV" <- function(data, partition, 
				       min.clusters = 1, shape, alpha = 1) 
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(shape))
    stop("no shape vector specified")
  shape <- sqrt(rev(sort(shape)))
  shape <- shape/shape[1]
  if(length(shape) != p)
    stop("length of shape vector must equal ncol(data)")
  if(any(shape <= 0)) stop("nonpositive shape")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = T)
  l <- length(unique(partition))
  attr(partition, "unique") <- l	
  ##===========================================================================
  m <- l - min.clusters
  if(m <= 0) {
    return(structure(NA, change = NA, initial.partition = partition,
		     shape = shape^2, dimension = dimx, call = match.call(), 
		     class = "mhtree"))
  }
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ld <- max(n, ll + 1, 3 * m)	
  ##	dp <- duplicated(partition)
  ##    x[c((1:n)[!dp],(1:n)[dp]), ], 
  ##    as.integer(c(partition[!dp], partition[dp])), 
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvfv",
		   data,
		   as.integer(n),
		   as.integer(p),
		   as.integer(partition),
		   as.integer(l),
		   as.integer(m),
		   as.integer(lwork <- max(4 * p, 5 * p - 4)),
		   as.double(alpha),
		   as.double(shape),
		   double(p),
		   double(lwork),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   as.integer(ld),
		   double(ld))[c(1, 7, 16, 17)]
  if(temp[[3]])
    stop("SVD does not converge")
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = F]
  temp[[4]] <- temp[[4]][1:m]
  structure(t(temp[[1]]), change = temp[[4]], initial.partition = 
	    partition, shape = shape^2, dimensions = dimx, workspace.svd = 
	    c(use = lwork, opt = temp[[2]]), call = match.call(), model = 
	    "VFV", class = "mhtree")
}
"mhtree.VI" <-
function(data, partition, min.clusters = 1, alpha = 1)
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = T)
  l <- length(unique(partition))
  attr(partition, "unique") <- l	
  ##===========================================================================
  m <- l - min.clusters
  if(m <= 0) {
    return(structure(NA, change = NA, initial.partition = partition,
		     dimension = dimx, call = match.call(), class = "mhtree"
		     ))
  }
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ld <- max(n, ll, 3 * m)
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvi",
		   data,
		   as.integer(n),
		   as.integer(p),
		   as.integer(partition),
		   as.integer(l),
		   as.integer(m),
		   as.double(alpha),
		   double(p),
		   as.integer(ld),
		   double(ld))[c(1, 10)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = F]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), change = temp[[2]], initial.partition = 
	    partition, dimensions = dimx, call = match.call(), model = "VI",
	    class = "mhtree")
}

"mhtree.VVV" <-function(data, partition, min.clusters = 1, alpha = 1, beta = 1)
{
  if(min.clusters < 1)
    stop("min.clusters must be positive")
  if(any(is.na(data)))
    stop("missing values not allowed in data")
  dimx <- dim(data)
  if(is.null(dimx))
    stop("data must be a matrix with  at least 2 columns")
  data <- as.matrix(data)
  n <- dimx[1]
  p <- dimx[2]
  if(p == 1)
    stop("data must be a matrix with  at least 2 columns")
  if(n <= p) warning("# of observations <= data dimension")	
  ##===========================================================================
  if(missing(partition))
    partition <- 1:n
  else if(length(partition) != n)
    stop("partition must assign a class to each observation")
  partition <- partconv(partition, consec = T)
  l <- length(unique(partition))
  attr(partition, "unique") <- l	
  ##===========================================================================
  m <- l - min.clusters
  if(m <= 0) {
    return(structure(NA, change = NA, initial.partition = partition,
		     dimension = dimx, call = match.call(), class = "mhtree"
		     ))
  }
  storage.mode(data) <- "double"
  ll <- (l * (l - 1))/2
  ld <- max(n, ll + 1, 3 * m)	
  ##	dp <- duplicated(partition)
  ##    x[c((1:n)[!dp],(1:n)[dp]), ], 
  ##    as.integer(c(partition[!dp], partition[dp])), 
  alpha <- alpha * traceW(data/sqrt(n * p))
  alpha <- max(alpha, .Machine$double.eps)
  temp <- .Fortran("hcvvv",
		   cbind(data, 0),
		   as.integer(n),
		   as.integer(p),
		   as.integer(partition),
		   as.integer(l),
		   as.integer(m),
		   as.double(alpha),
		   as.double(beta),
		   double(p),
		   double(p * p),
		   double(p * p),
		   double(p * p),
		   as.integer(ld),
		   double(ld))[c(1, 14)]
  temp[[1]] <- temp[[1]][1:m, 1:2, drop = F]
  temp[[2]] <- temp[[2]][1:m]
  structure(t(temp[[1]]), change = temp[[2]], initial.partition = 
	    partition, dimensions = dimx, call = match.call(), model = 
	    "VVV", class = "mhtree")
}

"mixproj" <- function(data, ms, partition, dimens, scale = F, 
                      k = 15, title, xlim, ylim, xlab, ylab,
                      col=partition, pch=partition, ...)
{
  if(missing(dimens))
    dimens <- sample(1:ncol(data), 2)
  data <- data[, dimens]
  if(dim(data)[2] != 2)
    stop("need two dimensions")
  p <- nrow(ms$mu)
  d <- dim(ms$sigma)
  ld <- length(d)
  m <- ncol(ms$mu)
  if(is.null(d)) {
    l <- length(ms$sigma)
    if(l == 1) {
      ## uniform spherical model (EI)
      sigma <- array(diag(c(ms$sigma, ms$sigma)), c(2, 2, m))
    }
    else if(l == m) {
      ## spherical model (VI)
      sigma <- array(sapply(ms$sigma, function(z) diag(c(z, z))), c(2, 2, m))
    }
    else stop("mu and sigma are incompatible")
  }
  else if(ld == 3 && all(d[1:2] == p) && d[3] == m) {
    ## nonuniform variances (VVV, etc)
    sigma <- ms$sigma[dimens, dimens,  ]
  }
  else if(ld == 2 && all(d == p)) {
    ## uniform variance (EEE)
    sigma <- array(ms$sigma[dimens, dimens], c(2, 2, m))
  }
  else stop("mu and sigma are incompatible")
  mu <- ms$mu[dimens,  ]

  if (missing(xlim))
    xlim <- range(data[, 1])
  if (missing(ylim))
    ylim <- range(data[, 2])
  
  if(scale) {
    d <- diff(xlim) - diff(ylim)
    if(d > 0) 
      ylim <- c(ylim[1] - d/2, ylim[2] + d/2)
    else 
      xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
  }

  if (missing(partition))
    partition _ rep(1, dim(data)[1])
  labs _ dimnames(data)[[2]]
  if (missing(xlab)) {
    if (is.null(labs))
      xlab _ ""
    else
      xlab _ labs[1]
  }
  if (missing(ylab)) {
    if (is.null(labs))
      ylab _ ""
    else
      ylab _ labs[2]
  }

  plot(data[, 1], data[, 2], col=col, pch=pch, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim, ...)
  if(!missing(title))
    title(title)
  l <- ncol(mu)
  for(i in 1:l) {
    mvn2plot(mu = mu[,i], sigma = sigma[,,i], k = k)
  }
  dimens
}

"mstep" <- function(data, modelid, z, ...)
{
  ## ... z, eps, tol, itmax, equal = F, noise = F, Vinv
  switch(as.character(modelid),
	 EI = mstep.EI(data, z, ...),
	 VI = mstep.VI(data, z, ...),
	 EEE = mstep.EEE(data, z, ...),
	 VVV = mstep.VVV(data, z, ...),
	 EEV = mstep.EEV(data, z, ...),
	 VEV = mstep.VEV(data, z, ...),
	 stop("invalid model id"))
}

"mstep.EEE" <-function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  np <- n * p
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]	# number of groups
  if(all(is.na(z))) {
    G <- if(noise) K - 1 else K
    return(
	   if (equal) list(mu = matrix(NA, p, G), 
			   sigma = matrix(NA, p, p)) 
	   else list(mu = matrix(NA, p, G), 
		     sigma = matrix(NA, p, p), prob = rep(NA, K)) 
	   )
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Machine$double.eps
  if(!noise) {
					# no noise assumed
    G <- K
    temp <- .Fortran("mseee",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(1),
		     double(p * G),
		     double(p * p),
		     double(if(equal) 1 else G))[c(6, 8:11)]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("msneee",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(1),
		     double(p * G),
		     double(p * p),
		     double(if(equal) 1 else K),
		     as.double(Vinv))[c(6, 8:11)]
  }
  rc <- temp[[1]]
  loglik <- temp[[2]]
  mu <- matrix(temp[[3]], p, G)
  Sigma <- matrix(temp[[4]], p, p)
  if(!equal)
    prob <- temp[[5]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  out <- if(equal) list(mu = mu, sigma = Sigma) else list(mu = mu, sigma
				   = Sigma, prob = prob)
  if(rc <= abs(eps)) {
    warning("reciprocal condition estimate falls below threshold")
    attr(out, "warn") <- 
      list("reciprocal condition estimate falls below threshold")
  }
  if(loglik == .Machine$double.xmax)
    loglik <- NA
  attr(out, "loglik") <- loglik
  attr(out, "rcond") <- rc
  out
}

"mstep.EEV" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  np <- n * p
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]	# number of groups
  if(all(is.na(z))) {
    G <- if(noise) K - 1 else K
    return(
	   if(equal) 
	   list(mu = matrix(NA, p, G),
		sigma = array(NA, c(p, p, G))) 
	   else 
	   list(mu = matrix(NA, p, G), 
		sigma = array(NA, c(p, p, G)), prob = rep(NA, K))
	   )
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1)) 
    stop("improper specification of z")	
  ##	shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(missing(eps))
    eps <- c(.Machine$double.eps, .Machine$double.eps)
  else if(length(eps) == 1)
    eps <- c(eps, .Machine$double.eps)
  if(!noise) {
					# no noise assumed
    G <- K
    lwork <- max(4 * p, 5 * p - 4, G)
    temp <- .Fortran("mseev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(p),
		     double(lwork),
		     as.integer(lwork),
		     double(1),
		     double(p * G),
		     double(p * p * G),
		     double(if(equal) 1 else G))[c(6, 7, 11:14)]
  } else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    lwork <- max(4 * p, 5 * p - 4, G)
    temp <- .Fortran("msneev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(p),
		     double(lwork),
		     as.integer(lwork),
		     double(1),
		     double(p * G),
		     double(p * p * G),
		     double(if(equal) 1 else K),
		     as.double(Vinv))[c(6, 7, 11:14)]
  }
  lambda <- temp[[1]][1]
  rc <- temp[[1]][2]
  shape <- temp[[2]]
  loglik <- temp[[3]]
  mu <- matrix(temp[[4]], p, G)
  sigma <- array(temp[[5]], c(p, p, G))
  if(!equal)
    prob <- temp[[6]]
  warn <- NULL
  if(loglik == .Machine$double.xmax) {
    loglik <- NA
    warning("volume falls below threshold")
    warn <- list("volume falls below threshold")
    sigma <- array(NA, c(p, p, G))
    shape <- rep(NA, p)
  }
  else if(loglik ==  - .Machine$double.xmax) {
    loglik <- NA
    warning("reciprocal condition estimate falls below threshold")
    warn <- list("reciprocal condition estimate falls below threshold"
		 )
    sigma <- array(NA, c(p, p, G))
  }
  else if(loglik == .Machine$double.xmin) {
    loglik <- NA
    warning("LAPACK DGESVD fails to converge")
    warn <- list("LAPACK DGESVD fails to converge")
    sigma <- array(NA, c(p, p, G))
    lambda <- NA
    shape <- rep(NA, p)
  }
  else if(loglik ==  - .Machine$double.xmin) {
    loglik <- NA
    warning("input error for LAPACK DGESVD")
    warn <- list("input error for LAPACK DGESVD")
    sigma <- array(NA, c(p, p, G))
    lambda <- NA
    shape <- rep(NA, p)
  }
  dimnames(mu) <- list(NULL, as.character(1:G))
  out <- if(equal) list(mu = mu, sigma = sigma) else list(mu = mu, sigma
				   = sigma, prob = prob)
  attr(out, "loglik") <- loglik
  attr(out, "rcond") <- rc
  attr(out, "shape") <- shape * shape
  attr(out, "lambda") <- lambda
  attr(out, "warn") <- warn
  out
}

"mstep.EI" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  np <- n * p
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]	# number of groups
  if(all(is.na(z))) {
    G <- if(noise) K - 1 else K
    return(if(equal) list(mu = matrix(NA, p, G), sigma = NA)
    else list(mu = matrix(NA, p, G), sigma = NA, prob = 
	      rep(NA, K)))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Machine$double.eps
  if(!noise) {
					# no noise assumed
    G <- K
    temp <- .Fortran("msei",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p * G),
		     double(1),
		     double(if(equal) 1 else G))[6:9]
  } else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("msnei",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p * G),
		     double(1),
		     double(if(equal) 1 else K),
		     as.double(Vinv))[6:9]
  }
  loglik <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  sigma <- temp[[3]]
  if(!equal)
    prob <- temp[[4]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  out <- if(equal) list(mu = mu, sigma = sigma) else list(mu = mu, 
				   sigma = sigma, prob = prob)
  if(sigma <= abs(eps)) {
    warning("sigma-squared falls below threshold")
    attr(out, "warn") <- list("sigma-squared falls below threshold"
			      )
  }
  if(loglik == .Machine$double.xmax)
    loglik <- NA
  attr(out, "loglik") <- loglik
  attr(out, "rcond") <- sigma/(1 + sigma)
  out
}

"mstep.VEV" <- function(data, z, eps, tol, itmax, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  np <- n * p
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]	# number of groups
  if(all(is.na(z))) {
    G <- if(noise) K - 1 else K
    return(
	   if (equal) 
	   list(mu = matrix(NA, p, G), 
		sigma = array(NA, c(p, p, G))) 
	   else 
	   list(mu = matrix(NA, p, G), 
		sigma = array(NA, c(p, p, G)), prob = rep(NA, K))
	   )
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1)) 
    stop("improper specification of z")	
  ##	shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(missing(eps))
    eps <- c(.Machine$double.eps, .Machine$double.eps)
  else if(length(eps) == 1)
    eps <- c(eps, .Machine$double.eps)
  if(missing(tol))
    tol <- sqrt(.Machine$double.eps)
  if(missing(itmax))
    itmax <- .Machine$integer.max
  if(!noise) {
					# no noise assumed
    G <- K
    lwork <- max(4 * p, 5 * p - 4, G)
    temp <- .Fortran("msvev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(G),
		     double(p),
		     double(p),
		     double(p * G),
		     double(lwork),
		     as.integer(lwork),
		     double(1),
		     double(p * G),
		     double(p * p * G),
		     double(G))[c(6:10, 15:18)]
  } else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    lwork <- max(4 * p, 5 * p - 4, G)
    temp <- .Fortran("msnvev",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     as.double(tol),
		     as.integer(itmax),
		     double(G),
		     double(p),
		     double(p),
		     double(p * G),
		     double(lwork),
		     as.integer(lwork),
		     double(1),
		     double(p * G),
		     double(p * p * G),
		     double(K),
		     as.double(Vinv))[c(6:10, 15:18)]
  }
  rc <- temp[[1]][2]
  inner <- temp[[3]]
  inerr <- if(inner > 0) temp[[2]] else NA
  lambda <- temp[[4]]
  shape <- temp[[5]]
  loglik <- temp[[6]]
  mu <- matrix(temp[[7]], p, G)
  sigma <- array(temp[[8]], c(p, p, G))
  if(!equal)
    prob <- temp[[9]]
  warn <- NULL
  if(inner >= itmax) {
    warning("inner iteration limit reached")
    warn <- list("inner iteration limit reached")
  }
  if(loglik == .Machine$double.xmax) {
    loglik <- NA
    warning("volume falls below threshold")
    warn <- c(warn, list("volume falls below threshold"))
    sigma <- array(NA, c(p, p, G))
    shape <- rep(NA, p)
  }
  else if(loglik ==  - .Machine$double.xmax) {
    loglik <- NA
    warning("condition estimate falls below threshold")
    warn <- c(warn, list("condition estimate falls below threshold"
			 ))
    sigma <- array(NA, c(p, p, G))
    lambda <- rep(NA, G)
    shape <- rep(NA, p)
  }
  else if(loglik == .Machine$double.xmin) {
    loglik <- NA
    warning("LAPACK DGESVD fails to converge")
    warn <- list("LAPACK DGESVD fails to converge")
    sigma <- array(NA, c(p, p, G))
    lambda <- rep(NA, G)
    shape <- rep(NA, p)
  }
  else if(loglik ==  - .Machine$double.xmin) {
    loglik <- NA
    warning("input error for LAPACK DGESVD")
    warn <- list("input error for LAPACK DGESVD")
    sigma <- array(NA, c(p, p, G))
    lambda <- rep(NA, G)
    shape <- rep(NA, p)
  }
  dimnames(mu) <- list(NULL, as.character(1:G))
  out <- 
    if (equal) 
      list(mu = mu, sigma = sigma) 
    else 
      list(mu = mu, sigma = sigma, prob = prob)
  attr(out, "loglik") <- loglik
  attr(out, "rcond") <- rc
  attr(out, "shape") <- shape * shape
  attr(out, "lambda") <- lambda
  attr(out, "warn") <- warn
  ## Next line gives 
  ## Error: attempt to set an attribute on NULL
  ##  attr(attr(out, "info"), "inner") <- 
  ##    c(iterations = inner, maxerr = inerr)
  out
}

"mstep.VI" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  np <- n * p
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]	# number of groups
  if(all(is.na(z))) {
    G <- if(noise) K - 1 else K
    return(
	   if(equal) 
	   list(mu = matrix(NA, p, G), 
		sigma = rep(NA, G)) 
	   else list(mu = matrix(NA, p, G), 
		     sigma = rep(NA, G), prob = rep(NA, K))
	   )
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Machine$double.eps
  if(!noise) {
					# no noise assumed
    G <- K
    temp <- .Fortran("msvi",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p * G),
		     double(G),
		     double(if(equal) 1 else G))[6:9]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("msnvi",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p * G),
		     double(G),
		     double(if(equal) 1 else K),
		     as.double(Vinv))[6:9]
  }
  loglik <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  sigma <- temp[[3]]
  if(!equal)
    prob <- temp[[4]]
  temp <- min(sigma)
  dimnames(mu) <- list(NULL, as.character(1:G))
  out <- if(equal) list(mu = mu, sigma = sigma) else list(mu = mu, 
				   sigma = sigma, prob = prob)
  if(temp <= abs(eps)) {
    warning("sigma-squared falls below threshold")
    attr(out, "warn") <- list("sigma-squared falls below threshold"
			      )
  }
  if(loglik == .Machine$double.xmax)
    loglik <- NA
  attr(out, "loglik") <- loglik
  attr(out, "rcond") <- temp/(1 + temp)
  out
}

"mstep.VVV" <- function(data, z, eps, equal = F, noise = F, Vinv)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  np <- n * p
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]	# number of groups
  if(all(is.na(z))) {
    G <- if(noise) K - 1 else K
    return(
	   if(equal) 
	   list(mu = matrix(NA, p, G), 
			  sigma = array(NA, c(p, p, G))) 
	   else list(mu = matrix(NA, p, G), 
		     sigma = array(NA, c(p, p, G)), prob = rep(NA, K))
	   )
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(missing(eps))
    eps <- .Machine$double.eps
  if(!noise) {
					# no noise assumed
    G <- K
    temp <- .Fortran("msvvv",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(1),
		     double(p * G),
		     double(p * p * G),
		     double(if(equal) 1 else G))[c(6, 8:11)]
  }
  else {
    if(missing(Vinv))
      Vinv <- hypvol(data, reciprocal = T)
    G <- K - 1
    temp <- .Fortran("msnvvv",
		     as.double(data),
		     as.double(z),
		     as.integer(n),
		     as.integer(p),
		     as.integer(if(equal)  - G else G),
		     as.double(eps),
		     double(p),
		     double(1),
		     double(p * G),
		     double(p * p * G),
		     double(if(equal) 1 else K),
		     as.double(Vinv))[c(6, 8:11)]
  }
  rc <- temp[[1]]
  loglik <- temp[[2]]
  mu <- matrix(temp[[3]], p, G)
  sigma <- array(temp[[4]], c(p, p, G))
  if(!equal)
    prob <- temp[[5]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  out <- if(equal) list(mu = mu, sigma = sigma) else list(mu = mu, sigma
				   = sigma, prob = prob)
  if(rc <= abs(eps)) {
    warning("reciprocal condition estimate falls below threshold")
    attr(out, "warn") <- 
      list("reciprocal condition estimate falls below threshold")
  }
  if(loglik == .Machine$double.xmax)
    loglik <- NA
  attr(out, "loglik") <- loglik
  attr(out, "rcond") <- rc
  out
}

"mvn2plot" <- function(mu, sigma, k = 15, alone = F)
{
  p <- length(mu)
  if(p != 2)
    stop("two-dimensional case only")
  if(any(unique(dim(sigma)) != p))
    stop("mu and sigma are incompatible")
  ev <- eigen(sigma, symmetric = T)
  s <- sqrt(rev(sort(ev$values)))	# need descending order
  V <- t(ev$vectors[, rev(order(ev$values))])
  theta <- (0:k) * (pi/(2 * k))
  x <- s[1] * cos(theta)
  y <- s[2] * sin(theta)
  xy <- cbind(c(x,  - x,  - x, x), c(y, y,  - y,  - y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
  if(alone) {
    xymin <- apply(xy, 2, FUN = "min")
    xymax <- apply(xy, 2, FUN = "max")
    r <- ceiling(max(xymax - xymin)/2)
    xymid <- (xymin + xymax)/2
    plot(xy[, 1], xy[, 2], xlim = c( - r, r) + xymid[1], 
	 ylim = c( - r, r) + xymid[2], xlab = "x", ylab = "y", type = "n")
  }
  l <- length(x)
  i <- 1:l
  for(k in 1:4) {
    lines(xy[i,  ])
    i <- i + l
  }
#  if(F) {
#    y <- seq(from = 0, to = s[2], by = s[2]/(2^k))
#    x <- s[1] * sqrt(1 - (y/s[2])^2)
#    xy <- cbind(c(x,  - x,  - x, x), c(y, y,  - y,  - y))
#    xy <- xy %*% V
#    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
#    l <- length(x)
#    i <- 1:l
#    for(k in 1:4) {
#      lines(xy[i,  ])
#      i <- i + l
#    }
#  }
					# semi-major axes
  P <- cbind(c( - s[1], s[1]), c(0, 0))
  P <- P %*% V
  P <- sweep(P, 2, mu, FUN = "+")
  lines(P, lty = 2)
  P <- cbind(c(0, 0), c( - s[2], s[2]))
  P <- P %*% V
  P <- sweep(P, 2, mu, FUN = "+")
  lines(P, lty = 2)
  points(mu[1], mu[2], pch = "*")
  invisible()
}

"one.XI" <- function(data)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  temp <- .Fortran("onexi",
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   double(p),
		   double(1),
		   double(1))[c(4:6)]
  loglik <- temp[[3]]
  temp[[3]] <- NULL
  names(temp) <- c("mu", "sigsq")
  sigsq <- temp$sigsq
  attr(temp, "loglik") <- loglik
  attr(temp, "rcond") <- sigsq/(1 + sigsq)
  temp
}

"one.XXX" <- function(data)
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  temp <- .Fortran("onexxx",
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   double(p),
		   double(p * p),
		   double(1),
		   double(1))[c(4:7)]
  loglik <- temp[[3]]
  rc <- temp[[4]]
  temp[[3]] <- temp[[4]] <- NULL
  temp[[2]] <- matrix(temp[[2]], p, p)
  names(temp) <- c("mu", "sigma")
  attr(temp, "loglik") <- loglik
  attr(temp, "rcond") <- rc
  temp
}

"partconv" <- function(x, consec = T)
{
  n <- length(x)
  y <- numeric(n)
  u <- unique(x)
  if(consec) {
    ## number groups in order of first row appearance
    l <- length(u)
    for(i in 1:l)
      y[x == u[i]] <- i
  }
  else {
    ## represent each group by its lowest-numbered member
    for(i in u) {
      l <- x == i
      y[l] <- (1:n)[l][1]
    }
  }
  y
}

"partuniq" <- function(x, sep = "001")
{
  ## finds the classification that removes duplicates from x
  n <- nrow(x)
  ##  x <- charconv(x)
  if(!is.data.frame(x))
    x <- data.frame(x)
  x <- do.call("paste", c(as.list(x), sep = sep))

  k <- duplicated(x)
  partition <- 1:n
  partition[k] <- match(x[k], x)
  partition
}

"pcvol" <- function(data, reciprocal = F)
{
  ## hypervolume of the data region via principal components
  data <- as.matrix(data)
  dimd <- dim(data)
  n <- dimd[1]
  p <- dimd[2]
#  if(F) {
#    vol1 <- prod(apply(data, 2, function(z)
#		       diff(range(z))))
#    V <- matrix(temp[[1]], p, p)
#    xbar <- apply(data, 2, mean)
#    X <- sweep(data, 2, xbar)
#    library(Matrix)
#    print(V)
#    print(eigen.Hermitian(crossprod(X))$vectors)
#    X <- X %*% V
#    vol <- prod(apply(X, 2, function(z)
#		      diff(range(z))))
#  }
  lwgesvd <- max(3 * min(n, p) + max(n, p), 5 * min(n, p) - 4)	# min
  lwsyevd <- p * (3 * p + 2 * ceiling(log(p, base = 2)) + 5) + 1	# min
  lisyevd <- 5 * p + 2	# minimum
  lwsyevx <- 8 * p	# minimum
  lisyevx <- 5 * p + p
  lwork <- max(lwsyevd, lwsyevx, n)
  liwork <- lisyevx
  temp <- .Fortran("mclvol",
		   as.double(data),
		   as.integer(n),
		   as.integer(p),
		   double(p),
		   double(p * p),
		   double(p * p),
		   double(lwork),
		   as.integer(lwork),
		   integer(liwork),
		   as.integer(liwork),
		   integer(1))[c(4, 11)]
  if(temp[[2]])
    stop("problem in computing principal components")
  if(reciprocal)
    prod(1/temp[[1]])
  else prod(temp[[1]])
}

"plot.emclust" <- function(x, xlab="number of clusters", ylab="BIC",
                           pch=symbols, ...)
{
  BIC <- as.matrix(x)
  n <- nrow(BIC)
  symbols <- if(n <= 9) as.character(1:n) else LETTERS[1:n]
  xrange <- if(!is.null(dn <- dimnames(BIC)[[2]])) as.numeric(dn) else 1:
    ncol(BIC)
###  plot(xrange, BIC[1,  ], type = "n", 
###       ylim = range(as.vector(BIC[!is.na(BIC)])), 
###       xlim = range(xrange), xlab = "number of clusters", 
###       ylab = "BIC")
###  for(i in 1:nrow(BIC)) {
###    points(xrange, BIC[i,  ], pch = symbols[i])
###    lines(xrange, BIC[i,  ], lty = i)
###  }
  matplot(xrange, t(BIC), type="b", xlab=xlab, ylab = ylab, pch=pch, ...)
  invisible()
}

"plot.emclust1" <- function(x, xlab = "number of clusters", ylab="BIC", ...)
{
  N <- as.numeric(names(x))
  plot(N, x, xlab=xlab, ylab = ylab, ...)
  invisible()
}

"print.bic" <- function(x, ...)
{
  print(as.vector(x), ...)	
  ##	cat("\n reciprocal condition estimate:\n")
  ##	print(attr(x, "rcond"))
  ##	cat("\n model:\n")
  ##	print(attr(x, "model"))
  invisible()
}

"print.emclust" <-
function(x, ...)
{
  modelid <- dimnames(x)[[1]]
  equal <- attr(x, "equal")
  noise <- !is.null(attr(x, "noise"))
  subset <- !is.null(attr(x, "subset"))
  attr(x, "tree") <- attr(x, "subset") <- NULL
  attr(x, "noise") <- attr(x, "Vinv") <- NULL
  attr(x, "equal") <- attr(x, "rcond") <- attr(x, "class") <- NULL
  cat("\n BIC:\n")
  NextMethod("print")
  cat("\n")
  print(c(sample = subset, noise = noise, equal = equal), ...)
  invisible()
}

"print.emclust1" <- function(x, ...)
{
  modelid <- attr(x, "modelid")
  equal <- attr(x, "equal")
  noise <- !is.null(attr(x, "noise"))
  subset <- !is.null(attr(x, "subset"))
  attr(x, "modelid") <- attr(x, "tree") <- attr(x, "subset") <- NULL
  attr(x, "noise") <- attr(x, "Vinv") <- NULL
  attr(x, "equal") <- attr(x, "rcond") <- attr(x, "class") <- NULL
  cat("\n BIC:\n")
  NextMethod("print")
  cat("\n")
  M <- c(HC = switch(modelid[1],
	   EI = "uniform spherical",
	   VI = "spherical",
	   EEE = "uniform variance",
	   VVV = "unconstrained variance",
	   stop("invalid model id for HC")), 
	 EM = switch(as.character(
	   modelid[2]),
	   EI = "uniform spherical",
	   VI = "spherical",
	   EEE = "uniform variance",
	   VVV = "unconstrained variance",
	   EEV = "uniform shape and volume",
	   VEV = "uniform shape",
	   stop("invalid model id for EM")))
  if(subset)
    M["HC"] <- paste(M["HC"], "(on a sample)")
  if(noise)
    M["EM"] <- paste(M["EM"], "(with noise)")
  if(equal)
    M["EM"] <- paste(M["EM"], "(equal mixing proportions)")
  print(M, ...)
  invisible()
}

#"print.mclust" <- function(x, ...)
#{
#  attributes(x) <- if(!is.null(dim(x))) attributes(x)[c("dim", "model")]
#  else NULL
#  NextMethod("print")
#}

"print.mhtree" <- function(x, ...)
{
  attributes(x) <- if(!is.null(dim(x))) attributes(x)[c("dim", "model")]
  else NULL
  NextMethod("print")
}

"print.summary.emclust" <- function(x, ...)
{
  bic <- attr(x, "bic")
  l <- length(bic) > 1
  if(l)
    cat("\n best classification:\n")
  else cat("\n classification:\n")
  print(x$classification, ...)
  cat("\n uncertainty (quantiles):\n")
  print(quantile(x$uncertainty))
  if(l)
    cat("\n best BIC values:\n")
  else cat("\n BIC value:\n")
  print(bic)	
  ##	cat("\n reciprocal condition estimates:\n")
  ##	print(attr(x, "rcond"))
  M <- switch(attr(x, "modelid"),
	      EI = "uniform spherical",
	      VI = "spherical",
	      EEE = "uniform variance",
	      VVV = "unconstrained variance",
	      EEV = "uniform shape and volume",
	      VEV = "uniform shape",
	      stop("invalid model id for EM"))
  cat("\n best model:", M, "\n\n")
  print(attr(x, "options"))
  invisible()
}

"print.summary.emclust1" <- function(x, ...)
{
  bic <- attr(x, "bic")
  l <- length(bic) > 1
  if(l)
    cat("\n best classification:\n")
  else cat("\n classification:\n")
  print(x$classification, ...)
  cat("\n uncertainty (quantiles):\n")
  print(quantile(x$uncertainty))
  if(l)
    cat("\n best BIC values:\n")
  else cat("\n BIC value:\n")
  print(bic)	
  ##	cat("\n reciprocal condition estimates:\n")
  ##	print(attr(x, "rcond"))
  cat("\n model:\n")
  M <- c(HC = switch(attr(x, "modelid")[1],
	   EI = "uniform spherical (EI)",
	   VI = "spherical (VI)",
	   EEE = "uniform variance (EEE)",
	   VVV = "unconstrained variance (VVV)",
	   stop("invalid model id for EM")), 
	 EM = switch(attr(x, "model")[2],
	   EI = "uniform spherical(EI)",
	   VI = "spherical (VI)",
	   EEE = "uniform variance (EEE)",
	   VVV = "unconstrained variance (VVV)",
	   EEV = "uniform shape and volume (EEV)",
	   VEV = "uniform shape (VEV)",
	   stop("invalid model id for EM")))
  print(M)
  cat("\n")
  print(attr(x, "options"))
  invisible()
}

"summary.emclust" <- function(x, data, nclus, modelid)
{
  rc <- attr(x, "rcond")
  tree <- attr(x, "tree")
  if(missing(modelid))
    modelid <- dimnames(x)[[1]]
  smpl <- attr(x, "subset")
  equal <- attr(x, "equal")
  noise <- attr(x, "noise")
  Vinv <- attr(x, "Vinv")
  attr(x, "tree") <- attr(x, "subset") <- attr(x, "noise") <- 
    attr(x, "Vinv") <- attr(x, "equal") <- attr(x, "rcond") <- 
      attr(x, "call") <- attr(x, "class") <- NULL
  n <- nrow(data)
  nclus <- 
    if(missing(nclus)) dimnames(x)[[2]] 
    else as.character(sort(unique(nclus)))
  x <- x[modelid, nclus, drop = F]
  rc <- rc[modelid, nclus, drop = F]
  if(all(is.na(x))) {
    warning("selected BIC values are all missing")
    return(structure(rep(NA, n), modelid = modelid, 
		     options = c(sample = !is.null(smpl), 
		       noise = !is.null(noise), equal = equal), 
		     class = "summary.emclust"))
  }
  x[is.na(x)] <-  - Inf #.Machine$double.xmax
  bicmax <- max(x)
  n <- nrow(x)
  l <- ncol(x)
  if(min(n, l) > 1) {
    I <- matrix(rep(1:n, l), n, l)
    i <- I[x == bicmax][1]
    j <- nclus[x[i,  ] == bicmax][1]
    other <- if(any(!as.numeric(nclus))) max(x[ - i, -1]) else max(x[- i,])
    best <- c(bicmax, max(x[i, nclus != j]), other)
    J <- matrix(rep(nclus, n - 1), n - 1, l, byrow = T)
    j <- J[x[ - i,  ] == other][1]
    K <- matrix(rep(1:n, l), n, l)[ - i,  ]
    k <- K[x[ - i,  ] == other][1]
    rows <- dimnames(x)[[1]][c(i, i, k)]
    cols <- c(nclus[x[i,  ] == best[1]], nclus[x[i,  ] == best[2]], j)
    rcond <- c(rc[rows[1], cols[1]], 
	       rc[rows[2], cols[2]], 
	       rc[rows[3], cols[3]])
    names(best) <- names(rcond) <- paste(rows, ",", cols[1:3], sep = "")
    modelid <- rows[1]
    k <- cols[1]
  }
  else if(l != 1) {
    dn <- dimnames(x)
    x <- as.vector(x)
    i <- (1:l)[x == bicmax][1]
    nextbest <- max(x[ - i])
    j <- 
      if(nextbest == bicmax) (1:l)[x == bicmax][2] 
      else (1:l)[x == nextbest][1]
    best <- c(bicmax, nextbest)
    rcond <- as.vector(rc)[c(i, j)]
    modelid <- dn[[1]]
    names(best) <- names(rcond) <- 
      paste(modelid, ",", dn[[2]][c(i, j)], sep = "")
    k <- as.character(i)
  }
  else if(n != 1) {
    dn <- dimnames(x)
    x <- as.vector(x)
    i <- (1:n)[x == bicmax][1]
    nextbest <- max(x[ - i])
    j <- 
      if(nextbest == bicmax) (1:n)[x == bicmax][2] 
      else (1:n)[x == nextbest][1]
    best <- c(bicmax, nextbest)
    rcond <- as.vector(rc)[c(i, j)]
    names(best) <- names(rcond) <-
      paste(dn[[1]][c(i, j)], ",", dn[[2]], sep = "")
    modelid <- dn[[1]][i]
    k <- nclus
  }
  else {
    dn <- dimnames(x)
    k <- nclus
    rcond <- rc
    best <- bicmax
    names(best) <- names(rcond) <- paste(dn[[1]], ",", dn[[2]], sep = "")
  }
  if(is.null(smpl)) {
    if(is.null(noise)) {
      if(k == "1") {
	z <- matrix(1, nrow(data), 1)
	mu <- apply(data, 2, mean)
	params <- c(mu = mu, sigma = crossprod(sweep(
			       data, 2, mu)))
      }
      else {
	cl <- mhclass(tree, as.numeric(k))
	z <- me(data, modelid = modelid, ctoz(cl), 
		equal = equal)
	params <- mstep(data, modelid = modelid, z, 
			equal = equal)[c("mu", "sigma", "prob")]
      }
    }
    else {               # noise
      if(k == "0") {
	z <- cbind(rep(0, n), rep(1, n))
	params <- NULL
      }
      else {
	cl <- numeric(n)
	cl[!noise] <- mhclass(tree, as.numeric(k))
	cl[noise] <- as.numeric(k) + 1
	z <- me(data, modelid = modelid, ctoz(cl), 
		equal = equal, noise = T, Vinv = Vinv)
	params <- mstep(data, modelid = modelid, z, 
			equal = equal, noise = T, Vinv = Vinv)
      }
    }
  }
  else {
    ## a sample was used in the hierarchical clustering phase
    if(is.null(noise)) {
      if(k == "1") {
	z <- matrix(1, nrow(data), 1)
	mu <- apply(data, 2, mean)
	params <- c(mu = mu, sigma = crossprod(sweep(data, 2, mu)))
      }
      else {
	cl <- mhclass(tree, as.numeric(k))
	params <- mstep(data[smpl,  ], modelid = modelid, ctoz(cl), 
			equal = equal)[c("mu", "sigma", "prob")]
	z <- do.call("estep", c(list(data, modelid = modelid), params))
	z <- me(data, modelid = modelid, z, equal = 
		equal)
	params <- mstep(data, modelid = modelid, z, 
			equal = equal)[c("mu", "sigma", "prob")]
      }
    }
    else {                # noise
      if(k == "0") {
	z <- cbind(rep(0, n), rep(1, n))
	params <- NULL
      }
      else {
	cl <- mhclass(tree, as.numeric(k))
	params <- mstep(data[smpl,  ], modelid = 
			modelid, ctoz(cl), equal = equal)[c("mu", 
					     "sigma", "prob")]
	z <- do.call("estep", c(list(data[!noise,  ], 
				     modelid = modelid), params))
	cl <- numeric(n)
	cl[!noise] <- z
	cl[noise] <- as.numeric(k) + 1
	z <- me(data, modelid = modelid, ctoz(cl), 
		equal = equal, noise = T, Vinv = Vinv)
	params <- mstep(data, modelid = modelid, z, 
			equal = equal, noise = T, Vinv = Vinv)
      }
    }
  }
  out <- list(classification = ztoc(z), uncertainty = 1 - apply(z, 1, max
					  ), parameters = params, z = z)
  attr(out, "modelid") <- modelid
  attr(out, "options") <- c(sample = !is.null(smpl), 
			    noise = !is.null(noise), equal = equal)	
  ##---------------------------------------------------------------------------
  attr(out, "bic") <- best
  attr(out, "rcond") <- rc
  class(out) <- "summary.emclust"
  out
}

"summary.emclust1" <- function(x, data, nclus)
{
  rc <- attr(x, "rcond")
  tree <- attr(x, "tree")
  modelid <- attr(x, "modelid")
  smpl <- attr(x, "subset")
  equal <- attr(x, "equal")
  noise <- attr(x, "noise")
  Vinv <- attr(x, "Vinv")
  attr(x, "modelid") <- attr(x, "tree") <- attr(x, "subset") <- 
    attr(x, "noise") <- attr(x, "Vinv") <- attr(x, "equal") <- 
      attr(x, "rcond") <- attr(x, "call") <- attr(x, "class") <- NULL 
  n <- nrow(data)
  nclus <- 
    if(missing(nclus)) names(x) 
    else as.character(sort(unique(nclus)))
  x <- x[nclus]
  if(all(is.na(x))) {
    warning("selected BIC values are all missing")
    return(structure(rep(NA, n), modelid = modelid, 
		     options = c(sample = !is.null(smpl), 
		       noise = !is.null(noise), equal = equal), 
		     class = "summary.emclust1"))
  }
  nclus <- nclus[!is.na(x)]
  x <- x[!is.na(x)]
  bicmax <- max(x)
  k <- nclus[j <- ((1:length(x))[x == bicmax][1])]
  if(is.null(smpl)) {
    if(is.null(noise)) {
      if(k == "1") {
	z <- matrix(1, nrow(data), 1)
	mu <- apply(data, 2, mean)
	params <- c(mu = mu, sigma = crossprod(sweep(
			       data, 2, mu)))
      }
      else {
	cl <- mhclass(tree, as.numeric(k))
	z <- me(data, modelid = modelid[2], ctoz(cl), 
		equal = equal)
	params <- mstep(data, modelid = modelid[2], z, 
			equal = equal)[c("mu", "sigma", "prob")]
      }
    }
    else {           # noise
      if(k == "0") {
	z <- cbind(rep(0, n), rep(1, n))
	params <- NULL
      }
      else {
	cl <- numeric(n)
	cl[!noise] <- mhclass(tree, as.numeric(k))
	cl[noise] <- as.numeric(k) + 1
	z <- me(data, modelid = modelid[2], ctoz(cl), 
		equal = equal, noise = T, Vinv = Vinv)
	params <- mstep(data, modelid = modelid[2], z, 
			equal = equal, noise = T, Vinv = Vinv)
      }
    }
  }
  else {
    ## a sample was used in the hierarchical clustering phase
    if(is.null(noise)) {
      if(k == "1") {
	z <- matrix(1, nrow(data), 1)
	mu <- apply(data, 2, mean)
	params <- c(mu = mu, sigma = crossprod(sweep(
			       data, 2, mu)))
      }
      else {
	cl <- mhclass(tree, as.numeric(k))
	params <- mstep(data[smpl,  ], modelid = 
			modelid[2], ctoz(cl), equal = equal)[c("mu", 
						"sigma", "prob")]
	z <- do.call("estep", c(list(data, modelid = 
				     modelid[2]), params))
	z <- me(data, modelid = modelid[2], z, equal = 
		equal)
	params <- mstep(data, modelid = modelid[2], z, 
			equal = equal)[c("mu", "sigma", "prob")]
      }
    }
    else {               # noise
      if(k == "0") {
	z <- cbind(rep(0, n), rep(1, n))
	params <- NULL
      }
      else {
	cl <- mhclass(tree, as.numeric(k))
	params <- mstep(data[smpl,  ], 
			modelid = modelid[2], ctoz(cl), 
			equal = equal)[c("mu", "sigma", "prob")]
	z <- do.call("estep", c(list(data[!noise,  ], 
				     modelid = modelid[2]), params))
	cl <- numeric(n)
	cl[!noise] <- z
	cl[noise] <- as.numeric(k) + 1
	z <- me(data, modelid = modelid[2], ctoz(cl), 
		equal = equal, noise = T, Vinv = Vinv)
	params <- mstep(data, modelid = modelid[2], z, 
			equal = equal, noise = T, Vinv = Vinv)
      }
    }
  }
  out <- list(classification = ztoc(z), 
	      uncertainty = 1 - apply(z, 1, max), 
	      parameters = params, z = z) 
  attr(out, "modelid") <- modelid
  attr(out, "options") <- c(sample = !is.null(smpl), noise =
			    !is.null(noise), equal = equal)	 
##-----------------------------------------------------------------------------
  if(length(nclus) > 1) {
    nextbest <- max(x[ - j])
    indx <- c(nclus[j], (nclus[x == nextbest])[1])
    rc <- rc[indx]
    best <- c(bicmax, nextbest)
    names(best) <- names(rc) <- indx
  }
  else {
    best <- bicmax
    rc <- rc[j]
    names(best) <- names(rc) <- nclus
  }
  attr(out, "bic") <- best
  attr(out, "rcond") <- rc
  class(out) <- "summary.emclust1"
  out
}

"traceW" <- function(x)
{
### sum(as.vector(sweep(x, 2, apply(x, 2, mean)))^2)
  dimx <- dim(x)
  n <- dimx[1]
  p <- dimx[2]
  .Fortran("mcltrw",
	   as.double(x),
	   as.integer(n),
	   as.integer(p),
	   double(p),
	   double(1))[[5]]
}

"ztoc" <- function(z)
{
### converts conditional probabilities to a classification
  cl <- numeric(nrow(z))
  for(i in 1:nrow(z)) {
    cl[i] <- (1:ncol(z))[z[i,  ] == max(z[i,  ])]
  }
  cl
}


