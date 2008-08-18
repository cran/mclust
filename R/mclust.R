".First.lib" <- function(lib, pkg) {
  library.dynam("mclust", pkg, lib)
}
"cdensE" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- drop(data)
	n <- length(data)
        mu <- drop(parameters$mean)
	G <- length(mu)
        if(any(is.na(unlist(parameters[c("mean", "variance")]))) ||
            any(is.null(parameters[c("mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "E",
                       WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
        if(is.null(sigmasq))
                stop("variance parameters are missing")
        if(length(sigmasq) > 1)
                warning("more than one sigma-squared given")
        if(sigmasq < 0)
                stop("sigma-squared is negative")
        if(!sigmasq) {
                WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "E",
                                 WARNING = WARNING, returnCode = 9))
        }
	temp <- .Fortran("es1e",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(-1),
		as.integer(n),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[8:9]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        } 
        dimnames(z) <- list(names(data),NULL)
	structure(z, logarithm = logarithm, modelName = "E", 
                     WARNING = WARNING, returnCode = ret) 
}

"emE" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepE(data, parameters = parameters, warn = warn)$z
	meE(data, z = z, prior = prior, control = control, 
            Vinv = parameters$Vinv, warn = warn)
}

"estepE" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- drop(data)
	n <- length(data)
        pro <- parameters$pro
        pro <- pro/sum(pro)
	l <- length(pro)

        mu <- drop(parameters$mean)
	G <- length(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
	if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
		WARNING <- "parameters are missing"
		if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
		return(structure(list(modelName = "E", n=n, d=1, G=G, z=z,
 		       parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
	}
        sigmasq <- parameters$variance$sigmasq
	if(is.null(sigmasq))
		stop("variance parameters are missing")
	if(length(sigmasq) > 1)
		warning("more than one sigma-squared specified")
	if(sigmasq < 0)
		stop("sigma-squared is negative")
	if(!sigmasq) {
		WARNING <- "sigma-squared vanishes"
		if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
		return(structure(list(modelName = "E", n=n, d=1, G=G, z=z,
 		       parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
	}
	temp <- .Fortran("es1e",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(pro),
		as.integer(n),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[8:9]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(names(data),NULL) 
	structure(list(modelName = "E", n = n, d = 1, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"hcE" <-
function(data, partition, minclus = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(minclus < 1) stop("minclus must be positive")
	if(any(is.na(data)))
		stop("missing values not allowed in data")
	#====================================================================
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- as.vector(data)
	n <- length(data)
	if(missing(partition))
		partition <- 1:n
	else if(length(partition) != n)
		stop("partition must assign a class to each observation")
	partition <- partconv(partition, consec = TRUE)
	l <- length(unique(partition))
	attr(partition, "unique") <- l
	m <- l - minclus
	if(m <= 0)
		stop("initial number of clusters is not greater than minclus")
	storage.mode(data) <- "double"
	ld <- max(c((l * (l - 1))/2, 3 * m))
	temp <- .Fortran("hc1e",
		data,
		as.integer(n),
		as.integer(partition),
		as.integer(l),
		as.integer(m),
		as.integer(ld),
		double(ld),
                PACKAGE = "mclust")[c(1, 3, 7)]
	temp[[1]] <- temp[[1]][1:m]
	temp[[2]] <- temp[[2]][1:m]
	temp[[3]] <- temp[[3]][1:m]
        change <- temp[[3]]
	structure(rbind(temp[[1]], temp[[2]]), 	initialPartition = partition, 
                  dimensions = n, modelName = "E",
		  call = match.call())
}

"meE" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be 1 dimensional")
	data <- as.vector(data)
	n <- length(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal length of data")
	K <- dimz[2]
        if (!is.null(Vinv)) {
		G <- K - 1
 		if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
	       WARNING <- "z is missing"
	       if (warn) warning(WARNING)
               variance <- list(modelName = "E", d = 1, G = G, sigmasq = NA)
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="E", prior=prior, n=n, d=1, G=G,
                      z=z, parameters=parameters, control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("me1e",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(G),
			double(1),
			double(K),
                        PACKAGE = "mclust")[6:12]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "E"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("me1ep",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(G),
			double(1),
			double(K),
                        PACKAGE = "mclust")[c(10:16, 9)]
	}
	mu <- temp[[5]]
	names(mu) <- as.character(1:G)
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	sigmasq <- temp[[6]]
	pro <- temp[[7]]
	## log post <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6) || 
                            sigmasq <= max(control$eps,0)) {
		WARNING <- "sigma-squared falls below threshold"
		if(warn) warning(WARNING)
		mu[] <- pro[] <- sigmasq <- z[] <- loglik <- logprior <- NA
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
		ret <- if(control$equalPro) -2 else -3
	}
	else if(its >= control$itmax[1]) {
		WARNING <- "iteration limit reached"
		warning(WARNING)
		its <-  - its
		ret <- 1
	}
	else ret <- 0
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(names(data), NULL) 
        variance <- list(modelName = "E", d = 1, G = G, sigmasq = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "E", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)

}

"mstepE" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- as.vector(data)
	n <- length(data)
	##
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
# number of groups 
        G <- dimz[2]
	##
	if(all(is.na(z))) {
               WARNING <- "z is missing"
	       if (warn) warning(WARNING)
               variance <- list(modelName="E", d=1, G=G, sigmasq=NA)
               parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="E", prior=prior, n=n, d=1, G=G,
                                     z = z, parameters=parameters), 
                                WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("ms1e",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(G),
			double(G),
			double(1),
			double(G),
                        PACKAGE = "mclust")[5:7]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "E"), prior[names(prior) !=
			"functionName"]))
		storage.mode(z) <- "double"
		temp <- .Fortran("ms1ep",
			as.double(data),
			z,
			as.integer(n),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(G),
			double(1),
			double(G),
                        PACKAGE = "mclust")[9:11]
	}
	mu <- temp[[1]]
	names(mu) <- as.character(1:G)
	sigmasq <- temp[[2]]
	pro <- temp[[3]]
	WARNING <- NULL
	if(sigmasq > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
                pro[] <- mu[] <- sigmasq <- NA
                ret <- -1
 
	}
        else ret <- 0
        dimnames(z) <- list(names(data), NULL) 
        variance <- list(modelName = "E", d = 1, G = G, sigmasq = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "E", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

"simE" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, 2), modelName = "E"))
  }
  if(!is.null(seed))
    set.seed(seed)
  mu <- parameters$mean
  G <- length(mu)
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- rep(0, n)
  sd <- sqrt(parameters$variance$sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd)
  }
  structure(cbind(group = clabels, "1" = x), modelName = "E")
}

"cdensEEE" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) > 2)
		stop("data must be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EEE", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholSigma))
           stop("variance parameters are missing")
	temp <- .Fortran("eseee",
		as.logical(1),
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$cholSigma),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(p),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[10:12]
	lapackCholInfo <- temp[[1]][1]
	loglik <- temp[[2]]
	z <- matrix(temp[[3]], n, G)
        WARNING <- NULL
	if(lapackCholInfo) {
		if(lapackCholInfo > 0) {
			WARNING <- "sigma is not positive definite"
			if (warn) warning(WARNING)
		}
		else {
			WARNING <- "input error for LAPACK DPOTRF"
			if (warn) warning(WARNING)
		}
		z[] <- NA
                ret <- -9 
	}
	else if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
	structure(z, logarithm = logarithm, modelName = "EEE",
                  WARNING = WARNING, retrunCode = ret)
}

"emEEE" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepEEE(data, parameters = parameters, warn = warn)$z  
	meEEE(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepEEE" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) > 2)
		stop("data must be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EEE", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholSigma))
           stop("variance parameters are missing")
	temp <- .Fortran("eseee",
		as.logical(1),
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$cholSigma),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(p),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[10:12]
	lapackCholInfo <- temp[[1]][1]
	loglik <- temp[[2]]
	z <- matrix(temp[[3]], n, K)
	WARNING <- NULL
	if(lapackCholInfo) {
		if(lapackCholInfo > 0) {
			WARNING <- "sigma is not positive definite"
			warning(WARNING)
                        ret <- -4 
		}
		else {
			WARNING <- "input error for LAPACK DPOTRF"
			warning(WARNING)
                        ret <- -5
		}
		z[] <- loglik <- NA
	}
	else if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EEE", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"hcEEE" <-
function(data, partition, minclus = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(minclus < 1) stop("minclus must be positive")
	if(any(is.na(data)))
		stop("missing values not allowed in data")
	#=====================================================================
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	dimnames(data) <- NULL
	n <- nrow(data)
	p <- ncol(data)
	if(n <= p)
		warning("# of observations <= data dimension")
	if(missing(partition))
		partition <- 1:n
	else if(length(partition) != n)
		stop("partition must assign a class to each observation")
	partition <- partconv(partition, consec = TRUE)
	l <- length(unique(partition))
	attr(partition, "unique") <- l
	m <- l - minclus
	if(m <= 0)
		stop("initial number of clusters is not greater than minclus")
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
		double(p * p),
                PACKAGE ="mclust")[c(1, 7:10)]
	#
	# currently temp[[5]] is not output
	temp[[4]] <- temp[[4]][1:2]
	temp[[5]] <- temp[[5]][1:2]
	names(temp[[5]]) <- c("determinant", "trace")
	temp[[1]] <- temp[[1]][1:(m + 1),  ]
	if(p < 3)
		tree <- rbind(temp[[2]], temp[[3]])
	else if(p < 4)
		tree <- rbind(temp[[1]][-1, 3], temp[[3]])
	else tree <- t(temp[[1]][-1, 3:4, drop = FALSE])
	determinant <- temp[[1]][, 1]
	attr(determinant, "breakpoints") <- temp[[4]]
        trace <- temp[[1]][, 2]
	structure(tree,	initialPartition = partition, 
                  dimensions = dimdat, modelName = "EEE", 
                  call = match.call())
}

"meEEE" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
        if (!is.null(Vinv)) {
                G <- K - 1
                if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
        }
        else G <- K
	if(all(is.na(z))) {
              WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEE", d = p, G = G, 
                     Sigma = matrix(NA, p, p), cholSigma = matrix(NA, p, p)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EEE", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("meeee",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(p * p),
			double(K),
			double(p),
                        PACKAGE = "mclust")[7:13]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EEE"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("meeeep",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(p * p),
			double(K),
			double(p),
                        PACKAGE = "mclust")[c(11:17, 10)]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
   	cholSigma <- matrix(temp[[6]], p, p)
	pro <- temp[[7]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
                Sigma <- matrix( NA, p, p)
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn)
				warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn)
				warning(WARNING)
		}
		mu[] <- pro[] <- z[] <- loglik <- logprior <- NA
		sigma <- array(NA, c(p, p, G))
                Sigma <- matrix(NA, p, p)
		ret <- if(control$equalPro) -2 else -3
	}
	else {
   	        Sigma <- unchol(cholSigma, upper = TRUE)
		sigma <- array(0, c(p, p, G))
		for(k in 1:G)
			sigma[,  , k] <- Sigma
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- dimnames(cholSigma) <- 
                       list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EEE", d = p, G = G,
                         sigma = sigma, Sigma = Sigma, cholSigma = cholSigma) 
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "EEE", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepEEE" <-
function(data, z, prior = NULL,  warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	##
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEE", d = p, G = G, 
                      sigma <- array(NA, c(p,p, G)), 
                      Sigma = matrix(NA, p, p), cholSigma = matrix(NA, p, p)) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EEE", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("mseee",
                        as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p),
			double(p * G),
			double(p * p),
			double(G),
                        PACKAGE = "mclust")[7:9]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EEE"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("mseeep",
                        as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			double(p),
			double(p * G),
			double(p * p),
			double(G),
                        PACKAGE = "mclust")[11:13]
	}
	mu <- matrix(temp[[1]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	cholSigma <- matrix(temp[[2]], p, p)
	pro <- temp[[3]]
	sigma <- array(0, c(p, p, G))
	Sigma <- unchol(cholSigma, upper = TRUE)
	for(k in 1:G)
		sigma[,  , k] <- Sigma
	WARNING <- NULL
	if(any(mu > signif(.Machine$double.xmax, 6))) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
		mu[] <- sigma[] <- Sigma[] <- cholSigma[] <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- dimnames(cholSigma) <- 
               list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EEE", d = p, G = G, 
                         sigma = sigma, Sigma = Sigma, cholSigma= cholSigma)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEE", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

"simEEE" <- 
function(parameters, n, seed = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(!is.null(seed)) set.seed(seed)
	mu <- as.matrix(parameters$mean)
	d <- nrow(mu)
	G <- ncol(mu)
	if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
		parameters[c("mean", "variance")]))) {
		warn <- "parameters are missing"
		warning("parameters are missing")
		return(structure(matrix(NA, n, d + 1), modelName = "EEE"))
	}
	pro <- parameters$pro
	if(is.null(pro))
		pro <- rep(1/G, G)
	clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
	ctabel <- table(clabels)
	x <- matrix(0, n, d)
	if(is.null(cholSigma <- parameters$variance$cholSigma)) {
		if(is.null(Sigma <- parameters$variance$Sigma)) {
			stop("variance parameters must inlcude either Sigma or cholSigma"
				)
		}
		cholSigma <- chol(Sigma)
	}
	for(k in 1:G) {
		m <- ctabel[k]
		x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
			ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
			FUN = "+")
	}
	dimnames(x) <- list(NULL, 1:d)
	structure(cbind(group = clabels, x), modelName = "EEE")
}

"cdensEEI" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EEI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
	temp <- .Fortran("eseei",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[10:11]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1 
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(z, logarithm = logarithm, modelName = "EEI",
                  WARNING = WARNING, returnCode = ret)
}

"emEEI" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepEEI(data, parameters = parameters, warn = warn)$z  
	meEEI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepEEI" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
                Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	} 
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EEI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
		   is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
	temp <- .Fortran("eseei",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EEI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"meEEI" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data  should be in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEI", d = p, G = G, 
                                 scale = NA, shape = rep(NA,p)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("meeei",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(1),
			double(p),
			double(K),
                        PACKAGE = "mclust")[7:14]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EEI"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("meeeip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(1),
			double(p),
			double(K),
                        PACKAGE = "mclust")[11:18]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	scale <- temp[[6]]
	shape <- temp[[7]]
	pro <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		sigma <- array(NA, c(p, p, G))
		Sigma <- matrix(NA, p, p)
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn)
				warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn)
				warning(WARNING)
		}
		sigma <- array(NA, c(p, p, G))
		Sigma <- matrix(NA, p, p)
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(0, c(p, p, G))
		Sigma <- diag(scale * shape)
		for(k in 1:G)
			sigma[,  , k] <- Sigma
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
	variance <- list(modelName = "EEI", d = p, G = G, 
                         sigma = sigma, Sigma = Sigma, 
                         scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
      	structure(list(modelName = "EEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
		  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepEEI" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEI", d = p, G = G, 
                                 scale = NA, shape = rep(NA,p)) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))

	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("mseei",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p * G),
			double(1),
			double(p),
			double(G),
                        PACKAGE = "mclust")[6:9]
	}
	else {
		storage.mode(z) <- "double"
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EEI"), prior[names(
			prior) != "functionName"]))
		temp <- .Fortran("mseeip",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p * G),
			double(1),
			double(p),
			double(G),
                        PACKAGE = "mclust")[10:13]
	}
	mu <- matrix(temp[[1]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	scale <- temp[[2]]
	shape <- temp[[3]]
	pro <- temp[[4]]
	WARNING <- NULL
	if(any(c(shape, scale) > signif(.Machine$double.xmax, 6)) || any(!c(
		scale, shape))) {
		WARNING <- "cannot compute M-step"
		if(warn)
			warning(WARNING)
		mu[] <- pro[] <- scale <- shape[] <- NA
		sigma <- Sigma <- array(NA, c(p, p, G))
                ret <- -1
	}
	else {
		sigma <- array(0, c(p, p, G))
		Sigma <- diag(scale * shape)
		for(k in 1:G)
			sigma[,  , k] <- Sigma
                ret <- 0
	}
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EEI", d = p, G = G, 
                         sigma = sigma, Sigma = Sigma, 
                         scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

"simEEI" <- 
function(parameters, n, seed = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(!is.null(seed)) set.seed(seed)
	mu <- as.matrix(parameters$mean)
	d <- nrow(mu)
	G <- ncol(mu)
	if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
		parameters[c("mean", "variance")]))) {
		warn <- "parameters are missing"
		warning("parameters are missing")
		return(structure(matrix(NA, n, d + 1), modelName = "EEI"))
	}
	pro <- parameters$pro
	if(is.null(pro))
		pro <- rep(1/G, G)
	clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
	ctabel <- table(clabels)
	x <- matrix(0, n, d)
	shape <- parameters$variance$shape
	if(length(shape) != d)
		stop("shape incompatible with mean")
	cholSigma <- diag(sqrt(parameters$variance$scale * shape))
	for(k in 1:G) {
		m <- ctabel[k]
		x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
			ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
			FUN = "+")
	}
	dimnames(x) <- list(NULL, 1:d)
	structure(cbind(group = clabels, x), modelName = "EEI")
}

"cdensEEV" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EEV", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape) ||
                   is.null(parameters$variance$orientation)) 
          stop("variance parameters are missing")
	temp <- .Fortran("eseev",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(parameters$variance$orientation),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(p),
		double(p),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[13:14]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
	structure(z, logarithm = logarithm, modelName = "EEV",
                  WARNING = WARNING, returnCode = ret)
}

"emEEV" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepEEV(data, parameters = parameters, warn = warn)$z  
	meEEV(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepEEV" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EEV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape) ||
                   is.null(parameters$variance$orientation)) 
          stop("variance parameters are missing")
	temp <- .Fortran("eseev",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(parameters$variance$orientation),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(p),
		double(p),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[13:14]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1 
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EEV", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"meEEV" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
              WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEV", d = p, G = G, 
              scale = NA, shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
               parameters <- list(Vinv= Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))

	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p))
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("meeev",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			as.integer(lwork),
			double(p * G),
			double(1),
			double(p),
			double(p * p * G),
			double(K),
			double(lwork),
			double(p),
                        PACKAGE = "mclust")[7:16]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EEV"), 
                        prior[names(prior) !="functionName"]))
		temp <- .Fortran("meeevp",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			as.integer(lwork),
			double(p * G),
			double(1),
			double(p),
			double(p * p * G),
			double(K),
			double(lwork),
			double(p),
                        PACKAGE = "mclust")[11:20]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	lapackSVDinfo <- temp[[5]]
	mu <- matrix(temp[[6]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	scale <- temp[[7]]
	shape <- temp[[8]]
	O <- array(temp[[9]], c(p, p, G))
	pro <- temp[[10]]
	WARNING <- NULL
	if(lapackSVDinfo) {
		if(lapackSVDinfo > 0) {
			WARNING <- "LAPACK DGESVD fails to converge"
		}
		else {
			WARNING <- "input error for LAPACK DGESVD"
		}
		z[] <- O[] <- shape[] <- NA
		scale <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -9
	}
	else if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		shape[] <- NA
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "a z column sum fell below threshold"
			if(warn)
				warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn)
				warning(WARNING)
		}
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- scale * shapeO(shape, O, transpose = TRUE)
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]], 
                            NULL)
## Sigma = scale * t(O) %*% diag(shape) %*% O
	variance <- list(modelName = "EEV", d = p, G = G, sigma = sigma,
                         scale = scale, shape = shape, orientation = O) 
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance) 
	structure(list(modelName = "EEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik),
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepEEV" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	##
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EEV", d = p, G = G, 
                scale = NA, shape = rep(NA,p), orientation=array(NA,c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))

	}
	#	shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
	if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop(
			"improper specification of z")
	lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), G)
	if(is.null(prior)) {
		temp <- .Fortran("mseev",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(lwork),
			as.integer(lwork),
			double(p * G),
			double(1),
			double(p),
			double(p * p * G),
			double(G),
                        PACKAGE = "mclust")[7:12]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EEV"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("mseevp",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			double(lwork),
			as.integer(lwork),
			double(p * G),
			double(1),
			double(p),
			double(p * p * G),
			double(G),
                        PACKAGE = "mclust")[11:16]
	}
	lapackSVDinfo <- temp[[1]]
	mu <- matrix(temp[[2]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	scale <- temp[[3]]
	shape <- temp[[4]]
	O <- array(temp[[5]], c(p, p, G))
	pro <- temp[[6]]
	WARNING <- NULL
	if(lapackSVDinfo) {
		if(lapackSVDinfo > 0) {
			WARNING <- "LAPACK DGESVD fails to converge"
			warning(WARNING)
                        ret <- -4
		}
		else {
			WARNING <- "input error for LAPACK DGESVD"
			warning(WARNING)
                        ret <- -5
		}
		O[] <- shape[] <- scale <- NA
		sigma <- array(NA, c(p, p, G))
	}
	else if(any(c(abs(scale), shape) > signif(.Machine$double.xmax, 6))) {
		WARNING <- "cannot compute M-step"
		if(warn)
			warning(WARNING)
		mu[] <- pro[] <- scale <- O[] <- shape[] <- NA
		sigma <- array(NA, c(p, p, G))
                ret <- -1
	}
	else {
		sigma <- scale * shapeO(shape, O, transpose = TRUE)
                ret <- 0
	}
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(O) <- 
              list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
        variance <- list(modelName = "EEV", d = p, G = G, sigma = sigma,
                       scale = scale, shape = shape, orientation = O)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

"simEEV" <- 
function(parameters, n, seed = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(!is.null(seed)) set.seed(seed)
	mu <- as.matrix(parameters$mean)
	d <- nrow(mu)
	G <- ncol(mu)
	if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
		parameters[c("mean", "variance")]))) {
		warn <- "parameters are missing"
		warning("parameters are missing")
		return(structure(matrix(NA, n, d + 1), modelName = "EEV"))
	}
	pro <- parameters$pro
	if(is.null(pro))
		pro <- rep(1/G, G)
	clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
	ctabel <- table(clabels)
	x <- matrix(0, n, d)
	shape <- parameters$variance$shape
	if(length(shape) != d)
		stop("shape incompatible with mean")
	sss <- sqrt(parameters$variance$scale * shape)
	for(k in 1:G) {
		m <- ctabel[k]
		cholSigma <- parameters$variance$orientation[,  , k] * sss
		x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
			ncol = d) %*% cholSigma, MARGIN = 2, STAT = mu[, k],
			FUN = "+")
	}
	dimnames(x) <- list(NULL, 1:d)
	structure(cbind(group = clabels, x), modelName = "EEV")
}

"cdensEII" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EII", 
                                 WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
	if(sigmasq < 0)
		stop("sigma-squared is negative")
	if(!sigmasq) {
		WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EII", 
                                 WARNING = WARNING, returnCode = 9))
	}
	temp <- .Fortran("eseii",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[9:10]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
	structure(z, logarithm = logarithm, modelName = "EII",
                  WARNING = WARNING, returnCode = ret)
}

"emEII" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepEII(data, parameters = parameters, warn = warn)$z
	meEII(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepEII" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	p <- ncol(data)
	n <- nrow(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
        mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
                Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
	if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
		WARNING <- "parameters are missing"
		if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
		return(structure(list(modelName = "EII", n=n, d=p, G=G, z=z,
  		                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
	}
        sigmasq <- parameters$variance$sigmasq
	if(is.null(sigmasq))
		warning("variance parameters are missing")
	if(sigmasq < 0)
		stop("sigma-squared is negative")
	if(!sigmasq) {
		WARNING <- "sigma-squared vanishes"
		if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
		return(structure(list(modelName = "EII", n=n, d=p, G=G, z=z,
 		                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
	}
	temp <- .Fortran("eseii",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[9:10]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "sigma-squared falls below threshold"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EII", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"hcEII" <-
function(data, partition, minclus = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(minclus < 1) stop("minclus must be positive")
	if(any(is.na(data)))
		stop("missing values not allowed in data")
	#====================================================================
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	dimnames(data) <- NULL
	n <- nrow(data)
	p <- ncol(data)
	if(missing(partition))
		partition <- 1:n
	else if(length(partition) != n)
		stop("partition must assign a class to each observation")
	partition <- partconv(partition, consec = TRUE)
	l <- length(unique(partition))
	attr(partition, "unique") <- l
	m <- l - minclus
	if(m <= 0)
		stop("initial number of clusters is not greater than minclus")
	if(n <= p)
		warning("# of observations <= data dimension")
	#=============================================================
	storage.mode(data) <- "double"
	ld <- max(c((l * (l - 1))/2, 3 * m))
	temp <- .Fortran("hceii",
		data,
		as.integer(n),
		as.integer(p),
		as.integer(partition),
		as.integer(l),
		as.integer(m),
		double(p),
		as.integer(ld),
		double(ld),
                PACKAGE = "mclust")[c(1, 9)]
	temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
	temp[[2]] <- temp[[2]][1:m]
        change <- temp[[2]]
	structure(t(temp[[1]]), initialPartition = partition, 
                  dimensions = dimdat, modelName = "EII", 
                  call =  match.call())
}

"meEII" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	# number of groups
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
	       WARNING <- "z is missing"
	       if (warn) warning(WARNING)
               variance <- list(modelName = "EII", d = p, G = G, sigmasq = NA)
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("meeii",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(1),
			double(K),
                        PACKAGE = "mclust")[7:13]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
		data, G = G, modelName = "EII"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("meeiip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(1),
			double(K),
                        PACKAGE = "mclust")[c(11:17, 10)]
	}
	mu <- matrix(temp[[5]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	sigmasq <- temp[[6]]
	Sigma <- diag(rep(sigmasq, p))
	pro <- temp[[7]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6) || 
           sigmasq <= max(control$eps,0)) {
		WARNING <- "sigma-squared falls below threshold"
		if (warn) warning(WARNING)
		mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if (warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if (warn) warning(WARNING)
		}
		mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(0, c(p, p, G))
		for(k in 1:G)
			sigma[,  , k] <- Sigma
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL) 
        dimnames(mu) <- list(dimnames(data)[[2]], NULL) 
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL) 
        variance <- list(modelName = "EII", d = p, G = G, sigma = sigma, 
                         Sigma = Sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance = variance)
	structure(list(modelName = "EII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik),
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepEII" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EII", d = p, G = G, sigmasq = NA)
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mseii",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p * G),
			double(1),
			double(G),
                        PACKAGE = "mclust")[6:8]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EII"), 
                        prior[names(prior) !="functionName"]))
		temp <- .Fortran("mseiip",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p * G),
			double(1),
			double(G),
                        PACKAGE = "mclust")[10:12]
	}
	mu <- matrix(temp[[1]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	sigmasq <- temp[[2]]
	pro <- temp[[3]]
	sigma <- array(0, c(p, p, G))
	Sigma <- diag(rep(sigmasq, p))
	for(k in 1:G)
		sigma[,  , k] <- Sigma
	WARNING <- NULL
	if(sigmasq > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]], NULL) 
        dimnames(mu) <- list(dimnames(data)[[2]], NULL) 
        dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL) 
        variance <- list(modelName = "EII", d = p, G = G, sigma = sigma, 
                         Sigma = Sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance = variance)
        structure(list(modelName = "EII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)

}

"simEII" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "EII"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  sigmasq <- parameters$variance$sigmasq
  cholSigma <- diag(rep(sqrt(sigmasq), d))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      cholSigma, MARGIN = 2, STAT = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EII")
}

"cdensEVI" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        mu <- parameters$mean
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "EVI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
	temp <- .Fortran("esevi",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
         }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
	structure(z, logarithm = logarithm, modelName = "EVI",
                  WARNING = WARNING, returnCode = ret)
}

"emEVI" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepEVI(data, parameters = parameters, warn = warn)$z  
	meEVI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepEVI" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
                Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	} 
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "EVI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
     		    is.null(parameters$variance$shape)) 
          stop("variance parameters are missing")
	temp <- .Fortran("esevi",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "EVI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"meEVI" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EVI", d = p, G = G, 
                                 scale = NA, shape = matrix(NA,p,G)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="EVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("meevi",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(1),
			double(p * G),
			double(K),
                        PACKAGE = "mclust")[7:14]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EVI"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("meevip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(1),
			double(p * G),
			double(K),
                        PACKAGE = "mclust")[11:18]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	scale <- temp[[6]]
	shape <- matrix(temp[[7]], p, G)
	dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
	pro <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn) warning(WARNING)
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			if(warn) warning("z column sum fell below threshold")
			WARNING <- "z column sum fell below threshold"
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "EVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepEVI" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "EVI", d = p, G = G, 
                                 scale = NA, shape = matrix(NA,p,G)) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="EVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("msevi",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p * G),
			double(1),
			double(p * G),
			double(G),
                        PACKAGE = "mclust")[6:9]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "EVI"), prior[names(
			prior) != "functionName"]))
		temp <- .Fortran("msevip",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p * G),
			double(1),
			double(p * G),
			double(G),
                        PACKAGE = "mclust")[10:13]
	}
	mu <- matrix(temp[[1]], p, G)
	scale <- temp[[2]]
	shape <- matrix(temp[[3]], p, G)
	dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
	pro <- temp[[4]]
	WARNING <- NULL
	if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!c(
		scale, shape))) {
		WARNING <- "cannot compute M-step"
		if(warn)
			warning(WARNING)
		mu[] <- pro[] <- scale <- shape[] <- NA
		sigma <- array(NA, c(p, p, G))
                ret <- -1 
	}
	else {
		sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
                ret <- 0
	}
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "EVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "EVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)
}

"simEVI" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "EVI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  shape <- as.matrix(parameters$variance$shape)
  if(!all(dim(shape) == dim(mean)))
    stop("shape incompatible with mean")
  sss <- sqrt(parameters$variance$scale * shape)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(sss[, k]), MARGIN = 2, STAT = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "EVI")
}

"cdensV" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- drop(data)
	n <- length(data)
        mu <- drop(parameters$mean)
	G <- length(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")])))
            || any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "V", 
                                 WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
        if(is.null(sigmasq))
                stop("variance parameters are missing")
        if(any(sigmasq < 0))
                stop("sigma-squared is negative")
	if(any(!sigmasq)) {
                WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(names(data), NULL)
                return(structure(z, logarithm = logarithm, modelName = "V", 
                                 WARNING = WARNING, returnCode = 9))
	}
        if (length(sigmasq) == 1) sigmasq <- rep(sigmasq,G)
	temp <- .Fortran("es1v",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(-1),
		as.integer(n),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[8:9]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(names(data),NULL)
	structure(z, logarithm = logarithm, modelName = "V",
                     WARNING = WARNING, returnCode = ret)
}

"emV" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepV(data, parameters = parameters, warn = warn)$z  
	meV(data, z = z, prior = prior, control = control, 
            Vinv = parameters$Vinv, warn = warn)
}

"estepV" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- drop(data)
	n <- length(data)
        pro <- parameters$pro 
        pro <- pro/sum(pro)
	l <- length(pro)
        mu <- drop(parameters$mean)
	G <- length(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv 
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")])))
            || any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
                return(structure(list(modelName = "V", n=n, d=1, G=G, z=z,
                       parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq
	if(is.null(sigmasq))
		stop("variance parameters are missing")
	if(any(sigmasq < 0))
		stop("sigma-squared is negative")
	if(any(!sigmasq)) {
		WARNING <- "sigma-squared vanishes"
		if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(names(data), NULL)
                return(structure(list(modelName = "V", n=n, d=1, G=G, z=z,
                       parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
	}
	temp <- .Fortran("es1v",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(pro),
		as.integer(n),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[8:9]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(names(data),NULL)
        structure(list(modelName = "V", n = n, d = 1, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"hcV" <-
function(data, partition, minclus = 1, alpha = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(minclus < 1) stop("minclus must be positive")
	if(any(is.na(data)))
		stop("missing values not allowed in data")
	#=====================================================================
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- as.vector(data)
	n <- length(data)
	if(missing(partition))
		partition <- 1:n
	else if(length(partition) != n)
		stop("partition must assign a class to each observation")
	partition <- partconv(partition, consec = TRUE)
	l <- length(unique(partition))
	attr(partition, "unique") <- l
	m <- l - minclus
	if(m <= 0)
		stop("initial number of clusters is not greater than minclus")
	storage.mode(data) <- "double"
	alpha <- alpha * (vecnorm(data - mean(data))^2/n)
	alpha <- min(alpha, .Machine$double.eps)
	ld <- max(c((l * (l - 1))/2, 3 * m))
	temp <- .Fortran("hc1v",
		data,
		as.integer(n),
		as.integer(partition),
		as.integer(l),
		as.integer(m),
		as.double(alpha),
		as.integer(ld),
		double(ld),
                PACKAGE = "mclust")[c(1, 3, 8)]
	temp[[1]] <- temp[[1]][1:m]
	temp[[2]] <- temp[[2]][1:m]
	temp[[3]] <- temp[[3]][1:m]
        change <- temp[[3]]
	structure(rbind(temp[[1]], temp[[2]]), 	initialPartition = partition, 
                  dimensions = n, modelName = "V",
		  call = match.call())
}

"meV" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- as.vector(data)
	n <- length(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal length of data")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
	       WARNING <- "z is missing"
	       if (warn) warning(WARNING)
               variance <- list(modelName = "V", d=1, G=G, sigmasq = rep(NA,G))
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="V", prior=prior, n=n, d=1, G=G,
                  z=z, parameters=parameters, control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("me1v",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(G),
			as.double(if(is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(G),
			double(G),
			double(K),
                        PACKAGE = "mclust")[6:12]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "V"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("me1vp",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(G),
			double(G),
			double(K),
                        PACKAGE = "mclust")[c(10:16, 9)]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- temp[[5]]
	names(mu) <- as.character(1:G)
	sigmasq <- temp[[6]]
	pro <- temp[[7]]
	## logpost <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6) || any(sigmasq <= max(
		control$eps, 0))) {
		WARNING <- "sigma-squared falls below threshold"
		if(warn) warning(WARNING)
		mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
		ret <- if(control$equalPro) -2 else -3
	}
	else if(its >= control$itmax[1]) {
		WARNING <- "iteration limit reached"
		if (warn) warning(WARNING)
		its <-  - its
		ret <- 1
	}
	else ret <- 0
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(names(data),NULL)
        variance = list(modelName = "V", d = 1, G = G, 
                        sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "V", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik),
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepV" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one-dimensional")
	data <- as.vector(data)
	n <- length(data)
	##
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
# number of groups 
	G <- dimz[2]
	##
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "V", d=1, G=G, sigmasq=rep(NA,G))
               parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                                  variance=variance)
               return(structure(list(modelName="V", prior=prior, n=n, d=1, G=G,
                                     z=z, parameters=parameters), 
                                WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("ms1v",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(G),
			double(G),
			double(G),
			double(G),
                        PACKAGE = "mclust")[5:7]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "V"), prior[names(prior) !=
			"functionName"]))
		storage.mode(z) <- "double"
		temp <- .Fortran("ms1vp",
			as.double(data),
			z,
			as.integer(n),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(G),
			double(G),
			double(G),
                        PACKAGE = "mclust")[9:11]
	}
	mu <- temp[[1]]
	names(mu) <- as.character(1:G)
	sigmasq <- temp[[2]]
	pro <- temp[[3]]
	WARNING <- NULL
	if(any(sigmasq > signif(.Machine$double.xmax, 6))) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
                mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
                print(G)
                print(sigmasq)
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(names(data),NULL)
        variance = list(modelName = "V", d = 1, G = G, 
                        sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "V", prior = prior, n = n, d = 1, G = G, 
                       z = z, parameters = parameters),
                  WARNING = WARNING, returnCode = ret)
}

"simV" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, 2), modelName = "V"))
  }
  if(!is.null(seed))
    set.seed(seed)
  mu <- parameters$mean
  G <- length(mu)
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- rep(0, n)
  sd <- sqrt(parameters$variance$sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd[k])
  }
  structure(cbind(group = clabels, "1" = x), modelName = "V")
}

"cdensVEI" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VEI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
	temp <- .Fortran("esvei",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],dimnames(mu)[[2]])
	structure(z, logarithm = logarithm, modelName = "VEI",
                  WARNING = WARNING, returnCode = ret)
}

"emVEI" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepVEI(data, parameters = parameters, warn = warn)$z  
	meVEI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepVEI" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
                Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	} 
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VEI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
		    is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
	temp <- .Fortran("esvei",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "VEI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"meVEI" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEI", d = p, G = G, 
                                 scale = rep(NA,G), shape = rep(NA,p)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mevei",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax),
			as.double(control$tol),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(p),
			double(K),
			double(G),
			double(p),
			double(p * G),
                        PACKAGE = "mclust")[7:14]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VEI"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("meveip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax),
			as.double(control$tol),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(p),
			double(K),
			double(G),
			double(p),
			double(p * G),
                        PACKAGE = "mclust")[11:18]
	}
	z <- temp[[1]]
	its <- temp[[2]][1]
	inner <- temp[[2]][2]
	err <- temp[[3]][1]
	inerr <- temp[[3]][2]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	scale <- temp[[6]]
	shape <- temp[[7]]
	dimnames(mu) <- list(NULL, as.character(1:G))
	pro <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		sigma <- array(NA, c(p, p, G))
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		sigma <- array(NA, c(p, p, G))
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(0, c(p, p, G))
		for(k in 1:G)
			sigma[,  , k] <- diag(scale[k] * shape)
		if(inner >= control$itmax[2]) {
			WARNING <- "inner iteration limit reached"
			warning(WARNING)
			inner <-  - inner
			ret <- 2
		}
		else if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
	attr(info, "inner") <- c(iterations = inner, error = inerr)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
	variance <- list(modelName = "VEI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "VEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepVEI" <-
function(data, z, prior = NULL, warn = NULL, control = NULL,...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
              WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEI", d = p, G = G, 
                                 scale = rep(NA,G), shape = rep(NA,p)) 
               parameters <- list(pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VEI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))

	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
        if (is.null(control)) control <- emControl()
	itmax <- if(length(control$itmax) == 1) control$itmax else control$
			itmax[2]
	tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
	if(is.null(prior)) {
		temp <- .Fortran("msvei",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.integer(itmax),
			as.double(tol),
			double(p * G),
			double(G),
			double(p),
			double(G),
			double(G),
			double(p),
			double(p * G),
                        PACKAGE = "mclust")[6:11]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VEI"), prior[names(
			prior) != "functionName"]))
		temp <- .Fortran("msveip",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			as.integer(itmax),
			as.double(tol),
			double(p * G),
			double(G),
			double(p),
			double(G),
			double(G),
			double(p),
			double(p * G),
                        PACKAGE = "mclust")[10:15]
	}
	inner <- temp[[1]]
	inerr <- temp[[2]]
	mu <- matrix(temp[[3]], p, G)
	scale <- temp[[4]]
	shape <- temp[[5]]
	dimnames(mu) <- list(NULL, as.character(1:G))
	pro <- temp[[6]]
	WARNING <- NULL
	if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!
		c(scale, shape))) {
		WARNING <- "cannot compute M-step"
		if(warn)
			warning(WARNING)
		mu[] <- pro[] <- shape <- scale[] <- NA
		sigma <- array(NA, c(p, p, G))
                ret <- -1
	}
	else {
                ret <- 0
		sigma <- array(0, c(p, p, G))
		for(k in 1:G)
			sigma[,  , k] <- diag(scale[k] * shape)
		if(inner >= itmax) {
			WARNING <- "inner iteration limit reached"
			warning(WARNING)
			inner <-  - inner
		}
	}
   	info <- c(iterations = inner, error = inerr)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "VEI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VEI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control),
                  info = info, WARNING = WARNING, returnCode = ret)

}

"simVEI" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VEI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(rtscale[k] * rtshape), MARGIN = 2, STAT = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VEI")
}

"cdensVEV" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VEV", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                   is.null(parameters$variance$shape) ||
                   is.null(parameters$variance$orientation)) 
          stop("variance parameters are missing")
	temp <- .Fortran("esvev",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(parameters$variance$orientation),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(p),
		double(p),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[13:14]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL) 
	structure(z, logarithm = logarithm, modelName = "VEV",
                  WARNING = WARNING, returnCode = ret)
}

"emVEV" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepVEV(data, parameters = parameters, warn = warn)$z  
	meVEV(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepVEV" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
       if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VEV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape) ||
                    is.null(parameters$variance$orientation)) 
           stop("variance parameters are missing")
	temp <- .Fortran("esvev",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(parameters$variance$orientation),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(p),
		double(p),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[13:14]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1  
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL) 
        structure(list(modelName = "VEV", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)

}

"meVEV" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEV", d = p, G = G, 
           scale=rep(NA,G), shape=rep(NA,p), orientation=array(NA,c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mevev",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax),
			as.double(control$tol),
			as.double(control$eps),
			as.integer(lwork),
			double(p * G),
			double(G),
			double(p),
			double(p * p * G),
			double(K),
			double(lwork),
			double(p),
                        PACKAGE = "mclust")[7:16]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VEV"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("mevevp",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax),
			as.double(control$tol),
			as.double(control$eps),
			as.integer(lwork),
			double(p * G),
			double(G),
			double(p),
			double(p * p * G),
			double(K),
			double(lwork),
			double(p),
                        PACKAGE = "mclust")[11:20]
	}
	z <- temp[[1]]
	its <- temp[[2]][1]
	inner <- temp[[2]][2]
	err <- temp[[3]][1]
	inerr <- temp[[3]][2]
	loglik <- temp[[4]]
	lapackSVDinfo <- temp[[5]]
	mu <- matrix(temp[[6]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	scale <- temp[[7]]
	shape <- temp[[8]]
	O <- array(temp[[9]], c(p, p, G))
	pro <- temp[[10]]
	WARNING <- NULL
	if(lapackSVDinfo) {
		if(lapackSVDinfo > 0) {
			WARNING <- "LAPACK DGESVD fails to converge"
		}
		else {
			WARNING <- "input error for LAPACK DGESVD"
		}
		warning(WARNING)
		O[] <- shape[] <- scale[] <- NA
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -9
	}
	else if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		O[] <- shape[] <- scale[] <- NA
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- shapeO(shape, O, transpose = TRUE)
		sigma <- sweep(sigma, MARGIN = 3, STATS = scale, FUN = "*")
		if(inner >= control$itmax[2]) {
			WARNING <- "inner iteration limit reached"
			warning(WARNING)
			inner <-  - inner
			ret <- 2
		}
		else if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- structure(c(iterations = its, error = err), inner = c(
		iterations = inner, error = inerr))
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(O) <- 
           list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
##  Sigma = scale * t(O) %*% diag(shape) %*% O
	variance <- list(modelName = "VEV", d = p, G = G, sigma = sigma, 
                        scale = scale, shape = shape, orientation = O)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance) 
	structure(list(modelName = "VEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepVEV" <-
function(data, z, prior = NULL, warn = NULL, control = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	##
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VEV", d = p, G = G, 
      scale = rep(NA,G), shape = rep(NA,p), orientation = array(NA,c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VEV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))

		WARNING <- "z is missing"
		warning(WARNING)
		return(structure(list(n = n, d = p, G = G, mu = matrix(NA,
			p, G), sigma = array(NA, c(p, p, G)), decomp = list(
			d = p, G = G, scale = rep(NA, G), shape = rep(NA, p),
			orientation = array(NA, c(p, p, G))), pro = rep(NA,
			G), modelName = "VEV", prior = prior), WARNING = 
			WARNING))
	}
	#	shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
	if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop(
			"improper specification of z")
        if (is.null(control)) control <- emControl()
	itmax <- if(length(control$itmax) == 1) control$itmax else control$
			itmax[2]
	tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
	lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
	if(is.null(prior)) {
		temp <- .Fortran("msvev",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(lwork),
			as.integer(lwork),
			as.integer(itmax),
			as.double(tol),
			double(p * G),
			double(G),
			double(p),
			double(p * p * G),
			double(G),
                        PACKAGE = "mclust")[7:14]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VEV"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("msvevp",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			double(lwork),
			as.integer(lwork),
			as.integer(itmax),
			as.double(tol),
			double(p * G),
			double(G),
			double(p),
			double(p * p * G),
			double(G),
                        PACKAGE = "mclust")[11:18]
	}
	lapackSVDinfo <- temp[[1]]
	inner <- temp[[2]]
	inerr <- temp[[3]]
	mu <- matrix(temp[[4]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	scale <- temp[[5]]
	shape <- temp[[6]]
	O <- array(temp[[7]], c(p, p, G))
        pro <- temp[[8]] 
	WARNING <- NULL
	if(lapackSVDinfo) {
		if(lapackSVDinfo > 0) {
			WARNING <- "LAPACK DGESVD fails to converge"
			warning(WARNING)
		}
		else {
			WARNING <- "input error for LAPACK DGESVD"
			warning(WARNING)
		}
		O[] <- shape[] <- scale[] <- NA
		sigma <- array(NA, c(p, p, G))
                ret <- -9
	}
	else if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(
		!c(scale, shape))) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
		mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
		sigma <- array(NA, c(p, p, G))
                ret <- -1
	}
	else {
		sigma <- sweep(shapeO(shape, O, transpose = TRUE), MARGIN = 3,
			STATS = scale, FUN = "*")
		if(inner >= itmax) {
			WARNING <- "inner iteration limit reached"
			warning(WARNING)
			inner <-  - inner
		}
                ret <- 2
	}
	info <- c(iteration = inner, error = inerr)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(O) <-
           list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
        variance <- list(modelName = "VEV", d = p, G = G, sigma = sigma, 
                         scale = scale, shape = shape, orientation = O)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VEV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control),
                  info = info, WARNING = WARNING, returnCode = ret)

}

"simVEV" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VEV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    sss <- rtscale[k] * rtshape
    cholSigma <- parameters$variance$orientation[,  , k] * sss
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      cholSigma, MARGIN = 2, STAT = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VEV")
}

"cdensVII" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VII", 
                                 WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq 
	if(any(sigmasq < 0))
		stop("sigma-squared is negative")
	if(any(!sigmasq)) {
		WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VII", 
                                 WARNING = WARNING, returnCode = 9))
	}
	temp <- .Fortran("esvii",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[9:10]
        loglik <- temp[[1]]
        z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
        if(loglik > signif(.Machine$double.xmax, 6)) {
                WARNING <- "sigma-squared falls below threshold"
                if (warn) warning(WARNING)
                z[] <- NA
                ret <- -1
        }
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0 
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
	structure(z, logarithm = logarithm, modelName = "VII",
                  WARNING = WARNING, returnCode = ret)
}

"emVII" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepVII(data, parameters = parameters, warn = warn)$z  
	meVII(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepVII" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
      	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VII", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        sigmasq <- parameters$variance$sigmasq 
	if(is.null(sigmasq))
		stop("variance parameters are missing")
	if(any(sigmasq < 0))
		stop("sigma-squared is negative")
	if(any(!sigmasq)) {
		WARNING <- "sigma-squared vanishes"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VII", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = -1))
	}
	temp <- .Fortran("esvii",
		as.double(data),
		as.double(mu),
		as.double(sigmasq),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[9:10]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "VII", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"hcVII" <-
function(data, partition, minclus = 1, alpha = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(minclus < 1) stop("minclus must be positive")
	if(any(is.na(data)))
		stop("missing values not allowed in data")
	#=====================================================================
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	dimnames(data) <- NULL
	n <- nrow(data)
	p <- ncol(data)
	if(n <= p)
		warning("# of observations <= data dimension")
	if(missing(partition))
		partition <- 1:n
	else if(length(partition) != n)
		stop("partition must assign a class to each observation")
	partition <- partconv(partition, consec = TRUE)
	l <- length(unique(partition))
	attr(partition, "unique") <- l
	m <- l - minclus
	if(m <= 0)
		stop("initial number of clusters is not greater than minclus")
	storage.mode(data) <- "double"
	ll <- (l * (l - 1))/2
	ld <- max(n, ll, 3 * m)
	alpha <- alpha * traceW(data/sqrt(n * p))
	alpha <- max(alpha, .Machine$double.eps)
	temp <- .Fortran("hcvii",
		data,
		as.integer(n),
		as.integer(p),
		as.integer(partition),
		as.integer(l),
		as.integer(m),
		as.double(alpha),
		double(p),
		as.integer(ld),
		double(ld),
                PACKAGE = "mclust")[c(1, 10)]
	temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
	temp[[2]] <- temp[[2]][1:m]
        change <- temp[[2]]
	structure(t(temp[[1]]), initialPartition = partition, 
                  dimensions = dimdat, modelName = "VII", 
                  call = match.call())
}

"meVII" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data must be in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mevii",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(K),
                        PACKAGE = "mclust")[7:13]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VII"), prior[names(prior) !=
			"functionName"]))
		storage.mode(z) <- "double"
		temp <- .Fortran("meviip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(K),
                        PACKAGE = "mclust")[c(11:17, 10)]
	}
	mu <- matrix(temp[[5]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	sigmasq <- temp[[6]]
	pro <- temp[[7]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6) || 
           any(sigmasq <= max(control$eps, 0))) {
		WARNING <- "sigma-squared falls below threshold"
		if (warn) warning(WARNING)
		mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if (warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if (warn) warning(WARNING)
		}
		mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(0, c(p, p, G))
		for(k in 1:G)
			sigma[,  , k] <- diag(rep(sigmasq[k], p))
		if(its >= control$itmax[1]) {
			warning("iteration limit reached")
			WARNING <- "iteration limit reached"
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
	variance <- list(modelName = "VII", d = p, G = G, 
                         sigma = sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "VII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control, 
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"meVVI" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
	       WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVI", d = p, G = G, 
                                 scale = rep(NA,G), shape = matrix(NA,p,G)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
        }
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mevvi",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(p * G),
			double(K),
                        PACKAGE = "mclust")[7:14]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VVI"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("mevvip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(p * G),
			double(K),
                        PACKAGE = "mclust")[11:18]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	scale <- temp[[6]]
	shape <- matrix(temp[[7]], p, G)
	dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
	pro <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn) warning(WARNING)
		sigma <- array(NA, c(p, p, G))
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if (warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		sigma <- array(NA, c(p, p, G))
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
			FUN = "*"), 2, diag), c(p, p, G))
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
	variance <- list(modelName = "VVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepVII" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal number of observations")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                                WARNING = WARNING, returnCode = 9))
        }
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("msvii",
			as.double(data),
			z,
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p * G),
			double(G),
			double(G),
                        PACKAGE = "mclust")[6:8]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VII"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("msviip",
			as.double(data),
			z,
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p * G),
			double(G),
			double(G),
                        PACKAGE = "mclust")[10:12]
	}
	mu <- matrix(temp[[1]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	sigmasq <- temp[[2]]
	pro <- temp[[3]]
	sigma <- array(0, c(p, p, G))
	for(k in 1:G)
		sigma[,  , k] <- diag(rep(sigmasq[k], p))
	WARNING <- NULL
	if(any(sigmasq > signif(.Machine$double.xmax, 6))) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "VII", d = p, G = G, 
                         sigma = sigma, sigmasq = sigmasq, scale = sigmasq)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VII", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters),
                  WARNING = WARNING, returnCode = ret)
}

"simVII" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d), modelName = "VII"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  sigmasq <- parameters$variance$sigmasq
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(rep(sqrt(sigmasq[k]), d)), MARGIN = 2, STAT = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VII")
}

"cdensVVI" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mu", "variance")]))) ||
            any(is.null(parameters[c("pro", "mu", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(z, logarithm = logarithm, modelName = "VVI", 
                                 WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
                    is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
	temp <- .Fortran("esvvi",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, G)
        WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)  
	structure(z, logarithm = logarithm, modelName = "VVI",
                  WARNING = WARNING, retrinCode = ret)
}

"emVVI" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepVVI(data, parameters = parameters, warn = warn)$z  
	meVVI(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepVVI" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if (is.null(Vinv) || Vinv <= 0) 
                  Vinv <- hypvol(data, reciprocal = TRUE)
	}
        if(any(is.na(unlist(parameters[c("pro", "mu", "variance")]))) ||
            any(is.null(parameters[c("pro", "mu", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VVI", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$scale) ||
		    is.null(parameters$variance$shape)) 
           stop("variance parameters are missing")
	temp <- .Fortran("esvvi",
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$scale),
		as.double(parameters$variance$shape),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[10:11]
	loglik <- temp[[1]]
	z <- matrix(temp[[2]], n, K)
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL) 
        structure(list(modelName = "VVI", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"meVVI" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
	if (!is.null(Vinv)) {
		G <- K - 1
		if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
	       WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVI", d = p, G = G, 
                                 scale = rep(NA,G), shape = matrix(NA,p,G)) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VVI", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
        }
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mevvi",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(p * G),
			double(K),
                        PACKAGE = "mclust")[7:14]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VVI"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("mevvip",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(G),
			double(p * G),
			double(K),
                        PACKAGE = "mclust")[11:18]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	scale <- temp[[6]]
	shape <- matrix(temp[[7]], p, G)
	dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
	pro <- temp[[8]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn) warning(WARNING)
		sigma <- array(NA, c(p, p, G))
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if (warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		sigma <- array(NA, c(p, p, G))
		mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
			FUN = "*"), 2, diag), c(p, p, G))
		if(its >= control$itmax[1]) {
			WARNING <- "iteration limit reached"
			warning(WARNING)
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = err)
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
	variance <- list(modelName = "VVI", d = p, G = G, 
                         sigma = sigma, scale = scale, shape = shape)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance)
	structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepVVI" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters), 
                          WARNING = WARNING, returnCode = 9))

	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("msvvi",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p * G),
			double(G),
			double(p * G),
			double(G),
                        PACKAGE = "mclust")[6:9]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VVI"), prior[names(
			prior) != "functionName"]))
		temp <- .Fortran("msvvip",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p * G),
			double(G),
			double(p * G),
			double(G),
                        PACKAGE = "mclust")[10:13]
	}
	mu <- matrix(temp[[1]], p, G)
	scale <- temp[[2]]
	shape <- matrix(temp[[3]], p, G)
	dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
	pro <- temp[[4]]
	WARNING <- NULL
	if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!
		c(scale, shape))) {
		WARNING <- "cannot compute M-step"
		if(warn)
			warning(WARNING)
		mu[] <- pro[] <- shape <- scale[] <- NA
		sigma <- array(NA, c(p, p, G))
                ret <- -1
	}
	else {
		sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
			FUN = "*"), 2, diag), c(p, p, G))
                ret <- 0
	}
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        variance <- list(modelName = "VVI", d = p, G = G, 
                         sigma = sigma, sigmasq = scale, 
                         scale = scale, shape = shape)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters),
                  WARNING = WARNING, returnCode = ret)
}

"simVVI" <- 
function(parameters, n, seed = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c(
    "mean", "variance")]))) {
    warn <- "parameters are missing"
    warning("parameters are missing")
    return(structure(matrix(NA, n, d + 1), modelName = "VVI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- table(clabels)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(!all(dim(rtshape) == dim(mu)))
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
      diag(rtscale[k] * rtshape[, k]), MARGIN = 2, STAT = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, 1:d)
  structure(cbind(group = clabels, x), modelName = "VVI")
}

"cdensVVV" <-
function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,G)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VVV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholsigma))
          stop("variance parameters are missing")
	temp <- .Fortran("esvvv",
		as.logical(1),
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$cholsigma),
		as.double(-1),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(-1),
		double(p),
		double(1),
		double(n * G),
                PACKAGE = "mclust")[10:12]
	lapackCholInfo <- temp[[1]][1]
	loglik <- temp[[2]]
	z <- matrix(temp[[3]], n, G)
        WARNING <- NULL
	if(lapackCholInfo) {
		if(lapackCholInfo > 0) {
			WARNING <- "sigma is not positive definite"
			if (warn) warning(WARNING)
		}
		else {
			WARNING <- "input error for LAPACK DPOTRF"
			if (warn) warning(WARNING)
		}
		z[] <- NA
                ret <- -9
	}
	else if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- NA
                ret <- -1
	}
        else {
          if (!logarithm) z <- exp(z)
          ret <- 0
        }
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
	structure(z, logarithm = logarithm, modelName = "VVV",
                  WARNING = WARNING, returnCode = ret)
}

"emVVV" <-
function(data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        z <- estepVVV(data, parameters = parameters, warn = warn)$z  
	meVVV(data, z = z, prior = prior, control = control, 
              Vinv = parameters$Vinv, warn = warn)
}

"estepVVV" <-
function(data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	if(is.null(dimdat) || length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
        pro <- parameters$pro
	pro <- pro/sum(pro)
	l <- length(pro)
	mu <- as.matrix(parameters$mean)
	G <- ncol(mu)
	noise <- l == G + 1
	if(!noise) {
		if(l != G)
			stop("pro improperly specified")
		K <- G
		Vinv <- NULL
	}
	else {
		K <- G + 1
                Vinv <- parameters$Vinv
		if(is.null(Vinv) || Vinv <= 0)
			Vinv <- hypvol(data, reciprocal = TRUE)
	}
        if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
            any(is.null(parameters[c("pro", "mean", "variance")]))) {
                WARNING <- "parameters are missing"
                if (warn) warning(WARNING)
                z <- matrix(NA,n,K)
                dimnames(z) <- list(dimnames(data)[[1]], NULL)
                return(structure(list(modelName = "VVV", n=n, d=p, G=G, z=z,
                                      parameters=parameters, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
        }
        if (is.null(parameters$variance$cholsigma))
          stop("variance parameters are missing")
	temp <- .Fortran("esvvv",
        	as.logical(1),
		as.double(data),
		as.double(mu),
		as.double(parameters$variance$cholsigma),
		as.double(pro),
		as.integer(n),
		as.integer(p),
		as.integer(G),
		as.double(if (is.null(Vinv)) -1 else Vinv),
		double(p),
		double(1),
		double(n * K),
                PACKAGE = "mclust")[10:12]
	lapackCholInfo <- temp[[1]][1]
	loglik <- temp[[2]]
	z <- matrix(temp[[3]], n, K)
	WARNING <- NULL
	if(lapackCholInfo) {
		if(lapackCholInfo > 0) {
			WARNING <- "sigma is not positive definite"
			warning(WARNING)
		}
		else {
			WARNING <- "input error for LAPACK DPOTRF"
			warning(WARNING)
		}
		z[] <- loglik <- NA
                ret <- -9
	}
	else if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "cannot compute E-step"
		if (warn) warning(WARNING)
		z[] <- loglik <- NA
                ret <- -1
	}
        else ret <- 0
        dimnames(z) <- list(dimnames(data)[[1]],NULL)
        structure(list(modelName = "VVV", n = n, d = p, G = G, 
                       z = z, parameters = parameters, loglik = loglik),
                   WARNING = WARNING, returnCode = ret)
}

"hcVVV" <-
function(data, partition, minclus = 1, alpha = 1, beta = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(minclus < 1) stop("minclus must be positive")
	if(any(is.na(data)))
		stop("missing values not allowed in data")
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) > 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	dimnames(data) <- NULL
	n <- nrow(data)
	p <- ncol(data)
	if(n <= p)
		warning("# of observations <= data dimension")
	if(missing(partition))
		partition <- 1:n
	else if(length(partition) != n)
		stop("partition must assign a class to each observation")
	partition <- partconv(partition, consec = TRUE)
	l <- length(unique(partition))
	attr(partition, "unique") <- l
	m <- l - minclus
	if(m <= 0)
		stop("initial number of clusters is not greater than minclus")
	storage.mode(data) <- "double"
	ll <- (l * (l - 1))/2
	#	dp <- duplicated(partition)
	#x[c((1:n)[!dp],(1:n)[dp]), ], 
	#as.integer(c(partition[!dp], partition[dp])), 
	ld <- max(n, ll + 1, 3 * m)
	alpha <- alpha * traceW(data/sqrt(n * p))
	alpha <- max(alpha, .Machine$double.eps)
	temp <- .Fortran("hcvvv",
		cbind(data, 0.),
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
		double(ld),
                PACKAGE = "mclust")[c(1, 14)]
	temp[[1]] <- temp[[1]][1:m, 1:2, drop = FALSE]
	temp[[2]] <- temp[[2]][1:m]
        change <- temp[[2]] 
	structure(t(temp[[1]]), initialPartition = partition, 
                  dimensions = dimdat, modelName = "VVV", 
                  call = match.call())
}

"meVVV" <-
function(data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should in the form of a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("data and z should have the same row dimension")
	K <- dimz[2]
        if (!is.null(Vinv))  {
		G <- K - 1
		if(Vinv <= 0)Vinv <- hypvol(data, reciprocal = TRUE)
	}
        else G <- K
	if(all(is.na(z))) {
            WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVV", d = p, G = G, 
                sigma = array(NA, c(p,p,G)), cholsigma = array(NA, c(p,p,G))) 
               parameters <- list(Vinv=Vinv, pro=rep(NA,G), 
                                  mean=matrix(NA,p,G), variance=variance)
               return(structure(list(modelName="VVV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     control=control, loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	storage.mode(z) <- "double"
	if(is.null(prior)) {
		temp <- .Fortran("mevvv",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(p * p * G),
			double(K),
			double(p),
                        PACKAGE = "mclust")[7:13]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VVV"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("mevvvp",
			as.logical(control$equalPro),
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(if (is.null(Vinv)) -1 else Vinv),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			z,
			as.integer(control$itmax[1]),
			as.double(control$tol[1]),
			as.double(control$eps),
			double(p * G),
			double(p * p * G),
			double(K),
			double(p),
                        PACKAGE = "mclust")[c(11:17, 10)]
	}
	z <- temp[[1]]
	its <- temp[[2]]
	err <- temp[[3]]
	loglik <- temp[[4]]
	mu <- matrix(temp[[5]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	cholsigma <- array(temp[[6]], c(p, p, G))
	pro <- temp[[7]]
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn) warning(WARNING)
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- -1
	}
	else if(loglik <  - signif(.Machine$double.xmax, 6)) {
		if(control$equalPro) {
			WARNING <- "z column sum fell below threshold"
			if(warn) warning(WARNING)
		}
		else {
			WARNING <- "mixing proportion fell below threshold"
			if(warn) warning(WARNING)
		}
		mu[] <- pro[] <- z[] <- loglik <- NA
		sigma <- array(NA, c(p, p, G))
		ret <- if(control$equalPro) -2 else -3
	}
	else {
		sigma <- array(apply(cholsigma, 3, unchol, upper = TRUE), 
                               c(p,p,G))
		if(its >= control$itmax[1]) {
			warning("iteration limit reached")
			WARNING <- "iteration limit reached"
			its <-  - its
			ret <- 1
		}
		else ret <- 0
	}
	info <- c(iterations = its, error = abs(err))
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- dimnames(cholsigma) <-
              list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
        variance <- list(modelName = "VVV", d = p, G = G,
                          sigma = sigma, cholsigma = cholsigma)
        parameters <- list(Vinv=Vinv, pro=pro, mean=mu, variance=variance) 
	structure(list(modelName = "VVV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters, control = control,
                       loglik = loglik), 
                  info = info, WARNING = WARNING, returnCode = ret)
}

"mstepVVV" <-
function(data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD || length(dimdat) != 2)
		stop("data should be a matrix or a vector")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	z <- as.matrix(z)
	dimz <- dim(z)
	if(dimz[1] != n)
		stop("row dimension of z should equal data length")
	G <- dimz[2]
	if(all(is.na(z))) {
               WARNING <- "z is missing"
               if (warn) warning(WARNING)
               variance <- list(modelName = "VVV", d = p, G = G, 
                      sigma <- array(NA, c(p,p, G)), 
                      cholsigma = array(NA, c(p,p,G))) 
               parameters <- list(pro=rep(NA,G), mean=matrix(NA,p,G), 
                                  variance=variance)
               return(structure(list(modelName="VVV", prior=prior, n=n, d=p, 
                                     G=G, z=z, parameters=parameters, 
                                     loglik=NA), 
                          WARNING = WARNING, returnCode = 9))
	}
	if(any(is.na(z)) || any(z < 0) || any(z > 1))
		stop("improper specification of z")
	if(is.null(prior)) {
		temp <- .Fortran("msvvv",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			double(p),
			double(p * G),
			double(p * p * G),
			double(G),
                        PACKAGE = "mclust")[7:9]
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = G, modelName = "VVV"), 
                        prior[names(prior) != "functionName"]))
		temp <- .Fortran("msvvvp",
			as.double(data),
			as.double(z),
			as.integer(n),
			as.integer(p),
			as.integer(G),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			double(p),
			double(p * G),
			double(p * p * G),
			double(G),
                        PACKAGE = "mclust")[11:13]
	}
	mu <- matrix(temp[[1]], p, G)
	dimnames(mu) <- list(NULL, as.character(1:G))
	cholsigma <- array(temp[[2]], c(p, p, G))
        pro <- temp[[3]]
	WARNING <- NULL
	if(any(c(mu, cholsigma) > signif(.Machine$double.xmax, 6))) {
		WARNING <- "cannot compute M-step"
		if(warn) warning(WARNING)
		mu[] <- sigma[] <- cholsigma[] <- NA
                ret <- -1
	}
	else {
		sigma <- array(apply(cholsigma, 3, unchol, upper = TRUE), 
                                     c(p,p,G))
                ret <- 0
	}
        dimnames(z) <- list(dimnames(data)[[1]], NULL)
        dimnames(mu) <- list(dimnames(data)[[2]], NULL)
        dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                                NULL)
        dimnames(cholsigma) <- list(dimnames(data)[[2]], 
                                    dimnames(data)[[2]], NULL)
        variance <- list(modelName = "VVV", d = p, G = G, 
                         sigma = sigma, cholsigma= cholsigma)
        parameters <- list(pro=pro, mean=mu, variance=variance)
        structure(list(modelName = "VVV", prior = prior, n = n, d = p, G = G, 
                       z = z, parameters = parameters), 
                  WARNING = WARNING, returnCode = ret)

}

"simVVV" <- 
function(parameters, n, seed = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(!is.null(seed)) set.seed(seed)
	mu <- as.matrix(parameters$mean)
	d <- nrow(mu)
	G <- ncol(mu)
	if(any(is.na(parameters[c("mean", "variance")])) || any(is.null(
		parameters[c("mean", "variance")]))) {
		warn <- "parameters are missing"
		warning("parameters are missing")
		return(structure(matrix(NA, n, d + 1), modelName = "VVV"))
	}
	pro <- parameters$pro
	if(is.null(pro))
		pro <- rep(1/G, G)
	clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
	ctabel <- table(clabels)
	x <- matrix(0, n, d)
	if(is.null(cholsigma <- parameters$variance$cholsigma)) {
		if(is.null(sigma <- parameters$variance$sigma)) {
			stop("variance parameters must inlcude either sigma or cholsigma"
				)
		}
		cholsigma <- apply(sigma, 3, chol)
		for(k in 1:ncol(cholsigma))
			sigma[,  , k] <- cholsigma[, k]
		cholsigma <- sigma
	}
	if(dim(cholsigma)[3] != G)
		stop("variance incompatible with mean")
	for(k in 1:G) {
		m <- ctabel[k]
		x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
			ncol = d) %*% cholsigma[,  , k], MARGIN = 2, STAT = mu[
			, k], FUN = "+")
	}
	dimnames(x) <- list(NULL, 1:d)
	structure(cbind(group = clabels, x), modelName = "VVV")
}

"mvnX" <-
function(data, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(!oneD)
		stop("data must be one dimensional")
	data <- as.vector(data)
	n <- length(data)
	if(is.null(prior)) {
		temp <- .Fortran("mvn1d",
			as.double(data),
			as.integer(n),
			double(1),
			double(1),
			double(1),
                        PACKAGE = "mclust")[3:5]
		logpost <- NULL
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = 1, modelName = "X"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("mvn1p",
			as.double(data),
			as.integer(n),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(1),
			double(1),
			double(1),
                        PACKAGE = "mclust")[c(7:9, 6)]
		logpost <- temp[[4]]
	}
	mu <- temp[[1]]
	sigmasq <- temp[[2]]
	loglik <- temp[[3]]
	ret <- 0
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "sigma-squared vanishes"
		if(warn) warning(WARNING)
		loglik <- NA
		ret <- -1
	}
        variance = list(modelName= "X", d = 1, G = 1, sigmasq = sigmasq)
        parameters <- list(mean = mu, variance = variance)
	structure(list(modelName = "X", prior = prior, n = n, d = 1, G = 1, 
                       parameters = parameters, loglik = loglik),
                  WARNING = WARNING, returnCode = ret) 
}

"mvnXII" <-
function(data, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD)
		stop("for multidimensional data only")
	if(length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	if(is.null(prior)) {
		temp <- .Fortran("mvnxii",
			as.double(data),
			as.integer(n),
			as.integer(p),
			double(p),
			double(1),
			double(1),
                        PACKAGE = "mclust")[4:6]
		logpost <- NULL
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = 1, modelName = "XII"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("mnxiip",
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p),
			double(1),
			double(1),
                        PACKAGE = "mclust")[c(8:10, 7)]
		logpost <- temp[[4]]
	}
	mu <- temp[[1]]
	sigmasq <- temp[[2]]
	loglik <- temp[[3]]
	Sigma <- sigmasq * diag(p)
	ret <- 0
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		loglik <- NA
		ret <- -1
	}
        variance <- list(modelName = "XII", d = p, G = 1,	
                         sigmasq	 = sigmasq, Sigma = Sigma, 
                         sigma = array(Sigma, c(p, p, 1)), scale = sigmasq)
        parameters <- list(mean = matrix(mu, ncol = 1), variance = variance) 
        structure(list(modelName = "XII", prior = prior, n = n, d = p, G = 1, 
                       parameters = parameters, loglik = loglik), 
                  WARNING = WARNING, returnCode = ret) 
}

"mvnXXI" <-
function(data, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD)
		stop("for multidimensional data only")
	if(length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	if(is.null(prior)) {
		temp <- .Fortran("mvnxxi",
			as.double(data),
			as.integer(n),
			as.integer(p),
			double(p),
			double(1),
			double(p),
			double(1),
                        PACKAGE = "mclust")[4:7]
		logpost <- NULL
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = 1, modelName = "XXI"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("mnxxip",
			as.double(data),
			as.integer(n),
			as.integer(p),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(priorParams$scale),
			as.double(priorParams$dof),
			double(p),
			double(1),
			double(p),
			double(1),
                        PACKAGE = "mclust")[c(8:11, 7)]
		logpost <- temp[[5]]
	}
	mu <- temp[[1]]
	scale <- temp[[2]]
	shape <- temp[[3]]
	loglik <- temp[[4]]
	Sigma <- diag(scale * shape)
	ret <- 0
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		loglik <- NA
		ret <- -1
	}
	variance <- list(modelName = "XXI", d = p, G = 1, 
                         Sigma = Sigma, sigma = array(Sigma, c(p, p, 1)),
                         scale = scale, shape = shape)
        parameters <- list(mean = matrix(mu, ncol = 1), variance = variance)
	structure(list(modelName = "XXI", prior = prior, n = n, d = p, G = 1, 
                       parameters = parameters, loglik = loglik),
		  WARNING = WARNING, returnCode = ret) 
}

"mvnXXX" <-
function(data, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(warn)) warn <- .Mclust$warn
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD)
		stop("for multidimensional data only")
	if(length(dimdat) != 2)
		stop("data must be a matrix")
	data <- as.matrix(data)
	n <- nrow(data)
	p <- ncol(data)
	if(is.null(prior)) {
		temp <- .Fortran("mvnxxx",
			as.double(data),
			as.integer(n),
			as.integer(p),
			double(p),
			double(p * p),
			double(1),
                        PACKAGE = "mclust")[c(4:6)]
		logpost <- NULL
	}
	else {
		priorParams <- do.call(prior$functionName, c(list(data = 
			data, G = 1, modelName = "XXX"), prior[names(prior) !=
			"functionName"]))
		temp <- .Fortran("mnxxxp",
			as.double(data),
			as.integer(n),
			as.integer(p),
			double(p),
			as.double(priorParams$shrinkage),
			as.double(priorParams$mean),
			as.double(if(any(priorParams$scale)) chol(priorParams$
					scale) else priorParams$scale),
			as.double(priorParams$dof),
			double(p),
			double(p * p),
			double(1),
                        PACKAGE = "mclust")[c(9:11, 8)]
		logpost <- temp[[4]]
	}
	mu <- temp[[1]]
	cholSigma <- matrix(temp[[2]], p, p)
	Sigma <- unchol(cholSigma, upper = TRUE)
	loglik <- temp[[3]]
## Sigma = t(cholSigma) %*% cholSigma
	ret <- 0
	WARNING <- NULL
	if(loglik > signif(.Machine$double.xmax, 6)) {
		WARNING <- "singular covariance"
		if(warn)
			warning(WARNING)
		loglik <- NA
		ret <- -1
	}
        variance <- list(modelName = "XXX", d = 1, G = 1,
                         Sigma = Sigma, cholSigma = cholSigma, 
                         sigma = array(Sigma, c(p, p, 1))) 
        parameters <- list(mean = matrix(mu, ncol = 1), variance = variance)
	structure(list(modelName = "XXX", prior = prior, n = n, d = p, G = 1,
                       parameters = parameters, loglik = loglik), 
                  WARNING = WARNING, returnCode = ret)
}

"adjustedRandIndex" <-
function(x, y)
{
        x <- as.vector(x)
        y <- as.vector(y)
	xx <- outer(x, x, "==")
	yy <- outer(y, y, "==")
	upper <- row(xx) < col(xx)
	xx <- xx[upper]
	yy <- yy[upper]
	a <- sum(as.numeric(xx & yy))
	b <- sum(as.numeric(xx & !yy))
	c <- sum(as.numeric(!xx & yy))
	d <- sum(as.numeric(!xx & !yy))
	ni <- (b + a)
	nj <- (c + a)
	abcd <- a + b + c + d
	q <- (ni * nj)/abcd
	(a - q)/((ni + nj)/2 - q)
}

"bicEMtrain" <-
function(data, labels, modelNames=NULL) 
{
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
            modelNames <- .Mclust$emModelNames
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

"classError" <-
function(classification, truth)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	q <- function(map, len, x)
	{
		x <- as.character(x)
		map <- lapply(map, as.character)
		y <- sapply(map, function(x)
		x[1])
		best <- y != x
		if(all(len) == 1)
			return(best)
		errmin <- sum(as.numeric(best))
		z <- sapply(map, function(x)
		x[length(x)])
		mask <- len != 1
		counter <- rep(0, length(len))
		k <- sum(as.numeric(mask))
		j <- 0
		while(y != z) {
			i <- k - j
			m <- mask[i]
			counter[m] <- (counter[m] %% len[m]) + 1
			y[x == names(map)[m]] <- map[[m]][counter[m]]
			temp <- y != x
			err <- sum(as.numeric(temp))
			if(err < errmin) {
				errmin <- err
				best <- temp
			}
			j <- (j + 1) %% k
		}
		best
	}
        if (any(isNA <- is.na(classification))) {
          classification <- as.character(classification)
          nachar <- paste(unique(classification[!isNA]),collapse="")
          classification[isNA] <- nachar
        }
	MAP <- mapClass(classification, truth)
	len <- sapply(MAP[[1]], length)
	if(all(len) == 1) {
		CtoT <- unlist(MAP[[1]])
		I <- match(as.character(classification), names(CtoT), nomatch= 0)               
		one <- CtoT[I] != truth
	}
	else {
		one <- q(MAP[[1]], len, truth)
	}
	len <- sapply(MAP[[2]], length)
	if(all(len) == 1) {
		TtoC <- unlist(MAP[[2]])
		I <- match(as.character(truth), names(TtoC), nomatch = 0)
		two <- TtoC[I] != classification
	}
	else {
		two <- q(MAP[[2]], len, classification)
	}
	err <- if(sum(as.numeric(one)) > sum(as.numeric(two)))
		as.vector(one)
	else as.vector(two)
        bad <- seq(along = classification)[err]
        list(misclassified = bad, errorRate = length(bad)/length(truth))
}

"cv1EMtrain" <-
function(data, labels, modelNames=NULL) 
{
    z <- unmap(as.numeric(labels))
    G <- ncol(z)
    dimDataset <- dim(data)
    oneD <- is.null(dimDataset) || length(dimDataset[dimDataset > 1]) == 
        1
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
                  z = z[-i, ], warn = FALSE)
                eStep <- do.call("estep", c(mStep, list(data = data[i], 
                                            warn = FALSE)))
                if (is.null(attr(eStep, "warn"))) {
                  k <- (1:G)[eStep$z == max(eStep$z)]
                  l <- (1:G)[z[i, ] == max(z[i, ])]
                  cv[i, m] <- as.numeric(!any(k == l))
                }
            }
        }
    }
    else {
        if (is.null(modelNames)) 
            modelNames <- .Mclust$emModelNames
        n <- nrow(data)
        cv <- matrix(1, nrow = n, ncol = length(modelNames))
        dimnames(cv) <- list(NULL, modelNames)
        for (m in modelNames) {
            for (i in 1:n) {
                mStep <- mstep(modelName = m, data = data[-i,],
                               z = z[-i, ], warn = FALSE)
                eStep <- do.call("estep", c(mStep, list(data = data[i, 
                                       , drop = FALSE], warn = FALSE)))
                if (is.null(attr(eStep, "warn"))) {
                  k <- (1:G)[eStep$z == max(eStep$z)]
                  l <- (1:G)[z[i, ] == max(z[i, ])]
                  cv[i, m] <- as.numeric(!any(k == l))
                }
            }
        }
    }
    errorRate <- apply(cv, 2, sum)
    errorRate/n
}

"hclass" <-
function(hcPairs, G)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	initial <- attributes(hcPairs)$init
	n <- length(initial)
	k <- length(unique(initial))
	G <- if(missing(G)) k:2 else rev(sort(unique(G)))
	select <- k - G
	if(length(select) == 1 && !select)
		return(matrix(initial, ncol = 1, dimnames = list(NULL, 
			as.character(G))))
	bad <- select < 0 | select >= k
	if(all(bad))
		stop("No classification with the specified number of clusters")
	if(any(bad))
		warning("Some selected classifications are inconsistent\n                          with mclust object"
			)
	L <- length(select)
	cl <- matrix(NA, nrow = n, ncol = L, dimnames = list(NULL, as.character(
		G)))
	if(select[1])
		m <- 1
	else {
		cl[, 1] <- initial
		m <- 2
	}
	for(l in 1:max(select)) {
		ij <- hcPairs[, l]
		i <- min(ij)
		j <- max(ij)
		initial[initial == j] <- i
		if(select[m] == l) {
			cl[, m] <- initial
			m <- m + 1
		}
	}
	apply(cl[, L:1, drop = FALSE], 2, partconv, consec = TRUE)
}

"mapClass" <-
function(a, b)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	l <- length(a)
	x <- y <- rep(NA, l)
	if(l != length(b)) {
		warning("unequal lengths")
		return(x)
	}
	aChar <- as.character(a)
	bChar <- as.character(b)
	Tab <- table(a, b)
	Ua <- dimnames(Tab)[[1]]
	Ub <- dimnames(Tab)[[2]]
	aTOb <- rep(list(Ub), length(Ua))
	names(aTOb) <- Ua
	bTOa <- rep(list(Ua), length(Ub))
	names(bTOa) <- Ub
	# -------------------------------------------------------------
	k <- nrow(Tab)
	Map <- rep(0, k)
	Max <- apply(Tab, 1, max)
	for(i in 1:k) {
		I <- match(Max[i], Tab[i,  ], nomatch = 0)
		aTOb[[i]] <- Ub[I]
	}
	if(is.numeric(b))
		aTOb <- lapply(aTOb, as.numeric)
	k <- ncol(Tab)
	Map <- rep(0, k)
	Max <- apply(Tab, 2, max)
	for(j in (1:k)) {
		J <- match(Max[j], Tab[, j])
		bTOa[[j]] <- Ua[J]
	}
	if(is.numeric(a))
		bTOa <- lapply(bTOa, as.numeric)
	list(aTOb = aTOb, bTOa = bTOa)
}

"mclustDA" <-
function(train, test = NULL, pro = NULL, 
         G = NULL, modelNames = NULL, prior = NULL, 
         control = emControl(), initialization = NULL, warn = FALSE, 
         verbose = FALSE, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        if (is.null(test)) test <- train
	if(verbose) cat("training ...\n")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name("mclustDAtrain")
        mc[[2]] <- train$data
        mc[[3]] <- train$labels
        names(mc)[c(2,3)] <- c("data", "labels")
	trainingModels <- eval(mc, parent.frame())
	if(verbose)
		cat("testing ...\n")
	S <- data.frame(trainClass = as.factor(unique(train$labels)), 
                    mclustModel = as.factor(sapply(trainingModels, function(x)
	x$modelName)), numGroups = sapply(trainingModels, function(x) x$G))
        if (!is.list(test)) test <- list(test)
	testFit <- mclustDAtest(test$data, trainingModels)
	testSumry <- summary(testFit)
	trainFit <- mclustDAtest(train$data, trainingModels)
	trainSumry <- summary(trainFit, pro = pro)
	structure(list(test = list(classification = testSumry$classification,
               uncertainty = 1-apply(testSumry$z,1,max), labels = test$labels),
                train = list(classification = trainSumry$classification,
           uncertainty = 1 - apply(trainSumry$z,1,max), labels = train$labels),
                summary = S), class = "mclustDA")
}

"mclustDAtest" <-
function(data, models)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	densfun <- function(model, data)
	{
		do.call("dens", c(list(data = data), model))
	}
	den <- as.matrix(data.frame(lapply(models, densfun, data = data)))
	dimnames(den) <- list(NULL, names(models))
	structure(den, class = "mclustDAtest")
}

"mclustDAtrain" <-
function(data, labels, G = NULL, modelNames = NULL, prior = NULL, 
    control = emControl(), initialization = NULL, warn = FALSE, 
    verbose = TRUE, ...) 
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name("mclustBIC")
        mc$labels <- mc$verbose <- NULL
  	dimData <- dim(data)
	oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
	if(!oneD && length(dimData) != 2)
		stop("data must be a vector or a matrix")
	if(oneD) {
		data <- as.vector(data)
		n <- length(data)
		p <- 1
		data <- as.matrix(data)
	}
	else {
		data <- as.matrix(data)
		n <- nrow(data)
		p <- ncol(data)
	}
	U <- unique(labels)
	L <- length(U)
	S <- rep(0, L)
	M <- rep("XXX", L)
        if(is.null(G)) {
                G <- 1:9
        }
        else {
                G <- sort(G)
        }
        if(any(G) <= 0)
                stop("G must be positive")
	R <- rep(list(matrix(0, n, max(G) + 2)), L)
	Glabels <- as.character(G)
	for(l in 1:L) {
		I <- labels == U[l]
                mc[[2]] <- data[I,]
		BIC <- eval(mc, parent.frame())
     	        SUMMARY <- summary(BIC, data[I,  ])
		S[l] <- SUMMARY$G
		M[l] <- SUMMARY$modelName
		R[[l]] <- c(SUMMARY, list(observations = (1:n)[I]))
	}
	names(S) <- M
	if(verbose)
		print(S)
	names(R) <- U
	if (FALSE) R <- lapply(R, function(x)
	{
		i <- charmatch("Vinv", names(x), nomatch = 0)
		if(i)
			x[ - i]
		else x
	}
	)
        R$Vinv <- NULL
	structure(R, G = G, modelNames = modelNames, prior = prior, 
                  control = control, initialization = initialization,
                  warn = warn, class = "mclustDAtrain")
}

"plot.mclustDA" <-
function(x, trainData, testData, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if (missing(trainData) || missing(testData))
          stop("data not supplied")
        par(ask = TRUE)
	p <- ncol(as.matrix(trainData))
	if(p > 2) {
                dimens <- c(1,2)
		Data <- rbind(testData, trainData)
		xlim <- range((Data[, dimens])[, 1])
		ylim <- range((Data[, dimens])[, 2])
		cl <- c(rep(1, nrow(testData)), rep(2, nrow(trainData)))
		clPairs(Data, cl = cl, symbols = c(1, 3), ...)
		coordProj(dataset = Data, 
                          classification = cl, what = "classification",
                          identify = FALSE, symbols = c(1,3), 
                          xlim=xlim, ylim=ylim, ...)
 	        title("Training and Test Data")
		coordProj(dataset = trainData, what ="classification",
                    classification = x$train$labels, identify = FALSE, 
                    xlim = xlim, ylim = ylim, ...)
  		title("Training Data: known Classification")
		coordProj(dataset = testData, what = "classification",
			  classification = x$train$classification,
                          identify = FALSE, xlim=xlim, ylim=ylim, ...)
		title("Training Data: mclustDA Classification")
		coordProj(dataset = trainData, 
			  classification = x$train$classification, 
                          truth	 = x$train$labels, what = "errors", 
                          identify = FALSE, xlim = xlim, ylim = ylim, ...)
		title("Training Error")
                if (!is.null(x$test$labels)) {
		   coordProj(dataset = testData, what = "classification",
                             classification = x$test$labels, identify = FALSE, 
                             xlim = xlim, ylim = ylim, ...)
                    title("Test Data: known Classification")
                }
		coordProj(dataset = testData, what = "classification",
			  classification = x$test$classification,
                          identify = FALSE, xlim=xlim, ylim=ylim, ...)
		title("Test Data: mclustDA Classification")
                if (!is.null(x$test$labels)) {
		   coordProj(dataset = testData, what = "errors",
                             classification = x$test$classification, 
                             truth = x$test$labels, identify = FALSE, 
                             xlim = xlim, ylim = ylim, ...)
  		    title("Test Error")
                }
	}
	else if (p == 2) {
                dimens <- c(1,2)
		Data <- rbind(testData, trainData)
		xlim <- range((Data[, dimens])[, 1])
		ylim <- range((Data[, dimens])[, 2])
		cl <- c(rep(1, nrow(testData)), rep(2, nrow(trainData)))
		mclust2Dplot(dataset = Data, classification = cl,
            what = "classification", identify = FALSE, symbols = c(1, 3), ...)
		title("Training and Test Data")
		mclust2Dplot(dataset = trainData, 
                             classification = x$train$labels, identify = FALSE,
                             xlim = xlim, ylim = ylim, ...)
		title("Training Data: known Classification")
		mclust2Dplot(dataset = trainData, 
			     classification = x$train$classification, 
                             identify = FALSE, xlim = xlim, ylim = ylim, ...)
		title("Training Data: mclustDA Classification")
		mclust2Dplot(dataset = trainData, 
			     classification = x$train$classification, 
                             truth = x$train$labels, what = "errors", 
                             identify = FALSE, xlim = xlim, ylim = ylim, ...)
		title("Training Error")
                if (!is.null(x$test$labels)) {
   	          mclust2Dplot(dataset = trainData, 
                               classification = x$train$labels, 
                               identify = FALSE, xlim = xlim, ylim = ylim, ...)
 		  title("Testing Data: known Classification")
                }
		mclust2Dplot(dataset = testData, 
			     classification = x$test$classification, 
                             identify = FALSE, xlim = xlim, ylim = ylim, ...)
		title("Test Data: mclustDA Classification")
                if (!is.null(x$test$labels)) {
		  mclust2Dplot(dataset = testData, 
			       classification = x$test$classification, 
                               truth = x$test$labels, what = "errors", 
                              identify = FALSE, xlim = xlim, ylim = ylim, ...)
		  title("Test Error")
                }
	}
	else {
		##
		## p == 1
		##
		Data <- c(testData, trainData)
		xlim <- range(Data)
		cl <- c(rep(1, length(testData)), rep(2, length(trainData)))
		mclust1Dplot(dataset = Data, classification = cl,
			     identify = FALSE, xlim = xlim, ...)
		title("Training and Test Data")
		mclust1Dplot(dataset = trainData, 
			     classification = x$train$labels, identify = FALSE,
                             xlim = xlim, ...)
		title("Training Data: known Classification")
		mclust1Dplot(dataset = trainData, 
                             classification = x$train$classification, 
                             identify = FALSE, xlim = xlim, ...)
		title("Training Data: mclustDA Classification")
		mclust1Dplot(dataset = trainData, 
			     classification = x$train$classification, 
                             truth = x$train$labels, what = "errors", 
                             identify = FALSE, xlim = xlim, ...)
		title("Training Error")
                if (!is.null(x$test$labels)) {
   		  mclust1Dplot(dataset = testData, 
			     classification = x$test$labels, 
                             what = "classification", 
                             identify = FALSE, xlim = xlim, ...)
		  title("Test Data: known classification")
                }
		mclust1Dplot(dataset = testData, 
                             classification = x$test$classification, 
                             identify = FALSE, xlim = xlim, ...)
		title("Test Data: mclustDA Classification")
                if (!is.null(x$test$labels)) {
   		  mclust1Dplot(dataset = testData, 
			     classification = x$test$classification, 
                             truth = x$test$labels, what = "errors", 
                             identify = FALSE, xlim = xlim, ...)
		  title("Test Error")
                }
	}
	invisible()
}

"print.mclustDA" <-
function (x, ndigits = options()$digits, ...) 
{
    cat("\nModeling Summary:\n")
    print(x$summary)
    cat("\nTest Classification Summary:\n")
    print(table(x$test$classification))
    cat("\nTraining Classification Summary:\n")
    print(table(x$train$classification))
    err <- classError(x$train$classification,x$train$labels)$errorRate
    cat("\nTraining Error:",err,"\n")
    if (!is.null(x$test$labels)) {
      err <- classError(x$test$classification,x$test$labels)$errorRate
      cat("\nTest Error:",err,"\n")
    }
    invisible()
}

"print.mclustDAtrain" <-
function(x, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	oldClass(x) <- attr(x, "control") <- NULL
	NextMethod("print")
	invisible()
}

"qclass" <-
function (x, k) 
{
 q <- quantile(x, seq(from = 0, to = 1, by = 1/k))
 cl <- rep(0, length(x))
 q[1]  <- q[1] - 1
 for (i in 1:k) cl[x > q[i] & x <= q[i+1]] <- i
 cl
}

"[.mclustDAtest" <-
function (x, i, j, drop = FALSE) 
{
    clx <- oldClass(x)
    oldClass(x) <- NULL
    NextMethod("[")
}
"summary.mclustDAtest" <-
function(object, pro=NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	clfun <- function(x)
	{
		cl <- names(x)[x == max(x)]
		if(length(cl) > 1)
			NA
		else cl
	}
	if(!is.null(pro)) {
		if(length(pro) != ncol(object))
			stop("wrong number of prior probabilities")
		if(any(pro) < 0)
			stop("pro must be nonnegative")
                object <- sweep(object, MARGIN = 1, FUN = "/", 
                                STATS = apply( object, 1, max))
		object <- sweep(object, MARGIN = 2, FUN = "*",
                                STATS = pro/sum(pro))
	}
	cl <- apply(object, 1, clfun)
	z <- sweep(object, MARGIN=1, STATS=apply(object, 1, sum), FUN="/")
	attr(z, "class") <- NULL
	list(classification = cl, z = z)
}

"summary.mclustDAtrain" <-
function(object, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	L <- length(object)
	M <- max(unlist(lapply(object, function(y)
	y$n)))
	N <- names(object)
	s <- rep(list(list(modelName = "XXX", classification = rep(0, M))),
		times = L)
	names(s) <- N
	for(l in 1:L) {
		s[[l]]$modelName <- object[[l]]$modelName
		cl <- if(!is.null(object[[l]]$z)) map(object[[l]]$z) else rep(1, object[[l]]$n)
		s[[l]]$classification <- cl
	}
	structure(s, class = "summary.mclustDAtrain")
}

`EMclust` <-
function(data, G = NULL, modelNames = NULL, prior = NULL, control = 
  emControl(), initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),
  Vinv = NULL, warn = FALSE, x = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if (!is.null(x)) {
    if (!missing(prior) || !missing(control) || 
        !missing(initialization) || !missing(Vinv)) 
 stop("only G and modelNames may be specified as arguments when x is supplied")
    prior <- attr(x,"prior") 
    control <- attr(x,"control")
    initialization <- attr(x,"initialization")
    Vinv <- attr(x,"Vinv")
    warn <- attr(x,"warn")
  }
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  if (is.null(x)) {
    if (is.null(modelNames)) {
      if (d == 1) {
        modelNames <- c("E", "V")
      }
      else {
        modelNames <- .Mclust$emModelNames
        if (n <= d) {
          m <- match(c("EEE","EEV","VEV","VVV"),.Mclust$emModelNames,nomatch=0)
          modelNames <- modelNames[-m]
        }
      }
    }
    if (is.null(G)) {
      G <- if (is.null(initialization$noise)) 1:9 else 0:9
    }
    else {
      G <- sort(as.numeric(G))
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if (is.null(G)) G <- Glabels
    if (is.null(modelNames)) modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if (all(Gmatch) && all(Mmatch)) {
      attr( x, "G") <- as.numeric(G)
      attr( x, "modelNames") <- modelNames
      attr( x, "returnCodes") <- 
           attr(x, "returnCodes")[as.character(G),modelNames,drop=FALSE]
      return(x[as.character(G),modelNames,drop=FALSE])
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  if (any(as.logical(as.numeric(G))) < 0) {
    if (is.null(initialization$noise)) {
      stop("G must be positive")
    }
    else {
      stop("G must be nonnegative")
    }
   }
  if (d == 1 && any(nchar(modelNames) > 1)) {
    Emodel <- any(sapply(modelNames, function(x)
    charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(modelNames, function(x)
    charmatch("V", x, nomatch = 0)[1]) == 1)
    modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  l <- length(Gall)
  m <- length(Mall)
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, 
                       dimnames = list(as.character(Gall), as.character(Mall)))
  if (!is.null(x)) {
    BIC[dimnames(x)[[1]],dimnames(x)[[2]]] <- x
    RET[dimnames(x)[[1]],dimnames(x)[[2]]] <- attr(x, "returnCodes")
    BIC <- BIC[as.character(G),modelNames,drop=FALSE]
    RET <- RET[as.character(G),modelNames,drop=FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1", ] == EMPTY]) {
         out <- mvn(modelName = mdl, data = data, prior = prior)
         BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                              n = n, d = d, G = 1, equalPro = FALSE)
         RET["1", mdl] <- attr(out, "returnCode")
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
    #######################################################
    # all data in initial hierarchical clustering phase
    #######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = "VVV", data = data)
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data)
          } 
        }
        else {
          hcPairs <- NULL 
      #   hcPairs <- hc(modelName = "E", data = data)
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs))  clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           z <- unmap(clss[, g]) 
         }
         else {
           z <- unmap(qclass( data, as.numeric(g)))
         }
         for (modelName in modelNames[BIC[g,] == EMPTY]) {
            out <- me(modelName = modelName, data = data, z = z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
        }
       }
    }
    else {
    ######################################################
    # initial hierarchical clustering phase on a subset
    ######################################################
      if (is.logical(initialization$subset)) 
        initialization$subset <- (1:n)[initialization$subset]
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
             hcPairs <- hc(modelName = "VVV", 
                   data = data[initialization$subset,  ])
          }
          else {
             hcPairs <- hc(modelName = "EII", 
                   data = data[initialization$subset,  ])
          }
         }
       else {
          hcPairs <- NULL
     #    hcPairs <- hc(modelName = "E", 
     #                  data = data[initialization$subset])
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           z <- unmap(clss[, g]) 
         }
         else {
           z <- unmap(qclass(data[initialization$subset], as.numeric(g)))
         }
         dimnames(z) <- list(as.character(initialization$subset), NULL)
         for (modelName in modelNames[!is.na(BIC[g,])]) {
            ms <- mstep(modelName = modelName, z = z, 
                        data = as.matrix(data)[initialization$subset,  ],
                        prior = prior, control = control, warn = warn)
#
#  ctrl <- control
#  ctrl$itmax[1] <- 1
#  ms <- me(modelName = modelName, data = as.matrix(data)[
#           initialization$subset,  ], z = z, prior = prior, control = ctrl)
#
            es <- do.call("estep", c(list(data = data, warn = warn), ms))
            out <- me(modelName = modelName, data = data, z = es$z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
         }
       }
     }
   }
 else {
    ######################################################
    # noise case
    ######################################################
    if (!is.null(initialization$subset)) 
      stop("subset option not implemented with noise")
    if (is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
    noise <- initialization$noise
    if (!is.logical(noise))
      noise <- as.logical(match(1:n, noise, nomatch = 0))
    if (!G[1]) {
      hood <- n * logb(Vinv)
      BIC["0",  ] <- 2 * hood - logb(n)
      if (l == 1) {
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, 
     initialization = list(hcPairs = hcPairs, subset = initialization$subset), 
                       warn = warn, n = n, d = d, oneD = oneD,
                       returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$hcPairs)) {
      if (d != 1) {
        if (n > d) {
          hcPairs <- hc(modelName = "VVV", data = data[!noise,  ])
        }
        else {
          hcPairs <- hc(modelName = "EII", data = data[!noise,  ])
        }
      }
      else {
        hcPairs <- NULL 
   #    hcPairs <- hc(modelName = "E", data = data[!noise])
      }
    }
    else hcPairs <- initialization$hcPairs
    if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
    z <- matrix(0, n, max(G) + 1)
    for (g in Glabels) {
       z[] <- 0
       k <- as.numeric(g)
       if (d > 1 || !is.null(hcPairs)) {
         z[!noise, 1:k] <- unmap(clss[, g])
       }
       else {
         z[!noise, 1:k] <- unmap(qclass(data[!noise]))
       }
       z[noise, k+1] <- 1
       K <- 1:(k+1) 
       for (modelName in modelNames[BIC[g,] == EMPTY]) {
         out <- me(modelName = modelName, data = data, z = z[, K], 
                   prior = prior, control = control, Vinv = Vinv, warn = warn)
         BIC[g, modelName] <- bic(modelName = modelName, loglik = out$loglik, 
                                  n = n, d = d, G = k, 
                                  noise = TRUE, 
                                  equalPro = control$equalPro)
         RET[g, modelName] <- attr(out, "returnCode")
       }
     }
  }
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            control = control, 
            initialization = list(hcPairs = hcPairs, 
                                  subset = initialization$subset,
                                  noise = initialization$noise), 
            Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD,
            returnCodes = RET, class = "mclustBIC")
}

"Mclust" <-
function(data, G = NULL, modelNames = NULL, prior = NULL, 
         control = emControl(), initialization = NULL, warn = FALSE, ...) 
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        mc <- match.call(expand.dots = FALSE)
        mc[[1]] <- as.name("mclustBIC")
  	Bic <- eval(mc, parent.frame())
        G <- attr(Bic, "G")
        modelNames <- attr(Bic,"modelNames") 
  	Sumry <- summary(Bic, data, G=G, modelNames = modelNames) 
	if(!(length(G) == 1)) {
		bestG <- length(unique(Sumry$cl))
		if (bestG == max(G))
   		  warning("optimal number of clusters occurs at max choice")
		else if (bestG == min(G))
  	          warning("optimal number of clusters occurs at min choice")
	}
        attr(Bic,"n") <- attr(Bic, "warn") <- NULL
        attr(Bic,"initialization") <- attr(Bic,"control") <- NULL
        attr(Bic,"d") <- attr(Bic,"returnCodes") <- attr(Bic,"class") <- NULL
	oldClass(Sumry) <- NULL
	Sumry$bic <- Sumry$bic[1]
 	ans <- c(list(BIC=Bic),Sumry)
        orderedNames <- c("modelName", "n", "d", "G", "BIC", "bic", "loglik",
                          "parameters", "z", "classification", "uncertainty") 
	structure(ans[orderedNames], class = "Mclust")
}

`mclustBIC` <-
function(data, G = NULL, modelNames = NULL, prior = NULL, control = 
  emControl(), initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),
  Vinv = NULL, warn = FALSE, x = NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if (!is.null(x)) {
    if (!missing(prior) || !missing(control) || 
        !missing(initialization) || !missing(Vinv)) 
 stop("only G and modelNames may be specified as arguments when x is supplied")
    prior <- attr(x,"prior") 
    control <- attr(x,"control")
    initialization <- attr(x,"initialization")
    Vinv <- attr(x,"Vinv")
    warn <- attr(x,"warn")
  }
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  if (is.null(x)) {
    if (is.null(modelNames)) {
      if (d == 1) {
        modelNames <- c("E", "V")
      }
      else {
        modelNames <- .Mclust$emModelNames
        if (n <= d) {
          m <- match(c("EEE","EEV","VEV","VVV"),.Mclust$emModelNames,nomatch=0)
          modelNames <- modelNames[-m]
        }
      }
    }
    if (is.null(G)) {
      G <- if (is.null(initialization$noise)) 1:9 else 0:9
    }
    else {
      G <- sort(as.numeric(G))
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if (is.null(G)) G <- Glabels
    if (is.null(modelNames)) modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if (all(Gmatch) && all(Mmatch)) {
      attr( x, "G") <- as.numeric(G)
      attr( x, "modelNames") <- modelNames
      attr( x, "returnCodes") <- 
           attr(x, "returnCodes")[as.character(G),modelNames,drop=FALSE]
      return(x[as.character(G),modelNames,drop=FALSE])
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  if (any(as.logical(as.numeric(G))) < 0) {
    if (is.null(initialization$noise)) {
      stop("G must be positive")
    }
    else {
      stop("G must be nonnegative")
    }
   }
  if (d == 1 && any(nchar(modelNames) > 1)) {
    Emodel <- any(sapply(modelNames, function(x)
    charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(modelNames, function(x)
    charmatch("V", x, nomatch = 0)[1]) == 1)
    modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  l <- length(Gall)
  m <- length(Mall)
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, 
                       dimnames = list(as.character(Gall), as.character(Mall)))
  if (!is.null(x)) {
    BIC[dimnames(x)[[1]],dimnames(x)[[2]]] <- x
    RET[dimnames(x)[[1]],dimnames(x)[[2]]] <- attr(x, "returnCodes")
    BIC <- BIC[as.character(G),modelNames,drop=FALSE]
    RET <- RET[as.character(G),modelNames,drop=FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1", ] == EMPTY]) {
         out <- mvn(modelName = mdl, data = data, prior = prior)
         BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                              n = n, d = d, G = 1, equalPro = FALSE)
         RET["1", mdl] <- attr(out, "returnCode")
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
    #######################################################
    # all data in initial hierarchical clustering phase
    #######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = "VVV", data = data)
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data)
          } 
        }
        else {
          hcPairs <- NULL 
      #   hcPairs <- hc(modelName = "E", data = data)
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs))  clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           z <- unmap(clss[, g]) 
         }
         else {
           z <- unmap(qclass( data, as.numeric(g)))
         }
         for (modelName in modelNames[BIC[g,] == EMPTY]) {
            out <- me(modelName = modelName, data = data, z = z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
        }
       }
    }
    else {
    ######################################################
    # initial hierarchical clustering phase on a subset
    ######################################################
      if (is.logical(initialization$subset)) 
        initialization$subset <- (1:n)[initialization$subset]
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
             hcPairs <- hc(modelName = "VVV", 
                   data = data[initialization$subset,  ])
          }
          else {
             hcPairs <- hc(modelName = "EII", 
                   data = data[initialization$subset,  ])
          }
         }
       else {
          hcPairs <- NULL
     #    hcPairs <- hc(modelName = "E", 
     #                  data = data[initialization$subset])
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
      for (g in Glabels) {
         if (d > 1 || !is.null(hcPairs)) {
           z <- unmap(clss[, g]) 
         }
         else {
           z <- unmap(qclass(data[initialization$subset], as.numeric(g)))
         }
         dimnames(z) <- list(as.character(initialization$subset), NULL)
         for (modelName in modelNames[!is.na(BIC[g,])]) {
            ms <- mstep(modelName = modelName, z = z, 
                        data = as.matrix(data)[initialization$subset,  ],
                        prior = prior, control = control, warn = warn)
#
#  ctrl <- control
#  ctrl$itmax[1] <- 1
#  ms <- me(modelName = modelName, data = as.matrix(data)[
#           initialization$subset,  ], z = z, prior = prior, control = ctrl)
#
            es <- do.call("estep", c(list(data = data, warn = warn), ms))
            out <- me(modelName = modelName, data = data, z = es$z, 
                      prior = prior, control = control, warn = warn)
            BIC[g, modelName] <- bic(modelName = modelName, 
                                     loglik = out$loglik,
                                     n = n, d = d, G = as.numeric(g), 
                                     equalPro = control$equalPro)
            RET[g, modelName] <- attr(out, "returnCode")
         }
       }
     }
   }
 else {
    ######################################################
    # noise case
    ######################################################
    if (!is.null(initialization$subset)) 
      stop("subset option not implemented with noise")
    if (is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
    noise <- initialization$noise
    if (!is.logical(noise))
      noise <- as.logical(match(1:n, noise, nomatch = 0))
    if (!G[1]) {
      hood <- n * logb(Vinv)
      BIC["0",  ] <- 2 * hood - logb(n)
      if (l == 1) {
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, 
     initialization = list(hcPairs = hcPairs, subset = initialization$subset), 
                       warn = warn, n = n, d = d, oneD = oneD,
                       returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$hcPairs)) {
      if (d != 1) {
        if (n > d) {
          hcPairs <- hc(modelName = "VVV", data = data[!noise,  ])
        }
        else {
          hcPairs <- hc(modelName = "EII", data = data[!noise,  ])
        }
      }
      else {
        hcPairs <- NULL 
   #    hcPairs <- hc(modelName = "E", data = data[!noise])
      }
    }
    else hcPairs <- initialization$hcPairs
    if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
    z <- matrix(0, n, max(G) + 1)
    for (g in Glabels) {
       z[] <- 0
       k <- as.numeric(g)
       if (d > 1 || !is.null(hcPairs)) {
         z[!noise, 1:k] <- unmap(clss[, g])
       }
       else {
         z[!noise, 1:k] <- unmap(qclass(data[!noise]))
       }
       z[noise, k+1] <- 1
       K <- 1:(k+1) 
       for (modelName in modelNames[BIC[g,] == EMPTY]) {
         out <- me(modelName = modelName, data = data, z = z[, K], 
                   prior = prior, control = control, Vinv = Vinv, warn = warn)
         BIC[g, modelName] <- bic(modelName = modelName, loglik = out$loglik, 
                                  n = n, d = d, G = k, 
                                  noise = TRUE, 
                                  equalPro = control$equalPro)
         RET[g, modelName] <- attr(out, "returnCode")
       }
     }
  }
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            control = control, 
            initialization = list(hcPairs = hcPairs, 
                                  subset = initialization$subset,
                                  noise = initialization$noise), 
            Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD,
            returnCodes = RET, class = "mclustBIC")
}

"mclustModel" <-
function(data, BICvalues, G=NULL, modelNames=NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        mc <- match.call(expand.dots = FALSE)
        if (is.null(attr(BICvalues,"initialization")$noise)) {
          mc[[1]] <- as.name("summaryMclustBIC")
        }
        else {
          mc[[1]] <- as.name("summaryMclustBICn")
        }
        nm <- names(mc)
        mc[1:3] <- mc[c(1,3,2)]
        nm[1:3] <- nm[c(1,3,2)]
        nm[nm == "BICvalues"] <- "object" 
        names(mc) <- nm
        ans <- eval(mc, parent.frame())
        ans$classification <- ans$uncertainty <- NULL
        attr( ans, "bestBICvalues") <- NULL
        attr( ans, "prior") <- NULL
        attr( ans, "control") <- NULL
        attr( ans, "initialization") <- NULL
        oldClass(ans) <- "mclustModel"
        ans
}

"pickBIC" <-
function(x, k = 3)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  Glabels <- dimnames(x)[[1]]
  modelNames <- dimnames(x)[[2]]
  mis <- is.na(x)
  if (all(mis)) {
    warning("none of the selected models could be fitted")
    return(rep(NA,k))
  }
  x[mis] <-  - .Machine$double.xmax
  x <- data.frame(as.vector(x), Glabels[as.vector(row(x))], 
                                modelNames[as.vector(col(x))])
  x <- x[rev(order(x[,1])),]
  namesx <- apply(x[,-1,drop = FALSE], 1, function(z) 
                  paste(as.character(z[2]), as.character(z[1]), sep = ","))
  k <- min(k, nrow(x))
  x <- x[1:k,1]
  x[x ==  - .Machine$double.xmax] <- NA
  namesx <- namesx[1:k]
  namesx[is.na(x)] <- " "
  names(x) <- namesx
  x
}

"print.Mclust" <-
function(x, ndigits = options()$digits, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	M <- switch(x$model,
		X = "univariate normal",
		E = "equal variance",
		V = "unequal variance",
		XII = "spherical multivariate normal",
		EII = "spherical, equal volume",
		VII = "spherical, varying volume",
		XXI = "diagonal multivariate normal",
		EEI = "diagonal, equal volume and shape",
		VEI = "diagonal, equal shape",
		EVI = "diagonal, equal volume",
		VVI = "diagonal, varying volume and shape",
		EEE = "elliposidal, equal variance",
		VEE = "elliposidal, equal shape and orientation",
		EVE = "elliposidal, equal volume and orientation",
		VVE = "ellipsoidal, equal orientation",
		EEV = "ellipsoidal, equal volume and shape",
		VEV = "ellipsoidal, equal shape",
		EVV = "elliposidal, equal volume",
		XXX = "elliposidal multivariate normal",
		VVV = "ellipsoidal, unconstrained",
		stop("invalid model id for EM"))
	G <- length(unique(x$classification))
	cat("\n best model:", M, "with", G, "components\n")
#	aveUncer <- signif(mean(x$uncertainty), 3)
#	medUncer <- signif(median(x$uncertainty), 3)
#	cat("\n average/median classification uncertainty:", aveUncer, "/",
#		medUncer, "\n\n")
	invisible()
}

"print.mclustBIC" <-
function(x, fill = FALSE, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	subset <- !is.null(attr(x, "subset"))
	oldClass(x) <- attr(x, "args") <- NULL
	attr(x, "control") <- attr(x, "initialization") <- NULL
	attr(x, "oneD") <- attr(x, "warn") <- attr(x, "Vinv") <- NULL
	attr(x, "prior") <- attr(x, "G") <- attr(x, "modelNames") <- NULL
	ret <- attr(x, "returnCodes") == -3
	n <- attr(x, "n")
	d <- attr(x, "d")
	attr(x, "returnCodes") <- attr(x, "n") <- attr(x, "d") <- NULL
        ##
	## if(!subset && any(ret) && fill) {
	##	x <- bicFill(x, ret, n, d)
	## }
        ##
	cat("\n BIC:\n")
	NextMethod("print")
	cat("\n")
	invisible()
}

"print.summary.mclustBIC" <-
function(x, ndigits = options()$digits, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	cat("\nclassification table:")
	print(table(x$classification), ...)
	#----------------------------------------------------------------------
	#	cat("\nuncertainty (quartiles):\n")
	#	print(quantile(x$uncertainty), digits = ndigits, ...)
	#----------------------------------------------------------------------
	bic <- attr(x,"bestBICvalues")
	l <- length(bic)
	if(l == 1) {
          cat("\nBIC value:\n")
  	  print(round(bic, ndigits))
        }
        else {
          cat("\nbest BIC values:\n")
  	  print(round(bic, ndigits))
        }
	M <- switch(x$model,
		X = "univariate normal",
		E = "equal variance",
		V = "unequal variance",
		XII = "spherical multivariate normal",
		EII = "spherical, equal volume",
		VII = "spherical, varying volume",
		XXI = "diagonal multivariate normal",
		EEI = "diagonal, equal volume and shape",
		VEI = "diagonal, equal shape",
		EVI = "diagonal, equal volume",
		VVI = "diagonal, varying volume and shape",
		EEE = "ellipsoidal, equal variance",
		VEE = "ellipsoidal, equal shape and orientation",
		EVE = "ellipsoidal, equal volume and orientation",
		VVE = "ellipsoidal, equal orientation",
		EEV = "ellipsoidal, equal volume and shape",
		VEV = "ellipsoidal, equal shape",
		EVV = "ellipsoidal, equal volume",
		XXX = "ellipsoidal multivariate normal",
		VVV = "ellipsoidal, unconstrained",
		stop("invalid model id for EM"))
##	cat("\nbest model:", M, "\n\n")
	##
	##	print(x$options)
	invisible()
}

"summary.mclustBIC" <-
function(object, dataset, G, modelNames, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        mc <- match.call(expand.dots = FALSE)
        if (is.null(attr(object,"initialization")$noise)) {
          mc[[1]] <- as.name("summaryMclustBIC")
        }
        else {
          mc[[1]] <- as.name("summaryMclustBICn")
        }
       ans <- eval(mc, parent.frame())
       Glabels <- dimnames(object)[[1]]
       if (length(Glabels) != 1 && (!missing(G) && length(G) > 1)) {
         Grange <- range(as.numeric(Glabels))
         if (match(ans$G, Grange, nomatch = 0))
  warning("best model occurs at the min or max # of components considered")
       }
       ans
}

"summary.mclustModel" <-
function(object, ...) 
{
  classification <- map(object$z)
  G <- dim(object$parameters$mean)[2]
  G1 <- dim(object$z)[2]
  if (G1 == (G + 1)) classification[classification == G1] <- 0
  uncertainty <- 1 - apply(object$z, 1, max)
  names(classification) <- names(uncertainty) <- dimnames(object$z)[[1]]
  data.frame(classification = classification, uncertainty = uncertainty)
}

"summaryMclustBIC" <-
function(object, data, G=NULL, modelNames=NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  initialization <- attr(object, "initialization")
  hcPairs <- initialization$hcPairs
  subset <- initialization$subset
  prior <- attr(object, "prior")
  control <- attr(object, "control")
  warn <- attr(object, "warn")
  oldClass(object) <- NULL
  attr(object, "prior") <- attr(object, "warn") <- NULL
  attr(object, "modelNames") <- attr(object, "oneD") <- NULL
  attr(object, "initialization") <- attr(object, "control") <- NULL
  d <- if(is.null(dim(data))) 1 else ncol(data)
  ##
  if (is.null(G))
    G <- dimnames(object)[[1]]
  if (is.null(modelNames))
    modelNames <- dimnames(object)[[2]]
  bestBICs <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
  if (all(is.na(bestBICs))) {
    return(structure(NULL, bestBICvalues = bestBICs, prior = prior, control
      = control, initialization = initialization, class = "summary.mclustBIC"))
  }
  temp <- unlist(strsplit(names(bestBICs)[1], ","))
  bestModel <- temp[1]
  G <- as.numeric(temp[2])
  if(G == 1) {
    out <- mvn(modelName = bestModel, data = data, prior = prior)
    ans <- c(list(bic = bestBICs, classification = rep(1, n), 
                  uncertainty = rep(0, n)), out)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", 
                      "parameters", "z", "classification", "uncertainty")
    return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                     prior = prior,  control = control, 
                     initialization = initialization, 
                     class = "summary.mclustBIC"))
  }
  if(is.null(subset)) {
    if(d > 1 || !is.null(hcPairs)) {
      z <- unmap(hclass(hcPairs, G))
    }
    else {
      z <- unmap(qclass(data, G))
    }
    out <- me(modelName = bestModel, data = data, z = z, prior = prior, 
      control = control, warn = warn)
  }
  else {
    if(d > 1 || !is.null(hcPairs)) {
      z <- unmap(hclass(hcPairs, G))
    }
    else {
      z <- unmap(qclass(data[subset], G))
    }
    ms <- mstep(modelName = bestModel, prior = prior, z = z, data = 
      as.matrix(data)[subset,  ], control = control, warn = warn)
    es <- do.call("estep", c(list(data = data), ms))
    out <- me(modelName = bestModel, data = data, z = es$z, 
              prior = prior, control = control, warn = warn)
  }
  obsNames <- if (is.null(dim(data))) {
                names(data) 
              }
              else {
                dimnames(data)[[1]]
              }
  classification <- map(out$z)
  uncertainty <- 1 - apply(out$z, 1, max)
  names(classification) <- names(uncertainty) <- obsNames
  ans <- c(list(bic = as.vector(bestBICs[1]), classification = classification, 
    uncertainty = uncertainty), out)
  orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", 
    "z", "classification", "uncertainty")
  structure(ans[orderedNames], bestBICvalues = bestBICs, prior = prior, control
     = control, initialization = initialization, class = "summary.mclustBIC")
}

"summaryMclustBICn" <-
function(object, data, G=NULL, modelNames=NULL, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  initialization <- attr(object, "initialization")
  hcPairs <- initialization$hcPairs
  noise <- initialization$noise
  if(!is.logical(noise))
    noise <- as.logical(match(1:n, noise, nomatch = 0))
  prior <- attr(object, "prior")
  control <- attr(object, "control")
  warn <- attr(object, "warn")
  Vinv <- attr(object, "Vinv")
  oldClass(object) <- NULL
  attr(object, "control") <- attr(object, "initialization") <- NULL
  attr(object, "prior") <- attr(object, "Vinv") <- NULL
  attr(object, "warn") <- NULL
  ##
  if (is.null(G)) G <- dimnames(object)[[1]]
  if (is.null(modelNames)) modelNames <- dimnames(object)[[2]]
  bestBICs <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
  if (all(is.na(bestBICs))) {
    return(structure(NULL, bestBICvalues = bestBICs, prior = prior, control
      = control, initialization = initialization, class = "summary.mclustBIC"))
  }
  temp <- unlist(strsplit(names(bestBICs)[1], ","))
  bestModel <- temp[1] 
  G <- as.numeric(temp[2])
  if(G == 0) {
    ans <- list(bic = bestBICs, classification = rep(1, n), 
                uncertainty = rep(0, n), n = n, d = ncol(data), 
                G = 0, loglik = n * logb(Vinv), Vinv = Vinv)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "Vinv", 
                      "classification", "uncertainty")
    return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                     prior = prior, control = control, 
                     initialization = initialization, 
                     class = "summary.mclustBIC"))
  }
  G1 <- G + 1
  z <- matrix(0, n, G1)
  if(d > 1 || !is.null(hcPairs)) {
    z[!noise, 1:G] <- unmap(hclass(hcPairs, G))
  }
  else {
    z[!noise, 1:G] <- unmap(qclass(data, G))
  }
  z[noise, G1] <- 1
  out <- me(modelName = bestModel, data = data, z = z, prior = prior, 
    control = control, warn = warn, Vinv = Vinv)
  obsNames <- if (is.null(dim(data))) {
                names(data)
              }  
              else {
                dimnames(data)[[1]]
              }
  classification <- map(out$z)
  classification[classification == G1] <- 0
  uncertainty <- 1 - apply(out$z, 1, max)
  names(classification) <- names(uncertainty) <- obsNames
  ans <- c(list(bic = as.vector(bestBICs[1]), classification = classification, 
    uncertainty = uncertainty, Vinv = Vinv), out)
  orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", 
    "z", "Vinv", "classification", "uncertainty")
  structure(ans[orderedNames], bestBICvalues = bestBICs, prior = prior, control
     = control, initialization = initialization, class = "summary.mclustBIC")
}

"defaultPrior" <-
function(data, G, modelName, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	aux <- list(...)
	if(is.null(aux$shrinkage)) {
		shrinkage <- 0.01
	}
	else if(is.na(aux$shrinkage) || !aux$shrinkage) {
		shrinkage <- 0
	}
	else if(aux$shrinkage < 0) {
		stop("negative value given for shrinkage")
	}
	else {
		shrinkage <- aux$shrinkage
	}
	if(is.null(aux$mean)) {
		mean <- if (is.null(dim(data))) 
                          mean(data) else colMeans(data)
	}
	else if(any(is.na(aux$mean))) {
		if(shrinkage)
			stop("positive shrinkage with no prior mean specified")
		mean <- if (is.null(dim(data))) 
                          mean(data) else colMeans(data)
	}
	else {
		if(!shrinkage)
			stop("prior mean specified but not shrinkage")
		mean <- aux$mean
	}
	switch(EXPR = modelName,
		E = ,
		V = ,
		X = {
			dof <- 3
			if(is.null(aux$scale)) {
				scale <- var(data)/G^2
			}
			else {
				scale <- aux$scale
			}
			list(shrinkage = shrinkage, mean = mean, dof = dof,
				scale = scale)
		}
		,
		EII = ,
		VII = ,
		XII = ,
		EEI = ,
		EVI = ,
		VEI = ,
		VVI = ,
		XXI = {
			n <- nrow(data)
			p <- ncol(data)
			dof <- p + 2
			if(is.null(aux$scale)) {
				fac <- (1/G)^(2/p)
				scale <- (fac * sum(apply(data, 2, var)))/
					p
			}
			else {
				scale <- aux$scale
			}
			list(shrinkage = shrinkage, mean = mean, dof = dof,
				scale = scale)
		}
		,
		EEE = ,
		EEV = ,
		VEV = ,
		VVV = ,
		XXX = {
			n <- nrow(data)
			p <- ncol(data)
			dof <- p + 2
			if(is.null(aux$scale)) {
				fac <- (1/G)^(2/p)
				if(n > p) {
					scale <- fac * var(data)
				}
				else {
					scale <- fac * diag(apply(data,
						2, var))
				}
			}
			else {
				scale <- aux$scale
			}
			list(shrinkage = shrinkage, mean = mean, dof = dof,
				scale = scale)
		}
		,
		stop("no default prior for this model"))
}

".Mclust" <-
structure(list(emModelNames = c("EII", "VII", "EEI", "VEI", "EVI", 
"VVI", "EEE", "EEV", "VEV", "VVV"), hcModelNames = c("EII", "VII", 
"EEE", "VVV"), bicPlotSymbols = structure(c(17, 2, 16, 10, 13, 
1, 15, 12, 7, 0, 17, 2), .Names = c("EII", "VII", "EEI", "EVI", 
"VEI", "VVI", "EEE", "EEV", "VEV", "VVV", "E", "V")), bicPlotColors = structure(c("gray", 
"black", "orange", "brown", "red", "magenta", "forestgreen", 
"green", "cyan", "blue", "gray", "black"), .Names = c("EII", 
"VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV", "VVV", 
"E", "V")), classPlotSymbols = c(17, 0, 10, 4, 11, 18, 6, 7, 
3, 16, 2, 12, 8, 15, 1, 9, 14, 13, 5), classPlotColors = c("blue", 
"red", "green", "cyan", "magenta", "forestgreen", "purple", "orange", 
"gray", "brown", "black"), warn = TRUE), .Names = c("emModelNames", 
"hcModelNames", "bicPlotSymbols", "bicPlotColors", "classPlotSymbols", 
"classPlotColors", "warn"))
"emControl" <-
function(eps = .Machine$double.eps, tol = c(1.0e-05, 
	sqrt(.Machine$double.eps)), itmax = c(.Machine$integer.max, .Machine$
	integer.max), equalPro = FALSE)
{
	##
	# argList <- list(eps=eps, tol=tol, itmax=itmax, equalPro=equalPro)
	#	nullArgs <- sapply(argList, is.null)
	#	argList[nullArgs] <- .Mclust[names(nullArgs[nullArgs])]
	#	argList$itmax[argList$itmax == Inf] <- .Machine$integer.max
	#	argList
	##
	if(any(eps < 0)) stop("eps is negative")
	if(any(eps >= 1))
		stop("eps is not less than 1")
	if(any(tol < 0))
		stop("tol is negative")
	if(any(tol >= 1))
		stop("tol is not less than 1")
	if(any(itmax < 0))
		stop("itmax is negative")
	if(length(tol) == 1)
		tol <- rep(tol, 2)
	if(length(itmax) == 1)
		itmax <- c(itmax, .Machine$integer.max)
	i <- is.infinite(itmax)
	if(any(i))
		itmax[i] <- .Machine$integer.max
	list(eps = eps, tol = tol, itmax = itmax, equalPro = equalPro)
}

"mclustOptions" <-
function(emModelNames = NULL, hcModelNames = NULL, bicPlotSymbols = NULL, 
         bicPlotColors = NULL, classPlotSymbols = NULL, classPlotColors = NULL,
         warn = TRUE)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if(is.null(emModelNames)) emModelNames <- c("EII", "VII", "EEI", "VEI",
			"EVI", "VVI", "EEE", "EEV", "VEV", "VVV")
	if(is.null(hcModelNames))
		hcModelNames <- c("EII", "VII", "EEE", "VVV")
	if(is.null(bicPlotSymbols)) 	bicPlotSymbols <- 
           c(EII = 17, VII = 2, EEI = 16, EVI = 10, VEI = 13, VVI = 1,
             EEE = 15, EEV = 12, VEV = 7, VVV = 0, E = 17, V = 2)
	if(is.null(bicPlotColors)) 	bicPlotColors <- 
           c(EII = "gray", VII = "black", EEI = "orange", VEI = "brown", 
             EVI = "red", VVI = "magenta", EEE = "forestgreen",
             EEV = "green", "VEV" = "cyan", VVV = "blue", 
             E = "gray", V = "black")
	if(is.null(classPlotSymbols)) 	classPlotSymbols <- 
           c(17, 0, 10, 4, 11, 18, 6, 7, 3, 16, 2, 12, 8, 15, 1, 9, 14, 13, 5)
	if(is.null(classPlotColors)) 	classPlotColors <- 
           c("blue", "red", "green", "cyan", "magenta", "forestgreen",
             "purple", "orange", "gray", "brown", "black")
	list(emModelNames = emModelNames, hcModelNames = hcModelNames, 
             bicPlotSymbols = bicPlotSymbols, bicPlotColors = bicPlotColors,
             classPlotSymbols = classPlotSymbols, classPlotColors = 
             classPlotColors, warn = warn)
}

"priorControl" <- 
function(functionName = "defaultPrior", ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	c(list(functionName = functionName), list(...))
}

`clPairs` <-
function (data, classification, symbols=NULL, colors=NULL, 
          labels = dimnames(data)[[2]], CEX = 1, ...) 
{
    data <- as.matrix(data)
    m <- nrow(data)
    n <- ncol(data)
    if (missing(classification)) 
        classification <- rep(1, m)
    if (!is.factor(classification)) 
        classification <- as.factor(classification)
    l <- length(levels(classification))
    if (missing(symbols)) {
        if (l == 1) {
            symbols <- "."
        }
        if (l <= length(.Mclust$classPlotSymbols)) {
            symbols <- .Mclust$classPlotSymbols
        }
        else if (l <= 9) {
            symbols <- as.character(1:9)
        }
        else if (l <= 26) {
            symbols <- LETTERS[1:l]
        }
        else symbols <- rep( 16,l)
    }
    if (length(symbols) == 1) symbols <- rep(symbols, l)
    if (length(symbols) < l) {
        symbols <- rep( 16, l)
        warning("more symbols needed")
    }
    if (is.null(colors)) {
        if (l <= length(.Mclust$classPlotColors)) 
          colors <- .Mclust$classPlotColors[1:l]
    }
    if (length(colors) == 1) colors <- rep(colors, l)
    if (length(colors) < l) {
       colors <- rep( "black", l)
       warning("more colors needed")
    }
    pairs(x = data, labels = labels, pch = symbols[classification], 
        cex = CEX, col = colors[classification], ...)
    invisible()
}

`coordProj` <-
function(data, dimens = c(1,2), parameters = NULL, z = NULL, classification = 
  NULL, truth = NULL, uncertainty = NULL, what = c("classification", "errors", 
  "uncertainty"), quantiles = c(0.75, 0.94999999999999996), symbols = NULL, 
  colors = NULL, scale = FALSE, xlim = NULL, ylim = NULL, CEX = 1, PCH = ".", 
  identify = FALSE, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(is.null(dimens)) dimens <- c(1, 2)
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mean)
    sigma <- parameters$variance$sigma
    haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any(
      is.na(sigma))
  }
  else haveParams <- FALSE
  data <- data[, dimens, drop = FALSE]
  if(dim(data)[2] != 2)
    stop("need two dimensions")
  if(is.null(xlim))
    xlim <- range(data[, 1])
  if(is.null(ylim))
    ylim <- range(data[, 2])
  if(scale) {
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
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    noise <- classification[1] == "0"
    if(is.null(symbols)) {
      if(L <= length(.Mclust$classPlotSymbols)) {
        symbols <- .Mclust$classPlotSymbols
        if(noise) {
          first <- symbols[1]
          symbols[symbols == 16] <- first
          symbols[1] <- 16
        }
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
      if(L <= length(.Mclust$classPlotColors)) {
        colors <- .Mclust$classPlotColors[1:L]
        if(noise) {
          first <- colors[1]
          colors[colors == "black"] <- first
          colors[1] <- "black"
        }
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
  if(length(what) > 1)
    what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
    badClass <- (what == "errors" && (length(unique(classification)) != length(
      unique(truth))))
    if(badClass && !bad)
      warning("classification and truth differ in number of groups")
    bad <- bad && badClass
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad)
    what <- "bad"
  switch(EXPR = what,
    classification = {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), 
          "Coordinate Projection showing Classification")
        title(main = TITLE)
      }
      for(k in 1:L) {
        I <- classification == U[k]
        points(data[I, 1], data[I, 2], pch = symbols[k], col = colors[k], cex
           = if(U[k] == "0") CEX/4 else CEX)
      }
    }
    ,
    errors = {
      ERRORS <- classError(classification, truth)$misclassified
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), 
          "Coordinate Projection showing Errors")
        title(main = TITLE)
      }
      CLASSES <- unique(as.character(truth))
      symOpen <- c(2, 0, 1, 5)
      symFill <- c(17, 15, 16, 18)
      good <- rep(TRUE, length(classification))
      good[ERRORS] <- FALSE
      if(L > 4) {
        points(data[good, 1], data[good, 2], pch = 1, col = colors, cex = CEX)
        points(data[!good, 1], data[!good, 2], pch = 16, cex = CEX)
      }
      else {
        for(k in 1:L) {
          K <- truth == CLASSES[k]
          if(any(I <- (K & good))) {
            points(data[I, 1], data[I, 2], pch = symOpen[k], col = colors[k], 
              cex = CEX)
          }
          if(any(I <- (K & !good))) {
            points(data[I, 1], data[I, 2], pch = symFill[k], cex = CEX)
          }
        }
      }
    }
    ,
    uncertainty = {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), 
          "Coordinate Projection showing Uncertainty")
        title(main = TITLE)
      }
      breaks <- quantile(uncertainty, probs = sort(quantiles))
      I <- uncertainty <= breaks[1]
      points(data[I, 1], data[I, 2], pch = 16, col = "gray75", cex = 0.5 * CEX)
      I <- uncertainty <= breaks[2] & !I
      points(data[I, 1], data[I, 2], pch = 16, col = "gray50", cex = 1 * CEX)
      I <- uncertainty > breaks[2] & !I
      points(data[I, 1], data[I, 2], pch = 16, col = "black", cex = 1.5 * CEX)
    }
    ,
    {
      plot(data[, 1], data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim = 
        xlim, ylim = ylim, main = "", ...)
      if(identify) {
        TITLE <- paste(paste(dimens, collapse = ","), "Coordinate Projection")
        title(main = TITLE)
      }
      points(data[, 1], data[, 2], pch = PCH, cex = CEX)
    }
    )
  if(haveParams) {
    ## plot ellipsoids
    for(k in 1:G)
      mvn2plot(mu = mu[, k], sigma = sigma[,  , k], k = 15)
  }
  invisible()
}

"mclust1Dplot" <-
function(data, parameters=NULL, z=NULL,
         classification=NULL, truth=NULL, uncertainty=NULL, 
         what = c("classification", "density", "errors", "uncertainty"), 
         symbols=NULL, ngrid = length(data),  xlab=NULL,  xlim=NULL, CEX =1, 
         identify= FALSE, ...) 
{
grid1 <- function (n, range = c(0, 1), edge = TRUE) 
{
    if (any(n < 0 | round(n) != n)) 
        stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if (edge) {
        G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n - 
            1))
    }
    else {
        lj <- abs(diff(range))
        incr <- lj/(2 * n)
        G <- seq(from = min(range) + incr, to = max(range) - 
            incr, by = 2 * incr)
    }
    G
}
    densNuncer <- function(data, parameters) {
        cden <- cdensV(data = data, parameters = parameters)
        if (parameters$variance$G != 1) {
          z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = parameters$pro)
          den <- apply(z, 1, sum)
          z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
          data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
        }
        else {
          data.frame(density = cden, uncertainty =  rep(NA, length(cden)))
        }  
    }
    if (is.null(xlab)) xlab <- " "
    p <- ncol(as.matrix(data))
    if (p != 1) 
        stop("for one-dimensional data only")
    data <- as.vector(data)
    n <- length(data)
    if(is.null(classification) && !is.null(z))
                classification <- map(z)
    if(is.null(uncertainty) && !is.null(z))
                uncertainty <- 1 - apply(z, 1, max)
    if (!is.null(parameters)) {
      mu <- parameters$mean
      L <- ncol(mu)
      sigmasq <- parameters$variance$sigmasq
      haveParams <- !is.null(mu) && !is.null(sigmasq) && 
                    !any(is.na(mu)) && !any(is.na(sigmasq)) 
    }
    else haveParams <- FALSE
    if (is.null(xlim)) xlim <- range(data)
    if (haveParams) {
      G <- length(mu)
      if ((l <- length(sigmasq)) == 1) {
            sigmasq <- rep(sigmasq, G)
      }
      else if (l != G) {
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
    if (!is.null(classification)) {
      classification <- as.character(classification)
      U <- sort(unique(classification))
      L <- length(U)
      if (is.null(symbols)) {
        symbols <- rep("|", L)
      }
      else if (length(symbols) == 1) {
        symbols <- rep(symbols, L)
      }
      else if(length(symbols) < L) {
        warning("more symbols needed to show classification")
        symbols <- rep("|", L)
      }
    }
    if (length(what) > 1) what <- what[1]
    choices <- c("classification", "density", "errors", "uncertainty")
    m <- charmatch(what, choices, nomatch = 0)
    if (m) { 
      type <- choices[m] 
      bad <- what == "classification" && is.null(classification)
      bad <- bad || (what == "uncertainty" && is.null(uncertainty))
      bad <- bad || (what == "errors" && 
                (is.null(classification) || is.null(truth)))
      if (bad) warning("insufficient input for specified plot")
    }
    else {
      bad <- !m
      warning("what improperly specified")
    }
    if (bad) what <- "bad"
    M <- L
    switch (EXPR = what,
            "classification" = {
            plot(data, seq(from = 0, to = M, length = n), type = "n", 
           xlab = xlab, ylab = "", xlim = xlim, yaxt = "n", main = "", ...)
            if (identify) title("Classification")
            for (k in 1:L) {
                I <- classification == U[k]
                points(data[I], rep(0, length(data[I])), pch = symbols[k], 
                  cex = CEX)
                points(data[I], rep(k, length(data[I])), pch = symbols[k], 
                  cex = CEX)
            }
        },
        "errors" = {
            ERRORS <- classError(classification, truth)$misclassified
            plot(data, seq(from = 0, to = M, length = n), type = "n", 
             xlab = xlab, ylab = "", xlim = xlim, yaxt = "n", main = "", ...)
            if (identify) title("Classification Errors")
            good <- rep(TRUE, length(classification))
            good[ERRORS] <- FALSE
            sym <- "|"
            for (k in 1:L) {
                K <- classification == U[k]
                I <- K & good
                if (any(I)) {
                  if (FALSE) {
                    sym <- if (L > 4) 
                      1
                    else if (k == 4) 
                      5
                    else k - 1
                  }
                  l <- sum(as.numeric(I))
                  points(data[I], rep(0, l), pch = sym, cex = CEX)
                }
                I <- K & !good
                if (any(I)) {
                  if (FALSE) 
                    sym <- if (L > 5) 
                      16
                    else k + 14
                  l <- sum(as.numeric(I))
                  points(data[I], rep(k, l), pch = sym, cex = CEX)
                  points(data[I], rep(0, l), pch = sym, cex = CEX)
                  points(data[I], rep(-0.5, l), pch = sym, cex = CEX)
                }
            }
        },
        "uncertainty" = {
            x <- grid1(n = ngrid, range = xlim, edge = TRUE)
            lx <- length(x)
            Z <- densNuncer(data = x, parameters = parameters)
            plot(x, Z$uncertainty, xlab = xlab, ylab = "uncertainty", 
                xlim = xlim, ylim = c(0,1), type = "l", main = "", ...)
            if (identify)  title("Uncertainty")
        },
        "density" = {
            if (is.null(parameters$pro) && parameters$variance$G != 1) 
                stop("mixing proportions missing")
            x <- grid1(n = ngrid, range = xlim, edge = TRUE)
            lx <- length(x)
            Z <- densNuncer(data = x, parameters = parameters)
            plot(x, Z$density, xlab = xlab, ylab = "density", xlim = xlim, 
                type = "l", main = "", ...)
            if (identify) title("Density")
        },
        {
            plot(data, rep(0, n), type = "n", xlab = "", ylab = "", 
                xlim = xlim, main = "", ...)
            points(data, rep(0, n), pch = "|", cex = CEX)
            if (identify) title("Point Plot")
            return(invisible())
        }
       )
    invisible()
}

`mclust2Dplot` <-
function (data, parameters=NULL, z=NULL, classification=NULL, truth=NULL,
          uncertainty=NULL, what=c("classification", "uncertainty", "errors"), 
          quantiles = c(0.75, 0.95), symbols = NULL, colors = NULL, 
          scale = FALSE, xlim=NULL, ylim=NULL, CEX = 1, PCH = ".", 
          identify = FALSE, swapAxes = FALSE, ...) 
{
        if(dim(data)[2] != 2)
                stop("data must be two dimensional")
        if(is.null(classification) && !is.null(z))
                classification <- map(z)
        if(is.null(uncertainty) && !is.null(z))
                uncertainty <- 1 - apply(z, 1, max)
        if (!is.null(parameters)) {
          mu <- parameters$mean
          L <- ncol(mu)
          sigma <- parameters$variance$sigma
          haveParams <- !is.null(mu) && !is.null(sigma) && 
                         !any(is.na(mu)) && !any(is.na(sigma)) 
        }
        else haveParams <- FALSE
        if (is.null(xlim))
                xlim <- range(data[, 1])
        if (is.null(ylim))
                ylim <- range(data[, 2])
        if (scale) {
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
       if (haveParams) { 
          G <- ncol(mu)
          dimpar <- dim(sigma)
          if (length(dimpar) != 3) {
            haveParams <- FALSE
            warning("covariance must be a 3D matrix")
          }
          if (G != dimpar[3]) {
            haveParams <- FALSE
            warning("means and variance parameters are incompatible")
          }
          mu <- array(mu, c(2, G))
          sigma <- array(sigma, c(2, 2, G))
        }
      if (swapAxes) {
        if (haveParams) {
            mu <- mu[2:1, ]
            sigma <- sigma[2:1, 2:1, ]
        }
        data <- data[, 2:1]
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
        if (charmatch("classification", what, nomatch = 0) && 
            is.null(classification) && !is.null(z))
            classification <- map(z)  
        if (!is.null(classification)) {
                classification <- as.character(classification)
                U <- sort(unique(classification))
                L <- length(U)
                noise <- classification[1] == "0"
                if(is.null(symbols)) {
                        if(L <= length(.Mclust$classPlotSymbols)) {
                                symbols <- .Mclust$classPlotSymbols
                                if (noise) {
                                   first <- symbols[1]
                                   symbols[symbols == 16] <- first
                                   symbols[1] <- 16
                                }
                        }
                        else if(L <= 9) {
                                symbols <- as.character(1:9)
                        }
                        else if(L <= 26) {
                                symbols <- LETTERS
                        }
                }
                if(is.null(colors)) {
                        if(L <= length(.Mclust$classPlotColors)) {
                                colors <- .Mclust$classPlotColors[1:L]
                                if (noise) {
                                   first <- colors[1]
                                   colors[colors == "black"] <- first
                                   colors[1] <- "black"
                                }
                        }
                }
                else if (length(colors) == 1) colors <- rep(colors, L)
            if(length(symbols) < L) {
                warning("more symbols needed to show classification ")
                symbols <- rep(16,L)
            }
            if(length(colors) < L) {
                warning("more colors needed to show classification ")
                colors <- rep("black",L)
            }
        }
        if (length(what) > 1) what <- what[1]
        choices <- c("classification", "errors", "uncertainty")
        m <- charmatch(what, choices, nomatch = 0)
        if (m) { 
          what <- choices[m] 
          bad <- what == "classification" && is.null(classification)
          bad <- bad || (what == "uncertainty" && is.null(uncertainty))
          bad <- bad || (what == "errors" && 
                (is.null(classification) || is.null(truth)))
          if (bad) warning("insufficient input for specified plot")
        }
        else {
          bad <- !m
          warning("what improperly specified")
        }
        if (bad) what <- "bad"
        switch(EXPR = what,
            "classification"= {
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Classification")
            for (k in 1:L) {
                I <- classification == U[k]
                points(data[I, 1], data[I, 2], pch = symbols[k], 
                       col = colors[k], cex = if (U[k] == "0") CEX/4 else CEX)
            }
            },
            "errors" = {
            ERRORS <- classError(classification, truth)$misclassified
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Classification Errors")
            CLASSES <- unique(as.character(truth))
            symOpen <- c(2, 0, 1, 5)
            symFill <- c(17, 15, 16, 18)
            good <- rep(TRUE,length(classification))
            good[ERRORS] <- FALSE
            if (L > 4) {
                points(data[good, 1], data[good, 2], pch = 1, 
                  col = colors, cex = CEX)
                points(data[!good, 1], data[!good, 2], pch = 16, 
                  cex = CEX)
            }
            else {
                for (k in 1:L) {
                  K <- truth == CLASSES[k]
                  points(data[K, 1], data[K, 2], pch = symOpen[k], 
                    col = colors[k], cex = CEX)
                  if (any(I <- (K & !good))) {
                    points(data[I, 1], data[I, 2], pch = symFill[k], 
                      cex = CEX)
                  }
                }
            }
        },
        "uncertainty" = {
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Classification Uncertainty")
            breaks <- quantile(uncertainty, probs = sort(quantiles))
            I <- uncertainty < breaks[1]
            points(data[I, 1], data[I, 2], pch = 16, col = "gray75",
                   cex = 0.5 * CEX)
            I <- uncertainty < breaks[2] & !I
            points(data[I, 1], data[I, 2], pch = 16, col = "gray50",
                   cex = 1 * CEX)
            I <- uncertainty >= breaks[2]
            points(data[I, 1], data[I, 2], pch = 16, col = "black",
                   cex = 1.5 * CEX)
        },
        {
            plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
            if (identify) title("Point Plot")
            points(data[, 1], data[, 2], pch = PCH, cex = CEX)
          }
      )
    if (haveParams) {
## plot ellipsoids
          for (k in 1:G) 
             mvn2plot(mu = mu[, k], sigma = sigma[ ,  , k], k = 15)
    }
    invisible()
}

"mvn2plot" <-
function (mu, sigma, k = 15, alone = FALSE, col = 1) 
{
    p <- length(mu)
    if (p != 2) 
        stop("two-dimensional case only")
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
    if (alone) {
        xymin <- apply(xy, 2, FUN = "min")
        xymax <- apply(xy, 2, FUN = "max")
        r <- ceiling(max(xymax - xymin)/2)
        xymid <- (xymin + xymax)/2
        plot(xy[, 1], xy[, 2], xlim = c(-r, r) + xymid[1], ylim = c(-r, 
            r) + xymid[2], xlab = "x", ylab = "y", type = "n")
    }
    l <- length(x)
    i <- 1:l
    for (k in 1:4) {
        lines(xy[i, ], col = col)
        i <- i + l
    }
    x <- s[1]
    y <- s[2]
    xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
    lines(xy[1:2, ], lty = 2, col = col)
    lines(xy[3:4, ], lty = 2, col = col)
    points(mu[1], mu[2], pch = 8, col = col)
    invisible()
}

"plot.Mclust" <-
function (x, data = NULL, 
          what = c("BIC","classification","uncertainty","density"), 
          dimens = c(1,2), xlab = NULL, ylim = NULL,  
          legendArgs = list(x="bottomright", ncol=2, cex=1),
          identify = TRUE, ...) 
{
    parSave <- par(no.readonly = TRUE)
    on.exit(par(parSave))
    par(ask = TRUE)
    if (any(match("BIC", what, nomatch = 0)))
      plot.mclustBIC(x$BIC, xlab=xlab, ylim=ylim, legendArgs=legendArgs,...)
    # title("BIC")
    if (is.null(data)) {
        warning("data not supplied")
        return(invisible())
    }
    p <- ncol(as.matrix(data))
    if (p > 2) {
        if (any(match("classification", what, nomatch = 0)))
          clPairs(data[,1:min(5,p)], classification = x$classification, ...)
        if (any(match("classification", what, nomatch = 0)))
          coordProj(data = data, parameters=x$parameters, z = x$z,
                    what = "classification", dimens = dimens, identify = identify, 
                    ...)
        if (any(match("uncertainty", what, nomatch = 0)))
          coordProj(data = data, parameters=x$parameters, z = x$z,
                    what = "uncertainty", dimens = dimens, identify = identify, 
                    ...)
    }
    else if (p == 2) {
        if (any(match("classification", what, nomatch = 0)))
          mclust2Dplot(data = data, parameters=x$parameters, z=x$z,
                       what = "classification", identify = identify, ...)
        if (any(match("uncertainty", what, nomatch = 0)))
          mclust2Dplot(data = data, parameters=x$parameters, z=x$z,
                       what = "uncertainty", identify = identify, ...)
        if (any(match("density", what, nomatch = 0)))
           surfacePlot(data = data, parameters=x$parameters,
                     what = "density", identify = identify, ...)
    }
    else {
        if (any(match("classification", what, nomatch = 0)))
           mclust1Dplot(data = data, parameters=x$parameters, z=x$z,
                       what = "classification", identify = identify, ...)
        if (any(match("uncertainty", what, nomatch = 0)))
          mclust1Dplot(data = data, parameters=x$parameters, z=x$z,
                       what = "uncertainty", identify = identify, ...)
        if (any(match("density", what, nomatch = 0)))
          mclust1Dplot(data = data, parameters=x$parameters, z=x$z,
                       what = "density", identify = identify, ...)
    }
    invisible()
}

`plot.mclustBIC` <-
function (x, G = NULL, modelNames = NULL, symbols = NULL, colors = NULL, 
    xlab = NULL, ylim = NULL, legendArgs = list(x = "bottomright", 
        ncol = 2, cex = 1), CEX = 1, ...) 
{
    if (is.null(xlab)) 
        xlab <- "number of components"
    fill <- FALSE
    subset <- !is.null(attr(x, "initialization")$subset)
    noise <- !is.null(attr(x, "initialization")$noise)
    ret <- attr(x, "returnCodes") == -3
##
##    if (!subset && any(ret) && fill) {
##        x <- bicFill(x, ret, n, d)
##    }
##
    n <- ncol(x)
    dnx <- dimnames(x)
    x <- matrix(as.vector(x), ncol = n)
    dimnames(x) <- dnx
    if (is.null(modelNames)) 
        modelNames <- dimnames(x)[[2]]
    if (is.null(G)) 
        G <- as.numeric(dimnames(x)[[1]])
    if (is.null(symbols)) {
        colNames <- dimnames(x)[[2]]
        m <- length(modelNames)
        if (is.null(colNames)) {
            symbols <- if (m > 9) 
                LETTERS[1:m]
            else as.character(1:m)
            names(symbols) <- modelNames
        }
        else {
            symbols <- .Mclust$bicPlotSymbols[modelNames]
        }
    }
    if (is.null(colors)) {
        colNames <- dimnames(x)[[2]]
        if (is.null(colNames)) {
            colors <- 1:m
            names(colors) <- modelNames
        }
        else {
            colors <- .Mclust$bicPlotColors[modelNames]
        }
    }
    if (is.null(ylim)) 
        ylim <- range(as.vector(x[!is.na(x)]))
    xx <- as.numeric(dimnames(x)[[1]])
    matplot(xx, x, type = "b", xlim = range(xx), ylim = ylim, pch = symbols, 
        col = colors, lty = 1, xlab = xlab, ylab = "BIC", main = "")
    if (!is.null(legendArgs)) 
        do.call("legend", c(list(legend = modelNames, col = colors, 
            pch = symbols), legendArgs))
    invisible(symbols)
}

"plot.mclustDA" <-
function(x, trainData, testData, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if (missing(trainData) || missing(testData))
          stop("data not supplied")
        par(ask = TRUE)
	p <- ncol(as.matrix(trainData))
	if(p > 2) {
                dimens <- c(1,2)
		Data <- rbind(testData, trainData)
		xlim <- range((Data[, dimens])[, 1])
		ylim <- range((Data[, dimens])[, 2])
		cl <- c(rep(1, nrow(testData)), rep(2, nrow(trainData)))
		clPairs(Data[,1:min(5,ncol(Data))], cl = cl, 
                              symbols = c(1, 3), ...)
		coordProj(data = Data, 
                          classification = cl, what = "classification",
                          identify = FALSE, symbols = c(1,3), 
                          xlim=xlim, ylim=ylim, ...)
 	        title("Training and Test Data")
		coordProj(data = trainData, what ="classification",
                    classification = x$train$labels, identify = FALSE, 
                    xlim = xlim, ylim = ylim, ...)
  		title("Training Data: known Classification")
#	        coordProj(data = testData, what = "classification",
#		          classification = x$train$classification,
#                         identify = FALSE, xlim=xlim, ylim=ylim, ...)
#		title("Training Data: mclustDA Classification")
#		coordProj(data = trainData, 
#			  classification = x$train$classification, 
#                          truth	 = x$train$labels, what = "errors", 
#                          identify = FALSE, xlim = xlim, ylim = ylim, ...)
#		title("Training Error")
                if (!is.null(x$test$labels)) {
#		   coordProj(data = testData, what = "classification",
#                            classification = x$test$labels, identify = FALSE, 
#                            xlim = xlim, ylim = ylim, ...)
#                  title("Test Data: known Classification")
                }
		coordProj(data = testData, what = "classification",
			  classification = x$test$classification,
                          identify = FALSE, xlim=xlim, ylim=ylim, ...)
		title("Test Data: mclustDA Classification")
                if (!is.null(x$test$labels)) {
		   coordProj(data = testData, what = "errors",
                             classification = x$test$classification, 
                             truth = x$test$labels, identify = FALSE, 
                             xlim = xlim, ylim = ylim, ...)
  		    title("Test Error")
                }
	}
	else if (p == 2) {
                dimens <- c(1,2)
		Data <- rbind(testData, trainData)
		xlim <- range((Data[, dimens])[, 1])
		ylim <- range((Data[, dimens])[, 2])
		cl <- c(rep(1, nrow(testData)), rep(2, nrow(trainData)))
		mclust2Dplot(data = Data, classification = cl,
            what = "classification", identify = FALSE, symbols = c(1, 3), ...)
		title("Training and Test Data")
		mclust2Dplot(data = trainData, what = "classification",
                             classification = x$train$labels, identify = FALSE,
                             xlim = xlim, ylim = ylim, ...)
		title("Training Data: known Classification")
#		mclust2Dplot(data = trainData, what = "classification",
#			     classification = x$train$classification, 
#                             identify = FALSE, xlim = xlim, ylim = ylim, ...)
#		title("Training Data: mclustDA Classification")
#		mclust2Dplot(data = trainData, 
#			     classification = x$train$classification, 
#                            truth = x$train$labels, what = "errors", 
#                            identify = FALSE, xlim = xlim, ylim = ylim, ...)
# 	        title("Training Error")
                if (!is.null(x$test$labels)) {
#  	          mclust2Dplot(data = testData, what = "classification",
#                              classification = x$test$labels, 
#                              identify = FALSE, xlim = xlim, ylim = ylim, ...)
#		  title("Test Data: known Classification")
                }
		mclust2Dplot(data = testData, what = "classification",
			     classification = x$test$classification, 
                             identify = FALSE, xlim = xlim, ylim = ylim, ...)
		title("Test Data: mclustDA Classification")
                if (!is.null(x$test$labels)) {
		  mclust2Dplot(data = testData,
			       classification = x$test$classification, 
                               truth = x$test$labels, what = "errors", 
                              identify = FALSE, xlim = xlim, ylim = ylim, ...)
		  title("Test Error")
                }
	}
	else {
		##
		## p == 1
		##
		Data <- c(testData, trainData)
		xlim <- range(Data)
		cl <- c(rep(1, length(testData)), rep(2, length(trainData)))
		mclust1Dplot(data = Data, what = "classification",
                             classification = cl, identify = FALSE, 
                             xlim = xlim, ...)
		title("Training and Test Data")
		mclust1Dplot(data = trainData, what = "classification",
			     classification = x$train$labels, identify = FALSE,
                             xlim = xlim, ...)
		title("Training Data: known Classification")
#		mclust1Dplot(data = trainData, what = "classification",
#                            classification = x$train$classification, 
#                            identify = FALSE, xlim = xlim, ...)
#		title("Training Data: mclustDA Classification")
#	        mclust1Dplot(data = trainData,
#		             classification = x$train$classification, 
#                            truth = x$train$labels, what = "errors", 
#                            identify = FALSE, xlim = xlim, ...)
#		title("Training Error")
                if (!is.null(x$test$labels)) {
#  		  mclust1Dplot(data = testData, what = "classification",
#		               classification = x$test$labels, 
#                               what = "classification", 
#                              identify = FALSE, xlim = xlim, ...)
#		  title("Test Data: known classification")
                }
		mclust1Dplot(data = testData, what = "classification",
                             classification = x$test$classification, 
                             identify = FALSE, xlim = xlim, ...)
		title("Test Data: mclustDA Classification")
                if (!is.null(x$test$labels)) {
   		  mclust1Dplot(data = testData, 
			     classification = x$test$classification, 
                             truth = x$test$labels, what = "errors", 
                             identify = FALSE, xlim = xlim, ...)
		  title("Test Error")
                }
	}
	invisible()
}

`plot.mclustDAtrain` <-
function(x, data, dimens = c(1,2), symbols = NULL, colors = NULL,
         scale = FALSE, xlim = NULL, ylim = NULL,  CEX = 1, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	if (missing(data))
          stop("data not supplied")
	p <- ncol(as.matrix(data))
        if (p >= 2) {
          data <- data[, dimens, drop = FALSE]
          x <- lapply(x, function(z, dimens) {
           z$parameters$mean <- array(z$parameters$mean[dimens,], c(2,z$G))
           z$parameters$variance$sigma <- 
             array(z$parameters$variance$sigma[dimens,dimens,], c(2,2,z$G))
           z   
           }, dimens = dimens)
         }
	if (p >= 2) {
          L <- length(x)
          if (is.null(symbols)) {
            if (L <= length(.Mclust$classPlotSymbols)) {
              symbols <- .Mclust$classPlotSymbols
             }
            else if (L <= 9) {
                symbols <- as.character(1:9)
            }
            else if (L <= 26) {
                symbols <- LETTERS
            }
            if (length(symbols) == 1) symbols <- rep(symbols,L)
            if (length(symbols) < L) {
              warning("more symbols needed to show classification")
              symbols <- rep(16, L)
            }
           if (is.null(colors)) 
             colors <- .Mclust$classPlotColors
            }
            if (length(colors) == 1) colors <- rep(colors,L)
            if (length(colors) < L) {
              warning("more colors needed to show classification")
              colors <- rep("black", L)
            }
          if (is.null(xlim))
            xlim <- range(data[, 1])
          if (is.null(ylim))
            ylim <- range(data[, 2])
          if (scale) {
            par(pty = "s")
            d <- diff(xlim) - diff(ylim)
            if (d > 0) {
              ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
            }
            else {
              xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
            }
          }
          if (is.null(dnames <- dimnames(data)[[2]])) {
            xlab <- ylab <- ""
          }
          else {
              xlab <- dnames[1]
              ylab <- dnames[2]
          }
          plot(data[, dimens[1]], data[, dimens[2]], type = "n", 
               xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, 
               main = "", ...)
          for (l in 1:length(x)) { 
             I <- x[[l]]$observations
             points(data[I, dimens[1]], data[I, dimens[2]], pch = symbols[l], 
                    col = colors[l], cex = CEX)
             for (k in 1:(x[[l]]$G)) {
                mvn2plot(mu = x[[l]]$parameters$mean[, k], 
                        sigma = x[[l]]$parameters$variance$sigma[ ,  , k], 
                         k = 15)
             }
          }
	}
	else {
            if (is.null(xlim)) xlim <- range(data)
            plot(data, seq(from = 0, to = length(x), length = length(data)), type = "n", 
                 xlab = "", ylab = "", xlim = xlim, yaxt = "n", main = "", ...)
            sym <- "|"
            for (l in 1:length(x)) { 
               I <- x[[l]]$observations
               i <- length(I)
               points(data[I], rep(0,i), pch = sym, col = colors[l], cex = CEX)
               points(data[I], rep(l,i), pch = sym, col = colors[l], cex = CEX)
            }
	}
	invisible()
}

"print.mclustDA" <-
function (x, ndigits = options()$digits, ...) 
{
    cat("\nModeling Summary:\n")
    print(x$summary)
    cat("\nTest Classification Summary:")
    print(table(x$test$classification))
    cat("\nTraining Classification Summary:")
    print(table(x$train$classification))
    err <- classError(x$train$classification,x$train$labels)$errorRate
    cat("\nTraining Error:",err,"\n")
    if (!is.null(x$test$labels)) {
      err <- classError(x$test$classification,x$test$labels)$errorRate
      cat("\nTest Error:",err,"\n")
    }
    invisible()
}

`randProj` <-
function(data, seeds = 0, parameters = NULL, z = NULL, classification = NULL, 
  truth = NULL, uncertainty = NULL, what = c("classification", "errors", 
  "uncertainty"), quantiles = c(0.75, 0.94999999999999996), symbols = NULL, 
  colors = NULL, scale = FALSE, xlim = NULL, ylim = NULL, CEX = 1, PCH = ".", 
  identify = FALSE, ...)
{
  ##
  # This function is part of the MCLUST software described at
  #       http://www.stat.washington.edu/mclust
  # Copyright information and conditions for use of MCLUST are given at
  #        http://www.stat.washington.edu/mclust/license.txt
  ##
  if(scale) par(pty = "s")
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mu)
    sigma <- parameters$variance$sigma
    haveParams <- !is.null(mu) && !is.null(sigma) && !any(is.na(mu)) && !any(
      is.na(sigma))
  }
  else haveParams <- FALSE
  xlab <- ylab <- ""
  p <- ncol(data)
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
    cho <- array(apply(sigma, 3, chol), c(p, p, G))
  }
  if(!is.null(truth)) {
    if(is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if(length(unique(truth)) != length(unique(classification)))
        truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    noise <- classification[1] == "0"
    if(is.null(symbols)) {
      if(L <= length(.Mclust$classPlotSymbols)) {
        symbols <- .Mclust$classPlotSymbols
        if(noise) {
          first <- symbols[1]
          symbols[symbols == 16] <- first
          symbols[1] <- 16
        }
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
      if(L <= length(.Mclust$classPlotColors)) {
        colors <- .Mclust$classPlotColors[1:L]
        if(noise) {
          first <- colors[1]
          colors[colors == "black"] <- first
          colors[1] <- "black"
        }
      }
    }
    else if(length(colors) == 1)
      colors <- rep(colors, L)
    if(length(symbols) < L) {
      warning("more symbols needed to show classification ")
      symbols <- rep(16,L)
    }
    if (length(colors) < L) {
      warning("more colors needed to show classification ")
      colors <- rep("black",L)
    }
  }
  if(length(what) > 1)
    what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) {
    what <- choices[m]
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && (is.null(classification) || is.null(
      truth)))
    if(bad)
      warning("insufficient input for specified plot")
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if(bad)
    what <- "bad"
  nullXlim <- is.null(xlim)
  nullYlim <- is.null(ylim)
  if(length(seeds) > 1)
    par(ask = TRUE)
  for(seed in seeds) {
    set.seed(seed)
    Q <- orth2(p)
    Data <- as.matrix(data) %*% Q
    if(dim(Data)[2] != 2)
      stop("need two dimensions")
    if(nullXlim)
      xlim <- range(Data[, 1])
    if(nullYlim)
      ylim <- range(Data[, 2])
    if(scale) {
      d <- diff(xlim) - diff(ylim)
      if(d > 0) {
        ylim <- c(ylim[1] - d/2, ylim[2] + d/2.)
      }
      else {
        xlim <- c(xlim[1] + d/2, xlim[2] - d/2)
      }
    }
    switch(EXPR = what,
      classification = {
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        for(k in 1:L) {
          I <- classification == U[k]
          points(Data[I, 1], Data[I, 2], pch = symbols[k], col = colors[k], cex
             = CEX)
        }
        if(identify) {
          TITLE <- paste("Random Projection showing Classification: seed = ", 
            seed)
          title(TITLE)
        }
      }
      ,
      errors = {
        ERRORS <- classError(classification, truth)$misclassifiedPoints
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        if(identify) {
          TITLE <- paste("Random Projection showing Errors: seed = ", seed)
          title(TITLE)
        }
        CLASSES <- unique(as.character(truth))
        symOpen <- c(2, 0, 1, 5)
        symFill <- c(17, 15, 16, 18)
        good <- !ERRORS
        if(L > 4) {
          points(Data[good, 1], Data[good, 2], pch = 1, col = colors, cex = CEX
            )
          points(Data[!good, 1], Data[!good, 2], pch = 16, cex = CEX)
        }
        else {
          for(k in 1:L) {
            K <- truth == CLASSES[k]
            points(Data[K, 1], Data[K, 2], pch = symOpen[k], col = colors[k], 
              cex = CEX)
            if(any(I <- (K & ERRORS))) {
              points(Data[I, 1], Data[I, 2], pch = symFill[k], cex = CEX)
            }
          }
        }
      }
      ,
      uncertainty = {
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        if(identify) {
          TITLE <- paste("Random Projection showing Uncertainty: seed = ", seed
            )
          title(TITLE)
        }
        breaks <- quantile(uncertainty, probs = sort(quantiles))
        I <- uncertainty <= breaks[1]
        points(Data[I, 1], Data[I, 2], pch = 16, col = "gray75", cex = 0.5 * 
          CEX)
        I <- uncertainty <= breaks[2] & !I
        points(Data[I, 1], Data[I, 2], pch = 16, col = "gray50", cex = 1 * CEX)
        I <- uncertainty > breaks[2] & !I
        points(Data[I, 1], Data[I, 2], pch = 16, col = "black", cex = 1.5 * CEX
          )
      }
      ,
      {
        plot(Data[, 1], Data[, 2], type = "n", xlab = xlab, ylab = ylab, xlim
           = xlim, ylim = ylim, main = "", ...)
        if(identify) {
          TITLE <- paste("Random Projection: seed = ", seed)
          title(TITLE)
        }
        points(Data[, 1], Data[, 2], pch = PCH, cex = CEX)
      }
      )
    if(haveParams) {
      ## plot ellipsoids
      muTrans <- crossprod(Q, mu)
      sigmaTrans <- array(apply(cho, 3, function(R, Q)
      crossprod(R %*% Q), Q = Q), c(2, 2, G))
      for(k in 1:G)
        mvn2plot(mu = muTrans[, k], sigma = sigmaTrans[,  , k], k = 15)
    }
  }
  invisible()
}

"summary.mclustBIC" <-
function(object, data, G, modelNames, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        mc <- match.call(expand.dots = FALSE)
        if (is.null(attr(object,"initialization")$noise)) {
          mc[[1]] <- as.name("summaryMclustBIC")
        }
        else {
          mc[[1]] <- as.name("summaryMclustBICn")
        }
       ans <- eval(mc, parent.frame())
       Glabels <- dimnames(object)[[1]]
       if (length(Glabels) != 1 && (!missing(G) && length(G) > 1)) {
         Grange <- range(as.numeric(Glabels))
         if (match(ans$G, Grange, nomatch = 0))
  warning("best model occurs at the min or max # of components considered")
       }
       ans
}

"surfacePlot" <-
function( data, parameters, type = c("contour", "image", "persp"), 
          what = c("density", "uncertainty"), 
          transformation = c("none", "log", "sqrt"), 
          grid = 50, nlevels = 20, scale = FALSE, 
          xlim=NULL, ylim=NULL, identify = FALSE, verbose = FALSE, 
          swapAxes = FALSE, ...)
{
grid1 <- function (n, range = c(0, 1), edge = TRUE) 
{
    if (any(n < 0 | round(n) != n)) 
        stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if (edge) {
        G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n - 
            1))
    }
    else {
        lj <- abs(diff(range))
        incr <- lj/(2 * n)
        G <- seq(from = min(range) + incr, to = max(range) - 
            incr, by = 2 * incr)
    }
    G
}
grid2 <- function (x, y) 
{
    lx <- length(x)
    ly <- length(y)
    xy <- matrix(0, nrow = lx * ly, ncol = 2)
    l <- 0
    for (j in 1:ly) {
        for (i in 1:lx) {
            l <- l + 1
            xy[l, ] <- c(x[i], y[j])
        }
    }
    xy
}

    data <- as.matrix(data)
    if(dim(data)[2] != 2)
                stop("data must be two dimensional")
     densNuncer <- function(modelName, data, parameters) {
        if (is.null(parameters$variance$cholsigma)) {
          parameters$variance$cholsigma <- parameters$variance$sigma
          G <- dim(parameters$variance$sigma)[3]
          for (k in 1:G)
             parameters$variance$cholsigma[,,k] <- 
                 chol(parameters$variance$sigma[,,k])
        }
        cden <- cdensVVV(data = data, parameters = parameters)
        z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = parameters$pro)
        den <- apply(z, 1, sum)
        z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
        data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
    }
          pro <- parameters$pro
          mu <- parameters$mean
          sigma <- parameters$variance$sigma
          haveParams <- !is.null(mu) && !is.null(sigma) && !is.null(pro) &&
              !any(is.na(mu)) && !any(is.na(sigma))  && !(any(is.na(pro)))
        if (haveParams) { 
          G <- ncol(mu)
          dimpar <- dim(sigma)
          if (length(dimpar) != 3) {
            haveParams <- FALSE
            warning("covariance must be a 3D matrix")
          }
          if (G != dimpar[3]) {
            haveParams <- FALSE
            warning("means and variance parameters are incompatible")
          }
          mu <- array(mu, c(2, G))
          sigma <- array(sigma, c(2, 2, G))
        }
     if (!haveParams) 
        stop("need parameters to compute density")
      if (swapAxes) {
        if (haveParams) {
            parameters$pro <- pro[2:1]
            parameters$mean <- mu[2:1, ]
            parameters$variance$sigma <- sigma[2:1, 2:1, ]
        }
        data <- data[, 2:1]
      }
        if (is.null(xlim))
                xlim <- range(data[, 1])
        if (is.null(ylim))
                ylim <- range(data[, 2])
        if (scale) {
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
    if (length(grid) == 1) 
        grid <- c(grid, grid)
    x <- grid1(n = grid[1], range = xlim, edge = TRUE)
    y <- grid1(n = grid[2], range = ylim, edge = TRUE)
    xy <- grid2(x, y)
    if (verbose) 
        cat("\n computing density and uncertainty over grid ...\n")
    Z <- densNuncer(modelName = "VVV", data = xy, parameters = parameters)
    lx <- length(x)
    ly <- length(y)
    CI <- type
    DU <- what
    TRANS <- transformation
    if (length(CI) > 1) CI <- CI[1]
    if (length(DU) > 1) DU <- DU[1]
    if (length(TRANS) > 1) TRANS <- TRANS[1]
    switch(DU, 
       "density" = {
            zz <- matrix(Z$density, lx, ly)
            title2 <- "Density"
        },
        "uncertainty" = {
            zz <- matrix(Z$uncertainty, lx, ly)
            title2 <- "Uncertainty"
        },
        stop("what improperly specified")
    )
    switch(TRANS, 
       "none" = {
            title1 <- ""
        },
        "log" = {
             zz <- logb(zz)
             title1 <- "log"
         },
        "sqrt" = {
              zz <- sqrt(zz)
              title1 <- "sqrt"
         },
        stop("transformation improperly specified")
    )
    switch(CI,
           "contour" = {
            title3 <- "Contour"
            contour(x = x, y = y, z = zz, nlevels = nlevels, xlab = xlab, 
                ylab = ylab, main = "", ...)
        },
        "image" = {
            title3 <- "Image"
            image(x = x, y = y, z = zz, xlab = xlab, ylab = ylab, main="", ...)
        },
        "persp" = {
            title3 <- "Perspective"
            persp(x = x, y = y, z = zz, xlab = xlab, ylab = ylab, 
                theta = 60, phi = 30, expand = 0.6, main = "", ...)
        },
       stop("type improperly specified")
    )
    if (identify) {
      TITLE <- paste(c(title1, title2, title3, "Plot"), collapse = " ")
            title(TITLE)
    }
    invisible(list(x = x, y = y, z = zz))
}

"uncerPlot" <-
function (z, truth=NULL, ...) 
{
    parSave <- par(no.readonly = TRUE)
    par(pty = "m")
    uncer <- 1 - apply(z, 1, max)
    ord <- order(uncer)
    M <- max(uncer)
    plot(uncer[ord], ylab = "uncertainty", ylim = c(-(M/32), M), 
         xaxt = "n", xlab = "observations in order of increasing uncertainty", 
         type = "n")
    points(uncer[ord], pch = 15, cex = 0.5)
    lines(uncer[ord])
    abline(h = c(0, 0), lty = 4)
    if (!is.null(truth)) {
        truth <- as.numeric(as.factor(truth))
        n <- length(truth)
        result <- map(z)
        bad <- classError(result, truth)$misclassified
        if (length(bad)) {
            for (i in bad) {
                x <- (1:n)[ord == i]
                lines(c(x, x), c(-(0.5/32), uncer[i]), lty = 1)
            }
        }
    }
    par(parSave)
    invisible()
}

"decomp2sigma" <- 
function(d, G, scale, shape, orientation = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	nod <- missing(d)
	noG <- missing(G)
	lenScale <- length(scale)
	if(lenScale != 1) {
		if(!noG && G != lenScale)
			stop("scale incompatibile with G")
		G <- lenScale
	}
	shape <- as.matrix(shape)
	p <- nrow(shape)
	if(!nod && p != d)
		stop("shape incompatible with d")
	d <- p
	g <- ncol(shape)
	if(g != 1) {
		if(!is.null(G) && g != G)
			stop("shape incompatible with scale")
		if(!noG && g != G)
			stop("shape incompatible with G")
		G <- g
	}
	if(is.null(orientation)) {
		orientName <- "I"
		if(is.null(G)) {
			G <- if(noG) 1 else G
		}
		orientation <- array(diag(d), c(d, d, G))
	}
	else {
		dimO <- dim(orientation)
		l <- length(dimO)
		if(is.null(dimO) || l < 2 || l > 3 || dimO[1] != dimO[2])
			stop("orientation improperly specified")
		if(dimO[1] != d)
			stop("orientation incompatible with shape")
		if(l == 3) {
			orientName <- "V"
			if(is.null(G)) {
				if(!noG && dimO[3] != G)
					stop("orientation incompatible with G")
				G <- dimO[3]
			}
			else if(G != dimO[3])
				stop("orientation incompatible with scale and/or shape"
					)
		}
		else {
			orientName <- "E"
			if(is.null(G)) {
				G <- if(noG) 1 else G
			}
			orientation <- array(orientation, c(d, d, G))
		}
	}
	if(G == 1) {
		scaleName <- shapeName <- "X"
	}
	else {
		scaleName <- if(lenScale == 1) "E" else "V"
		shapeName <- if(g == 1) "E" else "V"
		scale <- rep(scale, G)
		shape <- matrix(shape, nrow = d, ncol = G)
	}
	sigma <- array(0, c(d, d, G))
	for(k in 1:G) {
		sigma[,  , k] <- crossprod(orientation[,  , k] * sqrt(scale[
			k] * shape[, k]))
	}
	structure(sigma, modelName = paste(c(scaleName, shapeName, orientName),
		collapse = ""))
}

"grid1" <-
function (n, range = c(0, 1), edge = TRUE) 
{
    if (any(n < 0 | round(n) != n)) 
        stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if (edge) {
        G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n - 
            1))
    }
    else {
        lj <- abs(diff(range))
        incr <- lj/(2 * n)
        G <- seq(from = min(range) + incr, to = max(range) - 
            incr, by = 2 * incr)
    }
    G
}
"grid2" <-
function (x, y) 
{
    lx <- length(x)
    ly <- length(y)
    xy <- matrix(0, nrow = lx * ly, ncol = 2)
    l <- 0
    for (j in 1:ly) {
        for (i in 1:lx) {
            l <- l + 1
            xy[l, ] <- c(x[i], y[j])
        }
    }
    xy
}
"hypvol" <-
function(data, reciprocal = FALSE)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	## finds the minimum hypervolume between principal components and 
	## variable bounds
	##
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	if(oneD) {
		## 1D case
		n <- length(as.vector(data))
		if(reciprocal) {
			ans <- 1/(max(data) - min(data))
		}
		else {
			ans <- max(data) - min(data)
		}
		return(ans)
	}
	if(length(dimdat) != 2)
		stop("data must be a vector or a matrix")
	data <- as.matrix(data)
	dimd <- dim(data)
	n <- dimd[1]
	p <- dimd[2]
##
#		vol1 <- prod(apply(data, 2, function(z)
#		diff(range(z))))
#		V <- matrix(temp[[1]], p, p)
#		xbar <- apply(data, 2, mean)
#		X <- sweep(data, 2, xbar)
#		library(Matrix)
#		print(V)
#		print(eigen.Hermitian(crossprod(X))$vectors)
#		X <- X %*% V
#		vol <- prod(apply(X, 2, function(z)
#		diff(range(z))))
##
	lwgesvd <- max(3 * min(n, p) + max(n, p), 5 * min(n, p) - 4)
	# min
	lwsyevd <- p * (3 * p + 2 * ceiling(logb(p, base = 2)) + 5) + 1
	# minimum
	lisyevd <- 5 * p + 3
	# minimum
	lwsyevx <- 8 * p
	lisyevx <- 5 * p + p
	lwork <- max(lwsyevd, lwsyevx, n)
	liwork <- max(lisyevd, lisyevx)
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
		integer(1),
                PACKAGE = "mclust")[c(4, 11)]
	if(temp[[2]])
		stop("problem in computing principal components")
	if(reciprocal) {
		pcvol <- prod(1/temp[[1]])
		bdvol <- prod(1/(apply(data, 2, max) - apply(data, 2, min)))
		ans <- max(pcvol, bdvol)
	}
	else {
		pcvol <- sum(log(temp[[1]]))
		bdvol <- sum(log(apply(data, 2, max) - apply(data, 2, min)))
		maxlog <- log(.Machine$double.xmax)
		if(pcvol > maxlog || bdvol > maxlog) {
			warning("hypervolume greater than largest machine representable number"
				)
			ans <- Inf
		}
		else {
			ans <- min(pcvol, bdvol)
		}
	}
	ans
}

"map" <- 
function(z, warn = TRUE, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	nrowz <- nrow(z)
	cl <- numeric(nrowz)
	I <- 1:nrowz
	J <- 1:ncol(z)
	for(i in I) {
		cl[i] <- (J[z[i,  ] == max(z[i,  ])])[1]
	}
	if(warn) {
		K <- as.logical(match(J, sort(unique(cl)), nomatch = 0))
		if(any(!K))
			warning(paste("no assignment to", paste(J[!K], collapse
				 = ",")))
	}
	cl
}

"orth2" <-
function (n) 
{
    u <- rnorm(n)
    u <- u/vecnorm(u)
    v <- rnorm(n)
    v <- v/vecnorm(v)
    Q <- cbind(u, v - sum(u * v) * u)
    dimnames(Q) <- NULL
    Q
}
"partconv" <- 
function(x, consec = TRUE)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	n <- length(x)
	y <- numeric(n)
	u <- unique(x)
	if(consec) {
		# number groups in order of first row appearance
		l <- length(u)
		for(i in 1:l)
			y[x == u[i]] <- i
	}
	else {
		# represent each group by its lowest-numbered member
		for(i in u) {
			l <- x == i
			y[l] <- (1:n)[l][1]
		}
	}
	y
}

"partuniq" <-
function(x)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	# finds the classification that removes duplicates from x
"charconv" <- function(x, sep = "001")
{
        if(!is.data.frame(x)) x <- data.frame(x)
        do.call("paste", c(as.list(x), sep = sep))
}

	n <- nrow(x)
	x <- charconv(x)
	k <- duplicated(x)
	partition <- 1.:n
	partition[k] <- match(x[k], x)
	partition
}

"shapeO" <-
function(shape, O, transpose = FALSE)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
##
  dimO <- dim(O)
  if(dimO[1] != dimO[2])
    stop("leading dimensions of O are unequal")
  if((ldO <- length(dimO)) != 3) {
    if(ldO == 2) {
      dimO <- c(dimO, 1)
      O <- array(O, dimO)
    }
    else stop("O must be a matrix or an array")
  }
  l <- length(shape)
  if(l != dimO[1])
    stop("dimension of O and length s are unequal")
  storage.mode(O) <- "double"
  .Fortran("shapeo",
           as.logical(transpose),
           as.double(shape),
           O,
           as.integer(l),
           as.integer(dimO[3]),
           double(l * l),
           integer(1),
           PACKAGE="mclust")[[3]]
}

"sigma2decomp" <-
function (sigma, G=NULL, tol=NULL, ...) 
{
    dimSigma <- dim(sigma)
    if (is.null(dimSigma)) 
        stop("sigma improperly specified")
    d <- dimSigma[1]
    if (dimSigma[2] != d) 
        stop("sigma improperly specified")
    l <- length(dimSigma)
    if (l < 2 || l > 3) 
        stop("sigma improperly specified")
    if (is.null(G)) {
        if (l == 2) {
            G <- 1
            sigma <- array(sigma, c(dimSigma, 1))
        }
        else {
            G <- dimSigma[3]
        }
    }
    else {
        if (l == 3 && G != dimSigma[3]) 
            stop("sigma and G are incompatible")
        if (l == 2 && G != 1) 
            sigma <- array(sigma, c(d,d,G))
    }
    decomp <- list(d = d, G = G, scale = rep(0, G), shape = matrix(0, 
        d, G), orientation = array(0, c(d, d, G)))
    for (k in 1:G) {
        ev <- eigen(sigma[, , k], symmetric = TRUE)
        temp <- log(ev$values)
        logScale <- sum(temp)/d
        decomp$scale[k] <- exp(logScale)
        decomp$shape[, k] <- exp(temp - logScale)
        decomp$orientation[, , k] <- t(ev$vectors)
    }
    if (is.null(tol)) 
        tol <- sqrt(.Machine$double.eps)
    scaleName <- "V"
    shapeName <- "V"
    orientName <- "V"
    uniq <- function(x, tol = sqrt(.Machine$double.eps)) {
        abs(max(x) - min(x)) < tol
    }
    if (uniq(decomp$scale)) {
        decomp$scale <- decomp$scale[1]
        scaleName <- "E"
    }
    if (all(apply(decomp$shape, 1, uniq, tol = tol))) {
        decomp$shape <- decomp$shape[, 1]
        if (all(uniq(decomp$shape, tol = tol))) {
            shapeName <- "I"
            decomp$shape <- rep(1, d)
        }
        else {
            shapeName <- "E"
        }
    }
    if (all(apply(matrix(decomp$orientation, nrow = d * d, ncol = G), 
        1, uniq, tol = tol))) {
        decomp$orientation = decomp$orientation[, , 1]
        if (all(apply(cbind(decomp$orientation, diag(d)), 1, 
            uniq, tol = tol))) {
            orientName <- "I"
            decomp$orientation <- NULL
        }
        else {
            orientName <- "E"
        }
    }
    modelName <-  paste(c(scaleName, shapeName, orientName), collapse = "")
    c(list(modelName = modelName, decomp))
}

"traceW" <- function(x)
{
        ##
        # This function is part of the MCLUST software described at
        #       http://www.stat.washington.edu/mclust
        # Copyright information and conditions for use of MCLUST are given at
        #        http://www.stat.washington.edu/mclust/license.txt
        ##
        # sum(as.vector(sweep(x, 2, apply(x, 2, mean)))^2)
        dimx <- dim(x)
        n <- dimx[1]
        p <- dimx[2]
        .Fortran("mcltrw",
                as.double(x),
                as.integer(n),
                as.integer(p),
                double(p),
                double(1),
                PACKAGE = "mclust")[[5]]
}
"unchol" <-
function(x, upper = NULL)
{
##
# This function is part of the MCLUST software described at
#       http://www.stat.washington.edu/mclust
# Copyright information and conditions for use of MCLUST are given at
#        http://www.stat.washington.edu/mclust/license.txt
##
  if(is.null(upper)) {
    upper <- any(x[row(x) < col(x)])
    lower <- any(x[row(x) > col(x)])
    if(upper && lower)
      stop("not a triangular matrix")
    if(!(upper || lower)) {
      x <- diag(x)
      return(diag(x * x))
    }
  }
  dimx <- dim(x)
  storage.mode(x) <- "double"
  .Fortran("unchol",
           as.logical(upper),
           x,
           as.integer(nrow(x)),
           as.integer(ncol(x)),
           integer(1),
           PACKAGE="mclust")[[2]]
}

"unmap" <-
function(classification, noise=NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	# converts a classification to conditional probabilities
	# classes are arranged in sorted order
	# if a noise indicator is specified, that column is placed last
	n <- length(classification)
	u <- sort(unique(classification))
	labs <- as.character(u)
	k <- length(u)
	if(!is.null(noise)) {
		l <- u == noise
		if(any(l)) {
			m <- max(u) + 1
			u[l] <- m
			labs <- labs[order(u)]
			u <- sort(u)
			classification[classification == noise] <- m
		}
	}
	z <- matrix(0, n, k)
	for(j in 1:k)
		z[classification == u[j], j] <- 1
	##
	## z <- matrix(1, n, k)
	## for(j in 1:k) z[classification == u[j], j] <- k + 1
	## z <- sweep(z, 1, apply(z, 1, sum), "/")
	##
	dimnames(z) <- list(NULL, labs)
	z
}

"vecnorm" <-
function (x, p = 2) 
{
    if (is.character(p)) {
        if (charmatch(p, "maximum", nomatch = 0) == 1) 
            p <- Inf
        else if (charmatch(p, "euclidean", nomatch = 0) == 1) 
            p <- 2
        else stop("improper specification of p")
    }
    if (!is.numeric(x) && !is.complex(x)) 
        stop("mode of x must be either numeric or complex")
    if (!is.numeric(p)) 
        stop("improper specification of p")
    if (p < 1) 
        stop("p must be greater than or equal to 1")
    if (is.numeric(x)) 
        x <- abs(x)
    else x <- Mod(x)
    if (p == 2) 
        return(.Fortran("d2norm", as.integer(length(x)), as.double(x), 
            as.integer(1), double(1), PACKAGE = "mclust")[[4]])
    if (p == Inf) 
        return(max(x))
    if (p == 1) 
        return(sum(x))
    xmax <- max(x)
    if (!xmax) 
        xmax <- max(x)
    if (!xmax) 
        return(xmax)
    x <- x/xmax
    xmax * sum(x^p)^(1/p)
}
"[.mclustBIC" <-
function (x, i, j, drop = FALSE) 
{
    ATTR <- attributes(x)[c("G", "modelNames", "prior", "control", 
                            "initialization", "Vinv", "warn", "n", "d", 
                            "oneD", "returnCodes", "class")]
    oldClass(x) <- NULL
    x <- NextMethod("[")
    if (is.null(dim(x))) return(x)
    ATTR$G <- as.numeric(dimnames(x)[[1]])
    ATTR$modelNames <- dimnames(x)[[2]]
    ATTR$returnCodes <- ATTR$returnCodes[dimnames(x)[[1]],dimnames(x)[[2]],
                                          drop=FALSE]
    do.call("structure", c(list(.Data = x), ATTR))
}

"[.mclustDAtest" <-
function (x, i, j, drop = FALSE) 
{
    clx <- oldClass(x)
    oldClass(x) <- NULL
    NextMethod("[")
}
"bic" <-
function(modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	modelName <- switch(EXPR = modelName,
		X = "E",
		XII = "EII",
		XXI = "EEI",
		XXX = "EEE",
  		modelName)
        checkModelName(modelName)
        if (G == 0) {
                ## one cluster case
                if(!noise) stop("undefined model")
                nparams <- 1
        }
        else {
                nparams <- nVarParams(modelName, d, G) + G*d
                if(!equalPro)
                        nparams <- nparams + (G - 1)
                if(noise)
                        nparams <- nparams + 2
        }
        2 * loglik - nparams * logb(n)
}

"cdens" <-
function(modelName, data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
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

"checkModelName" <-
function(modelName)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##	
        switch(EXPR = modelName,
		E = ,
		V = ,
	        EII  = ,
		VII = ,
		EEI  = ,
		VEI = ,
		EVI = ,
		VVI = ,
		EEE = ,
		VEE = ,
		EVE = ,
		VVE = ,
		EEV = ,
		VEV = ,
		VVV = TRUE,
		stop("invalid model name"))
}

"dens" <-
function(modelName, data, logarithm = FALSE, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
	aux <- list(...)
	cden <- cdens(modelName = modelName, data = data,
                  logarithm = TRUE, parameters = parameters, warn = NULL)
	dimdat <- dim(data)
	oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
	G <- if (oneD) {
               length(parameters$mean) 
             }
             else {
               ncol(as.matrix(parameters$mean))
             }
	if(G > 1) {
		pro <- parameters$pro
		if(is.null(pro))
			stop("mixing proportions must be supplied")
		#
		##	if(is.null(pro <- aux$pro)) pro <- rep(1/G, G)
		##	apply(sweep(cden, 2, FUN = "*", STATS = pro), 1, sum)
		proz <- !pro
		pro <- pro[!proz]
		cden <- cden[, !proz, drop = FALSE]
		cden <- sweep(cden, 2, FUN = "+", STATS = log(pro))
	}
	maxlog <- apply(cden, 1, max)
	cden <- sweep(cden, 1, FUN = "-", STATS = maxlog)
	den <- logb(apply(exp(cden), 1, sum)) + maxlog
	if (!logarithm)	den <- exp(den)
	den
}

"em" <-
function(modelName, data, parameters, prior = NULL, control = emControl(), 
         warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        checkModelName(modelName)
	funcName <- paste("em", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

"estep" <-
function(modelName, data, parameters, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        checkModelName(modelName)
	funcName <- paste("estep", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

"hc" <-
function(modelName, data, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##        
        switch(EXPR = modelName,
                E = ,
                V = ,
                EII  = ,
                VII = ,
                EEE = ,
                VVV = TRUE,
                stop("invalid model name for hierarchical clustering"))
        funcName <- paste("hc", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc[[2]] <- NULL
        eval(mc, parent.frame())
}

"mclustModelNames" <-
structure(list(univariateMixture = c("E", "V"), multivariateMixture = structure(list(
    spherical = c("EII", "VII"), diagonal = c("EEI", "EVI", "VEI", 
    "VVI"), ellipsoidal = c("EEE", "EEV", "VEV", "VVV")), .Names = c("spherical", 
"diagonal", "ellipsoidal")), singleComponent = c("X", "XII", 
"XXI", "XXX")), .Names = c("univariateMixture", "multivariateMixture", 
"singleComponent"))
"mclustVariance" <-
function(modelName, d=NULL, G=2) 
{
 x <- -1
 if (nchar(modelName) == 1) {
   if (!is.null(d) && d != 1)  stop("modelName and d are incompatible")
   varList <- switch(EXPR = modelName,
                     "X" = list(sigmasq = x),
                     "E" = list(sigmasq = x),
                     "V" = list(sigmasq = rep(x,G)),
                     stop("modelName not recognized"))
 }
 else {
   if (nchar(modelName) != 3) stop("modelName is misspecified")
   if (is.null(d)) d <- 3
   varList <- switch(EXPR = modelName,
                     "XII" = list(sigmasq = x),
                     "EII" = list(sigmasq = x),
                     "VII" = list(sigmasq = rep(x,G)),
                     "XXI" = list(scale = x, shape = rep(x,d)),
                     "EEI" = list(scale = x, shape = rep(x,d)),
                     "EVI" = list(scale = x, shape = matrix(x,d,G)),
                     "VEI" = list(scale = rep(x,G), shape = rep(x,d)),
                     "VVI" = list(scale = rep(x,G), shape = matrix(x,d,G)),
                     "XXX" = {M <- matrix(x,d,d); M[row(M) > col(M)] <- 0;
                              list(cholSigma = M)},
                     "EEE" = {M <- matrix(x,d,d); M[row(M) > col(M)] <- 0;
                              list(cholSigma = M)},
                     "EEV" = list(scale = x, shape = rep(x,d),
                                  orientation = array(x,c(d,d,G))),
                     "VEV" = list(scale = x, shape = matrix(x,d,G),
                                  orientation = array(x,c(d,d,G))),
                     "VVV" = {A <- array(x,c(d,d,G)); 
                              I <- row(A[,,1]) > col(A[,,1])
                              for (k in 1:G) A[,,k][I] <- 0
                              list(cholsigma = A)},
                     stop("modelName not recognized"))
  }
 c(modelName = modelName, d=d, G=G, varList)
}

"me" <-
function(modelName, data, z, prior = NULL, control = emControl(), 
         Vinv = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        checkModelName(modelName)
	funcName <- paste("me", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

"mstep" <-
function(modelName, data, z, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        checkModelName(modelName)
	funcName <- paste("mstep", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

"mvn" <-
function(modelName, data, prior = NULL, warn = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##	
        modelName <- switch(EXPR = modelName,
		E = "X",
		V = "X",
		X  =  "X",
		Spherical = "XII",
	        EII  = "XII",
		VII = "XII",
		XII = "XII",
		Diagonal =  "XXI",
		EEI  = "XXI",
		VEI = "XXI",
		EVI = "XXI",
		VVI = "XXI",
		XXI = "XXI",
		Ellipsoidal = "XXX",
		EEE = "XXX",
		VEE = "XXX",
		EVE = "XXX",
		VVE = "XXX",
		EEV = "XXX",
		VEV = "XXX",
		VVV = "XXX",
		XXX = "XXX",
		stop("invalid model name"))
        funcName <- paste("mvn", modelName, sep = "")
        mc <- match.call()
        mc[[1]] <- as.name(funcName)
        mc[[2]] <- NULL
        eval(mc, parent.frame())
}

"nVarParams" <-
function(modelName, d, G)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        r <- (d*(d-1))/2
        s <- (d*(d+1))/2 # s = r + d
	switch(EXPR = modelName,
		E = 1,
		V = G,
		EII = 1,
		VII = G,
		EEI = d,
		VEI = G+(d-1),
		EVI = 1 + G*(d-1),
		VVI = G*d,
		EEE = s,
		VEE = G + (d-1) + r,
		EVE = 1 + (d-1) + r,
		VVE = G + G*(d-1) + r,
		EEV = 1 + (d-1) + G*r,
		VEV = G + (d-1) + G*r,
		VVV = G * s,
		stop("invalid model name"))
}

"sim" <-
function(modelName, parameters, n, seed = NULL, ...)
{
	##
	# This function is part of the MCLUST software described at
	#       http://www.stat.washington.edu/mclust
	# Copyright information and conditions for use of MCLUST are given at
	#        http://www.stat.washington.edu/mclust/license.txt
	##
        modelName <- switch(EXPR = modelName,
                X = "E",
                XII = "EII",
                XXI = "EEI",
                XXX = "EEE",
                modelName)
        checkModelName(modelName)
	funcName <- paste("sim", modelName, sep = "")
        mc <- match.call(expand.dots = TRUE)
        mc[[1]] <- as.name(funcName)
        mc$modelName <- NULL
        eval(mc, parent.frame())
}

"chevron" <- 
structure(.Data = list(structure(.Data = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2)
, levels = c("data", "noise")
, class = "factor"
)
, c(54.234892748296261, 48.235562283080071, 55.687680731993169, 
	42.526220000814646, 57.036088942084461, 38.920804302673787, 
	50.352763144299388, 22.635154556017369, 45.188236500136554, 
	50.147577424999326, 62.931766451802105, 47.640057321172208, 
	32.879657312761992, 24.837501123547554, 86.163869962561876, 
	86.837749335449189, 30.155830167233944, 85.877189304679632, 
	52.519144220743328, 57.793714497238398, 43.241628542076796, 
	36.721126145217568, 76.442784171551466, 64.163636977784336, 
	39.956034869886935, 69.729254734702408, 56.731933660805225, 
	47.354641291312873, 96.121449000202119, 82.79412979632616, 
	27.851137142628431, 24.93720111111179, 70.782046268228441, 
	48.560211469419301, 93.043078976683319, 32.209129346068949, 
	51.966823963914067, 86.239495593588799, 66.112952909898013, 
	58.24884533463046, 58.838225212879479, 56.636853196658194, 
	28.494492165744305, 57.653920940123498, 58.695394218433648, 
	86.625092410948128, 79.084716598968953, 78.508364506531507, 
	67.891362857073545, 73.143478830344975, 53.011723873205483, 
	30.373510002158582, 73.353310863021761, 29.196813782211393, 
	47.826015655882657, 95.469707596115768, 66.416002587880939, 
	36.85956151643768, 83.042443388141692, 54.88353905454278, 
	20.405374596826732, 56.198116395156831, 88.915851167403162, 
	59.562612885143608, 24.462514431215823, 31.448000574950129, 
	84.677187889814377, 87.364191023632884, 42.840826678089797, 
	92.397713405545801, 34.250688253901899, 69.683547101449221, 
	63.546513228211552, 43.672808268573135, 41.922648540697992, 
	56.473112017847598, 23.171860761009157, 49.953823303803802, 
	46.436614512931556, 38.59219811623916, 36.842365437187254, 
	26.409177205059677, 37.328651405405253, 34.927986776456237, 
	23.656566813588142, 54.236644501797855, 52.210703836753964, 
	28.157244871836156, 40.096746440976858, 31.42181457253173, 
	39.851036574691534, 31.424492257647216, 21.042729932814837, 
	29.523033264558762, 15.452032738830894, 30.015379111282527, 
	17.606630416121334, 57.101884072180837, 20.443842625245452, 
	28.15221022348851, 49.629977764561772, 51.881898459978402, 
	57.992464969865978, 41.258769547566772, 43.240769172552973, 
	29.272754858247936, 22.174956740345806, 15.878613241948187, 
	47.602522533852607, 15.149760278873146, 29.244020255282521, 
	10.318592390976846, 15.792438546195626, 58.153129788115621, 
	20.091362772509456, 17.016000929288566, 55.212599416263402, 
	34.973146319389343, 53.71798048261553, 42.072467103134841, 
	33.85334680788219, 42.941173664294183, 30.652374024502933, 
	24.599534531589597, 58.054738822393119, 49.216255359351635, 
	58.694438498932868, 33.398411432281137, 33.487216115463525, 
	56.60395851591602, 40.388811649754643, 45.681617558002472, 
	44.670640942640603, 47.806818750686944, 63.221794744022191, 
	21.259777629747987, 46.282551544718444, 28.109394870698452, 
	39.038281014654785, 26.695298594422638, 52.581904833205044, 
	38.105826680548489, 53.261818196624517, 22.38856595242396, 
	34.513444658368826, 52.018882452975959, 46.146133805159479, 
	27.751676877960563, 22.251768377609551, 23.039224445819855, 
	27.827243870124221, 47.182311916258186, 33.193872570991516, 
	33.675142102874815, 47.399089089594781, 29.215442202985287, 
	19.19504533521831, 39.283584766089916, 45.919479306321591, 
	24.246974098496139, 40.75248803012073, 18.43075955985114, 
	38.122169103007764, 18.892768386285752, 33.499488248489797, 
	28.507504318840802, 46.171165551058948, 55.270138657651842, 
	53.00451789284125, 26.972963795997202, 59.249994799029082, 
	21.945157558657229, 56.327855265699327, 40.006134547293186, 
	37.262823963537812, 25.617648903280497, 55.879329103045166, 
	24.799662989098579, 40.955960997380316, 43.778304925654083, 
	29.831844971049577, 56.157954586669803, 64.189819009043276, 
	39.882790127303451, 43.470695824362338, 41.053441364783794, 
	19.858526606112719, 27.550494831521064, 60.940338368527591, 
	55.592591632157564, 47.432575717102736, 38.313838448375463, 
	63.655514353886247, 23.842276630457491, 55.173727651126683, 
	36.56159303849563, 51.171865740325302, 14.076206854078919, 
	28.135600308887661, 40.992705882526934, 33.077532607130706, 
	27.613630839623511, 31.261784746311605, 54.491466227918863, 
	17.863415980245918, 45.307226125150919, 14.00991789996624, 
	47.494816437829286, 30.550595805980265, 22.459614060353488, 
	25.26177371153608, 38.942916460800916, 53.661402179859579, 
	56.124480329453945, 60.420438817236573, 74.254679393488914, 
	74.963943171314895, 76.253521561156958, 72.367728634271771, 
	73.941699247807264, 84.805357372388244, 46.228766483254731, 
	65.132539374753833, 86.463511495385319, 60.991641394793987, 
	71.728394830133766, 64.568556053563952, 70.021648365072906, 
	96.391001576557755, 78.706610652152449, 86.667805800680071, 
	71.702183936722577, 84.603021314833313, 59.024730327073485, 
	86.668122820556164, 82.75913733523339, 56.638916553929448, 
	51.76590992603451, 60.289750411175191, 70.936234463006258, 
	69.140980679076165, 82.953922473825514, 94.843533814419061, 
	95.314143276773393, 55.126623443793505, 92.814025529660285, 
	74.851102968677878, 89.899483653716743, 89.917693706229329, 
	56.741508990526199, 83.786134459078312, 96.84357285965234, 
	47.818922619335353, 73.453598842024803, 60.395266306586564, 
	74.826042328495532, 78.906501308083534, 60.655096885748208, 
	73.928653816692531, 79.443247432354838, 49.473912711255252, 
	60.718888929113746, 59.773725408595055, 67.382725598290563, 
	86.062307746615261, 59.336203166749328, 66.709851818159223, 
	64.263534042984247, 58.793802275322378, 68.099827175028622, 
	55.757865277118981, 91.672143219038844, 66.455129287205637, 
	79.538362873718143, 62.841279662679881, 78.063955535180867, 
	59.176526232622564, 70.252592810429633, 64.215119602158666, 
	92.512990746181458, 77.417090851813555, 61.613753747660667, 
	62.781374778132886, 84.352812813594937, 97.196518606506288, 
	87.331147640943527, 75.889226933941245, 60.499696799088269, 
	74.571702843531966, 70.27152115944773, 65.689407796598971, 
	76.708736689761281, 95.170646188780665, 63.154367236420512, 
	63.813885564450175, 95.110226231627166, 69.104425115510821, 
	86.728514351416379, 62.054581318516284, 94.769142719451338, 
	79.102635909803212, 87.765462999232113, 54.075649841688573, 
	55.873192274011672, 61.642944214399904, 88.682003612630069, 
	53.23315066518262, 82.474950519390404, 60.009071207605302, 
	77.290643295273185, 69.260751297697425, 74.735402576625347, 
	63.513273554854095, 86.321564957033843, 75.078163505531847, 
	58.225159978028387, 73.75684337457642, 61.371734300628304, 
	54.869736698456109, 69.747709559742361, 67.061143447645009, 
	78.089036212768406, 98.98887749761343, 79.364319217856973, 
	54.671592689119279, 55.085043162107468, 57.299722910393029, 
	77.168781748041511, 55.688775470480323, 85.181349173653871, 
	56.649081571958959, 79.601930344942957, 50.526558889541775, 
	87.531525117810816, 79.453271641395986, 64.015951626934111, 
	70.957273826934397, 73.659181636758149, 71.422531078569591, 
	52.113710390403867, 97.270094121340662, 60.998607259243727, 
	86.876463498920202, 68.728787831496447, 70.580466347746551, 
	82.039276428986341, 88.508475383277982, 62.857153776567429, 
	52.871364126913249, 18.260733188129961, 48.242243095766753, 
	104.85621612565592, 1.6468111644499004, 76.292511423118412, 
	82.846835286356509, 115.03907246515155, 8.4224839634262025, 
	110.26245710160583, 94.396102156490088, 98.428357252851129, 
	124.09517906373367, 2.5983997522853315, 28.010582229122519, 
	119.88305351790041, 37.371480740141124, 39.761460340116173, 
	120.57179889362305, 116.0976435332559, 12.473358022514731, 
	100.81065120967105, 23.259111578576267, 57.104839666280895, 
	41.63804620411247, 125.79464095178992, 126.63052967749536, 
	51.325166821014136, 105.67170533211902, 77.81598101882264, 
	15.22361704101786, 16.829490541480482, 116.5133208623156, 
	69.092953691259027, 100.55104284081608, 84.581971760839224, 
	125.43950975919142, 15.455725486855954, 49.709088711533695, 
	54.316982565913349, 98.122927467338741, 25.547892076894641, 
	54.029593718703836, 99.194448021240532, 95.67859789961949, 
	53.914788392838091, 12.175311088562012, 36.264842866454273, 
	126.89444229984656, 109.88714680075645, 3.6363264066167176, 
	61.86572879133746, 9.9464483158662915, 111.91885753022507, 
	102.0582111896947, 49.773802845273167, 33.060074502136558, 
	34.408025745302439, 8.7543834708631039, 6.7590066683478653, 
	118.4957725442946, 7.5494405799545348, 21.147015653550625, 
	13.96845042752102, 89.82311094738543, 96.392688883468509, 
	100.12316969316453, 68.562505691777915, 33.892352260649204, 
	42.616237458307296, 124.94237693725154, 51.896007396746427, 
	6.1892119063995779, 21.306909190956503, 100.04009515233338, 
	50.746443669311702, 62.667456711176783, 4.1698864176869392, 
	60.298391849268228, 26.091923768632114, 23.167692231480032, 
	116.79595653200522, 87.406535304617137, 102.55936067225412, 
	87.005080546252429, 71.44628586852923, 36.024712040089071, 
	114.48299978068098, 119.21788525814191, 6.5931212874129415, 
	73.500397203024477, 124.14050844404846, 17.019308641087264, 
	9.7509031379595399, 93.795005307998508, 88.002286639530212, 
	3.4815936917439103, 62.028843510895967, 5.9062781170941889, 
	112.7428733361885, 74.688497961964458, 124.67752790823579, 
	103.71663541346788, 114.77148305391893, 122.43437082925811, 
	58.048309536185116, 48.08694654982537, 23.48470267187804, 
	86.648323563858867, 104.64580322988331, 85.281665418297052, 
	2.5967782204970717, 122.72191030345857, 37.825218854472041, 
	73.602769385557622, 106.07903028512374, 80.29449285985902, 
	64.573924379888922, 3.7934052241034806, 14.411271899472922, 
	9.1130882548168302, 121.82028404390439, 64.086428255308419, 
	30.142027483321726, 9.2218365068547428, 57.943212111014873, 
	88.525624772999436, 73.029497192241251, 39.642025673761964, 
	47.75100915087387, 59.43292783619836, 86.025916869752109, 
	106.87721258634701, 74.76745021995157, 50.658329186029732, 
	99.179475331678987, 80.779239672236145, 31.022781821433455, 
	27.701226088218391, 70.397150602657348, 59.157695125322789, 
	52.216560912784189, 15.965767580084503, 4.5535093531943858, 
	40.315483149141073, 119.72016133833677, 70.820953687187284, 
	55.254787453450263, 107.04329040553421, 58.29693005932495, 
	60.950731799006462, 94.778017394710332, 16.855303465854377, 
	5.1008030329830945, 125.20530501380563, 19.123279816005379, 
	117.05310684815049, 30.921934833284467, 15.594612206798047, 
	36.559576065745205, 73.97708597406745, 9.1608772254548967, 
	60.407161391340196, 30.619238510727882, 39.713751916773617, 
	96.817864642944187, 86.451950440183282, 63.908302642405033, 
	8.8742578355595469, 60.937965170945972, 38.38533392501995, 
	70.242368211038411, 127.44251992600039, 97.17243448831141, 
	113.34552459372208, 45.340859668329358, 47.084083457943052, 
	123.78408719366416, 79.461673388257623, 76.060394993517548, 
	83.992182818241417, 100.19052438251674, 85.134706522803754, 
	32.949363013263792, 56.675030249170959, 4.1729824617505074, 
	82.929926445242018, 54.178752406965941, 111.72389138210565, 
	61.557107215747237, 100.60162192955613, 89.411813451442868, 
	1.5963736544363201, 77.989755954127759, 37.082840334624052, 
	122.86444625351578, 99.470816673710942, 73.503071467857808, 
	77.225126214325428, 6.604564858134836, 89.609252910129726, 
	78.241926767863333, 40.839448155835271, 68.662037074100226, 
	3.8544323379173875, 23.432667226530612, 119.54267649631947, 
	74.391020122449845, 66.262053582817316, 19.436904431320727, 
	98.661847421899438, 121.38916345639154, 3.3399431114085019, 
	19.744615266565233, 11.624164030887187, 25.567141579464078, 
	56.174077581148595, 77.977524415589869, 46.894562827423215, 
	67.996121652424335, 110.04782931134105, 123.87353478791192, 
	86.362637150567025, 26.791294113267213, 71.1956223892048, 
	1.7693870570510626, 8.8763646027073264, 113.83987859310582, 
	8.1820756569504738, 4.6885575391352177, 77.906489034183323, 
	15.124028530437499, 16.82679202966392, 41.167970235925168, 
	39.351054778322577, 92.120821432210505, 68.75347268441692, 
	84.129290046170354, 123.46324638044462, 11.986010336782783, 
	36.466133574023843, 63.905688048806041, 47.197478673886508, 
	65.146065894979984, 104.04342418909073, 51.899855570401996, 
	88.870257258415222, 16.011577168945223, 68.754992792848498, 
	93.34779418585822, 28.541898921132088, 98.468859472777694, 
	41.606486862991005, 62.116508715320379, 30.079980872571468, 
	14.848973709624261, 53.583268321119249, 51.051180050242692, 
	98.171318827662617, 96.896620204672217, 111.75066561065614, 
	121.21503516566008, 62.102458298206329, 117.31530580017716, 
	95.633917333092541, 81.658740204758942, 59.490287091117352, 
	118.55465504201129, 75.581577219069004, 6.9790993737988174, 
	55.863192815799266, 47.785586945712566, 44.529544668737799, 
	45.547190791461617, 35.004935494624078, 61.310466902330518, 
	26.708701082039624, 25.201865449082106, 19.590230466332287, 
	109.16138148587197, 19.669482204131782, 57.935748239047825, 
	39.934332637116313, 61.158326249103993, 93.898421884980053, 
	25.0810880321078, 1.1837924285791814, 78.590424932073802, 
	21.96231254003942, 113.68255837028846, 69.974941077176481, 
	120.12450284790248, 81.634451469406486, 5.2411539158783853, 
	110.4009982724674, 95.545921482611448, 117.5203898595646, 
	97.689060288015753, 23.467674252111465, 2.748198501765728, 
	6.6315603847615421, 71.713750023394823, 71.053446122445166, 
	60.598172856494784, 123.81770013598725, 5.4371770550496876, 
	104.46403920929879, 84.399017982184887, 95.456913054455072, 
	80.965890228748322, 109.38921471545473, 26.762251431588084, 
	37.649710851255804, 121.72102635679767, 126.67810095753521, 
	43.483440926298499, 121.78167054569349, 29.73739477340132, 
	88.853728622198105, 95.963571037631482, 36.41423639561981, 
	104.0495452512987, 78.288063809275627, 19.371748292818666, 
	1.5300123225897551, 8.3644131482578814, 104.9899572217837, 
	94.865336695220321, 72.51575589645654, 4.1911874515935779, 
	64.119399007409811, 13.819417232181877, 103.48210836062208, 
	68.552552068606019, 95.352681246586144, 81.033778943121433, 
	109.84260302269831, 15.163799851667136, 19.409917893819511, 
	1.6820975034497678, 86.994094415567815, 59.279218757990748, 
	123.63031288515776, 110.5170444524847, 2.6766205793246627, 
	12.684948138892651, 75.433121575042605, 122.52482425933704, 
	106.65533459326252, 111.56037309812382, 84.331093725282699, 
	102.10965347103775, 117.70885698823258, 20.62237522052601, 
	34.097953490447253, 101.26351714227349, 103.69528239639476, 
	37.422555768396705, 72.108009223360568, 121.90218550385907, 
	35.984306571539491, 15.215452017262578, 51.837688966654241, 
	103.8142516980879, 11.267733495682478, 123.40852295001969, 
	127.43777437740937, 37.49273966671899, 6.8147043543867767, 
	34.755856880918145, 69.906424306798726, 42.680059125646949, 
	36.846497123595327, 34.352577625773847, 67.987395754549652, 
	3.5804628948681056, 97.458876638673246, 62.403605297673494, 
	24.263261712156236, 53.560122328344733, 52.885797698050737, 
	111.79160557733849, 72.723893656395376, 54.686022360343486, 
	107.16449734894559, 55.267980591859668, 124.39840580383316, 
	45.689928760286421, 65.10774543043226, 78.407875664066523, 
	96.507699007634073, 125.61199676571414, 40.713500986341387, 
	117.64885599492118, 45.13085891213268, 39.144257285166532, 
	51.140286944806576, 68.564957298338413, 28.072209012228996, 
	33.947671634610742, 67.218429885338992, 100.59878586046398, 
	71.113587689120322, 76.079066116828471, 35.112949891015887, 
	116.12001450965181, 13.554063984192908, 41.743181419093162, 
	109.03411462996155, 93.742289054673165, 51.076567114796489, 
	63.012148394715041, 14.186432403512299, 121.90845459094271, 
	4.7878580428659916, 25.204745103605092, 126.70518590230495, 
	94.604904875159264, 15.336439732462168, 34.994519759435207, 
	54.699474409222603, 2.1904001417569816, 111.89486372098327, 
	23.308798969723284, 25.591879297979176, 57.072093228343874, 
	45.181771958712488, 22.431650512851775, 14.472638057079166, 
	21.756868501659483, 1.7477369848638773, 46.290041595231742, 
	112.94135494204238, 25.64110055565834, 44.971031255554408, 
	67.092865049373358, 10.63843137351796, 42.277234117034823, 
	48.672111710999161, 116.27852579159662, 110.81964662577957, 
	41.385295828338712, 38.072848681360483, 51.678193079773337, 
	36.094593974761665, 115.86631702305749, 18.309762539807707, 
	68.309295037295669, 89.182596073951572, 74.504600916523486, 
	17.972932083066553, 15.980133503675461, 126.9110757894814, 
	16.171689165756106, 117.12328152079135, 27.26665338082239, 
	53.044957084115595, 4.617964044213295, 101.03790410840884, 
	44.606407842133194, 88.691944470629096, 6.4330017799511552, 
	41.898747816216201, 47.175609255209565, 126.01479918323457, 
	66.758705249987543, 58.376265360508114, 11.086736034601927, 
	11.030157355125993, 41.216917916666716, 121.85634823795408, 
	6.8825527359731495, 65.465182333718985, 83.111177435144782, 
	52.452948359772563, 9.1498894388787448, 14.269479030743241, 
	47.595099097117782, 112.13049998460338, 116.07430108077824, 
	6.0370199205353856, 38.151698983740062, 92.737055672332644, 
	120.86318689631298, 57.287359778769314, 118.79266786482185, 
	49.239546886645257, 118.0956919984892, 114.77072282228619, 
	24.904520731419325, 99.054818219970912, 110.38396541727707, 
	95.475155420135707, 32.069896012544632, 13.677838505711406, 
	83.00076205516234, 70.97665379755199, 55.328131975606084, 
	7.9869723655283451, 11.523691274225712, 31.495840937364846, 
	20.604578053113073, 125.19503529276699, 47.596286844462156, 
	76.019509848672897, 48.157864126376808, 116.5377921005711, 
	49.164441619534045, 47.285962361842394, 84.865094321779907, 
	109.33376204222441, 14.316439345944673, 108.73915001982823, 
	44.551680212840438, 121.51971181621775, 87.250408151652664, 
	121.70577121945098, 47.787018937058747, 123.2182313837111, 
	35.028488304466009, 47.08707742812112, 60.103482178878039, 
	100.84473324241117, 96.709089009556919, 76.630582911428064, 
	77.592385038267821, 122.52542410604656, 72.508548037149012, 
	102.56203067908064, 46.310730540659279, 22.537071191705763, 
	4.1330515244044363, 102.90699747717008, 69.452451351564378, 
	14.738065110519528, 29.667451511602849, 35.88278167322278, 
	34.731176704633981, 7.5309895719401538, 86.908354658167809, 
	11.805155282374471, 82.665265182498842, 32.997436915524304, 
	53.888008841779083, 83.23318748595193, 69.892676030751318, 
	90.730076037812978, 37.473621216136962, 127.20602490520105, 
	29.535287474747747, 109.34459417499602, 7.859430986456573, 
	104.13339983532205, 93.895296862814575, 35.378559447359294, 
	17.274222793988883, 112.83955960115418, 79.045637499541044, 
	59.226491859648377, 100.95573297468945, 102.18857787456363, 
	116.51885751308873, 107.04330383008346, 74.539201479870826, 
	114.01300944574177, 127.49163325503469, 86.395918145775795, 
	56.897064783144742, 12.838676023297012, 6.8202548436820507, 
	11.682935288175941, 93.948530109133571, 61.740293410606682, 
	3.0686158798635006, 59.252136829309165, 62.359675263985991, 
	125.81377933220938, 68.792252954561263, 99.920908810570836, 
	85.574392623268068, 117.45267974492162, 80.532076333649457, 
	67.656891335733235, 27.311675179284066, 50.9422103818506, 
	76.833226186688989, 15.926986836828291, 25.256496267858893, 
	1.6712627685628831, 37.518133236560971, 90.936042805667967, 
	47.754298283718526, 99.053256300278008, 68.140343063510954, 
	100.42655545799062, 43.250397020019591, 53.686695897951722, 
	54.387864126823843, 113.96069220313802, 39.537609825376421, 
	35.504880375228822, 46.30884944787249, 17.218782066833228, 
	61.886631110217422, 104.71265198523179, 45.434518727008253, 
	29.029845034237951, 88.705212538130581, 103.89363685343415, 
	111.74549112701789, 85.929819569922984, 78.448190472554415, 
	112.35312785999849, 105.10448920726776, 71.397151604294777, 
	121.96557190828025, 56.93980978988111, 25.294881371315569, 
	3.6926564066670835, 34.967766820918769, 103.30494464328513, 
	77.03356321901083, 76.635911806952208, 72.96371441707015, 
	27.069402043242007, 10.617021760437638, 70.410816557239741, 
	113.8048245520331, 10.498747933190316, 126.10665698675439, 
	119.87614507926628, 76.764509881380945, 45.719567503780127, 
	18.124324036762118, 42.609524000901729, 29.660280318465084, 
	125.27662083227187, 123.40803913399577, 68.861474188510329, 
	27.371462917421013, 119.41296566277742, 21.388890074566007, 
	7.3319276352412999, 123.09525157185271, 91.953839023131877, 
	39.443951417226344, 125.84959058882669, 120.18227855954319, 
	38.929774562828243, 98.623340787831694, 69.045814539771527, 
	47.660000819712877, 53.13523676386103, 13.191091790329665, 
	122.00679414533079, 112.7695490391925, 45.675307001452893, 
	73.082477146293968, 79.783281884621829, 99.721362643875182, 
	8.3743955129757524, 49.882082361727953, 111.73017028626055, 
	66.770922181196511, 73.236168541945517, 7.1926889792084694, 
	48.820143505930901, 121.72270194161683, 8.7430376568809152, 
	101.28646525600925, 24.302910142578185, 89.778497385326773, 
	35.051496089436114, 63.365828080102801, 108.21273630438372, 
	65.046511330176145, 107.93921152735129, 22.546639879234135, 
	71.300482785329223, 45.963531245943159, 121.55102502088994, 
	26.931785392109305, 57.460004154127091, 29.97740755835548, 
	67.428827523253858, 126.03811366297305, 111.33692224463448, 
	33.805739079602063, 73.186330228112638, 95.211517202667892, 
	41.343767133075744, 18.037236738484353, 35.619029317516834, 
	95.217294844798744, 86.480003785807639, 110.24814067734405, 
	20.788849980570376, 97.057657016906887, 117.95269003324211, 
	11.816918676719069, 96.746365848463029, 43.315061063040048, 
	57.50252029299736, 15.764123576227576, 84.346445022616535, 
	95.762169148772955, 37.572667067404836, 41.099415498320013, 
	2.5994864902459085, 15.284549473319203, 1.6876453314907849, 
	7.6129189454950392, 3.1049187588505447, 47.758530328515917, 
	80.731851249933243, 125.15207892330363, 52.50024718651548, 
	109.18804352777079, 68.070850263349712, 119.50526943337172, 
	69.075300113297999, 118.62118628015742, 100.21520160185173, 
	108.88847565744072, 86.729670181870461, 64.036794798914343, 
	61.030501062050462, 90.623023596126586, 17.425544608850032, 
	85.020715777296573, 13.209796918556094, 19.354592605959624, 
	115.24600189737976, 64.559205160010606, 27.832398009486496, 
	4.8956242920830846, 19.512852606363595, 73.795337684452534, 
	45.51142584066838, 25.107413691468537, 74.064806059468538, 
	68.486655332613736, 50.806855205446482, 64.279663819354028, 
	112.59897831315175, 110.95997669175267, 1.4558871067129076, 
	16.062978762667626, 67.404507799074054, 76.246016886550933, 
	118.46167040430009, 73.428288508206606, 80.329179943539202, 
	57.935969418846071, 82.847283323295414, 58.561363932676613)
, c(36.597284914925694, 40.698521614074707, 50.629234757740051, 
	37.862097983714193, 39.875036876183003, 32.654850579798222, 
	52.743632549419999, 19.151394074782729, 33.131827239412814, 
	49.661798519082367, 50.189789352938533, 33.444164694752544, 
	15.609109182376415, 19.741804568096995, 14.192305614240468, 
	20.925381786189973, 28.668452797457576, 15.127970912493765, 
	53.624140632804483, 44.274768250528723, 34.713419219478965, 
	24.449981816112995, 35.270084713120013, 47.488972298800945, 
	38.275700125843287, 45.098942134995013, 54.668989978265017, 
	44.501173805911094, 23.196823325939476, 25.374713617376983, 
	12.261850818991661, 5.8538755369372666, 40.726678393781185, 
	52.831362747587264, 16.877141304314137, 23.788769717793912, 
	40.814581041224301, 27.046392838936299, 46.253369156736881, 
	61.208049123641104, 55.82941138651222, 40.56769672036171, 
	9.4526651259511709, 61.886471257545054, 41.35168637521565, 
	17.617923899553716, 29.457988624460995, 36.676529378630221, 
	46.932702858932316, 32.006159091368318, 42.240199570544064, 
	19.164341490715742, 38.592847236432135, 28.566003083717078, 
	38.204341431614012, 19.854599759913981, 35.622397341765463, 
	28.054100892040879, 31.422245074529201, 62.222894483245909, 
	16.682205139659345, 41.318396627437323, 11.22821647208184, 
	41.668097087182105, 14.231764977332205, 22.588176617398858, 
	19.057407031301409, 25.257744455244392, 37.717769826762378, 
	17.011566579341888, 22.600028838496655, 41.415341802872717, 
	47.879320282954723, 37.787089704535902, 39.556233814917505, 
	47.577056798618287, 6.8728428245522082, 43.875519670546055, 
	39.579561182763427, 25.745320606511086, 25.036056828685105, 
	23.746478438843042, 27.632271365728229, 26.058300752192736, 
	15.194642627611756, 53.771233516745269, 34.867460625246167, 
	13.536488504614681, 39.008358605206013, 28.271605169866234, 
	35.431443946436048, 29.978351634927094, 3.6089984234422445, 
	21.293389347847551, 13.332194199319929, 28.297816063277423, 
	15.396978685166687, 40.975269672926515, 13.331877179443836, 
	17.24086266476661, 43.361083446070552, 48.23409007396549, 
	39.710249588824809, 29.063765536993742, 30.859019320923835, 
	17.046077526174486, 5.1564661855809391, 4.6858567232266068, 
	44.873376556206495, 7.1859744703397155, 25.148897031322122, 
	10.100516346283257, 10.082306293770671, 43.258491009473801, 
	16.213865540921688, 3.1564271403476596, 52.181077380664647, 
	26.546401157975197, 39.604733693413436, 25.173957671504468, 
	21.093498691916466, 39.344903114251792, 26.071346183307469, 
	20.556752567645162, 50.526087288744748, 39.281111070886254, 
	40.226274591404945, 32.617274401709437, 13.937692253384739, 
	40.663796833250672, 33.290148181840777, 35.736465957015753, 
	41.206197724677622, 31.900172824971378, 44.242134722881019, 
	8.3278567809611559, 33.544870712794363, 20.461637126281857, 
	37.158720337320119, 21.936044464819133, 40.823473767377436, 
	29.747407189570367, 35.784880397841334, 7.4870092538185418, 
	22.582909148186445, 38.386246252339333, 37.218625221867114, 
	15.647187186405063, 2.8034813934937119, 12.668852359056473, 
	24.110773066058755, 39.500303200911731, 25.428297156468034, 
	29.72847884055227, 34.310592203401029, 23.291263310238719, 
	4.8293538112193346, 36.845632763579488, 36.186114435549825, 
	4.8897737683728337, 30.895574884489179, 13.271485648583621, 
	37.945418681483716, 5.2308572805486619, 20.897364090196788, 
	12.234537000767887, 45.924350158311427, 44.126807725988328, 
	38.357055785600096, 11.317996387369931, 46.76684933481738, 
	17.525049480609596, 39.990928792394698, 22.709356704726815, 
	30.739248702302575, 25.264597423374653, 36.486726445145905, 
	13.678435042966157, 24.921836494468153, 41.774840021971613, 
	26.24315662542358, 38.628265699371696, 45.130263301543891, 
	30.252290440257639, 32.938856552354991, 21.910963787231594, 
	1.01112250238657, 20.635680782143027, 45.328407310880721, 
	44.914956837892532, 42.700277089606971, 22.831218251958489, 
	44.311224529519677, 14.818650826346129, 43.350918428041041, 
	20.398069655057043, 49.473441110458225, 12.468474882189184, 
	20.546728358604014, 35.984048373065889, 29.042726173065603, 
	26.340818363241851, 28.577468921430409, 47.886289609596133, 
	2.7299058786593378, 39.001392740756273, 13.123536501079798, 
	31.271212168503553, 29.419533652253449, 17.960723571013659, 
	11.491524616722018, 37.142846223432571, 47.128635873086751, 
	49.953823303803802, 46.436614512931556, 38.59219811623916, 
	36.842365437187254, 26.409177205059677, 37.328651405405253, 
	34.927986776456237, 23.656566813588142, 54.236644501797855, 
	52.210703836753964, 28.157244871836156, 40.096746440976858, 
	31.42181457253173, 39.851036574691534, 31.424492257647216, 
	21.042729932814837, 29.523033264558762, 15.452032738830894, 
	30.015379111282527, 17.606630416121334, 57.101884072180837, 
	20.443842625245452, 28.15221022348851, 49.629977764561772, 
	51.881898459978402, 57.992464969865978, 41.258769547566772, 
	43.240769172552973, 29.272754858247936, 22.174956740345806, 
	15.878613241948187, 47.602522533852607, 15.149760278873146, 
	29.244020255282521, 10.318592390976846, 15.792438546195626, 
	58.153129788115621, 20.091362772509456, 17.016000929288566, 
	55.212599416263402, 34.973146319389343, 53.71798048261553, 
	42.072467103134841, 33.85334680788219, 42.941173664294183, 
	30.652374024502933, 24.599534531589597, 58.054738822393119, 
	49.216255359351635, 58.694438498932868, 33.398411432281137, 
	33.487216115463525, 56.60395851591602, 40.388811649754643, 
	45.681617558002472, 44.670640942640603, 47.806818750686944, 
	63.221794744022191, 21.259777629747987, 46.282551544718444, 
	28.109394870698452, 39.038281014654785, 26.695298594422638, 
	52.581904833205044, 38.105826680548489, 53.261818196624517, 
	22.38856595242396, 34.513444658368826, 52.018882452975959, 
	46.146133805159479, 27.751676877960563, 22.251768377609551, 
	23.039224445819855, 27.827243870124221, 47.182311916258186, 
	33.193872570991516, 33.675142102874815, 47.399089089594781, 
	29.215442202985287, 19.19504533521831, 39.283584766089916, 
	45.919479306321591, 24.246974098496139, 40.75248803012073, 
	18.43075955985114, 38.122169103007764, 18.892768386285752, 
	33.499488248489797, 28.507504318840802, 46.171165551058948, 
	55.270138657651842, 53.00451789284125, 26.972963795997202, 
	59.249994799029082, 21.945157558657229, 56.327855265699327, 
	40.006134547293186, 37.262823963537812, 25.617648903280497, 
	55.879329103045166, 24.799662989098579, 40.955960997380316, 
	43.778304925654083, 29.831844971049577, 56.157954586669803, 
	64.189819009043276, 39.882790127303451, 43.470695824362338, 
	41.053441364783794, 19.858526606112719, 27.550494831521064, 
	60.940338368527591, 55.592591632157564, 47.432575717102736, 
	38.313838448375463, 63.655514353886247, 23.842276630457491, 
	55.173727651126683, 36.56159303849563, 51.171865740325302, 
	14.076206854078919, 28.135600308887661, 40.992705882526934, 
	33.077532607130706, 27.613630839623511, 31.261784746311605, 
	54.491466227918863, 17.863415980245918, 45.307226125150919, 
	14.00991789996624, 47.494816437829286, 30.550595805980265, 
	22.459614060353488, 25.26177371153608, 38.942916460800916, 
	53.661402179859579, 1.1891810544766486, 13.974514538887888, 
	112.22612937679514, 108.91524361725897, 99.350027310661972, 
	16.902064719237387, 116.52611221047118, 118.44030319387093, 
	1.8005443713627756, 105.51905111595988, 10.464443478733301, 
	105.72949772095308, 70.603553343098611, 4.6445390335284173, 
	59.767248257063329, 115.36569624301046, 99.718200305011123, 
	96.584081499371678, 67.470765164587647, 106.02789150504395, 
	117.65104419644922, 120.4371798792854, 75.325613526627421, 
	10.272050822619349, 86.328772985842079, 22.265794642269611, 
	109.86511220689863, 47.378903591539711, 121.86315448442474, 
	50.738730230834335, 29.86259401217103, 123.6715569444932, 
	20.12369928508997, 89.050810230895877, 9.4124473645351827, 
	110.91702653188258, 50.259674362838268, 112.20425055501983, 
	14.098267829511315, 95.70136451581493, 106.60690627107397, 
	16.091330937575549, 56.603717565070838, 123.85576429218054, 
	13.681861021090299, 79.74437233665958, 2.6124773141928017, 
	92.783686938229948, 124.67039604252204, 10.88409340986982, 
	106.01856771204621, 97.532898479606956, 106.57870614249259, 
	3.7366205044090748, 67.274674311280251, 121.8603552589193, 
	108.17202668404207, 85.303708528168499, 98.969699303619564, 
	108.43108559260145, 55.16086700046435, 111.21830724459141, 
	27.87293240102008, 113.35955679602921, 63.842655413784087, 
	41.029215455055237, 15.900279967579991, 81.035312949214131, 
	108.59107061801478, 30.689855365548283, 89.412319208029658, 
	104.32100826594979, 73.73577278200537, 62.240019062999636, 
	122.7273592511192, 127.69526018900797, 53.266940758563578, 
	26.385443199425936, 31.852994742337614, 60.318720106035471, 
	106.75023480178788, 123.1608775104396, 67.136526592075825, 
	74.612110088579357, 50.014800885692239, 7.3641170430928469, 
	34.079227840993553, 22.432145269587636, 56.670952616259456, 
	51.80826749978587, 22.195442789234221, 30.981756458058953, 
	75.236678312532604, 9.5786990523338318, 123.50675575854257, 
	107.84753859229386, 26.502576827071607, 123.75724596064538, 
	59.773312723264098, 25.075215117074549, 8.0129718733951449, 
	75.072658884804696, 7.5858649895526469, 88.483324609696865, 
	62.445394677575678, 98.868400434497744, 34.609224846586585, 
	77.441016503609717, 94.043007473926991, 101.08413251955062, 
	124.31909865280613, 96.04773828946054, 82.360889608506113, 
	18.891335835680366, 68.46037698443979, 60.479840363375843, 
	67.249214683193713, 1.8475548955611885, 97.836682190652937, 
	122.86259609041736, 20.111821161117405, 76.971484495792538, 
	33.145381717011333, 43.719234322197735, 86.634274802636355, 
	84.64367246767506, 61.662859959993511, 70.248976637609303, 
	107.94150966824964, 9.0721677448600531, 13.673135004937649, 
	57.889076344668865, 88.426921632140875, 102.89679186278954, 
	122.23789864778519, 110.53688008151948, 93.966986203100532, 
	92.076594402547926, 27.435205929912627, 37.485007658600807, 
	109.78439133847132, 103.388507258147, 104.15006951801479, 
	121.72934733005241, 120.76649495400488, 99.527349165640771, 
	28.798621307592839, 48.62988393381238, 29.444443963933736, 
	97.544316857121885, 97.303895185235888, 78.383364802692086, 
	28.332693761680275, 64.851169708650559, 31.460166114382446, 
	94.253272936213762, 54.216616495978087, 103.43817389151081, 
	67.886691409628838, 61.006220073904842, 108.77033335529268, 
	107.87698537064716, 59.092851595487446, 57.969143963884562, 
	112.69650015281513, 125.05267392098904, 27.910256846807897, 
	11.642416686750948, 90.939782814122736, 45.338254536967725, 
	109.72203395236284, 16.452337529044598, 119.88565072463825, 
	79.953929136507213, 18.934011354111135, 19.302236509043723, 
	93.045928588602692, 41.733016965445131, 59.477046996355057, 
	94.16006256500259, 66.712428285740316, 70.038897995371372, 
	34.763222227338701, 6.9380326722748578, 59.907406228594482, 
	6.6968645481392741, 76.787452613469213, 58.724504555109888, 
	116.29605393577367, 86.348621093202382, 67.041910664644092, 
	33.243913532234728, 113.08473249897361, 83.9518963964656, 
	50.902247272897512, 49.836591531988233, 49.785904100630432, 
	4.0624160408042371, 59.907728772610426, 2.9345146864652634, 
	112.50538549013436, 63.047596833668649, 95.98052541539073, 
	30.548221048433334, 64.676462328992784, 35.60327540198341, 
	3.2184868971817195, 18.215725228190422, 93.682220642920583, 
	26.527585106436163, 99.641311642713845, 66.776509573217481, 
	94.776973059400916, 51.949432726483792, 58.044449002481997, 
	87.996135239023715, 103.49550795322284, 9.1604594085365534, 
	6.4134506094269454, 76.917759567964822, 50.933169982861727, 
	81.086652269586921, 78.237770894076675, 104.57264517294243, 
	35.713973348028958, 106.51287221210077, 111.55305506289005, 
	19.486252128146589, 28.542523724492639, 90.392913988791406, 
	105.26240702951327, 108.32411298854277, 122.32006587833166, 
	1.900063393637538, 56.434400526341051, 108.7762912530452, 
	64.264323640149087, 84.198362663853914, 115.90487617254257, 
	61.198874360881746, 115.97609204612672, 121.22355076530948, 
	30.060146544594318, 28.798388950526714, 10.580789789091796, 
	65.893074819818139, 31.833009195979685, 110.60343799134716, 
	8.0986612443812191, 3.6456328718923032, 21.536151939071715, 
	78.216087408363819, 96.394372925162315, 7.4307776181958616, 
	88.162504672538489, 23.735272297635674, 106.9431948964484, 
	51.052008883096278, 8.2300993502140045, 35.497167291585356, 
	36.840035066008568, 71.918207032140344, 9.5465767672285438, 
	26.514410655945539, 74.923283629585057, 119.80558476131409, 
	119.06185893714428, 118.85850019846112, 94.20657235942781, 
	90.793678595218807, 103.8247912703082, 84.981947926338762, 
	107.10667231539264, 52.06488355435431, 106.47066353680566, 
	127.79056827910244, 87.410060520283878, 57.615856547374278, 
	49.256717476528138, 98.481984601356089, 111.88392525631934, 
	90.868358179926872, 102.14387483522296, 56.013101504649967, 
	70.860478753689677, 3.7261798125691712, 96.49347253376618, 
	86.058997620362788, 18.674823581706733, 110.45522812893614, 
	46.838189892470837, 45.84736270736903, 53.604904909618199, 
	17.984899566043168, 119.51147867180407, 22.409610714763403, 
	31.344530477188528, 42.900306420400739, 106.08571967296302, 
	85.72173237381503, 69.538689414039254, 74.25811324082315, 
	7.9587333234958351, 85.125704918988049, 42.732920446898788, 
	105.65350128849968, 3.4852933674119413, 7.572783209849149, 
	1.6293824329040945, 40.280327862128615, 81.089571783784777, 
	17.1322060842067, 79.84219974745065, 106.98096791142598, 
	15.674030893947929, 89.798768336419016, 4.9515464697033167, 
	46.41052556456998, 75.937270705122501, 72.20791667047888, 
	90.773783176671714, 101.54822202585638, 95.53999622957781, 
	75.273653720039874, 39.987971813417971, 78.374638195149601, 
	73.607084757182747, 10.151776146143675, 127.12993296328932, 
	68.22310304408893, 90.020829128101468, 46.924328956287354, 
	4.9403550084680319, 112.23323167301714, 99.919972344767302, 
	82.608843398746103, 15.215564381331205, 20.396968877408653, 
	86.134836092591286, 74.460751311853528, 69.427891523111612, 
	101.64312560530379, 125.21133517939597, 78.160739174578339, 
	102.03863861085847, 91.572874325793236, 62.886631165631115, 
	88.28104078117758, 124.19953790213913, 47.304079293739051, 
	63.565674199722707, 110.08833324629813, 126.36571483081207, 
	66.721421373542398, 115.55380753241479, 95.336104530375451, 
	41.631134749390185, 12.820422243792564, 9.8095324616879225, 
	88.936884183436632, 39.278557303827256, 14.152523647993803, 
	54.676685379352421, 32.540594122838229, 41.9198366003111, 
	52.253170723095536, 103.38720667362213, 104.76547622634098, 
	52.902209120802581, 63.046359764412045, 76.231106647755951, 
	6.075936506036669, 63.996591053437442, 68.716641753446311, 
	23.302505458239466, 86.276771190576255, 123.98339603561908, 
	104.80975831439719, 8.1545471106655896, 66.004489412065595, 
	57.513177492655814, 6.4055330231785774, 121.82742046331987, 
	29.742494091391563, 70.17482131998986, 84.239852977916598, 
	56.553547718096524, 108.75535717653111, 21.586255491245538, 
	95.802963929250836, 89.3324133339338, 80.385567367542535, 
	73.94469679472968, 2.2613418470136821, 64.620563392993063, 
	13.214770566206425, 120.59811645094305, 112.87877164501697, 
	51.176443573087454, 91.807992177084088, 38.261226983275265, 
	110.11891567055136, 74.61393428966403, 35.157809175550938, 
	66.479859383776784, 52.404002275783569, 104.76833258010447, 
	59.61862442176789, 7.7243995368480682, 13.236612603534013, 
	70.044677293393761, 49.177383476402611, 23.774166173767298, 
	111.73661134950817, 43.404185994993895, 101.32801293488592, 
	124.99400278599933, 44.028677515685558, 35.33714713761583, 
	110.81632118159905, 52.04907239228487, 12.83269949676469, 
	84.27451155660674, 57.374981930013746, 23.296429105103016, 
	46.875313147436827, 37.285776998382062, 109.23717672750354, 
	121.88945010118186, 56.522361129987985, 111.71540050255135, 
	106.97861281968653, 4.5539058209396899, 7.5747934030368924, 
	72.112247655168176, 41.279151524882764, 99.756925043649971, 
	81.794931074138731, 91.569000899791718, 112.73338886257261, 
	60.91184998722747, 102.51552472868934, 15.251541463658214, 
	74.753272654023021, 114.06800881866366, 30.259589750319719, 
	63.520985472016037, 80.653930789791048, 82.853134238626808, 
	31.650747284293175, 19.569080473389477, 38.008373587392271, 
	123.8591841221787, 76.673401191364974, 73.833030211273581, 
	118.86042168317363, 48.145031321793795, 96.870499875862151, 
	67.528523016255349, 97.816512071527541, 64.173833011649549, 
	18.983323212713003, 100.66148648923263, 56.38826963538304, 
	77.177857962436974, 124.0922244801186, 47.958741453941911, 
	117.78981328895316, 82.253885779064149, 63.196024090982974, 
	102.14547560922801, 110.68913197517395, 51.206975551787764, 
	20.13991430727765, 121.18264113692567, 93.827585033141077, 
	76.760984015185386, 125.17565302317962, 87.219013394322246, 
	91.531667583156377, 82.462917720433325, 106.2987010339275, 
	58.164641830604523, 16.949706374667585, 84.825106670148671, 
	28.579223012086004, 123.9073228999041, 26.653621753212065, 
	125.34306519525126, 46.765706611331552, 109.67028876114637, 
	6.0192889296449721, 65.203772591426969, 35.70885327225551, 
	86.050833424553275, 67.4117102175951, 103.95022966712713, 
	119.03416379634291, 96.946273472625762, 4.8419269821606576, 
	40.486013010144234, 18.979430566541851, 110.72478745970875, 
	48.207682451233268, 126.16710424283519, 123.52224721526727, 
	14.536492428276688, 112.05245101545006, 46.1478112754412, 
	47.268193403724581, 84.305266844108701, 71.446274513844401, 
	127.62415939988568, 45.672099894378334, 104.23190941428766, 
	8.3533351155929267, 8.5989660448394716, 42.405761461239308, 
	44.323119573760778, 74.859891015570611, 121.19334972836077, 
	43.31347992317751, 68.015393213834614, 86.78873009653762, 
	5.6137824412435293, 127.02277342090383, 50.204687763936818, 
	17.314854410011321, 108.27879833383486, 100.80081604188308, 
	53.53342746431008, 91.133807824458927, 75.540389578323811, 
	116.48509985767305, 89.081820703111589, 14.043018062133342, 
	4.5370793016627431, 84.958107168786228, 97.215218290220946, 
	27.028158693574369, 86.199492447543889, 60.526595643721521, 
	101.28304495289922, 104.57364846579731, 11.916640723124146, 
	75.952226776629686, 39.546004130970687, 71.562157411128283, 
	43.458098393399268, 57.30590310273692, 68.217438475694507, 
	66.678635265212506, 39.758875848259777, 118.08292152639478, 
	126.47749726753682, 55.257680118549615, 99.708781002555043, 
	94.731066541746259, 72.071128615643829, 96.452631861437112, 
	56.997833293862641, 16.081848888657987, 80.914266804698855, 
	8.8807347370311618, 118.66675274865702, 126.39894271921366, 
	95.722731075715274, 16.864393482450396, 54.309792152605951, 
	1.7568814684636891, 49.853363643400371, 8.5705695748329163, 
	113.50561654241756, 65.635460913181305, 75.626774719450623, 
	123.10319819534197, 121.87661073217168, 54.48451384389773, 
	123.58858524635434, 24.629198247566819, 41.741277498658746, 
	98.80146994907409, 7.2318868297152221, 23.423305643722415, 
	124.61030433000997, 86.878359070047736, 43.310668869875371, 
	4.3746497663669288, 99.028860996477306, 124.79221531096846, 
	10.745748994406313, 104.13702765712515, 110.05223049363121, 
	18.314868777059019, 92.362476802431047, 93.092822609003633, 
	25.012198331300169, 84.278742359485477, 75.020786012522876, 
	112.06068469956517, 22.832379458472133, 25.181845897808671, 
	85.522897235117853, 1.4414876513183117, 97.95209037931636, 
	59.588973318226635, 43.177464341744781, 39.950869020540267, 
	24.973086350597441, 75.107367790769786, 5.5766812451183796, 
	96.805069095920771, 118.05704615125433, 15.807796828448772, 
	21.886923110112548, 109.03895024722442, 68.248912952374667, 
	112.43861042195931, 64.544909789226949, 58.437038531526923, 
	6.2943387827835977, 3.1839116453193128, 64.848130733706057, 
	104.44216665625572, 105.33919474156573, 83.255195940379053, 
	90.763933756388724, 97.233933945186436, 19.554570960346609, 
	117.39853410236537, 91.656395957339555, 66.455358043778688, 
	60.547145375981927, 97.822906001005322, 45.296213698107749, 
	42.66535137873143, 37.171961008571088, 14.104529642499983, 
	83.472290003672242, 38.700142088811845, 80.29359867842868, 
	20.836555565241724, 84.154159703757614, 39.717092382255942, 
	78.037606784142554, 66.561022730078548, 64.273198627401143, 
	121.12763413507491, 39.029513168148696, 107.32401417102665, 
	47.943311383482069, 107.8816796457395, 1.7632742743007839, 
	102.35621690377593, 64.919669387862086, 103.92418284807354, 
	7.939144067466259, 6.097014881670475, 98.069574464578182, 
	43.672605658415705, 76.44370283652097, 93.10445033852011, 
	119.15660763159394, 50.602636467665434, 71.490312298759818, 
	112.6646269611083, 122.59639674099162, 9.7024168004281819, 
	101.21953997481614, 89.578021908644587, 69.700624902267009, 
	107.69414502056316, 110.16764246718958, 5.5434311795979738, 
	120.86370660969988, 31.780507617164403, 98.370596858207136, 
	63.076896888203919, 109.84483699779958, 99.538263856433332, 
	28.709555869456381, 7.4020868726074696, 84.07446667086333, 
	94.466958997305483, 61.367887306958437, 15.723022455815226, 
	68.929464148823172, 58.658576830290258, 85.163864111527801, 
	82.912024306133389, 63.12157680420205, 44.757207164075226, 
	9.7332256678491831, 52.734322240110487, 6.2348521803505719, 
	71.943686944898218, 8.7115876046009362, 114.33903189841658, 
	1.5682977177202702, 66.358886867295951, 55.296287052333355, 
	82.870324817486107, 86.492600330151618, 94.927373252343386, 
	8.8380059343762696, 121.20886661484838, 46.79184781620279, 
	84.711086887400597, 117.10366559308022, 123.6823690882884, 
	98.274647405836731, 59.672629254404455, 1.2354587875306606, 
	120.18848010990769, 6.4091389044187963, 99.734471154399216, 
	110.07859435677528, 59.816529481671751, 32.157749338075519, 
	117.33822038210928, 65.888795936945826, 5.9734674538485706, 
	77.262249410152435, 110.40264801355079, 62.366291496902704, 
	127.85841086506844, 6.2249982063658535, 65.05058440938592, 
	3.6131852040998638, 62.87757172388956, 52.612359500490129, 
	49.526569426991045, 122.05007542436942, 36.209037665743381, 
	76.259744937065989, 24.717098647728562, 84.960403476376086, 
	72.866034144069999, 79.582048246171325, 87.848932099528611, 
	27.529080491047353, 78.354156295303255, 101.83271454134956, 
	27.593141198158264, 28.380562806501985, 99.836039224173874, 
	40.271735736634582, 76.794884727802128, 16.596009006723762, 
	79.151404023170471, 7.7159544308669865, 120.52306316699833, 
	104.49418163951486, 7.3184443609789014, 98.028348383493721)
)
, names = c("class", "x", "y")
, row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
	"15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26",
	"27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
	"39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",
	"51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62",
	"63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74",
	"75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86",
	"87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98",
	"99", "100", "101", "102", "103", "104", "105", "106", "107", "108",
	"109", "110", "111", "112", "113", "114", "115", "116", "117", "118",
	"119", "120", "121", "122", "123", "124", "125", "126", "127", "128",
	"129", "130", "131", "132", "133", "134", "135", "136", "137", "138",
	"139", "140", "141", "142", "143", "144", "145", "146", "147", "148",
	"149", "150", "151", "152", "153", "154", "155", "156", "157", "158",
	"159", "160", "161", "162", "163", "164", "165", "166", "167", "168",
	"169", "170", "171", "172", "173", "174", "175", "176", "177", "178",
	"179", "180", "181", "182", "183", "184", "185", "186", "187", "188",
	"189", "190", "191", "192", "193", "194", "195", "196", "197", "198",
	"199", "200", "201", "202", "203", "204", "205", "206", "207", "208",
	"209", "210", "211", "212", "213", "214", "215", "216", "217", "218",
	"219", "220", "221", "222", "223", "224", "225", "226", "227", "228",
	"229", "230", "231", "232", "233", "234", "235", "236", "237", "238",
	"239", "240", "241", "242", "243", "244", "245", "246", "247", "248",
	"249", "250", "251", "252", "253", "254", "255", "256", "257", "258",
	"259", "260", "261", "262", "263", "264", "265", "266", "267", "268",
	"269", "270", "271", "272", "273", "274", "275", "276", "277", "278",
	"279", "280", "281", "282", "283", "284", "285", "286", "287", "288",
	"289", "290", "291", "292", "293", "294", "295", "296", "297", "298",
	"299", "300", "301", "302", "303", "304", "305", "306", "307", "308",
	"309", "310", "311", "312", "313", "314", "315", "316", "317", "318",
	"319", "320", "321", "322", "323", "324", "325", "326", "327", "328",
	"329", "330", "331", "332", "333", "334", "335", "336", "337", "338",
	"339", "340", "341", "342", "343", "344", "345", "346", "347", "348",
	"349", "350", "351", "352", "353", "354", "355", "356", "357", "358",
	"359", "360", "361", "362", "363", "364", "365", "366", "367", "368",
	"369", "370", "371", "372", "373", "374", "375", "376", "377", "378",
	"379", "380", "381", "382", "383", "384", "385", "386", "387", "388",
	"389", "390", "391", "392", "393", "394", "395", "396", "397", "398",
	"399", "400", "401", "402", "403", "404", "405", "406", "407", "408",
	"409", "410", "411", "412", "413", "414", "415", "416", "417", "418",
	"419", "420", "421", "422", "423", "424", "425", "426", "427", "428",
	"429", "430", "431", "432", "433", "434", "435", "436", "437", "438",
	"439", "440", "441", "442", "443", "444", "445", "446", "447", "448",
	"449", "450", "451", "452", "453", "454", "455", "456", "457", "458",
	"459", "460", "461", "462", "463", "464", "465", "466", "467", "468",
	"469", "470", "471", "472", "473", "474", "475", "476", "477", "478",
	"479", "480", "481", "482", "483", "484", "485", "486", "487", "488",
	"489", "490", "491", "492", "493", "494", "495", "496", "497", "498",
	"499", "500", "501", "502", "503", "504", "505", "506", "507", "508",
	"509", "510", "511", "512", "513", "514", "515", "516", "517", "518",
	"519", "520", "521", "522", "523", "524", "525", "526", "527", "528",
	"529", "530", "531", "532", "533", "534", "535", "536", "537", "538",
	"539", "540", "541", "542", "543", "544", "545", "546", "547", "548",
	"549", "550", "551", "552", "553", "554", "555", "556", "557", "558",
	"559", "560", "561", "562", "563", "564", "565", "566", "567", "568",
	"569", "570", "571", "572", "573", "574", "575", "576", "577", "578",
	"579", "580", "581", "582", "583", "584", "585", "586", "587", "588",
	"589", "590", "591", "592", "593", "594", "595", "596", "597", "598",
	"599", "600", "601", "602", "603", "604", "605", "606", "607", "608",
	"609", "610", "611", "612", "613", "614", "615", "616", "617", "618",
	"619", "620", "621", "622", "623", "624", "625", "626", "627", "628",
	"629", "630", "631", "632", "633", "634", "635", "636", "637", "638",
	"639", "640", "641", "642", "643", "644", "645", "646", "647", "648",
	"649", "650", "651", "652", "653", "654", "655", "656", "657", "658",
	"659", "660", "661", "662", "663", "664", "665", "666", "667", "668",
	"669", "670", "671", "672", "673", "674", "675", "676", "677", "678",
	"679", "680", "681", "682", "683", "684", "685", "686", "687", "688",
	"689", "690", "691", "692", "693", "694", "695", "696", "697", "698",
	"699", "700", "701", "702", "703", "704", "705", "706", "707", "708",
	"709", "710", "711", "712", "713", "714", "715", "716", "717", "718",
	"719", "720", "721", "722", "723", "724", "725", "726", "727", "728",
	"729", "730", "731", "732", "733", "734", "735", "736", "737", "738",
	"739", "740", "741", "742", "743", "744", "745", "746", "747", "748",
	"749", "750", "751", "752", "753", "754", "755", "756", "757", "758",
	"759", "760", "761", "762", "763", "764", "765", "766", "767", "768",
	"769", "770", "771", "772", "773", "774", "775", "776", "777", "778",
	"779", "780", "781", "782", "783", "784", "785", "786", "787", "788",
	"789", "790", "791", "792", "793", "794", "795", "796", "797", "798",
	"799", "800", "801", "802", "803", "804", "805", "806", "807", "808",
	"809", "810", "811", "812", "813", "814", "815", "816", "817", "818",
	"819", "820", "821", "822", "823", "824", "825", "826", "827", "828",
	"829", "830", "831", "832", "833", "834", "835", "836", "837", "838",
	"839", "840", "841", "842", "843", "844", "845", "846", "847", "848",
	"849", "850", "851", "852", "853", "854", "855", "856", "857", "858",
	"859", "860", "861", "862", "863", "864", "865", "866", "867", "868",
	"869", "870", "871", "872", "873", "874", "875", "876", "877", "878",
	"879", "880", "881", "882", "883", "884", "885", "886", "887", "888",
	"889", "890", "891", "892", "893", "894", "895", "896", "897", "898",
	"899", "900", "901", "902", "903", "904", "905", "906", "907", "908",
	"909", "910", "911", "912", "913", "914", "915", "916", "917", "918",
	"919", "920", "921", "922", "923", "924", "925", "926", "927", "928",
	"929", "930", "931", "932", "933", "934", "935", "936", "937", "938",
	"939", "940", "941", "942", "943", "944", "945", "946", "947", "948",
	"949", "950", "951", "952", "953", "954", "955", "956", "957", "958",
	"959", "960", "961", "962", "963", "964", "965", "966", "967", "968",
	"969", "970", "971", "972", "973", "974", "975", "976", "977", "978",
	"979", "980", "981", "982", "983", "984", "985", "986", "987", "988",
	"989", "990", "991", "992", "993", "994", "995", "996", "997", "998",
	"999", "1000", "1001", "1002", "1003", "1004", "1005", "1006", "1007",
	"1008", "1009", "1010", "1011", "1012", "1013", "1014", "1015", "1016",
	"1017", "1018", "1019", "1020", "1021", "1022", "1023", "1024", "1025",
	"1026", "1027", "1028", "1029", "1030", "1031", "1032", "1033", "1034",
	"1035", "1036", "1037", "1038", "1039", "1040", "1041", "1042", "1043",
	"1044", "1045", "1046", "1047", "1048", "1049", "1050", "1051", "1052",
	"1053", "1054", "1055", "1056", "1057", "1058", "1059", "1060", "1061",
	"1062", "1063", "1064", "1065", "1066", "1067", "1068", "1069", "1070",
	"1071", "1072", "1073", "1074", "1075", "1076", "1077", "1078", "1079",
	"1080", "1081", "1082", "1083", "1084", "1085", "1086", "1087", "1088",
	"1089", "1090", "1091", "1092", "1093", "1094", "1095", "1096", "1097",
	"1098", "1099", "1100", "1101", "1102", "1103", "1104")
, class = "data.frame"
)

"cross" <-
structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
1.26295428488079, -0.326233360705649, 1.3297992629225, 1.2724293214294, 
0.414641434456408, -1.53995004190371, -0.928567034713538, -0.29472044679056, 
-0.00576717274753695, 2.40465338885795, 0.76359346114046, -0.799009248989368, 
-1.14765700923635, -0.289461573688223, -0.299215117897316, -0.411510832795067, 
0.252223448156132, -0.891921127284569, 0.435683299355719, -1.23753842192996, 
-0.224267885278309, 0.377395645981701, 0.133336360814841, 0.804189509744908, 
-0.0571067743838087, 0.503607972233726, 1.08576936214569, -0.69095383969683, 
-1.28459935387219, 0.046726172188352, -0.235706556439501, -0.542888255010254, 
-0.433310317456782, -0.649471646796233, 0.726750747385451, 1.1519117540872, 
0.992160365445798, -0.429513109491881, 1.23830410085338, -0.279346281854269, 
1.75790308981071, 0.560746090888056, -0.452783972553158, -0.832043296117832, 
-1.16657054708471, -1.0655905803883, -1.563782051071, 1.15653699715018, 
0.83204712857239, -0.227328691424755, 0.266137361672105, -0.376702718583628, 
2.44136462889459, -0.795339117255372, -0.0548774737115786, 0.250141322854153, 
0.618243293566247, -0.172623502645857, -2.22390027400994, -1.26361438497058, 
0.358728895971352, -0.0110454784656636, -0.940649162618608, -0.115825322156954, 
-0.814968708869917, 0.242263480859686, -1.4250983947325, 0.36594112304922, 
0.248412648872596, 0.0652881816716207, 0.0191563916602738, 0.257338377155533, 
-0.649010077708898, -0.119168762418038, 0.66413569989411, 1.10096910219409, 
0.14377148075807, -0.117753598165951, -0.912068366948338, -1.43758624082998, 
-0.797089525071965, 1.25408310644997, 0.77214218580453, -0.21951562675344, 
-0.424810283377287, -0.418980099421959, 0.996986860909106, -0.275778029088027, 
1.2560188173061, 0.646674390495345, 1.29931230256343, -0.873262111744435, 
0.00837095999603331, -0.880871723252545, 0.59625901661066, 0.119717641289536, 
-0.282173877322452, 1.45598840106634, 0.229019590694692, 0.996543928544126, 
0.781859184600258, -0.776776621764597, -0.615989907707918, 0.0465803028049967, 
-1.13038577760069, 0.576718781896486, -1.28074943178832, 1.62544730346494, 
-0.500696596002705, 1.67829720781629, -0.412519887482398, -0.97228683550556, 
0.0253828675878054, 0.0274753367451927, -1.68018272239593, 1.05375086302862, 
-1.11959910457218, 0.335617209968815, 0.494795767113158, 0.138052708711737, 
-0.118792025778828, 0.197684262345795, -1.06869271125479, -0.80321321736474, 
-1.11376513631953, 1.58009168370384, 1.49781876103841, 0.262645458662762, 
-1.23290119957126, -0.00372353379218051, 1.51167228281089, -0.475698284429534, 
0.79791643753108, -0.974002561112527, 0.689372697765473, -0.955839103276798, 
-1.2317070584141, -0.956891881325619, -0.869782873686493, -0.910680682493289, 
0.741276305260208, 0.0685115332771444, -0.323750754587962, -1.0865030469937, 
-1.01592894685489, -0.767790184730859, -1.11972006112698, -0.448174236603396, 
0.47173637445323, -1.18049068288428, 1.47025699708299, -1.31142059162488, 
-0.0965249227629818, 2.36971990795134, 0.890626476437422, -0.252183161390947, 
-0.865763754844263, 0.582585999886443, -0.0125293469957165, -0.374854762113895, 
0.317885735092473, -0.488805634836622, 2.65865802699847, 1.68027820468955, 
0.779584009075853, 0.71324052011282, -0.542881937111931, 0.885778373900932, 
-0.348594684510753, -1.00805457816954, 1.8831825423826, -0.928971079201175, 
-0.294196453606913, -0.614950270796355, -0.94707579174102, 0.598975149810313, 
-1.52361488237755, -0.206189002117285, -0.574295414301406, -1.39016603663508, 
-0.0704173825071551, -0.430879529604084, -0.592225372961588, 
0.981116159906615, 0.532409356670506, -0.0904561242409946, 0.15649049169609, 
-0.737311690839969, -0.20134120587218, 1.10217659507369, -0.0167482561518744, 
0.161788633557242, 2.02476139007773, -0.703694253898738, 0.960792383950527, 
1.79048505353782, -1.06416516333343, 0.0176365464246015, -0.38990862930585, 
-0.490832752366552, -1.04571765178674, -0.896211263960579, 1.26938716360704, 
0.593840948597503, 0.775634318883504, 1.55737037634497, -0.36540179688389, 
0.816556448748125, -0.0606347781205525, -0.501378317795089, 0.92606272533302, 
0.0369376908421676, -1.06620017351302, -0.238456353279062, 1.49522344384323, 
1.17215854698947, -1.45770721012238, 0.095056226988643, 0.847664963596026, 
-1.62436452976059, 1.40856335683396, -0.541760362101113, 0.278664724106355, 
-0.193972744670297, 1.57615818121873, -1.47554763526052, -0.144608207310541, 
-1.07501019081817, 0.406542731944927, 2.22926220164091, -1.51449700825033, 
-0.0617074220120377, -0.147270790038966, 1.54159306882696, -0.98185566880387, 
0.496578172661699, 1.69694788072309, -0.260736308568123, -0.705928585667662, 
-0.161178506172346, 0.501321827723727, -1.01353967049479, 1.61475223546813, 
0.0056419848524955, -2.90489906034557, -1.10716481896875, 1.54756693261827, 
-0.976830350346719, -0.101503447631726, 0.0426502497966979, -0.797552248911132, 
10.3229809229786, 0.781313159613491, -2.65504530767321, 4.61099644940903, 
2.54982271068308, 3.46089344603239, 9.04323905082022, -2.57369212146266, 
-2.13605901599485, -1.83306930522421, 9.6522480508404, 17.5525000708904, 
8.68008660159401, 9.34564129581753, 15.919862433606, -5.14485132304678, 
-13.2325460080511, -9.99329800111672, 2.28963050204598, 0.210102128452628, 
-24.4433276505523, 10.6676032628399, -1.97779677756853, -0.510054846089218, 
3.65908058495626, -11.7848695467608, -6.36022158189192, 9.3050250294412, 
17.7548487070343, 6.18772566143988, 6.64418294972592, -15.3744320863865, 
9.4735781165804, 10.1091519806059, -3.19132805175383, 1.69414967988762, 
8.31420024551634, 19.3743075878042, -9.98765408403015, 9.26557399868874, 
12.3943927477343, 8.23329946688387, 2.64261088308114, -1.36564762228545, 
1.48402295760307, 7.81719687498848, -9.70231044043853, -11.0009877142141, 
-6.40300226425223, -12.8162849967668, -15.0240997074182, 12.4131251912361, 
-8.27707125047001, -4.54041039418265, -10.2125863144471, 4.61495749483626, 
4.75189904427085, 10.694225915881, -4.15832721845055, 11.1696235926977, 
1.23319396655783, -11.6968071496312, 10.3626804826353, 13.6210273011827, 
7.21924088445048, -3.87483247589599, -12.1023084619283, 0.336978631624503, 
-12.1114961026169, -16.3260388434941, 6.33630178620148, -2.69117697432793, 
-3.59002166646955, -5.78985596921356, 11.8404454925165, -14.2647109326606, 
-4.00046782396941, 5.09564770116935, 2.61293107040096, -5.53989899086854, 
-0.558070368799121, -1.733492294586, 15.1483855371761, 0.923164614571009, 
-5.89351782073246, 10.6791767966101, -11.4513270461538, 4.84258122106827, 
11.0776535237012, -10.4913626711625, 2.9178850354192, 8.5578063748781, 
10.0569024168482, 13.7157214191727, 17.2624356871327, 1.80159575159047, 
-7.62471445285245, -6.81180785830493, 6.25594154731617, -14.1562608703647, 
4.50210419899281, -10.3942679356763, -8.32681105500008, 3.07714102841799, 
-5.30665166278555, -7.45594543269999, -5.97357188793219, 12.7479676735673, 
-12.3435210839711, 0.979746507959679, -4.20140027038469, 26.6556902696747, 
-1.42044462122646, 16.3748274602587, 7.773682959456, -14.097988978937, 
4.6391037193933, 6.26597560911691, 0.228789904986294, -17.5660497111041, 
20.2082286284336, -1.68963457912297, -12.6112607073466, 11.2352265517169, 
-2.22370407366639, 2.21003564042203, -7.0085120169909, -16.1179756273137, 
-1.50771470353655, 6.75881943156902, -0.153210630366342, 0.248210766257436, 
1.99931532387323, 3.58768642273758, -2.35366571473631, -2.72859861664959, 
-6.60214606316465, -12.4366924232844, -12.9751863593939, -7.0767914781292, 
-6.68234015626664, -3.51882222273542, -10.9417955632484, -8.01365032159383, 
13.4844536799115, 3.34655199853207, 2.34995475092581, -0.220894711088577, 
-8.25920319206087, -5.32695795960113, -3.33893750068784, 0.791318314544834, 
-0.312537102520872, 16.2573684211216, -3.06212458717065, -6.74066929990784, 
-14.7522315477412, -9.01068908154044, 23.3619465406393, 2.7287893364201, 
8.17872649370587, 1.870649530081, 1.60212601785311, -1.49188498353237, 
5.01393258023639, 12.9990096307516, 8.1122135461283, -1.99831539964344, 
0.955721840846537, -13.0429875011322, 10.2469655118127, 16.3244223972279, 
-13.599348594188, -0.162781728369867, -7.92169288525344, -10.7772938943597, 
9.61883062015695, 10.5006317263396, 18.26962646283, 4.50155466790179, 
-16.4047480086872, 4.40244954835033, -6.24500457366339, -1.79158060650603, 
2.17758668816648, 0.633172530255781, 7.1473145936452, -5.8984082653835, 
-6.09718587383515, 1.65436662902189, 7.89937440280074, -3.5300782137138, 
-6.02053474558027, -9.09119918804333, -6.75747536193682, -2.55197696992539, 
-10.492226832889, 6.22102487150906, -18.508838677455, -1.38757052528768, 
-13.004048795097, 10.9453577727738, -4.60054023288268, 10.1568488145141, 
-3.77678941449556, -1.91465761062821, -9.18594701960981, -27.741272558644, 
14.8152902674956, -12.7009260561152, 23.0862386776491, 2.93034342467055, 
-5.08984036846773, 6.29008808803855, 7.92643587614759, -2.04979473854538, 
1.22670673999928, 16.5906740964623, -12.6485254615721, -0.441747633791901, 
7.14739055980878, -12.726211639258, -16.1093462675347, 12.1895759920122, 
-6.72990034731278, -9.33707230085963, -0.325549223142774, -4.98915545796906, 
-5.02980076462877, -12.7788604477508, -8.17163868669447, 7.62869785065482, 
4.49516789673977, -5.74042965808372, 21.739546330377, 13.328173124688, 
-6.6508358280159, -0.879886195683914, -29.1274715911867, 8.61357292271955, 
5.03400757210134, -4.47930641568919, -1.74582278605548, 0.77701441181892, 
9.15756784351, -8.72926220574165, 9.59717161351004, 6.79958239263148, 
-6.32927064789725, -14.3704621286748, 4.41870635337353, 3.79443028846278, 
16.8651350873577, 9.31062891549901, 0.736292793186125, -0.742713858154477, 
5.45466087775891, -7.98678130785388, 0.948792511743886, 3.17587025986629, 
4.95354022609547, -10.2089787166516, 13.1611638487178, 6.31905039600816, 
22.5640003363503, -17.0102442926162, -5.30831511172054, -15.4305206716124, 
-3.78898108035027, 2.79127238854218, 15.3231352741305, -3.99046323927534, 
-10.7873737498628, -2.76642822802949, 5.58948782464363, 1.63711972263847, 
11.8656083807322, -2.69018382087596, -14.8339956794295, 8.56348622662952, 
-10.0181065779127, 5.55269828120174, 4.62144344836878, 3.32513193357923, 
15.5150471754534, -1.85530109975289, -11.8277562591, 0.571266866482798, 
-2.08779706049955, 5.71554294832304, 14.7117985693535, -16.2749482141941, 
-1.92497356346681, 0.633295271566325, 4.94736902627626, -6.27141191412903, 
3.51509348133188, 3.43270133189534, -0.11135492805367, -1.11991487778501, 
13.2007013386852, 6.06535808793815, 17.6078273692125, -2.42136911072255, 
-11.2009636651738, -3.56132630253873, 0.876569891304792, -2.14548255053791, 
-3.70645160894906, -14.1949624400534, -7.1754848582477, -9.86613104525906, 
2.77477873802221, 3.1031561129672, 13.8568327876692, -2.96562773149361, 
8.53550415667865, -4.31330028103257, -13.6339811550144, 3.91082990282695, 
-4.67583000313985, -7.51103129717156, -6.80982848418005, 9.80553145682072, 
14.1518962573123, 9.06627929084592, -2.45842193768092, -11.7871158752981, 
2.00316447773904, 10.0208986542001, 7.53654957630197, 2.83065091488878, 
1.99996264927645, -7.59253850928803, 3.99424793074109, 0.502173998364998, 
0.611751505697686, -1.81765568119765, -10.4221372863818, -5.33468943067439, 
6.89458845603604, 0.350353049869129, 0.131821312276605, -1.67685271246378, 
12.605317445574, 0.166370131024709, 2.24276410485421, 1.34298662106725, 
-8.66909862408359, -0.598200966697581, 11.5822983701919, 4.123127334642, 
-13.0682546113154, 0.696080492300982, 5.03905745679683, -0.67451912963356, 
7.04427927987277, -1.55401707551582, -9.46164400557883, 6.56506150914, 
2.36398718421689, 4.89292073238146, 9.36954268470611, 1.77755538331996, 
-14.6662045858619, 1.08936209143034, -14.7367975693665, -4.77938795401386, 
8.58311812774644, -15.4858559586947, 0.95688557873219, -5.47801646923058, 
-2.71077284357703, 8.78582272209678, 4.1040780243487, 11.6496726067414, 
-10.1988198552611, -7.82514315715003, -6.79473262674202, -1.16671815243822, 
-9.01621201929096, -7.37881451344166, -8.77097221842987, 5.43878467840335, 
4.939088293605, 8.24789391153992, 23.9540972998792, -1.62231358898104, 
6.16513289497398, 29.3977306698789, 5.04540413640492, -0.621155677483524, 
-8.75198642377827, -4.91927927170642, -15.1982309519111, -14.1513542914466, 
-3.64488442635934, 2.87357776159579, 0.363849134787787, -3.51008606457252, 
-16.373000070736, 5.93262639197573, 4.13659505405807, 14.5496370314655, 
-16.7057144527503, -2.5814149485095, 15.752896998019, 1.04772251438527, 
12.4582784242735, 5.16798821376334, 1.22841730825787, 8.22794390071317, 
-16.2074368615798, -3.05892575218397, 5.45638114974549, 12.070172795336, 
6.90558559055942, 1.74353100248949, 10.2651002148432, 0.124783240247643, 
-9.94775319135224, -0.226463751712447, -1.47306003874045, 3.33053773093185, 
-3.42742082496591, 5.87657132377338, 18.5520762546619, -16.1698044322373, 
5.25669412032774, -6.50477811548134, -5.66248198470931, -16.3458544926822, 
-2.33360188525512, 3.01169046129524, -12.8445085139894, 17.447654098285, 
-6.8357719187998, -20.5089852740129, -1.02645689822056, 21.1667099907255, 
14.3668872293677, 11.4970095993213, 7.10073772319519, 4.15318844050235, 
-3.94262220442507, -13.5702630834918, -20.006515182013, -10.6080519816946, 
-16.044842013448, -8.89127225331884, 6.55639819537952, -7.96216423084755, 
-13.8459727544267, -8.7500412233177, -15.4632257129114, 7.23381785292835, 
-13.516084529682, -1.30949164017766, 5.21512751076284, 10.8137301851799, 
17.0451970789251, -15.8420370914651, 8.3209212242194, -5.00887951195279, 
-1.62525724374587, 13.0269625677894, -5.46418307551073, 6.11426187254737, 
-0.842018763498, -4.4107766257393, 12.6959344496608, -2.02116410775661, 
-1.91245940578181, 6.26740626374344, 8.23664254694211, -8.31036884262024, 
10.3218593661557, -5.72278538833103, -7.97798977629545, -20.9982301623309, 
-1.30941724208386, 2.85267621427191, -6.36722989025543, 11.179400022582, 
5.58136975702211, 0.89912760830658, 16.2693313890745, -13.5218698096516, 
2.57076418690114, 7.61136267938843, -8.95808373839301, -2.31168712096488, 
-0.502704389686109, -4.00504715880425, 0.627255927160182, -1.39204535415135, 
-7.481371134681, 6.8538991711042, -5.18856239277115, -0.626368228189236, 
0.481335325559074, 1.69527108394262, -1.7612262939127, 0.198013015105549, 
0.397349098800811, 0.0292254949810759, 2.56027338872093, 1.25712771158455, 
-0.534537685615439, -0.625227428575037, 0.913848686849621, 1.00719953467783, 
0.719291823271818, -0.604711660854287, 0.53905440643653, -0.07683088553045, 
1.84991956032624, -0.854907551346389, 0.0326372948933622, -1.02505948075354, 
-0.982249075565237, 0.00410195695681333, -0.233427177912246, 
-0.498888219103613, 1.54971296288404, 0.087496916504864, 1.31870113056505, 
-0.981224119069448, -0.245622587644691, -1.40393383691901, 1.4408931473099, 
-0.981359990937249, 1.47424490392592, -0.99119724527207, -0.0944973445755386, 
-2.87514168382416, -0.246866100837098, 0.0147444897981717, -1.91908770329685, 
-0.287813743700098, -0.34663744540489, -1.83958858464688, 0.89858894092329, 
-1.21285501097807, -0.218964231724761, 0.564268160651169, -0.525434375699569, 
0.744374224954032, 0.128981753196475, 1.48827425685469, -0.662681951149737, 
-1.16065500147004, 0.358774234497368, -0.194846384192831, -0.295282008014163, 
0.496640423719942, 0.484912788166199, 0.0187845005297514, 0.634774564676918, 
0.75444408760507, 0.83358903485805, 0.965761296217652, 1.29387996791096, 
-0.136551020847416, -0.440138688134939, -1.2272839144599, -0.237653069522654, 
-0.92685815670431, 0.411235355800902, -0.198864576818102, -0.557449639360495, 
-0.977157122426193, 0.0607352487288863, -0.597229479590505, -1.25894869747501, 
-1.41004720156942, 1.10130407942216, -0.677242026280637, -0.762173507550694, 
-0.291748602024485, -0.575188529371732, -0.443941868186958, -0.312570443962516, 
-0.603004426394854, -1.09393472002918, 0.714706191321827, -0.10881226323069, 
-1.44379793533529, 0.806123263592213, -1.73983507830327, -0.401320309499965, 
-0.287582488751326, -0.938407400051363, 0.28766713879535, -1.50540112432801, 
1.51929701331257, 0.367409355950846, 1.69986240880968, 0.644196977276918, 
-1.6878010286367, 0.647645994369042, 0.448794226781398, 1.02630215826019, 
1.0749782228105, 0.458309613000345, 0.631586758086112, -0.580466398492627, 
1.58419214220344, -1.76399294007705, -1.88062197024023, -1.29171904402856, 
0.909670446302926, -1.10775568207131, -0.384123874715623, 0.082734833711361, 
-0.483882474015024, -2.08474112541176, 1.16758696392551, -0.07682577072231, 
0.53042139527307, 0.00490879581096575, -0.530231429870355, 0.102202409887892, 
0.815459635567165, -1.50652186300026, -1.15718017182801, 1.30155493728303, 
-0.905640775072037, 0.0155634809297453, -1.4280051965488, 0.687979701560331, 
0.594196501448839, -0.270315244004084, 1.55407607614196, -0.510742458166137, 
-0.291842735699669, 1.10143816503647, 0.458665214698786, 0.851908039344607, 
1.46619991886849, 0.243863294906007, 0.64942193391407, -0.641118156058402, 
0.464459770665112, -0.309200447143362, 1.27602070890132, -0.0476112199953048, 
0.0301011659686939, -0.0035802700058728, -2.03457264308331, -0.992819962588807, 
1.9638580560558, 0.706528773018853, 0.0396749797167387, -0.158518589629338, 
-0.162257551192431, 0.438821911652957, 0.78069929204977, -0.981205622027531, 
-0.14158573183272, -1.2562489529667, -0.429033675985112, 0.508776266227222, 
-1.44688974290189, 1.01951282869911, 1.17854697590699, -0.0102587652079676, 
0.26862487113597, 1.34202887245852, -0.583606442977517, -0.937000890812421, 
-1.1521570560791, -0.975183430828597, -0.711413598264887, 0.37191247618326, 
-0.943514417304648, -0.274367308067409, 0.154468494125185, -1.28865320152105, 
1.41910231880672, 1.3078277617162, -1.80497577381603, -0.484069526514124, 
-0.373211780065891, 0.169085902211342, -1.17298360282238, -0.376196736960944, 
0.783127639643547, 0.9911713044315, 0.229284474811702, 0.336850948675991, 
-0.64723199918968, -0.443206682066575, 1.16175176076309, 0.0650510103873427, 
-0.682467580302316, 1.63434165849319, -0.584849244406111, 1.12993600059353, 
0.581201043276578, 0.379293067076565, -0.310740875935373, 0.886390001309243, 
-1.64186474990831, -0.98856374150857, -0.244003429490954, 0.15605692572279, 
0.102051958284238, -0.28745817163283, -0.293194886378227, 0.469477058057219, 
-0.662588512287211, -2.19037266579489, 0.00385467168029454, 0.862942829061133, 
0.980223359085623, -0.291724771104236, -0.0610134758331342, 1.51186840118448, 
-0.643609064922924, -0.130950479435566, 0.485211453408803, -0.64518243454701, 
-0.463956219640376, 0.567249696972656, -0.723407833759501, 0.455884553708724, 
-1.19793462638391, -2.01818055007467, 0.61972340496952, -0.132857466100936, 
-0.434960739564345, -0.521754178741525, 0.992288425621206, -1.08542240785116, 
0.9598324137389, -1.03901054087568, -0.137841967398184, -0.214708418065892, 
0.57371858498562, -1.77644931491255, -0.256392295161705, -0.138550636059799, 
0.557684974345772, 1.12234628895228, -0.982147207446364, 0.326905555961652, 
0.447826608391904, 1.1837241202911), .Dim = as.integer(c(500, 
3)), .Dimnames = list(NULL, c("group", "1", "2")))
"diabetes" <- 
structure(.Data = list("class" = structure(.Data = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 2, 2, 1,
	2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3,
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
	3, 3, 3, 3, 3)
, class = "factor"
, levels = c("chemical", "normal", "overt")
)
, "glucose" = c(80., 97., 105., 90., 90., 86., 100., 85., 97., 97., 91., 87., 78., 90., 86.,
	80., 90., 99., 85., 90., 90., 88., 95., 90., 92., 74., 98., 100., 86.,
	98., 70., 99., 75., 90., 85., 99., 100., 78., 106., 98., 102., 90.,
	94., 80., 93., 86., 85., 96., 88., 87., 94., 93., 86., 86., 96., 86.,
	89., 83., 98., 100., 110., 88., 100., 80., 89., 91., 96., 95., 82.,
	84., 90., 100., 86., 93., 107., 112., 94., 93., 93., 90., 99., 93.,
	85., 89., 96., 111., 107., 114., 101., 108., 112., 105., 103., 99.,
	102., 110., 102., 96., 95., 112., 110., 92., 104., 75., 92., 92., 92.,
	93., 112., 88., 114., 103., 300., 303., 125., 280., 216., 190., 151.,
	303., 173., 203., 195., 140., 151., 275., 260., 149., 233., 146., 124.,
	213., 330., 123., 130., 120., 138., 188., 339., 265., 353., 180., 213.,
	328., 346.)
, "insulin" = c(356., 289., 319., 356., 323., 381., 350., 301., 379., 296., 353., 306., 290.,
	371., 312., 393., 364., 359., 296., 345., 378., 304., 347., 327., 386.,
	365., 365., 352., 325., 321., 360., 336., 352., 353., 373., 376., 367.,
	335., 396., 277., 378., 360., 291., 269., 318., 328., 334., 356., 291.,
	360., 313., 306., 319., 349., 332., 323., 323., 351., 478., 398., 426.,
	439., 429., 333., 472., 436., 418., 391., 390., 416., 413., 385., 393.,
	376., 403., 414., 426., 364., 391., 356., 398., 393., 425., 318., 465.,
	558., 503., 540., 469., 486., 568., 527., 537., 466., 599., 477., 472.,
	456., 517., 503., 522., 476., 472., 45., 442., 541., 580., 472., 562.,
	423., 643., 533., 1468., 1487., 714., 1470., 1113., 972., 854., 1364.,
	832., 967., 920., 613., 857., 1373., 1133., 849., 1183., 847., 538.,
	1001., 1520., 557., 670., 636., 741., 958., 1354., 1263., 1428., 923.,
	1025., 1246., 1568.)
, "sspg" = c(124., 117., 143., 199., 240., 157., 221., 186., 142., 131., 221., 178., 136.,
	200., 208., 202., 152., 185., 116., 123., 136., 134., 184., 192., 279.,
	228., 145., 172., 179., 222., 134., 143., 169., 263., 174., 134., 182.,
	241., 128., 222., 165., 282., 94., 121., 73., 106., 118., 112., 157.,
	292., 200., 220., 144., 109., 151., 158., 73., 81., 151., 122., 117.,
	208., 201., 131., 162., 148., 130., 137., 375., 146., 344., 192., 115.,
	195., 267., 281., 213., 156., 221., 199., 76., 490., 143., 73., 237.,
	748., 320., 188., 607., 297., 232., 480., 622., 287., 266., 124., 297.,
	326., 564., 408., 325., 433., 180., 392., 109., 313., 132., 285., 139.,
	212., 155., 120., 28., 23., 232., 54., 81., 87., 76., 42., 102., 138.,
	160., 131., 145., 45., 118., 159., 73., 103., 460., 42., 13., 130.,
	44., 314., 219., 100., 10., 83., 41., 77., 29., 124., 15.)
)
, class = "data.frame"
, names = c("class", "glucose", "insulin", "sspg")
, row.names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14",
	"15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26",
	"27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38",
	"39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50",
	"51", "52", "53", "54", "55", "56", "57", "58", "59", "60", "61", "62",
	"63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74",
	"75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86",
	"87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98",
	"99", "100", "101", "102", "103", "104", "105", "106", "107", "108",
	"109", "110", "111", "112", "113", "114", "115", "116", "117", "118",
	"119", "120", "121", "122", "123", "124", "125", "126", "127", "128",
	"129", "130", "131", "132", "133", "134", "135", "136", "137", "138",
	"139", "140", "141", "142", "143", "144", "145")
)

"wreath" <-
structure(c(-11.2386964979405, 1.83115966271336, -13.2755616065126, 
-6.08965018034, 8.61834966314144, -6.97126691060875, 11.1604720347314, 
2.41218223350522, -15.4800712913441, 2.34046274622134, -14.3445764668188, 
-7.36455307851838, 10.1004267320028, 10.4165990684670, 1.93657719402692, 
-15.0145168966254, 9.34894164793177, 9.617430913839, 9.85494084129942, 
-14.5921718923085, -15.2151910510441, -3.22071957140354, -2.51011766255831, 
-5.41931522240214, 10.6759689838825, 1.53783588213465, 8.74446230846057, 
-6.44958679680057, 2.09171095744346, -5.92407872565621, -9.06940056790329, 
9.32648536540075, -10.0121381615291, -14.5621730191282, -14.9035112386180, 
9.77432469128358, 8.98104329914836, -13.5967559183586, 10.6681754379799, 
2.34022914024075, -13.6919378768203, 10.0056493442224, 3.22804417757, 
10.7216681915226, 2.46685181575292, -6.47000879100595, -13.0484290588526, 
12.4872704323511, 8.83741149014849, -5.24311752881296, -11.6072197927575, 
10.8584639029657, 10.2083826324439, 10.4820575564566, 16.3558610873592, 
-15.1448825261699, -13.6588021101666, -14.8611905093341, 2.91194475203669, 
-9.97503242669595, -3.32452605437299, -15.9719083049796, -9.55247759056975, 
10.3112169749218, -15.1183211336617, -3.16595113062191, -5.58371301719195, 
-2.24458242102976, -6.72347495582868, 8.09041582568345, 11.2152182268657, 
3.80298444279985, -4.27501562821953, 10.4599340740174, 17.0267393668817, 
-6.35109182636597, -3.41932811291792, 10.0994768117288, 9.10865210076086, 
-7.16698037602642, 10.6474059829683, -7.76770486833855, -2.95084132694201, 
-5.74776036731826, 2.04311945879799, -7.19839696342708, 9.02487256886354, 
18.3050780348422, 9.23014153186874, 16.9074478849641, -14.5564930175761, 
-4.96272863323056, -16.3672509231369, 10.6538850529445, 10.6334878964175, 
9.54078126733262, -15.7686677881934, 11.6954982846992, 8.99993356917881, 
11.6816101290350, 10.9528785066284, 9.94542153763287, 2.45921768225351, 
10.5311961203549, -6.01147813059003, 10.0013807529396, 15.2981291353173, 
17.4941746923735, -11.4895221264813, 15.3856161180067, 9.95710647196474, 
1.80315865694241, 9.39828503170493, -4.89391088011224, -15.8045976018996, 
-10.7121647499478, 3.57292215632835, 10.1227279063471, 3.54056197674299, 
9.85313955858608, 15.7357873345509, 2.77317114433339, 3.3563828315034, 
-4.00341432032181, -5.36478570849816, 2.10764946513349, 9.72420249568023, 
-14.3738847597705, 15.9300084217070, 9.88880012660264, 16.4844810970787, 
-15.5877853544767, 2.76090436760367, -10.8856958466543, 4.00545755493249, 
10.7214169438977, -10.7815340453610, 16.1624577703711, -15.7182473716481, 
10.2270240825303, -14.2138645611473, -3.09931978786650, -5.26184296067376, 
-5.5697131323875, 2.19372658234258, 10.1122029601578, 9.0705449408925, 
-2.23768490865656, 2.76749930060041, -11.0350806558014, -12.5106002693895, 
-5.79467259867304, -15.6231928227805, 2.12822527560618, -4.2478594969284, 
10.0007820145153, 15.0110283138323, -14.3670421525776, 10.1443381879506, 
-6.81954550124883, -7.87641071076614, 15.3575275526899, -14.5117033451974, 
-3.07070436250906, -15.5009348665837, 9.18078096967811, -10.260555813068, 
-3.49468499700614, -11.798128505515, 10.4680297881393, 1.47157226693817, 
-10.722502164296, 1.92687483336205, 16.5518055005626, 9.9357571982357, 
10.7341066073410, -14.5784273429487, 1.51792633777752, -14.9892644322061, 
-14.5962207653211, 2.95724085660481, -11.1952390825595, 2.42545304305339, 
2.3587813003271, -3.75551695514248, -0.50401812206966, 15.1503528190449, 
-2.71133795844576, -9.33246311721238, 10.2235761212933, 2.35435788263379, 
-3.98927542100194, 15.4973821784626, -14.6924781613748, 1.83177322937805, 
-7.46312819521342, -9.97640673468906, -6.77295934939427, -6.37514589480303, 
-11.2646331989725, -11.8432632937618, 2.3976612608984, 10.6390198237299, 
3.99089055934635, -6.56086902670986, -5.78537425129797, -16.120233549222, 
-8.83776368660502, 2.79936954516133, -7.05311052206471, 9.77003169391222, 
4.360987946973, 11.3858956758756, 12.4427209991856, -7.94859440588261, 
10.0722474793232, -5.54378022566211, 10.8525391333768, 17.0211320214723, 
9.20775972130904, -6.56499951986763, -13.3244218088278, -7.15211172857176, 
-15.8599477843073, -3.00876443882427, 9.56551540142705, -3.28745374128052, 
1.90851637719091, -2.1273820400969, 3.46862662195685, -3.24153360930814, 
-6.23768899418298, 10.4928916430503, -10.3808222860904, -5.66265757211923, 
9.68208772825144, 11.7412661628536, 16.8179315693341, -2.3789642880603, 
-10.4888447005172, 15.0008147659872, -8.49878169945668, -4.42114889701737, 
9.37968533000916, 3.20801318842253, 15.9546900377614, 9.9502899480819, 
11.6436777959449, 3.44200583570006, 3.77949750037942, 9.85429472345102, 
2.65803561173403, -5.51189808277265, -8.61669285703423, -13.6658413160229, 
-11.7496838495353, -5.13950700083953, -0.442462352922628, -11.6924362868038, 
-13.9131567706766, -6.70559582379996, 1.82899657614080, 9.96927876496193, 
-2.86888382342606, 15.1059838512221, 10.0525567717615, -3.42871032773052, 
-12.2232772354555, 8.69878875482753, -14.6884885716608, -8.2130283026391, 
-16.9820826110057, 17.7432257993048, -12.0031065367718, -10.8773834401704, 
1.89508724816144, 4.75463423413484, -3.91417913651252, 10.0180166420854, 
-6.5823783932861, 8.9112502027601, 15.5321467022147, -15.3166664693606, 
4.89928936977969, 10.5963044232854, 11.4466118934799, -15.1144818794580, 
-7.01683006539529, 11.5147818771352, -15.0675738541593, 1.82788096984990, 
11.0038159237422, 9.9584866074009, -15.7540100615528, 14.3957533794042, 
-14.9504164688526, -13.0678042723483, -15.2687903500396, -8.34564613856797, 
9.9044392395632, -4.42543605141263, 14.6663587005292, -11.3132342526793, 
2.16019063013915, 2.78124883555558, -15.4972801604997, 10.9584102330689, 
-12.9970864840165, 10.3886188088644, -14.891889655147, -1.07005202909553, 
-11.7042666651863, -16.1138280241897, 10.4927403004998, -8.87755355377131, 
-5.40826346788424, -3.87031511858271, 2.84969554774167, -14.2879800851385, 
9.59899337264933, 9.89814870243902, -6.52574477808707, 10.8955068781074, 
15.1344945935759, 11.6047910110159, -15.1621308681804, -8.26775168250843, 
-3.77250453256282, 9.85773906676595, 1.72648781064607, -4.31194859490819, 
-3.26050963127663, -7.35651967817361, -13.4900142315358, -2.20700058389674, 
-8.36536299596102, -6.1522580863184, 8.95819600563956, 16.7521864086132, 
10.3099059463516, 9.0031706899185, 11.0693901681791, -10.9357801328478, 
-13.8669516371000, -0.488191546280416, 11.9223770659272, -2.37579791784729, 
-4.59266599925245, -4.18069506566463, -14.5424675478197, -8.92489974521917, 
-14.412569841921, 8.55673493194602, -17.0904238593959, 10.6264015315213, 
9.88870464326725, 10.7828473114595, -3.21563378924548, -14.6862204923150, 
-3.31531046799333, -11.7538214847493, -5.1295641148069, -6.04301654730991, 
17.9195451519294, -7.97657673111423, -4.86170320664296, -14.4525680946154, 
9.05759606650975, 2.12957314192643, 10.5138972269540, 3.23241069471021, 
-13.0263271769719, 11.9283091027524, -8.15507176734, 10.137359518901, 
-4.76317188736922, 17.3688222964071, 2.56720768324838, 10.3077940139887, 
9.78093858070392, -1.56720527202330, 4.04883213986221, -4.23872643663579, 
-10.4123388146418, -7.36914707402065, 3.21970930577393, 9.67994797878788, 
10.8628931862140, 17.4654183377448, 1.29250722410589, -6.11761857649584, 
2.96728995527122, -7.33727919467158, 10.0292775711781, -11.2937399191323, 
0.500580996149582, -7.57257273758771, -15.1664460046299, 5.59313148027172, 
-3.95961082606566, -10.5036991297925, 10.9207076534575, -3.5882222547685, 
2.30501147935406, 17.2431779563282, 12.2900658772424, -2.27746591030934, 
-14.5709993255700, 9.5960125272682, 9.96662625401746, -3.83204525942611, 
2.58784119128867, 10.3975693821373, 16.2745007834994, -11.1559951226741, 
8.93750542963907, -14.2572414202525, 17.1433074650989, -5.9650079500302, 
-11.7105415235415, 16.8992791992384, -3.48779011783783, -7.3853776842185, 
8.22571718712929, 2.85950025034937, 9.999462418492, -2.62943430107748, 
17.4278024039205, -11.1993552812158, -6.17597476791329, -4.56100543347593, 
16.7051843518175, -10.4818698868865, 15.5181055346196, 10.1612131116338, 
11.0306016941394, 17.2722317887119, 2.35617384545703, -2.86233373968986, 
2.79187339281202, 9.84216654549, -16.7516793498269, 10.1769972301452, 
-7.21826426078727, 9.667094723072, 10.8656961193246, -7.8435486971012, 
-15.0322224184389, -9.62073439217585, 9.36259393046984, -4.4307510929508, 
1.3190179379882, 1.37940699299537, 9.74811437332064, 9.3059426021506, 
-13.8861201186586, 18.0234586471180, -7.61604940436484, -2.52622041895487, 
10.2748332769311, 1.01017721775308, -8.43276252167552, -3.95712238180623, 
11.3095905041828, 17.5508720152637, -7.52865491639488, -14.0957194339371, 
10.0922276153715, 10.9180359755636, 10.1054530055654, 8.5343867631544, 
-9.25539083748225, -6.29682730083617, -13.3183430967517, -6.91988806623811, 
-16.2493616127969, -9.7238966526933, 9.09261370999002, -15.9994770347614, 
2.73038361792296, 9.6155067836161, -3.54308986435847, -6.96488031382164, 
-3.21004981532537, 11.2525309901907, -11.8435898181696, 1.13186659429983, 
11.4560085268317, 5.0061185958409, 16.3056605259224, 2.71904039913201, 
14.9233494560940, 16.1183141562360, -3.11766293513710, -15.4241912267593, 
17.319433197037, -14.3377670228763, 11.2573987796670, 2.35654007116905, 
7.91483519509626, 9.03613805642758, 1.79172766934661, -15.8395275541943, 
-3.51025480835507, 17.8987268570069, 10.5543542790803, 9.2436316465111, 
-10.7942534720045, 9.97442046719675, -6.69479457352644, 10.4259987297729, 
-15.1701069092044, 11.6309741263563, 9.39008153300713, 2.44012784899782, 
2.39313398810698, 9.75999479227037, -4.15690916729588, -15.4703152409982, 
-8.37433433855386, 10.9475822266629, -17.2744260274279, 1.81337286467103, 
-14.7855115558539, 9.07669268579858, 8.85577721948855, -8.19335821383445, 
1.88158525786932, -3.58392705350903, 9.47586680445077, -2.90284700820525, 
-5.67424467611289, 9.60484468556706, 0.638541880549261, 10.6264102940772, 
-14.1156688501051, 1.96375780835315, -2.38490044586928, -2.55584455975025, 
-13.8940735701212, 11.2639150747794, 1.34818529209421, -9.95900527732529, 
8.98638050446266, 2.75077558197308, -7.96368413485251, -16.3702425361679, 
11.0615080060331, 9.5265652263236, 12.0077599210567, 1.95327849946340, 
10.7039438584050, -3.57121790083117, 16.7166069690944, 10.5195245468196, 
1.75003642847380, -14.9106485127422, -14.0893670235633, 10.4381832010566, 
-16.2453607275592, -4.66296473076568, -4.43926621870093, 2.14154482139578, 
-4.96355279142349, 10.9395773797830, 16.6549869877003, -3.70391535590729, 
9.56685564495026, -4.20391213982419, 10.0241760894034, -3.58397156937559, 
9.72022425741758, -7.75900865391029, -4.88337797620157, 11.1275256678536, 
1.57091742460019, 1.53455859496511, 14.9835106342893, 9.10102870565992, 
1.31423993502662, 15.6475703302671, 9.74715857145425, 1.81347939565084, 
3.90815091302834, -11.7580163050377, -6.50566852215227, -8.35221698568819, 
-6.51994708736881, -10.0065732481801, -13.5047389286737, -7.19997593744441, 
-11.6096144763732, -3.56603256402707, -6.82134587002845, 9.90546552178776, 
-5.20134738534694, -7.08701911511189, 17.5906794987032, 4.26353529079196, 
-7.53684965286602, 2.39015001653393, 16.4671227242036, -3.92662263096532, 
10.4723821777678, 9.03671485560433, -7.06175560730069, -7.8843513365158, 
10.5351161967170, 2.60580490265640, -7.76811134403203, -3.45454625772472, 
-14.2434049078736, 10.0319959552353, 9.87952366558173, 9.00932816779263, 
9.23063052114108, 3.58844752372151, 8.45462488552199, 15.5416024306633, 
-4.86111427603746, 17.6092333159264, 10.8899097147208, -6.78821874951393, 
-8.12672256001221, -14.4940425719552, -13.5014387173946, 10.2942874424512, 
-5.42624669725885, -4.54003158976207, -6.95723090076794, -13.6387845977892, 
-7.366025198762, 3.54620376788196, 10.6353880657664, -5.30249417128239, 
10.0905253264699, 9.52849472016373, 2.49674824684802, -6.84580386063646, 
10.1639435792708, 2.16324076574626, -13.3974090681049, 16.9893800541817, 
-12.3966825746052, -6.43666626349215, -7.49960436630336, 8.92757441564112, 
-4.25161231086135, -3.70293467121071, -9.14577741948215, 1.12793476709212, 
-13.2730956035163, 9.52347190791168, -13.3485210089539, -11.651249128058, 
-10.9507927860490, 9.29139474642793, 10.2597437980508, 7.58804775333839, 
8.85969915370348, -11.8070195918728, 8.8580872386524, -10.9312383588121, 
3.37004657015535, 16.1610578867449, -6.20799120962959, -6.3379393868028, 
8.8862306150944, -4.02594932324065, -15.0512301844940, 11.2368314866606, 
-4.67218671316451, -8.39987182214553, -7.33919202456771, -2.99919131500728, 
17.3085561706313, -2.670729036938, -3.65457141219382, 17.7435082927323, 
7.88241322342226, -2.97545135334451, -6.16677021737042, 9.07881142450722, 
11.4894488687227, 10.4516013383809, 9.01141311302206, -7.79140704250658, 
16.7923428707389, 9.58097379169746, -15.4051615074937, -14.5545317372536, 
0.81905051067818, 14.9624734390047, -12.2819455463693, 2.87428241800367, 
-15.7164618493054, 1.88362177956713, -12.1292542688906, -5.88647533108693, 
-6.03630236228699, -1.96442686986989, 10.1676083857308, -14.833189431956, 
11.7044609614261, -12.1004602573813, 9.77190619685575, 16.6559127531448, 
-13.6509181925919, 9.21198198868027, 1.23922555417386, -12.6288430941127, 
10.4163740311041, 10.3996116785619, -11.6773547734866, 8.62349254580602, 
-3.92576648789726, -14.2797551136098, -7.50003806834817, 2.79120839724323, 
10.8463357195677, 10.8353116518891, 0.477994211816012, -15.6612178534413, 
-15.1064188221933, -7.6388726508189, -3.88430298986173, -11.1233341628078, 
10.7799400466072, 16.8646040992879, 16.8455875296525, -7.35176131919292, 
-7.12800293339285, -4.38234672090008, 10.1576550737418, -14.4659237892992, 
-10.5081040834174, 9.80938561839667, 16.6151051828565, -4.66761208553589, 
14.7997898770119, -4.5826156823598, 3.09059308284727, -15.8174242855252, 
10.2727364851880, 9.73117671522287, 10.5245127663555, -3.59562312918278, 
-13.0247240397470, 8.97144380779272, 14.9240397119324, 10.4590440731700, 
10.0503811600504, -11.6223420804872, 9.13100225679144, 11.0724202146024, 
16.7118695519528, -13.1161682773607, -6.57296902778447, -14.8789001878213, 
-4.15596555393055, -6.75885782574575, -5.79585274496402, -3.42869870280976, 
2.53168681475086, 2.13197592174722, -5.86806603408053, 10.8595034206253, 
-3.98648723227565, 9.3706162334928, -13.5996645221853, 2.93237311431414, 
0.102156265390823, -15.0604819798228, 16.6560863641913, -5.80903008307203, 
-11.0259182520448, -8.41132572886349, -16.6219606045212, -15.3446586386741, 
2.77299105572585, 1.92781869890196, 2.02066794538683, 2.17759020615433, 
-3.47595242139475, 10.0606644872596, -15.0850202127938, 10.7573284001417, 
9.83161591588974, -14.1189390230158, 10.5843410690955, 17.8230508603530, 
3.03699425651433, -3.86761392250578, -14.9774026413260, -15.4830245982538, 
8.1791527088244, 2.67764726620177, 10.5377530373156, 10.1372696909593, 
10.8281501254544, -5.54517752839591, 9.37959414139635, 2.13843000089106, 
11.0946369187401, 15.8414581207384, 8.4291173571179, -3.04404563386369, 
3.24462657544557, -2.25220830937196, -15.2835242934661, 10.4710120260902, 
-15.0393178817624, -11.3846098579582, 10.8899278632066, -15.5416178353660, 
-6.61671128800518, 2.22282889074544, 17.4631007865983, 2.34357592337344, 
-2.08879993800750, 9.61141094069922, 15.6218584932282, 14.2538903367713, 
-2.48188418011941, 9.14031632614819, 15.3273268324015, -13.6778798102939, 
9.85362945358164, 9.2036846852088, -2.82466341323237, -4.83736031883026, 
10.8404189947462, 3.37747608693880, -10.3865198955795, 10.8647451892570, 
-14.1332146080208, 2.83550058011166, 16.7973286289587, -10.9611234241571, 
0.316874610548145, -5.66048359452901, 9.63912618815169, 9.80032413148861, 
7.92406590770553, -10.1652644338858, -9.91757827082545, 8.98017913401257, 
-7.00774607942056, -2.52455877915015, 9.6140753718901, 9.07270308526547, 
16.7115546978795, -2.50021321868451, 3.28250968014336, 9.15439165298725, 
-3.97151016623668, 16.6210130668188, -13.4455286853795, 9.99793567089976, 
10.1105894732479, 11.1345065836635, -12.9409809160113, 9.7093242366545, 
-7.62473610580534, -15.3577438870573, 9.9545170150867, 9.8876681890375, 
-10.3048029218030, -5.41835626439192, 10.4107681831807, -11.8180015553745, 
9.72895128746265, 10.3309468884837, 9.71146702087681, 9.61007029980527, 
10.3903275709965, -4.72414078279521, -14.3104109309545, 11.0805389317456, 
10.3182731197429, -1.47396374614396, 10.7058440967968, 3.1384193469694, 
-3.70967248930602, 2.37701143094543, -5.68646473058713, -15.0204512594926, 
1.88956779935909, -15.3765690826572, 10.1194596174309, -5.11734582179228, 
-15.1231955700812, 2.6510644465053, 9.10942616776687, 1.94855619306723, 
8.69759748003102, 10.4436250891719, 15.6071901584842, 8.66855520440082, 
9.6619987703089, -3.84739546346505, -13.6982441566520, 10.4048210217573, 
9.22088133660368, -4.20828890220959, -15.1207154353537, -14.785935379741, 
6.56317582483018, 10.2147164283287, -6.97431088385615, -4.35309257538295, 
-11.6738061164525, 10.3170723523606, -16.5248698128151, -3.74357971274233, 
-3.33263193695214, 8.82805359435328, -14.1822348338325, 10.9016615148907, 
17.1636867838628, 9.29253911051622, 16.2829324525069, 2.12659331276857, 
16.2300012846835, 10.0380612081840, 9.32548569627807, 11.3685820836308, 
-16.0943171052415, -12.0029455081456, -7.78450546777121, 2.68270518761699, 
-6.50851946814591, -14.4403644617114, 10.2481193119614, -10.6861224788176, 
3.15202842527633, -4.53791754211828, 10.6422449656944, 10.8986563947834, 
10.6820480811321, 2.54578165115354, 10.5677513108456, 1.06350558476753, 
9.71833992875844, -6.04139041705728, 10.3065569568089, 1.98727812074518, 
-4.92624196878213, 9.87900158054387, -6.16021058877239, 10.8440859196655, 
-2.90421983109618, -16.4330383932040, 10.6902692081713, -7.7790449815649, 
10.5751635621164, 9.2679830345754, -14.5690444565269, 10.2054131332415, 
-4.85903529833702, 9.90630334808874, -16.0891913463696, -14.3782742758195, 
11.0579757720963, 2.78804988130309, 10.4935787119213, 8.9542906424322, 
-6.52628899189127, 9.17052189587863, -3.43151001795364, 9.27587587605473, 
10.6488092109470, 1.01659966158639, 5.71866298040969, -7.12900847178424, 
10.3616089965773, 7.59029793642208, 10.2704559148254, -4.53910342845411, 
9.61392658108726, -7.34564633007988, -15.9825973496907, 10.5286323851111, 
11.5985004009261, 9.77658696159269, 8.99930657351907, 10.2388724479805, 
-13.4777591794855, -5.10393701180477, 11.2072532123306, 16.2208158123691, 
-7.22571008070076, 11.0778138672083, -7.53944852654606, -7.11386098477548, 
2.27162321173453, -16.4384013993388, 10.7382776642573, 10.3108410195358, 
9.64110315319856, 0.669564498161358, -11.1224087641629, -1.21827947457496, 
9.85395688895126, 4.56985172643432, -9.52755354510191, 3.63818234541943, 
-11.6831692177871, 6.27617290028265, 9.894731906944, -7.0293013286745, 
-8.44993325324441, -11.5118947248609, -12.5584915932373, 9.52252955688902, 
6.61840290901133, 6.25522054806308, 12.947093907926, 12.3935262437724, 
-6.98922259013254, -6.90457536019144, -17.6221328894949, -18.1296671798680, 
-7.1194086964891, 3.66885832827975, 12.2624493224423, 4.89122140164327, 
-8.60449774629653, -10.6668593501492, 8.80159717445751, 2.23413396726873, 
-3.73086339223512, 1.65466045371694, -6.13411629932314, 7.64821095690774, 
4.58977659618752, 3.94936632301877, 5.29937818301589, -13.9540357611350, 
12.2285824943602, -5.81884691313834, -4.62776853325884, -10.9722978296856, 
13.2482273803242, 9.95733188369791, -9.37188717019116, -0.0374895442560701, 
-17.126569731316, 10.1296322773939, 15.4053770970427, 0.0678259170532346, 
-11.8784862963666, 4.95014888888384, -4.60421492792592, 0.319933248290643, 
6.79138814372422, 6.44379610476804, -7.40491228731225, 10.2642413321490, 
-0.511151223032745, -16.0043058061785, 6.79645018168384, 0.228579837193977, 
-12.2049964364673, -7.08458214552706, 17.4311436291543, 16.1366754266881, 
17.1839321245091, 9.0810006156266, -5.31376701805263, -13.9096790483278, 
-9.75781669889026, 15.7420219894337, 13.3688640160761, -1.08676506468952, 
8.73315206553107, -16.0220124283456, -13.7602275666167, -11.5794958876244, 
-8.6872747995579, -4.91675093366072, -9.8709402702754, 17.1692292224619, 
-6.97407703987337, 11.3501454093044, -8.37267974530858, -4.30767900385631, 
-0.888027519847287, -4.86489618203382, 0.343217273769061, 7.48524849660966, 
-7.35623227574026, -7.54798083385105, -3.92876039021238, -14.0818197718053, 
5.57006600936992, -8.84099364063592, -4.17077991201085, -5.3212074954239, 
-14.2733298076254, -13.3576525009957, 12.5641125762396, 11.0383303111395, 
13.7074663065070, -7.57369387457501, -7.29967074151825, 1.01493124156095, 
0.323411825260864, -0.0992039339008624, -2.53933246353360, 13.8604033506006, 
11.2345312635542, 5.05244108908745, 14.7318378030121, 6.67508984111829, 
0.172638650023211, -9.06636742088698, 11.9166164959574, -8.85711998342802, 
-11.330598386401, 1.19523407082394, -10.3346315333226, -10.6677906547045, 
-14.0623750668234, -7.925262984237, 11.2882311312835, -13.252702344747, 
-7.5451298152351, 1.14113313384533, 12.7113385220659, 0.950278563748496, 
-7.2606139979041, 11.6390333190491, 1.94496181851706, -10.1128966961782, 
12.5395092688817, 0.530995784359435, -0.250107990470434, 7.83022256116742, 
13.5139812699990, 6.90564444812948, -15.9756078852999, -7.02002537430106, 
15.2198682449111, -10.3314090994179, -5.35837774600508, -11.2105203436832, 
-16.7801037802198, 9.5740659271606, -2.06907595168055, -1.16203145172392, 
-7.51938411431563, -6.77051497692064, -11.1854593705623, -16.2492694887237, 
-11.8895661683720, -0.302574486700522, 7.22617851579874, -12.1381768392708, 
-8.05685560236947, -8.40765435279984, 0.488643651682308, -7.0658723623929, 
-16.7075763214283, -7.77265922756926, 4.02619991408887, 0.781017068363895, 
16.4811753269198, -1.82503298399345, -2.87088572644481, 12.1479978274480, 
-0.262674087723081, 11.6131630817751, 0.603572153973005, -2.58278658105523, 
-13.6872487061209, 7.79393097931622, 11.1255796393861, -7.6084323082668, 
8.02418342881206, 11.3575519971288, -0.838560605420274, -10.7783349494576, 
-11.113677189184, -17.5349291754641, 15.1368635873108, 0.493011567302912, 
14.6671321445162, 6.96939129602477, -12.9175267213546, 11.2787381598009, 
-14.5799387507504, 0.295514625017510, -6.35751963705468, -10.7758242610982, 
-8.50967433471131, -0.134915403038966, 9.66828443790014, -8.39756550190324, 
-0.478638718703674, -0.476002488468145, 8.71995293237802, -13.404100328763, 
9.8571029695055, -9.15263822649187, 10.4742712480317, 8.9139895604706, 
0.542372136686739, -10.7594639050601, 9.40968653914639, 11.3740432855404, 
-11.8772378960120, 13.6758980292966, 4.75481607593675, 8.64333614854457, 
4.06778289869539, 14.3504609940131, 4.16519944809867, -0.964242624808928, 
-5.83938936163377, 8.55920037128453, -6.95431302839074, 10.6012134642147, 
-8.39339299048576, 16.3594323032942, -4.05130046106252, -18.3098511673368, 
10.8470221293998, -14.6952671138587, 11.3634173784853, 16.5587484578978, 
9.35357361963851, 13.1881416021801, 0.0344767515041787, -7.32367604190981, 
-6.28422706888944, 14.8895939177338, 1.12243764641219, -15.8153727143943, 
0.707659832235996, 0.71416969582718, -9.46684671012853, 16.0013445281825, 
-3.98019654399036, -10.6420581571407, -0.538748431596675, -12.8780578355531, 
-15.4445661637694, 11.1852532988945, -10.2440345610826, 5.57550767141347, 
9.27148384649385, 14.8625391991403, 8.21719315575151, -6.9736431969179, 
0.727568692152246, -7.73588543856142, -10.1567358380362, -1.40485146262265, 
7.30905430448811, -7.74286630019117, -9.72291957349947, 5.10259750418355, 
16.4571511105101, 1.45195585958731, 3.8839675301166, -17.6266489295792, 
-0.839385736503454, 6.2207874533019, -6.00770404147201, -9.73794855143163, 
-7.19832392670321, 0.301804862040879, -1.47020960923405, 0.353419268499872, 
11.3801966943273, -9.2518183576395, -15.4882877534571, 14.613825588567, 
8.8858572776798, 4.72364433503849, 1.48039307596794, -5.7877969717161, 
8.2461692692884, -6.1179660298182, 14.8841728582108, 8.47803425374655, 
-9.36503904691948, 3.79121496564572, 7.74681947253907, -11.0347431198037, 
14.3493043782300, -13.4973742690661, -6.62827631774599, 2.49781144230973, 
5.29591051898038, 6.30755147909517, 7.03990150014745, -7.99225207282498, 
4.84969656216721, -15.8003187191447, 0.365027911544345, 0.845729124400444, 
-10.5310767663955, 9.03100851516943, -7.27055752285814, -5.74517627095755, 
-1.5447949156769, 12.8779924046420, -7.50603595660607, 15.9272276797929, 
-1.19503981249166, -6.40480238557624, 12.1126821701190, 6.99590048731604, 
-8.27679271244322, 16.5833597436398, -11.4147432866094, -5.85146134680287, 
-12.3743721075226, -13.1978167710023, 9.16177551239802, -12.9513730301134, 
0.571386644970654, -14.5217142135577, -6.52288816814937, 8.32490088251531, 
15.3605080461602, 12.5802928179715, 11.9674031136031, -17.3645218653737, 
-15.8792018452742, -8.96149036350123, 6.56924068332295, -17.2693743748136, 
7.29595177419094, -7.78778435046516, -7.76190643313986, 0.689175089488489, 
-3.75079580500959, -4.8326721290741, -12.8069975013006, 0.249711424593709, 
7.18298001775735, 11.9080423003454, 4.30092807031768, -15.2530360215426, 
15.7304926261928, -16.1950376277854, -8.5005004951268, -8.8964700487664, 
5.03615028810655, 4.88937739981627, 8.56842757298954, 13.7945232337316, 
-2.99833195578187, -12.7561519956817, 15.4881546926816, -6.32394721875537, 
15.2522149025948, -1.07554395867047, -6.90109881295476, 10.7738786064756, 
-1.06485110654353, -7.63488844943663, 9.86365194620722, 7.60789783580282, 
6.09496155250842, -10.1625366097617, 5.0215697482585, 9.87350819530297, 
5.94740929264849, -12.8757402430549, 6.9073509477314, 4.25402929768251, 
-12.1410994899870, -0.320793924006678, -11.7101519744935, -11.6607026385618, 
12.0311292184200, 17.2533841503797, 9.33177790546494, -15.9252150992959, 
0.443720531903639, 8.14874375250565, -10.4881684949123, 12.0852639012961, 
-4.55956032056736, -1.28954213120317, -12.7259668395799, 9.6202674504421, 
9.43631306471659, -7.710458010068, -12.1769869465287, 0.171428076198571, 
10.7909744596209, 6.81647294122911, 7.92922459331989, -9.15200057451904, 
-17.2185644200707, 0.502892755023551, -12.9704746408454, 14.9183038727402, 
-10.1033816928871, -0.852279624363491, 3.99457543283617, -17.5451186443964, 
7.28929576542792, 11.2390583769267, -11.7032100284982, 15.3916549612568, 
-11.5839087759133, -12.6045700979750, 0.325607833128144, 0.187005711617548, 
5.5171451316534, -1.02150447151481, 0.82043915505651, 8.39951795605516, 
-0.869093863717084, -0.403488780932644, 14.4252944105789, -8.45465316691037, 
5.75245212115319, 9.91475943249583, -2.37768999674011, -17.4317063999103, 
-0.749218819619638, 0.376409846202254, -8.80392254286776, -14.6689714969162, 
-0.708638699416182, -0.557201092158183, 1.41432228630067, -4.9696057936468, 
12.5476759626673, -1.15219027445493, 10.7816202245585, -16.5016406983670, 
-10.3566488311027, -3.62256869534719, -8.25440427024394, -11.8133909347809, 
-7.82262721643602, -11.5145037178435, -14.2003123547243, 8.89752168784053, 
-6.43417429390791, 5.27679356425655, -13.2191479993779, -14.9318529867288, 
-10.8856018287140, 10.4056012796619, 11.5928112020379, -4.95864216907872, 
5.58254316266925, -0.409983399902434, 8.62675530067775, 15.4597937159888, 
11.6366305602196, 11.8119368533821, 7.25190633114951, -16.3597730912949, 
4.66891159734818, -0.480392180967629, -8.6529987549181, -6.65787090834176, 
-11.9304185643888, 4.65339834530443, 11.154088083672, 12.9393093505660, 
0.291559529943779, 10.264753852518, -7.22241673379078, -8.82178876285456, 
7.71946969328553, 0.648859989202802, -6.24204831854886, 6.63224460243811, 
11.0650578276382, -13.7853195971773, -14.7976492594334, -8.98589223617527, 
-17.1156376657967, 4.773716196325, 0.286395843566046, -10.1735341902673, 
4.09471723317832, -11.1608645200006, -0.83381874016489, 11.2524105311559, 
-0.673599914898767, 0.427258611908803, -17.6796120485094, 6.95806202397006, 
-1.17335417297171, -7.60280938906679, 13.8944995446409, 9.86070616320226, 
5.87780882517215, 4.28812726002423, 11.1001298653007, 9.04217235114844, 
16.6148315563604, -0.662561813766104, 4.21770985898378, 11.776950624882, 
0.665002047696374, 10.3639777795745, -7.47238434782462, 13.0737584308980, 
-8.37631736762643, 13.7597583451637, -12.7958933506755, -11.2566055747386, 
11.5532229444652, 12.6406852928228, -15.5402461003191, -6.51710695427806, 
-8.4734291246741, 14.5316359216352, -7.46040493512675, -10.1799982442507, 
-5.49108553289755, 10.8969625777593, -4.96542606793784, 7.58176735346337, 
11.2328161422690, 14.9557757876323, -4.46026411927597, -17.364830618143, 
-8.0374867584922, -13.6988819788909, -10.1893250561728, 14.3811125472977, 
6.2363620513889, -10.4066578811096, -18.1919144777136, -15.5214687884638, 
8.96164770754233, -13.1999147839269, -11.8783340306753, -0.153726394533027, 
-7.42563741057596, -10.0413594499031, -8.74592066973293, -8.63888343097298, 
-3.04757878982657, -10.5374773950548, 4.52940024914457, 11.5686859167173, 
4.52858016976236, -16.6419046005611, 0.570663228935474, -11.8487287799365, 
-10.3432205111372, 7.14289135598123, 7.08045878629339, 12.9809019889313, 
7.57656466635153, -15.9523051055966, -16.1398839753712, -10.8528295159058, 
15.0307194396418, -4.34969643792615, -1.04337877643836, 16.0596702289004, 
5.67465337828717, 16.7928106807815, -13.3301368134429, 16.0251208960328, 
13.1816656405121, -9.0347447337774, 16.6793748015069, 4.34851154244023, 
12.2774474544538, -11.6825590073493, 0.989205375215975, 4.33455383978332, 
-10.2277344373732, 0.767803009566575, -11.8295291763760, -10.4213027890236, 
8.74733950643028, -4.03736803410318, 8.7560767028001, -9.31272850025152, 
9.29411065921633, -0.133116299453187, 6.7720352658845, 9.284712738966, 
0.0383967677973552, 16.9256727156522, -8.5084118425064, 13.2215682373925, 
15.8421141316462, 9.32885191643153, 0.281914251177875, -9.79026537694126, 
-9.23934524188988, -10.7498303647591, 0.624658756196318, 15.6900407366708, 
11.3714550742289, -12.5792156679881, 8.35362436480984, 8.48903451586498, 
13.5454951876721, 12.7901674488367, 9.40398446258466, -14.0222591567398, 
6.22941696488162, -4.93562113602407, -13.7303900402433, 5.6294730846409, 
12.4233759866629, 10.5022300777322, -5.3422958100677, 1.49462062111552, 
-15.5303790075732, -0.474853310399765, -15.0006801320148, 8.19492788966676, 
8.87986037966079, -7.69740269045705, 5.92421876616792, -13.3660831553906, 
-8.05777593109512, 15.7091844311861, 7.90992954207523, 6.51887828868774, 
7.59298601839383, 9.69062749015023, -12.7465046890236, -9.1392251206403, 
12.9402214538936, 12.6770693693889, 10.7763181942584, 8.46802009719553, 
3.61104841629957, 11.8034350282936, -6.55468864551859, -0.890897203411001, 
-0.488726070616372, -8.91298157137203, -7.17927459085432, -6.72983019974658, 
16.3551805836507, -15.9399447799024, 0.5141372596416, -10.3529110508638, 
-5.71913571108229, -4.92197058597302, -0.610559572740231, 6.2649384973737, 
-0.799571559595107, -6.31311359335425, 13.1479763691012, 5.60221640840798, 
5.18281325633641, -0.496505588904928, 9.89027068494463, 0.229691259895004, 
8.9608520673356, -0.675448921516744, -8.75450790895652, 8.09045413821645, 
-13.0457932663558, -14.5095908375917, -7.7883828237196, 3.43604734281269, 
15.3004137915534, -8.19613935920084, 8.65304469745153, -16.4016806012537, 
-0.461026541147531, 16.4097921390742, 17.2217513938662, -1.15638789083156, 
-5.93220828124551, -8.40342931523463, -8.5002462669983, 5.2412895988449, 
-2.25989731625749, 11.6426009615087, 10.2151840223568, 8.95487226735071, 
0.254207581351533, 4.92759122528375, -6.32854373893162, 5.79581220807568, 
-11.4856638209533, 0.396137608134576, -1.2653208267338, 9.746164202381, 
7.11825424887194, -11.5684822912435, 0.108008300224074, -7.03222231137959, 
7.73878094926483, 15.7452081105860, -3.55126770183500, 6.51184085532122, 
14.4480985898393, 1.28276229406795, 3.61767945044127, 0.338197035392082, 
7.861609671528, 12.7289562215536, -11.2216355441924, 5.40534843319926, 
-4.91937409814718, 13.7926306609608, -0.427539065851059, -5.29954902774542, 
-14.1317949818434, 6.57636611503613, -8.57451605272355, -8.97778844049929, 
-13.1418739581660, -12.2058214658171, -10.7374102503629, -7.88590948147416, 
-6.49201569321482, -8.97389070900898, -15.516871491674, 0.307256931428659, 
12.6593326183567, 0.207416018199114, 0.102233520767882, 9.78436207669583, 
7.89577224084169, -15.6765668222844, -4.18125722866093, 8.35214583137549, 
1.48493171311897, -13.1020256400900, -1.16686114590398, -14.3631296373892, 
0.97436304461095, -13.4574101048635, -11.7129428617852, 9.49348064563155, 
-12.4212732621118, -11.3745664368126, 13.6059783814440, -14.5480136584309, 
-6.37866967099029, 12.5021854527990, 1.8729018252887, 13.3680756044623, 
-5.87811424049503, 0.0652192342684699, 4.7606875105475, -12.5893828559028, 
-0.353154808725318, -5.37707225980112, 9.61807152767236, 8.85263901178718, 
-14.4027608399923, -9.5717439790055, 14.6143243790518, -16.0805018601439, 
10.3323693022545, 10.2680411632092, -7.91366037136095, -4.55476707598394, 
17.7618594848188, -11.7338302467223, 5.66478232754316, 10.0568720646262, 
11.7978560400290, 7.19358700835289, -1.27016863878705, 10.2097406115686, 
-0.303535393180968, -9.6455839039492, -7.85274533309063, -7.27459182366418, 
-9.9028140329518, -9.59359426847224, 9.88845600804822, 12.2515630054559, 
-17.3570245168556, 12.1889122472111, 7.35318237064061, -3.48701898352817, 
-14.6322947557962, 7.34533098528654, -4.12103051255532, 0.629396992935123, 
10.1800897915410, -15.6145460884448, -8.25316459969651, 7.47036965148855, 
-6.6854265200523, 10.6336147402601, -3.93096998951069, -12.1263284727554, 
-14.5582417175287, 15.5164631213845, -13.1332026287899, 11.8915706259184, 
12.7076475776909, -0.0515532077600599, -10.5460672160236, -14.7985729130601, 
-11.2217781486136, 16.5921628722896, 7.32390049820996, 5.40124875692726, 
-7.3489049541825, 0.331284711501529, -11.4180404478619, 6.587711906788, 
9.34635902614025, 11.1803712367057, -1.42897097329507, 12.3269842482458, 
16.6843117361577, 4.46140682922035, 0.894196479876898, 1.2292275527635, 
16.8527941139532, -4.02123419720994, 0.850539917346756, -6.04200102591898, 
4.52590167505428, 14.3777597046871, 14.7772997545932, -14.7036634655102, 
-4.48255716673165, 9.35766394950108, -0.0685420364837732, -4.82982468378193, 
5.7110259182823, -10.7453979266732, -0.348948595000408, -0.938665440759572, 
-10.6497764779988, 10.5768104993793, -5.81490795915895, 3.40895103686266, 
6.3714109907443, 0.135983345353485, 0.0910037614209204, 4.64981371057071, 
8.26899494727962, 17.2469565724694, 4.73052696094166, 4.73818441302043, 
0.439953200002556, 16.5843323366297, 10.3477010291782, -5.30653513387452, 
-15.1630104615242, -0.087162446298646, 5.72216246930388, 6.92705124092693, 
-3.9220685823848, -4.13085259127815, -6.02442702386006, 3.4893116430093, 
6.8643704793719, 7.25989644920063, 12.9501097984654, -3.1298683973973, 
1.40679088663524, -7.5106111643187, -3.81451457967619, -0.394177718216364, 
-5.83862432024646, 12.8370141681628, -4.02108322926193, 11.9562047196719, 
-4.18420239777577, 16.6882960539463, -7.50112018520996, -13.2409291424255, 
-11.8545679769953, 16.4845479686077, 12.8124525757916, -10.7630165043531, 
-15.3069955632818, 11.101606749794, 8.47066278752448, -4.589169301351, 
-11.4363810459258, 6.43727181054747, -13.6610435572057, 15.9178248187276, 
7.81765278473105, -10.6489258893580, 5.2813694308993, 10.1928168259915, 
-4.07202657801250, 11.1452046030866, 0.118739095381662, -6.33011507279038, 
4.52575941265368, -17.3139220266939, -7.50350143630174, -14.2184712290569, 
-12.6414464145081, -15.6583971350054, -6.79140599704174, 7.5458577195777, 
4.11407869859976, -3.28405982955024, 7.55511097017177, -15.5236213972406, 
0.75404866673091, 3.6971661544258, -8.28851134281774, 15.6714129541334, 
15.5883393823809, 4.96541604473415, 3.54436829441541, 13.6614360525626, 
-0.559223096880561, 5.0423105859329, -0.736871935756264, 10.5305044888496, 
-0.624104065849514, -3.76001397276165, -13.1044678254193, -3.00565248342931, 
-7.51027315717966, 0.281650861973847, -8.33155350227017, 11.0805308283235, 
-7.96352673389653, 5.0161494753301, -13.2449199465834, -0.974318849543622, 
10.9028961294168, 15.5337953006257, -4.95398315777852, -12.6522103803535, 
-13.7619214218485, 8.77364588488233, -13.0385500190376, -10.7148436538757, 
-4.52047810815942, 9.63138424539835, -5.60117575501695, -11.2835540840346, 
-13.3784509432975, 13.5996174590418, 7.90241822619217, 3.59529767787274, 
15.6609512324120, -8.0376550667201, 13.9573547807454, 8.0698269667341, 
4.74772545651483, -4.99241442147477, -5.91622519982884, -14.9735260927174, 
-14.5332213151386, 11.0057822511197, 9.79276552003525, 7.5597457654126, 
5.28498197720025, -10.3003310173208, -12.0206177559917, -13.186199668982, 
-7.28878249246488, 5.17301415465024, 14.8652952309111, 4.95177951253865, 
12.375006058162, -11.1777424180643, 4.28658382505351, -8.93793580540988, 
-14.0253465505712, 5.22050700708083, -6.19041976337743, -14.2342142587558, 
4.48391666082333, -8.49462761336629, -8.09326468783904, 11.8680579531135, 
-13.7733675300449, 12.6444702750609, -4.41750223691927, 12.9883139511511, 
-1.33167956883833, -7.72594822431541, -4.09641994520465, -0.287963260434826, 
8.00613599526203, -12.9265950595204, 7.42524681719568, 7.22231001422963, 
-10.6542776383927, -7.93257313086292, 14.2701678005275, -5.24179932597196, 
13.7822061864296), .Dim = as.integer(c(1000, 2)))
