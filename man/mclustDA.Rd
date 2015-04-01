\name{MclustDA}
\alias{MclustDA}
\alias{print.MclustDA}

\title{MclustDA discriminant analysis}

\description{
Discriminant analysis based on Gaussian finite mixture modeling.
}

\usage{
MclustDA(data, class, G = NULL, modelNames = NULL, 
         modelType = c("MclustDA", "EDDA"), 
         prior = NULL, 
         control = emControl(), 
         initialization = NULL, 
         warn = mclust.options("warn"), 
         \dots)
}

\arguments{
  \item{data}{
    A data frame or matrix giving the training data.
  }
  \item{class}{
    A vector giving the class labels for the observations in 
    the training data.}
    
 \item{G}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the BIC is to be calculated within each class. 
    The default is \code{G = 1:5}.\cr
    A different set of mixture components for each class can be specified
    by providing this argument with a list of integers for each class. 
    See the examples below.
  }
  
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    by EM within each class (see the description in 
    \code{\link{mclustModelNames}}).
    A different set of mixture models for each class can be specified
    by providing this argument with a list of character strings.
    See the examples below.
  }
   
  \item{modelType}{
    A character string specifying whether the models given in
    \code{modelNames} should fit a different number of mixture 
    components and covariance structures for each class 
    (\code{"MclustDA"}, the default) or should be constrained 
    to have a single component for each class with the same covariance 
    structure among classes (\code{"EDDA"}).
    See Details section and the examples below.
  }

  \item{prior}{
    The default assumes no prior, but this argument allows specification of a 
    conjugate prior on the means and variances through the function 
    \code{\link{priorControl}}.
  }
  
  \item{control}{
    A list of control parameters for EM. The defaults are set by the call
    \code{emControl()}. 
  }
  
  \item{initialization}{
    A list containing zero or more of the following components:
  \describe{
  \item{\code{hcPairs}}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{hc}. The default is to compute a hierarchical
    clustering tree by applying function \code{hc} with
    \code{modelName = "E"} to univariate data and
    \code{modelName = "VVV"} to multivariate data or a
    subset as indicated by the \code{subset} argument. 
    The hierarchical clustering results are used as starting values 
    for EM.}
  \item{\code{subset}}{
    A logical or numeric vector specifying a subset of the data
    to be used in the initial hierarchical clustering phase.}
  }
  }
  
  \item{warn}{
    A logical value indicating whether or not certain warnings
    (usually related to singularity) should be issued when
    estimation fails. 
    The default is controlled by \code{\link{mclust.options}}.
  }
  
  \item{\dots }{Further arguments passed to or from other methods.}
}

\value{
 An object of class \code{'MclustDA'} providing the optimal (according 
 to BIC) mixture model.
 
 The details of the output components are as follows:

  \item{call}{The matched call.} 

  \item{data}{The input data matrix.} 
  
  \item{class}{The input class labels.} 
    
  \item{models}{A list of \code{\link{Mclust}} objects containing information
  on fitted model for each class.} 
    
  \item{n}{The total number of observations in the data.}

  \item{d}{The dimension of the data.}

  \item{BIC}{All BIC values.}
  
  \item{bic}{Optimal BIC value.}
    
  \item{loglik}{Log-likelihood for the selected model.}
  
  \item{df}{Number of estimated parameters.}
}

\details{
The \code{"EDDA"} method for discriminant analysis is described in Bensmail and Celeux (1996), while \code{"MclustDA"} in Fraley and Raftery (2002).
}

\section{References}{
  Bensmail, H., and Celeux, G. (1996)
  Regularized Gaussian Discriminant Analysis Through Eigenvalue
  Decomposition.
  \emph{Journal of the American Statistical Association}, 91, 1743-1748.
  
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association}, 97, 611-631. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\author{Luca Scrucca}

\seealso{
  \code{\link{summary.MclustDA}}, 
  \code{\link{plot.MclustDA}}, 
  \code{\link{predict.MclustDA}}, 
  \code{\link{classError}}
}
\examples{
odd <- seq(from = 1, to = nrow(iris), by = 2)
even <- odd + 1
X.train <- iris[odd,-5]
Class.train <- iris[odd,5]
X.test <- iris[even,-5]
Class.test <- iris[even,5]

# common EEE covariance structure (which is essentially equivalent to linear discriminant analysis)
irisMclustDA <- MclustDA(X.train, Class.train, modelType = "EDDA", modelNames = "EEE")
summary(irisMclustDA, parameters = TRUE)
summary(irisMclustDA, newdata = X.test, newclass = Class.test)

# common covariance structure selected by BIC
irisMclustDA <- MclustDA(X.train, Class.train, modelType = "EDDA")
summary(irisMclustDA, parameters = TRUE)
summary(irisMclustDA, newdata = X.test, newclass = Class.test)

# general covariance structure selected by BIC
irisMclustDA <- MclustDA(X.train, Class.train)
summary(irisMclustDA, parameters = TRUE)
summary(irisMclustDA, newdata = X.test, newclass = Class.test)

plot(irisMclustDA)
plot(irisMclustDA, dimens = 3:4)
plot(irisMclustDA, dimens = 4)

plot(irisMclustDA, what = "classification")
plot(irisMclustDA, what = "classification", newdata = X.test)
plot(irisMclustDA, what = "classification", dimens = 3:4)
plot(irisMclustDA, what = "classification", newdata = X.test, dimens = 3:4)
plot(irisMclustDA, what = "classification", dimens = 4)
plot(irisMclustDA, what = "classification", dimens = 4, newdata = X.test)

plot(irisMclustDA, what = "train&test", newdata = X.test)
plot(irisMclustDA, what = "train&test", newdata = X.test, dimens = 3:4)
plot(irisMclustDA, what = "train&test", newdata = X.test, dimens = 4)

plot(irisMclustDA, what = "error")
plot(irisMclustDA, what = "error", dimens = 3:4)
plot(irisMclustDA, what = "error", dimens = 4)
plot(irisMclustDA, what = "error", newdata = X.test, newclass = Class.test)
plot(irisMclustDA, what = "error", newdata = X.test, newclass = Class.test, dimens = 3:4)
plot(irisMclustDA, what = "error", newdata = X.test, newclass = Class.test, dimens = 4)

\dontrun{
# simulated 1D data
n <- 250 
set.seed(1)
triModal <- c(rnorm(n,-5), rnorm(n,0), rnorm(n,5))
triClass <- c(rep(1,n), rep(2,n), rep(3,n))
odd <- seq(from = 1, to = length(triModal), by = 2)
even <- odd + 1
triMclustDA <- MclustDA(triModal[odd], triClass[odd])
summary(triMclustDA, parameters = TRUE)
summary(triMclustDA, newdata = triModal[even], newclass = triClass[even])
plot(triMclustDA, what = "scatterplot")
plot(triMclustDA, what = "classification")
plot(triMclustDA, what = "classification", newdata = triModal[even])
plot(triMclustDA, what = "train&test", newdata = triModal[even])
plot(triMclustDA, what = "error")
plot(triMclustDA, what = "error", newdata = triModal[even], newclass = triClass[even])

# simulated 2D cross data
data(cross)
odd <- seq(from = 1, to = nrow(cross), by = 2)
even <- odd + 1
crossMclustDA <- MclustDA(cross[odd,-1], cross[odd,1])
summary(crossMclustDA, parameters = TRUE)
summary(crossMclustDA, newdata = cross[even,-1], newclass = cross[even,1])
plot(crossMclustDA, what = "scatterplot")
plot(crossMclustDA, what = "classification")
plot(crossMclustDA, what = "classification", newdata = cross[even,-1])
plot(crossMclustDA, what = "train&test", newdata = cross[even,-1])
plot(crossMclustDA, what = "error")
plot(crossMclustDA, what = "error", newdata =cross[even,-1], newclass = cross[even,1])
}
}
\keyword{multivariate}
