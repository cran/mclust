\name{bicE}
\alias{bicE}
\alias{bicV}
\alias{bicEII}
\alias{bicVII}
\alias{bicEEI}
\alias{bicVEI}
\alias{bicEVI}
\alias{bicVVI}
\alias{bicEEE}
\alias{bicEEV}
\alias{bicVEV}
\alias{bicVVV}
\title{
  BIC for a Parameterized MVN Mixture Model
}
\description{
  Compute the BIC (Bayesian Information Criterion) for a parameterized
  mixture model given the loglikelihood, the dimension of the data, 
  and number of mixture components in the model.
}
\usage{
bicE(loglik, n, G, equalPro, noise = F, \dots)
bicV(loglik, n, G, equalPro, noise = F, \dots)
bicEII(loglik, n, d, G, equalPro, noise = F, \dots)
bicVII(loglik, n, d, G, equalPro, noise = F, \dots)
bicEEI(loglik, n, d, G, equalPro, noise = F, \dots)
bicVEI(loglik, n, d, G, equalPro, noise = F, \dots)
bicEVI(loglik, n, d, G, equalPro, noise = F, \dots)
bicVVI(loglik, n, d, G, equalPro, noise = F, \dots)
bicEEE(loglik, n, d, G, equalPro, noise = F, \dots)
bicEEV(loglik, n, d, G, equalPro, noise = F, \dots)
bicVEV(loglik, n, d, G, equalPro, noise = F, \dots)
bicVVV(loglik, n, d, G, equalPro, noise = F, \dots)
}
\arguments{
  \item{loglik}{
    The loglikelihood for a data set with respect to the MVN mixture model.
  }
  \item{n}{
    The number of observations in the data used to compute
    \code{loglik}.
  }
  \item{d}{
    The dimension of the data used to compute \code{loglik}.
  }
  \item{G}{
    The number of components in the MVN mixture model used to compute
    \code{loglik}.
  }
  \item{equalPro}{
    A logical variable indicating whether or not the components in the
    model are assumed to be present in equal proportion. The default is
    \code{.Mclust\$equalPro}.
  }
  \item{noise}{
    A logical variable indicating whether or not the model includes and
    optional Poisson noise component. The default is to assume that the
    model does not include a noise component. 
  }
  \item{\dots}{
    Catch unused arguments from a \code{do.call} call.
  }
}
\value{
  The BIC or Bayesian Information Criterion for the MVN mixture
  model and data set corresponding to the input arguments.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002a).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
  See \url{http://www.stat.washington.edu/mclust}. 
  
  
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and 
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{bic}},
  \code{\link{EMclust}},
  \code{\link{estepE}},
  \code{\link{mclustOptions}},
  \code{\link{do.call}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

n <- nrow(irisMatrix)
d <- ncol(irisMatrix)
G <- 3

emEst <- meVVI(data=irisMatrix, unmap(irisClass))
names(emEst)

bicVVI(loglik=emEst$loglik,n=n,d=d,G=G)
do.call("bicVVI", emEst)  ## alternative call
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
