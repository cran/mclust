\name{mclust-package}
\alias{mclust-package}
\alias{mclust}
\docType{package}

\title{
Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation
}
\description{
Finite Gaussian mixture modeling fitted via EM algorithm for model-based clustering, classification, and density Estimation, including Bayesian regularization and dimension reduction.
}

\details{Main functions to use:
\tabular{ll}{
  Clustering: \tab \code{\link{Mclust}} \cr
  Classification: \tab \code{\link{MclustDA}} \cr
  Density estimation: \tab \code{\link{densityMclust}} \cr
  Dimension reduction: \tab \code{\link{MclustDR}} 
  }
}

\author{
Chris Fraley, Adrian Raftery and Luca Scrucca.

Maintainer: Luca Scrucca \email{luca@stat.unipg.it}
}

\section{References}{
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
  
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611:631}. 
}

\examples{
# Clustering
mod1 = Mclust(iris[,1:4])
summary(mod1)
plot(mod1,  what = c("BIC", "classification"))

# Classification
data(banknote)
mod2 = MclustDA(banknote[,2:7], banknote$Status)
summary(mod2)
plot(mod2)

# Density estimation
mod3 = densityMclust(faithful$waiting)
summary(mod3)
plot(mod3, faithful$waiting)
}

\keyword{package}
