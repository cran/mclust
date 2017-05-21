\name{cvMclustDA}
\alias{cvMclustDA}

\title{MclustDA cross-validation}

\description{
K-fold cross-validation for discriminant analysis based on Gaussian finite mixture modeling.
}

\usage{
cvMclustDA(object, nfold = 10, verbose = interactive(), \dots) 
}

\arguments{
  \item{object}{
    An object of class \code{'MclustDA'} resulting from a call to \code{\link{MclustDA}}.
  }
  \item{nfold}{
    An integer specifying the number of folds.
  }
    
 \item{verbose}{
    A logical controlling if a text progress bar is displayed during the cross-validation procedure. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise.
  }
  
  \item{\dots }{Further arguments passed to or from other methods.}
}

\value{
The function returns a list with the following components:

  \item{classification}{a factor of cross-validated class labels.} 

  \item{error}{the cross-validation error.} 
    
  \item{se}{the standard error of cv error.}

}

%\details{}

\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

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
X <- iris[,-5]
Class <- iris[,5]

# common EEE covariance structure (which is essentially equivalent to linear discriminant analysis)
irisMclustDA <- MclustDA(X, Class, modelType = "EDDA", modelNames = "EEE")
cv <- cvMclustDA(irisMclustDA) # default 10-fold CV
cv[c("error", "se")]

cv <- cvMclustDA(irisMclustDA, nfold = length(Class)) # LOO-CV
cv[c("error", "se")] 
# compare with 
# cv1EMtrain(X, Class, "EEE")

# general covariance structure selected by BIC
irisMclustDA <- MclustDA(X, Class)
cv <- cvMclustDA(irisMclustDA) # default 10-fold CV
cv[c("error", "se")]
}

\keyword{multivariate}
