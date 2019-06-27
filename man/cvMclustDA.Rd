\name{cvMclustDA}
\alias{cvMclustDA}

\title{MclustDA cross-validation}

\description{
K-fold cross-validation for discriminant analysis based on Gaussian finite mixture modeling.
}

\usage{
cvMclustDA(object, nfold = 10, 
           metric = c("error", "brier"),
           prop = object$prop,
           verbose = interactive(), \dots) 
}

\arguments{
  \item{object}{
    An object of class \code{'MclustDA'} resulting from a call to \code{\link{MclustDA}}.
  }
  \item{nfold}{
    An integer specifying the number of folds.
  }
    
  \item{metric}{
    A character string specifying the statistic to be used in the 
    cross-validation resampling process. Possible values are \code{"error"} 
    for the classification error, and \code{"brier"} for the Brier score.
  }

  \item{prop}{
    A vector of class prior probabilities, which if not provided default
    to the class proportions in the training data. 
  }

 \item{verbose}{
    A logical controlling if a text progress bar is displayed during the cross-validation procedure. By default is \code{TRUE} if the session is interactive, and \code{FALSE} otherwise.
  }
  
  \item{\dots }{Further arguments passed to or from other methods.}
}

\value{
The function returns a list with the following components:

  \item{classification}{a factor of cross-validated class labels.} 

  \item{z}{a matrix containing the cross-validated probabilites for class assignement.} 

  \item{error}{the cross-validation classification error if \code{metric = "error"}, \code{NA} otherwise.}
  
  \item{brier}{the cross-validation Brier score if \code{metric = "brier"}, \code{NA} otherwise.} 
    
  \item{se}{the standard error of the cross-validated statistic.}

}

%\details{}

\author{Luca Scrucca}

\seealso{
  \code{\link{summary.MclustDA}}, 
  \code{\link{plot.MclustDA}}, 
  \code{\link{predict.MclustDA}}, 
  \code{\link{classError}},
  \code{\link{BrierScore}}
}

\examples{
\dontrun{
X <- iris[,-5]
Class <- iris[,5]

# common EEE covariance structure (which is essentially equivalent to linear discriminant analysis)
irisMclustDA <- MclustDA(X, Class, modelType = "EDDA", modelNames = "EEE")
cv <- cvMclustDA(irisMclustDA)  # default 10-fold CV
cv[c("error", "se")]

cv <- cvMclustDA(irisMclustDA, nfold = length(Class))  # LOO-CV
cv[c("error", "se")] 

cv <- cvMclustDA(irisMclustDA, metric = "brier")  # 10-fold CV with Brier score metric
cv[c("brier", "se")]

# general covariance structure selected by BIC
irisMclustDA <- MclustDA(X, Class)
cv <- cvMclustDA(irisMclustDA)  # default 10-fold CV
cv[c("error", "se")]
cv <- cvMclustDA(irisMclustDA, metric = "brier")  # 10-fold CV with Brier score metric
cv[c("brier", "se")]
}
}

\keyword{multivariate}
