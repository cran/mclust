\name{cvMclustDA}
\alias{cvMclustDA}

\title{MclustDA cross-validation}

\description{
V-fold cross-validation for classification models based on Gaussian 
finite mixture modelling.
}

\usage{
cvMclustDA(object, nfold = 10, 
           prop = object$prop,
           verbose = interactive(), 
           \dots) 
}

\arguments{
  \item{object}{
    An object of class \code{'MclustDA'} resulting from a call to 
    \code{\link{MclustDA}}.
  }
  \item{nfold}{
    An integer specifying the number of folds (by defaul 10-fold CV is 
    used).
  }

  \item{prop}{
    A vector of class prior probabilities, which if not provided default
    to the class proportions in the training data. 
  }

 \item{verbose}{
    A logical controlling if a text progress bar is displayed during 
    the cross-validation procedure. By default is \code{TRUE} if the 
    session is interactive, and \code{FALSE} otherwise.
  }
  
  \item{\dots }{Further arguments passed to or from other methods.}
}

\details{
The function implements V-fold cross-validation for classification 
models fitted by \code{\link{MclustDA}}. 
Classification error and Brier score are the metrics returned, but other
metrics can be computed using the output returned by this function
(see Examples section below).
}

\value{
The function returns a list with the following components:

  \item{classification}{a factor of cross-validated class labels.} 

  \item{z}{a matrix containing the cross-validated probabilites for class assignment.} 

  \item{error}{the cross-validation classification error.}
  
  \item{se.error}{the standard error of the cross-validated classification error.}

  \item{brier}{the cross-validation Brier score.} 
    
  \item{se.brier}{the standard error of the cross-validated Brier score.}
}

\author{Luca Scrucca}

\seealso{
  \code{\link{MclustDA}}, 
  \code{\link{predict.MclustDA}}, 
  \code{\link{classError}},
  \code{\link{BrierScore}}
}

\examples{
\donttest{
# Iris data
Class <- iris$Species
X <- iris[,1:4]

## EDDA model with common covariance (essentially equivalent to linear discriminant analysis)
irisEDDA <- MclustDA(X, Class, modelType = "EDDA", modelNames = "EEE")
cv <- cvMclustDA(irisEDDA)                         # 10-fold CV (default)
str(cv)
cv <- cvMclustDA(irisEDDA, nfold = length(Class))  # LOO-CV
str(cv)

## MclustDA model selected by BIC
irisMclustDA <- MclustDA(X, Class)
cv <- cvMclustDA(irisMclustDA)                     # 10-fold CV (default)
str(cv)

# Banknote data
data("banknote")
Class <- banknote$Status
X <- banknote[,2:7]

## EDDA model selected by BIC
banknoteEDDA <- MclustDA(X, Class, modelType = "EDDA")
cv <- cvMclustDA(banknoteEDDA)                     # 10-fold CV (default)
str(cv)

(ConfusionMatrix <- table(Pred = cv$classification, Class))
TP <- ConfusionMatrix[1,1]
FP <- ConfusionMatrix[1,2]
FN <- ConfusionMatrix[2,1]
TN <- ConfusionMatrix[2,2]
(Sensitivity <- TP/(TP+FN))
(Specificity <- TN/(FP+TN))
}
}

\keyword{multivariate}
