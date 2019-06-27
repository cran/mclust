\name{classError}
\alias{classError}

\title{Classification error}

\description{
Computes the errore rate of a given classification relative to the known classes, and the location of misclassified data points.}

\usage{
classError(classification, class)
}
\arguments{
  \item{classification}{
    A numeric, character vector or factor specifying the predicted class 
    labels. Must have the same length as \code{class}.
  }
  \item{class}{
    A numeric, character vector or factor of known true class labels. 
    Must have the same length as \code{classification}.
  }
  
}
\value{
  A list with the following two components:
  \item{misclassified}{
    The indexes of the misclassified data points in a minimum error
    mapping between the predicted classification and the known true classes.
  }
  \item{errorRate}{
    The error rate corresponding to a minimum error mapping 
    between the predicted classification and the known true classes.
  }
}

\details{
  If more than one mapping between predicted classification and the known 
  truth corresponds to the minimum number of classification errors,
  only one possible set of misclassified observations is returned.
}

\seealso{
  \code{\link{map}}
  \code{\link{mapClass}},
  \code{\link{table}}
}
\examples{
(a <- rep(1:3, 3))
(b <- rep(c("A", "B", "C"), 3))
classError(a, b)

(a <- sample(1:3, 9, replace = TRUE))
(b <- sample(c("A", "B", "C"), 9, replace = TRUE))
classError(a, b)

class <- factor(c(5,5,5,2,5,3,1,2,1,1), levels = 1:5)
probs <- matrix(c(0.15, 0.01, 0.08, 0.23, 0.01, 0.23, 0.59, 0.02, 0.38, 0.45, 
                  0.36, 0.05, 0.30, 0.46, 0.15, 0.13, 0.06, 0.19, 0.27, 0.17, 
                  0.40, 0.34, 0.18, 0.04, 0.47, 0.34, 0.32, 0.01, 0.03, 0.11, 
                  0.04, 0.04, 0.09, 0.05, 0.28, 0.27, 0.02, 0.03, 0.12, 0.25, 
                  0.05, 0.56, 0.35, 0.22, 0.09, 0.03, 0.01, 0.75, 0.20, 0.02),
                nrow = 10, ncol = 5)
cbind(class, probs, map = map(probs))
classError(map(probs), class)
}
\keyword{cluster}
