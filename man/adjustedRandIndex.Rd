\name{adjustedRandIndex}
\alias{adjustedRandIndex}
\title{
  Adjusted Rand Index
}
\description{
   Computes the adjusted Rand index comparing two classifications. 
}
\usage{
adjustedRandIndex(x, y)
}
\arguments{
  \item{x}{
    A numeric or character vector of class labels.
  }
  \item{y}{
    A numeric or character vector of class labels.
    The length of \code{y} should be the same as that of \code{x}.
  }
}
\value{
  The adjusted Rand index comparing the two partitions (a scalar).  
  This index has zero expected value in the case of random partition, and it is bounded above by 1 in the case of perfect agreement between two partitions.  
}

\references{
 L. Hubert and P. Arabie (1985) Comparing Partitions, \emph{Journal of the Classification}, 2, pp. 193-218.
}
\seealso{
  \code{\link{classError}},
  \code{\link{mapClass}},
  \code{\link{table}}
}
\examples{
a <- rep(1:3, 3)
a
b <- rep(c("A", "B", "C"), 3)
b
adjustedRandIndex(a, b)

a <- sample(1:3, 9, replace = TRUE)
a
b <- sample(c("A", "B", "C"), 9, replace = TRUE)
b
adjustedRandIndex(a, b)

a <- rep(1:3, 4)
a
b <- rep(c("A", "B", "C", "D"), 3)
b
adjustedRandIndex(a, b)

irisHCvvv <- hc(modelName = "VVV", data = iris[,-5])
cl3 <- hclass(irisHCvvv, 3)
adjustedRandIndex(cl3,iris[,5])

irisBIC <- mclustBIC(iris[,-5])
adjustedRandIndex(summary(irisBIC,iris[,-5])$classification,iris[,5])
adjustedRandIndex(summary(irisBIC,iris[,-5],G=3)$classification,iris[,5])
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
