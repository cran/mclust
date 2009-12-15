\name{imputeData}
\alias{imputeData}
\title{
  Missing Data Imputation via the mix package
}
\description{
  Imputes missing data using the mix package.
}
\usage{
imputeData(x, categorical = NULL, seed = NULL)
}
\arguments{
  \item{x}{
    A numeric vector, matrix, or data frame of observations containing
    missing values. Categorical variables are allowed. If a matrix
    or data frame, rows correspond to observations and columns
    correspond to variables. 
  }
  \item{categorical}{
    A logical vectors whose \emph{i}th entry is \code{TRUE} if the
    \emph{i}th variable or column of \code{x} is to be interpreted as
    categorical and \code{FALSE} otherwise. The default is to assume that a
    variable is to be interpreted as categorical only if it is a factor.
  }
  \item{seed}{
    A seed for the function \code{rngseed} that is used to initialize
    the random number generator in \code{mix}. By default, a seed is
    chosen uniformly in the interval \code{(.Machine$integer.max/1024,
    .Machine$integer.max)}.
  }
}
\value{
  A dataset of the same dimensions as \code{x} with missing values
  filled in.
}
\section{References}{
  J. L. Schafer, Analysis of Imcomplete Multivariate Data,
  Chapman and Hall, 1997.
}
\seealso{
  \code{\link{imputePairs}}
}
\examples{
library(mix)

# impute the continuos variables in the stlouis data
stlimp <- imputeData( stlouis[,-(1:3)])

# plot imputed values
imputePairs( stlouis[,-(1:3)], stlimp)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
