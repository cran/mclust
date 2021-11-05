\name{imputeData}
\alias{imputeData}
\alias{matchCluster}

\title{Missing data imputation via the \pkg{mix} package}

\description{
  Imputes missing data using the \pkg{mix} package.
}

\usage{
imputeData(data, categorical = NULL, seed = NULL, verbose = interactive())
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations containing
    missing values. Categorical variables are allowed. If a matrix
    or data frame, rows correspond to observations and columns
    correspond to variables. 
  }
  \item{categorical}{
    A logical vectors whose \emph{i}th entry is \code{TRUE} if the
    \emph{i}th variable or column of \code{data} is to be interpreted as
    categorical and \code{FALSE} otherwise. The default is to assume that a
    variable is to be interpreted as categorical only if it is a factor.
  }
  \item{seed}{
    A seed for the function \code{rngseed} that is used to initialize
    the random number generator in \pkg{mix}. By default, a seed is
    chosen uniformly in the interval \code{(.Machine$integer.max/1024,
    .Machine$integer.max)}.
  }
  \item{verbose}{
    A logical, if \code{TRUE} reports info about iterations of the algorithm.
  }
}
\value{
  A dataset of the same dimensions as \code{data} with missing values
  filled in.
}
\references{
Schafer J. L. (1997). Analysis of Imcomplete Multivariate Data, Chapman and Hall.
}
\seealso{
  \code{\link{imputePairs}}
}
\examples{
\donttest{
# Note that package 'mix' must be installed
data(stlouis, package = "mix")
 
# impute the continuos variables in the stlouis data
stlimp <- imputeData(stlouis[,-(1:3)])

# plot imputed values
imputePairs(stlouis[,-(1:3)], stlimp)
}
}
\keyword{cluster}
