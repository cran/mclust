\name{partuniq}
\alias{partuniq}
\title{
  Classifies Data According to Unique Observations
}
\description{
  Gives a one-to-one mapping from unique observations to rows of a data matrix.
}
\usage{
partuniq(x)
}
\arguments{
  \item{x}{Matrix of observations.}
}
\value{
  A vector of length \code{nrow(x)} with integer entries. An observation
  \code{k} is assigned an integer \code{i} whenever observation \code{i}
  is the first row of \code{x} that is identical to observation \code{k}
  (note that \code{i <= k}).
}
\seealso{
   \code{\link{partconv}}
}
\examples{
set.seed(0)

mat <- data.frame(lets = sample(LETTERS[1:2],9,TRUE), nums = sample(1:2,9,TRUE))
mat

ans <- partuniq(mat)
ans

partconv(ans,consec=TRUE)
}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.


