\name{bicEMtrain}
\alias{bicEMtrain}
\title{Select models in discriminant analysis using BIC}
\description{
For the ten available discriminant models the BIC is calulated. The
models for one-dimensional data are "E" and "V"; for higher dimensions
they are "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EEV", "VEV"
and "VVV". This function is much faster than \code{cv1EMtrain}.
}
\usage{
bicEMtrain(data, labels, modelNames)
}
\arguments{
  \item{data}{A data matrix}
  \item{labels}{Labels for each row in the data matrix}
  \item{modelNames}{Vector of model names that should be tested.}
}
\value{
  Returns a vector where each element is the BIC for the
  corresponding model.
}
\author{C. Fraley}
\seealso{\code{\link{cv1EMtrain}}}
\examples{
data(lansing)
odd <- seq(from=1, to=nrow(lansing), by=2)
round(bicEMtrain(lansing[odd,-3], labels=lansing[odd, 3]), 1)
}
\keyword{cluster}
