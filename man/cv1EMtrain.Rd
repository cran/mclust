\name{cv1EMtrain}
\alias{cv1EMtrain}
\title{Select discriminant models using cross validation}
\description{
For the ten available discriminant models the leave-one-out cross
validation error is calulated. The models for one-dimensional data are
"E" and "V"; for higher dimensions they are "EII", "VII", "EEI",
"VEI", "EVI", "VVI", "EEE", "EEV", "VEV" and "VVV".
}
\usage{
cv1EMtrain(data, labels, modelNames)
}
\arguments{
  \item{data}{A data matrix}
  \item{labels}{Labels for each row in the data matrix}
  \item{modelNames}{Vector of model names that should be tested.}
}
\value{
  Returns a vector where each element is the error rate for the
  corresponding model.
}
\author{C. Fraley}
\seealso{\code{\link{bicEMtrain}}}
\examples{
data(lansing)
odd <- seq(from=1, to=nrow(lansing), by=2)
round(cv1EMtrain(data=lansing[odd,-3], labels=lansing[odd,3]), 3)

cv1Modd <- mstepEEV(data=lansing[odd,-3], z=unmap(lansing[odd,3]))
cv1Zodd <- do.call("estepEEV", c(cv1Modd, list(data=lansing[odd,-3])))$z
compClass(map(cv1Zodd), lansing[odd,3])

even <- (1:nrow(lansing))[-odd]
cv1Zeven <- do.call("estepEEV", c(cv1Modd, list(data=lansing[even,-3])))$z
compClass(map(cv1Zodd), lansing[odd,3])$error
}
\keyword{cluster}

