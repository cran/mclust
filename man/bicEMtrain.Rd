\name{bicEMtrain}
\alias{bicEMtrain}
\title{Select models in discriminant analysis using BIC}
\description{
Computes the BIC given a dataset and labels for selected models.
}
\usage{
bicEMtrain(data, labels, modelNames=NULL)
}
\arguments{
  \item{data}{A numeric vector or matrix of observations.}
  \item{labels}{Labels for each element or row in the data.}
  \item{modelNames}{Vector of model names that should be tested.
            The default is to select all available model names.}
}
\value{
  Returns a vector where each element is the BIC for the
  dataset and labels corresponding to each model.
}
\section{References}{
  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling
  and Model-Based Clustering, Technical Report no. 504, 
  Department of Statistics, University of Washington.
}
\author{C. Fraley}
\seealso{\code{\link{cv1EMtrain}}}
\examples{
even <- seq(from=2, to=nrow(chickwts), by=2)
round(bicEMtrain(chickwts[even,1], labels=chickwts[even,2]), 1)
}
\keyword{cluster}
