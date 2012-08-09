\name{cv1EMtrain}
\alias{cv1EMtrain}
\title{Select discriminant models using cross validation}
\description{
Leave-one-out cross validation given a dataset and labels for selected models.
}
\usage{
cv1EMtrain(data, labels, modelNames=NULL)
}
\arguments{
  \item{data}{A numeric vector or matrix of observations.}
  \item{labels}{Labels for each element or row in the dataset.}
  \item{modelNames}{Vector of model names that should be tested.
            The default is to select all available model names.}
}
\value{
  Returns a vector where each element is the the crossvalidated error rate
  for the dataset and labels corresponding to each model.
}
\section{References}{
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\author{C. Fraley}
\seealso{\code{\link{bicEMtrain}}}
\examples{
even <- seq(from=2, to=nrow(chickwts), by=2)
round(cv1EMtrain(chickwts[even,1], labels=chickwts[even,2]), 1)
}
\keyword{cluster}

