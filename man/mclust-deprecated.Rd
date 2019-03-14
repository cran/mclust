\name{mclust-deprecated}
\alias{cv.MclustDA}
\alias{cv1EMtrain}
\alias{bicEMtrain}

\title{Deprecated Functions in mclust package}

\description{
  These functions are provided for compatibility with older versions of the \pkg{mclust}
  package only, and may be removed eventually. 
}

\usage{
cv.MclustDA(\dots)
cv1EMtrain(data, labels, modelNames=NULL)
bicEMtrain(data, labels, modelNames=NULL)
}

\arguments{
  \item{\dots}{pass arguments down.}
  \item{data}{A numeric vector or matrix of observations.}
  \item{labels}{Labels for each element or row in the dataset.}
  \item{modelNames}{Vector of model names that should be tested.
                    The default is to select all available model names.}

}

\seealso{\code{\link{deprecated}}}
