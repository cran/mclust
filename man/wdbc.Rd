\name{wdbc}
\alias{wdbc}
\docType{data}

\title{Wisconsin diagnostic breast cancer (WDBC) data}

\description{
The data set provides data for 569 patients on 30 features of the cell nuclei obtained from a digitized image of a fine needle aspirate (FNA) of a breast mass. For each patient the cancer was diagnosed as malignant or benign.}

\usage{data(wdbc)}

\format{A data frame with 569 observations on the following variables:
  \describe{
    \item{\code{ID}}{ID number}
    \item{\code{Diagnosis}}{cancer diagnosis: \code{M} = malignant, \code{B} = benign}
    \item{\code{Radius_mean}}{a numeric vector}
    \item{\code{Texture_mean}}{a numeric vector}
    \item{\code{Perimeter_mean}}{a numeric vector}
    \item{\code{Area_mean}}{a numeric vector}
    \item{\code{Smoothness_mean}}{a numeric vector}
    \item{\code{Compactness_mean}}{a numeric vector}
    \item{\code{Concavity_mean}}{a numeric vector}
    \item{\code{Nconcave_mean}}{a numeric vector}
    \item{\code{Symmetry_mean}}{a numeric vector}
    \item{\code{Fractaldim_mean}}{a numeric vector}
    \item{\code{Radius_se}}{a numeric vector}
    \item{\code{Texture_se}}{a numeric vector}
    \item{\code{Perimeter_se}}{a numeric vector}
    \item{\code{Area_se}}{a numeric vector}
    \item{\code{Smoothness_se}}{a numeric vector}
    \item{\code{Compactness_se}}{a numeric vector}
    \item{\code{Concavity_se}}{a numeric vector}
    \item{\code{Nconcave_se}}{a numeric vector}
    \item{\code{Symmetry_se}}{a numeric vector}
    \item{\code{Fractaldim_se}}{a numeric vector}
    \item{\code{Radius_extreme}}{a numeric vector}
    \item{\code{Texture_extreme}}{a numeric vector}
    \item{\code{Perimeter_extreme}}{a numeric vector}
    \item{\code{Area_extreme}}{a numeric vector}
    \item{\code{Smoothness_extreme}}{a numeric vector}
    \item{\code{Compactness_extreme}}{a numeric vector}
    \item{\code{Concavity_extreme}}{a numeric vector}
    \item{\code{Nconcave_extreme}}{a numeric vector}
    \item{\code{Symmetry_extreme}}{a numeric vector}
    \item{\code{Fractaldim_extreme}}{a numeric vector}
  }
}

\details{
The recorded features are:
\itemize{
\item \code{Radius} as mean of distances from center to points on the perimeter
\item \code{Texture} as standard deviation of gray-scale values
\item \code{Perimeter} as cell nucleus perimeter
\item \code{Area} as cell nucleus area
\item \code{Smoothness} as local variation in radius lengths
\item \code{Compactness} as cell nucleus compactness, perimeter^2 / area - 1
\item \code{Concavity} as severity of concave portions of the contour
\item \code{Nconcave} as number of concave portions of the contour
\item \code{Symmetry} as cell nucleus shape
\item \code{Fractaldim} as fractal dimension, "coastline approximation" - 1
}
For each feature the recorded values are computed from each image as \code{<feature_name>_mean}, \code{<feature_name>_se}, and \code{<feature_name>_extreme}, for the mean, the standard error, and the mean of the three largest values.
}

\source{Breast Cancer Wisconsin (Diagnostic) Data Set (\code{wdbc.data}, \code{wdbc.names}) is available at UCI Machine Learning Repository
\url{https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+(Diagnostic)}. Please note the UCI conditions of use.}

\references{
Mangasarian, O. L., Street, W. N., and Wolberg, W. H. (1995) Breast cancer diagnosis and prognosis via linear programming. \emph{Operations Research}, 43(4), pp. 570-577.
}

\keyword{datasets}
