\name{clPairs}
\alias{clPairs}
\title{
  Pairwise Scatter Plots showing Classification
}
\description{
  Creates a scatter plot for each pair of variables in given data.
  Observations in different classes are represented by different symbols.
}
\usage{
clPairs(data, classification, symbols, labels=dimnames(data)[[2]], 
        CEX=1, col, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{classification}{
    A numeric or character vector representing a classification of observations (rows) of \code{data}.
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique class \code{classification}. Elements in \code{symbols}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function
    \code{unique}). Default: If \emph{G} is the number of groups in the
    classification, the first \emph{G} symbols in
    \code{.Mclust\$symbols}, otherwise if \emph{G} is less than 27 then
    the first \emph{G} capital letters in the Roman alphabet. If no
    \code{classification} argument is given the default symbol is
    \code{"."}.
  }
  \item{labels}{
    A vector of character strings for labeling the variables. The default
    is to use the column dimension names of \code{data}.
  }
  \item{CEX}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{col}{
    Color vector to use. Default is one color per class. Splus default:
    all black.
  }
  \item{\dots}{
    Additional arguments to be passed to the graphics device.
  }
}
\section{Side Effects}{
  Scatter plots for each combination of variables in \code{data} are
  created on the current graphics device. Observations of different
  classifications are labeled with different symbols. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002b).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis. 
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{pairs}},
  \code{\link{coordProj}},
  \code{\link{mclustOptions}}
}
\examples{
data(iris)
irisMatrix <- as.matrix(iris[,1:4])
irisClass <- iris[,5]

clPairs(irisMatrix, cl=irisClass, symbols=as.character(1:3))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
