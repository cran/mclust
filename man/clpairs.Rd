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
clPairs(data, classification, symbols, colors, labels=dimnames(data)[[2]], 
        CEX=1, \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{classification}{
    A numeric or character vector representing a classification of observations
   (rows) of \code{data}.
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique class in \code{classification}. Elements in \code{symbols}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{unique}). 
    The default is given is \code{.Mclust\$classPlotSymbols}.
  }
  \item{colors}{
    Either an integer or character vector assigning a color to each
    unique class in \code{classification}. Elements in \code{colors}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{unique}). 
    The default is given is \code{.Mclust\$classPlotColors}.
  }
  \item{labels}{
    A vector of character strings for labeling the variables. The default
    is to use the column dimension names of \code{data}.
  }
  \item{CEX}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
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
  C. Fraley and A. E. Raftery (2006, revised 2010).
  MCLUST Version 3 for R: Normal Mixture Modeling
  and Model-Based Clustering, Technical Report no. 504, 
  Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{pairs}},
  \code{\link{coordProj}},
  \code{\link{mclustOptions}}
}
\examples{
clPairs(iris[,-5], cl=iris[,5], symbols=as.character(1:3))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
