\name{clPairs}
\alias{clPairs}
\alias{clPairsLegend}

\title{Pairwise Scatter Plots showing Classification}

\description{
Creates a scatter plot for each pair of variables in given data.
Observations in different classes are represented by different colors and symbols.
}

\usage{
clPairs(data, classification, symbols, colors, 
        labels = dimnames(data)[[2]], cex.labels = 1.5, 
        gap = 0.2, \dots)

clPairsLegend(x, y, class, col, pch, box = TRUE, \dots)
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
    The default is given by \code{mclust.options("classPlotSymbols")}.
  }
  \item{colors}{
    Either an integer or character vector assigning a color to each
    unique class in \code{classification}. Elements in \code{colors}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{unique}). 
    The default is given by \code{mclust.options("classPlotColors")}.
  }
  \item{labels}{
    A vector of character strings for labelling the variables. The default
    is to use the column dimension names of \code{data}.
  }
  \item{cex.labels}{
    An argument specifying the size of the text labels.
  }
  \item{gap}{
    An argument specifying the distance between subplots (see \code{\link{pairs}}).
  }
  \item{x,y}{
    The x and y co-ordinates with respect to a graphic device having 
    plotting region coordinates \code{par("usr" = c(0,1,0,1))}.
  }
  \item{class}{
    The class labels.
  }
  \item{box}{
    A logical, if \code{TRUE} then a box is drawn around the current plot figure. 
  }
  \item{col, pch}{
    The colors and plotting symbols appearing in the legend.
  }
  \item{\dots}{
    For a \code{clPairs} call may be additional arguments to be passed to \code{\link{pairs}}. 
    For a \code{clPairsLegend} call may be additional arguments to be passed to \code{\link{legend}}.
    }
}

\details{
  The function \code{clPairs()} draws scatter plots on the current graphics device for each combination of variables in \code{data}. Observations of different classifications are labeled with different symbols. 
  
  The function \code{clPairsLegend()} can be used to add a legend. See examples below.
}

\value{
  The function \code{clPairs()} invisibly returns a list with the following components:
  \item{class}{A character vector of class labels.}
  \item{col}{A vector of colors used for each class.}
  \item{pch}{A vector of plotting symbols used for each class.}
}

\seealso{
  \code{\link{pairs}},
  \code{\link{coordProj}},
  \code{\link{mclust.options}}
}
\examples{
clPairs(iris[,1:4], cl = iris$Species)

clp <- clPairs(iris[,1:4], cl = iris$Species, lower.panel = NULL)
clPairsLegend(0.1, 0.4, class = clp$class, 
              col = clp$col, pch = clp$pch, 
              title = "Iris data")

}
\keyword{cluster}
