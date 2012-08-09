\name{plot.mclustBIC}
\alias{plot.mclustBIC}
\title{
  BIC Plot
}
\description{
  Plots the BIC from mclust modeling via function \code{mclustBIC}.
}
\usage{
\method{plot}{mclustBIC}(x, G = NULL, modelNames = NULL, symbols = NULL,
                         colors = NULL, xlab = NULL, ylim = NULL,
                         legendArgs = list(x="bottomright", ncol=2, cex=1), 
                         CEX = 1, \dots)
}
\arguments{
  \item{x}{
    Output from \code{mclustBIC}.
  }
 \item{G}{
    One or more numbers of components corresponding to models fit in \code{x}.
    The default is to plot the BIC for all of the numbers of components fit.
  }
 \item{modelNames}{
    One or more model names corresponding to models fit in \code{x}.
    The default is to plot the BIC for all of the models fit.
  }
 \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique class in \code{classification}. Elements in \code{colors}
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
 \item{xlab}{
    Optional label for the horizontal axis of the BIC plot.
  }
 \item{ylim}{
    Optional limits for the vertical axis of the BIC plot.
  }
 \item{legendArgs}{
    Arguments to pass to the \code{legend} function. Set to \code{NULL}
    for no legend.
  }
  \item{CEX}{
    A scalar controling the size of the splot symbols.
  }
  \item{\dots}{
    Other graphics parameters.
  }
}
\value{
  A plot of the BIC values for the models specified in the 
 \code{modelNames} argument.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}

\seealso{
  \code{\link{mclustBIC}}
}
\examples{
\dontrun{
plot(mclustBIC(precip), legendArgs =  list(x = "bottomleft"))

plot(mclustBIC(faithful))

plot(mclustBIC(iris[,-5]))
}
}
\keyword{cluster}
% docclass is function
