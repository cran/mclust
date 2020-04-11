\name{plot.mclustBIC}
\alias{plot.mclustBIC}

\title{BIC Plot for Model-Based Clustering}

\description{
Plots the BIC values returned by the \code{\link{mclustBIC}} function.
}

\usage{
\method{plot}{mclustBIC}(x, G = NULL, modelNames = NULL, 
     symbols = NULL, colors = NULL, 
     xlab = NULL, ylab = "BIC", 
     legendArgs = list(x = "bottomright", ncol = 2, cex = 1, inset = 0.01), 
     \dots)
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
 \item{ylab}{
    Label for the vertical axis of the BIC plot.
  }
 \item{legendArgs}{
    Arguments to pass to the \code{legend} function. Set to \code{NULL}
    for no legend.
  }
  \item{\dots}{
    Other graphics parameters.
  }
}
\value{
  A plot of the BIC values.
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
