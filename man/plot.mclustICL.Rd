\name{plot.mclustICL}
\alias{plot.mclustICL}

\title{ICL Plot for Model-Based Clustering}

\description{
Plots the ICL values returned by the \code{\link{mclustICL}} function.
}

\usage{
\method{plot}{mclustICL}(x, ylab = "ICL", \dots)     
}

\arguments{
  \item{x}{
    Output from \code{\link{mclustICL}}.
  }

 \item{ylab}{
    Label for the vertical axis of the plot.
  }

  \item{\dots}{
    Further arguments passed to the \code{\link{plot.mclustBIC}} function.
  }
}

\value{
  A plot of the ICL values.
}

\seealso{
  \code{\link{mclustICL}}
}
\examples{
\dontrun{
data(faithful)
faithful.ICL = mclustICL(faithful)
plot(faithful.ICL)
}
}

\keyword{cluster}
% docclass is function
