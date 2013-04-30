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

\section{References}{
  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
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
