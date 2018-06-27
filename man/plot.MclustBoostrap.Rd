\name{plot.MclustBootstrap}
\alias{plot.MclustBootstrap}

\title{Plot of bootstrap distributions for mixture model parameters}

\description{
Plots the bootstrap distribution of parameters as returned by the \code{\link{MclustBootstrap}} function.
}

\usage{
\method{plot}{MclustBootstrap}(x, what = c("pro", "mean", "var"), 
     show.parest = TRUE, show.confint = FALSE,
     hist.col = "grey", hist.border = "lightgrey", breaks = "Sturges", 
     col = "forestgreen", lwd = 2, lty = 3, 
     xlab = NULL, xlim = NULL, ylim = NULL, \dots)
}

\arguments{
  \item{x}{Object returned by \code{MclustBootstrap}.}
  
  \item{what}{Character string specifying if mixing proportions (\code{"pro"}),
    component means (\code{"mean"}) or component variances (\code{"var"}) 
    should be drawn.}
    
  \item{show.parest}{A logical specifying if the parameter estimate should be drawn as vertical line.}
    
  \item{show.confint}{A logical specifying if resampling-based confidence interval should be drawn at the bottom of the graph. Confidence level can be provided as further argument \code{conf.level}; see \code{\link{summary.MclustBootstrap}}.}

  \item{hist.col}{The color to be used to fill the bars of the histograms.}

  \item{hist.border}{The color of the border around the bars of the histograms.}
  
  \item{breaks}{See the argument in function \code{\link[graphics]{hist}}.}

  \item{col, lwd, lty}{The color, line width and line type to be used to represent the estimated parameters and confidence intervals.}

 \item{xlab}{Optional label for the horizontal axis.}

 \item{xlim, ylim}{A two-values vector of axis range for, respectively, horizontal and
   vertical axis.}
   
  \item{\dots}{Other graphics parameters.}
}

\value{
  A plot for each variable/component of the selected parameters.
}

\seealso{
  \code{\link{MclustBootstrap}}
}
\examples{
\dontrun{
data(diabetes)
X <- diabetes[,-1]
modClust <- Mclust(X, G = 3, modelNames = "VVV")
bootClust <- MclustBootstrap(modClust, nboot = 99)
par(mfrow = c(1,3), mar = c(4,2,2,0.5))
plot(bootClust, what = "pro")
par(mfrow = c(3,3), mar = c(4,2,2,0.5))
plot(bootClust, what = "mean")
}
}

\keyword{cluster}
