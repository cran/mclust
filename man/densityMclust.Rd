\name{densityMclust}
\alias{densityMclust}

\title{Density Estimation via Model-Based Clustering}

\description{
  Produces a density estimate for each data point using a Gaussian finite 
  mixture model from \code{Mclust}.
}

\usage{
densityMclust(data, \dots)
}

\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{\dots }{
    Additional arguments for the \code{\link{Mclust}} function. 
    In particular, setting the arguments \code{G} and \code{modelNames} 
    allow to specify the number of mixture components and the type of
    model to be fitted. By default an "optimal" model is selected based
    on the BIC criterion. 
  }
}

\value{
  An object of class \code{densityMclust}, which inherits from 
  \code{Mclust}, is returned with the following slot added:
  \item{density}{The density evaluated at the input \code{data}
    computed from the estimated model.}
}

%\details{}
  
\references{
Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 205-233. 

Fraley C. and Raftery A. E. (2002) Model-based clustering, discriminant analysis and density estimation, \emph{Journal of the American Statistical Association}, 97/458, pp. 611-631.

Fraley C., Raftery A. E., Murphy T. B. and Scrucca L. (2012) mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation. \emph{Technical Report} No. 597, Department of Statistics, University of Washington.
}

\author{Revised version by Luca Scrucca based on 
  the original code by C. Fraley and A.E. Raftery.}

\seealso{
  \code{\link{plot.densityMclust}}, 
  \code{\link{Mclust}}, 
  \code{\link{summary.Mclust}},
  \code{\link{predict.densityMclust}}.
}

\examples{
dens <- densityMclust(faithful$waiting)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "BIC", legendArgs = list(x = "topright"))
plot(dens, what = "density", data = faithful$waiting)

dens <- densityMclust(faithful, modelNames = "EEE", G = 3)
summary(dens)
summary(dens, parameters = TRUE)
plot(dens, what = "density", data = faithful, 
     drawlabels = FALSE, points.pch = 20)
plot(dens, what = "density", type = "hdr")
plot(dens, what = "density", type = "hdr", prob = c(0.1, 0.9))
plot(dens, what = "density", type = "hdr", data = faithful)
plot(dens, what = "density", type = "persp")

\dontrun{
dens <- densityMclust(iris[,1:4], G = 2)
summary(dens, parameters = TRUE)
plot(dens, what = "density", data = iris[,1:4], 
     col = "slategrey", drawlabels = FALSE, nlevels = 7)
plot(dens, what = "density", type = "hdr", data = iris[,1:4])
plot(dens, what = "density", type = "persp", col = grey(0.9))
}
}

\keyword{cluster}
