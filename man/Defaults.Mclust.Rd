\name{Defaults.Mclust}
\alias{Defaults.Mclust}
\alias{.Mclust}
\title{
  List of values controlling defaults for some MCLUST functions.
}
\description{
  A named list of values including an enumeration of models used as defaults
  in MCLUST functions. 
}
\value{
  A list with the following components:

 \item{emModelNames}{
    A vector of character strings associated with multivariate models for
    which EM estimation is available in MCLUST. \cr
    The current default is the following list:\cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEI": diagonal, equal volume and shape\cr 
    "VEI": diagonal, varying volume, equal shape\cr 
    "EVI": diagonal, equal volume, varying shape \cr
    "VVI": diagonal, varying volume and shape \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "EEV": ellipsoidal, equal volume and equal shape\cr
    "VEV": ellipsoidal, equal shape \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation 
  }
  \item{hcModelNames}{
    A vector of character strings associated with multivariate models for
    which model-based hierarchical clustering is available in MCLUST. \cr
    The current default is the following list:\cr\cr
    "EII": spherical, equal volume \cr
    "VII": spherical, unequal volume \cr
    "EEE": ellipsoidal, equal volume, shape, and orientation \cr
    "VVV": ellipsoidal, varying volume, shape, and orientation 
  }
  \item{bicPlotSymbols}{
    A vector whose entries correspond to graphics symbols for plotting the 
    BIC values output from \code{Mclust} and \code{mclustBIC}. 
    These are displayed in the legend which appears at the lower right
    of the BIC plots.
  }
  \item{bicPlotColors}{
    A vector whose entries correspond to colors for plotting the 
    BIC curves from output from \code{Mclust} and \code{mclustBIC}. 
    These are displayed in the legend which appears at the lower right
    of the BIC plots.
  }
  \item{classPlotSymbols}{
    A vector whose entries are either integers corresponding to graphics 
    symbols or single characters for indicating classifications when
    plotting data. Classes are assigned symbols in the given order. 
  }
  \item{classPlotColors}{
    A vector whose entries correspond to colors for indicating 
    classifications when plotting data. Classes are assigned colors 
    in the given order. 
  }
\item{warn}{
    A logical value indicating whether or not to issue certain warnings
    (usually involving singularity). Default: \code{warn = TRUE}.
  }
}
\details{
  A function \code{mclustOptions} is supplied for assigning values to
  the \code{.Mclust} list. 
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association}. 
  
  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling  and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics, University of Washington. 
}
\seealso{
  \code{\link{mclustOptions}},
  \code{\link{Mclust}},
  \code{\link{mclustBIC}}
}
\examples{
irisBIC <- Mclust(iris[,-5])
summary(irisBIC, iris[-5])

.Mclust
.Mclust <- mclustOptions(emModelNames = c("VII", "VVI", "VVV"))
.Mclust

irisBIC <- Mclust(iris[,-5])
summary(irisBIC, iris[-5])

.Mclust <- mclustOptions() # restore defaults
.Mclust
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
