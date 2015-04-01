\name{mclust.options}
\alias{mclust.options}

\title{Default values for use with MCLUST package}

\description{Set or retrieve default values for use with MCLUST package.}

\usage{
mclust.options(\dots)
}

\arguments{
  \item{\dots}{
one or more arguments provided in the \code{name = value} form, or no argument at all may be given. \cr
Available arguments are described in the Details section below.}
}
  
\details{
\code{mclust.options} is provided for assigning values to the \code{.mclust}
variable list, which is used to supply default values to various functions 
in \code{MCLUST}. \cr

Available options are:
\describe{
  \item{\code{emModelNames}}{
    A vector of 3-character strings that are associated with multivariate 
    models for which EM estimation is available in MCLUST. \cr
    The current default is all of the multivariate mixture models
    supported in MCLUST.
    The help file for \code{\link{mclustModelNames}} describes the 
    available models.
  }
  \item{\code{hcModelNames}}{
    A vector of character strings associated with multivariate models for
    which model-based hierarchical clustering is available in MCLUST. \cr
    The available models are the following:\cr
    \code{"EII"} = spherical, equal volume \cr
    \code{"EEE"} = ellipsoidal, equal volume, shape, and orientation \cr
    \code{"VII"} = spherical, unequal volume \cr
    \code{"VVV"} = ellipsoidal, varying volume, shape, and orientation. \cr
    The last model in this list is used as default for initialization of 
    EM-algorithm.
  }
  \item{\code{hcUse}}{
    A string specifying the type of input variables to be used for
    model-based hierarchical clustering to start the EM-algorithm.
    Possible values are:\cr
    \code{"VARS"} = original variables (default);\cr
    \code{"STD"} = standardized variables;\cr
    \code{"PCR"} = principal components computed using SVD on standardized 
    (center and scaled) variables (i.e. using the correlation matrix);\cr
    \code{"PCS"} = principal components computed using SVD on centered 
    variables (i.e. using the covariance matrix);\cr
    \code{"SPH"} = sphered variables (centered, scaled, uncorrelated)  
    computed using SVD;\cr
    \code{"SVD"} = scaled SVD transformation;\cr
    \code{"RANDOM"} = random partition.
  }
  \item{\code{bicPlotSymbols}}{
    A vector whose entries correspond to graphics symbols for plotting the 
    BIC values output from \code{\link{Mclust}} and \code{\link{mclustBIC}}. 
    These are displayed in the legend which appears at the lower right
    of the BIC plots.
  }
  \item{\code{bicPlotColors}}{
    A vector whose entries correspond to colors for plotting the 
    BIC curves from output from \code{\link{Mclust}} and
    \code{\link{mclustBIC}}. 
    These are displayed in the legend which appears at the lower right
    of the BIC plots.
  }
  \item{\code{classPlotSymbols}}{
    A vector whose entries are either integers corresponding to graphics 
    symbols or single characters for indicating classifications when
    plotting data. Classes are assigned symbols in the given order. 
  }
  \item{\code{classPlotColors}}{
    A vector whose entries correspond to colors for indicating 
    classifications when plotting data. Classes are assigned colors 
    in the given order. 
  }
  \item{\code{warn}}{
    A logical value indicating whether or not to issue certain warnings.
    Most of these warnings have to do with situations in which 
    singularities are encountered. 
    The default is \code{warn = FALSE}. 
  }
}

The parameter values set via a call to this function will remain in effect for the rest of the session, affecting the subsequent behaviour of the functions for which the given parameters are relevant.
}

\value{
If the argument list is empty the function returns the current list of values.  
If the argument list is not empty, the returned list is invisible.
}

\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of
  Washington.
}

\seealso{
  \code{\link{Mclust}}, 
  \code{\link{MclustDA}}, 
  \code{\link{densityMclust}}, 
  \code{\link{emControl}}
}
\examples{
opt <- mclust.options() # save default values
irisBIC <- mclustBIC(iris[,-5])
summary(irisBIC, iris[,-5])

mclust.options(emModelNames = c("EII", "EEI", "EEE"))
irisBIC <- mclustBIC(iris[,-5])
summary(irisBIC, iris[,-5])

mclust.options(opt)    # restore default values
mclust.options()

oldpar <- par(mfrow = c(2,1), no.readonly = TRUE)
n <- with(mclust.options(), 
          max(sapply(list(bicPlotSymbols, bicPlotColors),length)))
plot(seq(n), rep(1,n), ylab = "", xlab = "", yaxt = "n", 
     pch = mclust.options("bicPlotSymbols"), 
     col = mclust.options("bicPlotColors"))
title("mclust.options(\"bicPlotSymbols\") \n mclust.options(\"bicPlotColors\")")
n <- with(mclust.options(), 
          max(sapply(list(classPlotSymbols, classPlotColors),length)))
plot(seq(n), rep(1,n), ylab = "", xlab = "", yaxt = "n", 
     pch = mclust.options("classPlotSymbols"), 
     col = mclust.options("classPlotColors"))
title("mclust.options(\"classPlotSymbols\") \n mclust.options(\"classPlotColors\")")
par(oldpar)
}
\keyword{cluster}
