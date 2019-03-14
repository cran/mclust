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
\code{mclust.options} is provided for assigning or retrieving default values used by various functions in \code{MCLUST}.\cr

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
  \item{\code{hcModelName}}{
    A string associated with multivariate models for which model-based 
    hierarchical clustering is available in MCLUST. \cr
    The available models are the following:
    \describe{
    \item{\code{"EII"}}{spherical, equal volume}
    \item{\code{"EEE"}}{ellipsoidal, equal volume, shape, and orientation}
    \item{\code{"VII"}}{spherical, unequal volume}
    \item{\code{"VVV"}}{ellipsoidal, varying volume, shape, and orientation.}
    }
    The \code{"VVV"} is used as default for initialization of EM algorithm.
  }
  \item{\code{hcUse}}{
    A string or a vector of character strings specifying the type of input 
    variables to be used in model-based hierarchical clustering to start 
    the EM algorithm.
    Possible values are:
    \describe{
    \item{\code{"VARS"}}{original variables;}
    \item{\code{"STD"}}{standardized variables;}
    \item{\code{"SPH"}}{sphered variables (centered, scaled, uncorrelated)  
    computed using SVD;}
    \item{\code{"PCS"}}{principal components computed using SVD on centered 
    variables (i.e. using the covariance matrix);}
    \item{\code{"PCR"}}{principal components computed using SVD on standardized 
    (center and scaled) variables (i.e. using the correlation matrix);}
    \item{\code{"SVD"}}{scaled SVD transformation;}
    \item{\code{"RND"}}{no transformation is applied but a random hierarchical structure is returned (see \code{\link{randomPairs}}).}
    }
    For further details see Scrucca and Raftery (2015), Scrucca et al. (2016).
  }
  \item{\code{subset}}{
    A value specifying the maximal sample size to be used in the model-based 
    hierarchical clustering to start the EM algorithm. 
    If data sample size exceeds this value, a random sample is drawn of size
    specified by \code{subset}.
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

\seealso{
  \code{\link{Mclust}}, 
  \code{\link{MclustDA}}, 
  \code{\link{densityMclust}}, 
  \code{\link{emControl}}
}

\references{
Scrucca L. and Raftery A. E. (2015) Improved initialisation of model-based clustering using Gaussian hierarchical partitions. \emph{Advances in Data Analysis and Classification}, 9/4, pp. 447-460.

Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 205-233.
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
