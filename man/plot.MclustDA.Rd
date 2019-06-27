\name{plot.MclustDA}
\alias{plot.MclustDA}

\title{Plotting method for MclustDA discriminant analysis}

\description{
Plots for model-based mixture discriminant analysis results, such as scatterplot of training and test data, classification of train and test data, and errors.
}

\usage{
\method{plot}{MclustDA}(x, what = c("scatterplot", "classification", "train&test", "error"), 
     newdata, newclass, dimens = NULL, 
     symbols, colors, main = NULL, \dots)
}

\arguments{
  \item{x}{
    An object of class \code{'MclustDA'} resulting from a call to \code{\link{MclustDA}}.
  }

 \item{what}{
    A string specifying the type of graph requested. Available choices are:
    \describe{
      \item{\code{"scatterplot"} =}{a plot of training data with points marked based on the known classification. Ellipses corresponding to covariances of mixture components are also drawn.}
      \item{\code{"classification"} =}{a plot of data with points marked on based the predicted classification; if \code{newdata} is provided then the test set is shown otherwise the training set.}
      \item{\code{"train&test"} =}{a plot of training and test data with points marked according to the type of set.}
      \item{\code{"error"} =}{a plot of training set (or test set if \code{newdata} and \code{newclass} are provided) with misclassified points marked.}
   }
   If not specified, in interactive sessions a menu of choices is proposed.
  }

  \item{newdata}{
    A data frame or matrix for test data.
  }

  \item{newclass}{
    A vector giving the class labels for the observations in 
    the test data (if known).
  }

  \item{dimens}{
    A vector of integers giving the dimensions of the desired coordinate
    projections for multivariate data. The default is to take all the
    the available dimensions for plotting.
  }

  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique class. Elements in \code{colors} correspond to classes in order of
    appearance in the sequence of observations (the order used by the 
    function \code{factor}). 
    The default is given by \code{mclust.options("classPlotSymbols")}.
  }

  \item{colors}{
    Either an integer or character vector assigning a color to each
    unique class in \code{classification}. Elements in \code{colors}
    correspond to classes in order of appearance in the sequence of
    observations (the order used by the function \code{factor}). 
    The default is given by \code{mclust.options("classPlotColors")}.
  }

  \item{main}{
    A logical, a character string, or \code{NULL} (default) for the main title. 
    If \code{NULL} or \code{FALSE} no title is added to a plot. 
    If \code{TRUE} a default title is added identifying the type of plot drawn.
    If a character string is provided, this is used for the title.
  }
  
  \item{\dots}{further arguments passed to or from other methods.}
}

%\value{}
  
\details{
 For more flexibility in plotting, use \code{mclust1Dplot}, 
 \code{mclust2Dplot}, \code{surfacePlot}, \code{coordProj}, or
 \code{randProj}. 
}

\author{Luca Scrucca}

\seealso{
  \code{\link{MclustDA}},
  \code{\link{surfacePlot}},
  \code{\link{coordProj}},
  \code{\link{randProj}}
}

\examples{
\dontrun{
odd <- seq(from = 1, to = nrow(iris), by = 2)
even <- odd + 1
X.train <- iris[odd,-5]
Class.train <- iris[odd,5]
X.test <- iris[even,-5]
Class.test <- iris[even,5]

# common EEE covariance structure (which is essentially equivalent to linear discriminant analysis)
irisMclustDA <- MclustDA(X.train, Class.train, modelType = "EDDA", modelNames = "EEE")
summary(irisMclustDA, parameters = TRUE)
summary(irisMclustDA, newdata = X.test, newclass = Class.test)

# common covariance structure selected by BIC
irisMclustDA <- MclustDA(X.train, Class.train, modelType = "EDDA")
summary(irisMclustDA, parameters = TRUE)
summary(irisMclustDA, newdata = X.test, newclass = Class.test)

# general covariance structure selected by BIC
irisMclustDA <- MclustDA(X.train, Class.train)
summary(irisMclustDA, parameters = TRUE)
summary(irisMclustDA, newdata = X.test, newclass = Class.test)

plot(irisMclustDA)
plot(irisMclustDA, dimens = 3:4)
plot(irisMclustDA, dimens = 4)

plot(irisMclustDA, what = "classification")
plot(irisMclustDA, what = "classification", newdata = X.test)
plot(irisMclustDA, what = "classification", dimens = 3:4)
plot(irisMclustDA, what = "classification", newdata = X.test, dimens = 3:4)
plot(irisMclustDA, what = "classification", dimens = 4)
plot(irisMclustDA, what = "classification", dimens = 4, newdata = X.test)

plot(irisMclustDA, what = "train&test", newdata = X.test)
plot(irisMclustDA, what = "train&test", newdata = X.test, dimens = 3:4)
plot(irisMclustDA, what = "train&test", newdata = X.test, dimens = 4)

plot(irisMclustDA, what = "error")
plot(irisMclustDA, what = "error", dimens = 3:4)
plot(irisMclustDA, what = "error", dimens = 4)
plot(irisMclustDA, what = "error", newdata = X.test, newclass = Class.test)
plot(irisMclustDA, what = "error", newdata = X.test, newclass = Class.test, dimens = 3:4)
plot(irisMclustDA, what = "error", newdata = X.test, newclass = Class.test, dimens = 4)

# simulated 1D data
n <- 250 
set.seed(1)
triModal <- c(rnorm(n,-5), rnorm(n,0), rnorm(n,5))
triClass <- c(rep(1,n), rep(2,n), rep(3,n))
odd <- seq(from = 1, to = length(triModal), by = 2)
even <- odd + 1
triMclustDA <- MclustDA(triModal[odd], triClass[odd])
summary(triMclustDA, parameters = TRUE)
summary(triMclustDA, newdata = triModal[even], newclass = triClass[even])
plot(triMclustDA)
plot(triMclustDA, what = "classification")
plot(triMclustDA, what = "classification", newdata = triModal[even])
plot(triMclustDA, what = "train&test", newdata = triModal[even])
plot(triMclustDA, what = "error")
plot(triMclustDA, what = "error", newdata = triModal[even], newclass = triClass[even])

# simulated 2D cross data
data(cross)
odd <- seq(from = 1, to = nrow(cross), by = 2)
even <- odd + 1
crossMclustDA <- MclustDA(cross[odd,-1], cross[odd,1])
summary(crossMclustDA, parameters = TRUE)
summary(crossMclustDA, newdata = cross[even,-1], newclass = cross[even,1])
plot(crossMclustDA)
plot(crossMclustDA, what = "classification")
plot(crossMclustDA, what = "classification", newdata = cross[even,-1])
plot(crossMclustDA, what = "train&test", newdata = cross[even,-1])
plot(crossMclustDA, what = "error")
plot(crossMclustDA, what = "error", newdata =cross[even,-1], newclass = cross[even,1])
}
}

\keyword{multivariate}
