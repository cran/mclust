\name{mclust1Dplot}
\alias{mclust1Dplot}
\title{
  Plot one-dimensional data modeled by an MVN mixture.
}
\description{
  Plot one-dimensional data given parameters of an MVN mixture model 
  for the data.
}
\usage{
mclust1Dplot(data, parameters = NULL, z = NULL, 
             classification = NULL, truth = NULL, uncertainty = NULL, 
             what = c("classification", "density", "error", "uncertainty"),
             symbols = NULL, colors = NULL, ngrid = length(data), 
             xlab = NULL, ylab = NULL, 
             xlim = NULL, ylim = NULL,
             cex = 1, main = FALSE, \dots)
}
\arguments{
  \item{data}{
    A numeric vector of observations.
    Categorical variables are not allowed.
  }
 \item{parameters}{
     A named list giving the parameters of an \emph{MCLUST} model,
     used to produce superimposing ellipses on the plot.
     The relevant components are as follows:
     \describe{
        \item{\code{pro}}{
              Mixing proportions for the components of the mixture. 
              There should one more mixing proportion than the number of 
              Gaussian components if the mixture model includes 
              a Poisson noise term.
        }
        \item{\code{mean}}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{\code{variance}}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{\link{mclustVariance}}
              for details.
        }
     }
  }
\item{z}{
        A matrix in which the \code{[i,k]}th entry gives the
        probability of observation \emph{i} belonging to the \emph{k}th class.
        Used to compute \code{classification} and
        \code{uncertainty} if those arguments aren't available.
  }
\item{classification}{
        A numeric or character vector representing a classification of
        observations (rows) of \code{data}. If present argument \code{z}
        will be ignored.
  }
  \item{truth}{
        A numeric or character vector giving a known
        classification of each data point.
        If \code{classification} or \code{z} is also present,
        this is used for displaying classification errors.
  }
  \item{uncertainty}{
        A numeric vector of values in \emph{(0,1)} giving the
        uncertainty of each data point. If present argument \code{z}
        will be ignored.
  }
  \item{what}{
    Choose from one of the following options: \code{"classification"}
    (default), \code{"density"}, \code{"error"}, \code{"uncertainty"}.
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to
    each unique class \code{classification}. Elements in \code{symbols}
    correspond to classes in \code{classification} in order of
    appearance in the observations (the order used by the 
    function \code{unique}). The default is to use a single plotting
    symbol \emph{|}. Classes are delineated by showing them in separate
    lines above the whole of the data.
  }
  \item{colors}{
    Either an integer or character vector assigning a color to each
    unique class \code{classification}. Elements in \code{colors}
    correspond to classes in order of appearance in the observations 
    (the order used by the function \code{unique}).
    The default is given is \code{mclust.options("classPlotColors")}.
  }
  \item{ngrid}{
    Number of grid points to use for density computation over the interval
    spanned by the data. The default is the length of the data set.
  }
  \item{xlab, ylab}{
    An argument specifying a label for the axes.
  }
  \item{xlim, ylim}{
    An argument specifying bounds of the plot.
    This may be useful for when comparing plots.
  }
  \item{cex}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{main}{
    A logical variable or \code{NULL} indicating whether or not to add a title
    to the plot identifying the dimensions used.
  }
\item{\dots}{
    Other graphics parameters.
  }
}
\value{
  A plot showing location of the mixture components, classification, uncertainty, density and/or classification errors. Points in the different classes are shown in separated levels above the whole of the data.
}

\seealso{
  \code{\link{mclust2Dplot}},
  \code{\link{clPairs}},
  \code{\link{coordProj}}
}
\examples{
\dontrun{
n <- 250 ## create artificial data
set.seed(1)
y <- c(rnorm(n,-5), rnorm(n,0), rnorm(n,5))
yclass <- c(rep(1,n), rep(2,n), rep(3,n))

yModel <- Mclust(y)

mclust1Dplot(y, parameters = yModel$parameters, z = yModel$z, 
             what = "classification")

mclust1Dplot(y, parameters = yModel$parameters, z = yModel$z, 
             what = "error", truth = yclass)

mclust1Dplot(y, parameters = yModel$parameters, z = yModel$z, 
             what = "density")

mclust1Dplot(y, z = yModel$z, parameters = yModel$parameters,
            what = "uncertainty")

}
}
\keyword{cluster}

