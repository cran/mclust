\name{mclust2Dplot}
\alias{mclust2Dplot}
\title{
  Plot two-dimensional data modelled by an MVN mixture.
}
\description{
  Plot two-dimensional data given parameters of an MVN mixture model 
  for the data.
}
\usage{
mclust2Dplot(data, parameters=NULL, z=NULL,
             classification=NULL, truth=NULL, uncertainty=NULL,
             what = c("classification","uncertainty","errors"), 
             quantiles = c(0.75, 0.95), symbols=NULL, scale = FALSE,
             xlim=NULL, ylim=NULL, CEX = 1, PCH = ".", identify = FALSE, 
             swapAxes = FALSE, \dots) 
}
\arguments{
  \item{data}{
    A numeric matrix or data frame of observations.
    Categorical variables are not allowed.
    If a matrix or data frame, rows correspond to observations and
    columns correspond to variables. 
    In this case the data are two dimensional, so there are two columns.
  }
  \item{parameters}{
     A named list giving the parameters of an \emph{MCLUST} model, 
     used to produce superimposing ellipses on the plot. 
     The relevant components are as follows:
     \describe{
        \item{mean}{
              The mean for each component. If there is more than one component,
              this is a matrix whose kth column is the mean of the \emph{k}th
              component of the mixture model.
        }
        \item{variance}{
              A list of variance parameters for the model.
              The components of this list depend on the model
              specification. See the help file for \code{mclustVariance}
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
        If \code{classification}
        or \code{z} is also present, 
        this is used for displaying classification errors.
  }
  \item{uncertainty}{
        A numeric vector of values in \emph{(0,1)} giving the
        uncertainty of each data point. If present argument \code{z}
        will be ignored.
  }
  \item{what}{
    Choose from one of the following three options: \code{"classification"}
    (default), \code{"errors"}, \code{"uncertainty"}. 
  }
  \item{quantiles}{
    A vector of length 2 giving quantiles used in plotting uncertainty.
    The smallest symbols correspond to the smallest quantile (lowest 
    uncertainty), medium-sized (open) symbols to points falling between 
    the given quantiles, and large (filled) symbols to those in the largest 
    quantile (highest uncertainty).
    The default is \emph{(0.75,0.95)}.
  }
  \item{symbols}{
    Either an integer or character vector assigning a plotting symbol to each
    unique class \code{classification}. 
    Elements in \code{symbols}
    correspond to classes in \code{classification}
    in order of appearance in the observations
    (the order used by the S-PLUS function 
    \code{unique}). 
    Default: If \emph{G} is the number of groups in the classification, 
    the first \emph{G} symbols in \code{.Mclust\$symbols},
    otherwise if \emph{G} is less than 27 then the first \emph{G} capital
    letters in the Roman alphabet.
  }
  \item{scale}{
    A logical variable indicating whether or not the two chosen
    dimensions should be plotted on the same scale, and
    thus preserve the shape of the distribution.
    Default: \code{scale=FALSE} 
  }
  \item{xlim, ylim}{
    An argument specifying bounds for the ordinate, abscissa of the plot.
    This may be useful for when comparing plots.
  }
  \item{CEX}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{PCH}{
    An argument specifying the symbol to be used when a classificatiion
    has not been specified for the data. The default value is a small dot ".".
  }
  \item{identify}{
    A logical variable indicating whether or not to add a title to the plot
    identifying the dimensions used.
  }
  \item{swapAxes}{
    A logical variable indicating whether or not the axes should be swapped
    for the plot.
  }
  \item{\dots}{
    Other graphics parameters.
  }
}
\section{Side Effects}{
  A plot showing the data, together with the location of the mixture 
  components, classification, uncertainty, and/or classification errors.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3: An R Package for Normal Mixture Modeling
  and Model-Based Clustering, Technical Report, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{surfacePlot}},
  \code{\link{clPairs}},
  \code{\link{coordProj}},
  \code{\link{mclustOptions}}
}
\examples{
faithfulModel <- mclustModel(faithful,mclustBIC(faithful))

mclust2Dplot(faithful, parameters=faithfulModel$parameters, 
             z=faithfulModel$z, what = "classification", identify = TRUE)

mclust2Dplot(faithful, parameters=faithfulModel$parameters, 
             z=faithfulModel$z, what = "uncertainty", identify = TRUE)
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
