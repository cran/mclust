\name{mclust1Dplot}
\alias{mclust1Dplot}
\title{
  Plot one-dimensional data modelled by an MVN mixture.
}
\description{
  Plot one-dimensional data given parameters of an MVN mixture model 
  for the data.
}
\usage{
mclust1Dplot(data, \dots,
             type = c("classification","uncertainty","density","errors"),
             ask = TRUE, symbols, grid = 100, identify = FALSE, CEX = 1, xlim) 
}
\arguments{
  \item{data}{
    A numeric vector of observations.
    Categorical variables are not allowed.
  }
  \item{\dots}{
    One or more of the following:
    \describe{
      \item{classification}{
	A numeric or character vector representing a classification of
	observations (rows) of \code{data}.
      }
      \item{uncertainty}{
	A numeric vector of values in \emph{(0,1)} giving the
	uncertainty of each data point.
      }
      \item{z}{
	A matrix in which the \emph{[i,k]}the entry gives the
	probability of observation \emph{i} belonging to the \emph{k}th class. 
	Used to compute \code{classification} and
	\code{uncertainty} if those arguments aren't available.
      }
      \item{truth}{
	A numeric or character vector giving a known
	classification of each data point.
	If \code{classification}
	or \code{z} is also present, this  
	is used for displaying classification errors.
      }
      \item{mu}{
	A vector whose entries are the means of each group. 
      }
      \item{sigma}{
	Either a vector whose entries are the variances for each group or a
	scalar giving a common variance for the groups.
      }
      \item{pro}{
	The vector of mixing proportions.
      }
    }
  }
  \item{type}{
    Any subset of
    \code{c("classification","uncertainty","density","errors")}.
    The function will produce the corresponding plot if
    it has been supplied sufficient information to do so. If more than
    one plot is possible then users will be asked to choose from a menu
    if \code{ask=TRUE}. 
  }
  \item{ask}{
    A logical variable indicating whether or not a menu should be produced
    when more than one plot is possible. The default is \code{ask=TRUE}.
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
  \item{grid}{
    Number of grid points to use.
  }
  \item{identify}{
    A logical variable indicating whether or not to add a title to the plot
    identifying the dimensions used.
  }
  \item{CEX}{
    An argument specifying the size of the plotting symbols. 
    The default value is 1.
  }
  \item{xlim}{
    An argument specifying bounds of the plot.
    This may be useful for when comparing plots.
  }
}
\section{Side Effects}{
  One or more plots showing location of the mixture components, classification,
  uncertainty, density and/or classification errors. Points in the different
  classes are shown in separate lines above the whole of the data.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 
  See \url{http://www.stat.washington.edu/mclust}.
  
  
  C. Fraley and A. E. Raftery (2002).
  MCLUST:Software for model-based clustering, density estimation and
  discriminant analysis.
  Technical Report, Department of Statistics, University of Washington. 
  See \url{http://www.stat.washington.edu/mclust}.
}
\seealso{
  \code{\link{mclust2Dplot}},
  \code{\link{clPairs}},
  \code{\link{coordProj}},
  \code{\link{do.call}}
}
\examples{
n <- 250 ## create artificial data
set.seed(0)
y <- c(rnorm(n,-5), rnorm(n,0), rnorm(n,5))
yclass <- c(rep(1,n), rep(2,n), rep(3,n))

yEMclust <- summary(EMclust(y),y)

mclust1Dplot(y, identify = TRUE, truth = yclass, z = yEMclust$z, ask=FALSE,
             mu = yEMclust$mu, sigma = yEMclust$sigma, pro = yEMclust$pro)

do.call("mclust1Dplot",
        c(list(data = y, identify = TRUE, truth = yclass, ask=FALSE), 
        yEMclust))
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.

