\name{mclustOptions}
\alias{mclustOptions}
\title{
  Set default values for use with MCLUST.
}
\description{
  Supplies a list of values an enumeration of models for use
  with \emph{MCLUST}.
}
\usage{
mclustOptions(emModelNames=NULL, hcModelNames=NULL, 
              bicPlotSymbols=NULL, bicPlotColors=NULL, 
              classPlotSymbols=NULL, classPlotColors=NULL, warn=TRUE) 
}
\arguments{
  \item{emModelNames}{
    A vector of 3-character strings that are associated with multivariate 
    models for which EM estimation is available in MCLUST. \cr
    The current default is  all of the multivariate mixture models
    supported in \code{MCLUST}.
    The help file for \code{mclustModelNames} describes the available models.
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
    A vector whose entries are either integers corresponding to graphics
    symbols or single characters for plotting BIC curves.
    The default is \cr
   \code{c(EII=17,VII=2,EEI=16,EVI=10,VEI=13,VVI=1, \cr
           EEE=15,EEV=12,VEV=7,VVV=0,E=17,V=2)}.
  }
  \item{bicPlotColors}{
    A vector whose entries are either integers corresponding to colors
    to BIC curves.
   \code{c(EII="gray",VII="black",\cr
           EEI="orange",EVI="brown",VEI="red",VVI="magenta", \cr
           EEE="forestgreen",EEV="green",VEV="cyan",VVV="blue",\cr
           E="gray",V="black")}.
  }
  \item{classPlotSymbols}{
    A vector whose entries are either integers corresponding to graphics
    symbols or single characters for plotting for
    classifications. Classes are assigned symbols in the given order.
    The default is \code{c(17,0,10,4,11,18,6,7,3,16,2,12,8,15,1,9,14,13,5)}. 
  }
  \item{classPlotColors}{
    A vector whose entries are either integers corresponding to graphics
    symbols or single characters for plotting for
    classifications. Classes are assigned symbols in the given order.
    The default is\cr
     \code{ "blue", "red", "green", "cyan", "magenta", \cr  
     "forestgreen", "purple", "orange", "gray", "brown", "black")}
  }
  \item{warn}{
    A logical value allowing some types of warnings to be turned on
    or off globally. Most of these warnings have to do with  
    situations in which singularities are encountered. 
    The default is \code{warn = TRUE}. 
  }
}
\value{
  A named list in which the names are the names of the arguments
  and the values are the values supplied to the arguments.
}
\details{
  \code{mclustOptions} is provided for assigning values
  to the \code{.Mclust} list, which is used
  to supply default values to various functions in MCLUST.
  
  Calls to \code{mclustOptions} do not in themselves affect
  the outcome of computations.
}
\section{References}{
  C. Fraley and A. E. Raftery (2002).
  Model-based clustering, discriminant analysis, and density estimation.
  \emph{Journal of the American Statistical Association 97:611-631}. 

  C. Fraley and A. E. Raftery (2006).
  MCLUST Version 3 for R: Normal Mixture Modeling and Model-Based Clustering, 
  Technical Report no. 504, Department of Statistics,
  University of Washington.
}
\seealso{
  \code{\link{.Mclust}},
  \code{\link{emControl}}
}
\examples{
irisBIC <- mclustBIC(iris[,-5])
summary(irisBIC, iris[,-5])

.Mclust
.Mclust <- mclustOptions(emModelNames = c("VII", "VVI", "VVV"))
.Mclust

irisBIC <- mclustBIC(iris[,-5])
summary(irisBIC, iris[,-5])

.Mclust <- mclustOptions() # restore default values
.Mclust
}
\keyword{cluster}
% docclass is function
% Converted by Sd2Rd version 1.21.
