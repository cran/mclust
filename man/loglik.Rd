\name{loglik}
\alias{loglik}
\title{
Loglikelihood for model-based hierarchical clustering.
}
\section{DESCRIPTION}{
Gives the loglikelihood for each stage of model-based hierarchical clustering.
}
\usage{
loglik(tree, data, Vinv)


}
\arguments{
\item{tree}{
an \code{"mhtree"} object.
}
\item{data}{
the data used to produce the \code{"mhtree"} object.
}
\item{Vinv}{
approximate recoprocal hypervolume of the region from which the data is drawn.
The default (for those methods that need this quantity) is determined by the
function \code{hypvol}.
}}
\value{
the loglikelihood corresponding to the initial partition and to each stage of 
merging in hierarchical clustering, together with the following attribute:

\item{nmerge}{
the number of clusters merged at each stage.
}}
\section{NOTES}{
The value given is equal to the loglikelihood up to an additive constant. 
For those models in which they arise, indeterminate terms are 
assigned the value \code{-k*log(volume)}, where \code{k} is the number of observations 
associated with the term. 
If you scaled your data before using \code{mhtree}, be 
sure to use the same scaling when supplying the data to \code{loglik}.
}
\seealso{
\code{\link{mhtree}}, \code{\link{awe}}
}
\examples{
data(iris)
loglik(mhtree(iris[,1:4]),iris[,1:4])
}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.




