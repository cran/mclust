\name{awe}
\alias{awe}
\title{
Approximate weight of evidence for model-based hierarchical clustering.
}
\description{
Computes a Bayesian criterion for assessing the number of clusters present
in the data.
}
\usage{
awe(tree, data)


}
\arguments{
\item{tree}{
an \code{"mhtree"} object.
}
\item{data}{
the data used to produce the \code{"mhtree"} object.
}}
\value{
the approximate weight of evidence for each possible stage of merging.
}
\section{NOTES}{
Since \code{"mhtree"} allows stopping and
starting at any stage, the result will contain NAs for those stages that have 
been eliminated. 
If you scaled your data before using \code{mhtree}, be sure
to use the same scaling when supplying the data to \code{awe}.
}
\references{
J. D. Banfield and A. E. Raftery, Model-based Gaussian and non-Gaussian
Clustering, \emph{Biometrics, }49:803-821 (1993).
}
\seealso{
\code{\link{mhtree}}, \code{\link{loglik}}
}
\examples{
data(iris)
iris.m _ iris[,1:4]
awe.val <- awe(mhtree(iris.m), iris.m)
plot(awe.val)


}
\keyword{cluster}
% Converted by Sd2Rd version 0.3-2.




