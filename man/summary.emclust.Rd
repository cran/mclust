\name{summary.emclust}
\alias{print.summary.emclust}
\title{
Summary method for `emclust' objects.
}
\description{
Optimal model parameters for \code{emclust} objects.
}
\usage{
summary.emclust(x, data, nclus, modelid)
print.summary.emclust(summary.emclust.object, \dots)
}
\arguments{
\item{x}{
an \code{emclust} object.
}
\item{data}{
the data used to obtain the \code{emclust} object.
}
\item{nclus}{
a vector of integers giving the number of clusters over which the
summary is to take place (must be a subset of the column names of \code{x}).
The default is to summarize over all numbers of clusters used in the original 
analysis.
}
\item{modelid}{
a vector of character strings denoting the models over which the
summary is to take place (must be a subset of the row names of \code{x}).
The default is to summarize over all models used in the original
analysis.
}}
\value{
A list giving the optimal (according to BIC) parameter and z values, 
together with the associated classification and its uncertainty.
}
\seealso{
\code{\link{emclust}}
}
\examples{
data(iris)
iris.bic <- emclust(iris[,1:4], nclus=3:5)
summary(iris.bic, iris[,1:4])


}
\keyword{clustering}
% Converted by Sd2Rd version 0.3-2.


