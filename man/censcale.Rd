\name{censcale}
\alias{censcale}
\title{
Centering and Scaling of Data
}
\section{DESCRIPTION}{
Transforms a data set by subtracting the column mean from each column and
dividing each column by the square root of its variance when the latter is
sufficiently large in magnitude.
}
\usage{
censcale(x, tol = 0.0001)


}
\arguments{
\item{x}{
Matrix of observations.
}
\item{tol}{
Threshold for column scaling. A column is divided by the square root of its
variance (standard deviation) when that quantity exceeds \code{tol}, The default
value is \code{0.0001}.
}}
\value{
The centered and scaled matrix of observations.
}
\examples{
data(iris)
iris.sc <- censcale(iris[,1:4])


}
\keyword{}
% Converted by Sd2Rd version 0.3-2.


