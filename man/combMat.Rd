\name{combMat}
\alias{combMat}
\title{
Combining Matrix
}
\description{
Create a combining matrix
}
\usage{
combMat(K, l1, l2)
}
\arguments{
  \item{K}{
	The original number of classes: the matrix will define a combining from K to (K-1) classes.
	}
  \item{l1}{
	Label of one of the two classes to be combined.
	}
  \item{l2}{
	Label of the other class to be combined.
	}
}
\value{
If \code{z} is a vector (length \emph{K}) whose \emph{k}th entry is the probability that an observation belongs to the \emph{k}th class in a \emph{K}-classes classification, then \code{combiM \%*\% z} is the vector (length \emph{K-1}) whose \emph{k}th entry is the probability that the observation belongs to the \emph{k}th class in the \emph{K-1}-classes classification obtained by merging classes \code{l1} and \code{l2} in the initial classification.
}
\author{
J.-P. Baudry, A. E. Raftery, L. Scrucca
}
\seealso{
\code{\link{clustCombi}}, \code{\link{combiPlot}}
}
%\examples{}
\keyword{ cluster }

