\name{ztoc}
\alias{ctoz}
\title{Conversion between conditional probabilities and a classification}
\description{The function ctoz converts a discrete classification into
  the corresponding z matrix. ztoc converts a z matrix to the nearest
  classification. 
}
\usage{
ztoc(z)
ctoz(c, noise)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{z}{Matrix with conditional probabilities; each object one row,
    and each cluster one column.}
  \item{c}{Classification vector.}
  \item{noise}{Noise indicator vector.}
}

\examples{
data(iris)
cltree _ mhtree(iris[,1:4], modelid="VVV")
cl _ mhclass(cltree, 3)
origz _ ctoz(cl)
z _ me(iris[,1:4], modelid="VVV", origz)
newcl _ ztoc(z)
}
\keyword{clustering}%-- one or more ...
