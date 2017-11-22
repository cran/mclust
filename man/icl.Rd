\name{icl}
\alias{icl}
\title{
  ICL for an estimated Gaussian Mixture Model
}
\description{
Computes the ICL (Integrated Complete-data Likelihood) for criterion for a Gaussian Mixture Model fitted by \code{\link{Mclust}}.
}

\usage{
icl(object, \dots)
}

\arguments{
  \item{object}{
    An object of class \code{'Mclust'} resulting from a call to \code{\link{Mclust}}.
  }
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
  The ICL for the given input MCLUST model.
}

\references{
  Biernacki, C., Celeux, G., Govaert, G. (2000). 
  Assessing a mixture model for clustering with the integrated completed likelihood.
  \emph{IEEE Trans. Pattern Analysis and Machine Intelligence}, 22 (7), 719-725.
}
\seealso{
  \code{\link{Mclust}},
  \code{\link{mclustBIC}},
  \code{\link{mclustICL}},
  \code{\link{bic}}.
}

\examples{
mod <- Mclust(iris[,1:4])
icl(mod)
}

\keyword{cluster}
