\name{mclust-package}
\alias{mclust-package}
\alias{mclust}
\docType{package}

\title{Gaussian Mixture Modelling for Model-Based Clustering, Classification, and Density Estimation}

\description{
Gaussian finite mixture models estimated via EM algorithm for model-based clustering, classification, and density estimation, including Bayesian regularization and dimension reduction.
}

\details{
For a quick introduction to \pkg{mclust} see the vignette \href{../doc/mclust.html}{A quick tour of mclust}.

See also:
\itemize{
  \item \code{\link{Mclust}} for clustering;
  \item \code{\link{MclustDA}} for supervised classification;
  \item \code{\link{MclustSSC}} for semi-supervised classification;
  \item  \code{\link{densityMclust}} for density estimation.
  }
}

\author{
Chris Fraley, Adrian Raftery and Luca Scrucca.

Maintainer: Luca Scrucca \email{luca.scrucca@unipg.it}
}

\references{
Scrucca L., Fraley C., Murphy T. B. and Raftery A. E. (2023) \emph{Model-Based Clustering, Classification, and Density Estimation Using mclust in R}. Chapman & Hall/CRC, ISBN: 978-1032234953, https://mclust-org.github.io/book/

Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 289-317. 

Fraley C. and Raftery A. E. (2002) Model-based clustering, discriminant analysis and density estimation, \emph{Journal of the American Statistical Association}, 97/458, pp. 611-631.
}

\examples{
\donttest{
# Clustering
mod1 <- Mclust(iris[,1:4])
summary(mod1)
plot(mod1,  what = c("BIC", "classification"))

# Classification
data(banknote)
mod2 <- MclustDA(banknote[,2:7], banknote$Status)
summary(mod2)
plot(mod2)

# Density estimation
mod3 <- densityMclust(faithful$waiting)
summary(mod3)
}
}

\keyword{package}
