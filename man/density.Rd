\name{density}
\alias{density}
\alias{print.density}
\title{Kernel Density Estimation}
\usage{
density(\dots, method, G)
}
\arguments{
  \item{\dots}{Arguments to the \code{density} function in the base
    package.}

  \item{method}{If equal to "mclust", \code{EMclust} is used to estimate
    the density.}

  \item{G}{The number of gaussians to consider in the model-based
    density estimation. Default: 1:9. Ignored if method is not equal to
    "mclust".}
}
\description{
  This is exactly the same function as in the \code{stats} package but for the
  \code{method} argument: if it is given and equals \code{"mclust"}, the
  \code{mclust} density estimation is used. Optionally, the number of
  gaussians to be considered can be given as well (\code{G}).}
}
\value{
  If \code{give.Rkern} is true, the number \eqn{R(K)}, otherwise
  an object with class \code{"density"} whose
  underlying structure is a list containing the following components.
  \item{x}{the \code{n} coordinates of the points where the density is
    estimated.}
  \item{y}{the estimated density values.}
  \item{bw}{the bandwidth used.}
  \item{N}{the sample size after elimination of missing values.}
  \item{call}{the call which produced the result.}
  \item{data.name}{the deparsed name of the \code{x} argument.}
  \item{has.na}{logical, for compatibility (always FALSE).}
}
\references{
  Fraley, C. and Raftery, A.E. (2002)
  MCLUST: software for model-based clustering, density estimation and
  discriminant analysis.
  Technical Report No. 415,
  Dept. of Statistics, University of Washington.
  
  Scott, D. W. (1992)
  \emph{Multivariate Density Estimation. Theory, Practice and Visualization}.
  New York: Wiley.

  Sheather, S. J. and Jones M. C. (1991)
  A reliable data-based bandwidth selection method for kernel density
  estimation.
  \emph{J. Roy. Statist. Soc.} \bold{B}, 683--690.

  Silverman, B. W. (1986)
  \emph{Density Estimation}.
  London: Chapman and Hall.

  Venables, W. N. and Ripley, B. D. (1999)
  \emph{Modern Applied Statistics with S-PLUS}.
  New York: Springer.
}
\seealso{
  \code{\link[base]{density}} (base package), \code{\link{bw.nrd}},
  \code{\link{plot.density}}, \code{\link{hist}}.
}
\examples{
plot(density(c(-20,rep(0,98),20)), xlim = c(-4,4))# IQR = 0

# The Old Faithful geyser data
data(faithful)
d <- density(faithful$eruptions, bw = "sj")
d
plot(d)
dmc <- density(faithful$eruptions, method="mclust")
plot(dmc, type = "n")
polygon(dmc, col = "wheat")
lines(d, col="red")

## Missing values:
x <- xx <- faithful$eruptions
x[i.out <- sample(length(x), 10)] <- NA
doRmc <- density(x=x, method="mclust", na.rm = TRUE)
lines(doRmc, col="blue")
doR <- density(x, bw = 0.15, na.rm = TRUE)
lines(doR, col = "green")
rug(x)
points(xx[i.out], rep(0.01, 10))

(kernels <- eval(formals(stats:::density.default)$kernel))

## show the kernels in the R parametrization
plot (density(0, bw = 1), xlab = "",
      main="R's density() kernels with bw = 1")
for(i in 2:length(kernels))
   lines(density(0, bw = 1, kern =  kernels[i]), col = i)
legend(1.5,.4, legend = kernels, col = seq(kernels),
       lty = 1, cex = .8, y.int = 1)

data(precip)
bw <- bw.SJ(precip) ## sensible automatic choice
plot(density(precip, bw = bw, n = 2^13))
lines(density(precip, G=2:5, method="mclust"), col="red")
rug(precip)
}
\keyword{cluster}
\keyword{distribution}
\keyword{smooth}
