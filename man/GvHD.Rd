\name{GvHD}
\alias{GvHD}
\alias{GvHD.pos}
\alias{GvHD.control}

\docType{data}
\title{GvHD Dataset}
\description{
GvHD (Graft-versus-Host Disease) data of Brinkman et al. (2007). Two samples of this flow cytometry data, one from a patient with the GvHD, and the other from a control patient. The GvHD positive and control samples consist of 9083 and 6809 observations, respectively. Both samples include four biomarker variables, namely, CD4, CD8b, CD3, and CD8. The objective of the analysis is to identify CD3+ CD4+ CD8b+ cell sub-populations present in the GvHD positive sample.

A treatment of this data by combining mixtures is proposed in Baudry et al. (2010).
}
\usage{data(GvHD)}	
\format{
  GvHD.pos (positive patient) is a data frame with 9083 observations on the following 4 variables, which are biomarker measurements.
  \describe{
    \item{CD4}{}
    \item{CD8b}{}
    \item{CD3}{}
    \item{CD8}{}
  }  
  GvHD.control (control patient) is a data frame with 6809 observations on the following 4 variables, which are biomarker measurements.
  \describe{
    \item{CD4}{}
    \item{CD8b}{}
    \item{CD3}{}
    \item{CD8}{}
  }
}
\references{
R. R. Brinkman, M. Gasparetto, S.-J. J. Lee, A. J. Ribickas, J. Perkins, W. Janssen, R. Smiley and C. Smith (2007). High-content flow cytometry and temporal data analysis for defining a cellular signature of Graft-versus-Host Disease. \emph{Biology of Blood and Marrow Transplantation, 13: 691-700.}

K. Lo, R. R. Brinkman, R. Gottardo (2008). Automated gating of flow cytometry data via robust model-based clustering. \emph{Cytometry A, 73: 321-332.}

J.-P. Baudry, A. E. Raftery, G. Celeux, K. Lo and R. Gottardo (2010). Combining mixture components for clustering. \emph{Journal of Computational and Graphical Statistics, 19(2):332-353.}
}
\examples{
\dontrun{
data(GvHD)

dat <- GvHD.pos[1:500,] # only a few lines for a quick example

output <- clustCombi(dat) 

output # is of class clustCombi

# plots the hierarchy of combined solutions, then some "entropy plots" which 
# may help one to select the number of classes
plot(output, dat) 
}

}
\keyword{datasets}
