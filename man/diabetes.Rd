\name{diabetes}
\alias{diabetes}
\docType{data}

\title{Diabetes Data (flawed)}

\description{The data set contains three measurements made on 145 non-obese adult patients classified into three groups.}

\usage{data(diabetes)}

\format{A data frame with the following variables:
\describe{
\item{class}{The type of diabete: \code{Normal}, \code{Overt}, and \code{Chemical}.}
\item{glucose}{Area under plasma glucose curve after a three hour oral glucose tolerance test (OGTT).}
\item{insulin}{Area under plasma insulin curve after a three hour oral glucose tolerance test (OGTT).}
\item{sspg}{Steady state plasma glucose.}
}
}

\details{This dataset is \emph{flawed} (compare with the reference) and it is provided here only for backward compatibility. A 5-variable version of the Reaven and Miller data is available in package \pkg{rrcov}. The \emph{glucose} and \emph{sspg} columns in this dataset are identical to the \emph{fpg} and \emph{insulin} columns, respectively in the \pkg{rrcov} version. The \emph{insulin} column in this dataset differs from the \emph{glucose} column in the \pkg{rrcov} version in one entry: observation 104 has the value 45 in the \emph{insulin} column in this data, and 455 in the corresponding \emph{glucose} column of the \pkg{rrcov} version.}

\source{Reaven, G. M. and Miller, R. G. (1979). An attempt to define the nature of chemical diabetes using a multidimensional analysis. \emph{Diabetologia} 16:17-24.}

\keyword{datasets}
