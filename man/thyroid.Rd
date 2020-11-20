\name{thyroid}
\alias{thyroid}
\docType{data}

\title{Thyroid gland data}

\description{
Data on five laboratory tests administered to a sample of 215 patients. The tests are used to predict whether a patient's thyroid can be classified as euthyroidism (normal thyroid gland function), hypothyroidism (underactive thyroid not producing enough thyroid hormone) or hyperthyroidism (overactive thyroid producing and secreting excessive amounts of the free thyroid hormones T3 and/or thyroxine T4). Diagnosis of thyroid operation was based on a complete medical record, including anamnesis, scan, etc.}

\usage{data(thyroid)}

\format{A data frame with the following variables:

\describe{
\item{Diagnosis}{Diagnosis of thyroid operation: \code{Hypo}, \code{Normal}, and \code{Hyper}.}
\item{RT3U}{T3-resin uptake test (percentage).}
\item{T4}{Total Serum thyroxin as measured by the isotopic displacement method.}
\item{T3}{Total serum triiodothyronine as measured by radioimmuno assay.}
\item{TSH}{Basal thyroid-stimulating hormone (TSH) as measured by radioimmuno assay.}
\item{DTSH}{Maximal absolute difference of TSH value after injection of 200 micro grams of thyrotropin-releasing hormone as compared to the basal value.}
}

}

\source{UCI \url{ftp://ftp.ics.uci.edu/pub/machine-learning-databases/thyroid-disease/}}

\references{
Coomans, D., Broeckaert, M. Jonckheer M. and Massart D.L. (1983)
Comparison of Multivariate Discriminant Techniques for Clinical Data - Application to the Thyroid Functional State, \emph{Meth. Inform. Med.} 22, pp. 93-101.

Coomans, D. and I. Broeckaert (1986) \emph{Potential Pattern Recognition in Cemical and Medical Decision Making}, Research Studies Press, Letchworth, England.
} 

\keyword{datasets}
