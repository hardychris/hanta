% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{pca.elbow}
\alias{pca.elbow}
\title{PCA Elbow}
\usage{
pca.elbow(pca_obj, elbow_thresh = 0.05)
}
\arguments{
\item{pca_obj}{A Principal Component Analysis object from the prcomp function in the 'stats' package.}

\item{elbow_thresh}{The absolute fractional difference to define the stopping point along the elbow. Defaults to 0.05 (ie. a < 5\% change from 1 PC to the next defines the elbow).}
}
\description{
This function is used to estimate the number of Principal Components used as input for tSNE.
}
\examples{
# Find the number of principal components at the elbow of a scree plot.
}
\keyword{Elbow}
\keyword{PCA}
