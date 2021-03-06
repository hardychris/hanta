% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pca.R
\name{gm.pca}
\alias{gm.pca}
\title{Gene Matrix Principal Component Analysis}
\usage{
gm.pca(gm, metadata, grouping_var = NULL, filt_method = FALSE,
  thresh_cut = 0, top_genes = ncol(gm), transform_type = "none",
  plot = FALSE, pcs = 4, file, ext = "pdf")
}
\arguments{
\item{gm}{Gene Matrix.}

\item{metadata}{Metadata table with columns of named features and rows of unique barcodes.}

\item{grouping_var}{Variables of interest from metadata file.}

\item{filt_method}{Filtering method. Choose from from 'quant_exp' (Filter genes below quantile value), 'top_var' (Retain genes with highest variance), or 'top_exp (Retain highest expressing genes).}

\item{thresh_cut}{Filtering threshold if using 'quant_exp' method. Defaults to 0 (all genes retained).}

\item{top_genes}{Number of genes to retain if using 'top_var' or 'top_exp' methods. Defaults to the total number of genes (all genes retained).}

\item{transform_type}{Transformation type from ("none", "ln" or any integer).  Defaults to "none" for no transformation.  Possible options are "ln" for natural log, or an integer indicating the base of the lograthim used in the transformation (ie. For log2 transformed data, enter 2)}

\item{plot}{A boolean defining whether or not to create a PCA plot. Defaults to FALSE.}

\item{pcs}{The number of Principal Components to plot (up to 4).}

\item{file}{A file destination/name for the plot output.}
}
\description{
This function is used to filter a gene matrix, calculate the prinicpal components, and create a PCA plot.  Returns a PCA object from the prcomp function in the 'stats' package.
}
\examples{
# Perform Principal Component Analysis on 500 highest expressing genes, plot
# the results of the first 4 Principal Components.
}
\keyword{Gene}
\keyword{Matrix}
\keyword{PCA}
