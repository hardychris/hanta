% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor.R
\name{cor.heatmap}
\alias{cor.heatmap}
\title{Correlation Matrix Heatmap}
\usage{
cor.heatmap(cm, gm = NULL, metadata, cor_type = "", grouping_var = NULL,
  gene_highlight = c(), gene_highlight_multi = FALSE, cluster = FALSE,
  label = FALSE, rel_scale = FALSE, palette = "RdGy",
  palette_rev = FALSE, file, ext = "pdf", plot_height = 6.69,
  plot_width = 10, plot_res = 600, legend_title = "group", ...)
}
\arguments{
\item{cm}{Correlation matrix.}

\item{metadata}{Metadata table with columns of named features and rows of unique barcodes.}

\item{cor_type}{Correlation type. Choose from 'pearson', 'spearman' or 'kendall'. Defaults to 'pearson'.}

\item{grouping_var}{Grouping variable(s) used to define plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))}

\item{cluster}{TRUE/FALSE Perform clustering.}

\item{label}{TRUE/FALSE Label rows and columns with sample names.}

\item{rel_scale}{TRUE/FALSE relative scale. If FALSE, color scale will range from 0:1 or -1:1 depending on the minimum value of the correlation matrix.  If TRUE color scale will range from the mimimum of the correlation matrix to 1.}

\item{palette}{A color palette selected from the brewer.pal function in the RColorBrewer package in R.  Defaults to "RdBu".  See 'http://www.datavis.ca/sasmac/brewerpal.html' for available options.}

\item{file}{A file destination/name for the plot output.  File format must be '.pdf'.  If not provided, plot will be rendered in console.}

\item{...}{Optional named parameters to be sent to the aheatmap function in the NMF package.}
}
\description{
Generates a heatmap from a correlation matrix.
}
\examples{

}
\keyword{Correlation}
\keyword{Gene}
\keyword{Matrix}
