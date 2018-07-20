#' Gene Matrix - Find Highly Variable Genes
#'
#' This function is used to find highly variable genes.  Genes are seperated into a user-defined number of bins based on their mean expression.  Variance is then scaled for each bin, and the genes with the highest scaled variance are returned.
#' @param gm Gene Matrix.
#' @param n_bin Number of bins from which to calculate scaled variance.  Defaults to 20.
#' @param transform_type Transformation type from ("none", "ln" or any integer).  Defaults to "none" for no transformation.  Possible options are "ln" for natural log, or an integer indicating the base of the lograthim used in the transformation (ie. For log2 transformed data, enter 2)
#' @param lo_mean_cut Only consider genes highly variable with a mean greater than this number.  Defaults to 0.10.
#' @param hi_mean_cut Only consider genes highly variable with a mean less than this number.  Defaults to Inf for Infinity.
#' @param lo_dsprs_scaled_cut Only consider genes highly variable if their scaled dispersion is greater than this number.  Defaults to 1.
#' @param hi_dsprs_scaled_cut Only consider genes highly variable if their scaled dispersion is less than this number.  Defaults to Inf for Infinity.
#' @param legend TRUE/FALSE, if generating the MA plot, include a legend?  Defaults to FALSE.
#' @param file A file destination/name for the plot output.  File format must be '.pdf'.  If not provided, plot will be rendered in console.
#' @keywords Gene Filtering Highly Variable
#' @export
#' @examples
#' # Find highly variable genes, that have mean expression greater than 0.25 and

gm.var_genes <- function(gm, n_bin, transform_type = "none", lo_mean_cut = 0.10, hi_mean_cut = Inf,
                        lo_dsprs_scaled_cut = 1, hi_dsprs_scaled_cut = Inf, legend = FALSE, file, ext = "pdf"){

  ext = hanta:::setup.figure_type(ext)

  df = pca.var_genes(gm, n_bin, transform_type)

  df_pass = pca.var_genes_sig(df, lo_mean_cut, hi_mean_cut, lo_dsprs_scaled_cut, hi_dsprs_scaled_cut)

  df_fail = df[-match(rownames(df_pass), rownames(df)), ]

  if(missing(file) == F){

    if(ext == "png"){

      png(pkgmaker::file_extension(file, ".png"), height = 5.69, width = 10.00, units = 'in', res = 600)

    } else {

      pdf(pkgmaker::file_extension(file, ".pdf"), height = 5.69, width = 10.00, onefile = FALSE)

    }

  }

  plot(x = df_fail$means, y = df_fail$dsprs_scaled,
       col = rgb(0, 0, 0, 0.5),
       main = "MA Plot",
       font.main = 3,
       xlab = "Mean",
       ylab = "Scaled Dispersion",
       xlim = c(min(c(df_fail$means, df_pass$means)), max(c(df_fail$means, df_pass$means))),
       ylim = c(min(c(df_fail$dsprs_scaled, df_pass$dsprs_scaled)), max(c(df_fail$dsprs_scaled, df_pass$dsprs_scaled))),
       pch = 16)

  points(x = df_pass$means, y = df_pass$dsprs_scaled, col = rgb(1, 0, 0, 0.5), pch = 16)

  if(legend == T){
    legend('topright',
           legend = "pass filter",
           col = 'red',
           pch = 16,
           bty = 'n')
  }

  if(missing(file) == F){

    graphics.off()

  }

  out = list(total = df, pass = df_pass)

  return(out)

}


# Calculate means and scaled dispersions for genes in a gene matrix.  Derived from Seurat.
pca.var_genes <- function(gm, n_bin, transform_type){

  if(missing(transform_type) == TRUE){transform_type = "ln"}

  means <- apply(gm, 1, function(x){

    if(transform_type == "ln"){

      log1p(mean(expm1(x)))

    } else if(transform_type == "none"){

      mean(x)

    } else if(is.numeric(as.numeric(transform_type)) == TRUE){

      log(mean(as.numeric(transform_type) ^ (x) - 1) + 1, base = as.numeric(transform_type))

    }

  })

  dsprs <- apply(gm, 1, function(x){

    if(transform_type == "ln"){

      mu = mean(expm1(x))

      var = sum((expm1(x) - mu) ^ 2)

      var = (var + (length(x) - length(x[x!=0])) * (mu ^ 2)) / (length(x) - 1)

      log(var / mu)

    } else if(transform_type == "none"){

      mu = mean(x)

      var = sum((x - mu) ^ 2)

      var = (var + (length(x) - length(x[x!=0])) * (mu ^ 2)) / (length(x) - 1)

      var / mu

    } else if(is.numeric(as.numeric(transform_type)) == TRUE){

      mu = mean((as.numeric(transform_type) ^ x) - 1)

      var = sum((((as.numeric(transform_type) ^ x) - 1) - mu) ^ 2)

      var = (var + (length(x) - length(x[x!=0])) * (mu ^ 2)) / (length(x) - 1)

      log(var / mu, base = transform_type)

    }

  })

  dsprs[is.na(dsprs)] <- 0
  means[is.na(means)] <- 0

  means_bins = cut(means, breaks = n_bin)
  mean_dsprs = tapply(dsprs, means_bins, mean)
  sd_dsprs = tapply(dsprs, means_bins, sd)

  dsprs_scaled = (dsprs - mean_dsprs[as.numeric(means_bins)])/sd_dsprs[as.numeric(means_bins)]
  dsprs_scaled[is.na(dsprs_scaled)] <- 0
  out_df = data.frame(means, dsprs, dsprs_scaled)
  rownames(out_df) = rownames(gm)

  out_df = out_df[order(-dsprs_scaled), ]

  return(out_df)

}


# Filter genes that are outside of user-defined thresholds for mean expression and scaled dispersion.
pca.var_genes_sig <- function(df, lo_mean_cut = 0.10, hi_mean_cut = Inf,
                          lo_dsprs_scaled_cut = 1, hi_dsprs_scaled_cut = Inf){

  pass = subset(df, df$means > lo_mean_cut & df$means < hi_mean_cut & df$dsprs_scaled > lo_dsprs_scaled_cut & df$dsprs_scaled < hi_dsprs_scaled_cut)

  return(pass)

}


#' PCA Gene Filter
#'
#' This function is used to filter genes from a gene matrix as input for Principal Component Analysis
#' @param gm Gene Matrix.
#' @param method Filtering method. Choose from from 'quant_exp' (Filter genes below quantile value), 'top_var' (Retain genes with highest scaled variance), or 'top_exp (Retain highest expressing genes).
#' @param thresh_cut Quantile filtering threshold if using 'quant_exp' method.  Value should be between 0.0 and 1.0.  Defaults to 0 (all genes retained).
#' @param top_genes Number of genes to retain if using 'top_var' or 'top_exp' methods.  Defaults to the number of columns in the gene matrix (all genes retained).
#' @keywords PCA Filter
#' @export
#' @examples
#'

pca.gene_filter <- function(gm, method, thresh_cut = 0.0, top_genes = ncol(gm), transform_type = "none"){

  if(missing(method) == TRUE | method == FALSE){
    return(gm)
    } else if(!(method %in% c("quant_exp", "top_var", "top_exp")) == TRUE){
        stop("Choose a gene filtering method from 'quant_exp' (Filter genes below quantile value), or 'top_var' (Pick genes with highest variance), or 'top_exp' (Pick highest expressing genes).")
  }

  if(method == 'quant_exp'){
    gene_cov = rowSums(gm)
    pca_genes = gene_cov > quantile(gene_cov, probs = c(thresh_cut))
    pca_genes = gm[pca_genes, ]
    return(pca_genes)
  }

  if(method == 'top_var'){
    hi_var_genes = head(pca.var_genes_sig(pca.var_genes(gm, 20, transform_type = transform_type)), top_genes)
    pca_genes = gm[rownames(hi_var_genes), ]
    return(pca_genes)
  }

  if(method == 'top_exp'){
    gene_cov = rowSums(gm)
    pca_genes = names(head(sort(gene_cov, decreasing = TRUE), top_genes))
    pca_genes = gm[pca_genes, ]
    return(pca_genes)
  }

}


#' PCA Elbow
#'
#' This function is used to estimate the number of Principal Components used as input for tSNE.
#' @param pca_obj A Principal Component Analysis object from the prcomp function in the 'stats' package.
#' @param elbow_thresh The absolute fractional difference to define the stopping point along the elbow. Defaults to 0.05 (ie. a < 5\% change from 1 PC to the next defines the elbow).
#' @keywords PCA Elbow
#' @export
#' @examples
#' # Find the number of principal components at the elbow of a scree plot.

pca.elbow <- function(pca_obj, elbow_thresh = 0.05){

  lag = 1

  x = summary(pca_obj)$importance[2,]

  pc_change = abs(c(diff(x, lag),rep(0, lag))/ x)

  pcs = min(which(pc_change < elbow_thresh))

  return(pcs)

}


#' PCA Count
#'
#' This function is used to estimate the number of Principal Components that explain a user-defined percentage of the variance in the data.
#' @param pca_obj A Principal Component Analysis object from the prcomp function in the 'stats' package.
#' @param perc_var The number of Principal Components required to explain 'per_var' (\% Variance) of the data.  Value should be between 0.0 and 1.0, defaults to 0.95 for 95\% of variance.
#' @keywords PCA Count
#' @export
#' @examples
#' # Find the number of principal components that explain 95\% of the data


pca.count <- function(pca_obj, perc_var = 0.95){

  frac_exp = summary(pca_obj)$importance[3,]

  pcs = min(which(frac_exp > perc_var))

  return(pcs)

}


#' PCA Plot
#'
#' This function is used to plot the first x Prinicipal Components (max of 4).
#' @param pca_obj A Principal Component Analysis object from the prcomp function in the 'stats' package.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Grouping variable(s) used to define plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @param pcs The number of Principal Components to plot.  Values must be between 2 & 4.  Defaults to 4.
#' @param file A file destination/name for the plot output.  File format must be '.pdf'.  If not provided, plot will be rendered in console.
#' @keywords PCA Plot
#' @export
#' @examples
#' # Create a 4 x 4 plot of the of the first 4 Principal Components


pca.plot <- function(pca_obj, gm = NULL, metadata, grouping_var = NULL, gene_highlight = c(), pcs = 4, file, ext = "pdf", marker_size = 0.4){

  ext = hanta:::setup.figure_type(ext)

  if(length(gene_highlight) > 0 & is.null(gm) == T){stop("Must enter gene matrix to highlight genes by mean expression.")}

  if(pcs < 2 | pcs > 4){pcs = 4}

  if(is.null(grouping_var) == TRUE & is.null(gene_highlight) == TRUE){

    f_type = "character"

    groups = as.factor(rep('undefined', nrow(metadata)))

    col_t = structure(hanta:::add.alpha("gray"), names = levels(groups))

  } else {

    # Sort barcodes in PCA data and Metadata by grouping variables.
    g = hanta:::gm.grouper(metadata, grouping_var)

    sort_order = g$sort_order

    groups = g$groups

    pca_obj$x = pca_obj$x[sort_order, ]

    metadata = metadata[sort_order, ]

    col_t = structure(c(rainbow(length(unique(groups)), alpha = 0.5)), names = levels(groups))

    if(length(gene_highlight) > 0){

      cols = match(gene_highlight, rownames(gm))

      col_val = apply(gm[, sort_order], 2, function(x) {

        mean(x[cols])

      })

      cuts = seq(0, ceiling(max(col_val)), length = ceiling(max(col_val)) + 1)

      col_cuts = cut(col_val, breaks = cuts, include.lowest = TRUE)

      groups = as.factor(as.numeric(col_cuts))

      col_t = structure(colorRampPalette(c('grey','red'), alpha = 0.5)(ceiling(max(col_val))), names = levels(groups))

    }

  }

  # Plot PCA
  if(missing(file) == FALSE){

    if(ext == "png"){

      png(pkgmaker::file_extension(file, ".png"), height = 5.69, width = 10.00, units = 'in', res = 600)

    } else {

      pdf(pkgmaker::file_extension(file, ".pdf"), height = 5.69, width = 10.00, onefile = FALSE)

    }

  }

  # Define keys for legend - list of lists containing plotting information
  if(length(gene_highlight) > 0){
    keys = list(space = "right",
                title = if(length(gene_highlight) == 1){
                  "expression"
                } else {
                  "mean expression"
                },
                cex.title = 1,
                points = list(pch = 16,
                              col = col_t),
                text = list(lab = sapply(levels(col_cuts), function(x) gsub(",", " - ", gsub("\\(|\\)|\\]|\\[", "", x)))))
  } else {
    keys = list(space = "right",
                points = list(pch = 16,
                              col = col_t),
                text = list(lab = names(col_t)))
  }

  # Create the PCA graph object
  t <- lattice::splom(as.data.frame(pca_obj$x[ , 1:pcs]),
          main = "Principal Component Analysis",
          key = keys,
          cex = marker_size,
          pch = 16,
          type = 'p',
          xlab = NULL,
          col = col_t[match(groups, names(col_t))],
          scales = list(cex = 0.5),
          par.settings = list(
            par.main.text = list(font = 3),
            layout.heights = list(bottom.padding = 0,
                                  top.padding = 0.5),
            layout.widths = list(left.padding = 0,
                                  right.padding = 0)
          )

  )

  plot(t)

  if(missing(file) == FALSE){

    graphics.off()

  }

}


#' Gene Matrix Principal Component Analysis
#'
#' This function is used to filter a gene matrix, calculate the prinicpal components, and create a PCA plot.  Returns a PCA object from the prcomp function in the 'stats' package.
#' @param gm Gene Matrix.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Variables of interest from metadata file.
#' @param filt_method Filtering method. Choose from from 'quant_exp' (Filter genes below quantile value), 'top_var' (Retain genes with highest variance), or 'top_exp (Retain highest expressing genes).
#' @param thresh_cut Filtering threshold if using 'quant_exp' method. Defaults to 0 (all genes retained).
#' @param top_genes Number of genes to retain if using 'top_var' or 'top_exp' methods. Defaults to the total number of genes (all genes retained).
#' @param transform_type Transformation type from ("none", "ln" or any integer).  Defaults to "none" for no transformation.  Possible options are "ln" for natural log, or an integer indicating the base of the lograthim used in the transformation (ie. For log2 transformed data, enter 2)
#' @param plot A boolean defining whether or not to create a PCA plot. Defaults to FALSE.
#' @param file A file destination/name for the plot output.
#' @param pcs The number of Principal Components to plot (up to 4).
#' @keywords Gene Matrix PCA
#' @export
#' @examples
#' # Perform Principal Component Analysis on 500 highest expressing genes, plot
#' # the results of the first 4 Principal Components.


gm.pca <- function(gm, metadata, grouping_var = NULL, filt_method = FALSE,
                   thresh_cut = 0, top_genes = ncol(gm), transform_type = "none",
                   plot = FALSE, pcs = 4, file, ext = "pdf"){

  ext = hanta:::setup.figure_type(ext)

  # Filter genes by expression, variance or quantile thresholds.
  if(filt_method != FALSE){

    gm = pca.gene_filter(gm = gm, method = filt_method, top_genes = top_genes,
                         thresh_cut = thresh_cut, transform_type = transform_type)

  }

  # Perform Principal Component Analysis
  gm_pca = stats::prcomp(t(gm), center = FALSE, scale. = FALSE)

  # Plot PCA
  if(plot == TRUE){

    pca.plot(pca_obj = gm_pca, metadata = metadata, grouping_var = grouping_var,
             pcs = pcs, file = file, ext = ext)

  }

  pca_res = list(pca_obj = gm_pca, pca_genes = gm)

  return(pca_res)

}
