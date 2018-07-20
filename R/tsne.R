#' tSNE Plot
#'
#' Plot the first two dimensions of a tSNE object.
#' @param tsne_obj A t-Distributed Stochastic Neighbor Embedding (tSNE) object from the Rtsne function in the 'Rtsne' package.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Grouping variable(s) used to define plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @param file A file destination/name for the plot output.  File format must be '.pdf'.  If not provided, plot will be rendered in console.
#' @keywords tSNE Plot
#' @export
#' @examples
#' # Create a tSNE plot


tsne.plot <- function(tsne_obj, gm = NULL, metadata, grouping_var = NULL, gene_highlight = NULL, file,
                      gene_highlight_col = c('gray', 'firebrick3'),
                      numeric_gradient_col = c('blue', 'cyan', 'yellow', 'red'),
                      col_override = NULL, ext = "pdf", marker_size = 1){

  ext = hanta:::setup.figure_type(ext)

  if(length(gene_highlight) > 0 & is.null(gm) == T){stop("Must enter gene matrix to highlight genes by mean expression.")}

  if(is.null(grouping_var) == TRUE & is.null(gene_highlight) == TRUE){

    f_type = "character"

    groups = as.factor(rep('undefined', nrow(metadata)))

    col_t = structure(hanta:::add.alpha("gray"), names = levels(groups))

  } else {

    # Sort barcodes in tSNE data and Metadata by grouping variables.
    g = hanta:::gm.grouper(metadata, grouping_var)

    sort_order = g$sort_order

    groups = g$groups

    tsne_obj$Y = tsne_obj$Y[sort_order, ]

    metadata = metadata[sort_order, ]

    f_type = g$f_type

    if(is.null(gene_highlight) == FALSE){

      cols = match(gene_highlight, rownames(gm))

      col_val = apply(gm[, sort_order], 2, function(x) {

        mean(x[cols])

      })

      cuts = seq(0, ceiling(max(col_val)), length = 101)

      col_cuts = cut(col_val, breaks = cuts, include.lowest = TRUE)

      groups = col_cuts

      col_t = structure(hanta:::add.alpha(colorRampPalette(gene_highlight_col)(100)), names = levels(groups))

    } else if(f_type == "numeric"){

      groups = as.numeric(levels(groups))[groups]

      cuts = seq(floor(min(groups)), ceiling(max(groups)), length = 101)

      col_cuts = cut(groups, breaks = cuts, include.lowest = TRUE)

      groups = col_cuts

      col_t = structure(hanta:::add.alpha(colorRampPalette(numeric_gradient_col)(100)), names = levels(groups))

    } else {

      col_t = structure(c(rainbow(length(unique(groups)), alpha = 0.5)), names = unique(levels(groups)))

    }

  }

  if(is.null(col_override) == FALSE){

    col_t = col_override

  }

  # Create tSNE plot.
  if(missing(file) == FALSE){

    if(ext == "png"){

      png(pkgmaker::file_extension(file, ".png"), height = 5.69, width = 10.00, units = 'in', res = 600)

    } else {

      pdf(pkgmaker::file_extension(file, ".pdf"), height = 5.69, width = 10.00, onefile = FALSE)

    }

  }

  par(pty = 's', xpd = TRUE)

  graphics::plot(tsne_obj$Y[,1], tsne_obj$Y[,2],
                 main = "tSNE Analysis",
                 font.main = 3,
                 col = col_t[match(groups, names(col_t))],
                 pch = 16,
                 cex = marker_size,
                 xlab = "tSNE 1",
                 ylab = "tSNE 2"
  )

  if(is.null(gene_highlight) == FALSE){

    hanta:::clust.legend(
      title = if(length(gene_highlight) == 1){
        gene_highlight
      } else {
        "mean exp"
      },
      col = gene_highlight_col,
      gradient = cuts)

  } else if(f_type == "numeric"){

    hanta:::clust.legend(title = paste0(grouping_var, collapse = "."),
                 col = numeric_gradient_col,
                 gradient = cuts)

  } else if(is.null(grouping_var) == FALSE){

    hanta:::clust.legend(title = paste0(grouping_var, collapse = "."), col = col_t)

  }

  if(missing(file) == F){

    graphics.off()

  }

}


#' Gene Matrix tSNE Analysis
#'
#' Perform t-Distributed Stochastic Neighbor Embedding (tSNE) on a gene matrix.
#' @param gm Gene Matrix
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Grouping variable(s) used to define perplexity options & plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @param dims The number of output dimensions. Defaults to 2.
#' @param initial_dims The number of dimensions that should be retained in the initial PCA step. Defaults to 50.
#' @param perplexity Perplexity parameter. Defaults to the average number of samples per group divided by the total number of groups (group defined by the 'grouping_var' option).
#' @param theta Speed/accuracy trade-off (increase for less accuracy), set to 0.0 for exact tSNE. Defaults to 0.25.
#' @param max_iter Number of iterations. Defaults to 1500.
#' @param pca_center TRUE/FALSE, Should data be centered before PCA is applied? Defaults to F.
#' @param pca_scale	TRUE/FALSE, Should data be scaled before PCA is applied? Defaults to F
#' @param plot A boolean defining whether or not to create a PCA plot. Defaults to FALSE.
#' @param file A file destination/name for the plot output (ie. '/path/to/tsne.pdf').
#' @note Defaults for the tSNE parameters have been optimized for general use cases of the ICELL8 platform. These defaults are different from the standard defaults in the Rtsne Package, and may need to be re-optimized for unique applications. Please refer to the original Rtsne package for more information on how to optimize parameters or visit https://distill.pub/2016/misread-tsne/ for more information on optimizing the tSNE analysis.
#' @keywords Gene Matrix tSNE
#' @export
#' @examples
#' # Perform tSNE on the 500 most variable genes


gm.tsne <- function(gm, metadata, grouping_var = NULL, dims = 2, initial_dims = 50,
                    perplexity, theta = 0.25, max_iter = 1500, pca_center = FALSE,
                    pca_scale = FALSE, plot = FALSE, file, ext = FALSE){

  if(missing(perplexity) == TRUE){

    if(is.null(grouping_var) == FALSE){

      cols = match(grouping_var, colnames(metadata))

      groups = apply(metadata, 1, function(x) {

        paste(x[cols], collapse = ".")

      })

      perplexity = mean(table(groups)) / length(table(groups))

    } else {

      perplexity = ifelse(ncol(gm) < 100, ncol(gm) / 4, 50)

    }

  }

  set.seed(1986)

  tsne_obj <- Rtsne::Rtsne(t(as.matrix(gm)),
                dims = dims,
                initital_dims = initial_dims,
                perplexity = perplexity,
                theta = theta,
                max_iter = max_iter,
                pca_center = pca_center,
                pca_scale = pca_scale,
                check_duplicates = FALSE
  )

  rownames(tsne_obj$Y) <- rownames(metadata)

  if(plot == T){

    tsne.plot(tsne_obj = tsne_obj, metadata = metadata, grouping_var = grouping_var, file = file, ext = ext)

  }

  return(tsne_obj)

}


# set number of initial dims for tsne
tsne.initial_dims <- function(..., count_perc_var = 0.95, elbow_thresh = 0.05){

  if(missing(...) == FALSE){

    pca_obj = list(...)[[1]]

    if(is(pca_obj, "prcomp") == TRUE){

      sig_pcs_count = icell8::pca.count(pca_obj, 0.95)

      sig_pcs_elbow = icell8::pca.elbow(pca_obj, 0.05)

      return(sig_pcs_count / sig_pcs_elbow)

    } else if(is.numeric(pca_obj)) {

      return(pca_obj)

    } else {

      return(50)

    }

  } else {

    return(50)

  }

}

clust.legend <- function(title = "", col, gradient = NULL, rounder = FALSE){

  usr <- par("usr")
  xlen <- (usr[2] - usr[1])
  ylen <- (usr[4] - usr[3])
  r <- usr[1] + (xlen * 1.075)
  t <- usr[3] + (ylen * 0.925)

  title <- as.expression(bquote(italic(.(title))))

  if(is.null(gradient) == FALSE){

    xl <- r + (xlen * 0.025)
    yb <- t - (ylen * 0.6)
    xr <- r + (xlen * 0.125)
    yt <- t - (ylen * 0.125)

    plotrix::gradient.rect(xl, yb, xr, yt,
                  col = hanta:::add.alpha(colorRampPalette(col)(100)),
                  border = T, gradient = "y")

    legend(r, t,
           title = title,
           legend = c(""),
           col = "",
           bty= 'n')

    tick_x <- xr
    tick_y <- seq(yb, yt, length.out = 5)
    tick_labels <- seq(min(gradient), max(gradient), length.out = length(tick_y))

    for(i in seq(1, length(tick_y), by = 1)){
      text(tick_x, tick_y[i],
           labels = if(rounder == TRUE){
             hanta:::rounder(tick_labels[i])
             } else {
               format(tick_labels[i], nsmall = 2)
             } , pos = 4)
    }

  } else {

    if(length(col) > 20){

      xl <- r + (xlen * 0.01)
      yb <- t - (ylen * 0.10)
      xr <- r + (xlen * 0.185)
      yt <- t - (ylen * 0.08)

      legend(r, t,
             title = title,
             legend = c(""),
             col = "",
             bty= 'n')

      plotrix::gradient.rect(xl, yb, xr, yt,
                    col = hanta:::add.alpha(colorRampPalette(col)(100)),
                    border = T, gradient = "x")

    } else {

      legend(r, t,
             title = title,
             legend = names(col),
             col = col,
             pch = 16,
             bty = 'n'
      )

    }

  }

}
