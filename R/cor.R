# Find the 'q' quantile expression of expressed genes in each barcode of a gene matrix 'gm'. Return expression vector 'ev'.
cor.quant <- function(gm, q){

  if(q != 0){

    ev = apply(gm, 2, function(x){

      quantile(x[x > 0], probs = c(q))[[1]]

    })

  } else {

      ev = rep(0, ncol(gm))

  }

  return(ev)

}


# Flag genes to be exluded from correlation analysis, whose expression is outside of desired thresholds.
cor.filt <- function(gm, lo, hi, quant_or_abs){

  if(missing(quant_or_abs) == TRUE || !(quant_or_abs %in% c("quant", "abs")) == TRUE){

    stop("Choose 'quant' or 'abs' threshold")

  }

  if(quant_or_abs == "quant"){

    lo = cor.quant(gm, lo)

    hi = cor.quant(gm, hi)

  }

  if(quant_or_abs == "abs"){

    lo = lo

    hi = hi

  }

  drops <- mapply(function(X, Y, Z){

    gm[, X] < Y | gm[, X] > Z

  }, X = seq(1, ncol(gm)), Y = lo, Z = hi)

  gm[drops] <- NA

  return(gm)

}


#' Gene Matrix to Correlation Matrix
#'
#' Generate a correlation matrix from a gene matrix.  Optional filtering parameters allow user to set low and high thresholds for correlation.  User can also return a 1-sided distribution of correlation values instead of the entire matrix.
#' @param gm Gene Matrix.
#' @param cor_type Correlation type. Choose from 'pearson', 'spearman' or 'kendall'. Defaults to 'pearson'.
#' @param quant_or_abs Choose 'quant' for quantile or 'abs' for absolute threshold filtering.
#' @param lo Low Filter, can be absolute or quantile based on the 'quant_or_abs' parameter.
#' @param hi High Filter, can be absolute or quantile based on the 'quant_or_abs' parameter.
#' @param cor_dist TRUE/FALSE. Return 1-sided correlation distribution instead of Correlation Matrix.
#' @keywords Gene Matrix Correlation Matrix
#' @export
#' @examples
#'

gm.cormatrix <- function(gm, cor_type, quant_or_abs = 'quant', lo = 0, hi = 1.0, cor_dist = FALSE){

  if(!(tolower(cor_type) %in% c("pearson", "spearman", "kendall"))){

    stop("Choose correlation type (cor_type) of 'pearson', 'spearman' or 'kendall'")

  }

  if(!(tolower(quant_or_abs) %in% c("quant", "abs"))){

    stop("Choose filter type (quant_or_abs) of 'quant' for filtering based on quantile threshold or 'abs' for filtering based on the absoulte expression levels.")

  }

  if(quant_or_abs == 'quant' & (lo > 1.0 | lo < 0.0 | hi > 1.0 | hi < 0.0 | lo > hi)){

    stop("Low Threshold 'lo' and High Threshold 'hi' must be between 0.0 & 1.0. Low Threshold must be lower than High Threshold.")

  }

  if(quant_or_abs == 'abs' & (lo < 0.0 | hi < 0.0 | lo > hi)){

    stop("Low Threshold 'lo' and High Threshold 'hi' must be greater than 0.0. Low Threshold must be lower than High Threshold.")

  }

  # Filter gene matrix, by user defined thresholds
  gm = cor.filt(gm, lo, hi, quant_or_abs)

  # Calculate correlation matrix (cm)
  cm = WGCNA::cor(gm, use = "pairwise.complete.obs", method = cor_type)

  if(cor_dist == F){

    return(cm)

  } else {

    corr_dist = cm[upper.tri(cm)]

    return(corr_dist)

  }

}


#' Correlation Matrix Heatmap
#'
#' Generates a heatmap from a correlation matrix.
#' @param cm Correlation matrix.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param cor_type Correlation type. Choose from 'pearson', 'spearman' or 'kendall'. Defaults to 'pearson'.
#' @param grouping_var Grouping variable(s) used to define plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @param cluster TRUE/FALSE Perform clustering.
#' @param label TRUE/FALSE Label rows and columns with sample names.
#' @param rel_scale TRUE/FALSE relative scale. If FALSE, color scale will range from 0:1 or -1:1 depending on the minimum value of the correlation matrix.  If TRUE color scale will range from the mimimum of the correlation matrix to 1.
#' @param palette A color palette selected from the brewer.pal function in the RColorBrewer package in R.  Defaults to "RdBu".  See 'http://www.datavis.ca/sasmac/brewerpal.html' for available options.
#' @param file A file destination/name for the plot output.  File format must be '.pdf'.  If not provided, plot will be rendered in console.
#' @param ... Optional named parameters to be sent to the aheatmap function in the NMF package.
#' @keywords Gene Matrix Correlation Matrix
#' @importFrom NMF aheatmap
#' @export
#' @examples
#'

cor.heatmap <- function(cm, gm = NULL, metadata, cor_type = "", grouping_var = NULL,
                        gene_highlight = c(), gene_highlight_multi = FALSE,
                        cluster = FALSE, label = FALSE, rel_scale = FALSE,
                        palette = "RdGy", palette_rev = FALSE, file, ext = "pdf",
                        plot_height = 6.69, plot_width = 10.00, plot_res = 600,
                        legend_title = "group", ...){

  ext = hanta:::setup.figure_type(ext)

  if(palette %in% rownames(RColorBrewer::brewer.pal.info) == FALSE){

    stop(paste("Must select heatmap palette from:", paste0(rownames(RColorBrewer::brewer.pal.info), collapse = ", "), sep = " "))

  }

  if(length(gene_highlight) > 0 & is.null(gm) == T){stop("Must enter gene matrix to highlight genes by mean expression.")}

  if(length(gene_highlight) > 4 & gene_highlight_multi == TRUE){stop("May only enter 4 or less genes/groups under the gene_highlight option.")}

  if(is.null(grouping_var) == FALSE){

    # Sort correlation matrix by grouping variable(s)
    g = hanta:::gm.grouper(metadata, grouping_var)

    sort_order = g$sort_order

    groups = g$groups

    if(is.null(gm) == FALSE) {

      gm = gm[, sort_order]

    }

    cm = cm[sort_order, sort_order]

    metadata = metadata[sort_order, ]

    if(typeof(legend_title) != "character"){legend_title = " "}
    if(legend_title == ""){legend_title = " "}

    heatmap_df <- setNames(data.frame(row.names = rownames(metadata), legend_title = groups), legend_title)

    col_t = structure(c(rainbow(length(unique(groups)), alpha = 0.5)), names=unique(levels(groups)))

    heatmap_colors <- setNames(list(legend_title = col_t), legend_title)

  }

  if(length(gene_highlight) > 0){

    if(gene_highlight_multi == TRUE){

      for(i in 1:length(gene_highlight)){

        cols = match(gene_highlight[[i]], rownames(gm))

        col_val = apply(gm, 2, function(x) {

          mean(x[cols])

        })

        heatmap_df[[ifelse(length(gene_highlight[[i]]) > 1, paste("mean marker exp. | group ", i, sep = ""), gene_highlight[[i]])]] <- col_val

      }

    } else {

      cols = match(gene_highlight, rownames(gm))

      col_val = apply(gm, 2, function(x) {

        mean(x[cols])

      })

      heatmap_df[[ifelse(length(gene_highlight) > 1, "mean marker exp.", gene_highlight)]] <- col_val

    }

  }

  # Define color ramp
  pal_length = RColorBrewer::brewer.pal.info[palette, 1]

  col_r = grDevices::colorRampPalette(RColorBrewer::brewer.pal(pal_length, palette))(50)

  if(palette_rev == TRUE){

    col_r = rev(col_r)

  }

  if(rel_scale == FALSE & min(cm) >= 0) {

    col_r = col_r[26:50]

  }

  if(missing(file) == F){

    if(ext == "png"){

      png(pkgmaker::file_extension(file, ".png"), height = plot_height, width = plot_width, units = 'in', res = plot_res)

    } else {

      pdf(pkgmaker::file_extension(file, ".pdf"), height = plot_height, width = plot_width, onefile = FALSE)

    }

  }

  # Default args for heatmap
  aheatmap_args = list(x <- cm,
                       color = col_r,
                       breaks = if(rel_scale == TRUE){
                         seq(min(cm), 1, length.out = 51)
                       } else if(min(cm) >= 0) {
                         seq(0, 1, length.out = 51)
                       } else {
                         seq(-1, 1, length.out = 51)
                       },
                       na.color = 'red',
                       border_color = NA,
                       cellwidth = (0.75 * 451) / nrow(cm), # Not perfect, but size ~ same regardless of cm dim.
                       cellheight = (0.75 * 451) / nrow(cm),
                       scale = "none",
                       Rowv = if(cluster == T){T} else {NA},
                       Colv = if(cluster == T){T} else {NA},
                       labRow = if(label == T){rownames(cm)} else {NA},
                       labCol = if(label == T){colnames(cm)} else {NA},
                       cexRow = min(0.2 + 1/log10(ncol(cm)), 0.3),
                       cexCol = min(0.2 + 1/log10(ncol(cm)), 0.3),
                       fontsize = 10,
                       annRow = if(is.null(grouping_var) == FALSE){heatmap_df} else {NA},
                       annCol = if(is.null(grouping_var) == FALSE){heatmap_df} else {NA},
                       annColors = if(is.null(grouping_var) == FALSE){heatmap_colors} else {NA},
                       annLegend = if(is.null(grouping_var) == FALSE){
                         if(length(unique(groups)) > 32){
                            warning("Number of groups exceeds plot area, legend set to FALSE")
                            FALSE
                         } else {
                            TRUE
                         }} else {FALSE},
                       main = if(cor_type == "spearman"){
                         expression(paste(italic("Spearman's Correlation Matrix Heatmap ("), italic(rho), italic(")")))
                       } else if (cor_type == "kendall"){
                         expression(paste(italic("Kendall's Correlation Matrix Heatmap ("), italic(tau), italic(")")))
                       } else if (cor_type == "pearson"){
                         expression(paste(italic("Pearson's Correlation Matrix Heatmap (r)")))
                       } else {expression(paste(italic("Correlation Matrix Heatmap")))
                       },
                       verbose = FALSE,
                       trace = FALSE
  )

  # Find named arguments from input and overwrite default arguments.
  aheatmap_custom_args = list(...)

  aheatmap_args[names(aheatmap_custom_args)] <- aheatmap_custom_args

  options(warn = -1)

  do.call(NMF::aheatmap, aheatmap_args)

  options(warn = 0)

  if(missing(file) == FALSE){

    graphics.off()

  }

}


#' Correlation Matrix Parser
#'
#' Parse a correlation matrix to find inter and intra-group distributions.  Returns a list with named elements intracor_list & intercor_list
#' @param cm Correlation Matrix.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Grouping variable(s) used to define groups within correlation matrix.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @keywords Gene Matrix Correlation Matrix
#' @export
#' @examples
#' # Parse a correlation matrix for inter and intra-group distributions.
#'

cor.parser <- function(cm, metadata = NULL, grouping_var = NULL){

  if(is.null(grouping_var) == TRUE){

    intracor_list = list()

    intracor_list[["total"]] = as.vector(cm[upper.tri(cm)])

    cor_data = list(intracor_list = intracor_list, intercor_list = NULL)

    return(cor_data)

  } else {

    # Sort the correlation matrix by grouping variable(s)
    g = hanta:::gm.grouper(metadata, grouping_var)

    sort_order = g$sort_order

    groups = g$groups

    cm = cm[sort_order, sort_order]

    metadata = metadata[sort_order, ]

    intracor_list = list()

    intercor_list = list()

    for(i in unique(groups)){

      i_labels <- rownames(subset(metadata, groups == toString(i)))

      for(j in unique(groups)){

        j_labels <- rownames(subset(metadata, groups == toString(j)))

        if(j >= 1){

          mm <- cm[i_labels, j_labels]

          if(i == j){

            if(length(mm) <= 1){

              intracor_list[[paste(i, j, sep = "_")]] <- NA

            } else {

              corr_dist = mm[upper.tri(mm)]

              intracor_list[[paste(i)]] <- corr_dist

            }

          } else {

            intercor_list[[paste(i, j, sep = "_")]] <- as.vector(mm)

          }

        }

      }

    }

    cor_data = list(intracor_list = intracor_list, intercor_list = intercor_list)

    return(cor_data)

  }

}


#' Gene Matrix Correlation Distribution Boxplot
#'
#' This function is used to create a boxplot of the well-to-well correlation across different groups.
#' @param cor_list List of well-to-well correaltion distributions.  Recommended to use the output list(s) from the cor.parser function.
#' @param cor_type Correlation type.  Choose from 'pearson', 'spearman' or 'kendall'.
#' @param file A file destination/name for the plot output.  File format must be '.pdf'.
#' @keywords Gene Matrix Correlation Boxplot
#' @export
#' @examples
#'

cor.dist_boxplot<-function(cor_list, cor_type, file, ext = "pdf"){

  ext = hanta:::setup.figure_type(ext)

  cor_list = hanta:::na.omit_list(cor_list)

  col_t <- c(rainbow(length(cor_list), alpha = 0.5))

  if(missing(file) == FALSE){

    if(ext == "png"){

      png(pkgmaker::file_extension(file, ".png"), height = 5.69, width = 10.00, units = 'in', res = 600)

    } else {

      pdf(pkgmaker::file_extension(file, ".pdf"), height = 5.69, width = 10.00, onefile = FALSE)

    }

  }

  par(mar = par("mar") + c(7, 0, 0, 0))

  boxplot(cor_list,
          main = "",
          outcol = rgb(0, 0, 0, 0.5),
          outcex = 0.5,
          col = col_t,
          ylab = if(cor_type == "spearman"){expression(paste("Spearman's Correlation (", rho, ")"))
          } else if (cor_type == "kendall"){expression(paste("Kendall's Correlation (", tau, ")"))
          } else if (cor_type == "pearson"){expression(paste("Pearson's Correlation (r)"))
          } else {"Correlation"},
          ylim = c(min(unlist(cor_list)) * 0.9, 1.0),
          las = 2)

  title(main = "Intragroup Well to Well Correlation", font.main = 3)

  if(missing(file) == FALSE){

    graphics.off()

  }

}


#' Gene Matrix Correlation Distribution Statistics Table
#'
#' This function is used to create a table of summary statistics for the distribution of well-to-well correlations across different groups.
#' @param cor_list List of well-to-well correaltion distributions. Recommended to use the output list(s) from the cm_parser function.
#' @param file Location to save table in .csv format (eg. "path/to/intracor_stats.csv")
#' @keywords Gene Matrix Correlation Boxplot
#' @export
#' @examples
#'

cor.stat_table <- function(cor_list, cor_dist_obj = FALSE, file){

  t = hanta:::na.omit_list(cor_list)

  if(cor_dist_obj == TRUE){

    t = unlist(t, recursive = F)

  }

  tbl = lapply(t, summary)

  if(length(tbl) == 1){

    tbl = data.frame(matrix(Reduce(rbind, tbl), byrow = T, nrow = 1), row.names = names(tbl))

  } else{

    tbl = data.frame(Reduce(rbind, tbl), row.names = names(tbl))

  }

  colnames(tbl) = c("Min", "1st Qu", "Median", "Mean", "3rd Qu", "Max")

  tbl[] = lapply(tbl, sprintf, fmt = "%.2f")

  if(missing(file) == FALSE){

    write.csv(tbl, file = pkgmaker::file_extension(file, ".csv"), row.names = T, quote = F)

  }

  return(tbl)

}


#' Gene Matrix Correlation Analysis
#'
#' Run correlation analysis.
#' @param gm Gene matrix.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param cor_type Correlation type. Choose from 'pearson', 'spearman' or 'kendall'. Defaults to 'pearson'.
#' @param quant_or_abs Choose 'quant' for quantile or 'abs' for absolute threshold filtering.
#' @param lo Low Filter, can be absolute or quantile based on the 'quant_or_abs' parameter.
#' @param hi High Filter, can be absolute or quantile based on the 'quant_or_abs' parameter.
#' @param cor_dist TRUE/FALSE. Return 1-sided correlation distribution instead of Correlation Matrix.
#' @param grouping_var Grouping variable(s) used to define plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @param heatmap TRUE/FALSE Generate heatmap.  Defaults to FALSE.
#' @param dist_boxplot TRUE/FALSE Generate correlation distribution boxplot.  Defaults to FALSE.
#' @param dist_table TRUE/FALSE Generate correlation distribution summary statistics table.  Defaults to FALSE.
#' @param heatmap_cluster TRUE/FALSE Perform clustering on heatmap.  Defaults to FALSE.
#' @param heatmap_label TRUE/FALSE Label rows and columns with sample names on heatmap.  Defaults to FALSE.
#' @param ... Optional named parameters to be sent to the aheatmap function in the NMF package.
#' @keywords Gene Matrix Correlation Boxplot
#' @export
#' @examples
#'
#'

gm.cor_analysis <- function(gm, metadata, cor_type, quant_or_abs = 'quant',
                            lo = 0.0, hi = 1.0, grouping_var = NULL,
                            heatmap = FALSE, dist_boxplot = FALSE,
                            dist_table = FALSE, heatmap_cluster = FALSE,
                            heatmap_label = FALSE, ext = "pdf",
                            output_dir = getwd(), ...){

  ext = hanta:::setup.figure_type(ext)

  cm = gm.cormatrix(gm = gm, cor_type = cor_type, quant_or_abs = quant_or_abs,
                    lo = lo, hi = hi, cor_dist = FALSE)

  if(heatmap == TRUE){

    cor.heatmap(cm = cm, metadata = metadata, cor_type = cor_type,
                grouping_var = grouping_var, label = heatmap_label,
                cluster = heatmap_cluster,
                file = file.path(output_dir, "heatmap.pdf"), ext = ext, ...)

  }

  cor_dist = cor.parser(cm = cm, metadata = metadata, grouping_var = grouping_var)

  if(dist_boxplot == TRUE){

    cor.dist_boxplot(cor_list = cor_dist$intracor_list, cor_type = cor_type,
                     file = file.path(output_dir, "boxplot.pdf"), ext = ext)

  }

  if(dist_table == TRUE){

    cor_stat_table = cor.stat_table(cor_dist, cor_dist_obj = TRUE,
                                    file = file.path(output_dir, "cor_stats.csv"))

  }

  cor_analysis = list(cm = cm, gm = gm, metadata = metadata, cor_dist = cor_dist,
                      cor_type = cor_type, grouping_var = grouping_var)

  return(cor_analysis)

}

