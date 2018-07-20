#' Gene Matrix Log Transform
#'
#' Log transforms expression values in a gene matrix.  Adds a pseudocount of 1 to each value in the matrix to avoid undefined log(0).
#' @param gm Gene matrix.
#' @param log_base Logarithmic base. Value must be numeric or 'ln' for natural log. Defaults to 'ln' for natural log.
#' @keywords Gene Matrix Log Transform
#' @export
#' @examples
#' # Log10 transform a gene matrix
#' gm_log <- gm.log(raw_data$gm, 10)

gm.log <- function(gm, log_base = 'ln'){

  if(tolower(log_base) == 'ln'){

    gm = log(gm + 1, base = exp(1))

  } else if(is.numeric(as.numeric(log_base)) == TRUE){

    gm = log(gm + 1, base = as.numeric(log_base))

  } else {

    stop("Unsupported log_base option")

  }

  return(gm)

}


#' QC Median Absolute Deviations
#'
#' Calculates lower and upper bounds of 'n' median absolute deviations from median of vector.  Used to find outlier cells or genes.  Returns named list with elements 'lo' & 'hi'.
#' @param v Vector of values (eg. vector of gene coverage per cell)
#' @param n Number of median absolute deviations from median.
#' @keywords Gene Matrix Quality Control Median Absolute Deviation
#' @export
#' @examples
#' # Calculate the threshold for dropping cells with lower than 3 Median
#' # Absolute Deviations from median overall coverage.


qc.mad <- function(v, n){

  mad = median(abs(v - median(v)))

  lo = median(v) - (n * mad)

  hi = median(v) + (n * mad)

  return(list(lo = lo, hi = hi))

}


#' QC Low Cell Coverage Flagger
#'
#' Flag cells to drop with low overall coverage.
#' @param gm Gene matrix.
#' @param n Number of median absolute deviations from median.
#' @param abs_thresh Absolute Threshold, set minimum read depth per cell. Defaults to -Inf reads.
#' @keywords Gene Matrix Quality Control Cell Coverage Dropper
#' @export
#' @examples
#' # Flag cells to drop that are lower than 3 Median Absolute Deviations from


# Flag cells to drop with low overall coverage (less than user-defined Median Absolute Deviations below the median, or the absolute threshold)
qc.flag_cells_lowcov <- function(gm, n, abs_thresh = -Inf){

  cell_cov = colSums(gm)

  cell_cov.drop = cell_cov < qc.mad(cell_cov, n)$lo | cell_cov < abs_thresh

  return(cell_cov.drop)

}


#' QC Low Total Gene Count Flagger
#'
#' This function is used to flag cells to drop with too few total genes expressed.
#' @param gm Gene matrix.
#' @param n Number of median absolute deviations from median.
#' @param abs_thresh Absolute Threshold, set minimum gene counts per cell. Defaults to -Inf genes.
#' @keywords Gene Matrix Quality Control Low Total Gene Count
#' @export
#' @examples
#' # Flag cells to drop that have less than 3 Median Absolute Deviations from
#' # median total gene count


# Flag cells to drop with too few total genes expressed (less than user-defined Median Absolute Deviations below the median, or the absolute threshold)
qc.flag_cells_lowgenecount <- function(gm, n, abs_thresh = -Inf){

  cell_genecount = colSums(gm != 0)

  cell_genecount.drop = cell_genecount < qc.mad(cell_genecount, n)$lo | cell_genecount < abs_thresh

  return(cell_genecount.drop)

}


#' QC Low Total Cell Count Flagger
#'
#' This function is used to flag genes to drop that are found in too few number of cells.
#' @param gm Gene matrix
#' @param cell_count The required number of cells in which a given gene has > 0 expression.  If a gene is expressed in less than the number of cells defined in cell_count, the gene is flagged to be dropped from the matrix.
#' @keywords Gene Matrix Quality Control Low Total Cell Count
#' @export
#' @examples
#'

# Flag genes to drop that are expressed in fewer than x number of cells
qc.flag_genes_lowcellcount <- function(gm, cell_count){

  gene_count = rowSums(gm != 0)

  gene_count.drop = gene_count <= cell_count

  return(gene_count.drop)

}


#' QC Low Gene Expression Flagger
#'
#' This function is used to flag genes to drop that have low overall expression.
#' @param gm Gene matrix.
#' @param exp_cutoff Flag genes to drop that have < than the defined expression cutoff.
#' @keywords Gene Matrix Quality Control Low Expression Genes
#' @export
#' @examples
#'

# Flag genes to drop with poor overall expression
qc.flag_genes_lowcov <- function(gm, exp_cutoff){

  gene_cov = rowSums(gm)

  gene_cov.drop = gene_cov <= exp_cutoff

  return(gene_cov.drop)

}


#' QC Normalization
#'
#' This function is used to normalize gene matrices.
#' @param gm Gene matrix.
#' @param method Normalize to 'tps' (Transcript Per Scale) or 'cc_mc' (cell coverage by median coverage). If 'tps' is chosen, the 'scale' paramter must be set. To normalize each cell to Transcripts Per Million choose method = 'tps', scale = 1e6.
#' @param scale The scale with which to normalize read counts when using the 'tps' method. Defaults to Transcripts per Million (TPM), scale = 1e6.
#' @keywords Gene Matrix Normalization
#' @export
#' @examples
#' # Normalize Gene Matrix to Transcripts per Million (TPM)


# Normalize matrix to specified scale (default Transcripts per Million)
gm.norm <- function(gm, method, scale = 1e6, genelengths = NULL){

  if(missing(method) == TRUE){

    stop("Specify normalization method. (method = 'tps' (transcripts per scale) or 'cc_mc' (cell coverage divided by median coverage)")

  }

  if(!(method %in% c("tps", "cc_mc", "rpkm")) == TRUE){

    stop("Set method to 'tps' (transcripts per scale), 'cc_mc' (cell coverage divided by median coverage) or 'rpkm' (reads per kilobase million)")

  }

  if(is.numeric(scale) == F){stop("Select a numeric value for scale.")}

  if(method == 'tps'){

    norm_factor = colSums(gm) / scale

    gm = sweep(gm, 2, norm_factor, FUN = '/')

    gm[is.na(gm)] <- 0

  }

  if(method == 'cc_mc'){

    cell_cov = colSums(gm)

    norm_factor = cell_cov / median(cell_cov)

    gm = sweep(gm, 2, norm_factor, FUN = '/')

    gm[is.na(gm)] <- 0

  }

  if(method == 'rpkm'){

    if(is.null(genelengths) == TRUE){

      stop("Must enter gene length table ('gm_norm_genelengths' option) when using 'rpkm' normalization.")

    }

    gene.keep = intersect(row.names(gm), row.names(genelengths))

    if(length(gene.keep) == 0){

      stop("Genes in gene matrix and gene length table ('gm_norm_genelengths' option) do not match.")

    }

    if(length(gene.keep) < nrow(gm)){

      warning(paste0(nrow(gm) - length(gene.keep),
                     " gene(s) removed from gene matrix as gene length(s) from",
                     " gene length table."))

    }

    glpk = matrix(genelengths[gene.keep, ] / 1000, nrow = 1)

    nf = matrix(colSums(gm) / 1e6, nrow = 1)

    t = crossprod(glpk, nf)

    gm = gm / t

  }

  return(gm)

}

# retain only specified cell types from gm & metadata
qc.cell_type_filter <- function(gm, metadata, cell_types = "all"){

  if("all" %in% tolower(cell_types)){

    res = list(gm = gm, metadata = metadata)

  } else {

    bc.flag = tolower(metadata$Type) %in% tolower(as.vector(cell_types))

    gm = gm[, bc.flag]

    metadata = metadata[bc.flag, ]

    res = list(gm = gm, metadata = metadata)

  }

  return(res)

}



#' Gene Matrix Quality Control (gmqc)
#'
#' This function reads in a gene matrix, filters based on quality control criteria and applies optional normalization and log transformation techniques
#' @param gm Gene matrix. Optionally, include the raw file location in the 'gm_loc' option.
#' @param gm_loc Gene matrix file location (eg. /path/to/genematrix.csv).
#' @param gm_qcfilt TRUE/FALSE, perform QC filtering.
#' @param gm_norm TRUE/FALSE, normalize gene matrix.
#' @param gm_norm_method Normalize to 'tps' (Transcript Per Scale) or 'cc_mc' (cell coverage by median coverage).  If 'tps' is chosen, the 'gm_norm_scale' paramter must be set.  To normalize each cell to Transcripts Per Million choose method = 'tps', gm_norm_scale = 1e6.
#' @param gm_norm_scale The scale with which to normalize read counts when using the 'tps' method.  Defaults to Transcripts per Million (TPM), gm_norm_scale = 1e6.
#' @param gm_log TRUE/FALSE, perform log transform on gene matrix.
#' @param gm_log_base Logrithmic base for transformation.  Choose integer or 'ln' for natural logrithm.  Default is gm_log_base = 'ln'.
#' @param mads Number of median absolute deviations from median for QC filtering. Default is mads = 3.
#' @param cell_abslowcov Absolute Threshold, set minimum read depth per cell. Defaults to 10,000 reads.
#' @param cells_abslowgenecount Absolute Threshold, set minimum gene counts per cell. Defaults to 3 genes
#' @param gene_cellcount The number of cells a gene must have > 0 expression to be kept in the gene matrix.  Default is gene_cellcount = 1.
#' @param gene_totcov The minimum number of reads for a gene across all cells. Default is gene_totcov = 5.
#' @keywords Gene Matrix Normalization
#' @export
#' @examples
#' # Perform quality control and normalization on a gene matrix file.

gm.qc <- function(gm, metadata = NULL, gm_loc = NULL, gm_qcfilt = TRUE,
                  gm_norm = FALSE, gm_norm_method = "tps", gm_norm_scale = 1e6,
                  gm_norm_genelengths = NULL, gm_log = FALSE, gm_log_base = 'ln',
                  mads = 3, cell_abslowcov = 1e4, cells_abslowgenecount = 3,
                  gene_cellcount = 1, gene_totcov = 5){

  if(missing(gm) == TRUE){

    if(is.null(gm_loc) == TRUE){

      stop("Must enter a gene matrix object in the 'gm' option, or the gene matrix file location in the 'gm_loc' option.")

    } else {

      gm = gm.reader(gm_loc)

    }

  }

  if(gm_qcfilt == TRUE){

    gm = gm[!(qc.flag_genes_lowcellcount(gm, cell_count = gene_cellcount) | qc.flag_genes_lowcov(gm, exp_cutoff = gene_totcov)),
            !(qc.flag_cells_lowcov(gm, n = mads, abs_thresh = cell_abslowcov) | qc.flag_cells_lowgenecount(gm, n = mads, abs_thresh = cells_abslowgenecount))]

  }

  if(gm_norm == TRUE){

    gm = gm.norm(gm, gm_norm_method, gm_norm_scale, gm_norm_genelengths)

  }

  if(gm_log == TRUE){

    gm = gm.log(gm, gm_log_base)

  }

  if(is.null(metadata) == FALSE){

    metadata = metadata[match(colnames(gm), rownames(metadata)), ]

  }

  qc_data = list(gm = gm, metadata = metadata)

  return(qc_data)

}

