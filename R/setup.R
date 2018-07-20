#' Gene Matrix Reader
#'
#' Reads gene matrix or metadata files from the ICELL8 Primary Pipeline into R
#' @param file_loc File location of gene matrix, or metadata file (eg. '/path/to/genematrix.csv' or '/path/to/genematrix_report.csv').
#' @param type File type.  Choose from 'gm' (gene matrix), 'metadata' or 'metadata_header'.
#' @keywords Gene Matrix Reader
#' @export

gm.reader <- function(file_loc, type, ICELL8 = TRUE){

  if(!(type %in% c("gm", "metadata", "metadata_header")) == T){

    stop("Must choose either 'gm', 'metadata' or 'metadata_header' for the 'type' parameter.")

  }

  if(file.exists(file_loc) == FALSE){

    stop(paste(file_loc, "does not exist!", sep = " "), call. = FALSE)

  }

  if(type == 'gm'){

    gm = as.matrix(read.csv(file_loc, sep = ","))

    return(gm)

  }

  if(type == 'metadata'){

    if(ICELL8 == TRUE){

      metadata = read.csv(file = file_loc, skip = 8, header = T, stringsAsFactors = F, check.names = F)

      metadata = metadata[-nrow(metadata), ]

    } else {

      metadata = read.csv(file = file_loc, header = T, stringsAsFactors = F, check.names = F)

    }

    return(metadata)

  }

  if(type == 'metadata_header'){

    if(ICELL8 == TRUE){

      metadata_header = read.csv(file = file_loc, nrows = 8, header = F, stringsAsFactors = F, check.names = F)

      colnames(metadata_header) = c("Description", "Value")

    } else {

      metadata_header = NULL

    }

    return(metadata_header)

  }

}


#' Gene Matrix Merger
#'
#' Merge gene matrix & metadata files into large matrices for downstream analysis.  Merge by either the intersection or union of expressed genes.  Returns a list with named elements "gm" and "metadata".
#' @param gms List of gene matrices.
#' @param metadata Matched list of metadata files (genematrix_report.csv files) for each gene matrix in 'gms'.
#' @param merge_type Method type to merge gene matrices.  Options are to merge by the 'union' or 'intersect' of expressed genes.  For the union method, non-observed genes will be replaced with zero read counts.
#' @keywords Gene Matrix Merge
#' @export
#' @examples
#' # Merge a list of gene matrices into one large gene matrix with the union of
#' # expressed genes across gene matrices.


gm.merger <- function(gms, metadata, merge_type){

  if(length(gms) != length(metadata)){

    stop("Number of gene matrix files must match the number of metadata files.")

  }

  if(is.null(names(gms)) == TRUE){

    name_list = as.list(seq(1, length(gms)))

  } else {

    name_list = as.list(names(gms))

  }

  # Filter mdata to bcs found in gm, add gm specific id to bcs in case of duplicates across chips.
  res = mapply(function(X, Y, Z){

    bc = colnames(X)

    mt = Y[match(bc, Y$Barcode), ]

    rownames(mt) = paste(mt$Barcode, Z, sep = "_")

    colnames(X) = paste(mt$Barcode, Z, sep = "_")

    return(list(X, mt))

  }, SIMPLIFY = FALSE, X = gms, Y = metadata, Z = name_list)

  gms = sapply(res, "[", 1)

  metadata = Reduce(rbind, sapply(res, "[", 2))

  if(merge_type == "intersect"){

    g_t = Reduce(intersect, lapply(gms, rownames))

    gms = Reduce(cbind, lapply(gms, function(x) x[g_t, ]))

  }

  if(merge_type == "union"){

    g_t = Reduce(union, lapply(gms, rownames))

    gms = Reduce(cbind, lapply(gms, function(x){

      dif = setdiff(g_t, rownames(x))

      rbind(x, matrix(0, length(dif), ncol(x), dimnames = list(dif)))

    }))

  }

  return(list(gm = gms, metadata = metadata))

}


#' Gene Matrix Setup
#'
#' Reads gene matrix or metadata files from the ICELL8 Primary Pipeline into R and merges multiple output files into a singular object.  Returns a list with named elements 'gm' and 'metadata'.
#' @param gms_loc List of gene matrices.
#' @param metadata_loc Matched list of metadata files (genematrix_report.csv files) for each gene matrix in 'gms'.
#' @param merge_type Method type to merge gene matrices.  Options are to merge by the 'union' or 'intersect' of expressed genes.  For the union method, non-observed genes will be replaced with zero read counts.  Defaults to 'union'.
#' @param transpose TRUE/FALSE Transpose gene matrix.  Gene matrix gm[i,j] must have 'i' rows of genes and 'j' columns of samples.
#' @keywords Gene Matrix Merge
#' @export
#' @examples
#' # Merge a list of gene matrices into one large gene matrix with the union of

gm.setup <- function(gms_loc, metadata_loc, gm_names = NULL, merge_type = "union", transpose = F, ICELL8 = T){

  if(missing(gms_loc) == TRUE){

    stop("Must enter list of gene matrix file(s).")

  }

  if(missing(metadata_loc) == TRUE){

    stop("Must enter list of report file(s).")

  }

  gms = lapply(gms_loc, function(x){

    gm.reader(x, type = "gm", ICELL8 = ICELL8)

  })

  names(gms) <- hanta:::setup.gms_namer(gms_loc = gms_loc, gm_names = gm_names)

  if(transpose == T){

    gms = lapply(gms, t)

  }

  mts = lapply(metadata_loc, function(x){

    gm.reader(x, type = "metadata", ICELL8 = ICELL8)

  })

  mts_headers = lapply(metadata_loc, function(x){

    gm.reader(x, type = "metadata_header", ICELL8 = ICELL8)

  })

  for(i in 1:length(mts)){

    mts[[i]][["group"]] <- rep(names(gms)[i], nrow(mts[[i]]))

  }

  merged = gm.merger(gms = gms, metadata = mts, merge_type = merge_type)

  merged[["mts_headers"]] <- mts_headers

  return(merged)

}


setup.fl_read <- function(ICELL8_data_loc = NULL, gm_file_list = NULL, md_file_list = NULL){

  if(all(is.null(c(ICELL8_data_loc, gm_file_list, md_file_list)))){

    stop("Must enter 'ICELL8_data_loc' option, or 'gm_file_list' & 'gm_report_list' options")

  }

  if(is.null(ICELL8_data_loc) == TRUE){

    gms = setup.fl_exists(file_list = gm_file_list, file_type = "gene matrix list file")

    mds = setup.fl_exists(file_list = md_file_list, file_type = "report list file")

  } else {

    if(file.exists(ICELL8_data_loc) == FALSE){

      stop("ICELL8 data folder not found!")

    }

    gms = setup.fl_exists(file_list = file.path(ICELL8_data_loc, "gm_list"), file_type = "gene matrix list file")

    gms = lapply(gms, function(x){file.path(ICELL8_data_loc, substring(x, 3))})

    mds = setup.fl_exists(file_list = file.path(ICELL8_data_loc, "md_list"), file_type = "metadata list file")

    mds = lapply(mds, function(x){file.path(ICELL8_data_loc, substring(x, 3))})

  }

  if(length(mds) != length(gms)){

    stop("There must be a results file specified for each gene matrix file!")

  }

  return(list(gms = gms, mds = mds))

}

setup.fl_exists <- function(file_list, file_type){

  # check that file list is specified
  if(is.null(file_list) == TRUE){

    stop(paste("Must enter", file_type, "!", sep = " "))

  }

  # check that file list exists
  if(file.exists(file_list) == FALSE){

    stop(paste("Cannot find", file_type, ":", file_list, sep = " "))

  }

  # pull files from file list
  fs = as.list(read.table(file_list, stringsAsFactors = FALSE)$V1)

  return(fs)

}


# parse unique cell types from metadata (report) files
setup.metadata_factor_list <- function(mds, factor, ICELL8 = TRUE){

  # read in metadata (results) files, define sample factors
  cell_types <- Reduce(union, lapply(mds, function(x) {

    f = gm.reader(file_loc = x, type = 'metadata', ICELL8 = ICELL8)

    if(factor %in% c("col", "cols", "colnames")){

      unique(colnames(f))

    } else {

      unique(f[[factor]])

    }

  }))

  return(cell_types)

}

setup.grouping_var_check <- function(mds, grouping_var, ICELL8 = TRUE){

  if(is.null(grouping_var) == FALSE){

    meta_factors = c(hanta:::setup.metadata_factor_list(mds, "colnames", ICELL8), "group")

    if(all(sapply(grouping_var, function(x) x %in% meta_factors) != TRUE)){

      stop(paste("One or more grouping variables\n",
                 paste0(grouping_var, collapse = ", "),
                 "\nnot found.\nPlease select from:\n",
                 paste0(meta_factors, collapse = ", "), sep = " ")
      )

    }

  }

}


# scrub whitespace and commas for input
setup.scrub <- function(v = NULL){

  if(is.null(v) == TRUE){

    return(NULL)

  } else {

    v = as.vector(unlist(strsplit(v, "\\, *|\\,")))

    return(v)

  }

}

setup.gms_namer <- function(gms_loc, gm_names){

  if(is.null(gm_names) == TRUE){

    gm_names = tools::file_path_sans_ext(basename(unlist(gms_loc)))

  } else if(length(gms_loc) != length(gm_names)){

    warning("Number of arguments in --names not equal to number of files in --gm_file_list. Group names defaulted to file names.")

    gm_names = tools::file_path_sans_ext(basename(unlist(gms_loc)))

  }

  return(gm_names)

}


# check & set figure_type
setup.figure_type <- function(ext){

  if(!(tolower(ext) %in% c('pdf', 'png'))){

    stop("Choose either 'pdf' or 'png' for figure type extension.")

  } else {

    return(tolower(ext))

  }

}


# check that cell types exist
setup.cell_type_check <- function(mds, cell_types, ICELL8 = TRUE){

  mts_types = hanta:::setup.metadata_factor_list(mds, "Type", ICELL8 = ICELL8)

  if("all" %in% cell_types){

    cell_types = mts_types

  } else {

    t = setdiff(tolower(cell_types), tolower(c(mts_types)))

    if(length(t) > 0){

      stop(paste("Could not find ", paste0(t, collapse = ", "), ".\nPlease choose from ",
                 toString(mts_types), sep = ""))

    }

  }

  return(cell_types)

}


# check that cor_type matches an acceptable option
setup.cor_check <- function(cor_type){

  if(!(tolower(cor_type) %in% c("pearson", "kendall", "spearman"))){

    stop("Correlation type must be pearson, kendall or spearman!")

  }

}


# check / set number of avaialble cores for parallel processing
setup.core_count <- function(core_count){

  if(is.numeric(core_count) == FALSE){

    stop("Must enter number of parallel cores to run.")

  }

  tot_cores = parallel::detectCores()

  if(core_count > tot_cores){

    stop(paste("The number of selected cores (", core_count,
               ") is greater than the number of available cores (", tot_cores,
               ")", sep = ""))

  } else if(core_count == tot_cores){

    core_count = core_count - 1

  }

  if(core_count == 0){

    stop("Must have at least 1 available core to run.")

  }

  if(core_count > 1){

    WGCNA::allowWGCNAThreads(core_count)

  }

}



# check threshold values for cutoffs
setup.thresh_check <- function(lo, hi, quant_or_abs){

  if(tolower(quant_or_abs) %in% c('quant', 'abs') == FALSE) {

    stop("--quant_or_abs must be either 'quant' or 'abs'")

  }

  if(quant_or_abs == "quant"){

    if(is.numeric(lo) == FALSE | lo < 0.0 | lo > 1.0){
      stop("--lo_thresh option must be between 0 & 1")
    }

    if(is.numeric(hi) == FALSE | hi < 0.0 | hi > 1.0){
      stop("--hi_thresh option must be between 0 & 1")
    }

  }

  if(quant_or_abs == "abs"){

    if(is.numeric(lo) == FALSE | lo < 0.0){

      stop("--lo_thresh option must be >= 0")

    }

    if(is.numeric(hi) == FALSE | hi <= 0.0){

      stop("--hi_thresh option must be > 0")

    }

  }

  if(hi < lo){

    stop("--hi_thresh option must be larger than --lo_thresh option.")

  }

}


# read in genelengths for rpkm normalization
setup.genelengths <- function(genelengths){

  if(is.null(genelengths) == FALSE){

    genelengths = data.frame(read.table(genelengths, header = T,
                                        stringsAsFactors = F, row.names = 1L))

  } else {

    genelengths = NULL

  }

  return(genelengths)

}


# normalization parameters check
setup.norm_check <- function(gm_norm, gm_norm_method, gm_norm_scale,
                             gm_norm_genelengths, opt){

  if(tolower(gm_norm_method) %in% c("tps", "cc_mc", "rpkm") == FALSE){

    stop(paste("Set --gm_norm_method to 'tps' (transcripts per scale)",
               "'cc_mc' (cell coverage divided by median coverage) or ",
               "'rpkm' (reads per kilobase million)", sep = ""))

  }

  if(is.numeric(gm_norm_scale) == FALSE){

    stop("--gm_norm_scale must be numeric.")

  }

  if(tolower(gm_norm_method) == "rpkm" & is.null(gm_norm_genelengths) == TRUE){

    stop(paste0("Table of gene names and gene lengths must be provided for",
                "'rpkm' normalization method"))

  }

  if(gm_norm == FALSE){

    opt$gm_norm_method <- ""

    opt$gm_norm_scale <- ""

    opt$gm_norm_genelengths <- ""

  }

  if(gm_norm_method == 'cc_mc'){

    opt$gm_norm_scale <- ""

    opt$gm_norm_genelengths <- ""

  }

  if(gm_norm_method == 'rpkm'){

    opt$gm_norm_scale <- ""

  }

  return(opt)

}


# transformation parameters check
setup.tform_check <- function(gm_log, gm_log_base, opt){

  if(is.numeric(gm_log_base) == FALSE){

    if(tolower(gm_log_base) != "ln"){

      stop("--gm_log_base option must be numeric or 'ln' for natural log")

    }

  }

  if(gm_log == FALSE){

    opt$gm_log_base = ""

  }

  return(opt)

}



# check gm.pca options
setup.pca_check <- function(pca_filt_method, pca_thresh_cut, pca_top_genes, opt){

  if(is.numeric(pca_top_genes) == FALSE){

    stop("--pca_top_genes must be numeric")

  }

  if(is.numeric(pca_thresh_cut) == FALSE | pca_thresh_cut < 0.0 | pca_thresh_cut > 1.0){

    stop("--pca_top_genes must be numeric & between 0 & 1")

  }

  if(pca_filt_method == ""){

    pca_filt_method = FALSE

    opt$pca_filt_method = FALSE

  } else if(tolower(pca_filt_method) %in% c("quant_exp", "top_var", "top_exp") == FALSE){

    stop(paste("Set --pca_filt_method to one of the following: 'quant_exp'",
               "(Filter genes below quantile value), or 'top_var'",
               "(Pick genes with highest variance), or 'top_exp'",
               "(Pick highest expressing genes).", sep  = " "))

  }

  if(pca_filt_method == FALSE){

    opt$pca_thresh_cut = ""

    opt$pca_top_genes = ""

  }

  if(pca_filt_method %in% c("top_var", "top_exp")){

    opt$pca_thresh_cut = ""

  }

  if(pca_filt_method %in% c("quant_exp")){

    opt$pca_top_genes = ""

  }

  return(list(pca_filt_method = pca_filt_method, opt = opt))

}

# check merge method
setup.merge_method_check <- function(merge_method){

  if(merge_method %in% c("intersect", "union") == FALSE){

    stop("Choose 'intersect' or 'union' for merge method argument.")

  }

}


# check and create output directory
setup.output_dir <- function(output_dir){

  if(dir.exists(output_dir) == FALSE){

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  } else {

    if(basename(output_dir) == "ICELL8_output"){

      warning("ICELL8_output already exists, overwriting current output.")

    } else {

      warning(paste(basename(output_dir), "already exists. Defaulting output to ICELL8_output"))

    }

    output_dir = file.path(getwd(), "ICELL8_output")

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  }

  return(output_dir)

}

# check report style
setup.report_check <- function(report = report){

  # if(tolower(report) %in% c("static", "interactive") == FALSE){
  #
  #   stop("Enter 'static' or 'interactive' for --report option.")
  #
  # }

}


# check pandoc installation and location
setup.pandoc_check <- function(pandoc_loc = Sys.which("pandoc")[1]){

  if(file.exists(pandoc_loc) == TRUE){

    Sys.setenv(RSTUDIO_PANDOC = pandoc_loc)

  } else {

    stop("Pandoc is either not installed or not placed in the system PATH!")

  }

}

