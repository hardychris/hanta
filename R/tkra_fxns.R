#' Library Load
#'
#' This function suppresses output when loading libraries.
#' @param libs List of libraries for loading.
#' @keywords General utility
#'

lib_load <- function(libs){

  os = .Platform$OS.type

  if(os == 'unix'){
    for(i in 1:length(libs)){
      capture.output(suppressWarnings(suppressMessages(base::library(libs[[i]], character.only = T))), file='/dev/null')
    }
  }

  if(os == 'windows'){
    for(i in 1:length(libs)){
      capture.output(suppressWarnings(suppressMessages(base::library(libs[[i]], character.only = T))), file='NUL')
    }
  }

}


#' gm.grouper
#'
#' Group variables in metadata file
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Grouping variable(s) used to define plot color scheme.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#'

gm.grouper <- function(metadata, grouping_var){

  if(is.null(grouping_var) == TRUE){

    grouping_var = "Barcode"

  }

  if(any(grouping_var %in% colnames(metadata)) == FALSE){

    stop("One or more grouping variables not found in metadata header.")

  }

  cols = match(grouping_var, colnames(metadata))

  f_type = "character"

  if(length(cols) == 1){

    ft = sapply(as.data.frame(metadata), is.numeric)

    if(ft[cols] == TRUE){

      f_type = "numeric"

    } else {

      f_type = "character"

    }

  }

  groups = do.call(paste, c(metadata[cols], list(sep = ".")))

  sort_order = naturalsort::naturalorder(groups)

  groups = as.factor(groups[sort_order])

  return(list(groups = groups, sort_order = sort_order, f_type = f_type))

}


#' Remove NA's from Lists
#'
#' This function removes NA's from lists
#' @param y List.
#' @keywords General utility
#'

na.omit_list <- function(y) {

  y <- y[!sapply(y, function(x) all(is.null(x)))]

  return(y[!sapply(y, function(x) all(is.na(x)))])

}

exp_summary <- function(metadata_loc = ''){

  metadata = gm.reader(metadata_loc, type = 'metadata')

  metadata_header = gm.reader(metadata_loc, type = 'metadata_header')

  sample_name = metadata_header[1, 2]
  sample_description = metadata_header[2, 2]
  total_reads = as.integer(metadata_header[3, 2])
  valid_reads = as.integer(metadata_header[6, 2])
  fraction_valid = valid_reads / total_reads
  barcode_count = nrow(metadata)
  valid_reads_per_bc = valid_reads / barcode_count

  res = data.frame(
    list(
      `Sample Name` = sample_name,
      `Sample Description` = sample_description,
      `Total Reads` = total_reads,
      `Valid Reads` = valid_reads,
      `Fraction Valid` = fraction_valid,
      `# of Barcodes` = barcode_count,
      `Valid Reads / Barcode` = valid_reads_per_bc
    ), check.names = FALSE
  )

  kable(t(res), format = 'latex', booktabs = TRUE)

}


#' Run tsne
#'
#' Run tsne shiny
#' @export
#'
tsne.shiny <- function(cluster_data = NULL,
                       gm = NULL, gm_top = NULL, metadata = NULL, tsne_obj = NULL,
                       app_loc = file.path(.libPaths(), "icell8", "shiny_apps", "tsne"), ...) {

  if(is.null(cluster_data) == FALSE){

    gm <<- gm
    gm_top <<- gm_top
    metadata <<- metadata
    tsne_obj <<- tsne_obj
    p_default <<- tsne_obj$perplexity

  } else {

    gm <<- cluster_data$gm
    gm_top <<- cluster_data$gms_top_genes
    metadata <<- cluster_data$metadata
    tsne_obj <<- cluster_data$tsne
    p_default <<- tsne_obj$perplexity

  }

  shiny::runApp(appDir = app_loc, ...)

}

add.alpha <- function(col, alpha = 0.5){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, function(x) {
    rgb(x[1], x[2], x[3], alpha = alpha)
  })
}


rounder <- function(number, rounding = F){

  lut <- c(1e-24, 1e-21, 1e-18, 1e-15, 1e-12, 1e-09, 1e-06,
           0.001, 1, 1000, 1e+06, 1e+09, 1e+12, 1e+15, 1e+18, 1e+21,
           1e+24)

  pre <- c("y", "z", "a", "f", "p", "n", "u", "m", "", "K",
           "M", "G", "T", "P", "E", "Z", "Y")

  ix <- findInterval(number, lut)

  print(lut[ix])

  if(lut[ix] != 1){
    if(rounding == T){
      out <- paste(sprintf("%.2f", round(number/lut[ix])), pre[ix], sep = "")
    } else {
      out <- paste(sprintf("%.2f", number/lut[ix]), pre[ix], sep = "")
    }
  }
  else {
    out <- as.character(number)
  }
  return(out)
}


find.file <- function(){

  list.files(path = "", pattern = "$pandoc^", recursive = TRUE, full.names= TRUE)

}

is.empty <- function(x){
  if(is.null(x) == TRUE){return(TRUE)}
  else if(is.na(x) == TRUE){return(TRUE)}
  else if(x == ""){return(TRUE)}
  else {return(FALSE)}
}





