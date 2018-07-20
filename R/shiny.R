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
