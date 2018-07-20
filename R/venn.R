#' Gene Matrix Venn Diagram
#'
#' Generate a venn diagram of gene counts for groups in a gene matrix.
#' @param gm Gene Matrix.
#' @param metadata Metadata table with columns of named features and rows of unique barcodes.
#' @param grouping_var Grouping variable(s) used to define groups within correlation matrix.  Must be vector where each element matches a column header in the metadata table (ie. c("Type", "Group"))
#' @param cell_no Cell Number.  The minimum number of cells a gene must be expressed in to be counted within the specified group.  Defaults to 1.
#' @param cell_thresh Cell threshold.  The minimum absolute expression level required for a gene to be counted within the specified group.  Defaults to 0 (ie. All expressed genes).
#' @param raw_or_percent Choose between c("raw"), c("percent"), or c("raw", "percent"), to display the given statistic type on the venn diagram.  Defaults to c("raw", "percent").
#' @param file A file destination/name for the plot output.  File format must be '.pdf'.
#' @keywords Gene Matrix Venn Diagram
#' @export
#' @examples
#' # Generate a venn diagram from a gene matrix, by cell type.
#' gm.venn(gm, metadata = metadata, grouping_var = c("Type"),
#'        file = 'gm_venn.pdf')
#'
#' # Generate a venn diagram from a gene matrix by cell type, only count genes
#' # that are expressed in 3 or more cells per group and have expression levels
#' # greater than or equal to 1.
#' gm.venn(gm, metadata = metadata, grouping_var = c("Type"), cell_no = 3,
#'         cell_thresh = 1, file = 'gm_venn.pdf')
#'

gm.venn <- function(gm, metadata, grouping_var, cell_no = 1, cell_thresh = 0,
                    raw_or_percent = c("raw", "percent"), file, ext = "pdf", col = NULL){

  # Sort gene matrix by grouping variable(s)
  sort_order = with(metadata, naturalsort::naturalorder(get(grouping_var)))

  gm = gm[sort_order, ]

  metadata = metadata[sort_order, ]

  cols = match(grouping_var, colnames(metadata))

  plot_groups = as.vector(apply(metadata, 1, function(x) {

    paste(x[cols], collapse = ".")

  }))

  if(length(unique(plot_groups)) > 5){stop("Venn Diagram can only handle up to 5 groups")}

  idx = c(0, cumsum(table(plot_groups)))

  out = list()

  for(i in 1:(length(idx) - 1)){
    
    gm_t = gm[(as.integer(idx[i]) + 1):as.integer(idx[i+1]), , drop = FALSE]
    
    gm_t = which((colSums(gm_t > cell_thresh) >= cell_no) == T)

    out[[names(idx)[i+1]]] = gm_t

  }

  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  if(is.null(col) == TRUE){
    
    col_t = RColorBrewer::brewer.pal(n = 5, name = "Set1")
    
  } else {
    
    col_t = col
    
  }
  
  if(missing(file) == F){

    if(tolower(ext) == "pdf"){
      
      pdf(pkgmaker::file_extension(file, ".pdf"), height = 5.69, width = 10.00, onefile = FALSE)
    
    } else if(tolower(ext) == "png"){
      
      png(pkgmaker::file_extension(file, ".png"), height = 5.69, width = 10.00, units = "in", res = 600)
      
    }

    vd = VennDiagram::venn.diagram(out, filename = NULL,
                              print.mode = raw_or_percent,
                              category.names = names(out),
                              cat.default.pos = c("outer"),
                              sigdigs = 2,
                              fill = col_t[1:length(unique(plot_groups))])

    g = grid::gTree(children = grid::gList(vd))

    blank = grid::grid.rect(gp = grid::gpar(col="white"))

    gridExtra::grid.arrange(blank, g, blank, ncol = 3, widths = c(1, 3, 1))

    graphics.off()

  }

  out = venn.gene_lists(out, gm)

  return(out)

}

venn.gene_lists <- function(gm.venn_obj, gm){

  combos <- unlist(lapply(1:length(gm.venn_obj), function(j) {
                        combn(names(gm.venn_obj), j, simplify = FALSE)
                      }), recursive = FALSE)

  names(combos) <- sapply(combos, function(i) paste0(i, collapse = ""))

  gene_lists <- lapply(combos, function(i) {
    set <- venn.setdiff(gm.venn_obj[i], gm.venn_obj[setdiff(names(gm.venn_obj), i)])
    colnames(gm)[set]
  })

}

venn.setdiff <- function (x, y) {
  i <- Reduce('intersect', x)
  u <- Reduce('union', y)
  setdiff(i, u)
}

