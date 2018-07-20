#' ICELL8 Analysis
#'
#' Perform quality control, normalization, transformation & filtering, correlation and clustering analysis on single cell RNAseq data from ICELL8.  Produce a standardized report.
#' @param ICELL8_data_loc Full path to ICELL8_data folder, obtained from primary pipeline (ie. /path/to/ICELL8_data).
#' @param gm_file_list Plain text file with a list of gene matrices for analysis.  Each line of the file must include the full path to the gene matrix (ie. /path/to/genematrix.csv).
#' @param gm_report_list Plain text file with a list of report files for analysis.  Each line of the must include the full path to the report file for each gene matrix in the --gm_file_list option (ie. /path/to/genematrix_report.csv).
#' @param cell_type Optional, comma separated list of cell types/groups of interest. (ie. –t “Positive_Control, K562”).  Defaults to all cell types found in the 'Type' Column of the report files under the --gm_report_list option.
#' @param correlation_type Correlation method type, pick from pearson, spearman or kendall.  Defaults to ‘pearson’.
#' @param parallel_cores The number of cores to run in parallel.  If running large datasets on a server, increasing the number of parallel cores can dramatically increase the speed of the analysis.
#' @param names Optional, comma separated list of names for each gene matrix in the file supplied in the --gm_file_list option (ie. –-names “group1, group2, group3”).  If not specified the basename of each gene matrix file will be used as the name.  For example, /analyses/raw_data/human/A375.100K.withUMI.csv would be shortened to A375.100K.withUMI.
#' @param cor_low_thresh For correlation analysis, filter genes from each cell with expression less than this value.
#' @param cor_high_thresh For correlation analysis, filter genes from each cell with expression higher than this value.
#' @param quant_or_abs Choose 'quant' or 'abs' to define whether the values for the --lo_thresh and --hi_thresh options are quantile or absolute cutoffs.
#' @param gm_norm TRUE/FALSE, should gene matrix be normalized?  Defaults to FALSE.  If TRUE, defaults to transcripts per million unless the --gm_norm_method option is also specified.
#' @param gm_norm_method What normalization method to use? Normalize to 'tps' (Transcript Per Scale), 'cc_mc' (cell coverage by median coverage) or 'rpkm' (reads per kilobase million). If 'tps' is chosen, the --gm_norm_scale option must be set.  To normalize each cell to Transcripts Per Million set --gm_norm_method to 'tps' and --gm_norm_scale to 1e6.
#' @param gm_norm_scale What scaling factor should be used for normalization?  Defaults to 1e6 for Transcripts per Million.
#' @param gm_norm_genelengths Provide a table of gene lengths used for rpkm normalization.  Defaults to NULL.
#' @param gm_log TRUE/FALSE, should gene matrix be log transformed?  If TRUE, defaults to natural log 'ln', unless the -–gm_log_base option is specified.
#' @param gm_log_base What log base should be used for log transformation?  For natural log select ‘ln’.  Defaults to 'ln'.
#' @param qc_cell_abslowcov Define minimum read depth per cell. Defaults to 10,000 reads.
#' @param qc_cells_abslowgenecount Define minimum gene counts per cell. Defaults to 3 genes.
#' @param qc_gene_cellcount The number of cells a gene must have > 0 expression to be kept in the gene matrix.  Defaults to 1.
#' @param qc_gene_totcov The minimum number of reads for a gene across all cells. Defaults to 5.
#' @param grouping_var Comma separated list of factors from report files to create groups for analysis. Defaults to ‘Group’, which classifies groups by each gene matrix that is entered.  For example, if the analysis should observe differences between cell types and groups enter \"Group, Type\".
#' @param pca_filt_method Method to choose which genes will be used for clustering analysis.  Choose from 'quant_exp' (Filter genes below a quantile value specified in the --pca_thresh_cut option), 'top_var' (Select the genes with the highest variance- must also set the --pca_top_genes option), or 'top_exp' (Select the highest expressing genes- must also set the --pca_top_genes option).  If the option is not specified, all genes are used in the analysis.
#' @param pca_thresh_cut For clustering analysis, filter genes with expression less than this quantile value (values must be between 0 & 1). Defaults to 0.0 (No genes filtered).
#' @param pca_top_genes Select the number of genes to be used in the clustering analysis?  Defaults to 500.
#' @param merge_method How should multiple gene matrices be merged?  Choose from ‘intersect’ or ‘union’.  Defaults to ‘union’, where missing genes are filled in with zero read counts.
#' @param correlation_analysis TRUE/FALSE, perform correlation analysis. Defaults to TRUE.
#' @param cluster_analysis TRUE/FALSE, perform clustering analysis (PCA & tSNE). Defaults to TRUE.
#' @param report TRUE/FALSE, generate final report.
#' @param ICELL8_data TRUE/FALSE, is data derived from ICELL8 primary pipeline?
#' @param figure_type Enter 'pdf' or 'png' for figure outputs.  Defaults to 'pdf'
#' @param transpose TRUE/FALSE, transpose gene matrix. Should be columns of samples, rows of genes. Defaults to FALSE.
#' @param verbose TRUE/FALSE, print status to console?  Defaults to TRUE.
#' @param output_dir Name of the output directory.  Defaults to ICELL8_output in current working directory.  The script will not overwrite existing folders, with the exception of the ICELL8_output folder.
#' @param title Title of experiment.
#' @param author Name of indvidual(s) performing analysis.
#' @keywords tertiary analysis single-cell single cell rna seq rnaseq rna-seq clustering pca tsne correlation ICELL8 takara
#' @export

ICELL8.analysis <- function(ICELL8_data_loc = NULL,
                            gm_file_list = NULL,
                            gm_report_list = NULL,
                            cell_type = "all",
                            correlation_type = "pearson",
                            parallel_cores = 1,
                            names = NULL,
                            cor_low_thresh = 0,
                            cor_high_thresh = 1,
                            quant_or_abs = "quant",
                            gm_norm = FALSE,
                            gm_norm_method = "tps",
                            gm_norm_scale = 1e6,
                            gm_norm_genelengths = NULL,
                            gm_log = FALSE,
                            gm_log_base = "ln",
                            qc_cell_abslowcov = 1e4,
                            qc_cells_abslowgenecount = 3,
                            qc_gene_cellcount = 1,
                            qc_gene_totcov = 5,
                            grouping_var = NULL,
                            pca_filt_method = "top_var",
                            pca_thresh_cut = 0,
                            pca_top_genes = 500,
                            merge_method = "union",
                            correlation_analysis = TRUE,
                            cluster_analysis = TRUE,
                            report = TRUE,
                            ICELL8_data = TRUE,
                            figure_type = "png",
                            transpose = FALSE,
                            verbose = FALSE,
                            output_dir = "ICELL8_output",
                            title = "ICELL8 Analysis Report",
                            author = "")
{

  # ---------- start script ---------- #

  start_time <- Sys.time()


  # ---------- parse user options ---------- #

  opt <- list(
    ICELL8_data_loc <- ICELL8_data_loc,
    gm_loc <- gm_file_list,
    md_loc <- gm_report_list,
    gm_names <- hanta:::setup.scrub(v = names),
    grouping_var <- hanta:::setup.scrub(v = grouping_var),
    cell_types <- hanta:::setup.scrub(v = cell_type),
    cor_type <- tolower(correlation_type),
    transpose <- transpose,
    core_count <- parallel_cores,
    lo <- cor_low_thresh,
    hi <- cor_high_thresh,
    quant_or_abs <- tolower(quant_or_abs),
    gm_norm <- gm_norm,
    gm_norm_method <- tolower(gm_norm_method),
    gm_norm_scale <- gm_norm_scale,
    gm_norm_genelengths <- gm_norm_genelengths,
    gm_log <- gm_log,
    gm_log_base <- gm_log_base,
    cell_abslowcov <- qc_cell_abslowcov,
    cells_abslowgenecount <- qc_cells_abslowgenecount,
    gene_cellcount <- qc_gene_cellcount,
    gene_totcov <- qc_gene_totcov,
    pca_filt_method <- tolower(pca_filt_method),
    pca_thresh_cut <- pca_thresh_cut,
    pca_top_genes <- pca_top_genes,
    perform_correlation_analysis <- correlation_analysis,
    perform_cluster_analysis <- cluster_analysis,
    report <- report,
    ICELL8_data <- ICELL8_data,
    figure_type <- tolower(figure_type),
    merge_method <- tolower(merge_method),
    verbose <- verbose,
    output_dir <- file.path(output_dir),
    title <- title,
    author <- author
  )

  names(opt) <- c("ICELL8_data_loc",
                  "gm_file_list",
                  "gm_report_list",
                  "names",
                  "grouping_var",
                  "cell_type",
                  "correlation_type",
                  "transpose",
                  "parallel_cores",
                  "cor_low_thresh",
                  "cor_high_thresh",
                  "quant_or_abs",
                  "gm_norm",
                  "gm_norm_method",
                  "gm_norm_scale",
                  "gm_norm_genelengths",
                  "gm_log",
                  "gm_log_base",
                  "qc_cell_abslowcov",
                  "qc_cells_abslowgenecount",
                  "qc_gene_cellcount",
                  "qc_gene_totcov",
                  "pca_filt_method",
                  "pca_thresh_cut",
                  "pca_top_genes",
                  "correlation_analysis",
                  "cluster_analysis",
                  "report",
                  "ICELL8_data",
                  "figure_type",
                  "merge_method",
                  "verbose",
                  "output_dir",
                  "title",
                  "author")


  # check gene matrix & metadata (report) list files are specified / exists & return locations
  res = hanta:::setup.fl_read(ICELL8_data_loc = ICELL8_data_loc,
                               gm_file_list = gm_loc, md_file_list = md_loc)
  gms = res$gms
  mds = res$mds


  # check that grouping variables exist in metadata (report) files
  hanta:::setup.grouping_var_check(mds = mds, grouping_var = grouping_var,
                                    ICELL8 = ICELL8_data)


  # check that cell types exist in metadata (report) files
  cell_types = hanta:::setup.cell_type_check(mds = mds, cell_types = cell_types,
                                              ICELL8 = ICELL8_data)


  # check that correlation type matches an acceptable option
  hanta:::setup.cor_check(cor_type = cor_type)


  # check / set number of avaialble cores for parallel processing
  hanta:::setup.core_count(core_count = core_count)


  # check threshold values for cutoffs
  hanta:::setup.thresh_check(lo = lo, hi = hi, quant_or_abs = quant_or_abs)


  # check qm.qc options / parameters
  # normalization
  opt = hanta:::setup.norm_check(gm_norm = gm_norm,
                                  gm_norm_method = gm_norm_method,
                                  gm_norm_scale = gm_norm_scale,
                                  gm_norm_genelengths = gm_norm_genelengths,
                                  opt = opt)

  genelengths = hanta:::setup.genelengths(gm_norm_genelengths)

  # transformation
  opt = hanta:::setup.tform_check(gm_log = gm_log, gm_log_base = gm_log_base,
                                   opt = opt)


  # check gm.pca options
  res = hanta:::setup.pca_check(pca_filt_method = pca_filt_method,
                                 pca_thresh_cut = pca_thresh_cut,
                                 pca_top_genes = pca_top_genes, opt = opt)
  pca_filt_method = res$pca_filt_method
  opt = res$opt


  # check report style
  hanta:::setup.report_check(report = report)


  # check & set figure_type
  ext = hanta:::setup.figure_type(figure_type)


  # check merge method
  hanta:::setup.merge_method_check(merge_method = merge_method)


  # check and create output directory
  output_dir = hanta:::setup.output_dir(output_dir = output_dir)

  # create log file, save option information
  log_file = file.path(output_dir, paste("log", basename(output_dir), sep = "_"))

  hanta:::log.options_summary(log_file = log_file, start_time = start_time,
                               gms = gms, mds = mds, opt = opt)

  # cleanup objects
  rm(list = c("opt", "res"))

  # ---------- main ---------- #


  # ---------- load & merge genematrix and metadata files ---------- #

  hanta:::log.out(paste0("Setting up gene matrix and metadata file(s) for ICELL8 ",
                          "tertiary analysis."), file = log_file, append = T, echo = verbose)

  raw_data <- gm.setup(gms_loc = gms, metadata_loc = mds, gm_names = gm_names,
                       merge_type = merge_method, transpose = transpose,
                       ICELL8 = ICELL8_data)

  save(raw_data, file = file.path(output_dir, "ICELL8_raw_data.rda"))


  # ---------- qc, normalziation & transformation ----------#

  hanta:::log.out("Performing QC Procedures.", file = log_file, append = T,
                   echo = verbose)

  qc_data <- gm.qc(gm = raw_data$gm, metadata = raw_data$metadata, gm_loc = NULL,
                   gm_qcfilt = TRUE, gm_norm = gm_norm,
                   gm_norm_method = gm_norm_method, gm_norm_scale = gm_norm_scale,
                   gm_norm_genelengths = genelengths,
                   gm_log = gm_log, gm_log_base = gm_log_base, mads = 3,
                   cell_abslowcov = cell_abslowcov,
                   cells_abslowgenecount = cells_abslowgenecount,
                   gene_cellcount = gene_cellcount, gene_totcov = gene_totcov)

  save(qc_data, file = file.path(output_dir, "ICELL8_qc_data.rda"))

  rm(raw_data)


  # ---------- subset by cell type parameter ---------- #

  res <- hanta:::qc.cell_type_filter(gm = qc_data$gm,
                                      metadata = qc_data$metadata,
                                      cell_types = cell_types)

  rm(qc_data)


  # ---------- correlation analysis ---------- #


  if(perform_correlation_analysis == TRUE){

    hanta:::log.out("Performing Correlation Analysis", file = log_file, append = T,
                     echo = verbose)

    correlation_data <- gm.cor_analysis(gm = res$gm, metadata = res$metadata,
                                        cor_type = cor_type, quant_or_abs = quant_or_abs,
                                        lo = lo, hi = hi, grouping_var = grouping_var,
                                        heatmap = TRUE, dist_boxplot = TRUE,
                                        dist_table = TRUE, heatmap_cluster = FALSE,
                                        ext = ext, output_dir = output_dir)


    save(correlation_data, file = file.path(output_dir, "ICELL8_correlation_data.rda"))

    rm(correlation_data)

  }


  # ---------- cluster analysis ---------- #

  if(perform_cluster_analysis == TRUE){


    # ---------- pca analysis ----------#

    hanta:::log.out("Performing PCA Analysis", file = log_file, append = T,
                     echo = verbose)

    pca_data <- gm.pca(gm = res$gm, metadata = res$metadata,
                       grouping_var = grouping_var, filt_method = pca_filt_method,
                       thresh_cut = pca_thresh_cut, top_genes = pca_top_genes,
                       transform_type = if(gm_log == TRUE){gm_log_base} else {"none"},
                       plot = TRUE, file = file.path(output_dir, "pca.pdf"), ext = ext)


    # ---------- tsne analysis ----------#

    hanta:::log.out("Performing tSNE Analysis", n = 2, file = log_file, append = T,
                     echo = verbose)

    tsne_data <- gm.tsne(gm = pca_data$pca_genes, metadata = res$metadata,
                         grouping_var = grouping_var, dims = 2, theta = 0.25,
                         max_iter = 1500, pca_center = FALSE, pca_scale = FALSE,
                         initial_dims = hanta:::tsne.intial_dims(pca_data$pca_obj),
                         plot = TRUE, file = file.path(output_dir, "tSNE.pdf"), ext = ext)

    cluster_data <- list(gm = res$gm,
                         gms_top_genes = pca_data$pca_genes,
                         metadata = res$metadata,
                         pca = pca_data$pca_obj, tsne = tsne_data,
                         grouping_var = grouping_var,
                         transformation_type = if(gm_log == TRUE){gm_log_base} else {"none"})

    save(cluster_data, file = file.path(output_dir, "ICELL8_cluster_data.rda"))

    rm(list = c("res", "pca_data", "tsne_data", "cluster_data"))

  }


  # ---------- write out time statistics ---------- #

  hanta:::log.end_time(log_file = log_file, start_time = start_time, end_time = Sys.time())


  # ---------- end ---------- #

  # if(report == "interactive"){
  #
  #   rmarkdown::run(system.file(file.path("reports", "ICELL8_interactive_report.Rmd"), package = "hanta"),
  #                  shiny_args = list(launch.browser = TRUE),
  #                  render_args = list(params = list(set_title = title,
  #                                                   set_author = author,
  #                                                   data_dir = output_dir)),
  #                  quiet = !(verbose)
  #   )
  #
  # }

  if(report == TRUE){

    rmarkdown::render(system.file(file.path("reports", "ICELL8_static_report.Rmd"), package = "hanta"),
                      params = list(set_title = title,
                                    set_author = author,
                                    data_dir = output_dir),
                      output_dir = output_dir,
                      output_file = "ICELL8_report.html",
                      quiet = !(verbose)
    )

  }

}












