#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# library(icell8, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# library(plotly, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# library(shiny, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# library(shinyBS, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# library(shinythemes, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# library(colourpicker, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# library(crosstalk, lib.loc = '/opt/shiny-server/samples/sample-apps/tsne/tsne_libs/')
# pca_tsne_all_loc <<- '/opt/shiny-server/samples/sample-apps/tsne/pca_tsne_all.rda'
# pca_tsne_all_loc <- '/Users/chardy/Desktop/icell8_dev/scripts/tertiary_analysis/tsne_shiny/ICELL8_cluster_data.rda'
library(icell8)
library(plotly)
library(shiny)
library(shinyBS)
library(shinythemes)
library(colourpicker)
library(crosstalk)
library(shinycssloaders)
set.seed(1986)



# A375 Marker MIA (ENSG00000261857)
# KU812 Markers HBB (ENSG00000244734), HBG1 (ENSG00000213934), HBG2 (ENSG00000196565)
# gene_id <- "ENSG00000244734"
#
# sapply(unique(metadata$Type), function(x){
#   t <- rownames(subset(metadata, metadata$Type == x))
#   mean(gm[gene_id,t])
# })
#
# sapply(unique(metadata$Type), function(x){
#   t <- rownames(subset(metadata, metadata$Type == x))
#   max(gm[gene_id,t])
# })


# pca_tsne_all_loc <<- '/Users/chardy/Documents/internal_support/ishminder/sce-1446/gm.cluster_data.rda'
# load(pca_tsne_all_loc)
# tsne_obj <<- pca_tsne$tsne
# metadata <<- pca_tsne$metadata
# gm_top <<- pca_tsne$gms_top_genes
# gm <- pca_tsne$gm
#
# p_default <<- tsne_obj$perplexity

#tsne_shiny(pca_tsne_all_loc)
#pca_tsne_all_loc <<- '/Users/chardy/Desktop/icell8_dev/scripts/tertiary_analysis/example/all/pca_tsne_all.rda'
# load(pca_tsne_all_loc)
# tsne_obj <<- pca_tsne$tsne
# metadata <<- pca_tsne$metadata
# gm_top <<- pca_tsne$gms_top_genes
# gm <- pca_tsne$gms_merged$gms_merged
# gm <- gm[rownames(metadata),]
# gm <<- log(gm + 1, base = exp(1))
# p_default <<- tsne_obj$perplexity



# Define UI for application that draws a histogram
ui <- fluidPage(

  # theme = shinytheme("spacelab"),

  # div(style="",
  #   div(style="padding: 1px 0px; width: '100%'",
  #       titlePanel(
  #         title="", windowTitle=HTML("ICELL8&reg; Analysis")
  #       )
  #   ),
  #   navbarPage(
  #     title=div(style = "width:100%;", HTML("ICELL8<sup>&reg;</sup> Analysis"), div(style = "position:absolute; right:0; top:7px; padding: 0px 10px",
  #       a(href="http://takarabio.com/", target="_blank", img(src="takara_noback.png", width = 100))
  #       )),
  #       navbarMenu("More",
  #              tabPanel("Correlation Analysis", ""),
  #              tabPanel("PCA Analysis", "")
  #   ))
  # ),

  # tags$br(),

    # inputPanel(

      # actionButton("ud_input_data", label = "Upload New Dataset", icon("upload"), width = "100%"),

      fluidRow(
        column(4, uiOutput("grouping_var")),
        column(4, uiOutput("gene_id")),
        column(4, uiOutput("perplexity"))

      ),

      fluidRow(
        column(4, actionButton("ud_grouping_var", label = "Color by group(s)")),
        column(4, actionButton("ud_gene_id", label = "Highlight gene(s)")),
        column(4, actionButton("ud_perp", label = "Update perplexity")),

        bsPopover(id = "q1", title = "Note",
                  content = paste("Perplexity defaults have been optimized for general use cases of the ICELL8 platform. ",
                                  "These values are different from the standard defaults in the Rtsne Package, and may need ",
                                  "to be re-optimized for unique applications. For more information please refer to the following resources: ",
                                  "https://github.com/jkrijthe/Rtsne & ",
                                  "distill.pub/2016/misread-tsne", sep = ""),
                  placement = "right",
                  trigger = "click"
        )
      ),

    # ),

    fluidRow(
    mainPanel(

      div(style = "width: 750px; height: 525px;",
        div(style = "position:absolute; top: 50px; width: 750px; height: 500px;",
          plotlyOutput("Plot1", inline = TRUE)
        ),
      div(style = "position:absolute; top:25px; left:295px;",
          p(style = "color:black; font-family:helvetica; font-size:20px;",tags$em("tSNE Analysis"))
        )
      )

      # fluidRow(style = "padding-bottom: 20px;",
      #
      #          # column(4, align = 'left', downloadButton("downloadData", "Download")),
      #          # column(4, fileInput("input_data", "Choose ICELL8_cluster_data file.")),
      #
      #
      # ))

  )),

  fluidRow(style = "padding-top: 25px;",

    column(4, actionButton("trigger_download", label = "Download", icon("download"))),

    column(4,
        conditionalPanel(
          condition = "input.grouping_var == 'custom'",
          actionButton("clear_custom", "Clear custom selections")
        ), offset = 4)

  ),

  bsModal("startupModal", title = "Interactive tSNE App", trigger = '',
          actionButton("input_data", "Start App", width = "100%"),
          tags$style("#startupModal .modal-footer{ display:none}"),
          size = "small"),

  bsModal("downloadModal", title = "Download Data", trigger = '',
          fluidRow(
            column(4, downloadButton("downloadPlot", label = "Plot", icon("download"))),
            column(4, downloadButton("downloadData", label = "Data", icon("download")))
          ),
          tags$style("#downloadModal .modal-footer{ display:none}"),
          size = "small")


  # fluidRow(style = "padding-bottom: 20px;",
  #
  #   plotOutput('tsne_plot')
  #
  # ),
  #
  # fluidRow(style = "padding-bottom: 20px;",
  #
  #   column(4, downloadButton("downloadData", "Download"))
  #
  # )
  #
  #   ,

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  options(shiny.maxRequestSize=30*1024^2)

  # output$tsne_plot <- renderPlot({
  #   tsne.plot(tsne_obj, metadata = metadata, grouping_var = if(input$ud_perp == 0){c("Type")})
  # })
  #
  # observeEvent(input$ud_perp,{
  #   output$tsne_plot <- renderPlot({
  #     withProgress(message = "Please wait: updating tSNE plot", {
  #       incProgress(1/2, detail = "Calculating tSNE")
  #       tsne_obj <<- gm.tsne(gm_top, metadata = metadata, grouping_var = c("Type"),
  #                            perplexity = as.numeric(input$perplexity), plot = F)
  #     })
  #     tsne.plot(tsne_obj, gm = gm, metadata = metadata, grouping_var = c("Type"))
  #   })
  #   updateSelectizeInput(session, 'gene_id', selected = "")
  #   updateSelectizeInput(session, 'grouping_var', selected = "")
  # })
  #
  # observeEvent(input$ud_grouping_var,{
  #   output$tsne_plot <- renderPlot({
  #     tsne.plot(tsne_obj, gm = gm, metadata = metadata,
  #               grouping_var = if(length(input$grouping_var) == 0){
  #                                 c("Type")
  #                              } else {as.vector(input$grouping_var)}
  #              )
  #   })
  #   updateSelectizeInput(session, 'gene_id', selected = "")
  # })
  #
  # observeEvent(input$ud_gene_id,{
  #   output$tsne_plot <- renderPlot({
  #     tsne.plot(tsne_obj, gm = gm, metadata = metadata, grouping_var = c("Type"),
  #               gene_highlight = as.vector(input$gene_id))
  #   })
  #   updateSelectizeInput(session, 'grouping_var', selected = "")
  # })
  #
  output$downloadPlot <- downloadHandler(
    filename = function() {

      gv <- if(length(input$grouping_var) == 0){
              c("")
            } else {paste(input$grouping_var, collapse = ".")}

      gh <- if(length(input$gene_id) == 0){
              c("")
            } else {paste(as.vector(input$gene_id), collapse = ".")}

      p <- paste("p", input$perplexity, sep = "")

      t <- if(gv == "" & gh == ""){
              p
           } else if(gv == ""){
              paste(gh, p, sep = "_")
           } else if(gv == ""){
              paste(gh, p, sep = "_")
           }

      if(length(t) > 0){
        paste("tsne_", t, ".pdf", sep="")
      } else {
        "tsne.pdf"
      }

    },

    content = function(file) {

      pdf(pkgmaker::file_extension(file, ".pdf"), height = 5.69, width = 10.00)

      add.alpha <- function(col, alpha=1){
        if(missing(col))
          stop("Please provide a vector of colours.")
        apply(sapply(col, col2rgb)/255, 2, function(x) {
          rgb(x[1], x[2], x[3], alpha=alpha)
        })
      }

      tsne.plot(d$tsne_obj, gm = d$gm, metadata = d$metadata,
                grouping_var = if(is.null(input$grouping_var) == FALSE){
                                  as.vector(input$grouping_var)
                               },
                gene_highlight = if(is.null(input$gene_id) == FALSE){
                                    as.vector(input$gene_id)
                                 },
                col_override = if(is.null(input$grouping_var) == FALSE){
                                  if(input$grouping_var == "custom"){
                                    icell8:::add.alpha(d$custom_col)
                                  }})
      graphics.off()

      }
  )


  Dataset <- reactive({
    if (is.null(input$input_data)) {
      # User has not uploaded a file yet
      return("NULL")
    }

    cluster_data <- list(gm = d$gm,
                         gms_top_genes = d$gm_top,
                         metadata = d$metadata,
                         pca = d$pca_obj, tsne = d$tsne_obj,
                         grouping_var = d$grouping_var,
                         transformation_type = d$transformation_type)
    return(cluster_data)

  })

  ud_data <- function(d. = d){

    cluster_data <- list(gm = d.$gm,
                         gms_top_genes = d.$gm_top,
                         metadata = d.$metadata,
                         pca = d.$pca_obj, tsne = d.$tsne_obj,
                         grouping_var = d.$grouping_var,
                         transformation_type = d.$transformation_type)

    cluster_data <- serialize(cluster_data, NULL)

    return(cluster_data)

  }

  output$downloadData <- downloadHandler(
    filename = function() {

      "ICELL8_cluster_data.rda"

    },

    content = function(file) {

      cluster_data <- ud_data()

      save(cluster_data, file = file)

    },

    contentType = '.rda'
  )

  # output$downloadData <- plotly_IMAGE(p, format = "png", out_file = "output.png")

  #### TESTING #####

  print(Sys.time())

  toggleModal(session, "startupModal", toggle = "open")

  inFile <- reactive({
    if(input$input_data == 0){
      return(NULL)
    } else {
      return(out = list(datapath = cluster_data_loc))
    }
    })

  loader <- reactive({
    if(is.null(inFile())){
      return(NULL)
    } else {
      return(inFile()$datapath)
    }
  })

  d <- reactiveValues(gm = NULL, gm_top = NULL, metadata = NULL,
                      pca_obj = NULL, tsne_obj = NULL,
                      grouping_var = NULL,
                      transformation_type = NULL,
                      plot.df = NULL)

  observeEvent(loader(),{
    load(loader())
    withProgress(message = "Please wait:", {
      incProgress(1/4, detail = "Loading Gene Matrix")
      d$gm <- cluster_data$gm
      d$gm_top <- cluster_data$gms_top_genes
      d$pca_obj <- cluster_data$pca
      d$grouping_var <- cluster_data$grouping_var
      d$transformation_type <- cluster_data$transformation_type

      incProgress(1/4, detail = "Loading Metadata")
      d$metadata <- cluster_data$metadata
      d$tsne_obj <- cluster_data$tsne
      p_default <- d$tsne_obj$perplexity
      d$metadata$custom <- rep("undefined", nrow(d$metadata))
      d$custom_col <- c(undefined = "gray")
      d$plot.df <- data.frame(d$tsne_obj$Y[,1], d$tsne_obj$Y[,2], Class = d$metadata$custom)
      colnames(d$plot.df) <- c("x", "y", "Class")
      rownames(d$tsne_obj$Y) <- rownames(d$metadata)

      incProgress(1/4, detail = "Rendering tSNE Plot")
      updateSelectizeInput(session, "grouping_var", choices = colnames(d$metadata))
      updateSelectizeInput(session, "gene_id", choices = rownames(d$gm))
      updateNumericInput(session, "perplexity", value = p_default)
      incProgress(1/4, detail = "Rendering tSNE Plot")
    })
    toggleModal(session, "startupModal", toggle = "close")
  })

  observeEvent(input$ud_input_data, {
    toggleModal(session, "startupModal", toggle = "open")
  })

  observeEvent(input$trigger_download, {
    toggleModal(session, "downloadModal", toggle = "open")
  })

  print("Finished loading user data")

  output$grouping_var <- renderUI({
    selectInput("grouping_var", label = p("Select Group(s)"), choices = NULL, multiple = T, selectize = T)
  })

  output$gene_id <- renderUI({
    selectInput("gene_id", label = p("Select Gene(s)"), choices = NULL, multiple = T, selectize = T)
  })

  output$perplexity <- renderUI({
    numericInput("perplexity", label = p("Perplexity",
                                         tags$style(type = "text/css", "#q1 {vertical-align: top; margin-left: 10px;}"),
                                         bsButton("q1", label = "", icon = icon("info-circle"), style = "info", size = "extra-small")

    ), value = NULL, min = 0, step = .01)
  })

  gm.plotly <- function(plot.df, col_t, cols_t){
    print('actual plot start')

    par(pty = 's')

    t <- list(
      family = "Helvetica",
      color = 'black')

    p <- plot_ly(d$plot.df, x = ~x, y = ~y,
          source = "A", type = "scatter",
          mode = "markers", alpha = 0.5,
          text = if(is.numeric(col_t) == TRUE){
            round(col_t, 2)
          } else {
            col_t
          },
          hoverinfo = 'text',
          color = col_t,
          colors = cols_t,
          height = 500, width = 750,
          marker = list(size = 8)) %>%
    layout(
           # title = "<i>tSNE Analysis</i>",
           titlefont = 16,
           font = t,
           showlegend = TRUE, legend = list(font = list(size = 14)),
           autosize = FALSE,
           margin = list(
             l=125,
             r=225,
             b=75,
             t=25,
             pad=2
           ),
           xaxis = list(title = "tSNE 1",
                        showgrid = FALSE,
                        zeroline = FALSE,
                        hoverformat = '.2f',
                        linecolor = '#7f7f7f',
                        mirror = TRUE,
                        tick = "outside",
                        ticklen = 10,
                        tickfont = list(family = "Helvetica", size = 14),
                        tickcolor = "#7f7f7f"
           ),
           yaxis = list(title = "tSNE 2",
                        showgrid = FALSE,
                        zeroline = FALSE,
                        hoverformat = '.2f',
                        linecolor = '#7f7f7f',
                        mirror = TRUE,
                        tick = "outside",
                        ticklen = 10,
                        tickfont = list(family = "Helvetica", size = 14),
                        tickcolor = "#7f7f7f"
           ),
           dragmode = "lasso",
           d$plot.df) %>% config(displaylogo = FALSE,
                               collaborate = FALSE,
                               modeBarButtonsToRemove = list(
                                 'sendDataToCloud',
                                 'toImage',
                                 'autoScale2d',
                                 'resetScale2d',
                                 'hoverClosestCartesian',
                                 'hoverCompareCartesian',
                                 'zoom2d',
                                 'pan2d','select2d','lasso2d','zoomIn2d','zoomOut2d','toggleSpikelines')
           )

    return(p)

  }

  observeEvent(input$ud_perp,{
    print("updating tsne obj")
    withProgress(message = "Please wait: updating tSNE plot", {
      incProgress(1/2, detail = "Calculating tSNE")
      d$tsne_obj <- gm.tsne(d$gm_top, metadata = d$metadata, grouping_var = c("Type"),
                           perplexity = as.numeric(input$perplexity), plot = F)
      })
      d$plot.df <- data.frame(d$tsne_obj$Y[,1], d$tsne_obj$Y[,2], Class = d$metadata$Type)
      colnames(d$plot.df) <- c("x", "y", "Class")
      rownames(d$tsne_obj$Y) <- rownames(d$metadata)
  }, ignoreInit = TRUE)

  add_selection <- function(v, metadata, custom_col, tsne_obj, label, new_col){
    print('adding selection')
    sel <- apply(v[, c(3:4)], 1, paste, collapse = ".")
    ref <- apply(tsne_obj$Y, 1, paste, collapse = ".")
    i <- which(ref %in% sel)
    metadata[["custom"]][i] <- label
    custom_col[[label]] <- new_col
    print(custom_col)
    return(list(metadata = metadata, custom_col = custom_col))
  }

  observeEvent(event_data("plotly_selected", source = "A"), {
    print('data selected')
    col_r <- c(rainbow(10, alpha = 0.5))
    col_r <- col_r[c(1,3,5,7,9,2,4,6,8,10)]
    s <- length(unique(d$metadata$custom))
    print(s)
    d <- isolate(event_data("plotly_selected", source = "A"))
    if(length(d) > 0){
      showModal(modalDialog(
        size = "s",
        textInput("new_group_label", "Enter Label for Selected Group."),
        colourInput("col", "Select colour", col_r[s],
                    showColour = c("background")),
        easyClose = TRUE,
        footer = list(actionButton("set_group", "Set Custom Group"),
                      modalButton("Dismiss"))
      ))
    }
  }, ignoreInit = TRUE)

  gene_color <- function(gm, metadata, gene_id){
    print("Gene Colors")
    sort_order = rownames(metadata)
    cols = match(gene_id, rownames(gm))
    col_val = apply(gm[, sort_order], 2, function(x) {mean(x[cols])})
    return(col_val)
  }

  observeEvent(input$set_group, {
    print('set group activated')
    removeModal()
    res <- add_selection(v = event_data("plotly_selected"), d$metadata, d$custom_col,
                               d$tsne_obj, label = input$new_group_label,
                               new_col = input$col)
    d$metadata <- res$metadata
    d$custom_col <- res$custom_col
    updateSelectizeInput(session, 'grouping_var', selected = "custom")
    updateSelectizeInput(session, 'gene_id', selected = "")
  }, ignoreInit = TRUE)

  observeEvent(input$clear_custom, {
    print('clear custom activated')
    d$metadata$custom <- rep("undefined", nrow(d$metadata))
    d$custom_col <- c(undefined = "gray")
    updateSelectizeInput(session, 'grouping_var', selected = "")
  }, ignoreInit = TRUE)

  observeEvent(input$ud_gene_id, {
    print('ud gene id activated')
    if(is.null(isolate(input$grouping_var)) == FALSE){
      updateSelectizeInput(session, 'grouping_var', selected = "")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$ud_grouping_var, {
    print('ud grouping var activated')
    if(is.null(isolate(input$gene_id)) == FALSE){
      updateSelectizeInput(session, 'gene_id', selected = "")
    }
  }, ignoreInit = TRUE)

  plot_style <- eventReactive(list(input$grouping_var,
                                   input$gene_id,
                                   input$ud_grouping_var,
                                   input$ud_gene_id,
                                   input$set_group,
                                   input$clear_custom,
                                   input$ud_perp,
                                   input$perplexity
                                     ),{

    t1 <- input$grouping_var
    t2 <- input$gene_id
    t3 <- input$ud_grouping_var
    t4 <- input$ud_gene_id
    t5 <- input$set_group
    t6 <- input$clear_custom
    t7 <- input$ud_perp
    t8 <- input$perplexity

    col_t <- rep("undefined", nrow(d$metadata))
    cols_t <- colorRampPalette(c('gray'), alpha = 0.5)(1)
    print('plot style activated')
    if(is.null(t1) == FALSE){
      if(t1 == "custom"){
        print("set group")
        col_t <- d$metadata[["custom"]]
        cols_t <- d$custom_col
      }
    }
    if(t3){
      print("t3 activated")
      if(is.null(t1) == FALSE & is.null(t2) == TRUE){
        print("Plot grouping Var!")
        if(length(t1) > 1){
          cols <- match(t1, colnames(d$metadata))
          col_t <- do.call(paste, c(d$metadata[cols], list(sep = ".")))
          cols_t <- rainbow(length(unique(col_t)), alpha = 0.5)
        } else if(t1 == "custom"){
          col_t <- d$metadata[["custom"]]
          cols_t <- d$custom_col
        } else {
          col_t <- d$metadata[[t1]]
          if(is.numeric(col_t) == TRUE){
            cols_t <- colorRampPalette(c('blue', 'cyan', 'yellow', 'red'), alpha = 0.5)(100)
          } else {
            cols_t <- rainbow(length(unique(col_t)), alpha = 0.5)
          }
        }
      }
    }
    if(t4){
      print("t4 activated")
      if(is.null(t1) == TRUE & is.null(t2) == FALSE){
        print("Plot gene id!")
        col_t <- gene_color(d$gm, d$metadata, isolate(t2))
        cols_t <- colorRampPalette(c('antiquewhite', 'firebrick3'), alpha = 0.5)(ceiling(max(col_t) + 0.1))
      }
    }
    return(list(col_t = col_t, cols_t = cols_t))
  }, ignoreInit = TRUE)

  # output$Plot1 <- renderPlotly({
  #   if(is.null(input$input_data) == TRUE){
  #     return(NULL)
  #   } else {
  #     t <- plot_style()
  #     gm.plotly(d$plot.df, col_t = t$col_t, cols_t = t$cols_t)
  #   }
  # })

  thePlot <- reactive({
    if(input$input_data == 0){
      return(NULL)
    } else {
      t <- plot_style()
      p <- gm.plotly(d$plot.df, col_t = t$col_t, cols_t = t$cols_t)
    }
  })

  output$Plot1 <- renderPlotly({
    thePlot()
  })

  session$onSessionEnded(stopApp)

}


# Run the application
shinyApp(ui, server)




