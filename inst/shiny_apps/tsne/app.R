#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(icell8)
library(plotly)
library(shiny)
library(shinyBS)
library(shinythemes)
library(colourpicker)
library(crosstalk)
library(shinycssloaders)
set.seed(1986)


is.empty <- function(x){
  if(is.null(x) == TRUE){return(TRUE)}
  else if(is.na(x) == TRUE){return(TRUE)}
  else if(x == ""){return(TRUE)}
  else {return(FALSE)}
}


# ---------- ui function ---------- #

ui <- fluidPage(

  theme = shinytheme("spacelab"),

  div(style="",
    div(style="padding: 1px 0px; width: '100%'",
        titlePanel(
          title="", windowTitle=HTML("ICELL8&reg; Analysis")
        )
    ),
    navbarPage(
      title=div(style = "width:100%;", HTML("ICELL8<sup>&reg;</sup> Analysis"), div(style = "position:absolute; right:0; top:7px; padding: 0px 10px",
        a(href="http://takarabio.com/", target="_blank", img(src="takara_noback.png", width = 100))
        )),
        navbarMenu("More",
               tabPanel("Correlation Analysis", ""),
               tabPanel("PCA Analysis", "")
    ))
  ),

  tags$br(),

  sidebarLayout(
    sidebarPanel(

      actionButton("ud_input_data", label = "Upload New Dataset", icon("upload"), width = "100%"),

      tags$br(),

      tags$br(),

      tags$br(),

      uiOutput("grouping_var"),

      actionButton("ud_grouping_var", label = "Color by group(s)", width = '100%'),

      tags$br(),

      tags$br(),

      tags$br(),

      uiOutput("gene_id"),

      actionButton("ud_gene_id", label = "Highlight gene(s)", width = '100%'),

      tags$br(),

      tags$br(),

      tags$br(),

      uiOutput("perplexity"),

      actionButton("ud_perp", label = "Update perplexity", width = '100%'),

      bsPopover(id = "q1", title = "Note",
                content = paste("Perplexity defaults have been optimized for general use cases of the ICELL8 platform. ",
                                "These values are different from the standard defaults in the Rtsne Package, and may need ",
                                "to be re-optimized for unique applications. For more information please refer to the following resources: ",
                                "https://github.com/jkrijthe/Rtsne & ",
                                "distill.pub/2016/misread-tsne", sep = ""),
                placement = "right",
                trigger = "click"
      ),

      tags$br(),

      tags$br(),

      tags$br(),

      tags$br(),

      tags$br(),

      actionButton("trigger_download", label = "Download", icon("download"), width = "100%"),

      width = 3

    ),

    mainPanel(

      tags$br(),

      tags$br(),

      div(style = "width: 750px; height: 625px;",
        div(style = "position:absolute; top: 75px; width: 750px; height: 500px;",
          plotlyOutput("Plot1", inline = TRUE)
        ),
        div(style = "position:absolute; left:295px;",
            p(style = "color:black; font-family:helvetica; font-size:20px;",tags$em("tSNE Analysis"))
        ),
        div(style = "position:absolute; bottom:50px; left: 500px",
          conditionalPanel(
            condition = "input.grouping_var == 'custom'",
            actionButton("clear_custom", "Clear custom selections")
          )
        )
      )
  )),

  bsModal("startupModal", title = "Select Input Data", trigger = '',
          fileInput("input_data", "Choose ICELL8_cluster_data file."),
          tags$style("#startupModal .modal-footer{ display:none}"),
          size = "small"),

  bsModal("downloadModal", title = "Download Data", trigger = '',
          fluidRow(downloadButton("downloadPlot", label = "Download Plot"),
          downloadButton("downloadData", label = "Download Data")),
          tags$style("#downloadModal .modal-footer{ display:none}"),
          size = "small")

)



# ---------- server function ---------- #

server <- function(input, output, session) {

  options(shiny.maxRequestSize=1000*1024^2)


  # ---------- start app ---------- #

  print(Sys.time())


  # ---------- prompt user to load data | load data ---------- #

  toggleModal(session, "startupModal", toggle = "open")

  inFile <- reactive({
    if (is.null(input$input_data)) {
      return(NULL)
    } else {
      input$input_data
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
      d$cluster_data <- cluster_data
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
      updateNumericInput(session, "perplexity", value = p_default)
      incProgress(1/8, detail = "Loading Gene IDs")
      updateSelectizeInput(session, "gene_id", choices = rownames(d$gm))
    })
    toggleModal(session, "startupModal", toggle = "close")
  })
  
  observeEvent(input$ud_input_data, {
    toggleModal(session, "startupModal", toggle = "open")
  })

  observeEvent(input$trigger_download, {
    toggleModal(session, "downloadModal", toggle = "open")
  })

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


  # ---------- plotly fxn ---------- #

  gm.plotly <- function(plot.df, col_t, cols_t){

    print('actual plot start')

    par(pty = 's')

    t <- list(
      family = "Helvetica",
      color = 'black')

    p <- plotly::plot_ly(d$plot.df, x = ~x, y = ~y, source = "A", type = "scatter",
          mode = "markers", alpha = 0.5, hoverinfo = 'text', color = col_t,
          text = ifelse(is.numeric(col_t) == TRUE, round(col_t, 2), col_t),
          colors = cols_t, height = 500, width = 750, marker = list(size = 8)
          ) %>%
      layout(titlefont = 16, font = t, showlegend = TRUE, autosize = FALSE,
             legend = list(font = list(size = 14)),
             margin = list(l = 125, r = 225, b = 75, t = 25, pad = 2),
             xaxis = list(title = "tSNE 1", showgrid = FALSE, zeroline = FALSE,
                          hoverformat = '.2f', linecolor = '#7f7f7f',
                          mirror = TRUE, tick = "outside", ticklen = 10,
                          tickfont = list(family = "Helvetica", size = 14),
                          tickcolor = "#7f7f7f"),
             yaxis = list(title = "tSNE 2", showgrid = FALSE, zeroline = FALSE,
                          hoverformat = '.2f', linecolor = '#7f7f7f',
                          mirror = TRUE, tick = "outside", ticklen = 10,
                          tickfont = list(family = "Helvetica", size = 14),
                          tickcolor = "#7f7f7f"),
             dragmode = "lasso", d$plot.df
             ) %>%
      config(displaylogo = FALSE, collaborate = FALSE,
             modeBarButtonsToRemove = list('sendDataToCloud', 'toImage',
                                           'autoScale2d', 'resetScale2d',
                                           'hoverClosestCartesian', 'pan2d',
                                           'hoverCompareCartesian', 'zoom2d',
                                           'select2d','lasso2d','zoomIn2d',
                                           'zoomOut2d','toggleSpikelines'))

    return(p)

  }


  # ---------- perplexity update ---------- #

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


  # ---------- handle custom (lasso) selection ---------- #

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


  # ---------- trigger plot ---------- #

  observeEvent(input$ud_grouping_var, {
    print('trigger')
    updateSelectizeInput(session, 'gene_id', selected = "")
    if(length(input$grouping_var) > 1){
      cols <- match(input$grouping_var, colnames(d$metadata))
      col_t <- do.call(paste, c(d$metadata[cols], list(sep = ".")))
      cols_t <- rainbow(length(unique(col_t)), alpha = 0.5)
    } else if(input$grouping_var == "custom"){
      col_t <- d$metadata[["custom"]]
      cols_t <- d$custom_col
    } else {
      col_t <- d$metadata[[input$grouping_var]]
      if(is.numeric(col_t) == TRUE){
        cols_t <- colorRampPalette(c('blue', 'cyan', 'yellow', 'red'), alpha = 0.5)(100)
      } else {
        cols_t <- rainbow(length(unique(col_t)), alpha = 0.5)
      }
    }
    output$Plot1 <- renderPlotly({
      if(is.empty(isolate(input$input_data)[[1]]) == TRUE){
        plotly::plotly_empty(type = "scatter") %>% plotly::config(staticPlot = TRUE)
      } else {
        gm.plotly(d$plot.df, col_t = col_t, cols_t = cols_t)
      }
    })
  }, ignoreInit = TRUE)

  observeEvent(input$ud_gene_id, {
    print('trigger')
    updateSelectizeInput(session, 'grouping_var', selected = "")
    col_t <- gene_color(d$gm, d$metadata, isolate(input$gene_id))
    cols_t <- colorRampPalette(c('gray', 'firebrick3'), alpha = 0.5)(ceiling(max(col_t) + 0.1))
    output$Plot1 <- renderPlotly({
      if(is.empty(isolate(input$input_data)[[1]]) == TRUE){
        plotly::plotly_empty(type = "scatter") %>% plotly::config(staticPlot = TRUE)
      } else {
        gm.plotly(d$plot.df, col_t = col_t, cols_t = cols_t)
      }
    })
  }, ignoreInit = TRUE)

  observe({
    print("no data")
    if(is.empty(input$input_data[[1]]) == TRUE){
      output$Plot1 <- renderPlotly({
        plotly::plotly_empty(type = "scatter") %>% plotly::config(staticPlot = TRUE)
      })
    }
  })

  observe({
    print("no input | default")
    if(is.empty(input$input_data[[1]]) == FALSE){
      if(is.empty(input$grouping_var[[1]]) == TRUE & is.empty(input$gene_id[[1]]) == TRUE){
        output$Plot1 <- renderPlotly({
          col_t <- rep("undefined", nrow(d$metadata))
          cols_t <- colorRampPalette(c('gray'), alpha = 0.5)(1)
          gm.plotly(isolate(d$plot.df), col_t = col_t, cols_t = cols_t)
        })
      }
      if(is.empty(input$grouping_var[[1]]) == FALSE & is.empty(input$gene_id[[1]]) == TRUE){
        if(input$grouping_var == "custom"){
          output$Plot1 <- renderPlotly({
            col_t <- d$metadata[["custom"]]
            cols_t <- d$custom_col
            gm.plotly(isolate(d$plot.df), col_t = col_t, cols_t = cols_t)
          })
        }
      }
    }
  })


  # ---------- download data ---------- #

  dset <- reactive({
    if (is.null(input$input_data)) {
      return(NULL)
    } else {
      d$cluster_data <- list(gm = d$gm,
                             gms_top_genes = d$gm_top,
                             metadata = d$metadata,
                             pca = d$pca_obj, tsne = d$tsne_obj,
                             grouping_var = d$grouping_var,
                             transformation_type = d$transformation_type)
      return(d$cluster_data)
    }
  })

  output$downloadData <- downloadHandler(

    filename = function() {

      "ICELL8_cluster_data.rda"

    },

    content = function(file) {

      cluster_data <- dset()

      save(cluster_data, file = file)

    }

  )


  # ---------- download plot ---------- #

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
        paste("tsne_", t, ".png", sep="")
      } else {
        "tsne.png"
      }

    },

    content = function(file) {

      png(pkgmaker::file_extension(file, ".png"), height = 5.69, width = 10.00, units = "in", res = 600)
      
      tsne.plot(d$tsne_obj, gm = d$gm, metadata = d$metadata,
                grouping_var = if(is.empty(input$grouping_var)[1] == TRUE){
                  NULL
                } else {
                  as.vector(input$grouping_var)
                },
                gene_highlight = if(is.empty(input$gene_id)[1] == TRUE){
                  NULL
                } else {
                  as.vector(input$gene_id)
                },
                col_override = if(is.empty(input$grouping_var)[1] == FALSE){
                  if(input$grouping_var == "custom"){
                    icell8:::add.alpha(d$custom_col)
                  }
                }
      )
      
      graphics.off()

    }
    
  )


  # ---------- end session ---------- #

  session$onSessionEnded(stopApp)

}


# Run the application
shinyApp(ui, server)

