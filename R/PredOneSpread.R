#' Launch an interactive Shiny app to explore predicted spread
#'
#' This function launches a Shiny application that allows users to interactively explore the predicted spread
#' of unmixed signals based on a given userm object (\code{Userm}) and a selected population.
#' The app provides multiple tabs for visualizing the unmixing matrix, prediction plots, slope and intercept matrices,
#' and exporting results.
#'
#' @param Userm A userm object created by \code{\link{CreateUserm}}, containing unmixing matrix, detector and fluor information, and intensity data.
#' @param population_id A character string specifying the column name in \code{Userm$Intensity_mtx} representing the population to be visualized.
#'
#' @return A Shiny UI object that can be passed to \code{shinyApp()} or used within a Shiny server function.
#'
#' @details
#' The app includes the following features:
#' \itemize{
#'   \item Fluor and detector selection panels
#'   \item Interactive prediction plot with scale and axis controls
#'   \item Visualization of unmixing, pseudoinverse, slope, and intercept matrices
#'   \item Export options for HTML reports and NxN/Nx1 plots
#' }
#'
#' @examples
#' \dontrun{
#' PredOneSpread(Userm, population_id = "P1")
#' }
#' @export
#' @import shiny
#' @importFrom MASS ginv
#' @importFrom stats quantile cov dnorm
#' @importFrom rmarkdown render
#' @importFrom lsa cosine

PredOneSpread = function(Userm,population_id){

  if (!population_id %in% colnames(Userm$Intensity_mtx)) {
    stop(paste0("The ",population_id," column is missing in Userm$Intensity_mtx."))
  }

  ui <- fluidPage(
    tags$head(

      tags$style(HTML("
      body {
        background-color: #f7f9fc;
        font-family: 'Segoe UI', 'Helvetica Neue', Arial, sans-serif;
        color: #333;
        padding-right: 10px;
      }
      #input_sidebarPanel {
        background-color: #ffffff;
        border-radius: 8px;
        box-shadow: 0px 2px 4px rgba(0,0,0,0.15);
        padding: 20px;
        margin-left: 10px;
      }
      #output_mainPanel {
        background-color: #ffffff;
        border-radius: 8px;
        box-shadow: 0px 6px 16px rgba(0,0,0,0.15);
        padding: 20px;
        /*height: 800px;*/
      }
      .scrollable-table {
        overflow-x: auto;
        overflow-y: auto;
        max-height: 1500px;
      }

      #titlePanel {
        font-size: 40px;
        font-weight: 600;
        color: #2e4e7e;
        margin-bottom: 0px;
        padding-bottom: 0px;
        padding-left: 15px;
      }
      #subtitle{
        font-size: 16px;
        color: #2e4e7e;
        margin-top: 0;
        padding-top: 0px;
        padding-left: 15px;
        padding-bottom: 10px;
      }
      th.rotate {
        /* Something you can count on */
        height: 180px;
        white-space: nowrap;
      }

      th.rotate > div {
        transform:
          /* Magic Numbers */
          translate(3px, 81px)
          /* 45 is really 360 - 45 */
          rotate(315deg);
        width: 30px;
      }
      th.rotate > div > span {
        padding: 5px 10px;
      }

    "))
    ),

    div(id = "titlePanel","USERM"),
    div(id = "subtitle","Unmixing Spread Estimation based on Residual Model"),
    sidebarLayout(
      sidebarPanel(
        id = "input_sidebarPanel",
        width = 3,
        tabsetPanel(id = "input_tabs",
                    tabPanel("Fluors",
                             tagList(
                               lapply(Userm$fluors, function(fluor) {
                                 checkboxInput(inputId = fluor, label = fluor, value = TRUE)
                               })
                             )
                    ),
                    tabPanel("Detectors",
                             tagList(
                               lapply(Userm$detectors, function(detector) {
                                 checkboxInput(inputId = detector, label = detector, value = TRUE)
                               })
                             )
                    ),
                    tabPanel("Intensity",
                             uiOutput("intensity_inputs")
                    )

        )
      ),
      mainPanel(
        id = "output_mainPanel",
        width = 9,
        tabsetPanel(

          tabPanel("Signature matrix",
                   div(class = "scrollable-table",
                       uiOutput("unmixing_matrix_html"),
                       p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                       fluidRow(
                         column(4, numericInput("unmixing_matrix_color_min", "Min", value = -1)),
                         column(4, numericInput("unmixing_matrix_color_mid", "Mid", value = 0)),
                         column(4, numericInput("unmixing_matrix_color_max", "Max", value = 1))
                       ),
                       actionButton(inputId = "unmixing_matrix_set_color_threshold", label = "Set")
                   )

          ),
          tabPanel("Prediction",
                   div(class = "scrollable-table",
                     selectInput("prediction_plot_mode", "Select display mode：", choices = NULL),

                     # X axis controls in one row
                     fluidRow(
                       column(4, selectInput("fluor_selector_x_prediction", "Select x axis：", choices = NULL)),
                       column(4, numericInput("x_min_input", "Min", value = 0)),
                       column(4, numericInput("x_max_input", "Max", value = 1000))
                     ),
                     fluidRow(
                       column(4, selectInput("x_scale", "Scale", choices = NULL)),
                       column(8, numericInput("x_cofactor", "Cofactor(only for Arcsinh)", value = 100, min = 1))
                     ),

                     # Y axis controls in one row
                     fluidRow(
                       column(4, selectInput("fluor_selector_y_prediction", "Select y axis：", choices = NULL)),
                       column(4, numericInput("y_min_input", "Min", value = 0)),
                       column(4, numericInput("y_max_input", "Max", value = 1000))
                     ),
                     fluidRow(
                       column(4, selectInput("y_scale", "Scale", choices = NULL)),
                       column(8, numericInput("y_cofactor", "Cofactor(only for Arcsinh)", value = 100, min = 1))
                     ),
                     numericInput(
                       inputId = paste0("intercept_factor"),
                       label = paste("Factor for default Autofluorescence (100: baseline; 0: no AF)"),
                       value = 100
                     ),
                     plotOutput("prediction_plot", width = "1000px", height = "800px")
                   )
          ),

          tabPanel("pinv matrix",
                   div(class = "scrollable-table",
                       uiOutput("pinv_matrix_html"),
                       p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                       fluidRow(
                         column(4, numericInput("pinv_matrix_color_min", "Min", value = -1)),
                         column(4, numericInput("pinv_matrix_color_mid", "Mid", value = 0)),
                         column(4, numericInput("pinv_matrix_color_max", "Max", value = 1))
                       ),
                       actionButton(inputId = "pinv_matrix_set_color_threshold", label = "Set")
                   )
          ),
          tabPanel("Intercept matrix",
                   tabsetPanel(
                     tabPanel("raw",
                              selectInput("fluor_selector_intercept_matrix_raw", "Select fluor SCC file：", choices = NULL),
                              div(class = "scrollable-table",
                                  uiOutput("intercept_matrix_raw_html"),
                                  p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                                  fluidRow(
                                    column(4, numericInput("intercept_matrix_raw_color_min", "Min", value = -10)),
                                    column(4, numericInput("intercept_matrix_raw_color_mid", "Mid", value = 0)),
                                    column(4, numericInput("intercept_matrix_raw_color_max", "Max", value = 10))
                                  ),
                                  actionButton(inputId = "intercept_matrix_raw_set_color_threshold", label = "Set")
                              )),
                     tabPanel("summary",
                              div(class = "scrollable-table",
                                  uiOutput("intercept_matrix_weighted_html"),
                                  p("*spread from row into column."),
                                  p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                                  fluidRow(
                                    column(4, numericInput("intercept_matrix_weighted_color_min", "Min", value = -10)),
                                    column(4, numericInput("intercept_matrix_weighted_color_mid", "Mid", value = 0)),
                                    column(4, numericInput("intercept_matrix_weighted_color_max", "Max", value = 10))
                                  ),
                                  actionButton(inputId = "intercept_matrix_weighted_set_color_threshold", label = "Set")
                              ))
                   )
          ),
          tabPanel("Slop matrix",
                   tabsetPanel(
                     tabPanel("raw",
                              selectInput("fluor_selector_slop_matrix_raw", "Select fluor SCC file：", choices = NULL),
                              div(class = "scrollable-table",
                                  uiOutput("slop_matrix_raw_html"),
                                  p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                                  fluidRow(
                                    column(4, numericInput("slop_matrix_raw_color_min", "Min", value = -1)),
                                    column(4, numericInput("slop_matrix_raw_color_mid", "Mid", value = 0)),
                                    column(4, numericInput("slop_matrix_raw_color_max", "Max", value = 1))
                                  ),
                                  actionButton(inputId = "slop_matrix_raw_set_color_threshold", label = "Set")
                              )),
                     tabPanel("weighted",
                              selectInput("fluor_selector_slop_matrix_weighted_from", "Select fluor SCC file：", choices = NULL),
                              selectInput("fluor_selector_slop_matrix_weighted_to", "Select spread channel：", choices = NULL),
                              div(class = "scrollable-table",
                                  uiOutput("slop_matrix_weighted_html"),
                                  p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                                  fluidRow(
                                    column(4, numericInput("slop_matrix_weighted_color_min", "Min", value = -1)),
                                    column(4, numericInput("slop_matrix_weighted_color_mid", "Mid", value = 0)),
                                    column(4, numericInput("slop_matrix_weighted_color_max", "Max", value = 1))
                                  ),
                                  actionButton(inputId = "slop_matrix_weighted_set_color_threshold", label = "Set")
                              )),
                     tabPanel("summary",
                              div(class = "scrollable-table",
                                  uiOutput("slop_matrix_summary_html"),
                                  p("*spread from row into column."),
                                  p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                                  fluidRow(
                                    column(4, numericInput("slop_matrix_summary_color_min", "Min", value = -1)),
                                    column(4, numericInput("slop_matrix_summary_color_mid", "Mid", value = 0)),
                                    column(4, numericInput("slop_matrix_summary_color_max", "Max", value = 1))
                                  ),
                                  actionButton(inputId = "slop_matrix_summary_set_color_threshold", label = "Set")
                              ))
                   )
          ),
          tabPanel("Coefficient matrix",
                   div(class = "scrollable-table",
                       uiOutput("coef_matrix_html"),
                       p("*spread from row into column."),
                       p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                       fluidRow(
                         column(4, numericInput("coef_matrix_color_min", "Min", value = -1)),
                         column(4, numericInput("coef_matrix_color_mid", "Mid", value = 0)),
                         column(4, numericInput("coef_matrix_color_max", "Max", value = 1))
                       ),
                       actionButton(inputId = "coef_matrix_set_color_threshold", label = "Set")
                   )
          ),
          tabPanel("Similarity matrix",
                   div(class = "scrollable-table",
                       uiOutput("similarity_matrix_html"),
                       p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                       fluidRow(
                         column(4, numericInput("similarity_matrix_color_min", "Min", value = -1)),
                         column(4, numericInput("similarity_matrix_color_mid", "Mid", value = 0)),
                         column(4, numericInput("similarity_matrix_color_max", "Max", value = 1))
                       ),
                       actionButton(inputId = "similarity_matrix_set_color_threshold", label = "Set")
                   )
          ),
          tabPanel("Hotspot matrix",
                   div(class = "scrollable-table",
                       uiOutput("hotspot_matrix_html"),
                       p("Set color threshold:", style = "font-size:18px; font-weight: bold; margin-top:20px;"),
                       fluidRow(
                         column(4, numericInput("hotspot_matrix_color_min", "Min", value = -1)),
                         column(4, numericInput("hotspot_matrix_color_mid", "Mid", value = 0)),
                         column(4, numericInput("hotspot_matrix_color_max", "Max", value = 1))
                       ),
                       actionButton(inputId = "hotspot_matrix_set_color_threshold", label = "Set")
                   )
          ),
          tabPanel("Export",
                   selectInput("prediction_plot_mode_report", "Select display mode：", choices = c("Pseudo-color","Contour line"), selected = "Pseudo-color"),
                   downloadButton("downloadHTML", "Export HTML report"),
                   downloadButton("downloadNxN", "Export N x N plot"),
                   p("Note: all parameters are from input Userm object."),
                   selectInput("prediction_Nx1_fluor", "Select fluor for Nx1 plot：", choices = NULL),
                   downloadButton("downloadNx1", "Export N x 1 plot")
          ),tabPanel(
            "About",

            fluidRow(
              column(12,
                     h3("About USERM"),
                     tags$hr()
              )
            ),

            fluidRow(
              column(12,
                     p("The USERM package provides an out-of-box tool to apply the residual model approach,
              which characterizes and predicts the spread of unmixed spectral flow cytometry data,
              which arises from instrumental noise or deviations between actual cellular emission
              and the average fluorescence signatures.")
              )
            ),
            fluidRow(
              column(12,
                     p("You can predict spread and interactively adjust various parameters.")
              )
            ),

            fluidRow(
              column(12,
                     p("To support panel design and spread interpretation,
                       the USERM also supports computing various matrices, including:"),
                     tags$ul(
                       tags$li("Coefficient Matrix"),
                       tags$li("Similarity Matrix"),
                       tags$li("Hotspot Matrix"),
                       tags$li("Signature Matrix"),
                       tags$li("pseudo-inverse Matrix"),
                       tags$li("Slop Matrix"),
                       tags$li("Intercept Matrix")

                     )
              )
            ),

            fluidRow(
              column(12,
                     p(
                       "Here is the link to ",
                       tags$a(href = "https://github.com/xiangmingcai/USERM",
                              "GitHub Repository", target = "_blank",
                              style = "color:#2c3e50; font-weight:bold;"),
                       ". You can find more instructions on the repository."
                     )
              )
            )
          )


        )
      )

    )
  )


  # Server
  server <- function(input, output, session) {

    # Initialize cache object

    #detector_cache
    detector_cache <- reactiveValues(selected = Userm$detectors)
    observe({
      selected <- Userm$detectors[unlist(lapply(Userm$detectors, function(det) input[[det]]))]
      detector_cache$selected <- selected
    })

    #fluor_cache
    fluor_cache <- reactiveValues(selected = Userm$fluors)
    observe({
      selected <- Userm$fluors[unlist(lapply(Userm$fluors, function(fluor) input[[fluor]]))]
      fluor_cache$selected <- selected
    })

    #intensity_cache
    intensity_cache <- reactiveValues()
    intensity_matrix = Userm$Intensity_mtx[,population_id,drop = FALSE]
    intensity_cache$matrix = intensity_matrix

    #intercept_matrix_cache
    intercept_matrix_cache <- reactiveValues()

    #intensity_cache
    lapply(Userm$fluors, function(fluor) {
      observeEvent(input[[paste0("intensity_", fluor)]], {
        intensity_cache$matrix[fluor,1] = input[[paste0("intensity_", fluor)]]
      }, ignoreNULL = TRUE)
    })
    observeEvent(input[["set_all_intensity"]], {
      intensity_cache$matrix[,1] = input[["intensity_value_for_all"]]
    })

    #observe and calculate intercept_matrix
    observe({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)

      weighted_matrix <- array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)

      for (i in 1:ncol(A)) {
        fluor = colnames(A)[i]
        intercept_matrix = Userm$Res[[fluor]]$interceptMtx
        intercept_matrix = intercept_matrix[detector_cache$selected,detector_cache$selected]
        weighted_matrix[i,] = diag((A_pinv %*% intercept_matrix) %*% t(A_pinv))
      }

      intercept_matrix_cache$matrix <- weighted_matrix
    })

    # Dynamically create numericInput components, using cached values as default when available.
    output$intensity_inputs <- renderUI({

      req(fluor_cache$selected)
      tagList(
        numericInput(
          inputId = "intensity_value_for_all",
          label = "Set Intensity for All",
          value = 0
        ),
        actionButton(inputId = "set_all_intensity", label = "Set"),
        lapply(fluor_cache$selected, function(fluor) {
          numericInput(
            inputId = paste0("intensity_", fluor),
            label = paste("Intensity for", fluor),
            value = intensity_cache$matrix[fluor,1]

          )
        })
      )
    })

    #color_threshold
    color_thre_mtx_init = matrix(nrow = 10,ncol = 3)
    rownames(color_thre_mtx_init) = c("signature_mtx","pinv_mtx","intercept_mtx_raw","intercept_mtx_weighted",
                                     "slop_mtx_raw","slop_mtx_weighted","slop_mtx_summary",
                                     "coef_mtx","similarity_mtx","hotspot_mtx")
    colnames(color_thre_mtx_init) = c("min","mid","max")
    color_thre_mtx_init["signature_mtx",] = c(-1, 0, 1)
    color_thre_mtx_init["pinv_mtx",] = c(-1, 0, 1)
    color_thre_mtx_init["intercept_mtx_raw",] = c(-10, 0, 10)
    color_thre_mtx_init["intercept_mtx_weighted",] = c(-10, 0, 10)
    color_thre_mtx_init["slop_mtx_raw",] = c(-1, 0, 1)
    color_thre_mtx_init["slop_mtx_weighted",] = c(-1, 0, 1)
    color_thre_mtx_init["slop_mtx_summary",] = c(-1, 0, 1)
    color_thre_mtx_init["coef_mtx",] = c(-1, 0, 1)
    color_thre_mtx_init["similarity_mtx",] = c(-1, 0, 1)
    color_thre_mtx_init["hotspot_mtx",] = c(-1, 0, 1)

    color_cache <- reactiveValues()
    color_cache$matrix = color_thre_mtx_init

    observeEvent(input[["unmixing_matrix_set_color_threshold"]], {
      color_cache$matrix["signature_mtx","min"] = input[["unmixing_matrix_color_min"]]
      color_cache$matrix["signature_mtx","mid"] = input[["unmixing_matrix_color_mid"]]
      color_cache$matrix["signature_mtx","max"] = input[["unmixing_matrix_color_max"]]
    })
    observeEvent(input[["pinv_matrix_set_color_threshold"]], {
      color_cache$matrix["pinv_mtx","min"] = input[["pinv_matrix_color_min"]]
      color_cache$matrix["pinv_mtx","mid"] = input[["pinv_matrix_color_mid"]]
      color_cache$matrix["pinv_mtx","max"] = input[["pinv_matrix_color_max"]]
    })
    observeEvent(input[["intercept_matrix_raw_set_color_threshold"]], {
      color_cache$matrix["intercept_mtx_raw","min"] = input[["intercept_matrix_raw_color_min"]]
      color_cache$matrix["intercept_mtx_raw","mid"] = input[["intercept_matrix_raw_color_mid"]]
      color_cache$matrix["intercept_mtx_raw","max"] = input[["intercept_matrix_raw_color_max"]]
    })
    observeEvent(input[["intercept_matrix_weighted_set_color_threshold"]], {
      color_cache$matrix["intercept_mtx_weighted","min"] = input[["intercept_matrix_weighted_color_min"]]
      color_cache$matrix["intercept_mtx_weighted","mid"] = input[["intercept_matrix_weighted_color_mid"]]
      color_cache$matrix["intercept_mtx_weighted","max"] = input[["intercept_matrix_weighted_color_max"]]
    })
    observeEvent(input[["slop_matrix_raw_set_color_threshold"]], {
      color_cache$matrix["slop_mtx_raw","min"] = input[["slop_matrix_raw_color_min"]]
      color_cache$matrix["slop_mtx_raw","mid"] = input[["slop_matrix_raw_color_mid"]]
      color_cache$matrix["slop_mtx_raw","max"] = input[["slop_matrix_raw_color_max"]]
    })
    observeEvent(input[["slop_matrix_weighted_set_color_threshold"]], {
      color_cache$matrix["slop_mtx_weighted","min"] = input[["slop_matrix_weighted_color_min"]]
      color_cache$matrix["slop_mtx_weighted","mid"] = input[["slop_matrix_weighted_color_mid"]]
      color_cache$matrix["slop_mtx_weighted","max"] = input[["slop_matrix_weighted_color_max"]]
    })
    observeEvent(input[["slop_matrix_summary_set_color_threshold"]], {
      color_cache$matrix["slop_mtx_summary","min"] = input[["slop_matrix_summary_color_min"]]
      color_cache$matrix["slop_mtx_summary","mid"] = input[["slop_matrix_summary_color_mid"]]
      color_cache$matrix["slop_mtx_summary","max"] = input[["slop_matrix_summary_color_max"]]
    })
    observeEvent(input[["coef_matrix_set_color_threshold"]], {
      color_cache$matrix["coef_mtx","min"] = input[["coef_matrix_color_min"]]
      color_cache$matrix["coef_mtx","mid"] = input[["coef_matrix_color_mid"]]
      color_cache$matrix["coef_mtx","max"] = input[["coef_matrix_color_max"]]
    })
    observeEvent(input[["similarity_matrix_set_color_threshold"]], {
      color_cache$matrix["similarity_mtx","min"] = input[["similarity_matrix_color_min"]]
      color_cache$matrix["similarity_mtx","mid"] = input[["similarity_matrix_color_mid"]]
      color_cache$matrix["similarity_mtx","max"] = input[["similarity_matrix_color_max"]]
    })
    observeEvent(input[["hotspot_matrix_set_color_threshold"]], {
      color_cache$matrix["hotspot_mtx","min"] = input[["hotspot_matrix_color_min"]]
      color_cache$matrix["hotspot_mtx","mid"] = input[["hotspot_matrix_color_mid"]]
      color_cache$matrix["hotspot_mtx","max"] = input[["hotspot_matrix_color_max"]]
    })

    colormin = "#4473c5"#098ebb  2166ac
    colormid = "#f7f7f7"
    colormax = "#ef7e30" #e96449  b2182b



    # unmixing_matrix_table
    output$unmixing_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      A = Userm$A
      mat = A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      HTML(Userm_html_table(mat = mat,
                            val_min = color_cache$matrix["signature_mtx","min"],
                            val_mid = color_cache$matrix["signature_mtx","mid"],
                            val_max = color_cache$matrix["signature_mtx","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })

    # pinv_matrix_table
    output$pinv_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      A = Userm$A
      mat = A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      mat_pinv = ginv(mat)
      colnames(mat_pinv) = detector_cache$selected
      rownames(mat_pinv) = fluor_cache$selected
      HTML(Userm_html_table(mat = mat_pinv,
                            val_min = color_cache$matrix["pinv_mtx","min"],
                            val_mid = color_cache$matrix["pinv_mtx","mid"],
                            val_max = color_cache$matrix["pinv_mtx","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })

    # intercept_matrix_table
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "fluor_selector_intercept_matrix_raw",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
    })
    output$intercept_matrix_raw_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      intercept_matrix = Userm$Res[[input$fluor_selector_intercept_matrix_raw]]$interceptMtx
      intercept_matrix = intercept_matrix[detector_cache$selected,detector_cache$selected]
      HTML(Userm_html_table(mat = intercept_matrix,
                            val_min = color_cache$matrix["intercept_mtx_raw","min"],
                            val_mid = color_cache$matrix["intercept_mtx_raw","mid"],
                            val_max = color_cache$matrix["intercept_mtx_raw","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })
    output$intercept_matrix_weighted_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]

      req(intercept_matrix_cache$matrix)
      weighted_matrix = intercept_matrix_cache$matrix

      final_row = matrix(apply(weighted_matrix,1,median),nrow = 1)
      sqrt_row = sqrt(abs(final_row))
      final_matrix = rbind(weighted_matrix, final_row)
      final_matrix = rbind(final_matrix, sqrt_row)
      rownames(final_matrix) = c(colnames(A),"median(Sigma^2)","abs(Sigma)")
      colnames(final_matrix) = colnames(A)

      HTML(Userm_html_table(mat = final_matrix,
                            val_min = color_cache$matrix["intercept_mtx_weighted","min"],
                            val_mid = color_cache$matrix["intercept_mtx_weighted","mid"],
                            val_max = color_cache$matrix["intercept_mtx_weighted","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })

    # slop_matrix_table
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "fluor_selector_slop_matrix_raw",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
      updateSelectInput(session, "fluor_selector_slop_matrix_weighted_from",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
      updateSelectInput(session, "fluor_selector_slop_matrix_weighted_to",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
    })
    output$slop_matrix_raw_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)
      slop_matrix = Userm$Res[[input$fluor_selector_slop_matrix_raw]]$slopMtx
      slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]
      HTML(Userm_html_table(mat = slop_matrix,
                            val_min = color_cache$matrix["slop_mtx_raw","min"],
                            val_mid = color_cache$matrix["slop_mtx_raw","mid"],
                            val_max = color_cache$matrix["slop_mtx_raw","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })
    output$slop_matrix_weighted_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      slop_matrix = Userm$Res[[input$fluor_selector_slop_matrix_weighted_from]]$slopMtx
      slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)
      A_weight_mtx = A_pinv[input$fluor_selector_slop_matrix_weighted_to,]%o%A_pinv[input$fluor_selector_slop_matrix_weighted_to,]
      weighted_slop_matrix = A_weight_mtx * slop_matrix * intensity_cache$matrix[input$fluor_selector_slop_matrix_weighted_from,1]

      HTML(Userm_html_table(mat = weighted_slop_matrix,
                            val_min = color_cache$matrix["slop_mtx_weighted","min"],
                            val_mid = color_cache$matrix["slop_mtx_weighted","mid"],
                            val_max = color_cache$matrix["slop_mtx_weighted","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })
    output$slop_matrix_summary_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)

      weighted_matrix = array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)

      for (i in 1:ncol(A)) {
        fluor = colnames(A)[i]
        slop_matrix = Userm$Res[[fluor]]$slopMtx
        slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]
        weighted_matrix[i,] = diag((A_pinv %*% slop_matrix) %*% t(A_pinv)) * intensity_cache$matrix[fluor,1]
      }

      final_row = matrix(apply(weighted_matrix,1,sum),nrow = 1)
      sqrt_row = sqrt(abs(final_row))
      final_matrix = rbind(weighted_matrix, final_row)
      final_matrix = rbind(final_matrix, sqrt_row)
      rownames(final_matrix) = c(colnames(A),"sum(Sigma^2)","abs(Sigma)")
      colnames(final_matrix) = colnames(A)

      HTML(Userm_html_table(mat = final_matrix,
                            val_min = color_cache$matrix["slop_mtx_summary","min"],
                            val_mid = color_cache$matrix["slop_mtx_summary","mid"],
                            val_max = color_cache$matrix["slop_mtx_summary","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })

    # estimation mtx
    # coef_matrix_html
    output$coef_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)

      #calculate weighted slop
      Coef_matrix = array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)
      for (i in 1:ncol(A)) {
        fluor = colnames(A)[i]
        slop_matrix = Userm$Res[[fluor]]$slopMtx
        slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]
        Coef_matrix[,i] = diag((A_pinv %*% slop_matrix) %*% t(A_pinv))
      }
      colnames(Coef_matrix) = colnames(A)
      rownames(Coef_matrix) = colnames(A)
      Coef_matrix = t(Coef_matrix)

      for (row in 1:nrow(Coef_matrix)) {
        for (col in 1:ncol(Coef_matrix)) {
          if(Coef_matrix[row,col] < 0){
            Coef_matrix[row,col] = -sqrt(abs(Coef_matrix[row,col]))
          }else if(Coef_matrix[row,col] > 0){
            Coef_matrix[row,col] = sqrt(Coef_matrix[row,col])
          }
        }
      }
      diag(Coef_matrix) = 0# residual model cannot predict spread at the positive axis

      HTML(Userm_html_table(mat = Coef_matrix,
                            val_min = color_cache$matrix["coef_mtx","min"],
                            val_mid = color_cache$matrix["coef_mtx","mid"],
                            val_max = color_cache$matrix["coef_mtx","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })
    #similarity_matrix_html
    output$similarity_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      cos_sim_matrix <- cosine(A)
      colnames(cos_sim_matrix) = colnames(A)
      rownames(cos_sim_matrix) = colnames(A)

      HTML(Userm_html_table(mat = cos_sim_matrix,
                            val_min = color_cache$matrix["similarity_mtx","min"],
                            val_mid = color_cache$matrix["similarity_mtx","mid"],
                            val_max = color_cache$matrix["similarity_mtx","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })
    #hotspot_matrix_html
    output$hotspot_matrix_html <- renderUI({
      req(detector_cache$selected, fluor_cache$selected)

      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)
      H_mtx = (t(A) %*% A)
      H_mtx = ginv(H_mtx)
      H_mtx = sqrt(abs(H_mtx))
      colnames(H_mtx) = colnames(A)
      rownames(H_mtx) = colnames(A)

      HTML(Userm_html_table(mat = H_mtx,
                            val_min = color_cache$matrix["hotspot_mtx","min"],
                            val_mid = color_cache$matrix["hotspot_mtx","mid"],
                            val_max = color_cache$matrix["hotspot_mtx","max"],
                            colormin = colormin, colormid = colormid, colormax = colormax))
    })

    #prediction
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "prediction_plot_mode",
                        choices = c("Pseudo-color","Contour line"),
                        selected = "Pseudo-color")
      updateSelectInput(session, "fluor_selector_x_prediction",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
      updateSelectInput(session, "fluor_selector_y_prediction",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
      updateSelectInput(session, "x_scale",
                        choices = c("Linear","Log10","Arcsinh"),
                        selected = "Linear")
      updateSelectInput(session, "y_scale",
                        choices = c("Linear","Log10","Arcsinh"),
                        selected = "Linear")

    })
    output$prediction_plot <- renderPlot({
      req(detector_cache$selected, fluor_cache$selected)
      req(intercept_matrix_cache$matrix)
      intercept_matrix = intercept_matrix_cache$matrix
      intercept_matrix = t(intercept_matrix)
      intercept_col = matrix(apply(intercept_matrix,1,median),ncol = 1)

      # A = Userm$A
      A = Userm$A[detector_cache$selected, fluor_cache$selected, drop = FALSE]
      A_pinv = ginv(A)
      colnames(A_pinv) = rownames(A)
      rownames(A_pinv) = colnames(A)

      weighted_matrix = array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)
      intensity_col = array(0, dim = c(ncol(A), 1)) #(channel, 1)
      for (i in 1:ncol(A)) {
        fluor = colnames(A)[i]
        slop_matrix = Userm$Res[[fluor]]$slopMtx
        slop_matrix = slop_matrix[detector_cache$selected,detector_cache$selected]
        weighted_matrix[,i] = diag((A_pinv %*% slop_matrix) %*% t(A_pinv)) * intensity_cache$matrix[fluor,1]
        intensity_col[i,1] = intensity_cache$matrix[fluor,1]
      }
      spread_col = matrix(apply(weighted_matrix,1,sum),ncol = 1)

      final_col = spread_col + intercept_col * (input$intercept_factor / 100)

      sqrt_col = sqrt(abs(final_col))
      rownames(sqrt_col) = colnames(A)
      rownames(intensity_col) = colnames(A)
      #Calculate iso band data

      mu_x = intensity_col[input$fluor_selector_x_prediction,1]
      mu_y = intensity_col[input$fluor_selector_y_prediction,1]
      sigma_x = sqrt_col[input$fluor_selector_x_prediction,1]
      sigma_y = sqrt_col[input$fluor_selector_y_prediction,1]

      # create grid
      x <- create_seq(min = input$x_min_input,
                      max = input$x_max_input, len=100,
                      scale = input$x_scale, cofactor = input$x_cofactor)
      y <- create_seq(min = input$y_min_input,
                      max = input$y_max_input, len=100,
                      scale = input$y_scale, cofactor = input$y_cofactor)
      grid <- expand.grid(x = x, y = y)

      # calculate intensity for x and y
      fx <- dnorm(grid$x, mean = mu_x, sd = sigma_x)
      fy <- dnorm(grid$y, mean = mu_y, sd = sigma_y)

      # Joint probability density (due to independence, it equals the product of individual densities)
      grid$z <- fx * fy

      #visualize iso band data
      p = PredIsobandPlot(grid = grid,
                          x_scale = input$x_scale,
                          y_scale = input$y_scale,
                          x_cofactor = input$x_cofactor,
                          y_cofactor = input$y_cofactor,
                          x_label = input$fluor_selector_x_prediction,
                          y_label = input$fluor_selector_y_prediction,
                          x_min = input$x_min_input,
                          x_max = input$x_max_input,
                          y_min = input$y_min_input,
                          y_max = input$y_max_input,
                          mode = input$prediction_plot_mode)
      p

    })

    #Export
    observe({
      req(fluor_cache$selected)
      updateSelectInput(session, "prediction_Nx1_fluor",
                        choices = fluor_cache$selected,
                        selected = fluor_cache$selected[1])
    })
    #downloadHTML
    output$downloadHTML <- downloadHandler(
      filename = function() {
        paste0("UsermOne_report-", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"), ".html")
      },
      content = function(file) {
        tempReport <- system.file("report_tmp", "UsermOnereport.Rmd", package = "USERM")
        # tempReport <- file.path(report_tmp_dir, "UsermOnereport.Rmd")
        file.copy("UsermOnereport.Rmd", tempReport, overwrite = TRUE)

        params <- list(Userm = Userm,
                       detector_selected = detector_cache$selected,
                       fluor_selected = fluor_cache$selected,
                       prediction_plot_mode_report = input$prediction_plot_mode_report,
                       population_id = population_id)

        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          output_format = "html_document",
                          envir = new.env(parent = globalenv()))
      }
    )
    #downloadNxN
    output$downloadNxN <- downloadHandler(
      filename = function() {
        paste0("UsermOne_NxN-", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"), ".html")
      },
      content = function(file) {
        tempReport <- system.file("report_tmp", "UsermOneNxN.Rmd", package = "USERM")
        # tempReport <- file.path(report_tmp_dir, "UsermOneNxN.Rmd")
        file.copy("UsermOneNxN.Rmd", tempReport, overwrite = TRUE)

        params <- list(Userm = Userm,
                       detector_selected = detector_cache$selected,
                       fluor_selected = fluor_cache$selected,
                       prediction_plot_mode_report = input$prediction_plot_mode_report,
                       population_id = population_id)

        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          output_format = "html_document",
                          envir = new.env(parent = globalenv()))
      }
    )
    #downloadNx1
    output$downloadNx1 <- downloadHandler(
      filename = function() {
        paste0("UsermOne_Nx1-", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"), ".html")
      },
      content = function(file) {
        tempReport <- system.file("report_tmp", "UsermOneNx1.Rmd", package = "USERM")
        # tempReport <- file.path(report_tmp_dir, "UsermOneNx1.Rmd")
        file.copy("UsermOneNx1.Rmd", tempReport, overwrite = TRUE)

        params <- list(Userm = Userm,
                       detector_selected = detector_cache$selected,
                       fluor_selected = fluor_cache$selected,
                       prediction_plot_mode_report = input$prediction_plot_mode_report,
                       x_fluor = input$prediction_Nx1_fluor,
                       population_id = population_id)

        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          output_format = "html_document",
                          envir = new.env(parent = globalenv()))
      }
    )

  }

  # running App
  shinyApp(ui = ui, server = server)

}
