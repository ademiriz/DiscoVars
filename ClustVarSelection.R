
library(shiny)
library(datamods)
library(shinythemes)
library(DT)
require(shinydashboard)
require(shinyjs)
library(doParallel)
library(parallel)
library(sna)
library(tidygraph)
library(tidyverse)
library(ggraph)
library(igraph)
library(mlbench)
library(dplyr)
library(mclust)
if (!require(factoextra)) install.packages("factoextra")
library(factoextra)
if (!require(rlang)) install.packages("rlang")
library(rlang)
if (!require(clustvarsel)) install.packages("clustvarsel")
library(clustvarsel)
if (!require(sparcl)) install.packages("sparcl")
library(sparcl)
if (!require(RSKC)) install.packages("RSKC")
library(RSKC)

#setwd("C:\\Users\\AyhanDemiriz\\Documents\\MyRFolder")

source("helpers.R", local = TRUE)



ui <-dashboardPage(
  title='Variable Select',
  dashboardHeader(
    title = 'Variable Select'
  ),
  dashboardSidebar(
   # wellPanel(
      div(id='tab1_sidebar',
                         checkboxGroupInput(
                           inputId = "from",
                           label = "From",
                           choices = c("env", "file", "copypaste", "googlesheets"),
                           selected = c("env", "file", "copypaste", "googlesheets")
                         ),
                         actionButton("imprtData", "Import Data"),
          conditionalPanel(condition="output.dataLoaded",
          radioButtons("rb1", "Direction Method",
                       choiceNames = list(
                         HTML("<p style='color:blue;'>Backward</p>"),
                         HTML("<p style='color:blue;'>Forward</p>")
                       ),
                       choiceValues = list(
                         "backward", "forward"
                       )),

            sliderInput("numClus", "Number of Clusters", 3, 15,
                        value = 9, step = 1

          ),

                           withBusyIndicatorUI(
                             actionButton("runclustVarSel", "Run clustVarSel",
                                          class = "btn-primary")
                           )


          ),
          conditionalPanel(condition="output.dataLoaded",
          sliderInput("numClus2", "Number of Clusters - RSKC", 3, 15,
                      value = 9, step = 1

          ),
          sliderInput("alpha", "Alpha - RSKC", 0, 0.3,
                      value = 0.05, step = 0.01

          ),
          sliderInput("L1Param", "L1 - RSKC", 1, 10,
                      value = 1.25, step = 0.01

          ),


                           withBusyIndicatorUI(
                             actionButton("runSparcl", "Run RSKC",
                                          class = "btn-primary")
                           )


          )




          )

  ),
  dashboardBody(
    useShinyjs(),
    tabsetPanel(
      id = "navbar",
      tabPanel(title="Data",id="tab1",value='tab1_val',
               column(
                 width = 8,
                 tags$b("Imported data:"),
                 verbatimTextOutput(outputId = "name"),
                 DT::dataTableOutput('data')
               )),
      tabPanel(title="ClusterVarSel",id="tab2",value='tab2_val',
               conditionalPanel(condition="output.dataLoaded",

                                column(
                                  width = 8,
                                  tags$b("Finalized for clustVarSel"),
                                  verbatimTextOutput(outputId = "txtClustVars1"),
                                  verbatimTextOutput("clustVarSelOut"),
                                  # DT::dataTableOutput('DTclustData')
                                  wellPanel(plotOutput('ClustVarSelPlot',width = "100%"))
                                )


               )


               ),
      tabPanel(title="RSKCSel",id="tab3",value='tab3_val',
               conditionalPanel(condition="output.dataLoaded",

                                column(
                                  width = 8,
                                  tags$b("Finalized for RSKC"),
                                  verbatimTextOutput(outputId = "txtClustVars2"),
                                  verbatimTextOutput("sparclOut"),
                                  # DT::dataTableOutput('DTclustData')
                                  wellPanel(plotOutput('sparclPlot',width = "100%"))
                                )


               )


      )
    )
  )
)


server <- function(input, output, session) {

  observeEvent(input$imprtData, {
    req(input$from)
    import_modal(
      id = "myid",
      from = input$from,
      title = "Import data to be used in application"
    )
  })

  imported <- import_server("myid", return_class = "data.frame")

  #constructDepNet <- reactiveVal()

  output$dataLoaded<-reactive(!is.null(imported$data()))
  outputOptions(output, 'dataLoaded', suspendWhenHidden = FALSE)

  output$name <- renderPrint({
    req(imported$name())
    imported$name()
  })

  output$data = DT::renderDataTable(imported$data(), server = TRUE)


  observeEvent(input$runclustVarSel, { withBusyIndicatorServer("runclustVarSel",{

    output$ClustVarSelPlot <- renderPlot({
      tbdata<-req(imported$data())
      tbdata<-select_if(tbdata, is.numeric)
      input$numclus
      if(!is.null(tbdata)){


          out<-clustvarsel(as.data.frame(tbdata),G = input$numClus,direction=input$rb1, parallel = TRUE)


        output$txtClustVars1 <- renderPrint({
          out$subset

        })
        output$clustVarSelOut <- renderPrint({
          summary(out$model)
        })

        fviz_mclust(out$model, "classification", geom = "point")

      }
    }, height=500)
  })
  })

  observeEvent(input$runSparcl, { withBusyIndicatorServer("runSparcl",{

    output$sparclPlot <- renderPlot({
      tbdata<-req(imported$data())
      tbdata<-select_if(tbdata, is.numeric)
      input$numclus2
      if(!is.null(tbdata)){


        r2<-RSKC(tbdata, ncl=input$numClus2, alpha=input$alpha, L1 = input$L1Param)

        selvars<-r2$weight[r2$weights>0.0000001]

        r2data<-tbdata[names(selvars)]


        output$txtClustVars2 <- renderPrint({
          names(selvars)

        })
        output$sparclOut <- renderPrint({
          print(r2)
        })



        #plot(fit,what="classification")
        fviz_cluster(object=list(data = r2data, cluster = as.vector(r2$labels)),geom = "point")

      }
    }, height=500)
  })
  })


  }


if (interactive())
  shinyApp(ui, server)
