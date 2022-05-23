
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
library(netrankr)
library(statnet)
library(MASS)
library(igraph)
library(mlbench)
library(dplyr)
library(olsrr)
library(grid)
library(plotly)
library(visNetwork)
library(mclust)
if (!require(factoextra)) install.packages("factoextra")
library(factoextra)
if (!require(rlang)) install.packages("rlang")
library(rlang)

if (!require(clustvarsel)) install.packages("clustvarsel")
library(clustvarsel)

if (!require(glmnet)) install.packages("glmnet")
library(glmnet)

#setwd("C:\\Users\\AyhanDemiriz\\Documents\\MyRFolder")

source("helpers.R", local = TRUE)

model.MyolsSTW <- function(x) {

  library(olsrr)
  response <- names(tbdata)[x]
  nmsregress<-names(tbdata)[-x]
  regressors <- paste0(nmsregress,"+")
  regressors[d-1]=nmsregress[d-1]
  regformula <- formula(paste(c(response, "~", regressors,"-1"),collapse=" "))
  fit <- lm(regformula,data=tbdata)
  k <- ols_step_both_p(fit, pent = 0.10, prem = 0.25)
  return(k$model$coefficient)
}

model.MyolsFRW <- function(x) {

  library(olsrr)
  response <- names(tbdata)[x]
  nmsregress<-names(tbdata)[-x]
  regressors <- paste0(nmsregress,"+")
  regressors[d-1]=nmsregress[d-1]
  regformula <- formula(paste(c(response, "~", regressors,"-1"),collapse=" "))
  fit <- lm(regformula,data=tbdata)
  k <- ols_step_forward_p(fit, pent = 0.1)
  return(k$predictors[1:length(k$rsquare)])
}

model.MyolsAIC <- function(x) {

  library(MASS)
  response <- names(tbdata)[x]
  nmsregress<-names(tbdata)[-x]
  regressors <- paste0(nmsregress,"+")
  regressors[d-1]=nmsregress[d-1]
  regformula <- formula(paste(c(response, "~", regressors,"-1"),collapse=" "))
  fit <- lm(regformula,data=tbdata)
  k <- stepAIC(fit,direction="both",trace=FALSE)
  return(k$coefficients)
}

model.MyolsLasso <- function(x) {
  library(glmnet)
  response <- names(tbdata)[x]
  nmsregress<-names(tbdata)[-x]
  x<-makeX(tbdata[nmsregress])
  y<-makeX(tbdata[response])
  nrw<-nrow(tbdata)
 # ncl<-ncol(tbdata)
  fit <- glmnet(x, y,intercept = TRUE)
  mycoef<-as.data.frame(summary(coef(fit, s = 32.0 / (2*nrw))))
  mycoef<-mycoef[-1,] #remove intercept
  mydf<-data.frame(t(mycoef$x))
  names(mydf)<-nmsregress[mycoef$i-1] #intercept is the first variable 
  return(mydf)
}

centChoices<-c('alpha','authority','betweenness','closeness','degree','eigen','hub','pagerank','power')


wssplot <- function(data, nc = 15, seed = 6754) {
  wss <- (nrow(data) - 1) * sum(apply(data, 2, var))
  for (i in 2:nc) {
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers = i)$withinss)
  }
  plot(1:nc,
       wss,
       type = "b",
       xlab = "Number of Clusters",
       ylab = "Within groups sum of squares")
}



ui <-dashboardPage(
  title='DiscoVars',
  dashboardHeader(
    title = 'DiscoVars'
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
          radioButtons("rb", "Variable Selection Method",
                       choiceNames = list(
                         HTML("<p style='color:blue;'>Stepwise</p>"),
                         HTML("<p style='color:blue;'>Forward</p>"),
                         HTML("<p style='color:blue;'>StepAIC</p>"),
                         HTML("<p style='color:blue;'>Lasso</p>")
                       ),
                       choiceValues = list(
                         "stw", "fwd", "aic","las"
                       )),
          textOutput("txt"),
          conditionalPanel(condition="output.dataLoaded",
                           withBusyIndicatorUI(
                             actionButton("runAlgo", "Run Dependency Discoverer",
                                          class = "btn-primary")
                           )


          ),
          conditionalPanel(condition="output.graphReady",
                           withBusyIndicatorUI(
                             actionButton("finalizeVars", "Finalize Selected Variables",
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
      tabPanel(title="Dependencies",id="tab2",value='tab2_val',

               column(2,
                      fluidRow(
                        wellPanel(
                          selectInput(inputId = "CentralMeasure", label="Centrality Measure", choices = centChoices, selected = "alpha", multiple = FALSE, selectize = TRUE)
                        )
                      ),
                      fluidRow(
                        wellPanel(
                          sliderInput("topVars", "Top Variables", 1, 20,
                                      value = 5, step = 1
                          )
                      )
                      ),
                      fluidRow(
                        wellPanel(
                           dataTableOutput('resTopVars')
                        )
                      )
               ),
               column(10,
                      wellPanel(plotOutput('DepNet',width = "100%")
                      ))


               ),
      tabPanel(title="Cluster Data I",id="tab3",value='tab3_val',
               column(
                 width = 8,
                 tags$b("Finalized for Clustering Analysis - mclust"),
                 verbatimTextOutput(outputId = "txtClustVars"),
                 verbatimTextOutput("mclustModelOut"),
                # DT::dataTableOutput('DTclustData')
                wellPanel(plotOutput('ClustPlot',width = "100%")
                )
               )),
      tabPanel(title="Cluster Data II",id="tab4",value='tab4_val',
               column(
                 width = 8,
                 tags$b("Finalized for Clustering Analysis - kmeans"),
                 verbatimTextOutput(outputId = "txtClustVars2"),
                 wellPanel(plotOutput('kmeansClustPlot',width = "100%"))


               ),
               column(
                 width = 4,
                 numericInput('ksolution', 'Choose k for kmeans', 5),
                 wellPanel(plotOutput('elbowPlot',width = "100%"))

               )
               ),

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

  constructDepNet <- reactiveVal()

  output$dataLoaded<-reactive(!is.null(imported$data()))
  outputOptions(output, 'dataLoaded', suspendWhenHidden = FALSE)

  output$name <- renderPrint({
    req(imported$name())
    imported$name()
  })

  output$data = DT::renderDataTable(imported$data(), server = TRUE)


  output$txt <- renderText({
    paste("You chose", input$rb)
  })


  observeEvent(input$runAlgo, { withBusyIndicatorServer("runAlgo",{
    tbdata<-req(imported$data())
    tbdata<-select_if(tbdata, is.numeric)

    d<-length(tbdata)
    if(d-2>5)
      {updateSliderInput(session, "topVars",min=2,max=d-2,value=5)}
    else
      {updateSliderInput(session, "topVars", max=d-2,value=2)}


    no_cores <- detectCores(logical = TRUE)

    if(input$rb=="stw"){
      system.time({
        clust <- makeCluster( no_cores-1, 'PSOCK')
        registerDoParallel(clust)
        clusterExport(clust, varlist=c("tbdata","d"),envir=environment())
        reglist <- parSapply(clust, 1:d, model.MyolsSTW)})

      stopCluster(clust)
      pnet<-matrix(data=0,ncol=d,nrow=d)
      colnames(pnet)<-names(tbdata)
      rownames(pnet)<-names(tbdata)
      for(i in 1:d){
        jj<-unlist(reglist[i])
        nms_jj<-names(jj)
        pnet[nms_jj[-1],i]<-1
        }

    }

    if(input$rb=="fwd"){
      system.time({
        clust <- makeCluster( no_cores-1, 'PSOCK')
        registerDoParallel(clust)
        clusterExport(clust, varlist=c("tbdata","d"),envir=environment())
        reglist <- parSapply(clust, 1:d, model.MyolsFRW)})

      stopCluster(clust)
      pnet<-matrix(data=0,ncol=d,nrow=d)
      colnames(pnet)<-names(tbdata)
      rownames(pnet)<-names(tbdata)
      for(i in 1:d){
        jj<-unlist(reglist[i])

        pnet[jj,i]<-1
      }

    }
    if(input$rb=="aic"){
      system.time({
        clust <- makeCluster( no_cores-1, 'PSOCK')
        registerDoParallel(clust)
        clusterExport(clust, varlist=c("tbdata","d"),envir=environment())
        reglist <- parSapply(clust, 1:d, model.MyolsAIC)})

      stopCluster(clust)

      pnet<-matrix(data=0,ncol=d,nrow=d)
      colnames(pnet)<-names(tbdata)
      rownames(pnet)<-names(tbdata)
      for(i in 1:d){
        jj<-unlist(reglist[i])
        nms_jj<-names(jj)
        pnet[nms_jj,i]<-1
      }

    }

    if(input$rb=="las"){
      system.time({
        clust <- makeCluster( no_cores-1, 'PSOCK')
        registerDoParallel(clust)
        clusterExport(clust, varlist=c("tbdata","d"),envir=environment())
        reglist <- parSapply(clust, 1:d, model.MyolsLasso)})

      stopCluster(clust)

      pnet<-matrix(data=0,ncol=d,nrow=d)
      colnames(pnet)<-names(tbdata)
      rownames(pnet)<-names(tbdata)
      for(i in 1:d){
        jj<-unlist(reglist[i])
        nms_jj<-names(jj)
        pnet[nms_jj,i]<-1
      }

    }

    constructDepNet(pnet)

  })
  })



  output$DepNet <- renderPlot({
    if(!is.null(constructDepNet())){
      methodChoice <- switch(input$CentralMeasure,
                             'alpha' = centrality_alpha,
                             'authority' = centrality_authority,
                             'betweenness' = centrality_betweenness,
                             'closeness' = centrality_closeness,
                             'degree' = centrality_degree,
                             'eigen'=centrality_eigen,
                             'hub'=centrality_hub,
                             'pagerank'=centrality_pagerank,
                             'power'=centrality_power

      )
    set.seed(1)
    pnet.graph <- as_tbl_graph(constructDepNet(), directed = TRUE)
    ggraph(pnet.graph) +
      geom_edge_link() +
      geom_node_point() +
      geom_node_text(
        aes(label = name), size = 3, repel = TRUE
      ) +
      theme_graph()


    set.seed(123)
    pnet.graph %>%
      activate(nodes) %>%
      mutate(centrality = methodChoice()) %>%
      ggraph(layout = "graphopt") +
      geom_edge_link(width = 1, colour = "lightgray") +
      geom_node_point(aes(size = centrality, colour = centrality)) +
      geom_node_text(aes(label = name), repel = TRUE)+
      scale_color_gradient(low = "yellow", high = "red")+
      theme_graph()
    }

  }, height=600)


  TopVars<-reactive({
    if(is.null(constructDepNet())) return()
    methodChoice <- switch(input$CentralMeasure,
                    'alpha' = centrality_alpha,
                    'authority' = centrality_authority,
                    'betweenness' = centrality_betweenness,
                    'closeness' = centrality_closeness,
                    'degree' = centrality_degree,
                    'eigen'=centrality_eigen,
                    'hub'=centrality_hub,
                    'pagerank'=centrality_pagerank,
                    'power'=centrality_power

                    )
    set.seed(123)
    pnet.graph <- as_tbl_graph(constructDepNet(), directed = TRUE)
    ggraph(pnet.graph) +
      geom_edge_link() +
      geom_node_point() +
      geom_node_text(
        aes(label = name), size = 3, repel = TRUE
      ) +
      theme_graph()
    set.seed(123)
    pnet.graph %>%
      activate(nodes) %>%
      mutate(centrality = centrality_authority()) %>%
      ggraph(layout = "graphopt") +
      geom_edge_link(width = 1, colour = "lightgray") +
      geom_node_point(aes(size = centrality, colour = centrality)) +
      geom_node_text(aes(label = name), repel = TRUE)+
      scale_color_gradient(low = "yellow", high = "red")+
      theme_graph()

    xy<-pnet.graph %>%
      activate(nodes) %>%
      mutate(centrality = methodChoice())%>%
      arrange(desc(centrality))%>%
      as_tibble()
    return(as.data.frame(xy)[1:input$topVars,])
  })

  output$resTopVars = renderDataTable({if(!is.null(constructDepNet()))
    TopVars()}, server = TRUE)

  output$graphReady<-reactive(!is.null(TopVars()))
  outputOptions(output, 'graphReady', suspendWhenHidden = FALSE)

  clustVarsSelected<-reactiveVal()
  clustData<-reactiveVal()

  observeEvent(input$finalizeVars, { withBusyIndicatorServer("finalizeVars",{
    tmp<-TopVars()[,1]
    clustVarsSelected(tmp)
    tbdata<-req(imported$data())
    tbdata<-select_if(tbdata, is.numeric)
    clustData(tbdata[,clustVarsSelected()])
    output$txtClustVars <- renderPrint({
      req(clustVarsSelected())
      clustVarsSelected()
    })
    output$txtClustVars2 <- renderPrint({
      req(clustVarsSelected())
      clustVarsSelected()
    })
   output$ClustPlot <- renderPlot({
      if(!is.null(clustData())){

        suppressWarnings({
          fit <- Mclust(as.data.frame(clustData()))
        })

        output$mclustModelOut <- renderPrint({
          summary(fit)
        })

        fviz_mclust(fit, "classification", geom = "point")

      }
    }, height=500)

    rval_kclust<-reactive({
      req(input$ksolution)
      centers <- as.numeric(eval_tidy(input$ksolution))
      kmeans(clustData(), centers = centers)
    })

    output$kmeansClustPlot <- renderPlot({
      if(!is.null(clustData())){

        fviz_cluster(rval_kclust(), data = clustData())

      }
    }, height=600)

    output$elbowPlot<-renderPlot({
      if(!is.null(clustData())){
      wssplot(clustData())
      }
    })

  })
  })


  }


if (interactive())
  shinyApp(ui, server)
