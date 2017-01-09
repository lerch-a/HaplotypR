#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# display.mode = "showcase"
# runApp('~/switchdrive/Repos/haplotypR/inst/shinyApp', display.mode = "showcase")
# options(shiny.trace=TRUE)getShinyOption("myOption")

library(shiny)
library(shinyFiles)
library(HaplotypR)
library(Biostrings)
library(ShortRead)
source("guiHelper.R")


# Define UI for application 
ui <- shinyUI(
  navbarPage("HaplotypR", theme='haplotypRstyle.css',
             tabPanel("Home", homepage()),
             tabPanel("Input", setinput()),
             navbarMenu("Views",
                        tabPanel("Sample List", viewSampleList()),
                        tabPanel("Marker List", viewMarkerList())
                        ),
             #tabPanel("Project Information"),
             #tabPanel("Quality Control"),
             navbarMenu("Process Reads",
                        tabPanel("De-multiplex by Sample", demultiplexSample()),
                        tabPanel("De-multiplex by Marker", demultiplexMarker()),
                        tabPanel("Concatenate Paired Reads", concatreads())
                        ),
             tabPanel("Call Genotypes"),
             tabPanel("Call Haplotypes"),
             #tabPanel("Settings")
             tabPanel("Help")
  )
)

# Define server logic 
server <- shinyServer(function(input, output, session) {
  source("serverHelper.R", local=TRUE)
  volumes <- c(Home='/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/') # getVolumes()
  
  change <- reactiveValues(config=0)
  config <- list(cTest="empty", fnReadsR1=NULL, fnReadsR2=NULL, fnBarcodeF=NULL, fnBarcodeR=NULL, 
                           linkerF=NULL, linkerR=NULL, fnPrimer=NULL, rootVolumes=volumes, outputDir=NULL)
  
  ## Set Input Tab
  # output$inputDyn <- renderUI({
  #   textInput(inputId = "dynamic", label ="Dynlab" ,value="3")
  # })
  # observeEvent(input$iTest, { 
  #   config$cTest <- isolate(input$iTest)
  #   browser()
  #   change$config <- change$config+1
  #   })
  # output$oTest <- renderText({if(!is.null(change$config))
  #   return(config$cTest) })
  shinyFileSave(input, 'saveConfig', roots=volumes, session=session)
  observeEvent(input$saveConfig, { 
    configVar <- c("fileR1", "fileR2", "fileB1", "fileB2", "linkerF", "linkerR", "fileP1", "dirOut")
    configLst <- lapply(configVar, function(n) { isolate(input[[n]]) })
    names(configLst) <- configVar
    configVar$rootVolumes <- volumes
    saveRDS(configLst, file=as.character(parseSavePath(volumes, input$saveConfig)$datapath))
  })
  shinyFileChoose(input, 'loadConfig', roots=volumes, session=session)
  observeEvent(input$loadConfig, {
    configLst <- readRDS(as.character(parseFilePaths(volumes, input$loadConfig)$datapath))
    lapply(names(configLst), function(n) { input[[n]] <- configLst[[n]] })
  })
  shinyFileChoose(input, 'fileR1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileR2', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB2', roots=volumes, session=session) 
  shinyFileChoose(input, 'fileP1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileSample', roots=volumes, session=session)
  shinyDirChoose(input, 'dirOut', roots=volumes, session=session)
  # config$fnReadR1 <- renderText({browser()
  #   paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep="")})
  # config$fnReadR2 <- observeEvent(parseFilePaths(volumes, input$fileR2)$datapath)
  # config$fnBarcodeF <- observeEvent(parseFilePaths(volumes, input$fileB1)$datapath)
  # config$fnBarcodeR <- observeEvent(parseFilePaths(volumes, input$fileB2)$datapath)
  # config$outputDir <- observeEvent(parseDirPath(volumes, input$dirOut))
  # config$linkerF <- input$linkerF
  # config$linkerR <- input$linkerF
  # config$fnPrimer <- observeEvent(parseFilePaths(volumes, input$fileP1)$datapath)
  output$filepathsR1 <- renderText({ ifelse(is.null(input$fileR1), "Select forward read file", 
                                            paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep=""))})  
  output$filepathsR2 <- renderText({ ifelse(is.null(input$fileR2), "Select reverse read file", 
                                            paste(".", paste(input$fileR2$files[[1]], collapse = "/"), sep=""))})
  output$filepathsB1 <- renderText({ ifelse(is.null(input$fileB1), "Select forward barcode file", 
                                            paste(".", paste(input$fileB1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsB2 <- renderText({ ifelse(is.null(input$fileB2), "Select reverse barcode file", 
                                            paste(".", paste(input$fileB2$files[[1]], collapse = "/"), sep=""))})
  output$filepathsP1 <- renderText({ ifelse(is.null(input$fileP1), "Select primer file", 
                                            paste(".", paste(input$fileP1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsSample <- renderText({ ifelse(is.null(input$fileSample), "Select sample file", 
                                            paste(".", paste(input$fileSample$files[[1]], collapse = "/"), sep=""))})
  output$linkerF <- renderText(input$linkerF)
  output$linkerR <- renderText(input$linkerR)
  output$dirpathOut <- renderText({ ifelse(is.null(input$dirOut), "Select output directory",
                                                 paste(".", paste(input$dirOut$path, collapse = "/"), sep=""))})
  
  shinyFileSave(input, 'createSampleFile', roots=volumes, session=session)
  observeEvent(input$createSampleFile, { 
    lst <- createSampleFile(input, output, session, volumes)
    write.table(lst, file=as.character(parseSavePath(volumes, input$createSampleFile)$datapath), sep="\t", row.names=F)
  })
  
  ###Views
  sample <- reactive({
    fnSample <- input$fileSample
    if(is.null(fnSample)) return(data.frame())
    data <- read.delim(as.character(parseFilePaths(volumes, input$fileSample)$datapath))
    return(data)
  })
  output$sampleTable <- renderTable({ sample() })
  primer <- reactive({
    fnPrimer <- input$fileP1
    if(is.null(fnPrimer)) return(data.frame())
    data <- read.delim(as.character(parseFilePaths(volumes, input$fileP1)$datapath))
    return(data)
  })
  output$markerTable <- renderTable({ primer() })
  
  
  ### Demultiplex by sample Tab 
  output$fnR1_deplex <- renderText(paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep=""))
  output$fnR2_deplex <- renderText(paste(".", paste(input$fileR2$files[[1]], collapse = "/"), sep=""))
  output$fnB1_deplex <- renderText(paste(".", paste(input$fileB1$files[[1]], collapse = "/"), sep=""))
  output$fnB2_deplex <- renderText(paste(".", paste(input$fileB2$files[[1]], collapse = "/"), sep=""))
  output$fnS_deplex <- renderText(paste(".", paste(input$fileSample$files[[1]], collapse = "/"), sep=""))
  shinyDirChoose(input, 'dirOutDePlex', roots=volumes, session=session)
  output$dirpathOutDePlex <- renderText({ ifelse(is.null(input$dirOutDePlex), "Select output directory",
                                                  paste(".", paste(input$dirOutDePlex$path, collapse = "/"), sep=""))})
  # Action
  resDePlex <- eventReactive(input$startDePlex, {
    runDemultiplex(input, output, session, volumes)
  })
  outDePlexSample <- reactive({
    data <- isolate(resDePlex())
    if(!is.null(outDePlexSample)) return(data)
    data <- read.delim(file.path(out, "demultiplexSampleSummary.txt"))
    return(data)
  })
  output$table_dePlexS <- renderTable({ resDePlex() })
  
  ### Demultiplex by Marker Tab 
  # shinyDirChoose(input, 'dirOutTruncPrim', roots=volumes, session=session)  
  output$dirOutTruncPrim <- renderText({ ifelse(is.null(input$dirTruncPrim), "Select output directory",
                                                  paste(".", paste(input$dirTruncPrim$path, collapse = "/"), sep=""))})
  output$primerTableDPM <- renderTable({ primer() })
  # Action
  truncPrim <- eventReactive(input$startTruncPrim, {
    runTruncatePrimer(input, output, session, volumes)
  })
  output$table_dePlexM <- renderTable(truncPrim())

  ### Demultiplex by Marker Tab 
  # Action
  concatReads <- eventReactive(input$startConcatReads, {
    runConcatReads(input, output, session, volumes)
  })
  output$table_concat <- renderTable(concatReads())
  
  
})

# Run the application 
shinyApp(ui = ui, server = server)

