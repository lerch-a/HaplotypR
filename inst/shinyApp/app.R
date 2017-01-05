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
             tabPanel("Views"),
             #tabPanel("Project Information"),
             #tabPanel("Quality Control"),
             navbarMenu("Process Reads",
                        tabPanel("De-multiplex Sample", demultiplex()),
                        tabPanel("De-multiplex Marker", demultiplexMarker()),
                        tabPanel("Join paired read")
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
  
  config <- reactiveValues(fnReadsR1=NULL, fnReadsR2=NULL, fnBarcodeF=NULL, fnBarcodeR=NULL, 
                           linkerF=NULL, linkerR=NULL, fnPrimerF=NULL, fnPrimerR=NULL, resDePlex=NULL,
                           rootVolumes=NULL, outputDir=NULL)
  
  ## Set Input Tab
  #saveConfig <- eventReactive(input$saveConfig, { saveRDS(input, "")})
  shinyFileSave(input, 'saveConfig', roots=volumes, session=session)
  observeEvent(input$saveConfig, {
    saveRDS(isolate(config), file=as.character(parseSavePath(volumes, input$saveConfig)$datapath))
  })
  shinyFileChoose(input, 'loadConfig', roots=volumes, session=session)
  observeEvent(input$loadConfig, { browser()
    config <- readRDS(as.character(parseFilePaths(volumes, input$loadConfig)$datapath))
  })
  shinyFileChoose(input, 'fileR1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileR2', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB2', roots=volumes, session=session) 
  shinyFileChoose(input, 'fileP1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileP2', roots=volumes, session=session) 
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
  output$linkerF <- renderText(input$linkerF)
  output$linkerR <- renderText(input$linkerR)
  shinyDirChoose(input, 'dirOut', roots=volumes, session=session)
  output$dirpathOut <- renderText({ ifelse(is.null(input$dirOut), "Select output directory",
                                                 paste(".", paste(input$dirOut$path, collapse = "/"), sep=""))})
  
  ## Demultiplex Tab 
  output$fnR1_deplex <- renderText(paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep=""))
  output$fnR2_deplex <- renderText(paste(".", paste(input$fileR2$files[[1]], collapse = "/"), sep=""))
  output$fnB1_deplex <- renderText(paste(".", paste(input$fileB1$files[[1]], collapse = "/"), sep=""))
  output$fnB2_deplex <- renderText(paste(".", paste(input$fileB2$files[[1]], collapse = "/"), sep=""))
  shinyDirChoose(input, 'dirOutDePlex', roots=volumes, session=session)
  output$dirpathOutDePlex <- renderText({ ifelse(is.null(input$dirOutDePlex), "Select output directory",
                                                  paste(".", paste(input$dirOutDePlex$path, collapse = "/"), sep=""))})
  # Action
  config$resDePlex <- eventReactive(input$startDePlex, {
    runDemultiplex(input, output, session, volumes)
  })
  output$sampleTable_deplex <- renderDataTable({ config$resDePlex() })
  
  ##
  # shinyDirChoose(input, 'dirOutTruncPrim', roots=volumes, session=session)  
  # output$dirpathOutTruncPrim <- renderText({ ifelse(is.null(input$dirTruncPrim), "Select output directory",
  #                                                paste(".", paste(input$dirTruncPrim$path, collapse = "/"), sep=""))})
  output$linkerF_rmprim <- renderText(input$linkerF)
  output$linkerR_rmprim <- renderText(input$linkerR)
  # Action
  truncPrim <- eventReactive(input$startTruncPrim, {
    runTruncatePrimer(input, output, session, volumes)
  })
  output$runTruncPrim <- renderText(truncPrim())


  # 
  # output$runTruncPrim <- eventReactive(input$startTruncPrim, { #observeEvent
  #   runTruncatePrimer(input, output, session, volumes)
  # })
  
})

# Run the application 
shinyApp(ui = ui, server = server)

