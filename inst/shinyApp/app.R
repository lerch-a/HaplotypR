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
             tabPanel("Input"),
             tabPanel("View"),
             #tabPanel("Project Information"),
             #tabPanel("Quality Control"),
             navbarMenu("Process Reads",
                        tabPanel("De-multiplex", demultiplex()),
                        tabPanel("Truncate Primer Sequence", rmprimer()),
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
  # reactiveValues()
  volumes <- c(Home='/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/') # getVolumes()
  
  #Buttons
  shinyFileChoose(input, 'fileR1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileR2', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB2', roots=volumes, session=session) 
  shinyDirChoose(input, 'dirOutDePlex', roots=volumes, session=session)
  shinyFileChoose(input, 'fileP1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileP2', roots=volumes, session=session) 
  shinyDirChoose(input, 'dirOutTruncPrim', roots=volumes, session=session)  
  output$filepathsR1 <- renderText({ ifelse(is.null(input$fileR1), "Select forward read file", 
                                            paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsR2 <- renderText({ ifelse(is.null(input$fileR2), "Select reverse read file", 
                                            paste(".", paste(input$fileR2$files[[1]], collapse = "/"), sep=""))})
  output$filepathsB1 <- renderText({ ifelse(is.null(input$fileB1), "Select forward barcode file", 
                                            paste(".", paste(input$fileB1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsB2 <- renderText({ ifelse(is.null(input$fileB2), "Select reverse barcode file", 
                                            paste(".", paste(input$fileB2$files[[1]], collapse = "/"), sep=""))})
  output$dirpathOutDePlex <- renderText({ ifelse(is.null(input$dirOutDePlex), "Select output directory",
                                                  paste(".", paste(input$dirOutDePlex$path, collapse = "/"), sep=""))})
  output$filepathsP1 <- renderText({ ifelse(is.null(input$fileP1), "Select forward primer file", 
                                            paste(".", paste(input$fileP1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsP2 <- renderText({ ifelse(is.null(input$fileP2), "Select reverse primer file", 
                                            paste(".", paste(input$fileP2$files[[1]], collapse = "/"), sep=""))})
  output$dirpathOutTruncPrim <- renderText({ ifelse(is.null(input$dirTruncPrim), "Select output directory",
                                                 paste(".", paste(input$dirTruncPrim$path, collapse = "/"), sep=""))})
  
  # Action demultiplex
  resDePlex <- eventReactive(input$startDePlex, {
    runDemultiplex(input, output, session, volumes)
  })
  #output$runDePlex <- renderText({ ifelse(!is.null, c("De-Multiplexing finished"), c("")) })
  output$sampleTable <- renderDataTable({
    resDePlex()
  })
  
  output$runTruncPrim <- eventReactive(input$startTruncPrim, { #observeEvent
    runTruncatePrimer(input, output, session, volumes)
  })
  
})

# Run the application 
shinyApp(ui = ui, server = server)

