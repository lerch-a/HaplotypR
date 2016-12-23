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
library(haplotypR)

# Define UI for application 
ui <- shinyUI(
  navbarPage("haplotypR", theme='haplotypRstyle.css',
             ##
             tabPanel("Home"),
             #tabPanel("Project Information"),
             tabPanel("Quality Control"),
             ##De-multiplex Tab
             tabPanel("De-multiplex",
                      fluidPage(
                        titlePanel("De-multiplex Amplicon Sequence Run"),
                        hr(),
                        fluidRow(
                        column(12,
                               h5("Choose sequence run files (fastq format):"),
                               fluidRow(
                                 column(2, shinyFilesButton('fileR1', 'Browse...', 'Please select the forward read file', FALSE)),
                                 column(4, textOutput('filepathsR1')),
                                 column(2, shinyFilesButton('fileR2', 'Browse...', 'Please select the reverse read file', FALSE)),
                                 column(4, textOutput('filepathsR2'))
                               ),
                               h5("Choose barcode sequence files (fasta format):"),
                               fluidRow(
                                   column(2, shinyFilesButton('fileB1', 'Browse...', 'Please select the forward barcode file', FALSE)),
                                   column(4, textOutput('filepathsB1')),
                                   column(2, shinyFilesButton('fileB2', 'Browse...', 'Please select the reverse barcode file', FALSE)),
                                   column(4, textOutput('filepathsB2'))
                               ),
                               h5("Add adapter sequence:"),
                               fluidRow(
                                   column(6, textInput(inputId='AdapterF', label=NULL, value="", placeholder="Input forward adapter sequence")),
                                   column(6, textInput(inputId='AdapterR', label=NULL, value="", placeholder="Input reverse adapter sequence"))
                               ),
                               h5("Choose xxx:"),
                               fluidRow(
                                 column(2, shinyDirButton('dirOutDePlex', 'Browse...', 'Please select the output directory', FALSE)),
                                 column(4, textOutput("directorypath")),
                                 column(2, numericInput("numMisMatch", "Mis-match", 0, min=0, step=1)),
                                 column(4, checkboxInput("withIndel", "with-indel"),
                                        actionButton("startDePlex", "Start ..."),
                                        verbatimTextOutput("runDePlex"))
                               ),
                               h5("Choose marker sequence files (fasta format):"),
                               fluidRow(
                                 column(2, shinyFilesButton('fileB1', 'Browse...', 'Please select the forward marker file', FALSE)),
                                 column(4, helpText("Select forward marker adapter file")),
                                 column(2, shinyFilesButton('fileB2', 'Browse...', 'Please select the reverse marker file', FALSE)),
                                 column(4, helpText("Select reverse marker adapter file"))
                               )
                          )),
                        hr()
                      )
             ),
             tabPanel("Join paired read"),
             tabPanel("Call Genotypes"),
             tabPanel("Call Haplotypes"),
             tabPanel("Settings")
  )
)

# Define server logic 
server <- shinyServer(function(input, output, session) {
  # reactiveValues()
  #output$filepathsR1 <- "select file"
  volumes <- c(Home='/Users/anitalerch/CloudStation/') # getVolumes()
  #Buttons
  shinyFileChoose(input, 'fileR1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileR2', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB2', roots=volumes, session=session) 
  shinyDirChoose(input, 'dirOutDePlex', roots=volumes, session=session)  
  output$filepathsR1 <- renderText({ ifelse(is.null(input$fileR1), "Select forward read file", 
                                            paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsR2 <- renderText({ ifelse(is.null(input$fileR2), "Select reverse read file", 
                                            paste(".", paste(input$fileR2$files[[1]], collapse = "/"), sep=""))})
  output$filepathsB1 <- renderText({ ifelse(is.null(input$fileB1), "Select forward read file", 
                                            paste(".", paste(input$fileB1$files[[1]], collapse = "/"), sep=""))})
  output$filepathsB2 <- renderText({ ifelse(is.null(input$fileB2), "Select reverse read file", 
                                            paste(".", paste(input$fileB2$files[[1]], collapse = "/"), sep=""))})
  # demultiplexReads(fastqFileFwd, fastqFileRev, barcodeFileFwd, barcodeFileRev, 
  #                              outputDir, adapterFwd=NULL, adapterRev=NULL, max.mismatch=0, with.indels=F){
  output$directorypath <- renderPrint({parseDirPath(volumes, input$dirOutDePlex)})
  output$runDePlex <- eventReactive(input$startDePlex, { #observeEvent
    progress <- Progress$new(session, min=0, max=1)
    on.exit(progress$close())
    progress$set(message = 'Calculation in progress',
                 detail = 'This may take a while...', value=1)
    #progress$set(value = 1)
    Sys.sleep(1)
    isolate(parseDirPath(volumes, input$dirOutDePlex))
    withProgress(Sys.sleep(3), min = 0, max = 1, value =0.5,
               message = "Haallo", detail ="deetail")
  })
})

# Run the application 
shinyApp(ui = ui, server = server)

