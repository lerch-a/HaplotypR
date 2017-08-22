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
             # tabPanel("Home", homepage()), 
             tabPanel("Input", setinput()),
             navbarMenu("Views",
                        tabPanel("Sample List", viewSampleList()),
                        tabPanel("Marker List", viewMarkerList())
                        ),
             # tabPanel("Quality Control"),
             navbarMenu("Process Reads",
                        tabPanel("De-multiplex by Sample", demultiplexSample()),
                        tabPanel("De-multiplex by Marker", demultiplexMarker()),
                        tabPanel("Concatenate Paired Reads", concatreads())
                        ),
             tabPanel("Call Genotypes", callgenotype()),
             tabPanel("Call Haplotypes", callhaplotype())
             # tabPanel("Settings")
             # tabPanel("Help")
  )
)

# Define server logic 
server <- shinyServer(function(input, output, session) {
  source("serverHelper.R", local=TRUE)
  #volumes <- getVolumes()
  volumes <- c(Home='~')
  #runvolumes <- c(Home=getwd())
  # volumes <- reactive({
  #   if(change$config>0){
  #     vol <- config$rootVolumes
  #     browser()
  #     names(vol) <- names(config$rootVolumes)
  #     return(vol)
  #   }
  #   return(c(Home='~'))
  # })

  change <- reactiveValues(config=0)
  config <- list(fnReadR1=NULL, fnReadR2=NULL, fnBarcodeF=NULL, fnBarcodeR=NULL,
                           fnPrimer=NULL, rootVolumes=volumes, outputDir=NULL)


  ### Input Tab
  shinyFileSave(input, 'saveConfig', roots=volumes, session=session)
  observeEvent(input$saveConfig, {
    configLst <- list(fileR1=isolate(input$fileR1), fileR2=isolate(input$fileR2),
                      fileB1=isolate(input$fileB1), fileB2=isolate(input$fileB2),
                      fileP1=isolate(input$fileP1), fileSample=isolate(input$fileSample),
                      dirOut=isolate(input$dirOut))
    configLst$rootVolumes <- volumes
    saveRDS(configLst, file=as.character(parseSavePath(volumes, input$saveConfig)$datapath))
  })
  shinyFileChoose(input, 'loadConfig', roots=volumes, session=session)
  observeEvent(input$loadConfig, {
    configLst <- readRDS(as.character(parseFilePaths(volumes, input$loadConfig)$datapath))
    lapply(names(configLst), function(n) { config[[n]] <<- configLst[[n]] })
    change$config <- change$config+1
  })
  shinyFileChoose(input, 'fileR1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileR2', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileB2', roots=volumes, session=session)
  shinyFileChoose(input, 'fileP1', roots=volumes, session=session)
  shinyFileChoose(input, 'fileSample', roots=volumes, session=session)
  shinyDirChoose(input, 'dirOut', roots=volumes, session=session)
  shinyFileSave(input, 'createSampleFile', roots=volumes, session=session)
  # config$fnReadR1 <- observeEvent(input$fileR1, { return(isolate(input$fileR1)) })
  fileR1 <- reactive({
    if(change$config>0) return(paste(".", paste(config$fileR1$files[[1]], collapse = "/"), sep=""))
    if(is.null(input$fileR1)) return(NULL)
    return(paste(".", paste(input$fileR1$files[[1]], collapse = "/"), sep=""))
  })
  fileR2 <- reactive({
    if(change$config>0) return(paste(".", paste(config$fileR2$files[[1]], collapse = "/"), sep=""))
    if(is.null(input$fileR2)) return(NULL)
    return(paste(".", paste(input$fileR2$files[[1]], collapse = "/"), sep=""))
  })
  fileB1 <- reactive({
    if(change$config>0) return(paste(".", paste(config$fileB1$files[[1]], collapse = "/"), sep=""))
    if(is.null(input$fileB1)) return(NULL)
    return(paste(".", paste(input$fileB1$files[[1]], collapse = "/"), sep=""))
  })
  fileB2 <- reactive({
    if(change$config>0) return(paste(".", paste(config$fileB2$files[[1]], collapse = "/"), sep=""))
    if(is.null(input$fileB2)) return(NULL)
    return(paste(".", paste(input$fileB2$files[[1]], collapse = "/"), sep=""))
  })
  filePrimer <- reactive({
    if(change$config>0) return(paste(".", paste(config$fileP1$files[[1]], collapse = "/"), sep=""))
    if(is.null(input$fileP1)) return(NULL)
    return(paste(".", paste(input$fileP1$files[[1]], collapse = "/"), sep=""))
  })
  fileSample <- reactive({
    if(change$config>0) return(paste(".", paste(config$fileSample$files[[1]], collapse = "/"), sep=""))
    if(is.null(input$fileSample)) return(NULL)
    return(paste(".", paste(input$fileSample$files[[1]], collapse = "/"), sep=""))
  })
  dirOut <- reactive({
    if(change$config>0) return(paste(".", paste(config$dirOut$path, collapse = "/"), sep=""))
    if(is.null(input$dirOut)) return(NULL)
    return(paste(".", paste(input$dirOut$path, collapse = "/"), sep=""))
  })
  observeEvent(input$createSampleFile, {
    lst <- createSampleFile(input, output, session, volumes)
    write.table(lst, file=as.character(parseSavePath(volumes, input$createSampleFile)$datapath), sep="\t", row.names=F)
  })
  output$filepathsR1 <- renderText( ifelse(is.null(fileR1()), "Select forward read file", fileR1()))
  output$filepathsR2 <- renderText( ifelse(is.null(fileR2()), "Select reverse read file", fileR2()))
  output$filepathsB1 <- renderText( ifelse(is.null(fileB1()), "Select forward barcode file", fileB1()))
  output$filepathsB2 <- renderText( ifelse(is.null(fileB2()), "Select reverse barcode file", fileB2()))
  output$filepathsP1 <- renderText( ifelse(is.null(filePrimer()), "Select primer file", filePrimer()))
  output$filepathsSample <- renderText( ifelse(is.null(fileSample()), "Select sample file", fileSample()))
  output$dirpathOut <- renderText( ifelse(is.null(dirOut()), "Select output directory", dirOut()))

  baseOutDir <- reactive({
    if(is.null(dirOut())){
      showNotification("Output Directory missing", closeButton = T, type="error")
      return(NULL)
    }
    #out <- as.character(parseDirPath(volumes(), input$dirOut))
    out <- sub(pattern = "/\\./", "/", file.path(volumes, dirOut()))
    return(out)
  })


  ### Views
  sample <- reactive({
    fnSample <- fileSample()
    #fnSample <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/tmpShinyApp/sampleLst.txt"
    if(is.null(fnSample)){
      showNotification("Sample file missing", closeButton = T, type="error")
      return(NULL)
    }
    data <- read.delim(file.path(volumes, fileSample()), stringsAsFactors=F)
    return(data)
  })
  output$sampleTable <- renderTable({ sample() })

  primer <- reactive({
    fnPrimer <- filePrimer()
    if(is.null(fnPrimer)){
      showNotification("Primer sequence missing", closeButton = T, type="error")
      return(NULL) #return(data.frame()
    }
    fn <- sub("/\\./", "/", file.path(volumes, filePrimer()))
    if(!file.exists(fn)){
      showNotification("Primer sequence file not found", closeButton = T, type="error")
      return(NULL)
    }
    return(read.delim(fn, stringsAsFactors=F))
  })
  output$markerTable <- renderTable({ primer() })


  ### Demultiplex by sample Tab
  output$fnR1_deplex <- renderText(fileR1()) # ifelse(is.null(fileR1()), "Go to menu 'Input' and select forward read file", fileR1())
  output$fnR2_deplex <- renderText(fileR2())
  output$fnB1_deplex <- renderText(fileB1())
  output$fnB2_deplex <- renderText(fileB2())
  output$fnS_deplex <- renderText(fileSample())
  resDePlexS <- eventReactive(input$startDePlex, {
    runDemultiplex(input, output, session, volumes)
  })
  # TODO check this function
  outDePlexSample <- reactive({
    #data <- isolate(resDePlexS())
    fn <- file.path(baseOutDir(), "demultiplexSampleSummary.txt")
    if(!file.exists(fn)){
      showNotification("Summary file of de-multiplex by sample is missing. Did you run this step?", closeButton = T, type="error")
      return(NULL)
    }
    return(read.delim(fn))
  })
  output$table_dePlexS <- renderTable({ resDePlexS() })

  
  ### Demultiplex by Marker Tab
  output$primerTableDPM <- renderTable({ primer() })
  resDePlexM <- eventReactive(input$startTruncPrim, {
    runTruncatePrimer(input, output, session, volumes)
  })
  outDePlexMarker <- reactive({
    #data <- isolate(resDePlexM())
    fn <- file.path(baseOutDir(), "demultiplexMarkerSummary.txt")
    if(!file.exists(fn)){
      showNotification("Summary file of de-multiplex by marker is missing. Did you run this step?", closeButton = T, type="error")
      return(NULL)
    }
    return(read.delim(fn))
  })
  output$table_dePlexM <- renderTable(resDePlexM())

  
  ### Concatenate read Tab
  concatReads <- eventReactive(input$startConcatReads, {
    runConcatReads(input, output, session, volumes)
  })
  outConcatReads <- reactive({
    browser()
    # This code should not be used anymore
  #   out <- baseOutDir()
  #   if(is.null(out)) return(NULL)
  #   fn <- list.files(out, "processedReadSummary")
  #   browser()
  #   prj <- isolate(input$projectSelectedG)
  #   fn <- fn[grep(prj, fn)]
  #   if(!file.exists(file.path(out, fn))){
  #     showNotification("Summary file of concatenate paired reads is missing. Did you run this step?", closeButton = T, type="error")
  #     return(NULL)
  #   }
  #   return(read.delim(file.path(out, fn)))
  })
  output$table_concat <- renderTable(concatReads())

  ### project selection
  projectSelection <- reactive({
    out <- baseOutDir()
    if(is.null(out)) return(NULL)
    fn <- list.files(out, "processedReadSummary")
    prj <- do.call(rbind, strsplit(fn, "_bind"))
    prj <- paste("_bind", sub(".txt", "", prj[,2]), sep="")
    choiseProject <- as.list(prj)
    names(choiseProject) <- prj
    return(choiseProject)
  })
  
  ### marker selection
  markerSelection <- reactive({
    primerTable <- primer()
    #sampleTable <- outConcatReads()
    if(is.null(primerTable)) return(NULL)
    mID <- c("none", as.character(unique(primerTable$MarkerID)))
    choiseMarker <- as.list(mID)
    names(choiseMarker) <- mID
    return(choiseMarker)
  })

  
  ### Call Genotyp Tab
  output$uiProjectSelectionG <- renderUI({
    choiseProject <- projectSelection()
    if(is.null(choiseProject)){ 
      return(helpText("Summary file of concatenate paired reads not found. Did you run this step or is the ouput directory set?")) 
    }
    selectInput("projectSelectedG", label=NULL, choices=choiseProject)
  })
  
  output$uiMarkerSelectionG <- renderUI({
    choiseMarker <- markerSelection()
    if(is.null(choiseMarker)){ 
      return(helpText("Summary file of concatenate paired reads not found. Did you run this step or is the ouput directory set?")) 
    }
    selectInput("markerSelectedG", label=NULL, choices=choiseMarker)
  })
  callGen <- eventReactive(input$startCallGenotyp, {
    runCallGenotype(input, output, session, volumes)
  })
  outPotSNP <- reactive({
    out <- isolate(baseOutDir())
    if(is.null(out)) return(NULL)
    marker <- input$markerSelectedH
    proj <- input$projectSelectedH
    if(is.null(marker)) return(NULL)
    if(marker=="none") return(NULL)
    fnSNP <- file.path(out, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                    isolate(input$minMMrate)*100, 
                                    isolate(input$minOccGen), 
                                    marker,
                                    proj))
    if(!file.exists(fnSNP)){
      showNotification("Input 'Trim Parameter' is missing.", closeButton = T, type="error")
      return(NULL)
    }
    potSNP <- read.delim(fnSNP)
    return(potSNP)
  })
  outCallGen <- reactive({
    #isolate(callGen())
    out <- baseOutDir()
    if(is.null(out)) return(NULL)
    marker <- input$markerSelectedG
    if(is.null(marker)) return(NULL)
    if(marker=="none") return(NULL)
    fnErr <- file.path(out, sprintf("mismatchRate_%s.txt", marker))
    if(!file.exists(fnErr)){
      #showNotification("Input 'TODO' is missing. Did you run this step?", closeButton = T, type="error")
      return(NULL)
    }
    seqErr <- read.delim(fnErr)
    potSNP <- isolate(outPotSNP())
    return(list(seqErrors=seqErr, potentialSNP=potSNP))
  })
  #output$plotMisMatch <- renderPlot({ plotGenotype(input, output, session, volumes) })
  output$plotMisMatch <- renderPlot({ callGen() })
  
  
  ### Call Haplotyp Tab
  output$uiProjectSelectionH <- renderUI({
    choiseProject <- projectSelection()
    if(is.null(choiseProject)){ 
      return(helpText("Summary file of concatenate paired reads not found. Did you run this step or is the ouput directory set?")) 
    }
    selectInput("projectSelectedH", label=NULL, choices=choiseProject)
  })
  
  output$uiMarkerSelectionH <- renderUI({
    choiseMarker <- markerSelection()
    if(is.null(choiseMarker)){ 
      return(helpText("Summary file of concatenate paired reads is missing. Did you run this step?")) 
    }
    selectInput("markerSelectedH", label=NULL, choices=choiseMarker)
  })
  createHapTab <- eventReactive(input$startCallHaplotype, {
    runCreateHaplotypOverview(input, output, session, volumes)
  })
  #runCreateHaplotypOverview(input, output, session, volumes)
  #runCallHaplotype(input, output, session, volumes)

  # outHapOview <- reactive({
  #   browser()
  #   #overviewHap <- isolate(createHapTab())
  #   #if(!is.null(overviewHap)) return(overviewHap)
  #   out <- baseOutDir()
  #   if(is.null(out)) return(NULL)
  #   marker <- input$markerSelectedH
  #   if(is.null(marker)) return(NULL)
  #   if(marker=="none") return(NULL)
  #   fnHab <- file.path(out, sprintf("HaplotypeOverviewTable_%s.txt", marker))
  #   if(!file.exists(fnHab)){
  #     #showNotification("Input 'TODO' is missing. Did you run this step?", closeButton = T, type="error")
  #     overviewHap <- runCreateHaplotypOverview(input, output, session, volumes)
  #   }else{
  #     overviewHap <- read.delim(fnHab)
  #   }
  #   return(overviewHap)
  # })
  # outHapTab <- reactive({
  #   browser()
  #   #isolate(createHapTab())
  #   out <- baseOutDir()
  #   if(is.null(out)) return(NULL)
  #   marker <- input$markerSelectedH
  #   if(is.null(marker)) return(NULL)
  #   if(marker=="none") return(NULL)
  #   fnHab <- file.path(out, sprintf("finalHaplotypeTable_%s.txt", marker))
  #   if(!file.exists(fnHab)){
  #     #haplotypesSample <- runCallHaplotype(input, output, session, volumes)
  #     #showNotification("Input 'TODO' is missing. Did you run this step?", closeButton = T, type="error")
  #     return(NULL)
  #   }else{
  #     haplotypesSample <- read.delim(fnHab)
  #   }
  #   return(haplotypesSample)
  # })
  output$tableHaplotyps <- renderTable({ createHapTab() })

  # ### Views Haplotyp
  # output$overviewHap <- renderTable({ outHapOview() })
  # output$haplotypesSample <- renderTable({ outHapTab() })
  
  
})

# Run the application 
shinyApp(ui = ui, server = server)

