
runDemultiplex <- function(input, output, session, volumes){
  
  # read and check input values
  # if(is.null(input$fileR1)){
  #   showNotification("Input missing R1", closeButton = T, type="error")
  #   return(NULL)
  # }
  # if(is.null(input$fileR2)){
  #   showNotification("Input missing R2", closeButton = T, type="error")
  #   return(NULL)
  # }
  # if(is.null(input$fileB1)){
  #   showNotification("Input missing B1", closeButton = T, type="error")
  #   return(NULL)
  # }
  # if(is.null(input$fileB2)){
  #   showNotification("Input missing B2", closeButton = T, type="error")
  #   return(NULL)
  # }
  # if(is.null(input$dirOutDePlex)){
  #   showNotification("Output Directory missing", closeButton = T, type="error")
  #   return(NULL)
  # }
  fnR1 <- isolate(parseFilePaths(volumes, input$fileR1)$datapath)
  fnR2 <- isolate(parseFilePaths(volumes, input$fileR2)$datapath)
  fnB1 <- isolate(parseFilePaths(volumes, input$fileB1)$datapath)
  fnB2 <- isolate(parseFilePaths(volumes, input$fileB2)$datapath)
  dirOut <- isolate(parseDirPath(volumes, input$dirOutDePlex))
  numMM <- isolate(input$numMisMatch)
  wIndel <- isolate(input$withIndel)
  adapF <- isolate(input$AdapterF)
  adapR <- isolate(input$AdapterR)
  if(nchar(adapF)==0) 
    adapF <- NULL
  if(nchar(adapR)==0) 
    adapR <- NULL

  fnR1 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_raw_data/Fwd_01-Rev_01_R1.fastq"
  fnR2 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_raw_data/Fwd_01-Rev_01_R2.fastq"
  fnB1 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/AmpliconMultiplexIndexFwd.fasta"
  fnB2 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/AmpliconMultiplexIndexRev.fasta"
  dirOut <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_raw_data_run2"

  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=1)
  on.exit(progress$close())
  progress$set(message='De-Multiplexing in progress. This may take a while...', value=1)
  updateProgress <- function(detail=NULL){
    progress$set(detail=detail)
  }
  
  # run demultiplexing
  res <- demultiplexReads(as.character(fnR1), as.character(fnR2), as.character(fnB1), as.character(fnB2), dirOut,
                   adapterFwd=NULL, adapterRev=NULL,
                   max.mismatch=numMM, with.indels=wIndel,
                   progressReport=updateProgress)
  
  return(res) # "De-Multiplexing finished"
}

runTruncatePrimer <- function(input, output, session, volumes){
  return("Run truncation")
}
