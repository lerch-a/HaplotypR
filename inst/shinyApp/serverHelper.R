
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
                   progressReport=updateProgress)
  write.table(res, file.path(dirOut, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)
  
  return(res) 
}

runTruncatePrimer <- function(input, output, session, volumes){

  numMM <- isolate(input$numMisMatch)
  wIndel <- isolate(input$withIndel)
  # if(length(isolate(parseDirPath(volumes, input$dirOut)))){
  #     showNotification("Output Directory missing", closeButton = T, type="error")
  #     return(NULL)
  # }
  # out <- isolate(parseDirPath(volumes, input$dirOut))
  out <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/tmpShinyApp/"
  outDePlex <- file.path(out, "demultiplexSampleSummary.txt")
  if(file.exists(outDePlex)){
    dePlexFiles <- read.delim(file.path(out, "demultiplexSampleSummary.txt"))
  }else{
    showNotification("Input missing", closeButton = T, type="error")
    return(NULL)
  }
  
  dirOut <- file.path(out, "demultiplexByMarker")
  if(!file.exists(dirOut)){
    dir.create(dirOut)
  }else{
    showNotification("Output dir exists", closeButton = T, type="error")
    return(NULL)
  }

  markerTab <- read.delim("/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/markerSeqLst.txt")
  # cspprimerF <-    "AAATGACCCAAACCGAAATGT"
  # cspprimerR <-    "GGAACAAGAAGGATAATACCA"
  # linkerF <- DNAString("GTGACCTATGAACTCAGGAGTC")
  # linkerR <- DNAString("CTGAGACTTGCACATCGCAGC")
  # adapterF <- paste(linkerF, cspprimerF, sep="")
  # adapterR <- paste(linkerR, cspprimerR, sep="")
  
  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=length(dePlexFiles$FileR1))
  on.exit(progress$close())
  updateProgress <- function(detail=NULL){
    progress$set(detail=detail)
  }
  
  # Run
  res <- lapply(1:dim(markerTab)[1], function(j){
    mID <- as.character(markerTab[j, "MarkerID"])
    adapterF <- as.character(markerTab[j, "Forward"])
    adapterR <- as.character(markerTab[j, "Reverse"])
    dir.create(file.path(dirOut, mID), showWarnings=F)

    res <- lapply(seq_along(dePlexFiles$FileR1), function(i){
      progress$set(message=sprintf('De-Multiplexing of marker %s in progress. This may take a while...', mID), detail=NULL, value=i)
      outputFile <- file.path(dirOut, mID, sub("R1\\.fastq.gz", mID, basename(as.character(dePlexFiles$FileR1)[i])))
      removePrimer(as.character(dePlexFiles$FileR1)[i], as.character(dePlexFiles$FileR2)[i], outputFile, 
                   adapterF, adapterR, max.mismatch=numMM, with.indels=wIndel, progressReport=updateProgress)
    })
    cbind(BarcodePair=dePlexFiles$BarcodePair, MarkerID=mID, do.call(rbind, res))
  })
  res <- do.call(rbind, res)
  write.table(res, file.path(dirOut, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
  
  return(res)
}
