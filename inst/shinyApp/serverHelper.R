
createSampleFile <- function(input, output, session, volumes){
  fnB1 <- as.character(isolate(parseFilePaths(volumes, input$fileB1)$datapath))
  fnB2 <- as.character(isolate(parseFilePaths(volumes, input$fileB2)$datapath))
  fnB1 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/AmpliconMultiplexIndexFwd.fasta"
  fnB2 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/AmpliconMultiplexIndexRev.fasta"
  if(is.null(fnB1)){
    showNotification("Forwarde barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnB2)){
    showNotification("Reverse barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  barcodeF <- readFasta(fnB1)
  barcodeR <- readFasta(fnB2)
  lst <- do.call(rbind, lapply(as.character(id(barcodeR)), function(br){
    data.frame(BarcodeID_F=as.character(id(barcodeF)), BarcodeID_R=br)
  }))
  lst <- as.data.frame(cbind(SampleID="", lst, AddInfo1="", AddInfo2=""))
  return(lst)
}

runDemultiplex <- function(input, output, session, volumes){
  
  fnR1 <- as.character(isolate(parseFilePaths(volumes, input$fileR1)$datapath))
  fnR2 <- as.character(isolate(parseFilePaths(volumes, input$fileR2)$datapath))
  fnB1 <- as.character(isolate(parseFilePaths(volumes, input$fileB1)$datapath))
  fnB2 <- as.character(isolate(parseFilePaths(volumes, input$fileB2)$datapath))
  wSampleName <- isolate(input$withSampleName)
  #dirOut <- as.character(isolate(parseDirPath(volumes, input$dirOutDePlex)))
  
  fnR1 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_raw_data/Fwd_01-Rev_01_R1.fastq"
  fnR2 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_raw_data/Fwd_01-Rev_01_R2.fastq"
  fnB1 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/AmpliconMultiplexIndexFwd.fasta"
  fnB2 <- "/Users/anitalerch/CloudStation/Documents/PhD/Data/files_duplicate/amplicon/_meta_data/AmpliconMultiplexIndexRev.fasta"
  
  # read and check input values
  if(is.null(fnR1)){
    showNotification("Forward seuquence file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnR2)){
    showNotification("Reverse seuquence file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnB1)){
    showNotification("Forwarde barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(is.null(fnB2)){
    showNotification("Reverse barcode file missing.", closeButton = T, type="error")
    return(NULL)
  }
  if(wSampleName){
    sampleTab <- isolate(sample())
    if(length(sampleTab)==0){
      showNotification("Sample file missing.", closeButton = T, type="error")
      return(NULL)
    }
  }
  

  
  # Demultiplexed sample files
  out <- isolate(parseDirPath(volumes, input$dirOut))
  if(is.null(out)){
    showNotification("Output Directory missing", closeButton = T, type="error")
    return(NULL)
  }
  # Output dir
  outDePlexSample <- file.path(out, "demultiplexBySample")
  if(!file.exists(outDePlexSample)){
    dir.create(outDePlexSample)
  }else{
    showNotification("Output dir exists", closeButton = T, type="error")
    return(NULL)
  }

  # Configure shiny progress
  progress <- Progress$new(session, min=0, max=1)
  on.exit(progress$close())
  progress$set(message='De-Multiplexing in progress. This may take a while...', value=1)
  updateProgress <- function(detail=NULL){
    progress$set(detail=detail)
  }
  
  # run demultiplexing
  res <- demultiplexReads(fnR1, fnR2, fnB1, fnB2, outDePlexSample, progressReport=updateProgress)

  if(wSampleName && length(sampleTab)!=0){
    res <- renameDemultiplexedFiles(sampleTab, res)
  }
  
  write.table(res, file.path(out, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)
  
  return(res) 
}


renameDemultiplexedFiles <- function(sampleTab, resTab){

  sampleTab$BarcodePair <- paste(sampleTab$BarcodeID_F, sampleTab$BarcodeID_R, sep="-")
  resTab <- merge.data.frame(sampleTab, resTab, by="BarcodePair", all.y=T)

  resTab <- lapply(1:dim(resTab)[1], function(i){
    if(!is.na(resTab$SampleID[i])){
      old <- resTab[i, c("FileR1","FileR2")]
      new <- sub(resTab$BarcodePair[i], resTab$SampleID[i], old)
      resTab[i, c("FileR1","FileR2")] <- new
      file.rename(as.character(old), as.character(new))
    }
    return(resTab[i,])
  })
  resTab <- do.call(rbind, resTab)
  
  return(resTab[, c("SampleID","BarcodePair","ReadNumbers","FileR1","FileR2")])
}


runTruncatePrimer <- function(input, output, session, volumes){
  
  # Options
  numMM <- isolate(input$numMisMatch)
  wIndel <- isolate(input$withIndel)
  
  # Primer sequence
  markerTab <- isolate(primer())
  if(length(markerTab)==0){
    showNotification("Primer sequence missing", closeButton = T, type="error")
    return(NULL)
  }
  
  # Demultiplexed sample files
  out <- isolate(parseDirPath(volumes, input$dirOut))
  if(is.null(out)){
      showNotification("Output Directory missing", closeButton = T, type="error")
      return(NULL)
  }
  
  outDePlexSample <- file.path(out, "demultiplexSampleSummary.txt")
  if(file.exists(outDePlexSample)){
    dePlexFiles <- read.delim(outDePlexSample)
  }else{
    showNotification("Input missing", closeButton = T, type="error")
    return(NULL)
  }
  
  # Output dir
  outDePlexMarker <- file.path(out, "demultiplexByMarker")
  if(!file.exists(outDePlexMarker)){
    dir.create(outDePlexMarker)
  }else{
    # showModal(modalDialog(
    #     helpText(sprintf("Output directory '%s' exists! Press OK to overwrite.", outDePlexMarker)),
    #     title="Warning!",
    #     size="l",
    #     footer = tagList(
    #       modalButton("Cancel"),
    #       actionButton("ok", "OK")
    #     )
    #   ))
    showNotification("Output dir exists", closeButton = T, type="error")
    return(NULL)
  }

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
    dir.create(file.path(outDePlexMarker, mID), showWarnings=F)

    res <- lapply(seq_along(dePlexFiles$FileR1), function(i){
      progress$set(message=sprintf('De-Multiplexing of marker %s in progress. This may take a while...', mID), detail=NULL, value=i)
      outputFile <- file.path(outDePlexMarker, mID, sub("R1\\.fastq.gz", mID, basename(as.character(dePlexFiles$FileR1)[i])))
      removePrimer(as.character(dePlexFiles$FileR1)[i], as.character(dePlexFiles$FileR2)[i], outputFile, 
                   adapterF, adapterR, max.mismatch=numMM, with.indels=wIndel, progressReport=updateProgress)
    })
    cbind(BarcodePair=as.character(dePlexFiles$BarcodePair), MarkerID=mID, do.call(rbind, res))
  })
  res <- do.call(rbind, res)
  write.table(res, file.path(out, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
  
  return(res)
}
