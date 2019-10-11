
demultiplexReads <- function(fastqFileFwd, fastqFileRev, barcodeFileFwd, barcodeFileRev, outputDir, adapterFwd=NULL, adapterRev=NULL, max.mismatch=0, with.indels=F, progressReport=message){
  require(Biostrings)
  require(ShortRead)
  
  # Read barcodes
  barcodesFwd <- Biostrings::readDNAStringSet(barcodeFileFwd)
  barcodesFwdLength <- unique(width(barcodesFwd))
  if(length(barcodesFwdLength)>1)
    stop("Barcodes length must have equal length.")
  barcodesRev <- Biostrings::readDNAStringSet(barcodeFileRev)
  barcodesRevLength <- unique(width(barcodesRev))
  if(length(barcodesRevLength)>1)
    stop("Barcodes length must have equal length.")
  
  # check existing output files
  of <- list.files(path=outputDir, pattern=names(barcodesFwd))
  if(length(of)>0)
    stop("Output directory must be empty. Found existing files: ", paste(of, collapse=", "))
  
  # check and set progress report function
  if(!is.function(progressReport))
    progressReport <- message
  
  # Start demultiplexing
  f1 <- FastqStreamer(fastqFileFwd)
  f2 <- FastqStreamer(fastqFileRev)
  mode <- "w"
  countReads <- 0
  sumDemultiplex <- character(0)
  msg <- paste("Processing file", basename(fastqFileFwd), "and", basename(fastqFileRev), ":")
  progressReport(detail=msg)
  while(length(sr1 <- yield(f1)) > 0){
    sr2 <- yield(f2)
    countReads <- countReads + length(sr1)
    
    if(is.null(adapterFwd)){
      sr1_trim <- narrow(sr1, start=1, width=barcodesFwdLength)
      sr2_trim <- narrow(sr2, start=1, width=barcodesRevLength)
    } else {
      
      #      matchPDict( whichPDict
      #       index <- vmatchPattern(adapterFwd, sread(sr1), max.mismatch=max.mismatch, with.indels=with.indels)
      #
      #       fPrimEnd <- vmatchPattern(primerFw, sread(sr1), max.mismatch=max.mismatch, with.indels=with.indels)
      #       fPrimEnd <- unlist(lapply(fPrimEnd, end))
      #index <- unlist(lapply(index, length))==1
      #sr1 <- sr1[index]
      #sr2 <- sr2[index]
      #sr1_trim <- narrow(sr1, start=fPrimEnd, width=ifelse(width(sr1)-fPrimEnd>=read1Length, read1Length, NA))
    }
    
    idxFwd <- Biostrings::match(sread(sr1_trim), barcodesFwd)
    idxRev <- Biostrings::match(sread(sr2_trim), barcodesRev)
    
    idx <- paste(names(barcodesFwd)[idxFwd], names(barcodesRev)[idxRev], sep="-")
    sumDemultiplex <- c(sumDemultiplex, idx)
    sr1_lst <- split(sr1, idx)
    sr2_lst <- split(sr2, idx)
    
    lapply(seq_along(sr1_lst), function(i){
      outFile <- file.path(outputDir, names(sr1_lst[i]))
      writeFastq(sr1_lst[[i]], file=paste(outFile, "_R1.fastq.gz", sep=""), mode=mode, compress=T)
      writeFastq(sr2_lst[[i]], file=paste(outFile, "_R2.fastq.gz", sep=""), mode=mode, compress=T)
    })
    
    mode <- "a"
    progressReport(detail=paste(msg, countReads, "reads done ..."))
  }
  progressReport(detail=paste(msg, "finished, total" , countReads, "demultiplexed reads."))
  close(f1)
  close(f2)
  
  # Format output
  tab <- as.data.frame(table(sumDemultiplex), stringsAsFactors=F)
  colnames(tab) <- c("BarcodePair", "ReadNumbers")
  
  outfiles <- list.files(outputDir)
  names(outfiles) <- sub("_R..fastq.gz$", "", outfiles)
  outFileR1 <- outfiles[grep("_R1.fastq.gz$", outfiles)][tab$BarcodePair]
  outFileR2 <- outfiles[grep("_R2.fastq.gz$", outfiles)][tab$BarcodePair]
  tab$FileR1 <- file.path(outputDir, outFileR1)
  tab$FileR2 <- file.path(outputDir, outFileR2)
  
  return(invisible(tab))
}


removePrimer <- function(fastqFileR1, fastqFileR2, outputFile, primerFwd, primerRev,
                         max.mismatch=0, with.indels=F, outputPrimerSequence=F, progressReport=message){
  
  # check and set progress report function
  if(!is.function(progressReport))
    progressReport <- message
  msg <- paste("Processing file", basename(fastqFileR1), "and", basename(fastqFileR2))
  progressReport(detail=msg)
  
  if(length(fastqFileR1) != length(fastqFileR2))
    stop("Vector length of fastqFileR1 and fastqFileR2 not identical.")
  
  f1 <- FastqStreamer(fastqFileR1)
  f2 <- FastqStreamer(fastqFileR2)
  mode <- "w"
  totalReads <- 0
  filteredReads <- 0
  while(length(sr1 <- yield(f1)) > 0){
    sr2 <- yield(f2)
    totalReads <- totalReads+length(sr1)
    
    index <- which(vcountPattern(primerFwd, sread(sr1), max.mismatch=max.mismatch, with.indels=with.indels) == 1)
    sr1 <- sr1[index]
    sr2 <- sr2[index]
    
    index <- which(vcountPattern(primerRev, sread(sr2), max.mismatch=max.mismatch, with.indels=with.indels) == 1)
    sr2 <- sr2[index]
    sr1 <- sr1[index]
    
    if(length(sr1)>1){
      fPrimEnd <- vmatchPattern(primerFwd, sread(sr1), max.mismatch=max.mismatch, with.indels=with.indels)
      fPrimEnd <- unlist(endIndex(fPrimEnd))
      sr1_trim <- narrow(sr1, start=fPrimEnd+1, width=NA)
      start <- fPrimEnd-nchar(primerFwd)+1
      sr1_prim <- narrow(sr1, start=ifelse(start<1, 1, start), end=fPrimEnd)
      
      rPrimEnd <- vmatchPattern(primerRev, sread(sr2), max.mismatch=max.mismatch, with.indels=with.indels)
      rPrimEnd <- unlist(endIndex(rPrimEnd))
      sr2_trim <- narrow(sr2, start=rPrimEnd+1, width=NA)
      start <-rPrimEnd-nchar(primerRev)+1
      sr2_prim <- narrow(sr2, start=ifelse(start<1, 1, start), end=rPrimEnd)
      
      # ensure no zero length reads present
      index <- which((width(sr1_trim) > 0) & (width(sr2_trim) > 0))
      sr1_trim <- sr1_trim[index]
      sr2_trim <- sr2_trim[index]
      sr1_prim <- sr1_prim[index]
      sr2_prim <- sr2_prim[index]
      
      if(outputPrimerSequence){
        writeFastq(sr1_prim,
                   file=paste(outputFile, "_primerF.fastq.gz", sep=""), mode=mode, compress=T)
        writeFastq(sr2_prim,
                   file=paste(outputFile, "_primerR.fastq.gz", sep=""), mode=mode, compress=T)
      }
      writeFastq(sr1_trim, file=paste(outputFile, "_F.fastq.gz", sep=""), mode=mode, compress=T)
      writeFastq(sr2_trim, file=paste(outputFile, "_R.fastq.gz", sep=""), mode=mode, compress=T)
      
      mode <- "a"
      filteredReads <- filteredReads + length(sr1_trim)
    }
  }
  close(f1)
  close(f2)
  suppressWarnings(rm(sr1, sr2, sr1_trim, sr2_trim))
  gc()
  gc()
  if(filteredReads==0)
    return(data.frame(numReadIn=totalReads, numReadOut=filteredReads, FileR1=NA_character_, FileR2=NA_character_, stringsAsFactors=F))
  else
    return(data.frame(numReadIn=totalReads, numReadOut=filteredReads,
                      FileR1=paste(outputFile, "_F.fastq.gz", sep=""), FileR2=paste(outputFile, "_R.fastq.gz", sep=""),
                      stringsAsFactors=F))
}

demultiplexByMarker <- function(sampleTable, markerTable, trimFilenameExt="R1\\.fastq.gz", progressReport=message){
  resM <- do.call(rbind, lapply(1:dim(markerTable)[1], function(j){
    mID <- as.character(markerTable[j, "MarkerID"])
    adapterF <- as.character(markerTable[j, "Forward"])
    adapterR <- as.character(markerTable[j, "Reverse"])
    
    # check and set progress report function
    if(!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing marker", mID, sep=" ")
    progressReport(detail=msg)
    
    # process each file demultiplexed by sample
    res <- do.call(rbind.data.frame, lapply(seq_along(sampleTable$FileR1), function(i){
      outputFile <- file.path(outDeplexMarker, sub(trimFilenameExt, mID, basename(as.character(sampleTable$FileR1)[i])))
      # demultiplex by marker and truncate primer sequence
      removePrimer(as.character(sampleTable$FileR1)[i], as.character(sampleTable$FileR2)[i], outputFile,
                   adapterF, adapterR, max.mismatch=2, with.indels=F)
    }))
    cbind.data.frame(BarcodePair=as.character(sampleTable$BarcodePair), MarkerID=mID, res, stringsAsFactors=F)
  }))
  resM <- merge.data.frame(sampleTable[,c("SampleID", "SampleName", "BarcodePair")], resM, by="BarcodePair")
  return(resM)
}

mergeAmpliconReads <- function(fastqFileR1, fastqFileR2, outputDir, mergePrefix="_merge", progressReport=message){
  
  if(length(fastqFileR1) != length(fastqFileR2))
    stop("Vector length of fastqFileR1 and fastqFileR2 not identical.")
  
  
  tab <- lapply(seq_along(fastqFileR1), function(i){
    
    # check and set progress report function
    if(!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing file", basename(fastqFileR1[i]), "and", basename(fastqFileR2[i]), "...")
    progressReport(detail=msg, value=i)
    
    outputFile <- file.path(outputDir, sub("_F\\.fastq.gz", "", basename(fastqFileR1[i])))
    outputFile <- paste(outputFile, mergePrefix, ".fastq.gz", sep="")
    args <- paste("--fastq_mergepairs", fastqFileR1[i], "--reverse", fastqFileR2[i],
                  "--fastqout", outputFile, "--fastq_truncqual", 1, "--fastq_maxns", 0)
    Rvsearch:::.vsearchBin(args=args)
    return(data.frame(numRead=NA_integer_, ReadFile=outputFile, stringsAsFactors=F))
  })
  tab <- do.call(rbind.data.frame, tab)
  return(tab)
}


bindAmpliconReads <- function(fastqFileR1, fastqFileR2, outputDir, read1Length=NULL, read2Length=read1Length, mergePrefix="_bind", progressReport=message){
  
  if(length(fastqFileR1) != length(fastqFileR2))
    stop("Vector length of fastqFileR1 and fastqFileR2 not identical.")
  
  
  tab <- lapply(seq_along(fastqFileR1), function(i){
    
    # check and set progress report function
    if(!is.function(progressReport))
      progressReport <- message
    msg <- paste("Processing file", basename(fastqFileR1[i]), "and", basename(fastqFileR2[i]), "...")
    progressReport(detail=msg, value=i)
    
    outputFile <- file.path(outputDir, sub("_F\\.fastq.gz", "", basename(fastqFileR1[i])))
    outputFile <- paste(outputFile, mergePrefix, read1Length, "_", read2Length, ".fastq.gz", sep="")
    
    f1 <- FastqStreamer(fastqFileR1[i])
    f2 <- FastqStreamer(fastqFileR2[i])
    mode <- "w"
    numReads <- 0
    
    while(length(sr1 <- yield(f1)) > 0){
      sr2 <- yield(f2)
      numReads <- numReads+length(sr1)
      
      if(!is.null(read1Length))
        sr1 <- narrow(sr1, start=1, width=ifelse(width(sr1)>=read1Length, read1Length, NA))
      
      if(!is.null(read2Length))
        sr2 <- reverseComplement(narrow(sr2, start=1, width=ifelse(width(sr2)>=read2Length, read2Length, NA)))
      else
        sr2 <- reverseComplement(sr2)
      
      
      writeFastq(ShortReadQ(sread=xscat(sread(sr1), sread(sr2)),
                            qual=xscat(quality(quality(sr1)), quality(quality(sr2))),
                            id=id(sr1)),
                 file=outputFile, mode=mode, compress=T)
      mode <- "a"
    }
    close(f1)
    close(f2)
    rm(sr1, sr2)
    gc()
    gc()
    if(numReads==0)
      return(NULL)
    else
      return(data.frame(numRead=numReads, ReadFile=outputFile, stringsAsFactors=F))
  })
  tab <- do.call(rbind.data.frame, tab)
  return(tab)
}
