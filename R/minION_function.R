
mergeMinIONfiles <- function(inDir, outDir="dePlexSample", sampleTab){
  rownames(sampleTab) <- sampleTab$BarcodePair
  dirLst <- file.path(inDir, sampleTab$BarcodePair)
  tab <- do.call(rbind, lapply(dirLst, function(dir){
    prefix <- sampleTab[basename(dir),"SampleID"]
    fnLst <- list.files(dir, "fastq", full.names=T)
    reads <- unlist(lapply(fnLst, readLines))
    outfilename <- file.path(outDir, paste(prefix, "_R1.fastq", sep = ""))
    writeLines(text=reads, con=outfilename)
    return(data.frame(BarcodePair=basename(dir), FileR1=outfilename))
  }))
  tab$SampleID <- sampleTab[tab$BarcodePair, "SampleID"]
  tab$SampleName <- sampleTab[tab$BarcodePair, "SampleName"]
  tab$FileR2 <- NA
  return(tab)
}


demultiplexByMarkerMinION <- function (sampleTable, markerTable, outputDir, trimFilenameExt = "R1\\.fastq$", 
          progressReport = message, ...) 
{
  stopifnot(is.data.frame(sampleTable), nrow(sampleTable) > 
              0, all(c("FileR1", "SampleID", "SampleName", 
                       "BarcodePair") %in% colnames(sampleTable)), all(file.exists(c(sampleTable$FileR1))), is.data.frame(markerTable), nrow(markerTable) > 
              0, all(c("MarkerID", "Forward", "Reverse") %in% colnames(markerTable)), 
            is.character(outputDir), length(outputDir) == 1, file.exists(outputDir))
  resM <- do.call(rbind, lapply(1:dim(markerTable)[1], function(j) {
    mID <- as.character(markerTable[j, "MarkerID"])
    adapterF <- as.character(markerTable[j, "Forward"])
    adapterR <- as.character(markerTable[j, "Reverse"])
    if (!is.function(progressReport)) 
      progressReport <- message
    msg <- paste("Processing marker", mID, sep = " ")
    progressReport(detail = msg)
    res <- do.call(rbind.data.frame, lapply(seq_along(sampleTable$FileR1), 
                                            function(i) {
                                              outputFile <- file.path(outputDir, sub(trimFilenameExt, 
                                                                                     mID, basename(as.character(sampleTable$FileR1)[i])))
                                              removePrimerMinION(as.character(sampleTable$FileR1)[i], outputFile, 
                                                           adapterF, adapterR, with.indels = F, ...)
                                            }))
    cbind.data.frame(SampleID = as.character(sampleTable$SampleID), 
                     BarcodePair = as.character(sampleTable$BarcodePair), 
                     MarkerID = mID, res, stringsAsFactors = F)
  }))
  resM <- merge.data.frame(sampleTable[, c("SampleID", "SampleName", 
                                           "BarcodePair")], resM, by = c("SampleID", "BarcodePair"))
  return(resM)
}




removePrimerMinION <- function (fastqFileR1, outputFile, primerFwd, primerRev, 
          max.mismatch = 0, with.indels = F, outputPrimerSequence = F, 
          progressReport = message) 
{
  if (!is.function(progressReport)) 
    progressReport <- message
  msg <- paste("Processing file", basename(fastqFileR1))
  progressReport(detail = msg)
  f1 <- FastqStreamer(fastqFileR1)

  primerRevRC <- as.character(reverseComplement(DNAString(primerRev)))
  primerFwdRC <- as.character(reverseComplement(DNAString(primerFwd)))

  mode <- "w"
  totalReads <- 0
  filteredReads <- 0
  while (length(sr <- yield(f1)) > 0) {

    totalReads <- totalReads + length(sr)
    index <- vcountPattern(primerFwd, sread(sr), max.mismatch = max.mismatch, 
                           with.indels = with.indels)
    srF <- sr[index==1]
    sr <- sr[index!=1]
    index <- vcountPattern(primerRevRC, sread(srF), max.mismatch = max.mismatch, 
                           with.indels = with.indels)
    srF <- srF[index==1]
    #
    index <- vcountPattern(primerRev, sread(sr), max.mismatch = max.mismatch, 
                           with.indels = with.indels)
    srR <- sr[index==1]
    index <- vcountPattern(primerFwdRC, sread(srR), max.mismatch = max.mismatch, 
                           with.indels = with.indels)
    srR <- srR[index==1]
    srR <- reverseComplement(srR)
    sr <- append(srF,srR)
    writeFastq(sr, file = paste(outputFile, "_untrimF.fastq.gz", sep = ""), mode = mode, 
               compress = T)
    rm(srF, srR)
    if (length(sr) > 0) {
      fPrimCoord <- vmatchPattern(primerFwd, sread(sr),
                                max.mismatch = max.mismatch, with.indels = with.indels)
      rPrimCoord <- vmatchPattern(primerRevRC, sread(sr),
                                  max.mismatch = max.mismatch, with.indels = with.indels)
      # remove reads without match of both primer
      index <- which((lengths(fPrimCoord) > 0) & (lengths(rPrimCoord) > 0))
      fPrimCoord <- fPrimCoord[index]
      rPrimCoord <- rPrimCoord[index]
      sr <- sr[index]
      # remove reads with primer in wrong order
      index <- (unlist(startIndex(rPrimCoord)) - unlist(endIndex(fPrimCoord))) > 0
      fPrimCoord <- unlist(fPrimCoord[index])
      rPrimCoord <- unlist(rPrimCoord[index])
      sr <- sr[index]
      # trim
      start <- end(fPrimCoord) + 1
      end <- start(rPrimCoord) - 1
      sr_trim <- narrow(sr, start = start, end=end)
      srF_prim <- narrow(sr, 
                         start = ifelse(start(fPrimCoord)<1,1,start(fPrimCoord)), 
                         end = end(fPrimCoord))
      srR_prim <- narrow(sr, 
                         start = start(rPrimCoord), 
                         end = ifelse(width(sr)<end(rPrimCoord),width(sr),end(rPrimCoord)))
      if (outputPrimerSequence) {
        writeFastq(srF_prim, file = paste(outputFile, 
                                          "_primerF.fastq.gz", sep = ""), mode = mode, 
                   compress = T)
        writeFastq(srR_prim, file = paste(outputFile, 
                                          "_primerR.fastq.gz", sep = ""), mode = mode, 
                   compress = T)
      }
      writeFastq(sr_trim, file = paste(outputFile, "_F.fastq.gz", 
                                        sep = ""), mode = mode, compress = T)
      mode <- "a"
      filteredReads <- filteredReads + length(sr_trim)
    }
  }
  close(f1)
  suppressWarnings(rm(sr, sr_trim, srF_trim, srR_trim))
  gc()
  gc()
  if (filteredReads == 0) 
    return(data.frame(numReadIn = totalReads, numReadOut = filteredReads, 
                      FileR1 = NA_character_, stringsAsFactors = F))
  else return(data.frame(numReadIn = totalReads, numReadOut = filteredReads, 
                         FileR1 = paste(outputFile, "_F.fastq.gz", sep = ""), 
                         stringsAsFactors = F))
}
