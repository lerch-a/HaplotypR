
clusterReads <- function(freqFile, outputDir, prefixOutputFiles, method="Swarm", progressReport=message){
  outfile <- file.path(outputDir, paste(prefixOutputFiles, ".swarms", sep=""))
  structfile <- file.path(outputDir, paste(prefixOutputFiles, ".struct", sep=""))
  logfile <- file.path(outputDir, paste(prefixOutputFiles, ".log.txt", sep=""))
  statsfile <- file.path(outputDir, paste(prefixOutputFiles, ".stats.txt", sep=""))
  repfile <- file.path(outputDir, paste(prefixOutputFiles, ".representatives.fasta", sep=""))
  syscall <- paste("-f -b 3 -w", repfile, "-l", logfile, "-s", statsfile, "-i", structfile, "-o", outfile, freqFile, sep=" ")
  
  # check and set progress report function
  if(!is.function(progressReport))
    progressReport <- message
  
  # require(Rswarm)
  lapply(seq_along(syscall), function(i){
    msg <- paste("Processing file", basename(freqFile[i]), "...", sep=" ")
    progressReport(detail=msg, value=i)
    .swarmBin(syscall[i])
  })
  return(cbind(RepresentativeFile=repfile, SwarmFile=outfile, StatisticsFile=statsfile, StructureFile=structfile, LogFile=logfile))
}


extractSwarmClusterReads <- function(swarmFile, swarmRepresentativeFile, fastafile, outputDir){
  require(Biostrings)
  require(ShortRead)
  prefix <- sub(".swarms","", basename(swarmFile))
  outDir <- file.path(outputDir, prefix)
  # i <- 1
  lapply(seq_along(prefix), function(i){
    dir.create(outDir[i])
    rep <- readFasta(repfile[i])
    hapInput <- readFasta(fastafile[i])
    swarm <- readLines(swarmFile[i])
    # j <- 1
    mode <- "w"
    lapply(seq_along(swarm), function(j){
      hap <- strsplit(swarm[[j]], " ")[[1]]
      idx <- id(hapInput) %in% hap
      #as.character(id(rep[j]))
      repFreq <- as.integer(strsplit(as.character(id(rep[j])), "_")[[1]][2])
      if(repFreq>1){ # check frequencie
        #aln <- pairwiseAlignment(hapInput[idx], sread(rep[j])) # patternQuality, subjectQuality
        #nmismatch(aln)
        #mismatchSummary(aln)
        #nedit(aln)
        #coverage(aln)
        #consensusMatrix(aln)
        #freq <- as.integer(do.call(rbind, strsplit(as.character(id(hapInput[idx])), "_"))[,2])
        #consensusString(aln, ambiguityMap="N", threshold=0.51) == sread(rep[j])
        writeFasta(hapInput[idx], file.path(outDir[i], paste(as.character(id(rep[j])), "fasta", sep=".")))
      }else{
        writeFasta(hapInput[idx], file.path(outDir[i], "singelton.fasta"), mode=mode)
        mode <- "a"
      }
    })
  })
}
