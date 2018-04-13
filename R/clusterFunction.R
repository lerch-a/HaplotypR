
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
    
    require(Rswarm)
    files = c(outfile, structfile, logfile, statsfile, repfile)
    lapply(seq_along(syscall), function(i) {
        if (!all(sapply(files, file.exists))) {
            msg <- paste("Processing file", basename(freqFile[i]), "...", sep=" ")
            progressReport(detail=msg, value=i)
            Rswarm:::.swarmBin(syscall[i])
        } else { 
            msg <- paste("Retrieving existing files for", basename(freqFile[i]), "...", sep=" ")
            progressReport(detail=msg, value=i)        
        }
    })
    return(cbind(RepresentativeFile=repfile, SwarmFile=outfile, StatisticsFile=statsfile, StructureFile=structfile, LogFile=logfile))
}

