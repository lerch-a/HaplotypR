# About

HaplotypR is a program for analysis of Amplicon-Seq genotyping experiments. HaplotypR provides a Shiny interface for simpler data analysis.

The HaplotypR project was developed by Anita Lerch. A paper with more details about the program is available from:

Lerch, A. et al. Development Of Amplicon Deep Sequencing Markers And Data Analysis Pipeline For Genotyping Multi-Clonal Malaria Infections. BMC Genomics (2017), 18(1), p.864, http://dx.doi.org/10.1186/s12864-017-4260-y

# License

HaplotypR is distributed under the GNU General Public License, version 3.

# Installation

To install HaplotypR start R and first install ShortRead by typing:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
```

Then install devtools by typing

```R
install.packages("devtools")
install.packages("git2r")
```

and install Rswarm, Rvsearch and HaplotypR by typing

```R
library(devtools)
library(git2r)
path <- file.path(tempfile(pattern="Rswarm-"), "Rswarm")
dir.create(path, recursive=TRUE)
repo <- clone("https://github.com/lerch-a/Rswarm.git", path)
clone("https://github.com/torognes/swarm.git", file.path(path, "src", "swarm"))
install(path)

path <- file.path(tempfile(pattern="Rvsearch-"), "Rvsearch")
dir.create(path, recursive=TRUE)
repo <- clone("https://github.com/lerch-a/Rvsearch.git", path)
clone("https://github.com/torognes/vsearch.git", file.path(path, "src", "vsearch"))
install(path)

devtools::install_github("lerch-a/HaplotypR")
```

# Run HaplotypR as Shiny App

Load HaplotypR package:
```R
library("HaplotypR")
```

Copy Example Files to a working directory 'outputDir':
```R
# Define output directory 
outputDir <- "~/exampleHaplotypR"  
# Create output directoy
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)
# Set working directory to output directory
setwd(outputDir)

# Copy example files to output directory
file.copy(from=system.file(package="HaplotypR", "extdata"), to=outputDir, recursive = T)
# List files example files in output directory
dir(file.path(outputDir, "extdata"))
```
The listed file can be used as example input files in the shiny app. The following files should be listed with the last R command: "barcode_Fwd.fasta", "barcode_Rev.fasta", "markerFile.txt", "readsF.fastq.gz", "readsR.fastq.gz", "sampleFile.txt". 


Run HaplotypR GUI:
```R
runShinyApp()
```

# Run HaplotypR on R command line (without shiny app)

```R
library("HaplotypR")
```

Copy example files to a working directory 'outputDir':
```R
# Define output directory 
outputDir <- "~/exampleHaplotypR"  
# Create output directoy
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)
  
# Set working directory to output directory
setwd(outputDir)

# Copy example files to output directory
file.copy(from=system.file(package="HaplotypR", "extdata"), to=outputDir, recursive = T)

# List files example files in output direcoty
dir(file.path(outputDir, "extdata"))
```
The following files should be listed with the last R command: "barcode_Fwd.fasta", "barcode_Rev.fasta", "markerFile.txt", "readsF.fastq.gz", "readsR.fastq.gz", "sampleFile.txt". 

Run demultiplexing by sample and rename output files
```R
# set input file path
primerFile <- "extdata/markerFile.txt"
sampleFile <- "extdata/sampleFile.txt"
fnBarcodeF <- "extdata/barcode_Fwd.fasta"
fnBarcodeR <- "extdata/barcode_Rev.fasta"
reads <- list.files("extdata", pattern="reads", full.names = T)

# create output subdirectory 
outDeplexSample <- file.path(outputDir, "dePlexSample")
dir.create(outDeplexSample)

# demultiplex by samples
dePlexSample <- demultiplexReads(reads[1], reads[2], fnBarcodeF, fnBarcodeR, outDeplexSample)

# rename output files to sample files
sampleTab <- read.delim(sampleFile, stringsAsFactors=F)
dePlexSample <- renameDemultiplexedFiles(sampleTab, dePlexSample)

# save summary table
write.table(dePlexSample, file.path(outputDir, "demultiplexSampleSummary.txt"), sep="\t", row.names=F)
```

Run demultiplex by marker and truncate primer sequence
```R
# create output subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)
  
# process each marker
markerTab <- read.delim(primerFile, stringsAsFactors=F)
dePlexMarker <- lapply(1:dim(markerTab)[1], function(j){
  mID <- as.character(markerTab[j, "MarkerID"])
  adapterF <- as.character(markerTab[j, "Forward"])
  adapterR <- as.character(markerTab[j, "Reverse"])
  
  # process each file demultiplexed by sample
  res <- lapply(seq_along(dePlexSample$FileR1), function(i){
    outputFile <- file.path(outDeplexMarker, sub("R1\\.fastq.gz", mID, basename(as.character(dePlexSample$FileR1)[i])))
    
    # demultiplex by marker and truncate primer sequence
    removePrimer(as.character(dePlexSample$FileR1)[i], as.character(dePlexSample$FileR2)[i], outputFile, 
                 adapterF, adapterR, max.mismatch=2, with.indels=F)
  })
  cbind(BarcodePair=as.character(dePlexSample$BarcodePair), MarkerID=mID, do.call(rbind, res))
})
dePlexMarker <- do.call(rbind, dePlexMarker)
dePlexMarker <- merge.data.frame(dePlexSample[,c("SampleID", "SampleName", "BarcodePair")], dePlexMarker, by="BarcodePair")

# save summary table
write.table(dePlexMarker, file.path(out, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
```

Trim and fuse paired reads
```R
# create output subdirectory 
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

# Trim options
numNtF <- 190
numNtR <- 120
concatSuffix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

# Adjust reference to trim options and save as fasta file
refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], concatSuffix, ".fasta", sep="")))
})

# Fuse paired read
procReads <- bindAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles, 
                         read1Length=numNtF, read2Length=numNtR)
procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", concatSuffix)), sep="\t", row.names=F)
```

Calculate mismatch rate and call SNPs
```R
# Options
minMMrate <- 0.5
minOccGen <- 2

# process each marker
snpLst <- lapply(markerTab$MarkerID, function(marker){
  # Calculate mismatch rate
  seqErrLst <- calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]), 
                                            refSeq[marker], 
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                            minCoverage=100L)
  names(seqErrLst) <- procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, concatSuffix)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  snps <- cbind(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N")
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                  minMMrate*100, minOccGen, marker, concatSuffix)), 
              row.names=F, col.names=T, sep="\t", quote=F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                             minMMrate*100, minOccGen, marker, concatSuffix)), 
      width=1500 , height=600)
  matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
          ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
  abline(v=snps[,"Pos"], lty=2, col="grey")
  abline(h=minMMrate, lty=1, col="red")
  dev.off()
  
  return(snps)
})
names(snpLst) <- markerTab$MarkerID
```


Call Haplotypes (Work in progress.)
```R
```
