# Workflow: Merge by overlapping sequence read

## Load R libraries
```R
library("HaplotypR")
library("ShortRead")
```

## Load example data
Copy example files to a working directory 'outputDir':
```R
# Define output directory 
outputDir <- "exampleHaplotypR"  
# Create output directoy
if(!dir.exists(outputDir))
  dir.create(outputDir, recursive=T)

# Copy example files to output directory
file.copy(from=system.file(package="HaplotypR", "extdata/ex2"), to=".", recursive = T)

# List files example files in output direcoty
dir(file.path(outputDir, "ex2"))
```
The following files should be listed with the last R command: "barcode_F.fasta", "barcode_R.fasta", "marker_file.txt", "reads2_F.fastq.gz", "reads2_R.fastq.gz", "sample_file.txt". 

## Demultiplex sequence reads by sample
Run demultiplexing by sample and rename output files
```R
# set input file path
primerFile <- "ex2/marker_file.txt"
sampleFile <- "ex2/sample_file.txt"
fnBarcodeF <- "ex2/barcode_F.fasta"
fnBarcodeR <- "ex2/barcode_R.fasta"
reads <- list.files("ex2", pattern="reads", full.names = T)

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

## Demulitplex by marker
Run demultiplex by marker and truncate primer sequence
```R
# create output subdirectory 
outDeplexMarker <- file.path(outputDir, "dePlexMarker")
dir.create(outDeplexMarker)
  
# process each marker
markerTab <- read.delim(primerFile, stringsAsFactors=F)
dePlexMarker <- demultiplexByMarker(dePlexSample, markerTab, outDeplexMarker)

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
```

## Merge by overlapping paired sequence read
Fuse paired reads. Method work only for overlapping sequence read pair by merging the overlap of the forward and reverse read (using vsearch wrapper).
```R
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

postfix <- "_merge"
refSeq <- DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

procReads <- mergeAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles)
procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F, quote=F)
```

## Call SNPs
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
  write.table(seqErr, file.path(outputDir, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
  snpRef <- unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  snps <- data.frame(Chr=marker, Pos=potSNP, Ref=snpRef, Alt="N", stringsAsFactors=F)
  write.table(snps, file=file.path(outputDir, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                  minMMrate*100, minOccGen, marker, postfix)), 
              row.names=F, col.names=T, sep="\t", quote=F)
  
  # Plot mismatch rate and SNP calls
  png(file.path(outputDir, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png", 
                             minMMrate*100, minOccGen, marker, postfix)), 
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


## Call Haplotypes
```R
# call haplotype options
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 2
minCovSample <- 25

# call final haplotypes
finalTab <- createFinalHaplotypTable(
  outputDir = outputDir, sampleTable = procReads, markerTable = markerTab, referenceSeq = refSeq,
  snpList = snpLst, postfix = postfix, minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample)
```

## Get haplotyp list and sequence
Todo
