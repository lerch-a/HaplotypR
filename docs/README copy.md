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
dePlexMarker <- merge.data.frame(dePlexFiles[,c("SampleID", "BarcodePair")], dePlexMarker, by="BarcodePair")

# save summary table
write.table(dePlexMarker, file.path(out, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
```

Call SNPs 
```R
```

Call SNPs 
```R
# Options
minMMrate <- 0.5  # minimum MismatchRate
minOccGen <- 2    # mininimum Replicate

# Load Sample
sampleTable <- read.delim(file.path(outputDir, typerPrj$pr ojects[prj, "filename"]), stringsAsFactors=T)
sampleTable <- sampleTable[sampleTable$MarkerID == marker,]
sampleTable <- sampleTable[!(is.na(sampleTable$ReadFile) | is.na(sampleTable$SampleID)),]

# Load Reference
fn <- file.path(outputDir, paste(marker, prj, ".fasta", sep=""))
refSeq <- readDNAStringSet(fn)


# Calculate mismatch rate
seqErrLst <- calculateMismatchFrequencies(as.character(sampleTable$ReadFile), refSeq, minCoverage=100L)
names(seqErrLst) <- sampleTable$SampleID
seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
  l[,"MisMatch"]/l[,"Coverage"]
}))
write.table(seqErr, file.path(out, sprintf("mismatchRate_rate%.0f_occ%i_%s%s.txt", 
                                           minMMrate*100, minOccGen, marker, prj)), sep="\t", row.names=F)

potSNP <- callGenotype(seqErr, minMismatchRate=minMMrate, minReplicate=minOccGen)
snpRef <- unlist(lapply(potSNP, function(snp){
  as.character(subseq(refSeq, start=snp, width=1))
}))
snpLst <- cbind(Chr=names(refSeq), Pos=potSNP, Ref=snpRef, Alt="N")
write.table(snpLst, file=file.path(out, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt", 
                                                minMMrate*100, minOccGen, marker, prj)), 
            row.names=F, col.names=T, sep="\t", quote=F)

matplot(seqErr, type="p", pch=16, cex=0.4, col="#00000088", ylim=c(0, 1),
        ylab="Mismatch Rate", xlab="Base Position", main=marker, cex.axis=2, cex.lab=2)
abline(v=snpLst[,"Pos"], lty=2, col="grey")
abline(h=minMMrate, lty=1, col="red")
```

Call Haplotypes
```R
```
