# About

HaplotypR is a program for analysis of Amplicon-Seq genotyping experiments. HaplotypR provides a Shiny interface for simpler data analysis.

The HaplotypR project was developed by Anita Lerch. A paper with more details about the program is available from:

  * Lerch, A. et al. Development Of Amplicon Deep Sequencing Markers And Data Analysis Pipeline For Genotyping Multi-Clonal Malaria Infections. BMC Genomics (2017), 18(1), p.864, http://dx.doi.org/10.1186/s12864-017-4260-y.
  * Lerch, A. et al. Longitudinal tracking and quantification of individual Plasmodium falciparum clones in complex infections. Sci. Rep. 9, 3333 (2019), http://dx.doi.org/10.1038/s41598-019-39656-7.

# License

HaplotypR is distributed under the GNU General Public License, version 3.

# Installation

Note:
* For Windows computer: install the precomplied NGmergeR package from https://github.com/lerch-a/NGmergeR/releases) and  install.packages("PathToPackage/NGmergeR_0.99.1-win_R-v3.6.zip", repos = NULL).

To install HaplotypR start R and first install ShortRead and dada2 by typing:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ShortRead","dada2")
```

Then install devtools by typing

```R
install.packages("devtools")
install.packages("git2r")
```

and install NGmergeR and HaplotypR by typing

```R
library(devtools)
library(git2r)

path <- file.path(tempfile(pattern="NGmergeR-"), "NGmergeR")
dir.create(path, recursive=TRUE)
repo <- clone("https://github.com/lerch-a/NGmergeR.git", path)
clone("https://github.com/jsh58/NGmerge", file.path(path, "src", "NGmerge"))
install(path)

detach("package:HaplotypR", unload=TRUE)
devtools::install_github("lerch-a/HaplotypR")
```

# Run HaplotypR on R command line

```R
library("HaplotypR")
library("ShortRead")
```

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
dir(file.path("ex2"))
```
The following files should be listed with the last R command: "barcode_F.fasta", "barcode_R.fasta", "marker_file.txt", "reads2_F.fastq.gz", "reads2_R.fastq.gz", "sample_file.txt". 

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

Fuse paired reads. Two methods are provided for fusing paired reads. 
First method works for non overlapping sequence read pairs. Trim to a fixed length (removes low quality bases) and then concatenate forward and reverse read.
```R
# create output subdirectory 
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

# Trim options
numNtF <- 190
numNtR <- 120
postfix <- sprintf("_bind%.0f_%.0f", numNtF, numNtR)

# Adjust reference to trim options and save as fasta file
refSeq <- as.character(markerTab$ReferenceSequence)
refSeq <- DNAStringSet(paste(substr(refSeq, 1,numNtF), substr(refSeq, nchar(refSeq)+1-numNtR, nchar(refSeq)), sep=""))
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

# Fuse paired read
procReads <- bindAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles, 
                         read1Length=numNtF, read2Length=numNtR)
procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F)
```

Second method work only for overlapping sequence read pair by merging the overlap of the forward and reverse read (using vsearch wrapper).
```R
# create output subdirectory 
outProcFiles <- file.path(outputDir, "processedReads")
dir.create(outProcFiles)

postfix <- "_merge"
refSeq <- DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) <- markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(outputDir, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

procReadsMerge <- mergeAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), outProcFiles, method="NGmerge")
procReads <- cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReadsMerge)
write.table(procReads, file.path(outputDir, sprintf("processedReadSummary%s.txt", postfix)), sep="\t", row.names=F, quote=F)

# subset: remove markers without reads 
procReads <- procReads[procReads$numRead>10,]

```

Call Haplotypes
```R
# call haplotype options
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 2
minCovSample <- 25

# remove samples without reads
procReads <- procReads[procReads$numRead>0,]

# call final haplotypes
finalTab <- createFinalHaplotypTableDADA2(
  outputDir = outputDir, sampleTable = procReads, markerTable = markerTab,
  referenceSequence=refSeq, filterIndel=T,
  minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample,
  multithread=FALSE, pool="pseudo", OMEGA_A=1e-120)
```

