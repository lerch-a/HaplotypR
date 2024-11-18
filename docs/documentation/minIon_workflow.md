# About

HaplotypR is a program for analysis of Amplicon-Seq genotyping experiments. 

The HaplotypR project was developed by Anita Lerch. A paper with more details about the program is available from:

  * Lerch, A. et al. Development Of Amplicon Deep Sequencing Markers And Data Analysis Pipeline For Genotyping Multi-Clonal Malaria Infections. BMC Genomics (2017), 18(1), p.864, http://dx.doi.org/10.1186/s12864-017-4260-y.
  * Lerch, A. et al. Longitudinal tracking and quantification of individual Plasmodium falciparum clones in complex infections. Sci. Rep. 9, 3333 (2019), http://dx.doi.org/10.1038/s41598-019-39656-7.
  * xx
  * xx

# License

HaplotypR is distributed under the GNU General Public License, version 3.

# Installation

To install HaplotypR start R and first install ShortRead and dada2 by typing:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ShortRead","dada2")
```

Then install devtools by typing

```R
install.packages("devtools")
```

and install HaplotypR by typing

```R
library(devtools)
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
dir(file.path("ex3"))
```
The following files should be listed with the last R command: "marker_file.txt", "reads2_F.fastq.gz", "sample_file.txt". 

Run demultiplexing by sample and rename output files
```R
# set input file path
primerFile <- "ex2/marker_file.txt"
sampleFile <- "ex2/sample_file.txt"
reads <- list.files("ex3", pattern="reads", full.names = T)

# create output subdirectory 
outDeplexSample <- file.path(outputDir, "dePlexSample")
dir.create(outDeplexSample)

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

subset: remove markers without sufficient reads 
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
