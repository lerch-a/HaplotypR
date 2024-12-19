# About

HaplotypR is a program for analysis of Amplicon-Seq genotyping experiments. 

The HaplotypR project was developed by Anita Lerch. A paper with more details about the program is available from:

  * Lerch, A. et al. Development Of Amplicon Deep Sequencing Markers And Data Analysis Pipeline For Genotyping Multi-Clonal Malaria Infections. BMC Genomics (2017), 18(1), p.864, http://dx.doi.org/10.1186/s12864-017-4260-y.
  * Lerch, A. et al. Longitudinal tracking and quantification of individual Plasmodium falciparum clones in complex infections. Sci. Rep. 9, 3333 (2019), http://dx.doi.org/10.1038/s41598-019-39656-7.
  * Holzschuh A, et al. (2024) Using a mobile nanopore sequencing lab for end-to-end genomic surveillance of Plasmodium falciparum: A feasibility study. PLOS Glob Public Health 4(2): e0002743, https://doi.org/10.1371/journal.pgph.0002743.
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

# Copy example files to working directory
file.copy(from=system.file(package="HaplotypR", "extdata/ex3"), to=".", recursive = T)

# List files example files in output direcoty
dir(file.path("ex3"))
```
The following files should be listed with the last R command: "marker_file_ex3.txt", "sample_file_ex3.txt" and others. 

Run demultiplexing by sample and rename output files
```R
# set input file path
primerFile <- "ex3/marker_file_ex3.txt"
sampleFile <- "ex3/sample_file_ex3.txt"
readsDir <- "ex3/read_dir_ex3"

# create output subdirectory 
outDeplexSample <- file.path(outputDir, "dePlexSample")
dir.create(outDeplexSample, recursive=T)

# rename sample files and merge to a single file per sample if needed
sampleTab <- read.delim(sampleFile, stringsAsFactors=F)
dePlexSample <- mergeMinIONfiles(inDir=readsDir, outDir=outDeplexSample, sampleTab=sampleTab)

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
# shorten primer sequence to same length for demultiplexing
markerTab$Forward <- substr(markerTab$Forward, start=nchar(markerTab$Forward)-20, stop=nchar(markerTab$Forward))
markerTab$Reverse <- substr(markerTab$Reverse, start=nchar(markerTab$Reverse)-20, stop=nchar(markerTab$Reverse))
dePlexMarker <- demultiplexByMarkerMinION(dePlexSample, markerTab, outDeplexMarker, max.mismatch=2)

# save summary table
write.table(dePlexMarker, file.path(outputDir, "demultiplexMarkerSummary.txt"), sep="\t", row.names=F)
```

Call Haplotypes
```R
# call haplotype options
minCov <- 3
detectionLimit <- 1/100
minOccHap <- 1
minCovSample <- 100

# call final haplotypes
finalTab <- createFinalHaplotypTableDADA2(
  outputDir = outputDir, sampleTable = dePlexMarker, markerTable = markerTab,
  referenceSequence=NULL, filterIndel=T,
  minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample,
  multithread=FALSE, pool="pseudo", OMEGA_A=1e-120)
  
write.csv(finalTab, file=file.path(outputDir, "finalHaplotypList_vMinION.csv"), row.names=F)
```

Calculate mismatch rate and call SNPs
```R
dePlexMarker <- dePlexMarker[dePlexMarker$]

refSeq <- DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) <- markerTab$MarkerID
snpLst <- createSNPsList(outputDir, sampleTable=dePlexMarker, markerTable=markerTab, refSeq=refSeq, postfix=postfix)
```


```R
# Assign Strain
strainLst <- unlist(lapply(markerTab$MarkerID, function(marker){
  haplotypes <- readDNAStringSet(file.path(outputDir, sprintf("HaplotypeList_%s.fasta", marker)))
  strain <- assignStrain(haplotypes, refFn=sprintf("ex3/%s.fasta", marker))
  return(strain)
}))
finalTab$Strain <- strainLst[finalTab$Haplotype]
```


```R
# Assign SNP
mutationTab <- read.delim("ex3/mutationList.txt")
mutationLst <- unlist(lapply(markerTab$MarkerID, function(marker){
  haplotypes <- readDNAStringSet(file.path(outputDir, sprintf("HaplotypeList_%s.fasta", marker)))
  mutTab <- mutationTab[mutationTab$MarkerID==marker,]
  mut <- assignSNP(haplotypes, snpTab=mutTab)
  return(mut)
}))
finalTab$Mutations <- mutationLst[finalTab$Haplotype]
```