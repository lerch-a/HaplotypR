## Table of contents

- [Installation](##Install HaplotypR)
- [Run HaplotypR on R command line](##Run HaplotypR on R command line)
- [Run HaplotypR as Shiny App (currently dysfunctional)](##Run HaplotypR as Shiny App (currently dysfunctional))


## Install HaplotypR
To install HaplotypR start R and first install ShortRead by typing:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ShortRead")
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

detach("package:HaplotypR", unload=TRUE)
devtools::install_gith
```

## Run HaplotypR on R command line
Workflow:
* [Merge sequence read by fixed length](/HaplotypR/documentation/bindReads)
* [Merge by overlapping sequence read](/HaplotypR/documentation/mergeReads)

## Run HaplotypR as Shiny App (currently dysfunctional)
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
install.packages("shiny")
install.packages("shinyFiles")
runShinyApp()
```
