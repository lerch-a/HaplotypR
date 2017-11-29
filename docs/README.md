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
# List files example files in output direcoty
dir(file.path(outputDir, "extdata"))
```

Run HaplotypR GUI:
```R
runShinyApp()
```

# Run HaplotypR in R command line (without shiny app)

Work in progress.
