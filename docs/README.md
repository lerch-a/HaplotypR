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

```bash
Usage: ~/tools/HaplotypR/R/run_haplotypr.R [options]

Options:
    -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
        Directory to save the output.

    -p AMPLICONS_FILE, --amplicons_file=AMPLICONS_FILE
        File with fwd/rev primers, reference seqs, fwd/rev read lengths and max indel threshold listed by amplicon.

    -s SAMPLES_DIR, --samples_dir=SAMPLES_DIR
        Directory with demultiplexed sample files.

    -t, --trim_reads
        If passed will look for fwd/rev read lengths in amplicons_file and trim reads to those lengths.

    --min_mismatch=MIN_MISMATCH
        Minimum rate of mismatch between haplotype and reference sequence.

    --min_genotype_occurrence=MIN_GENOTYPE_OCCURRENCE
        Minimum # of samples for a valid genotype to be called in.

    --detection_limit=DETECTION_LIMIT
        Minimum frequency for detecting a haplotype.

    --min_haplotype_coverage=MIN_HAPLOTYPE_COVERAGE
        Minimum coverage for haplotype to be recognized as valid.

    --min_haplotype_occurrence=MIN_HAPLOTYPE_OCCURRENCE
        Minimum # of samples for a valid haplotype to be called in.

    --min_sample_coverage=MIN_SAMPLE_COVERAGE
        Minimum coverage for a sample to be recognized as valid.

    -v, --verbose
        Run verbosely.

    -h, --help
        Show this help message and exit

```