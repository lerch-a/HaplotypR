# About

#### Still under development

HaplotypR is a program for analysis of Amplicon-Seq genotyping experiments. HaplotypR provides a Shiny interface for simpler data analysis.

The HaplotypR project was developed by Anita Lerch. A paper with more details about the program is available from:

Lerch A, et al. Development Of Amplicon Deep Sequencing Markers And Data Analysis Pipeline For Genotyping Multi-Clonal Malaria Infections. bioRxiv (2017), http://dx.doi.org/10.1101/121426. 


# Installation

To install HaplotypR start R and first install ShortRead by typing:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("ShortRead")
```

Then install devtools by typing

```R
install.packages("devtools")
```

and install HaplotypR by typing

```R
devtools::install_github("lerch-a/HaplotypR")
```

Next load HaplotypR with

```R
library("HaplotypR")
```


# License

HaplotypR is distributed under the GNU General Public License, version 3.