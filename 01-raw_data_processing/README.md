# Raw data processing

## 1. Install R packages (v3.6.1)

```R
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("dada2")
# BiocManager::install("ShortRead")
# BiocManager::install("Biostrings")
# BiocManager::install("phyloseq")
# install.packages("tidyverse")
# install.packages("stringr")
# install.packages("data.table")
# install.packages("broom")
# install.packages("qualpalr")
# install.packages("viridis")
# install.packages("seqinr")
# install.packages("ape")
# install.packages("phytools")
```
## 2. Load required libraries

```R
library(dada2)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(phyloseq)
library(ape)
library(phytools)
```

