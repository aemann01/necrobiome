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

## 3. File path setup

```R
rawpath <- "/Volumes/histolytica/necrobiome/raw"
wdpath <- "~/github/necrobiome/01-raw_data_processing/"
fnFs <- sort(list.files(rawpath, pattern="_R1_001.fastq.gz", full.names=T))
fnRs <- sort(list.files(rawpath, pattern="_R2_001.fastq.gz", full.names=T))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
message("sample names:")
sample.names
```

```text
sample names:
    'Blank''BlankE''Negctrl''NegCtrl''S01A''S02E''S04A''S05E''S06A''S07E''S08A''S09E''S10A''S11E''S12A''S13E''S14A''S15E''S16A''S17E''S18A''S19E''S20A''S21E''S24A''S25E''S26A''S27E''S30A''S31E''S32A''S33E''S34A''S35E''S36A''S37E''S48A''S49E''W03A''W04E''W06E''W07A''W09A''W10E''W11A''W12E''W13A''W14E''W15A''W16E''W17A''W18E''W19A''W20E''W23A''W24E''W25A''W26E''W27A''W28E''W29A''W30E''W31A''W32E'
```

## 4. Plot quality scores

```R
system("mkdir img")
pdf(paste(wdpath, "img/", "forward_quality_plot.pdf", sep=""))
fq <- plotQualityProfile(fnFs[5:15])
print(fq)
dev.off()
pdf(paste(wdpath, "img/", "reverse_quality_plot.pdf", sep=""))
rq <- plotQualityProfile(fnRs[5:15])
print(rq)
dev.off()
print(fq)
print(rq)
```

![forward quality plot](01-raw_data_processing/img/forward_quality_plot.png)
![reverse quality plot](img/reverse_quality_plot.png)

## 5. Preliminary filter (removes sequences with N's)

```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
```
