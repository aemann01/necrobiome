# Raw data processing

## Setup

Download raw data from ////

```bash
mkdir raw
cd raw
wget ////
cd ..
```

### 1. Install R packages (v3.6.1)

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
### 2. Load required libraries

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

### 3. File path setup

```R
rawpath <- "raw"
wdpath <- "~/necrobiome/01-raw_data_processing/" # if you didn't save to home folder, change to correct path
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

### 4. Plot quality scores

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

![forward quality plot](https://github.com/aemann01/necrobiome/blob/master/01-raw_data_processing/imgs/forward_quality_plot.png)
![reverse quality plot](https://github.com/aemann01/necrobiome/blob/master/01-raw_data_processing/imgs/reverse_quality_plot.png)

### 5. Preliminary filter (removes sequences with N's)

```R
fnFs.filtN <- file.path(rawpath, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(rawpath, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress = TRUE)
```

### 6. Primer removal

```R
cutadapt <- as.character(system("which cutadapt", intern=T))
system("cutadapt --version")
path.cut <- file.path(rawpath, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc("GTGYCAGCMGCCGCGGTAA")
REV.RC <- dada2:::rc("GGACTACNVGGGTWTCTAAT")
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", "GTGYCAGCMGCCGCGGTAA", "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", "GGACTACNVGGGTWTCTAAT", "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,"-o", fnFs.cut[i], "-p", fnRs.cut[i], fnFs.filtN[i], fnRs.filtN[i]))
}
cutFs <- sort(list.files(path.cut, pattern = "R1", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "R2", full.names = TRUE))
```

### 7. Filter and trim reads

```R
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, trimRight=50, minLen = c(150,120), maxN=c(0,0), maxEE=c(2,2), truncQ=c(2,2), rm.phix=TRUE, matchIDs=TRUE, compress=TRUE, multithread=TRUE)
retained <- as.data.frame(out)
retained$percentage_retained <- retained$reads.out/retained$reads.in*100
retained
```

```text
                                 reads.in reads.out percentage_retained
Blank_S35_L001_R1_001.fastq.gz     481523    436419            90.63305
BlankE_S94_L001_R1_001.fastq.gz       958       487            50.83507
Negctrl_S34_L001_R1_001.fastq.gz   293284    262418            89.47573
NegCtrl_S95_L001_R1_001.fastq.gz     4056      3078            75.88757
S01A_S68_L001_R1_001.fastq.gz      282815    257473            91.03937
S02E_S69_L001_R1_001.fastq.gz      162192    148004            91.25234
S04A_S1_L001_R1_001.fastq.gz       159433    141282            88.61528
S05E_S2_L001_R1_001.fastq.gz       268711    240115            89.35808
S06A_S3_L001_R1_001.fastq.gz       294200    262951            89.37831
S07E_S70_L001_R1_001.fastq.gz      184791    161918            87.62223
S08A_S71_L001_R1_001.fastq.gz      204775    180186            87.99219
S09E_S4_L001_R1_001.fastq.gz       267348    238461            89.19498
S10A_S5_L001_R1_001.fastq.gz       247957    232296            93.68399
S11E_S6_L001_R1_001.fastq.gz       325705    300339            92.21197
S12A_S72_L001_R1_001.fastq.gz      273152    251139            91.94112
S13E_S73_L001_R1_001.fastq.gz      143669    131726            91.68714
S14A_S7_L001_R1_001.fastq.gz       323750    292292            90.28324
S15E_S8_L001_R1_001.fastq.gz       323700    291564            90.07229
S16A_S74_L001_R1_001.fastq.gz      189223    169938            89.80832
S17E_S53_L001_R1_001.fastq.gz      436983    413395            94.60208
S18A_S9_L001_R1_001.fastq.gz       330394    296941            89.87482
S19E_S10_L001_R1_001.fastq.gz      319303    289852            90.77647
S20A_S11_L001_R1_001.fastq.gz      255848    234347            91.59618
S21E_S12_L001_R1_001.fastq.gz      312902    287657            91.93198
S24A_S75_L001_R1_001.fastq.gz      128753    118165            91.77650
S25E_S76_L001_R1_001.fastq.gz      185999    171827            92.38060
S26A_S13_L001_R1_001.fastq.gz      245292    226229            92.22845
S27E_S14_L001_R1_001.fastq.gz      340998    303188            88.91196
S30A_S15_L001_R1_001.fastq.gz      291497    257305            88.27021
S31E_S16_L001_R1_001.fastq.gz      222782    183182            82.22478
S32A_S77_L001_R1_001.fastq.gz      165190    148785            90.06901
S33E_S78_L001_R1_001.fastq.gz      125758    115564            91.89396
S34A_S17_L001_R1_001.fastq.gz      485443    413588            85.19806
S35E_S18_L001_R1_001.fastq.gz      392980    335547            85.38526
S36A_S19_L001_R1_001.fastq.gz      393499    342614            87.06858
S37E_S20_L001_R1_001.fastq.gz      372373    323852            86.96979
S48A_S79_L001_R1_001.fastq.gz       62542     53600            85.70241
S49E_S80_L001_R1_001.fastq.gz      235341    201294            85.53291
W03A_S81_L001_R1_001.fastq.gz      219571    199004            90.63310
W04E_S82_L001_R1_001.fastq.gz      182648    168858            92.44996
W06E_S83_L001_R1_001.fastq.gz      261847    239522            91.47403
W07A_S84_L001_R1_001.fastq.gz      111638    102056            91.41690
W09A_S21_L001_R1_001.fastq.gz      321086    297306            92.59388
W10E_S22_L001_R1_001.fastq.gz      472254    431200            91.30680
W11A_S85_L001_R1_001.fastq.gz      158644    141245            89.03268
W12E_S86_L001_R1_001.fastq.gz      124480    113300            91.01864
W13A_S87_L001_R1_001.fastq.gz      149514    134372            89.87252
W14E_S88_L001_R1_001.fastq.gz      125636    113540            90.37219
W15A_S23_L001_R1_001.fastq.gz      417570    376836            90.24499
W16E_S24_L001_R1_001.fastq.gz      277670    249288            89.77851
W17A_S25_L001_R1_001.fastq.gz      299571    263847            88.07495
W18E_S26_L001_R1_001.fastq.gz        4529      1787            39.45683
W19A_S27_L001_R1_001.fastq.gz        1808       976            53.98230
W20E_S28_L001_R1_001.fastq.gz      318742    272394            85.45909
W23A_S89_L001_R1_001.fastq.gz      227815    205604            90.25042
W24E_S90_L001_R1_001.fastq.gz      130542    118786            90.99447
W25A_S91_L001_R1_001.fastq.gz      184623    162714            88.13311
W26E_S92_L001_R1_001.fastq.gz      175793    155599            88.51263
W27A_S93_L001_R1_001.fastq.gz      181337    154429            85.16133
W28E_S29_L001_R1_001.fastq.gz      296214    233156            78.71201
W29A_S30_L001_R1_001.fastq.gz      336092    283659            84.39921
W30E_S31_L001_R1_001.fastq.gz      385568    292358            75.82528
W31A_S32_L001_R1_001.fastq.gz      344667    232620            67.49123
W32E_S33_L001_R1_001.fastq.gz      459447    412766            89.83974
```

### 8. Learn and plot error rates

```R
errF <- learnErrors(filtFs, multithread=T, random=T)
errR <- learnErrors(filtRs, multithread=T, random=T)
png(paste(wdpath, "imgs/", "error_plot.png", sep=""))
ep <- plotErrors(errF, nominalQ=TRUE) 
plot(ep)
dev.off()
print(ep)
```

![reverse quality plot](https://github.com/aemann01/necrobiome/blob/master/01-raw_data_processing/imgs/error_plot.png)

### 9. Dereplication

```R
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
```

### 10. Sample inference

```R
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
```

### 11. Filter out samples with fewer than 1000 reads

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN))
samples_to_keep <- track[,4] > 1000
samples_to_remove <- names(samples_to_keep)[which(samples_to_keep == FALSE)]
paste("number of samples removed: ", samples_to_remove)
```

### 12. Merge paired end reads

```R
mergers <- mergePairs(dadaFs[samples_to_keep], derepFs[samples_to_keep], dadaRs[samples_to_keep], derepRs[samples_to_keep], verbose=T, )
```

### 13. Construct sequence table

```R
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

```text
```

### 14. Sequence length distribution plot

```R
length.histogram <- as.data.frame(table(nchar(getSequences(seqtab))))
png(paste(wdpath, "img/", "length_hist.png", sep=""))
plot(x=length.histogram[,1], y=length.histogram[,2])
dev.off()
```

![reverse quality plot](https://github.com/aemann01/necrobiome/blob/master/01-raw_data_processing/imgs/length_hist.png)

### 15. Remove chimeras

```R
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=T, verbose=T)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
```

### 16. Processing summary

```R
getN <- function(x) sum(getUniques(x))
track <- cbind(out[samples_to_keep,], sapply(dadaFs[samples_to_keep], getN), sapply(dadaRs[samples_to_keep], getN), sapply(mergers, getN), rowSums(seqtab.nochim))
track <- cbind(track, 100-track[,6]/track[,5]*100, 100-track[,7]/track[,6]*100)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nochimeras", "percent_singletons", "percent_chimeras")
rownames(track) <- sample.names[samples_to_keep]
```

### 17. Save output

```R
write.table(data.frame("row_names"=rownames(track),track),"read_retention.16s.txt", row.names=FALSE, quote=F, sep="\t")
write.table(data.frame("row_names"=rownames(seqtab.nochim),seqtab.nochim),"sequence_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")
uniquesToFasta(seqtab.nochim, "rep_set.fa")
system("awk '/^>/{print \">ASV\" ++i; next}{print}' < rep_set.fa > rep_set_fix.fa")
system("mv rep_set_fix.fa rep_set.fa")
```

### 18. Clean up ASV names

```R
my_otu_table <- t(as.data.frame(seqtab.nochim)) 
ASV.seq <- as.character(unclass(row.names(my_otu_table))) 
ASV.num <- paste0("ASV", seq(ASV.seq), sep='') 
colnames(seqtab.nochim) <- ASV.num 
```

### 19. Assign taxonomy using VSEARCH and QIIME2

```R
# get executable paths
qiimepath <- system("conda env list | grep 'qiime' | awk '{print $2}'", intern=T)
qiimepathFix <- paste(qiimepath, "/bin/qiime", sep="")
qiimepathFix
# reference files
BACREF <- "ezbiocloud_qiime_full.fasta"
BACTAX <- "ezbiocloud_id_taxonomy.txt"
# format and run
system(paste(qiimepathFix, "tools import --input-path rep_set.fa --output-path rep_set.qza --type 'FeatureData[Sequence]'")
system("rm -r assigntax")
system(sprintf(paste(qiimepathFix, "feature-classifier classify-consensus-vsearch --i-query rep_set.qza --i-reference-reads %s --i-reference-taxonomy %s --output-dir assigntax"), BACREF, BACTAX)
system("unzip assigntax/classification.qza -d assigntax/")
#get file path for taxonomy file
tempfile <- subset(dir(path="assigntax"), !grepl("classification.qza", dir(path="assigntax/")))
newpath <- paste("assigntax/", tempfile, "/data/taxonomy.tsv", sep="")
#make 7 level taxonomy file -- test this out on tortoise data
system(sprintf("awk '{print $2}' %s | sed '1d' > taxonomy_strings.txt", newpath))
system(sprintf("awk '{print $1}' %s | sed '1d' > asv_ids.txt", newpath))
system("python fix_taxonomy_L7.py taxonomy_strings.txt > fix_string.txt")
system("paste asv_ids.txt fix_string.txt > taxonomy_L7.txt")
system("rm taxonomy_strings.txt fix_string.txt asv_ids.txt")
system("head taxonomy_L7.txt")
```

```text
```

### 20. Combine sequence and taxonomy tables

```R
#taxa will be the rows, columns will be samples, followed by each rank of taxonomy assignment, from rank1 (domain-level) to rank7/8 (species-level), followed by accession (if applicable)
#first check if the row names of the taxonomy table match the column headers of the sequence table
taxa <- read.table(newpath, header=T, sep="\t", row.names=1)
length(which(row.names(taxa) %in% colnames(seqtab.nosingletons.nochim)))
dim(taxa)
dim(seqtab.nosingletons.nochim)
#the number of taxa from the last three commands should match
#now ensure that the taxa in the tables are in the same order #this should be true if you haven't reordered one or the other of these matrices inadvertently
order.col <- row.names(taxa)
seqtab.nosingletons.nochim <- seqtab.nosingletons.nochim[,order.col]
row.names(taxa) == colnames(seqtab.nosingletons.nochim) #IMPORTANT: only proceed if this evaluation is true for every element. if it isn't you'll need to re-order your data. I'd suggest sorting both matrices by their rows after transposing the sequence table.
#as long as the ordering of taxa is the same, you can combine like this (note you need to transpose the sequence table so that the taxa are in the rows)
sequence_taxonomy_table <- cbind(t(seqtab.nosingletons.nochim), taxa)
#now write to file
write.table(data.frame("row_names"=rownames(sequence_taxonomy_table),sequence_taxonomy_table),"sequence_taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")
```

### 21. Filter out unwanted taxonomic groups

```R
system("grep -v -E 'Unassigned|Bacteria;unknown' sequence_taxonomy_table.16s.merged.txt | awk '{print $1}' | grep 'A' > wanted.ids")
wanted <- read.table("wanted.ids", header=F)
seqtab.filtered <- seqtab.nosingletons.nochim[, which(colnames(seqtab.nosingletons.nochim) %in% wanted$V1)]
write.table(data.frame("row_names"=rownames(seqtab.filtered),seqtab.filtered),"sequence_table.16s.filtered.txt", row.names=FALSE, quote=F, sep="\t")
```

### 22. Generate representative sequence tree

```R
system("seqtk subseq rep_set.fa wanted.ids > rep_set.filt.fa")
system("mafft --auto rep_set.filt.fa > rep_set.filt.align.fa")
system("fasttree -nt rep_set.filt.align.fa > rep_set.filt.tre")
#midpoint root tree
tre <- read.tree("rep_set.filt.tre")
root <- midpoint.root(tre)
write.tree(root, "rep_set.root.tre")
```
