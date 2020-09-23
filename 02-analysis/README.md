# Visualization (R v3.6.1 and QIIME2 2020.8)

### Install required libraries

```R
# install.packages("ggplot2")
# install.packages("gridExtra")
# install.packages("ggdendro")
# install.packages("ape")
# install.packages("RColorBrewer")
# install.packages("UpSetR")
# install.packages("vegan")
# install.packages("ggfortify")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("philr")
# BiocManager::install("phyloseq")
# BiocManager::install("ggtree")
# BiocManager::install("ALDEx2")
# devtools::install_github('reptalex/phylofactor')
```

### Load required libraries

```R
library(ape)
library(phyloseq)
library(philr)
library(ggdendro)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(ggfortify)
library(UpSetR)
library(plyr)
library(vegan)
library(phylofactor)
library(ggtree)
library(ALDEx2)
```

### Load data into R

Load metadata (map.txt), sequence table with no taxonomy information (sequence_table.16s.filtered.txt), taxonomy assignments (taxonomy_L7.txt), and representative sequence tree (rep_set.root.tre)

```R
rawmetadata <- read.table("map.txt", sep="\t", header=T)
seqtab.filtered <- read.table("../01-raw_data_processing/sequence_table.16s.filtered.txt", header=T, row.names=1)
system("sed 's/;/\t/g' ../01-raw_data_processing/taxonomy_L7.txt > tax_for_phyloseq.txt")
taxa <- read.table("tax_for_phyloseq.txt", header=F, sep="\t", row.names=1)
tree <- read.tree("../01-raw_data_processing/rep_set.root.tre")
```

Check to see if same samples in metadata and sequence table

```R
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID)
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))
notinmeta
notinraw
```

```text
character(0)
character(0)
```

Now create a phyloseq object from different files

```R
rownames(rawmetadata) <- rawmetadata$SampleID
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=F), sample_data(rawmetadata), tax_table(as.matrix(taxa)), tree)
ps.dat
```
```text
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2989 taxa and 61 samples ]
sample_data() Sample Data:       [ 61 samples by 17 sample variables ]
tax_table()   Taxonomy Table:    [ 2989 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 2989 tips and 2988 internal nodes ]
```

### PHILR transformation

```R
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
# [1] TRUE
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
# [1] TRUE
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
```

### Heirarchical cluster dendrogram from transformed data

```R
# system("mkdir imgs")
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,5), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
write.table(df2, "philr_cluster.txt", quote=F, sep="\t", col.names=NA)
hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
cols <- brewer.pal(6, "Set2")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$Season))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
png("imgs/philr_dendrogram_season.png")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()
```

### PCA of PHILR distances

```R
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
png("imgs/philr_screeplot.png")
screeplot(pca)
dev.off()
png("imgs/pca_season.png")
autoplot(pca, data=rawmetadata, colour="Season") + theme_minimal() + xlim(c(-0.25, 0.31)) + ylim(c(-0.25, 0.31))
dev.off()
png("imgs/pca_matrix.png")
autoplot(pca, data=rawmetadata, colour="Matrix") + theme_minimal() + xlim(c(-0.25, 0.31)) + ylim(c(-0.25, 0.31))
dev.off()
png("imgs/pca_temp.png")
autoplot(pca, data=rawmetadata, colour="Temperature_C") + theme_minimal() + xlim(c(-0.25, 0.31)) + ylim(c(-0.25, 0.31))
dev.off()
```
![screeplot](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/philr_screeplot.png)

![pca season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_season.png)
![pca matrix](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_matrix.png)
![pca temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_temp.png)

### Upset plot

How many ASVs are shared between groups?

```R
map <- as.matrix(read.table("map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merge$Season), FUN=sum) # first by season
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table -- ignore warnining message, still works
agg[agg>1] <- 1
#transpose again
agg <- data.frame(t(agg[,-1]))
png("imgs/upset_seqson.png")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs per species", mb.ratio = c(0.55, 0.45))
dev.off()
agg <- aggregate(merged[,2:n], by=list(merge$Insects), FUN=sum) # by insect status
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table -- ignore warnining message, still works
agg[agg>1] <- 1
#transpose again
agg <- data.frame(t(agg[,-1]))
png("imgs/upset_insects.png")
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs per species", mb.ratio = c(0.55, 0.45))
dev.off()
```
![upset season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/upset_seqson.png)
![upset insects](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/upset_insects.png)

### Taxonomic summary

Abundance and taxonomic composition

```R
png("imgs/abundance_tax_barplot.png")
plot_bar(ps.dat, "Season", fill=rank_names(ps.dat)[2])
dev.off()
```

![ab tax barpot](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/abundance_tax_barplot.png)


Hard to see, collapse low abundant phyla into "other" category

```R
physeq.2 <- filter_taxa(ps.dat, function(x) mean(x) > 0.1, TRUE) # remove low freq ASVs
physeq.3 <- transform_sample_counts(physeq.2, function(x) x/sum(x)) # get relative abundance
glom <- tax_glom(physeq.3, taxrank=rank_names(physeq.3)[2]) # collapse at phylum level
data <- psmelt(glom) # create dataframe from phyloseq object
data$V3 <- as.character(data$V3) # convert to character
data$V3[data$Abundance < 0.01] <- "< 1% abund" # rename low freq phyla
medians <- ddply(data, ~V3, function(x) c(median=median(x$Abundance)))
medians
```

```text
               V3     median
1      < 1% abund 0.00000000
2  Actinobacteria 0.12674612
3 Armatimonadetes 0.04565219
4   Bacteroidetes 0.07443391
5 Deferribacteres 0.01040156
6      Firmicutes 0.48802825
7  Proteobacteria 0.22158059
8     Tenericutes 0.01152360
```

Plot

```R
png("imgs/taxonomy_barchart.png")
ggplot(data, aes(x=SampleID, y=Abundance, fill=V3)) + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink")) + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
dev.off()
```

![ab tax barpot](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/taxonomy_barchart.png)

Alpha diversity

```R
png("imgs/adiv_allsamp.png")
plot_richness(ps.dat, measures=c("Observed", "Shannon"), color="Season") + theme_minimal()
dev.off()
png("imgs/adiv_insect_season.png")
plot_richness(ps.dat, x="Insects", color="Season", measures=c("Observed", "Shannon")) + theme_minimal()
dev.off()
```

![adiv all](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/adiv_allsamp.png)
![adiv cat](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/adiv_insect_season.png)

### PERMANOVA

Is there a difference in microbial diversity across samples by some metadata category?

```R
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ Season, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ Season, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Season     2    2290.8  1145.4  7.0658 0.19591  0.001 ***
Residuals 58    9401.9   162.1         0.80409
Total     60   11692.7                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ Matrix, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ Matrix, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Matrix     2     196.2  98.107 0.49495 0.01678  0.964
Residuals 58   11496.5 198.215         0.98322
Total     60   11692.7                 1.00000
```

```R
adonis(philr.dist ~ Insects, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ Insects, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Insects    3    2976.4  992.15  6.4881 0.25456  0.001 ***
Residuals 57    8716.3  152.92         0.74544
Total     60   11692.7                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ Temp_group, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ Temp_group, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Temp_group  4    3141.0  785.26  5.1422 0.26863  0.001 ***
Residuals  56    8551.7  152.71         0.73137
Total      60   11692.7                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Phylofactor

Differentially abundant taxa between groups

```R
OTUTable <- as.matrix(t(seqtab.filtered))
filt.list <- colnames(OTUTable)
filt.list <- filt.list[-1:-3] # remove blanks
filtmap <- rawmetadata[rawmetadata$SampleID %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$SampleID),]
filtmap$Season <- droplevels(filtmap$Season) # drop blank level
x <- as.factor(filtmap$Season) 
tree <- phy_tree(philr.dat)
tax <- read.table("tax_for_phyloseq.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
OTUTable <- OTUTable[,filt.list] # filter out blanks
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
png("imgs/phylofactor_tree.png")
gtree$ggplot + geom_tiplab()
dev.off()
```

By temperature group

```R
OTUTable <- as.matrix(t(seqtab.filtered))
filt.list <- colnames(OTUTable)
filt.list <- filt.list[-1:-3] # remove blanks
filtmap <- rawmetadata[rawmetadata$SampleID %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$SampleID),]
filtmap$Season <- droplevels(filtmap$Temp_group) # drop blank level
x <- as.factor(filtmap$Temp_group) 
tree <- phy_tree(philr.dat)
tax <- read.table("tax_for_phyloseq.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
OTUTable <- OTUTable[,filt.list] # filter out blanks
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
png("imgs/phylofactor_tree_tempG.png")
gtree$ggplot + geom_tiplab()
dev.off()
pdf("imgs/phylofactor_tree_tempG.pdf")
gtree$ggplot + geom_tiplab()
dev.off()
```

![phylo tree](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/phylofactor_tree.png)
![phylo tree tg](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/phylofactor_tree_tempG.png)


Boxplots and significance levels for each factor

Factor 1:

```R
y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor1_boxp.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
wilcox.test(dat[dat$V1 == "summer",]$V2, dat[dat$V1 == "winter",]$V2)
```

```text
	Wilcoxon rank sum test

data:  dat[dat$V1 == "summer", ]$V2 and dat[dat$V1 == "winter", ]$V2
W = 145, p-value = 1.425e-05
alternative hypothesis: true location shift is not equal to 0
```

Factor 2:

```R
y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor2_boxp.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2') + ylim(c(-3.5,9.5))
dev.off()
wilcox.test(dat[dat$V1 == "summer",]$V2, dat[dat$V1 == "winter",]$V2)
```

```text
	Wilcoxon rank sum test

data:  dat[dat$V1 == "summer", ]$V2 and dat[dat$V1 == "winter", ]$V2
W = 692, p-value = 2.125e-06
alternative hypothesis: true location shift is not equal to 0
```

Factor 3:

```R
y <- t(PF$basis[,3]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor3_boxp.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3') + ylim(c(-3.5,9.5))
dev.off()
wilcox.test(dat[dat$V1 == "summer",]$V2, dat[dat$V1 == "winter",]$V2)
```

```text
	Wilcoxon rank sum test

data:  dat[dat$V1 == "summer", ]$V2 and dat[dat$V1 == "winter", ]$V2
W = 654, p-value = 5.661e-05
alternative hypothesis: true location shift is not equal to 0
```

![factor1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor1_boxp.png)
![factor2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor2_boxp.png)
![factor3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor3_boxp.png)

### Differential abundance 

First need to format for qiime2 (transpose sequence table beforehand, columns = samples, rows = ASV IDs)

```R
system("biom convert -i ../01-raw_data_processing/sequence_table.16s.filtered.tr.txt -o sequence_table.16s.filtered.biom --table-type="OTU table" --to-hdf5")
system("biom summarize-table -i sequence_table.16s.filtered.biom")
```
```text
Num samples: 61
Num observations: 2,989
Total count: 7,928,033
Table density (fraction of non-zero values): 0.047

Counts/sample summary:
 Min: 2,394.000
 Max: 328,813.000
 Median: 119,915.000
 Mean: 129,967.754
 Std. dev.: 66,955.954
 Sample Metadata Categories: None provided
 Observation Metadata Categories: None provided

Counts/sample detail:
NegCtrl: 2,394.000
S48A: 16,000.000
W07A: 44,641.000
W12E: 47,275.000
W11A: 50,547.000
S13E: 52,958.000
S07E: 59,115.000
S24A: 61,119.000
S33E: 64,722.000
W04E: 66,235.000
S09E: 67,009.000
W14E: 71,210.000
S31E: 77,538.000
S25E: 80,651.000
W24E: 80,996.000
W03A: 83,347.000
W13A: 85,035.000
S08A: 87,589.000
S10A: 92,886.000
S15E: 95,415.000
W28E: 97,112.000
S32A: 99,606.000
S49E: 100,502.000
Negctrl: 101,417.000
W30E: 107,583.000
S06A: 109,901.000
S16A: 110,355.000
S12A: 110,985.000
W09A: 116,494.000
S11E: 117,369.000
S26A: 119,915.000
W27A: 122,649.000
S04A: 123,121.000
W17A: 124,667.000
W20E: 126,772.000
W16E: 129,583.000
W31A: 131,180.000
S02E: 132,168.000
W26E: 135,470.000
W06E: 135,878.000
W29A: 138,578.000
W25A: 147,346.000
W10E: 147,675.000
S27E: 152,616.000
S19E: 162,642.000
S05E: 163,124.000
S30A: 166,900.000
W23A: 173,031.000
S20A: 182,490.000
S14A: 189,895.000
W32E: 190,175.000
S21E: 201,080.000
S17E: 210,433.000
S01A: 211,940.000
S18A: 219,883.000
W15A: 229,232.000
S37E: 249,545.000
Blank: 249,745.000
S35E: 278,400.000
S36A: 295,081.000
S34A: 328,813.000
```

```R
system("qiime tools import --input-path sequence_table.16s.filtered.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path sequence_table.16s.filtered.qza")
```

Now can run ALDEx plugin through qiime

```R
system("qiime feature-table filter-samples --i-table sequence_table.16s.filtered.qza --m-metadata-file map.txt --p-where "[Sample-type]='swab'" --o-filtered-table swab-feature-table.qza")
system("qiime aldex2 aldex2 --i-table swab-feature-table.qza --m-metadata-file map.txt --m-metadata-column Season --output-dir season_aldex")
system("qiime aldex2 effect-plot --i-table season_aldex/differentials.qza --o-visualization season_aldex/season")
system("qiime tools view season_aldex/season.qzv")
```

![aldex season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/season_aldex/effect_plot.png)


```R
system("qiime aldex2 extract-differences --i-table season_aldex/differentials.qza --o-differentials season_aldex/season --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0")
system("qiime tools export --input-path season_aldex/season.qza --output-path season_aldex/")
system("awk '{print $1}' season_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > season_aldex/differentials.taxonomy.txt")
system("head season_aldex/differentials.taxonomy.txt")
```

```text
ASV4	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV5	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV6	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV7	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV10	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV16	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Burkholderiaceae	Paraburkholderia	Paraburkholderia_unknown
ASV17	Bacteria	Firmicutes	Clostridia	Clostridiales	Peptostreptococcaceae	Peptostreptococcaceae_unknown	Peptostreptococcaceae_unknown
ASV20	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV22	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV23	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
```

Low vs high temperature -- first need to filter out the na sample from mapping file (delete in excel) and then filter from qza

```R
system("qiime feature-table filter-samples --i-table swab-feature-table.qza --m-metadata-file map.filt.txt --o-filtered-table swab-feature-table.filt.qza")
system("qiime aldex2 aldex2 --i-table swab-feature-table.filt.qza --m-metadata-file map.filt.txt --m-metadata-column Temp_group_binary --output-dir temp-low-hi_aldex")
system("qiime aldex2 effect-plot --i-table temp-low-hi_aldex/differentials.qza --o-visualization temp-low-hi_aldex/temp-low-hi")
system("qiime tools view temp-low-hi_aldex/temp-low-hi.qzv")
```

![aldex season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/temp-low-hi_aldex/effect_plot.png)

```R
system("qiime aldex2 extract-differences --i-table temp-low-hi_aldex/differentials.qza --o-differentials temp-low-hi_aldex/temp-low-hi --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0")
system("qiime tools export --input-path temp-low-hi_aldex/temp-low-hi.qza --output-path temp-low-hi_aldex/")
system("awk '{print $1}' temp-low-hi_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > temp-low-hi_aldex/differentials.taxonomy.txt")
system("head temp-low-hi_aldex/differentials.taxonomy.txt")
```

```text
ASV4	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV5	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV6	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV10	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV16	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Burkholderiaceae	Paraburkholderia	Paraburkholderia_unknown
ASV20	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV22	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV23	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV26	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV32	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Morganellaceae	Providencia	Providencia_unknown
```

Clostridium seems to be very high in high temp, very low in low, plot these along temperature gradient. First merge with metadata and taxonomy.

```R
merged1 <- merge(seqtab.filtered, taxa, by=0)
merged.seqtab <- merge(merged1, taxa, by=0)





