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
# install.packages("ggrepel")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("philr")
# BiocManager::install("phyloseq")
# BiocManager::install("ggtree")
# BiocManager::install("ALDEx2")
# BiocManager::install("microbiome")
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
library(microbiome)
library(randomForest)
library(rfUtilities)
library(tidyverse)
library(ggrepel)
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
[1] "BlankE"  "NegCtrl"
```

Now create a phyloseq object from different files

```R
rownames(rawmetadata) <- rawmetadata$SampleID
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=F), sample_data(rawmetadata), tax_table(as.matrix(taxa)), tree)
ps.dat
```
```text
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 2879 taxa and 58 samples ]
sample_data() Sample Data:       [ 58 samples by 31 sample variables ]
tax_table()   Taxonomy Table:    [ 2879 taxa by 7 taxonomic ranks ]
phy_tree()    Phylogenetic Tree: [ 2879 tips and 2878 internal nodes ]
```

### PHILR transformation

```R
ps.dat.nocont <- subset_samples(ps.dat, Sample.type=="swab")
philr.dat <- transform_sample_counts(ps.dat.nocont, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
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
tip_labels <- as.vector(dend_data$labels$label)
png("imgs/philr_dendrogram_season.png")
ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
dev.off()
pdf("imgs/philr_dendrogram_season.pdf")
ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
dev.off()
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
```

### PCA of PHILR distances

```R
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))
png("imgs/philr_screeplot.png")
screeplot(pca)
dev.off()
png("imgs/pca_season.png")
autoplot(pca, data=sample_data(ps.dat.nocont), colour="Season") + theme_minimal() 
dev.off()
png("imgs/pca_matrix.png")
autoplot(pca, data=sample_data(ps.dat.nocont), colour="Matrix") + theme_minimal() 
dev.off()
png("imgs/pca_temp.png")
autoplot(pca, data=sample_data(ps.dat.nocont), colour="Temperature_C") + theme_minimal() 
dev.off()
png("imgs/pca_bmi.png")
autoplot(pca, data=sample_data(ps.dat.nocont), colour="BMI_classification") + theme_minimal() 
dev.off()
png("imgs/pca_sex.png")
autoplot(pca, data=sample_data(ps.dat.nocont), colour="Sex") + theme_minimal() 
dev.off()
png("imgs/pca_bodyID.png")
autoplot(pca, data=sample_data(ps.dat.nocont)) + theme_minimal() + geom_label_repel(aes(label=Body_id, colour=Body_id), box.padding=0.35, point.padding=0.5, segment.color='grey50', max.overlaps=15, cex=1.5)
dev.off()
pdf("imgs/pca_bodyID.pdf")
autoplot(pca, data=sample_data(ps.dat.nocont)) + theme_minimal() + geom_label_repel(aes(label=Body_id, colour=Body_id), box.padding=0.35, point.padding=0.5, segment.color='grey50', max.overlaps=15, cex=1.5)
dev.off()
```
![screeplot](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/philr_screeplot.png)
![pca season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_season.png)
![pca matrix](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_matrix.png)
![pca temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_temp.png)
![pca bmi](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_bmi.png)
![pca sex](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_sex.png)
![pca bodyID](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_bodyID.png)


Colored by surface temperature?

```R
temp <- sample_data(ps.dat.nocont)
temp$Temperature_C <- as.numeric(as.character(temp$Temperature_C))
png("imgs/pca_temperature_cont.png")
autoplot(pca, data=temp, colour="Temperature_C") + theme_minimal() + scale_color_gradient(low="blue",high="red")
dev.off()
```

![pca temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/pca_temperature_cont.png)

### Upset plot

How many ASVs are shared between groups?

```R
map <- as.matrix(read.table("map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merged$Season), FUN=sum) # first by season
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
agg <- aggregate(merged[,2:n], by=list(merged$Insects), FUN=sum) # by insect status
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
physeq.2 <- filter_taxa(ps.dat.nocont, function(x) mean(x) > 0.1, TRUE) # remove low freq ASVs
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
2  Actinobacteria 0.11948640
3 Armatimonadetes 0.04146412
4   Bacteroidetes 0.07630162
5 Deferribacteres 0.01073772
6      Firmicutes 0.52327533
7  Proteobacteria 0.19213305
8     Tenericutes 0.01108594
```

Plot

```R
data$SampleID <- factor(data$SampleID, levels=unique(data$SampleID))
png("imgs/taxonomy_barchart.png")
ggplot(data, aes(x=SampleID, y=Abundance, fill=V3)) + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink")) + theme_minimal() + theme(axis.text.x = element_text(angle = 90))
dev.off()
```

![ab tax barpot](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/taxonomy_barchart.png)

Alpha diversity

```R
png("imgs/adiv_allsamp.png")
plot_richness(ps.dat.nocont, measures=c("Observed", "Shannon"), color="Season") + theme_minimal()
dev.off()
png("imgs/adiv_insect_season.png")
plot_richness(ps.dat.nocont, x="Insects", color="Season", measures=c("Observed", "Shannon")) + theme_minimal()
dev.off()
adiv <- estimate_richness(ps.dat.nocont)
wilcox.test(adiv[grepl("W", rownames(adiv)),]$Observed, adiv[grepl("S", rownames(adiv)),]$Observed)
```

```text
	Wilcoxon rank sum test with continuity correction

data:  adiv[grepl("W", rownames(adiv)), ]$Observed and adiv[grepl("S", rownames(adiv)), ]$Observed
W = 407, p-value = 0.9937
alternative hypothesis: true location shift is not equal to 0
```

![adiv all](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/adiv_allsamp.png)
![adiv cat](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/adiv_insect_season.png)

### PERMANOVA

Is there a difference in microbial diversity across samples by some metadata category?

```R
metadata <- as(sample_data(ps.dat.nocont), "data.frame")
adonis(philr.dist ~ Season, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ Season, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Season     1    1827.0  1827.0  12.076 0.17739  0.001 ***
Residuals 56    8472.7   151.3         0.82261
Total     57   10299.8                 1.00000
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
Matrix     1      59.8  59.836 0.32723 0.00581  0.991
Residuals 56   10239.9 182.856         0.99419
Total     57   10299.8                 1.00000
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
Insects    2    2410.8 1205.38  8.4035 0.23406  0.001 ***
Residuals 55    7889.0  143.44         0.76594
Total     57   10299.8                 1.00000
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
Temp_group  4    2804.0  700.99  4.9564 0.27224  0.001 ***
Residuals  53    7495.8  141.43         0.72776
Total      57   10299.8                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ Sex, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ Sex, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sex        1     885.2  885.24  5.2656 0.08595  0.002 **
Residuals 56    9414.5  168.12         0.91405
Total     57   10299.8                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist ~ BMI_classification, data=metadata)
```

```text
Call:
adonis(formula = philr.dist ~ BMI_classification, data = metadata)

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
BMI_classification  4    2187.3  546.82  3.5725 0.21236  0.001 ***
Residuals          53    8112.5  153.07         0.78764
Total              57   10299.8                 1.00000
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Effect of metadata within season

```R
# filter on season
ps.dat.summer <- subset_samples(ps.dat, Season=="summer")
ps.dat.winter <- subset_samples(ps.dat, Season=="winter")
# get distance matrix for each





philr.dat.summer <- transform_sample_counts(ps.dat.summer, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
phy_tree(philr.dat.summer) <- makeNodeLabel(phy_tree(philr.dat.summer), method="number", prefix="n")
otu.table.summer <- otu_table(philr.dat.summer)
tree.summer <- phy_tree(philr.dat.summer)
metadata.summer <- sample_data(philr.dat.summer)
tax.summer <- tax_table(philr.dat.summer)
philr.t.summer <- philr(otu.table.summer, tree.summer, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist.summer <- dist(philr.t.summer, method="euclidean")












### Phylofactor

Differentially abundant taxa between groups

```R
OTUTable <- as.matrix(t(seqtab.filtered))
remove <- c("BlankE", "NegCtrl")
filtmap <- rawmetadata[!rawmetadata$SampleID %in% remove,]
filtmap$Season <- droplevels(filtmap$Season) # drop blank level
x <- as.factor(filtmap$Season) 
tree <- phy_tree(philr.dat)
tax <- read.table("tax_for_phyloseq.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
OTUTable <- OTUTable[rownames(OTUTable) %in% tree$tip.label,]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
png("imgs/phylofactor_tree.png")
gtree$ggplot + geom_tiplab()
dev.off()
pdf("imgs/phylofactor_tree.pdf")
gtree$ggplot + geom_tiplab()
dev.off()
```

![phylo tree](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/phylofactor_tree.png)

Boxplots and significance levels for each factor

Factor 1:

```R
y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor1_boxp.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor1_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
wilcox.test(dat[dat$V1 == "summer",]$V2, dat[dat$V1 == "winter",]$V2)
```

```text
	Wilcoxon rank sum test

data:  dat[dat$V1 == "summer", ]$V2 and dat[dat$V1 == "winter", ]$V2
W = 52, p-value = 2.72e-10
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
pdf("imgs/factor2_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2') + ylim(c(-3.5,9.5))
dev.off()
wilcox.test(dat[dat$V1 == "summer",]$V2, dat[dat$V1 == "winter",]$V2)
```

```text
	Wilcoxon rank sum test

data:  dat[dat$V1 == "summer", ]$V2 and dat[dat$V1 == "winter", ]$V2
W = 771, p-value = 8.248e-11
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
pdf("imgs/factor3_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3') + ylim(c(-3.5,9.5))
dev.off()
wilcox.test(dat[dat$V1 == "summer",]$V2, dat[dat$V1 == "winter",]$V2)
```

```text
	Wilcoxon rank sum test

data:  dat[dat$V1 == "summer", ]$V2 and dat[dat$V1 == "winter", ]$V2
W = 712, p-value = 2.731e-07
alternative hypothesis: true location shift is not equal to 0
```

![factor1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor1_boxp.png)
![factor2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor2_boxp.png)
![factor3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor3_boxp.png)

By temperature group

```R
OTUTable <- as.matrix(t(seqtab.filtered))
remove <- c("BlankE", "NegCtrl", "S48A", "S49E")
filtmap <- rawmetadata[!rawmetadata$SampleID %in% remove,]
filtmap$Temp_group <- droplevels(filtmap$Temp_group) # drop blank level
x <- as.factor(filtmap$Temp_group) 
tree <- phy_tree(philr.dat)
tax <- read.table("tax_for_phyloseq.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
OTUTable <- OTUTable[rownames(OTUTable) %in% tree$tip.label,]
OTUTable <- OTUTable[, !colnames(OTUTable) %in% remove]
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

![phylo tree tg](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/phylofactor_tree_tempG.png)

Factor 1:

```R
y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor1_boxp_tempG.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor1_boxp_tempG.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
```

Factor 2:

```R
y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor2_boxp_tempG.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor2_boxp_tempG.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2') + ylim(c(-3.5,9.5))
dev.off()
```

Factor 3:

```R
y <- t(PF$basis[,3]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor3_boxp_tempG.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor3_boxp_tempG.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3') + ylim(c(-3.5,9.5))
dev.off()
```

![factor1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor1_boxp_tempG.png)
![factor2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor2_boxp_tempG.png)
![factor3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor3_boxp_tempG.png)

By Insect activity

```R
OTUTable <- as.matrix(t(seqtab.filtered))
remove <- c("BlankE", "NegCtrl")
filtmap <- rawmetadata[!rawmetadata$SampleID %in% remove,]
filtmap$Insects <- droplevels(filtmap$Insects) # drop blank level
x <- as.factor(filtmap$Insects) 
tree <- phy_tree(philr.dat)
tax <- read.table("tax_for_phyloseq.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
OTUTable <- OTUTable[rownames(OTUTable) %in% tree$tip.label,]
OTUTable <- OTUTable[, !colnames(OTUTable) %in% remove]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
png("imgs/phylofactor_tree_insects.png")
gtree$ggplot + geom_tiplab()
dev.off()
pdf("imgs/phylofactor_tree_insects.pdf")
gtree$ggplot + geom_tiplab()
dev.off()
```

![phylo tree tg](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/phylofactor_tree_insects.png)

Factor 1:

```R
y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor1_boxp_insects.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor1_boxp_insects.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()
```

Factor 2:

```R
y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor2_boxp_insects.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor2_boxp_insects.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 2') + ylim(c(-3.5,9.5))
dev.off()
```

Factor 3:

```R
y <- t(PF$basis[,3]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
png("imgs/factor3_boxp_insects.png")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3') + ylim(c(-3.5,9.5))
dev.off()
pdf("imgs/factor3_boxp_insects.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 3') + ylim(c(-3.5,9.5))
dev.off()
```

![factor1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor1_boxp_insects.png)
![factor2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor2_boxp_insects.png)
![factor3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/factor3_boxp_insects.png)

### Differential abundance 

```bash
conda activate qiime2-2020.8
```

First need to format for qiime2 (transpose sequence table beforehand, columns = samples, rows = ASV IDs)

```R
biom convert -i ../01-raw_data_processing/sequence_table.16s.filtered.tr.txt -o sequence_table.16s.filtered.biom --table-type='OTU table' --to-hdf5
biom summarize-table -i sequence_table.16s.filtered.biom
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

```bash
qiime tools import --input-path sequence_table.16s.filtered.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path sequence_table.16s.filtered.qza
```

Now can run ALDEx plugin through qiime

```R
qiime feature-table filter-samples --i-table sequence_table.16s.filtered.qza --m-metadata-file map.txt --p-where "[Sample-type]='swab'" --o-filtered-table swab-feature-table.qza
qiime aldex2 aldex2 --i-table swab-feature-table.qza --m-metadata-file map.txt --m-metadata-column Season --output-dir season_aldex
qiime aldex2 effect-plot --i-table season_aldex/differentials.qza --o-visualization season_aldex/season
qiime tools view season_aldex/season.qzv
```

![aldex season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/season_aldex/effect_plot.png)


```R
qiime aldex2 extract-differences --i-table season_aldex/differentials.qza --o-differentials season_aldex/season --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0
qiime tools export --input-path season_aldex/season.qza --output-path season_aldex/
awk '{print $1}' season_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > season_aldex/differentials.taxonomy.txt
head season_aldex/differentials.taxonomy.txt
```

```text
ASV4	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV5	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV6	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV7	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV10	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV16	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Burkholderiaceae	Paraburkholderia	Paraburkholderia_unknown
ASV17	Bacteria	Firmicutes	Clostridia	Clostridiales	Peptostreptococcaceae	Peptostreptococcaceae_unknown	Peptostreptococcaceae_unknown
ASV20	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV22	Bacteria	Firmicutes	Erysipelotrichi	Erysipelotrichales	Erysipelotrichaceae	Erysipelotrichaceae_unknown	Erysipelotrichaceae_unknown
ASV23	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
```

Low vs high temperature -- first need to filter out the na sample from mapping file (delete in excel) and then filter from qza

```R
qiime feature-table filter-samples --i-table swab-feature-table.qza --m-metadata-file map.filt.txt --o-filtered-table swab-feature-table.filt.qza
qiime aldex2 aldex2 --i-table swab-feature-table.filt.qza --m-metadata-file map.filt.txt --m-metadata-column Temp_group_binary --output-dir temp-low-hi_aldex
qiime aldex2 effect-plot --i-table temp-low-hi_aldex/differentials.qza --o-visualization temp-low-hi_aldex/temp-low-hi
qiime tools view temp-low-hi_aldex/temp-low-hi.qzv
```

![aldex season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/temp-low-hi_aldex/effect_plot.png)

```R
qiime aldex2 extract-differences --i-table temp-low-hi_aldex/differentials.qza --o-differentials temp-low-hi_aldex/temp-low-hi --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0
qiime tools export --input-path temp-low-hi_aldex/temp-low-hi.qza --output-path temp-low-hi_aldex/
awk '{print $1}' temp-low-hi_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > temp-low-hi_aldex/differentials.taxonomy.txt
head temp-low-hi_aldex/differentials.taxonomy.txt
```

```text
ASV4	Bacteria	Actinobacteria	Actinobacteria_c	Corynebacteriales	Corynebacteriaceae	Corynebacterium	Corynebacterium_unknown
ASV5	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV6	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV10	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV16	Bacteria	Proteobacteria	Betaproteobacteria	Burkholderiales	Burkholderiaceae	Paraburkholderia	Paraburkholderia_unknown
ASV20	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV22	Bacteria	Firmicutes	Erysipelotrichi	Erysipelotrichales	Erysipelotrichaceae	Erysipelotrichaceae_unknown	Erysipelotrichaceae_unknown
ASV23	Bacteria	Firmicutes	Bacilli	Bacillales	Planococcaceae	Planococcaceae_unknown	Planococcaceae_unknown
ASV26	Bacteria	Firmicutes	Clostridia	Clostridiales	Clostridiaceae	Clostridium	Clostridium_unknown
ASV32	Bacteria	Proteobacteria	Gammaproteobacteria	Enterobacterales	Morganellaceae	Providencia	Providencia_unknown
```

### Beta dispersion

```R
dispr <- vegan::betadisper(philr.dist, phyloseq::sample_data(ps.dat.nocont)$Season)
dispr
```

```text
	Homogeneity of multivariate dispersions

Call: vegan::betadisper(d = philr.dist, group =
phyloseq::sample_data(ps.dat.nocont)$Season)

No. of Positive Eigenvalues: 57
No. of Negative Eigenvalues: 0

Average distance to median:
summer winter
 11.49  11.71

Eigenvalues for PCoA axes:
(Showing 8 of 57 eigenvalues)
 PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8
2884.0 2465.0 1063.1  546.6  439.1  374.7  319.5  245.1
```

```R
permutest(dispr)
```

```text
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq Mean Sq      F N.Perm Pr(>F)
Groups     1   0.67  0.6732 0.0394    999  0.872
Residuals 56 956.64 17.0828
```

```R
dispr <- vegan::betadisper(philr.dist, phyloseq::sample_data(ps.dat.nocont)$Sex)
dispr
```

```text
	Homogeneity of multivariate dispersions

Call: vegan::betadisper(d = philr.dist, group =
phyloseq::sample_data(ps.dat.nocont)$Sex)

No. of Positive Eigenvalues: 57
No. of Negative Eigenvalues: 0

Average distance to median:
Female   Male
 13.09  11.77

Eigenvalues for PCoA axes:
(Showing 8 of 57 eigenvalues)
 PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8
2884.0 2465.0 1063.1  546.6  439.1  374.7  319.5  245.1
```

```R
permutest(dispr)
```

```text
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq Mean Sq      F N.Perm Pr(>F)
Groups     1  24.74  24.745 2.2261    999  0.128
Residuals 56 622.49  11.116
```

```R
dispr <- vegan::betadisper(philr.dist, phyloseq::sample_data(ps.dat.nocont)$BMI_classification)
dispr
```

```text
	Homogeneity of multivariate dispersions

Call: vegan::betadisper(d = philr.dist, group =
phyloseq::sample_data(ps.dat.nocont)$BMI_classification)

No. of Positive Eigenvalues: 57
No. of Negative Eigenvalues: 0

Average distance to median:
         na      normal       obese  overweight underweight
      1.416      13.722      10.054      10.800       2.306

Eigenvalues for PCoA axes:
(Showing 8 of 57 eigenvalues)
 PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8
2884.0 2465.0 1063.1  546.6  439.1  374.7  319.5  245.1
```

```R
permutest(dispr)
```

```text
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq Mean Sq      F N.Perm Pr(>F)
Groups     4 705.51 176.378 27.396    999  0.001 ***
Residuals 53 341.22   6.438
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

NOTE: This is probably driven by higher number of "normal" individuals, only one underweight individual


```R
png("imgs/betadispr_bmi_type.png")
boxplot(dispr)
dev.off()
```

![bdisp_bmi](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/betadispr_bmi_type.png)

```R
dispr <- vegan::betadisper(philr.dist, phyloseq::sample_data(ps.dat.nocont)$Temp_group)
dispr
```

```text
	Homogeneity of multivariate dispersions

Call: vegan::betadisper(d = philr.dist, group =
phyloseq::sample_data(ps.dat.nocont)$Temp_group)

No. of Positive Eigenvalues: 57
No. of Negative Eigenvalues: 0

Average distance to median:
   20C    30C    40C    50C     na
12.662 10.914 11.008  2.713  5.931

Eigenvalues for PCoA axes:
(Showing 8 of 57 eigenvalues)
 PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8
2884.0 2465.0 1063.1  546.6  439.1  374.7  319.5  245.1
```

```R
permutest(dispr)
```

```text
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq Mean Sq      F N.Perm Pr(>F)
Groups     4 206.78  51.696 3.4759    999  0.024 *
Residuals 53 788.25  14.873
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
png("imgs/betadispr_temp_type.png")
boxplot(dispr)
dev.off()
```

![bdisp_bmi](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/betadispr_temp_type.png)

```R
dispr <- vegan::betadisper(philr.dist, phyloseq::sample_data(ps.dat.nocont)$Insects)
dispr
```

```text
	Homogeneity of multivariate dispersions

Call: vegan::betadisper(d = philr.dist, group =
phyloseq::sample_data(ps.dat.nocont)$Insects)

No. of Positive Eigenvalues: 57
No. of Negative Eigenvalues: 0

Average distance to median:
     n      p      y
14.024  7.980  9.652

Eigenvalues for PCoA axes:
(Showing 8 of 57 eigenvalues)
 PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8
2884.0 2465.0 1063.1  546.6  439.1  374.7  319.5  245.1
```

```R
permutest(dispr)
```

```text
Permutation test for homogeneity of multivariate dispersions
Permutation: free
Number of permutations: 999

Response: Distances
          Df Sum Sq Mean Sq     F N.Perm Pr(>F)
Groups     2 352.55 176.277 25.93    999  0.001 ***
Residuals 55 373.89   6.798
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
png("imgs/betadispr_insect_type.png")
boxplot(dispr)
dev.off()
```

![bdisp_bmi](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/betadispr_insect_type.png)







### OLD BELOW

Clostridium seems to be very high in high temp, very low in low, plot these along temperature gradient. First merge with taxonomy. Open in excel and get values from clostridium.

```R
seqtab.nochim <- read.table("../01-raw_data_processing/sequence_table.16s.filtered.tr.txt", header=T, row.names=1)
taxa <- read.table("../01-raw_data_processing/taxonomy_L7.txt", header=F, row.names=1)
merged <- merge(seqtab.nochim, taxa, by=0)
write.table(data.frame("row_names"=rownames(merged),merged),"sequence_taxonomy_table.16s.merged.txt", row.names=FALSE, quote=F, sep="\t")
```

Test plot differentially abundant ASVs by temperature

```R
test <- otu_table(philr.dat)[,"ASV4"]
test.m <- merge(rawmetadata, test, by=0)
test.m$Temperature_C <- as.numeric(as.character(test.m$Temperature_C))
test.m <- test.m[!is.na(test.m$Temperature_C),]
test.m$Temp_C_bin <- cut(test.m$Temperature_C, breaks=10)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}
df2 <- data_summary(test.m, varname="ASV4", groupnames=c("Temperature_C"))
png(paste("imgs/", "test_ASV4.png", sep=""))
ggplot(df2, aes(x=as.factor(Temperature_C), y=ASV4)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("IRL Transformed Read Counts") + geom_errorbar(aes(ymin=ASV4-sd, ymax=ASV4+sd), width=.2) + geom_point(test.m, mapping=aes(x=as.factor(Temperature_C), y=ASV4)) + geom_jitter()
dev.off()
```

![asv4 IRL counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/test_ASV4.png)

Now do a real plot based off of differentially abundant clades in phylofactor by temp group

```R
# get list of asvs to filter data
factor1 <- c("ASV5", "ASV19", "ASV2", "ASV80", "ASV23", "ASV85", "ASV212", "ASV122")
factor2 <- c("ASV7", "ASV125", "ASV20", "ASV335", "ASV471", "ASV10", "ASV194", "ASV6", "ASV93", "ASV26", "ASV35", "ASV38", "ASV34")
factor3 <- c("ASV14", "ASV9", "ASV1", "ASV56", "ASV25", "ASV3", "ASV67")
# convert philr.dat transformed otu table to dataframe
df <- otu_table(philr.dat)
df <- as.data.frame(df)
# filter data
df.fac1 <- df[,which((names(df) %in% factor1)==TRUE)]
df.fac2 <- df[,which((names(df) %in% factor2)==TRUE)]
df.fac3 <- df[,which((names(df) %in% factor3)==TRUE)]
# sum across all rows
df.fac1 <- as.data.frame(rowSums(df.fac1))
df.fac2 <- as.data.frame(rowSums(df.fac2))
df.fac3 <- as.data.frame(rowSums(df.fac3))
# merge data with metadata
df.fac1.m <- merge(rawmetadata, df.fac1, by=0)
df.fac2.m <- merge(rawmetadata, df.fac2, by=0)
df.fac3.m <- merge(rawmetadata, df.fac3, by=0)
# prep temperature category
df.fac1.m$Temperature_C <- as.numeric(as.character(df.fac1.m$Temperature_C))
df.fac2.m$Temperature_C <- as.numeric(as.character(df.fac2.m$Temperature_C))
df.fac3.m$Temperature_C <- as.numeric(as.character(df.fac3.m$Temperature_C))
df.fac1.m <- df.fac1.m[!is.na(df.fac1.m$Temperature_C),]
df.fac2.m <- df.fac2.m[!is.na(df.fac2.m$Temperature_C),]
df.fac3.m <- df.fac3.m[!is.na(df.fac3.m$Temperature_C),]
df.fac1.m$Temp_C_bin <- cut(df.fac1.m$Temperature_C, breaks=10)
df.fac2.m$Temp_C_bin <- cut(df.fac2.m$Temperature_C, breaks=10)
df.fac3.m$Temp_C_bin <- cut(df.fac3.m$Temperature_C, breaks=10)
# summarize data into bins
df.fac1.sum <- data_summary(df.fac1.m, varname="rowSums(df.fac1)", groupnames=c("Temperature_C"))
df.fac2.sum <- data_summary(df.fac2.m, varname="rowSums(df.fac2)", groupnames=c("Temperature_C"))
df.fac3.sum <- data_summary(df.fac3.m, varname="rowSums(df.fac3)", groupnames=c("Temperature_C"))
# clean up
colnames(df.fac1.sum) <- c("Temperature_C", "factor", "sd")
colnames(df.fac2.sum) <- c("Temperature_C", "factor", "sd")
colnames(df.fac3.sum) <- c("Temperature_C", "factor", "sd")
# plot
png(paste("imgs/", "fact1_barplot_tempG.png", sep=""))
ggplot(df.fac1.sum, aes(x=as.factor(Temperature_C), y=factor)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("ILR Transformed Read Counts") + geom_errorbar(aes(ymin=factor-sd, ymax=factor+sd), width=.2) + geom_point(df.fac1.sum, mapping=aes(x=as.factor(Temperature_C), y=factor)) + geom_jitter()
dev.off()
png(paste("imgs/", "fact2_barplot_tempG.png", sep=""))
ggplot(df.fac2.sum, aes(x=as.factor(Temperature_C), y=factor)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("ILR Transformed Read Counts") + geom_errorbar(aes(ymin=factor-sd, ymax=factor+sd), width=.2) + geom_point(df.fac2.sum, mapping=aes(x=as.factor(Temperature_C), y=factor)) + geom_jitter()
dev.off()
png(paste("imgs/", "fact3_barplot_tempG.png", sep=""))
ggplot(df.fac3.sum, aes(x=as.factor(Temperature_C), y=factor)) + geom_bar(stat="identity", color="black", fill="white") + theme_minimal() + xlab("Temperature C") + ylab("ILR Transformed Read Counts") + geom_errorbar(aes(ymin=factor-sd, ymax=factor+sd), width=.2) + geom_point(df.fac3.sum, mapping=aes(x=as.factor(Temperature_C), y=factor)) + geom_jitter()
dev.off()
```

![fact1 ILR counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/fact1_barplot_tempG.png)
![fact2 ILR counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/fact2_barplot_tempG.png)
![fact3 ILR counts temp](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/fact3_barplot_tempG.png)

### Random forest

```R
otu_table <- read.table("../01-raw_data_processing/sequence_table.16s.filtered.txt", sep="\t", header=T, row.names=1, stringsAsFactors=F, comment.char="")
otu_table <- t(otu_table)
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=T, comment.char="")
metadata <- metadata[metadata$Season %in% c("winter", "summer"),]
metadata$Season <- factor(metadata$Season)
otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.1)
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed), '/')*100
otu_table_scaled <- scale(otu_table_rare_removed_norm, center=T, scale=T)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Season"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_season <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_season
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T) 
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 6.9%
Confusion matrix:
       summer winter class.error
summer     32      2  0.05882353
winter      2     22  0.08333333
```

Lysing matrix

```R
metadata <- metadata[metadata$Matrix %in% c("E", "A"),]
metadata$Matrix <- factor(metadata$Matrix)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Matrix"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_matrix <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_matrix
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T) 
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 55.93%
Confusion matrix:
   A  E class.error
A 12 18   0.6000000
E 15 14   0.5172414
```

Insects

```R
metadata <- metadata[metadata$Insects %in% c("y", "n"),]
metadata$Insects <- factor(metadata$Insects)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Insects"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_insects <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_insects
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T) 
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 8.33%
Confusion matrix:
   n  y class.error
n 21  3  0.12500000
y  1 23  0.04166667
```

Temp group

```R
metadata <- metadata[metadata$Temp_group %in% c("20C", "30C", "40C", "50C"),]
metadata$Temp_group <- factor(metadata$Temp_group)
otu_table_scaled_var <- data.frame(t(otu_table_scaled))
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Temp_group"]
set.seed(151)
otu_table_scaled_var <- otu_table_scaled_var %>% filter(!is.na(var))
rf_tgroup <- randomForest(x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], y=otu_table_scaled_var[, ncol(otu_table_scaled_var)], ntree=10000, importance=T, proximity=T)
rf_tgroup
```

```text
Call:
 randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = T, proximity = T) 
               Type of random forest: classification
                     Number of trees: 10000
No. of variables tried at each split: 14

        OOB estimate of  error rate: 17.39%
Confusion matrix:
    20C 30C 40C class.error
20C   2   5   1  0.75000000
30C   0  27   1  0.03571429
40C   0   1   9  0.10000000
```

What top 5 taxa are most important in the random forest models?

```R
png(paste("imgs/", "rf_season_importance.png", sep=""))
varImpPlot(rf_season)
dev.off()
pdf(paste("imgs/", "rf_season_importance.pdf", sep=""))
varImpPlot(rf_season)
dev.off()
```

![rf_season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/rf_season_importance.png)

Season:
1. ASV6	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium_unknown
2. ASV10	Bacteria;Firmicutes;Clostridia;Clostridiales;Clostridiaceae;Clostridium;Clostridium_unknown
3. ASV2	Bacteria;Firmicutes;Bacilli;Bacillales;Planococcaceae;Sporosarcina;Sporosarcina_unknown
4. ASV4	Bacteria;Actinobacteria;Actinobacteria_c;Corynebacteriales;Corynebacteriaceae;Corynebacterium;Corynebacterium_unknown
5. ASV25	Bacteria;Proteobacteria;Gammaproteobacteria;Ignatzschineria_o;Ignatzschineria_f;Ignatzschineria_f_unknown;Ignatzschineria_f_unknown

```R
png(paste("imgs/", "rf_insects_importance.png", sep=""))
varImpPlot(rf_insects)
dev.off()
pdf(paste("imgs/", "rf_insects_importance.pdf", sep=""))
varImpPlot(rf_insects)
dev.off()
```

![rf insects](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/rf_insects_importance.png)

Insects:
1. ASV3	Bacteria;Proteobacteria;Gammaproteobacteria;Ignatzschineria_o;Ignatzschineria_f;Ignatzschineria_f_unknown;Ignatzschineria_f_unknown
2. ASV1	Bacteria;Proteobacteria;Gammaproteobacteria;Ignatzschineria_o;Ignatzschineria_f;Ignatzschineria_f_unknown;Ignatzschineria_f_unknown
3. ASV73	Bacteria;Firmicutes;Bacilli;Lactobacillales;Lactobacillaceae;Lactobacillus;Lactobacillus_unknown
4. ASV89	Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Morganellaceae;Providencia;Providencia_unknown
5. ASV53	Bacteria;Firmicutes;Bacilli;Lactobacillales;Carnobacteriaceae;Trichococcus;Trichococcus_unknown

Plot for figure

```R
impToPlot.season <- importance(rf_season, scale=F)[,3]
impToPlot.season <- sort(impToPlot.season)
short.imp <- tail(impToPlot.season, 10)
pdf(paste("imgs/", "rf_season_importance.filt.pdf", sep=""))
dotchart(short.imp, xlim=c(0.00, 0.06))
dev.off()
impToPlot.season <- importance(rf_season, scale=F)[,4]
impToPlot.season <- sort(impToPlot.season)
short.imp <- tail(impToPlot.season, 10)
pdf(paste("imgs/", "rf_season_importance.filt.gini.pdf", sep=""))
dotchart(short.imp, xlim=c(0.0, 2.5))
dev.off()
impToPlot.insect <- importance(rf_insects, scale=F)[,3]
impToPlot.insect <- sort(impToPlot.insect)
short.imp <- tail(impToPlot.insect, 10)
pdf(paste("imgs/", "rf_insects_importance.filt.pdf", sep=""))
dotchart(short.imp, xlim=c(0.00, 0.06))
dev.off()
impToPlot.insect <- importance(rf_insects, scale=F)[,4]
impToPlot.insect <- sort(impToPlot.insect)
short.imp <- tail(impToPlot.insect, 10)
pdf(paste("imgs/", "rf_insects_importance.filt.gini.pdf", sep=""))
dotchart(short.imp, xlim=c(0.0, 2.5))
dev.off()
```

### Correlation between temp/days in field and ratio between major phyla 

added a pseudo count of 1 to get around divide by zero issues

```R
png(paste("imgs/", "corr_actino-firmi_temp.png", sep=""))
qplot(as.numeric(as.character(seqtab.phylum$Temperature_C)), seqtab.phylum$actino.firmi, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_actino-proteo_temp.png", sep=""))
qplot(as.numeric(as.character(seqtab.phylum$Temperature_C)), seqtab.phylum$actino.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_firmi-proteo_temp.png", sep=""))
qplot(as.numeric(as.character(seqtab.phylum$Temperature_C)), seqtab.phylum$firmi.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
```

![corr_temp1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-firmi_temp.png)
![corr_temp2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-proteo_temp.png)
![corr_temp3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_firmi-proteo_temp.png)


```R
png(paste("imgs/", "corr_actino-firmi_days.png", sep=""))
qplot(seqtab.phylum$Days_in_Field, seqtab.phylum$actino.firmi, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_actino-proteo_days.png", sep=""))
qplot(seqtab.phylum$Days_in_Field, seqtab.phylum$actino.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
png(paste("imgs/", "corr_firmi-proteo_days.png", sep=""))
qplot(seqtab.phylum$Days_in_Field, seqtab.phylum$firmi.proteo, geom=c("point", "smooth")) + theme_minimal()
dev.off()
```

![corr_days1](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-firmi_days.png)
![corr_days2](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_actino-proteo_days.png)
![corr_days3](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/corr_firmi-proteo_days.png)
