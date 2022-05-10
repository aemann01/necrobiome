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

### Bray curtis biplot

```R
phylum.sum <- tapply(taxa_sums(ps.dat.nocont), tax_table(ps.dat.nocont)[,"V3"], sum, na.rm=T)
top6phyla <- names(sort(phylum.sum, T))[1:5]
ps1 <- prune_taxa((tax_table(ps.dat.nocont)[,"V3"] %in% top6phyla), ps.dat.nocont)
ps1.ord <- ordinate(ps1, "NMDS", "bray")
pdf("biplot.pdf")
plot_ordination(ps1, ps1.ord, type="split", color="V3", label="Body_id", shape="Season") + theme_minimal() + coord_flip() + scale_x_reverse()
dev.off()
png("biplot.png")
plot_ordination(ps1, ps1.ord, type="split", color="V3", label="Body_id", shape="Season") + theme_minimal() + coord_flip() + scale_x_reverse()
dev.off()
```

![biplot](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/biplot.png)

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
metadata.summer <- as(sample_data(ps.dat.summer), "data.frame")
tax.summer <- tax_table(philr.dat.summer)
philr.t.summer <- philr(otu.table.summer, tree.summer, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist.summer <- dist(philr.t.summer, method="euclidean")
# winter
philr.dat.winter <- transform_sample_counts(ps.dat.winter, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
phy_tree(philr.dat.winter) <- makeNodeLabel(phy_tree(philr.dat.winter), method="number", prefix="n")
otu.table.winter <- otu_table(philr.dat.winter)
tree.winter <- phy_tree(philr.dat.winter)
metadata.winter <- as(sample_data(ps.dat.winter), "data.frame")
tax.winter <- tax_table(philr.dat.winter)
philr.t.winter <- philr(otu.table.winter, tree.winter, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist.winter <- dist(philr.t.winter, method="euclidean")
```

Testing insects, temp_group, sex, and bmi_classification within season, only significant results reported below

```R
adonis(philr.dist.summer ~ Insects, data = metadata.summer)
```

```text
Call:
adonis(formula = philr.dist.summer ~ Insects, data = metadata.summer) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Insects    2    2908.4 1454.21  3.2852 0.17488  0.002 **
Residuals 31   13722.2  442.65         0.82512          
Total     33   16630.6                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist.summer ~ Temp_group, data=metadata.summer)
```

```text
Call:
adonis(formula = philr.dist.summer ~ Temp_group, data = metadata.summer) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Temp_group  3    4086.2 1362.08  3.2574 0.24571  0.001 ***
Residuals  30   12544.4  418.15         0.75429           
Total      33   16630.6                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist.summer ~ Sex, data = metadata.summer)
```

```text
Call:
adonis(formula = philr.dist.summer ~ Sex, data = metadata.summer) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sex        1    3362.3  3362.3  8.1092 0.20218  0.001 ***
Residuals 32   13268.3   414.6         0.79782           
Total     33   16630.6                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist.summer ~ BMI_classification, data = metadata.summer)
```

```text
Call:
adonis(formula = philr.dist.summer ~ BMI_classification, data = metadata.summer) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
BMI_classification  2      2943 1471.51  3.3327 0.17696  0.003 **
Residuals          31     13688  441.53         0.82304          
Total              33     16631                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Winter

```R
adonis(philr.dist.winter ~ Insects, data = metadata.winter)
```

```text
Call:
adonis(formula = philr.dist.winter ~ Insects, data = metadata.winter) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Insects    1    3184.3  3184.3   20.79 0.48586  0.001 ***
Residuals 22    3369.7   153.2         0.51414           
Total     23    6554.0                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
adonis(philr.dist.winter ~ Temp_group, data=metadata.winter)
```

```text
Call:
adonis(formula = philr.dist.winter ~ Temp_group, data = metadata.winter) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Temp_group  1    1226.0 1225.96  5.0621 0.18705  0.011 *
Residuals  22    5328.1  242.19         0.81295         
Total      23    6554.0                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

Sex not statistically significant in winter

```R
adonis(philr.dist.winter ~ BMI_classification, data = metadata.winter)
```

```text
Call:
adonis(formula = philr.dist.winter ~ BMI_classification, data = metadata.winter) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

                   Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
BMI_classification  4    2818.0  704.49  3.5827 0.42996  0.004 **
Residuals          19    3736.1  196.64         0.57004          
Total              23    6554.0                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Differential abundance 

Activate qiime2 environment

```bash
conda activate qiime2-2020.2
```

Install ALDEx2 through R

```R
install.packages("BiocManager")
BiocManager::install("ALDEx2")
```

Install qiime ALDEx2 plugin through conda

```bash
conda install -c dgiguere q2-aldex2
```

Transpose biom table

```bash
awk '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' sequence_table.16s.filtered.txt | sed 's/ /\t/g' > sequence_table.16s.filtered.tr.txt
```

```bash
biom convert -i sequence_table.16s.filtered.tr.txt -o sequence_table.16s.filtered.biom --table-type='OTU table' --to-hdf5
biom summarize-table -i sequence_table.16s.filtered.biom
```
```text
Num samples: 58
Num observations: 2,954
Total count: 7,694,223
Table density (fraction of non-zero values): 0.048

Counts/sample summary:
 Min: 16,020.000
 Max: 336,515.000
 Median: 124,689.000
 Mean: 132,659.017
 Std. dev.: 66,606.669
 Sample Metadata Categories: None provided
 Observation Metadata Categories: None provided

Counts/sample detail:
S48A: 16,020.000
W07A: 44,672.000
W12E: 47,493.000
W11A: 50,879.000
S13E: 51,681.000
S07E: 59,190.000
S24A: 61,351.000
S33E: 64,224.000
W04E: 66,186.000
S09E: 69,232.000
W14E: 71,479.000
S31E: 76,573.000
S25E: 81,041.000
W03A: 82,321.000
W24E: 82,708.000
W13A: 84,908.000
S08A: 88,837.000
S10A: 92,993.000
S15E: 96,206.000
S32A: 97,258.000
W28E: 97,368.000
S49E: 100,765.000
S12A: 107,710.000
W30E: 107,912.000
S16A: 109,860.000
S06A: 111,400.000
W09A: 119,388.000
S11E: 119,455.000
W27A: 124,238.000
S26A: 125,140.000
W17A: 125,184.000
W20E: 127,919.000
S04A: 128,248.000
W16E: 130,601.000
S02E: 132,464.000
W31A: 133,023.000
W06E: 134,716.000
W26E: 135,620.000
W29A: 139,534.000
W25A: 147,663.000
W10E: 150,912.000
S27E: 159,327.000
S30A: 167,832.000
S19E: 173,258.000
W23A: 174,604.000
S05E: 175,855.000
S14A: 189,582.000
W32E: 190,243.000
S20A: 194,536.000
S17E: 208,105.000
S01A: 210,705.000
S21E: 216,643.000
W15A: 229,778.000
S18A: 238,232.000
S37E: 251,537.000
S35E: 286,416.000
S36A: 296,683.000
S34A: 336,515.000
```

```bash
qiime tools import --input-path sequence_table.16s.filtered.biom --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path sequence_table.16s.filtered.qza
```

Now can run ALDEx plugin through qiime

Differential abundance by season

```bash
qiime feature-table filter-samples --i-table sequence_table.16s.filtered.qza --m-metadata-file map.txt --p-where "[Sample-type]='swab'" --o-filtered-table swab-feature-table.qza
qiime aldex2 aldex2 --i-table swab-feature-table.qza --m-metadata-file map.txt --m-metadata-column Season --output-dir season_aldex
qiime aldex2 effect-plot --i-table season_aldex/differentials.qza --o-visualization season_aldex/season
qiime aldex2 extract-differences --i-table season_aldex/differentials.qza --o-differentials season_aldex/season --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0
qiime tools export --input-path season_aldex/season.qza --output-path season_aldex/
awk '{print $1}' season_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > season_aldex/differentials.taxonomy.txt
head season_aldex/differentials.taxonomy.txt
```

```text
ASV4    Bacteria    Actinobacteria  Actinobacteria_c    Corynebacteriales   Corynebacteriaceae  Corynebacterium Corynebacterium_unknown
ASV5    Bacteria    Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown
ASV6    Bacteria    Firmicutes  Bacilli Bacillales  Planococcaceae  Planococcaceae_unknown  Planococcaceae_unknown
ASV10   Bacteria    Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown
ASV16   Bacteria    Proteobacteria  Betaproteobacteria  Burkholderiales Burkholderiaceae    Paraburkholderia    Paraburkholderia_unknown
ASV17   Bacteria    Firmicutes  Clostridia  Clostridiales   PeptostreptococcaceaePeptostreptococcaceae_unknown  Peptostreptococcaceae_unknown
ASV20   Bacteria    Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown
ASV21   Bacteria    Actinobacteria  Actinobacteria_c    Corynebacteriales   Corynebacteriaceae  Corynebacterium Corynebacterium_unknown
ASV23   Bacteria    Firmicutes  Bacilli Bacillales  Planococcaceae  Planococcaceae_unknown  Planococcaceae_unknown
ASV26   Bacteria    Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown
```

Differential abundance by sex

```bash
qiime aldex2 aldex2 --i-table swab-feature-table.qza --m-metadata-file map.txt --m-metadata-column Sex --output-dir sex_aldex
qiime aldex2 effect-plot --i-table sex_aldex/differentials.qza --o-visualization sex_aldex/season
qiime aldex2 extract-differences --i-table sex_aldex/differentials.qza --o-differentials sex_aldex/sex --p-sig-threshold 0.1 --p-effect-threshold 0 --p-difference-threshold 0
qiime tools export --input-path sex_aldex/sex.qza --output-path sex_aldex/
awk '{print $1}' sex_aldex/differentials.tsv | grep "ASV" | while read line; do grep -w $line tax_for_phyloseq.txt ; done > sex_aldex/differentials.taxonomy.txt
head sex_aldex/differentials.taxonomy.txt
```

```text
ASV2    Bacteria    Firmicutes  Bacilli Bacillales  Planococcaceae  Sporosarcina    Sporosarcina_unknown
ASV6    Bacteria    Firmicutes  Bacilli Bacillales  Planococcaceae  Planococcaceae_unknown  Planococcaceae_unknown
ASV17   Bacteria    Firmicutes  Clostridia  Clostridiales   PeptostreptococcaceaePeptostreptococcaceae_unknown  Peptostreptococcaceae_unknown
ASV28   Bacteria    Firmicutes  Tissierellia    Tissierellales  Tissierellaceae Tissierellaceae_unknown Tissierellaceae_unknown
ASV33   Bacteria    Firmicutes  Clostridia  Clostridiales   Ruminococcaceae Eubacterium_g23 Eubacterium_g23_unknown
ASV42   Bacteria    Firmicutes  Tissierellia    Tissierellales  Peptoniphilaceae    Peptoniphilaceae_unknown    Peptoniphilaceae_unknown
ASV43   Bacteria    Firmicutes  Tissierellia    Tissierellales  Tissierellaceae Tissierellaceae_unknown Tissierellaceae_unknown
ASV102  Bacteria    Firmicutes  Tissierellia    Tissierellales  Peptoniphilaceae    Peptoniphilaceae_unknown    Peptoniphilaceae_unknown
ASV105  Bacteria    Firmicutes  Erysipelotrichi Erysipelotrichales  Erysipelotrichaceae Erysipelothrix  Erysipelothrix_unknown
ASV174  Bacteria    Firmicutes  Tissierellia    Tissierellales  Tissierellaceae Tissierella Tissierella_unknown
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
Groups     1   0.67  0.6732 0.0394    999   0.84
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
Groups     1  24.74  24.745 2.2261    999   0.13
Residuals 56 622.49  11.116
```

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
Groups     4 206.78  51.696 3.4759    999   0.02 *
Residuals 53 788.25  14.873
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
png("betadispr_temp_type.png")
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
png("betadispr_insect_type.png")
boxplot(dispr)
dev.off()
```

![bdisp_bmi](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/betadispr_insect_type.png)


Abundance of different taxa across temperature

```R
map <- read.table("map.txt", sep="\t", header=TRUE)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
    library(doBy)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # Collapse the data
    formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
    datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)

    # Rename columns
    names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
    names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
    names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
}

tgc <- summarySE(map, measurevar="log10_abund", groupvars=c("genera","Temp_range"))
ggplot(tgc, aes(x=as.factor(Temp_range), y=log10_abund, colour=genera, group=genera)) + geom_line(position=pd) + geom_errorbar(aes(ymin=log10_abund-se, ymax=log10_abund+se), width=.2, position=pd, colour="black") + geom_point(position=pd, size=3, shape=21, fill="white") + theme_minimal() + scale_color_manual(values=c("#66c2a5", "#fc8d62", "#8da0cb"), labels=c("Ignatzschineria", "Proteus", "Providencia"), name="Genus") + xlab("Temperature range (celcius)") + ylab("Log10(abundance)")
```

### Random forest

```R
otu_table <- read.table("sequence_table.16s.filtered.txt", sep="\t", header=T, row.names=1, stringsAsFactors=F, comment.char="")
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

Insects

```R
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=T, comment.char="")
metadata <- metadata[metadata$Insects %in% c("y", "n", "p"),]
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

        OOB estimate of  error rate: 8.62%
Confusion matrix:
   n p  y class.error
n 21 0  3       0.125
p  0 8  2       0.200
y  0 0 24       0.000
```

Temp group

```R
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=T, comment.char="")
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

        OOB estimate of  error rate: 19.64%
Confusion matrix:
    20C 30C 40C 50C class.error
20C   2   5   1   0  0.75000000
30C   0  26   2   0  0.07142857
40C   0   1  17   0  0.05555556
50C   0   2   0   0  1.00000000
```

What top 5 taxa are most important in the random forest models?

```R
png("rf_season_importance.png")
varImpPlot(rf_season)
dev.off()
pdf("rf_season_importance.pdf")
varImpPlot(rf_season)
dev.off()
```

![rf_season](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/rf_season_importance.png)

Season:
1. ASV5 Bacteria    Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown
2. ASV10 Bacteria   Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown	
3. ASV2	Bacteria    Firmicutes  Bacilli Bacillales  Planococcaceae  Sporosarcina    Sporosarcina_unknown
4. ASV4	Bacteria    Actinobacteria  Actinobacteria_c    Corynebacteriales   Corynebacteriaceae  Corynebacterium Corynebacterium_unknown
5. ASV29 Bacteria   Proteobacteria  Gammaproteobacteria Ignatzschineria_o   Ignatzschineria_f   Ignatzschineria_f_unknown   Ignatzschineria_f_unknown

```R
png("rf_insects_importance.png")
varImpPlot(rf_insects)
dev.off()
pdf("rf_insects_importance.pdf")
varImpPlot(rf_insects)
dev.off()
```

![rf insects](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/rf_insects_importance.png)

Insects:
1. ASV3	Bacteria    Proteobacteria  Gammaproteobacteria Ignatzschineria_o   Ignatzschineria_f   Ignatzschineria_f_unknown   Ignatzschineria_f_unknown
2. ASV103 Proteobacteria    Gammaproteobacteria Enterobacterales    Morganellaceae  Providencia Providencia_unknown
3. ASV34 Bacteria   Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown
4. ASV71 Bacteria   Firmicutes  Bacilli Lactobacillales Lactobacillaceae    Lactobacillus   Lactobacillus_unknown
5. ASV57 Bacteria   Firmicutes  Clostridia  Clostridiales   Clostridiaceae  Clostridium Clostridium_unknown

Plot for figure

```R
impToPlot.season <- importance(rf_season, scale=F)[,3]
impToPlot.season <- sort(impToPlot.season)
short.imp <- tail(impToPlot.season, 10)
pdf("rf_season_importance.filt.pdf")
dotchart(short.imp, xlim=c(0.00, 0.06))
dev.off()
impToPlot.season <- importance(rf_season, scale=F)[,4]
impToPlot.season <- sort(impToPlot.season)
short.imp <- tail(impToPlot.season, 10)
pdf("rf_season_importance.filt.gini.pdf")
dotchart(short.imp, xlim=c(0.0, 2.5))
dev.off()
impToPlot.insect <- importance(rf_insects, scale=F)[,3]
impToPlot.insect <- sort(impToPlot.insect)
short.imp <- tail(impToPlot.insect, 10)
pdf("rf_insects_importance.filt.pdf")
dotchart(short.imp, xlim=c(0.00, 0.06))
dev.off()
impToPlot.insect <- importance(rf_insects, scale=F)[,4]
impToPlot.insect <- sort(impToPlot.insect)
short.imp <- tail(impToPlot.insect, 10)
pdf("rf_insects_importance.filt.gini.pdf")
dotchart(short.imp, xlim=c(0.0, 2.5))
dev.off()
```
Abundance of important taxa by temperature range

```R
library(plyr)
library(ggplot2)
getwd()
# [1] "/Users/mann/github/necrobiome/02-analysis"
tempdat <- read.table("map_abund.csv", sep=",", header=T)
# function to summarize data
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
# summarize the data
sumdat <- data_summary(tempdat, varname="log10Abund", groupnames=c("Genus", "Temp_group"))
# convert temp group to factor
sumdat$Temp_group <- as.factor(sumdat$Temp_group)
# graph 
pd <- position_dodge(0.1)
pdf("imgs/Fig5_temp_gradient_abundance.pdf")
ggplot(sumdat, aes(x=Temp_group, y=log10Abund, group=Genus, color=Genus)) + geom_line(position=pd, size=2) + geom_errorbar(aes(ymin=log10Abund-sd,ymax=log10Abund+sd), color="darkgrey", width=.1, position=pd) + geom_point(position=pd, size=3, shape=21, fill="white") + theme_minimal()
dev.off()
png("imgs/Fig5_temp_gradient_abundance.png")
ggplot(sumdat, aes(x=Temp_group, y=log10Abund, group=Genus, color=Genus)) + geom_line(position=pd, size=2) + geom_errorbar(aes(ymin=log10Abund-sd,ymax=log10Abund+sd), color="darkgrey", width=.1, position=pd) + geom_point(position=pd, size=3, shape=21, fill="white") + theme_minimal()
dev.off()
```

![rf insects](https://github.com/aemann01/necrobiome/blob/master/02-analysis/imgs/Fig5_temp_gradient_abundance.png)

