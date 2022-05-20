## Statistical analyses

All R scripts can be found in the analysis jupyter notebook in this directory. To run:

```bash
jupyter-notebook analyses.ipynb
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

