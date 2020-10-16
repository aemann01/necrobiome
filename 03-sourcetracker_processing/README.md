### Want to see if we have significant contamination in our samples

First need to convert transposed table to biom format for sourcetracker 2

```bash
biom convert -i ../01-raw_data_processing/sequence_table.16s.filtered.tr.txt --table-type="OTU table" --to-hdf5 -o asv_table.biom
```
Activate sourcetracker conda environment

```bash
conda activate st2
```

Run sourcetracker

```bash
rm -r results/
sourcetracker2 gibbs -i asv_table.biom -m ../02-analysis/map.txt -o results
```