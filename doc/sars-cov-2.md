# Sars-Cov-2

## Alignment with [minimap2](https://github.com/lh3/minimap2)
```
# default output
minimap2 -t 8 ~/projects/sars-cov-2/reference/NC_045512.2.fa ~/projects/sars-cov-2/data/fastq_pass/barcode01/all_barcode01.fastq -o ~/projects/sars-cov-2/alignment/barcode01.paf

# SAM outout
minimap2 -t 8 -a ~/projects/sars-cov-2/reference/NC_045512.2.fa ~/projects/sars-cov-2/data/fastq_pass/barcode01/all_barcode01.fastq -o ~/projects/sars-cov-2/alignment/barcode01.sam
```
Count how many reads are in input fastq and compare with minimap2 run
```
awk '(NR-2)%4==0 { print $1}' ~/projects/sars-cov-2/data/fastq_pass/barcode01/all_barcode01.fastq | wc -l
```

### Length distribution of nanopore reads

Count length of each read from fastq file

```
awk '(NR-2)%4==0 { print length($1)}' all_barcode01.fastq > length.txt
```

Plot lenght distribution with R

```
library(tidyverse)

myd = read.csv("~/projects/sars-cov-2/data/fastq_pass/barcode01/length.txt", header=F)

ggplot(myd, aes(V1)) + geom_density()
ggsave("~/projects/sars-cov-2/plots/barcode01_sars_length.pdf")

ggplot(myd, aes(V1)) + geom_density() + xlim (0,1000)
ggsave("~/projects/sars-cov-2/plots/barcode01_sars_length_Max1k.pdf")

```

## Coverage

Distribution of coverage along the sequence

```
samtools sort ~/projects/sars-cov-2/alignment/barcode01.sam > ~/projects/sars-cov-2/alignment/barcode01_sorted.sam

# tabular output
samtools coverage ~/projects/sars-cov-2/alignment/barcode01_sorted.sam > ~/projects/sars-cov-2/coverage/barcode01_coverage.txt

# histogram output
samtools coverage -m ~/projects/sars-cov-2/alignment/barcode01_sorted.sam > ~/projects/sars-cov-2/coverage/barcode01_coverageHisto.txt
```
Convert SAM to BAM
```
samtools view -b ~/projects/sars-cov-2/alignment/barcode01_sorted.sam > ~/projects/sars-cov-2/alignment/barcode01_sorted.bam
```
Print average of total coverage
```
samtools depth  ~/projects/sars-cov-2/alignment/barcode01_sorted.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}'
```
Coverage for each bases
```
bedtools genomecov -d -ibam ~/projects/sars-cov-2/alignment/barcode01_sorted.bam > ~/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.txt
```
Plot coverage for each bases with R
```
library(tidyverse)

myd=read.csv("~/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.txt", header=F , sep="\t")

ggplot(myd, aes(V2,V3)) + geom_area( fill="#69b3a2", alpha=0.4) + geom_line(color="#69b3a2", size=0.4) + xlab("length") + ylab("coverage") + ggtitle("Barcode01")

ggsave("~/projects/sars-cov-2/plots/barcode01_coverage_eachBases.pdf")
```
Plot Barcode01/02 together with R
```
library(tidyverse)

myd01=read.csv("~/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.txt", header=F , sep="\t")
myd02=read.csv("~/projects/sars-cov-2/coverage/barcode02_coverage_eachBases.txt", header=F , sep="\t")

myd01$id = "Barcode01"
myd02$id = "Barcode02"

myd_final = rbind(myd01,myd02)

ggplot(myd_final, aes(V2,V3, color=id)) + geom_line() + xlab("length") + ylab("coverage") + ggtitle("Coverage Distribution")

ggsave("~/projects/sars-cov-2/plots/all_coverage_eachBases.pdf")
```

## Assembly using [ArticNetwork](https://github.com/artic-network/artic-ncov2019)
Demultiplex
```
~/src/ont-guppy-cpu/bin/./guppy_barcoder --require_barcodes_both_ends -i /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/fastq_pass/ -s /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/demultiplex/ --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
```
Read filtering
```
artic guppyplex --min-length 400 --max-length 700 --directory /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/demultiplex/barcode02 --prefix 20200804_1408_MN25488_FAO00129_31063f42

#move output
```
Run the MinION pipeline
```
artic minion --normalise 200 --threads 4 --scheme-directory /lustrehome/gianluca/github/artic-ncov2019/primer_schemes --read-file /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/demultiplex/barcode01/20200804_1408_MN25488_FAO00129_31063f42_barcode01.fastq --fast5-directory /lustre/home/enza/sars-cov-2/ceinge/data/Covid-19_Run1/Pool_covid/20200804_1408_MN25488_FAO00129_31063f42/fast5_pass --sequencing-summary /lustre/home/enza/sars-cov-2/ceinge/data/Covid-19_Run1/Pool_covid/20200804_1408_MN25488_FAO00129_31063f42/sequencing_summary_FAO00129_be760d9d.txt nCoV-2019/V3 20200804_1408_MN25488_FAO00129_31063f42
```
