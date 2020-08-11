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

## Assembly using correctify script from [viral-assembly](https://github.com/ekg/viral-assembly) repo

```
~/github/viral-assembly/scripts/correctify -f /lustre/home/enza/sars-cov-2/ceinge/data/Covid-19_Run1/Pool_covid/20200804_1408_MN25488_FAO00129_31063f42/fastq_pass/barcode01/all_barcode01.fastq -o /lustre/home/enza/sars-cov-2/ceinge/assembly/barcode01_assembly -k 15 -a 100 -G 350 -L 550 -p barcode01 -t 48
```
## Validate assembly using [MashMap](https://github.com/marbl/MashMap)
```
# move into WorkingDirectory
cd /lustre/home/enza/sars-cov-2/ceinge/assembly/mashmap/barcode01/

lustrehome/gianluca/github/MashMap/./mashmap -s 500 -r /lustre/home/enza/sars-cov-2/ceinge/data/reference/NC_045512.2.fa -q /lustre/home/enza/sars-cov-2/ceinge/assembly/barcode01_assembly/minia.k15.a100.contigs.fa
```
DotPlot
```
~/github/MashMap/scripts/./generateDotPlot png large mashmap.out
```
