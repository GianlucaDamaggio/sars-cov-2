# Sars-Cov-2

## Alignment with [minimap2](https://github.com/lh3/minimap2)
```
# default output
minimap2 -t 8 ~/projects/sars-cov-2/reference/NC_045512.2.fa ~/projects/sars-cov-2/data/fastq_pass/barcode01/all_barcode01.fastq -o ~/projects/sars-cov-2/alignment/barcode01.paf

for id in 03 04 05 06 07 08 09 10 11 12 ; do echo qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/minimapPAF_barcode$id.out -e /lustrehome/gianluca/junk/minimapPAF_barcode$id.err -v id="$id" -N miniPAF_barcode$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-minimap_paf.sh ; done

# SAM outout
minimap2 -t 8 -a ~/projects/sars-cov-2/reference/NC_045512.2.fa ~/projects/sars-cov-2/data/fastq_pass/barcode01/all_barcode01.fastq -o ~/projects/sars-cov-2/alignment/barcode01.sam

for id in 03 04 05 06 07 08 09 10 11 12 ; do echo qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/minimapSAM_barcode$id.out -e /lustrehome/gianluca/junk/minimapSAM_barcode$id.err -v id="$id" -N miniSAM_barcode$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-minimap_sam.sh ; done

```
Count how many reads are in input fastq and compare with minimap2 run
```
awk '(NR-2)%4==0 { print $1}' ~/projects/sars-cov-2/data/fastq_pass/barcode01/all_barcode01.fastq | wc -l
```
<!-- ### Valuate with Dotploty alignment with reference

Not work !!
```
~/src/dotPlotly/./pafCoordsDotPlotly.R -i /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode02/barcode02.paf -o barcode02_dotplotly -s -t -m 10 -q 10 -s -p 15

for id in 03 04 05 06 07 08 09 10 11 12 ; do echo qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/dotplotly_barcode$id.out -e /lustrehome/gianluca/junk/dotplotly_barcode$id.err -v id="$id" -N dotplotly_barcode$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-dotPlotly.sh ; done

``` -->

<!-- ### Length distribution of nanopore reads

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

``` -->

## Coverage

Distribution of coverage along the sequence

<!-- ```
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
``` -->
Coverage for each bases
```
bedtools genomecov -d -ibam ~/projects/sars-cov-2/alignment/barcode01_sorted.bam > ~/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.txt

for id in 03 04 05 06 07 08 09 10 11 12 ; do echo qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/genomecov_barcode$id.out -e /lustrehome/gianluca/junk/genomecov_barcode$id.err -v id="$id" -N genomecov_barcode$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-genomecov.sh ; done

echo reference position cov | tr " " "\t" > header.txt

for id in 03 04 05 06 07 08 09 10 11 12 ; do cat header.txt /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20201130_1527_1F_PAG19577_0a1c6f33/coverage/barcode$id/barcode$id.coverage_eachBases.txt > /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20201130_1527_1F_PAG19577_0a1c6f33/coverage/barcode$id/barcode$id.coverage_eachBases.tsv ; done
```
## Coverage analysis with different ID (seminario IGB 18/12/2020)

```
Rscript /Users/gianlucadamaggio/projects/sars-cov-2/doc/comparison_differentID.R
```

<!-- Plot coverage for each bases with R
```
library(tidyverse)
library (grid)
library(lattice)
library(gridExtra)

bar01=read.csv("~/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.txt", header=F , sep="\t")
bar02=read.csv("~/projects/sars-cov-2/coverage/barcode02_coverage_eachBases.txt", header=F , sep="\t")

plot01=ggplot(bar01, aes(V2,V3)) + geom_area( fill="#69b3a2", alpha=0.4) + geom_line(color="#69b3a2", size=0.4) + xlab("length") + ylab("coverage") + ggtitle("Barcode01 Per-base coverage distribution")

plot02=ggplot(bar02, aes(V2,V3)) + geom_area( fill="#69b3a2", alpha=0.4) + geom_line(color="#69b3a2", size=0.4) + xlab("length") + ylab("coverage") + ggtitle("Barcode02 Per-base coverage distribution")

lay <- rbind(c(1,1),
             c(2,2))

pgrid=grid.arrange(plot01, plot02, nrow = 2, layout_matrix = lay)

ggsave(plot=pgrid,filename="~/projects/sars-cov-2/coverage/coverage_perBase.png",height=10, width=20)

ggsave(plot=plot01,filename="~/projects/sars-cov-2/plots/barcode01_coverage_eachBases.png")
ggsave(plot=plot02,filename="~/projects/sars-cov-2/plots/barcode02_coverage_eachBases.png")

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
``` -->

## Assembly using [ArticNetwork](https://github.com/artic-network/artic-ncov2019)
Demultiplex
```
~/src/ont-guppy-cpu/bin/./guppy_barcoder --require_barcodes_both_ends -i /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/fastq_pass/ -s /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/demultiplex/ --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"
```
Read filtering

.conda/envs/artic-ncov2019/bin/artic

```
/lustrehome/gianluca/.conda/envs/artic-ncov2019/bin/artic guppyplex --min-length 400 --max-length 700 --directory /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/demultiplex/barcode02 --prefix 20200804_1408_MN25488_FAO00129_31063f42

for id in 03 04 05 06 07 08 09 10 11 12 ; do echo qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/guppyplex_barcode$id.out -e /lustrehome/gianluca/junk/guppyplex_barcode$id.err -v id="$id" -N guppyplex_barcode$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-guppyplex.sh ; done

#move output
```
Run the MinION pipeline
```
artic minion --normalise 200 --threads 4 --scheme-directory /lustrehome/gianluca/github/artic-ncov2019/primer_schemes --read-file /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/demultiplex/barcode01/20200804_1408_MN25488_FAO00129_31063f42_barcode01.fastq --fast5-directory /lustre/home/enza/sars-cov-2/ceinge/data/Covid-19_Run1/Pool_covid/20200804_1408_MN25488_FAO00129_31063f42/fast5_pass --sequencing-summary /lustre/home/enza/sars-cov-2/ceinge/data/Covid-19_Run1/Pool_covid/20200804_1408_MN25488_FAO00129_31063f42/sequencing_summary_FAO00129_be760d9d.txt nCoV-2019/V3 20200804_1408_MN25488_FAO00129_31063f42

for id in 03 04 05 06 07 08 09 10 11 12 ; do echo qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/artic_barcode$id.out -e /lustrehome/gianluca/junk/artic_barcode$id.err -v id="$id" -N artic_barcode$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-artic.sh ; done
```


<!-- ## Obtain regions with coverage < threshold
```

library(tidyverse)
library (grid)
library(lattice)
library(gridExtra)

bar01=read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.tsv", sep="\t" , header=T)
bar02=read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/barcode02_coverage_eachBases.tsv", sep="\t" , header=T)

bar01$name="bar01"
bar02$name="bar02"

bar=rbind(bar01,bar02)

# write table for obtain interval
min20_cov= bar %>% filter(cov <= 20)
min50_cov= bar %>% filter(cov <= 50)

write.table(min20_cov, "/Users/gianlucadamaggio/projects/sars-cov-2/coverage/barcode_Min20_gap.txt" , sep="\t", quote=F , row.names=F)
write.table(min50_cov, "/Users/gianlucadamaggio/projects/sars-cov-2/coverage/barcode_Min50_gap.txt" , sep="\t", quote=F , row.names=F)

# plot
cov20 = bar %>% filter ( cov < 20) %>% ggplot( aes(position, name, fill=name )) + geom_tile()+ labs(colour="samples",x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("cov < 20") + theme(text = element_text(size = 20))

cov50 = bar %>% filter ( cov < 50) %>% ggplot( aes(position, name, fill=name  )) + geom_tile()+ labs(colour="samples",x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("cov < 50") + theme(text = element_text(size = 20))

lay <- rbind(c(1,1),
             c(2,2))

pgrid=grid.arrange(cov20, cov50, nrow = 2, layout_matrix = lay)

ggsave(plot=pgrid,filename="~/projects/sars-cov-2/coverage/gap_cov20cov50.png",height=10, width=20)



``` -->

## Obtain table for GAP interval (position with coverage < threshold)

delete header and filter by barcode
```
cat ~/projects/sars-cov-2/coverage/barcode_Min50_gap.txt | cut -f 2,4 | grep -v "position" |grep "bar01" | cut -f 1 > bar01_gapInterval_Min50.txt
cat ~/projects/sars-cov-2/coverage/barcode_Min50_gap.txt | cut -f 2,4 | grep -v "position" |grep "bar02" | cut -f 1 > bar02_gapInterval_Min50.txt
cat ~/projects/sars-cov-2/coverage/barcode_Min20_gap.txt | cut -f 2,4 | grep -v "position" |grep "bar02" | cut -f 1 > bar02_gapInterval_Min20.txt
cat ~/projects/sars-cov-2/coverage/barcode_Min20_gap.txt | cut -f 2,4 | grep -v "position" |grep "bar01" | cut -f 1 > bar01_gapInterval_Min20.txt
```

```
python3.7 gap_interval.py
```
or manually
```
import pandas as pd
from operator import itemgetter
import itertools

num_list=[]
with open('bar01_gapInterval_Min20.txt', 'r') as fh:
    for line in fh:
        num_list.append(int(line))

first_num=[]
for k, g in itertools.groupby( enumerate(num_list), lambda x: x[1]-x[0] ) :
  first_num.append(list(map(itemgetter(1), g))[0])


last_num=[]
for k, g in itertools.groupby( enumerate(num_list), lambda x: x[1]-x[0] ) :
   last_num.append(list(map(itemgetter(1), g))[-1])

df = pd.DataFrame(list(zip(first_num, last_num)), columns =['start', 'end'])

df

df.to_csv('~/projects/sars-cov-2/coverage/bar01_gap.csv')
```

## Assembly GAP from Consensus fasta (denovo assemby) https://bioinf.shenwei.me/seqkit/

```
for n in 03 04 05 06 07 08 09 10 11 12; do seqkit locate --ignore-case --only-positive-strand -r --pattern "N+" /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20201130_1527_1F_PAG19577_0a1c6f33/artic/barcode$n/20201130_1527_1F_PAG19577_0a1c6f33.consensus.fasta > /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20201130_1527_1F_PAG19577_0a1c6f33/artic/barcode$n/assemblyGap_barcode$n.tsv ; done
```

<!-- ```
library(tidyverse)

bar01=read.table("~/projects/sars-cov-2/assembly/gap_interval/assemblyGap_barcode01.tsv", header=T, sep="\t")
bar02=read.table("~/projects/sars-cov-2/assembly/gap_interval/assemblyGap_barcode02.tsv", header=T, sep="\t")

bar01$ID = "bar01"
bar02$ID = "bar02"

bar=rbind(bar01,bar02)

gap_interval = bar %>% ggplot( aes(start, ID, fill=ID )) + geom_tile()+ labs(colour="ID",x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("assembly gap") + theme(text = element_text(size = 20))

ggsave(plot=gap_interval,filename="~/projects/sars-cov-2/assemby/assembyGap_interval.png",height=10, width=20)

``` -->



## Variant Effect Predictor [custom](http://www.ensembl.info/2020/08/28/cool-stuff-the-ensembl-vep-can-do-annotating-sars-cov-2-variants/) for Sars-Cov-2

### VEP on run1
```
singularity exec -B /lustre/home/enza/ /lustre/home/enza/biocontainers/vep-101.sif vep -i /lustre/home/enza/sars-cov-2/ceinge/assembly/artic_analysis/20200804_1408_MN25488_FAO00129_31063f42/analysis/artic/sars_merged.vcf.gz -gff /lustre/home/enza/sars-cov-2/ceinge/variantCalling/data/ensembl/gff3/sars_cov_2_vep100/Sars_cov_2.ASM985889v3.100.primary_assembly.MN908947.3.gff.gz -fasta /lustre/home/enza/sars-cov-2/ceinge/variantCalling/data/ensembl/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna_sm.toplevel.fa.gz --force_overwrite
```

### VEP on run2 Promethaion
```
for id in 03 04 05 06 07 08 09 10 11 12; do qsub -q testqueue -l nodes=1:ppn=40 -o /lustrehome/gianluca/junk/vep_$id.out -e /lustrehome/gianluca/junk/vep_$id.err -v id="$id" -N vep_$id /lustrehome/gianluca/jobs/sars-cov-2/pbs-sarsCov2_VEP.sh ; done
```

#### Count number of reads for each Barcode ()

(base) [gianluca@ui02 barcode03]$ for c in 03 04 05 06 07 08 09 10 11 12; do awk '(NR-2)%4==0 { print $1}' ../barcode$c/20201130_1527_1F_PAG19577_0a1c6f33_barcode$c.fastq | wc -l ; done
3641115
1679310
1431084
215080
3899796
2195511
5361839
2126633
3822106
1954960
