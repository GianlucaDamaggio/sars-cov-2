
for n in 01 02 ; do minimap2 -t 8 -a /lustre/home/enza/sars-cov-2/ceinge/data/reference/NC_045512.2.fa /lustre/home/enza/sars-cov-2/ceinge/data/Covid-19_Run1/Pool_covid/20200804_1408_MN25488_FAO00129_31063f42/fastq_pass/barcode$n/all_barcode$n.fastq -o /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.sam ; done

for n in 01 02 ; do samtools view -S -b /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.sam > /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.bam ; done

for n in 01 02 ; do samtools sort  /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.bam > /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.sort.bam ; done

# Mpileup produce a VCF from a fasta and bam
for n in 01 02 ; do samtools mpileup -ABuf /lustre/home/enza/sars-cov-2/ceinge/data/reference/NC_045512.2.fa /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.sort.bam | bcftools call -cOz --pval-threshold 0.99 > /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/mapped.vcf.gz ; done

for n in 01 02 ; do tabix /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/mapped.vcf.gz ; done

for n in 01 02 ; do cat /lustre/home/enza/sars-cov-2/ceinge/data/reference/NC_045512.2.fa | bcftools consensus /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/mapped.vcf.gz > /lustre/home/enza/sars-cov-2/ceinge/alignment/barcode$n/barcode$n.mapped.fastq ; done

#Fasta from gff3
gff3_to_fasta -g 2019-nCoV_EPI_ISL_514432_variants.gff3 -f NC_045512.2.fa -st all -d simple -o 2019-nCoV_EPI_ISL_514432.fa



##### Multiple Alignment with Clustal-OMEGA http://www.clustal.org/omega/
/lustrehome/gianluca/src/clustalo --infile ~/All_Sars-Cov-2.fa --threads 8 --MAC-RAM 2000 --verbose --guidetree-out clustalo-I20201007-164636-0273-35651654-p2m.dnd --outfmt clustal --resno --outfile clustalo-I20201007-164636-0273-35651654-p2m.clustal_num --output-order tree-order --seqtype rna




```
library(ggtree)
sars=read.tree("clustalo-I20201007-164636-0273-35651654-p2m.dnd")
ggtree(sars) + geom_tiplab()
```

python -m venv /opt/conda/envs/venv_medaka --prompt "(fast5mod) " && source /opt/conda/envs/venv_medaka/bin/activate && pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir fast5mod
