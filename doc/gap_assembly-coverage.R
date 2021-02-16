
library(tidyverse)
library (grid)
library(lattice)
library(gridExtra)

bar01=read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/barcode01_coverage_eachBases.tsv", sep="\t" , header=T)
bar02=read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/barcode02_coverage_eachBases.tsv", sep="\t" , header=T)

bar01$ID="bar01"
bar02$ID="bar02"

bar=rbind(bar01,bar02)

# plot
bar01cov20 = bar %>% filter ( cov < 20, ID =="bar01") %>% ggplot( aes(position, ID )) + geom_tile()+ labs(x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("cov < 20") + theme(text = element_text(size = 15))+ theme_bw()
bar02cov20 = bar %>% filter ( cov < 20, ID =="bar02") %>% ggplot( aes(position, ID )) + geom_tile()+ labs(x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("cov < 20") + theme(text = element_text(size = 15))+ theme_bw()


bar01cov50 = bar %>% filter ( cov < 50, ID =="bar01") %>% ggplot( aes(position, ID )) + geom_tile()+ labs(x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("cov < 50") + theme(text = element_text(size = 15))+ theme_bw()
bar02cov50 = bar %>% filter ( cov < 50, ID =="bar02") %>% ggplot( aes(position, ID )) + geom_tile()+ labs(x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("cov < 50") + theme(text = element_text(size = 15))+ theme_bw()

gapbar01=read.table("~/projects/sars-cov-2/assembly/gap_interval/assemblyGap_barcode01.tsv", header=T, sep="\t")
gapbar02=read.table("~/projects/sars-cov-2/assembly/gap_interval/assemblyGap_barcode02.tsv", header=T, sep="\t")

gapbar01$ID="bar01"
gapbar02$ID="bar02"

gapbar=rbind(gapbar01,gapbar02)

bar01gap_interval = gapbar %>% filter(ID == "bar01" ) %>% ggplot( aes(start, ID )) + geom_tile()+ labs(x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("assembly gap") + theme(text = element_text(size = 15))+ theme_bw()
bar02gap_interval = gapbar %>% filter(ID == "bar02" ) %>% ggplot( aes(start, ID )) + geom_tile()+ labs(x="Genomic position (bp)") + theme(legend.position="top")+ ggtitle("assembly gap") + theme(text = element_text(size = 15))+ theme_bw()

lay <- rbind(c(1,1),
             c(2,2),
             c(3,3))

bar01pgrid=grid.arrange(bar01cov20, bar01cov50, bar01gap_interval, nrow = 3, layout_matrix = lay)
bar02pgrid=grid.arrange(bar02cov20, bar02cov50, bar02gap_interval, nrow = 3, layout_matrix = lay)

ggsave(plot=bar01pgrid,filename="~/projects/sars-cov-2/assembly/bar01assembyGap&cov_interval.png",height=8, width=10)
ggsave(plot=bar02pgrid,filename="~/projects/sars-cov-2/assembly/bar02assembyGap&cov_interval.png",height=8, width=10)
