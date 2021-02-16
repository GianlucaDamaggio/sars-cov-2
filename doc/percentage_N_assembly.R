library(tidyverse)

SC_01r2 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode03.tsv", header=T , sep="\t")
SC_02r2 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode04.tsv", header=T , sep="\t")
SC_03r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode05.tsv", header=T , sep="\t")
SC_04r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode06.tsv", header=T , sep="\t")
SC_06r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode07.tsv", header=T , sep="\t")
SC_07r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode08.tsv", header=T , sep="\t")
SC_08r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode09.tsv", header=T , sep="\t")
SC_09r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode10.tsv", header=T , sep="\t")
SC_10r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode11.tsv", header=T , sep="\t")
SC_11r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/run2/assemblyGap_barcode12.tsv", header=T , sep="\t")

SC_01r2$ID = "SC_01"
SC_02r2$ID = "SC_02"
SC_03r1$ID = "SC_03"
SC_04r1$ID = "SC_04"
SC_06r1$ID = "SC_06"
SC_07r1$ID = "SC_07"
SC_08r1$ID = "SC_08"
SC_09r1$ID = "SC_09"
SC_10r1$ID = "SC_10"
SC_11r1$ID = "SC_11"

SC_01r2$dCT = 18.01
SC_02r2$dCT = 24.5
SC_03r1$dCT = 31.8
SC_04r1$dCT = 35
SC_06r1$dCT = 15.87
SC_07r1$dCT = 30
SC_08r1$dCT = 20.06
SC_09r1$dCT = 24.56
SC_10r1$dCT = 27
SC_11r1$dCT = 30

SC_01r2$end_start = SC_01r2$end - SC_01r2$start + 1
SC_01r2$NinAssembly = sum(SC_01r2$end_start) / 29900

SC_02r2$end_start = SC_02r2$end - SC_02r2$start + 1
SC_02r2$NinAssembly = sum(SC_02r2$end_start) / 29900

SC_03r1$end_start = SC_03r1$end - SC_03r1$start + 1
SC_03r1$NinAssembly = sum(SC_03r1$end_start) / 29900

SC_04r1$end_start = SC_04r1$end - SC_04r1$start + 1
SC_04r1$NinAssembly = sum(SC_04r1$end_start) / 29900

SC_06r1$end_start = SC_06r1$end - SC_06r1$start + 1
SC_06r1$NinAssembly = sum(SC_06r1$end_start) / 29900

SC_07r1$end_start = SC_07r1$end - SC_07r1$start + 1
SC_07r1$NinAssembly = sum(SC_07r1$end_start) / 29900

SC_08r1$end_start = SC_08r1$end - SC_08r1$start + 1
SC_08r1$NinAssembly = sum(SC_08r1$end_start) / 29900

SC_09r1$end_start = SC_09r1$end - SC_09r1$start + 1
SC_09r1$NinAssembly = sum(SC_09r1$end_start) / 29900

SC_10r1$end_start = SC_10r1$end - SC_10r1$start + 1
SC_10r1$NinAssembly = sum(SC_10r1$end_start) / 29900

SC_11r1$end_start = SC_11r1$end - SC_11r1$start + 1
SC_11r1$NinAssembly = sum(SC_11r1$end_start) / 29900

final_df=rbind(SC_01r2,SC_02r2,SC_03r1,SC_04r1,SC_06r1,SC_07r1,SC_08r1,SC_09r1,SC_10r1,SC_11r1)

final_df %>% group_by(ID, NinAssembly,dCT) %>% tally() %>% ggplot(aes(x=ID,y=NinAssembly)) + geom_bar(stat="identity",aes(fill = dCT)) + ylab("% of N in Assembly") + coord_flip() + theme_minimal() #+ theme(legend.position="none")

ggsave(filename="/Users/gianlucadamaggio/projects/sars-cov-2/assembly/gap_interval/plots/percent_N_assembly.png",height=4, width=6,dpi=300)
