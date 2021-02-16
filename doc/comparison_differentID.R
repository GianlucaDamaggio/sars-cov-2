library(tidyverse)
library (grid)
library(lattice)
library(gridExtra)

SC_01r2 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode03.coverage_eachBases.tsv", header=T , sep="\t")
SC_02r2 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode04.coverage_eachBases.tsv", header=T , sep="\t")

SC_01r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run1/barcode01.coverage_eachBases.tsv", header=T , sep="\t")
SC_02r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run1/barcode02.coverage_eachBases.tsv", header=T , sep="\t")

SC_01r2$reference = as.factor("NC_045512.2")
SC_02r2$reference = as.factor("NC_045512.2")

SC_01_final=merge(SC_01r1, SC_01r2, by=c("reference","position"))
SC_02_final=merge(SC_02r1, SC_02r2, by=c("reference","position"))

SC_01_final <- SC_01_final[order(as.integer(SC_01_final$position),decreasing = FALSE),]
SC_02_final <- SC_02_final[order(as.integer(SC_02_final$position),decreasing = FALSE),]

lay <- rbind(c(1,1),
             c(2,2),
					   c(3,3))

coeff <- 10

plot01=ggplot(SC_01_final, aes(x=position)) + geom_area( aes(y=cov.x),fill="#69b3a2", alpha=0.7) + geom_area( aes(y=cov.y / coeff),fill="#fa9579", alpha=0.7) + scale_y_continuous(name = "Run1", sec.axis = sec_axis(~.*coeff, name="Run2 (x10)")) + theme(axis.title.y = element_text(color = "#69b3a2", size=13),axis.title.y.right = element_text(color = "#fa9579", size=13)) + ggtitle("SC_01 Per-base coverage distribution")

#ggsave(filename="~/projects/sars-cov-2/coverage/plots_compare/SC_01-comparison.coverage_eachBases.png")


plot02=ggplot(SC_02_final, aes(x=position)) + geom_area( aes(y=cov.x),fill="#69b3a2", alpha=0.7) + geom_area( aes(y=cov.y / coeff),fill="#fa9579", alpha=0.7) + scale_y_continuous(name = "Run1", sec.axis = sec_axis(~.*coeff, name="Run2 (x10)")) + theme(axis.title.y = element_text(color = "#69b3a2", size=13),axis.title.y.right = element_text(color = "#fa9579", size=13)) + ggtitle("SC_02 Per-base coverage distribution")

#ggsave(filename="~/projects/sars-cov-2/coverage/plots_compare/SC_02-comparison.coverage_eachBases.png")

# pgrid=grid.arrange(plot01, plot02, nrow = 2, layout_matrix = lay)
#
# ggsave(plot=pgrid,filename="~/projects/sars-cov-2/coverage/plots_compare/coverage_perBase_panel.png",height=10, width=20)

SC_03r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode05.coverage_eachBases.tsv", header=T , sep="\t")
SC_04r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode06.coverage_eachBases.tsv", header=T , sep="\t")
SC_06r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode07.coverage_eachBases.tsv", header=T , sep="\t")
SC_07r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode08.coverage_eachBases.tsv", header=T , sep="\t")
SC_08r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode09.coverage_eachBases.tsv", header=T , sep="\t")
SC_09r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode10.coverage_eachBases.tsv", header=T , sep="\t")
SC_10r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode11.coverage_eachBases.tsv", header=T , sep="\t")
SC_11r1 = read.table("/Users/gianlucadamaggio/projects/sars-cov-2/coverage/run2/barcode12.coverage_eachBases.tsv", header=T , sep="\t")

SC_03r1$ID = "SC_03"
SC_04r1$ID = "SC_04"
SC_06r1$ID = "SC_06"
SC_07r1$ID = "SC_07"
SC_08r1$ID = "SC_08"
SC_09r1$ID = "SC_09"
SC_10r1$ID = "SC_10"
SC_11r1$ID = "SC_11"

final_df=rbind(SC_03r1,SC_04r1,SC_06r1,SC_07r1,SC_08r1,SC_09r1,SC_10r1,SC_11r1)

plot03=ggplot(final_df, aes(position,cov, fill=ID)) + geom_area(alpha=0.6) + xlab("position") + ylab("coverage") + ggtitle("Multiples ID Per-base coverage distribution") 

# ggsave(filename="~/projects/sars-cov-2/coverage/plots_compare/all_barcode_coverage_perBase.png",height=10, width=20)


pgrid=grid.arrange(plot01, plot02, plot03, nrow = 3, layout_matrix = lay)

ggsave(plot=pgrid,filename="~/projects/sars-cov-2/coverage/plots_compare/coverage_perBase_panel.png",height=10, width=20,dpi=300)
