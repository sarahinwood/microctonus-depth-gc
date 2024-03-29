#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

gc <- snakemake@input[['gc']]
gc_hist_file <- snakemake@input[['gc_hist']]
viral_contig_list <- snakemake@input[['viral_contig_list']]
hic_scaffold_list <- snakemake@input[['hic_scaffold_list']]

########
# MAIN #
########

##data
gc_content_table <- fread(gc)
gc_hist <- fread(gc_hist_file)
##scaffold label files
viral_contigs <- fread(viral_contig_list)
hic_scaffolds <- fread(hic_scaffold_list, header=FALSE)
viral_hic <- intersect(viral_contigs$contig_id, hic_scaffolds$V1)

##GC histogram
pdf(snakemake@output[["gc_histogram"]])
ggplot(gc_hist, aes(x=gc_hist$`#GC`, y=scaffolds))+
  geom_col(alpha=0.7, fill="#440154FF", colour="#440154FF")+
  xlim(0.12, 0.58)+
  xlab("GC content")+
  ylab("Number of contigs")+
  theme_classic()
dev.off()

gc_content_table$plot_label <- ifelse(gc_content_table$`#Name` %in% viral_hic, 'Hi-C scaffold and viral',
                                   ifelse(gc_content_table$`#Name` %in% viral_contigs$contig_id, 'Viral contig',
                                          ifelse(gc_content_table$`#Name` %in% hic_scaffolds$V1, 'Hi-C scaffold', 'Other contig')))

gc_content_table$plot_label <- factor(gc_content_table$plot_label, levels=c("Hi-C scaffold", "Hi-C scaffold and viral",
                                                                            "Viral contig", "Other contig"))

##make plot without others for hic assemblies
hic_plot_table <- subset(gc_content_table, !(plot_label=="Other contig"))
pdf(snakemake@output[["hic_only"]])
ggplot(hic_plot_table, aes(x=plot_label, y=GC, colour=plot_label))+
  geom_boxplot()+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw()+
  xlab("Scaffold identity")+
  ylab("% GC content")+
  theme(legend.position = "none")
dev.off()

fwrite(gc_content_table, snakemake@output[["gc_table"]])


#write log
sessionInfo()