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
busco_results_file <- snakemake@input[['busco']]

########
# MAIN #
########

##data
gc_content_table <- fread(gc)
gc_hist <- fread(gc_hist_file)
##scaffold label files
viral_contigs <- fread(viral_contig_list, header=FALSE)
##busco labels
busco_res <- fread(busco_results_file, header=TRUE, fill=TRUE, skip = 4)
busco_contigs <- busco$Contig

##GC histogram
pdf(snakemake@output[["gc_histogram"]])
ggplot(gc_hist, aes(x=gc_hist$`#GC`, y=scaffolds))+
  geom_col(alpha=0.7, fill="#440154FF", colour="#440154FF")+
  xlim(0.12, 0.58)+
  xlab("GC content")+
  ylab("Number of contigs")+
  theme_classic()
dev.off()

gc_content_table$plot_label <- ifelse(gc_content_table$scaffold_id %in% viral_contigs$V1, 'Viral contig',
                                             ifelse(gc_content_table$scaffold_id %in% busco_contigs, 'BUSCO contig', 'Other contig'))

##plot GC content
pdf(snakemake@output[["all_contigs"]])
ggplot(gc_content_table, aes(x=plot_label, y=GC))+
  geom_boxplot()+
  theme_classic()+
  xlab("Scaffold identity")+
  ylab("GC content")+
  ylim(0,1)
dev.off()

##make plot without others for hic assemblies
busco_plot_table <- subset(gc_content_table, !(plot_label=="Other contig"))
pdf(snakemake@output[["busco_only"]])
ggplot(busco_plot_table, aes(x=plot_label, y=GC))+
    geom_boxplot()+
  theme_classic()+
  xlab("Scaffold identity")+
  ylab("GC content")
dev.off()

#write log
sessionInfo()