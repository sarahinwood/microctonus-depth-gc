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

gc_table_file <- snakemake@input[["gc_table"]]
coverage_file <- snakemake@input[["coverage_file"]]
viral_genome_table <- snakemake@input[["viral_genome_table"]]

########
# MAIN #
########

gc_table <- fread(gc_table_file)
coverage <- fread(coverage_file)
viral_genome_type <- fread(viral_genome_table)

gc_table_virus <- merge(gc_table, viral_genome_type, by.x="#Name", by.y="contig_id", all.x=TRUE)

depth_table <- coverage[,c(1,7)]
gc_depth <- merge(gc_table_virus, depth_table, by.x="#Name", by.y="#rname")
gc_depth_plot <- gc_depth[order(-plot_label)]

##remove other
gc_depth_plot <- subset(gc_depth_plot, !(plot_label=="Other contig"))
##plot mean depth
pdf(snakemake@output[["depth_boxplot"]])
ggplot(gc_depth_plot, aes(x=plot_label, y=meandepth, colour=plot_label))+
  geom_boxplot()+
  scale_colour_viridis(discrete=TRUE)+
  theme_bw()+
  xlab("Scaffold identity")+
  ylab("Mean sequencing depth")+
  theme(legend.position = "none")
dev.off()

##plot gc vs depth
pdf(snakemake@output[["gc_vs_depth"]])
ggplot(gc_depth_plot %>% arrange(plot_label), aes(x=GC, y=meandepth, colour=plot_label)) +
  geom_point()+
  scale_colour_viridis(discrete=TRUE, direction = -1)+
  theme_bw()+
  xlab("GC content")+
  ylab("Mean sequencing depth")+
  labs(colour="Contig identity")
dev.off()

fwrite(gc_depth_plot, snakemake@output[['gc_depth_table']])
  
#write log
sessionInfo()