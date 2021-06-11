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
library(ggplot2)
library(viridis)

###########
# GLOBALS #
###########

samtools_depth <- snakemake@input[["samtools_depth"]]
##scaffold labels
scaffold_id_table <- snakemake@input[["gc_depth_table"]]

########
# MAIN #
########

##scaffold ID table
scaffold_table <- fread(scaffold_id_table, header=TRUE)
##depth table
st_depth_names <- c("#Name", "BP", "depth")
st_depth <- fread(samtools_depth, col.names=st_depth_names)
##merge with table of scaffold ids for plotting
st_depth_labels <- merge(st_depth, scaffold_table, by='#Name', all.y=TRUE)
st_depth_boxpl <- subset(st_depth_labels, !(plot_label=="Other contig"))

#####################
## making box plot ##
#####################

pdf(snakemake@output[["boxplot_y_zoom"]])
ggplot(st_depth_boxpl, aes(x=reorder(`#Name`, plot_label), y=depth, colour=plot_label))+
  geom_boxplot(outlier.shape=NA)+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE)+
  coord_cartesian(ylim = c(0,1000))
dev.off()

jpeg(snakemake@output[["boxplot"]])
ggplot(st_depth_boxpl, aes(x=reorder(`#Name`, plot_label), y=depth, colour=plot_label))+
  geom_boxplot()+
  theme_bw(base_size=18)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  xlab("")+
  ylab("Depth")+
  stat_summary(fun=mean, geom="point", colour="grey35")+
  scale_colour_viridis(discrete=TRUE)
dev.off()

#write log
sessionInfo()
