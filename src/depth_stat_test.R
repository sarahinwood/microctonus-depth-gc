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

###########
# GLOBALS #
###########

gc_depth_table <- snakemake@input[["gc_depth_table"]]

####################
## test normality ##
####################

gc_depth <- fread(gc_depth_table, na.strings="")
gc_depth_dna <- subset(gc_depth, (viral_genome=="DNA" | is.na(viral_genome)))

##test normality of meandepth
shapiro <- shapiro.test(gc_depth_dna$meandepth)
shapiro_chars <- capture.output(print(shapiro))
writeLines(shapiro_chars, con=file(snakemake@output[["shapiro_res"]]))

##cannot use parametric tests

##########################
## non-parametric tests ##
##########################

summary <- group_by(gc_depth_dna, plot_label) %>%
  summarise(
    count = n(),
    mean = mean(meandepth, na.rm = TRUE),
    sd = sd(meandepth, na.rm = TRUE),
    median = median(meandepth, na.rm = TRUE),
    IQR = IQR(meandepth, na.rm = TRUE)
  )
summary_chars <- capture.output(print(summary))
writeLines(summary_chars, con=file(snakemake@output[["summary_stats"]]))

##does meandepth content between any groups differ significantly
kruskal_res <- kruskal.test(meandepth ~ plot_label, data = gc_depth_dna)
kruskal_chars <- capture.output(print(kruskal_res))
writeLines(kruskal_chars, con=file(snakemake@output[["kruskal_res"]]))

##which groups differ significantly
pair_wilcox_res <- pairwise.wilcox.test(gc_depth_dna$meandepth, gc_depth_dna$plot_label,
                                        p.adjust.method = "BH")
wilcox_chars <- capture.output(print(pair_wilcox_res))
writeLines(wilcox_chars, con=file(snakemake@output[["wilcox_res"]]))

#write log
sessionInfo()