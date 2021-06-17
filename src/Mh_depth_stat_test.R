library(data.table)
library(dplyr)

gc_depth <- fread('output/gc_depth/Mh/gc_vs_depth_table.csv')
high_outliers <- list("scaffold_90", "scaffold_995")
gc_depth_fil <- subset(gc_depth, !(`#Name` %in% high_outliers))

##test normality of meandepth
shapiro <- shapiro.test(gc_depth_fil$meandepth)
shapiro_chars <- capture.output(print(shapiro))

##cannot use parametric tests

##########################
## non-parametric tests ##
##########################

summary <- group_by(gc_depth_fil, plot_label) %>%
  summarise(
    count = n(),
    mean = mean(meandepth, na.rm = TRUE),
    sd = sd(meandepth, na.rm = TRUE),
    median = median(meandepth, na.rm = TRUE),
    IQR = IQR(meandepth, na.rm = TRUE)
  )
summary_chars <- capture.output(print(summary))

##does meandepth content between any groups differ significantly
kruskal_res <- kruskal.test(meandepth ~ plot_label, data = gc_depth_fil)
kruskal_chars <- capture.output(print(kruskal_res))

##which groups differ significantly
pair_wilcox_res <- pairwise.wilcox.test(gc_depth_fil$meandepth, gc_depth_fil$plot_label,
                                        p.adjust.method = "BH")
wilcox_chars <- capture.output(print(pair_wilcox_res))
