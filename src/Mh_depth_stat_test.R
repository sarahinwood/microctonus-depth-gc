library(data.table)
library(dplyr)

gc_depth <- fread('output/gc_depth/Mh/gc_vs_depth_table.csv', na.strings="")
##remove high outliers to show still sig without them
high_outliers <- list("scaffold_90", "scaffold_995")
gc_depth$plot_label <- tstrsplit(gc_depth$plot_label, " and viral", keep=c(1))
gc_depth$plot_label <- ifelse(gc_depth$`#Name` %in% high_outliers, 'viral high-depth outlier', gc_depth$plot_label)

##test normality of meandepth
shapiro.test(gc_depth$meandepth)
shapiro_chars <- capture.output(print(shapiro))

##cannot use parametric tests

##########################
## non-parametric tests ##
##########################

summary <- group_by(gc_depth, plot_label) %>%
  summarise(
    count = n(),
    mean = mean(meandepth, na.rm = TRUE),
    sd = sd(meandepth, na.rm = TRUE),
    median = median(meandepth, na.rm = TRUE),
    IQR = IQR(meandepth, na.rm = TRUE)
  )
summary_chars <- capture.output(print(summary))

##does meandepth content between any groups differ significantly
kruskal.test(meandepth ~ plot_label, data = gc_depth)
kruskal_chars <- capture.output(print(kruskal_res))

##which groups differ significantly
pairwise.wilcox.test(gc_depth$meandepth, gc_depth$plot_label,
                                        p.adjust.method = "BH")
wilcox_chars <- capture.output(print(pair_wilcox_res))
