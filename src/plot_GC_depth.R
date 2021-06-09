library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(gridExtra)
library(EnvStats)

mh_gc_depth_table <- fread('output/gc_depth/Mh/gc_vs_depth_table.csv')
mo_gc_depth_table <- fread('output/gc_depth/MO/gc_vs_depth_table.csv')
fr_gc_depth_table <- fread('output/gc_depth/FR/gc_vs_depth_table.csv')

##remove other
gc_depth_plot <- subset(gc_depth_table, !(plot_label=="Other contig"))

##mean depth boxplot
depth_plot <- ggplot(gc_depth_plot, aes(x=plot_label, y=meandepth, colour=plot_label))+
  geom_boxplot()+
  scale_colour_viridis(discrete=TRUE, direction = -1)+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())

##gc boxplot
gc_plot <- ggplot(gc_depth_plot, aes(x=plot_label, y=GC, colour=plot_label))+
  geom_boxplot()+
  stat_n_text(y.pos = 0.19, hjust=0)+
  scale_colour_viridis(discrete=TRUE, direction = -1)+
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())+
  ylim(0.19, 0.42)+
  coord_flip()

##scatterplot
scatterplot <- ggplot(gc_depth_plot %>% arrange(plot_label), aes(x=GC, y=meandepth, colour=plot_label)) +
  geom_point()+
  scale_colour_viridis(discrete=TRUE, direction = -1)+
  theme_bw()+
  xlab("GC content")+
  ylab("Mean sequencing depth")+
  labs(colour="Contig identity")+
  xlim(0.19, 0.42)+
  theme(legend.position = "none")

##blank
blankplot <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

##arrange in grid
all_plots <- grid.arrange(gc_plot, blankplot, scatterplot, depth_plot, 
             ncol=2, nrow=2, widths=c(2.5, 0.5), heights=c(0.5, 2.5))

##save as svg and fix in inkscape
ggsave(file="output/gc_depth/svg/fr_gc_depth.svg", plot=all_plots, width=5, height=5)
