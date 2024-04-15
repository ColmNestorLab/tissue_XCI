library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(ggplotify)
library(karyoploteR)
library(UpSetR)
library(grid)
library(gridExtra)
library(VennDiagram)
library(raster)
library(rstatix)
library(ggpubr)
library(ggbeeswarm)
source("r_sources/plot_parameters.R")

##### Read in data ##### 

# Supplementary table 4
df_tissues_plot <- unique(fread("supplementary_tables/Supplementary_Table4.tsv")[,c("individual", "tissue_id")])

df_tissues_plot$valz <- 1

tissue_order <- unique(sort(df_tissues_plot$tissue_id))

df_tissues_plot_complete <- tidyr::complete(df_tissues_plot, tissue_id, individual)
df_tissues_plot_complete$tissue_availability <- ifelse(is.na(df_tissues_plot_complete$valz), no = "available", yes = "missing")

tissues_with_overlap <- 
  unique(fread("supplementary_tables/Supplementary_Table4.tsv")[,c("individual", "tissue_id")]) %>% 
  dplyr::group_by(tissue_id) %>% 
  dplyr::count() %>% 
  subset(n > 1)

levlz <- c(tissues_with_overlap$tissue_id, setdiff(df_tissues_plot_complete$tissue_id, tissues_with_overlap$tissue_id))

ggsave2(filename = "plots/figure_3B.pdf", width = 6.5,
        ggplot(df_tissues_plot_complete, aes(x=factor(individual, levels = c("nmXCI-1", "nmXCI-2", "UPIC")), y=factor(tissue_id, levels =  rev(levlz)), fill = tissue_availability)) + geom_tile(col="black") + coord_fixed(1/1) + scale_fill_manual(values = c( "missing" = "grey", "available" = "red")) + labs(x="individual", y="tissue")  + guides(fill=guide_legend(ncol =1)) + theme_AL_box_rotX(legend.position = "top")
)
