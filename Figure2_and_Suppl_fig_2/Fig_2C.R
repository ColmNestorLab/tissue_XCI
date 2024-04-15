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
suppl_table4 <- subset(fread("supplementary_tables/Supplementary_Table4.tsv"), individual != "UPIC")

# select known escape, inactive, variable and PAR genes from the Tukiainen landscape paper.
goiz_ZZPU <- c("AKAP17A", "DDX3X", "APOOL", "PRKX")
goiz_13PLJ <- c("AKAP17A", "PUDP", "APOOL", "PRKX")

########## fig 2C ###############
ggsave2(filename = paste0(plot_dir, "Fig_2C.pdf"), width = 4, height = 4,
        
     ggplot(data=suppl_table4[(suppl_table4$gene %in% goiz_ZZPU & suppl_table4$individual == "nmXCI-2") | suppl_table4$gene %in% goiz_13PLJ & suppl_table4$individual == "nmXCI-1",], aes(x=factor(gene, levels = c("AKAP17A", "DDX3X", "PUDP", "APOOL", "PRKX")), y=allelic_expression)) + 
  geom_boxplot() +
  geom_quasirandom(aes(col=tissue_id)) +
  facet_wrap(~individual, scales = "free_x") + 
  theme_AL_box_rotX() + 
  theme(axis.text.x = element_text(size=6), 
        axis.text.y = element_text(size=6), 
        legend.position = "top", 
        legend.text = element_text(size=6), 
        strip.text.x = element_text(size=6, margin=margin(t=0, b=0))) +
  geom_hline(yintercept = 0.4, lty=2, size = 0.5)+
  scale_y_continuous(limits = c(0,0.5), breaks = c(0,0.1,0.2,0.3,0.4,0.5))+
  labs(x="", y="")+
    scale_color_manual(values = color_valzz)

)


