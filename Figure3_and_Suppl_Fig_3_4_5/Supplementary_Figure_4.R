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
source("r_sources/plot_parameters.R")

##### Read in data ##### 

# Supplementary table 4
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")

####### Suppl figure 1c #######
# get overlapping genes in overlapping tissues
# find overlapping genes
overlap1 <- intersect(unique(suppl_table4[suppl_table4$individual == "UPIC",]$gene), unique(suppl_table4[suppl_table4$individual == "nmXCI-1",]$gene))
overlap2 <- intersect(unique(suppl_table4[suppl_table4$individual == "UPIC",]$gene), unique(suppl_table4[suppl_table4$individual == "nmXCI-2",]$gene))
overlap3 <- intersect(unique(suppl_table4[suppl_table4$individual == "nmXCI-1",]$gene), unique(suppl_table4[suppl_table4$individual == "nmXCI-2",]$gene))

genes_overlap <- unique(c(overlap1, overlap2, overlap3))

# find overlapping tissues
overlap4 <- intersect(unique(suppl_table4[suppl_table4$individual == "UPIC",]$tissue_id), unique(suppl_table4[suppl_table4$individual == "nmXCI-1",]$tissue_id))
overlap5 <- intersect(unique(suppl_table4[suppl_table4$individual == "UPIC",]$tissue_id), unique(suppl_table4[suppl_table4$individual == "nmXCI-2",]$tissue_id))
overlap6 <- intersect(unique(suppl_table4[suppl_table4$individual == "nmXCI-1",]$tissue_id), unique(suppl_table4[suppl_table4$individual == "nmXCI-2",]$tissue_id))

UPIC_nmXCI1 <- na.omit(dcast(suppl_table4[suppl_table4$gene %in% overlap1 & suppl_table4$tissue_id %in% overlap4,], gene + tissue_id ~ individual, value.var = "allelic_expression"))
UPIC_nmXCI2 <- na.omit(dcast(suppl_table4[suppl_table4$gene %in% overlap2 & suppl_table4$tissue_id %in% overlap5,], gene + tissue_id ~ individual, value.var = "allelic_expression"))
nmXCI1_nmXCI2 <- na.omit(dcast(suppl_table4[suppl_table4$gene %in% overlap3 & suppl_table4$tissue_id %in% overlap6,], gene + tissue_id ~ individual, value.var = "allelic_expression"))

UPIC_nmXCI1$overlaping_between <- "UPIC vs nmXCI-1"
UPIC_nmXCI2$overlaping_between <- "UPIC vs nmXCI-2"
nmXCI1_nmXCI2$overlaping_between <- "nmXCI-1 vs nmXCI-2"

colnames(UPIC_nmXCI1) <- c("gene", "tissue_id", "participant1", "participant2", "overlaping_between")
colnames(UPIC_nmXCI2) <- c("gene", "tissue_id", "participant1", "participant2", "overlaping_between")
colnames(nmXCI1_nmXCI2) <- c("gene", "tissue_id", "participant1", "participant2", "overlaping_between")

df_overlaps <- rbind(UPIC_nmXCI1, UPIC_nmXCI2, nmXCI1_nmXCI2)

cor_test(data = df_overlaps, participant1, participant2, method = "spearman")

df_overlaps_with_class <- merge(df_overlaps, unique(dplyr::select(suppl_table4, new_category, gene)), by = "gene")

overlaping_between_order <- c("UPIC vs nmXCI-1", "UPIC vs nmXCI-2", "nmXCI-1 vs nmXCI-2")
df_overlaps_with_class$overlaping_between <- factor(df_overlaps_with_class$overlaping_between, levels = overlaping_between_order)

ggsave2(filename = paste(plot_dir, "Suppl_Fig_4.pdf"), width = 8, height = 8,
ggplot(df_overlaps_with_class, aes(x=participant1, y=participant2)) + 
  geom_point(aes(col = overlaping_between, shape = new_category)) + 
  geom_smooth(method = "lm", level=0.95, aes(col = overlaping_between)) +
  stat_cor(method = "spearman", aes(col = overlaping_between)) + 
  theme_AL_box(legend.position = "bottom")+
  guides(col=guide_legend(ncol = 1), shape=guide_legend(ncol = 1))
)

