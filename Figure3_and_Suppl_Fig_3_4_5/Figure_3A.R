library(data.table)
library(cowplot)
library(VennDiagram)
source("r_sources/plot_parameters.R")

##### Read in data ##### 
# Supplementary table 5
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")

# keep only highest covered hSNP.
suppl_table4_filtered <- data.table(suppl_table4)

####### figure 2a and figure 2b ######
# Generate 3 sets
set1_snps <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "nmXCI-1",]$position)
set2_snps <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "nmXCI-2",]$position)
set3_snps <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$position)

# Chart
fig3a <- venn.diagram(
  x = list(set1_snps, set2_snps, set3_snps),
  category.names = c("nmXCI-1" , "nmXCI-2", "UPIC"),
  filename = NULL,
  imagetype="tiff",
  main = paste("All hSNPs =", length(unique(suppl_table4_filtered$position))),
  fill = c("red", "yellow", "blue"), alpha = c(0.33, 0.33, 0.33), cex = 1, cat.fontface = 2,
  lty =2
)

set1_genes <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "nmXCI-1",]$gene)
set2_genes <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "nmXCI-2",]$gene)
set3_genes <- unique(suppl_table4_filtered[suppl_table4_filtered$individual == "UPIC",]$gene)

# Chart
fig3b <- venn.diagram(
  x = list(set1_genes, set2_genes, set3_genes),
  category.names = c("nmXCI-1 " , "nmXCI-2", "UPIC"),
  filename = NULL,
  imagetype="tiff",
  main = paste("All genes =", length(unique(suppl_table4_filtered$gene))),
  fill = c("red", "yellow", "blue"), alpha = c(0.33, 0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2
)


ggsave2(filename = paste0(plot_dir, "Fig_3A.pdf"), height = 8,  width=4,
        plot_grid(fig3a, fig3b, ncol=1)
)


