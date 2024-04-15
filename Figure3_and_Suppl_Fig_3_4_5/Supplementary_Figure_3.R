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

suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")
suppl_table4$sig <- ifelse(suppl_table4$binomP.adj < 0.01, yes = "sig", no = "ns")
suppl_table5 <- fread("supplementary_tables/Supplementary_Table5.tsv")

data_wide <- merge(suppl_table4, suppl_table5, by = "gene")
data_wide$tissue_individual <- paste0(data_wide$individual, "_", data_wide$tissue_id)
data_wide$minor_allele <- pmin(data_wide$refCount, data_wide$altCount)
data_wide$major_allele <- pmax(data_wide$refCount, data_wide$altCount)

data_long <- melt(data_wide,
                  # ID variables - all the variables to keep but not split apart on
                  id.vars=c("gene", "tissue_individual", "curation", "curation_reason","curation_reason_short", "tissue_id", "sig"),
                  # The source columns
                  measure.vars=c("minor_allele", "major_allele", "totalCount"),
                  # Name of the destination column that will identify the original
                  # column that the measurement came from
                  variable.name="condition",
                  value.name="measurement"
)

data_long$sig_lab <- ifelse(data_long$sig == "ns", yes = NA, no = "*")


ggsave2(filename = paste0(plot_dir, "Suppl_figure_3.pdf"), height = 12, width = 40, limitsize = F,
        ggplot(data_long[data_long$condition != "totalCount" & data_long$gene %in% suppl_table5$gene,],
               aes(x=tissue_individual, y=measurement, fill = condition)) + 
          geom_col()+ 
          facet_wrap(interaction(curation_reason_short, curation)~gene, scales = "free", ncol = 16) + 
          theme_AL_box_rotX() +
          geom_text(data=data_long[data_long$condition == "totalCount" & data_long$gene %in% suppl_table5$gene,], aes(label=sig_lab), position = "stack", size = 8)+
          labs(title = "low_power")
)
