library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(ggplotify)
library(grid)
library(gridExtra)
library(raster)
library(rstatix)
library(ggpubr)
library(karyoploteR)
library(ggbeeswarm)
library(VennDiagram)

source("r_sources/plot_parameters.R")

# Supplementary table 4
suppl_table_4 <- fread("supplementary_tables/Supplementary_Table4.tsv")
suppl_table_4$sig <- ifelse(suppl_table_4$binomP.adj < 0.01, yes = "sig", no = "ns")

##### Figure 3B #####
# example_genes
example_genes <- c("CSF2RA", 
                   "DDX3X", "PUDP", "ZFX",
                   "IDS", "TCEAL4", 
                   "ANOS1", "ARSD", "PRKX",
                   "ACOT9",
                   "AKAP17A", "GTPBP6"
)



df_fig3b <- suppl_table_4[suppl_table_4$gene %in% example_genes,]

df_fig3b_split1 <- df_fig3b[df_fig3b$gene %in% c("ACOT9", "AKAP17A", "GTPBP6", "ANOS1", "ARSD", "PRKX"),]
df_fig3b_split2 <- df_fig3b[df_fig3b$gene %in% c("CSF2RA", "IDS", "TCEAL4", "DDX3X", "PUDP", "ZFX"),]

df_fig3b_split1$category <- "kek"
df_fig3b_split1[df_fig3b_split1$gene %in% "ACOT9",]$category <- "inactive"
df_fig3b_split1[df_fig3b_split1$gene %in% "AKAP17A",]$category <- "escape"
df_fig3b_split1[df_fig3b_split1$gene %in% "GTPBP6",]$category <- "escape"
df_fig3b_split1[df_fig3b_split1$gene %in% "ANOS1",]$category <- "variable"
df_fig3b_split1[df_fig3b_split1$gene %in% "ARSD",]$category <- "variable"
df_fig3b_split1[df_fig3b_split1$gene %in% "PRKX",]$category <- "variable"

df_fig3b_split2$category <- "kek"
df_fig3b_split2[df_fig3b_split2$gene %in% "CSF2RA",]$category <- "variable"
df_fig3b_split2[df_fig3b_split2$gene %in% "IDS",]$category <- "variable"
df_fig3b_split2[df_fig3b_split2$gene %in% "TCEAL4",]$category <- "variable"
df_fig3b_split2[df_fig3b_split2$gene %in% "DDX3X",]$category <- "variable"
df_fig3b_split2[df_fig3b_split2$gene %in% "PUDP",]$category <- "variable"
df_fig3b_split2[df_fig3b_split2$gene %in% "ZFX",]$category <- "variable"


ggsave2(filename = paste0(plot_dir, "Fig_3C.pdf"), height = 20, width = 10, units = "cm",
        plot_grid(ncol=1,
                  ggplot(df_fig3b_split1, aes(x=factor(gene, levels = c("ACOT9", "AKAP17A", "GTPBP6", "ANOS1", "ARSD", "PRKX")), y=allelic_expression, shape=sig, col=sig)) + 
                    geom_point() +
                    facet_wrap(factor(category, levels = c("inactive","escape", "variable"))~factor(new_category, levels = c("PAR","escape_across_tissues", "inactive_across_tissues", "variable_across_tissues")), scales = "free_x", nrow = 1)+
                    labs(title="top facet = binomial test\nbottom facet = manual curation")+
                    theme_AL_box_rotX()+
                    scale_color_manual(values = c("ns" = "black", sig = "red"))+
                    geom_hline(yintercept = 0.4, lty=2)+
                    labs(x=""),
                  
                  ggplot(df_fig3b_split2, aes(x=factor(gene, levels = c("CSF2RA", "IDS", "TCEAL4", "DDX3X", "PUDP", "ZFX")), y=allelic_expression, shape=sig, col=sig)) + 
                    geom_point() + 
                    facet_wrap(factor(category, levels = c("variable","inactive","escape"))~factor(new_category, levels = c("PAR", "inactive_across_tissues","escape_across_tissues", "variable_across_tissues")), scales = "free_x", nrow = 1)+
                    labs(title="top facet = binomial test\nbottom facet = manual curation")+
                    theme_AL_box_rotX()+
                    scale_color_manual(values = c("ns" = "black", sig = "red"))+
                    geom_hline(yintercept = 0.4, lty=2)+
                    labs(x="")
        )
)
