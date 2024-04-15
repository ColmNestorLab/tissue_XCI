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
unique(suppl_table_4$new_category)


tuki_table1 <- fread("annotation_data/landscape.Suppl.Table.13.csv", header = T)
tuki_table2 <- fread("annotation_data/Suppl.Table.5.csv", header = T)

tuki_table2_count <- tuki_table2 %>% dplyr::group_by(`Gene name`) %>% dplyr::count()
tuki_table2_count <- tuki_table2_count[order(tuki_table2_count$n, decreasing = F),]

plot_suppl <- suppl_table_4[suppl_table_4$gene %in% tuki_table2_count[tuki_table2_count$n == 1,]$`Gene name`,]
plot_suppl_with_their_classification <- merge(plot_suppl, tuki_table2[,c("Gene name", "Incomplete XCI")],by.x = c("gene") ,by.y = c("Gene name"))


plot_order <- plot_suppl_with_their_classification %>% dplyr::group_by(gene, new_category) %>% dplyr::summarise(medianz = median(allelic_expression))
plot_order <- plot_order[order(plot_order$medianz, decreasing = T),]



ggsave2(filename = paste0(plot_dir, "Figure_4A.pdf"), height = 3,  width=6,
        plot_grid(
          ggplot(plot_suppl_with_their_classification[plot_suppl_with_their_classification$new_category != "inactive_in_single_tissue",],
                 aes(x=factor(gene, levels = unique(plot_order$gene)), y=allelic_expression)) + 
            geom_boxplot(outlier.shape = NA) + 
            geom_quasirandom(aes(col=tissue_id))+ 
            theme_AL_box_rotX()+
            scale_color_manual(values=color_valzz)+
            scale_y_continuous(limits = c(0,NA))+
            facet_grid(~new_category, space = "free", scales = "free_x")+
            labs(x="", y="allelic expression")+
            geom_hline(yintercept = 0.4, lty=2)+
            guides(col=guide_legend(ncol = 1))
            
        )
)
