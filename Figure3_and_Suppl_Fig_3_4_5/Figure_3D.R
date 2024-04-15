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

##### Figure 3D #####
suppl_table_4_stats <- suppl_table_4 %>% dplyr::group_by(gene, new_category) %>% rstatix::get_summary_stats(allelic_expression, type = "common")

ggsave2(filename = paste0(plot_dir, "Fig_3D.pdf"), height = 12, width = 6, units = "cm",
        plot_grid(ncol=1,rel_heights = c(0.5,1),
                  ggplot(unique(dplyr::select(suppl_table_4,new_category, gene)), 
                         aes(x=factor(new_category, 
                                      levels = c("PAR","escape_across_tissues", "escape_in_single_tissue","inactive_across_tissues", "inactive_in_single_tissue",  "variable_across_tissues")), 
                             fill = factor(new_category,
                                           levels = c("PAR","escape_across_tissues", "escape_in_single_tissue","inactive_across_tissues", "inactive_in_single_tissue",  "variable_across_tissues")))) + 
                    geom_bar(stat="count", position = "dodge") + 
                    scale_fill_manual(values=c("PAR" = "#00A651", "escape_across_tissues" = "red", "escape_in_single_tissue" = "pink", "inactive_across_tissues" = "black", "inactive_in_single_tissue"="grey90", "variable_across_tissues" = "blue")) + 
                    geom_text(stat="count", aes(label=after_stat(count)), vjust= -.1, size = 1.5) + 
                    theme_AL_box_rotX()+ 
                    theme(strip.text.x = element_text(size=6, margin=margin(t=0, b=0)),                                                                                             
                          axis.text.x = element_text(size=6), 
                          axis.text.y = element_text(size=6), 
                          axis.title.x = element_text(size=6), 
                          axis.title.y = element_text(size=6),
                          legend.position = "none")+
                    labs(x="",y="gene count"),
                  
                  
                  ggplot(suppl_table_4_stats, 
                         aes(x=factor(new_category, 
                                      levels = c("PAR","escape_across_tissues", "escape_in_single_tissue","inactive_across_tissues", "inactive_in_single_tissue",  "variable_across_tissues")), 
                             y=median, 
                             fill = new_category)) + 
                    theme_AL_box_rotX()+ 
                    geom_boxplot(outlier.shape = NA)+
                    geom_quasirandom(dodge.width = 0.75, alpha = 1,size=0.5, 
                                     aes(col=new_category)) +
                    scale_fill_manual(values=c("PAR" = "#00A651", "escape_across_tissues" = "red", "escape_in_single_tissue" = "pink", "inactive_across_tissues" = "black", "inactive_in_single_tissue"="grey90", "variable_across_tissues" = "blue"))+
                    scale_color_manual(values=c("PAR" = "#00A651", "escape_across_tissues" = "red", "escape_in_single_tissue" = "pink", "inactive_across_tissues" = "black", "inactive_in_single_tissue"="grey90", "variable_across_tissues" = "blue"))+
                    theme(strip.text.x = element_text(size=6, margin=margin(t=0, b=0)),
                          axis.text.x = element_text(size=6), 
                          axis.text.y = element_text(size=6), 
                          axis.title.x = element_text(size=6), 
                          axis.title.y = element_text(size=6),
                          legend.title = element_blank(),
                          legend.text = element_text(size=5),
                          legend.key.height = unit(2, 'mm'),
                          legend.key.width = unit(2, 'mm'),
                          legend.position = "none")+
                    labs(x="",y="allelic expression per gene across tissues (median)")
                  
                  
        )
)
