library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(rstatix)
library(ggplotify)
library(karyoploteR)

# plot parameters
source("r_sources/plot_parameters.R")

# read in data #
suppl_table_1 <- fread(file = "supplementary_tables/Supplementary_Table1.tsv")
suppl_table_2 <- fread(file = "supplementary_tables/Supplementary_Table2.tsv")
suppl_table_3 <- fread(file = "supplementary_tables/Supplementary_Table3.tsv")

##### Figure 1b
# plot median allelic expression of X-linked genes excluding PAR and variable genes.
#selected_samples <- fread("supplementary_tables/selected_samples.tsv")

suppl_table_2_metrics <- suppl_table_2[suppl_table_2$gene %in% suppl_table_1$gene,] %>% dplyr::group_by(individual, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")
suppl_table_2_metrics$tissue_id <- factor(suppl_table_2_metrics$tissue_id)
suppl_table_2_metrics <- suppl_table_2_metrics[order(suppl_table_2_metrics$median, decreasing = F),]

###### Figure 1c ########
# calculate metrics for all non-hits (i.e. samples with a median chrX nonPAR nonHET allele-specific expression lower than 0.475)
df_non_hits <- suppl_table_2[suppl_table_2$gene %in% suppl_table_1$gene & !suppl_table_2$individual %in% c("13PLJ", "ZZPU", "UPIC")]

df_hits <- suppl_table_3

# bind 'em
df_smash <- rbind(df_hits, df_non_hits)

# add tag to separate screen hits and non-hits
df_smash$participant <- "mosaics"
df_smash[df_smash$individual %in% c("ZZPU","13PLJ","UPIC"),]$participant <- df_smash[df_smash$individual %in% c("ZZPU","13PLJ","UPIC"),]$individual


ggsave(filename = paste0(plot_dir, "Fig1B_1C.pdf"), width = 300, height = 240, unit = "mm",
       plot_grid(ncol=1, 
                 plot_grid(ncol=3,rel_widths = c(0.5,1), labels = c("A","B"),
                           NULL,
                           plot_grid(
                             ggplot(suppl_table_2_metrics, 
                                    aes(x=reorder(individual, median), y=median)) + 
                               geom_ribbon(aes(ymin = median-se, ymax = median+se, group = 1), alpha = 0.15)+
                               geom_point(aes(col = tissue_id, fill = tissue_id), size = 0.75) + 
                               geom_hline(yintercept = 0.475, lty = 2) + 
                               theme_AL_box_rotX() + 
                               theme(legend.title = element_blank(), 
                                     axis.text.x = element_blank(), 
                                     axis.text.y = element_text(size=6), 
                                     axis.title.x = element_text(size=6), 
                                     axis.title.y = element_text(size=6),
                                     legend.text = element_text(size=4),
                                     legend.key.height = unit(3, 'mm'))+
                               labs(x="participant", y="allelic expression (median)")+
                               scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits = c(0,0.5))+
                               scale_color_manual(values=color_valzz, limits = force)+
                               geom_text_repel(data=suppl_table_2_metrics[suppl_table_2_metrics$individual %in% c("UPIC", "ZZPU", "13PLJ"),], 
                                               aes(x=individual, y=median, label = individual), size = 2, nudge_x = -50, direction = "y")),
                           
                           plot_grid(
                             ggplot(df_smash, aes(x=factor(participant, levels = c("mosaics","13PLJ", "ZZPU", "UPIC")), y=allelic_expression)) + 
                               stat_summary(geom="crossbar", 
                                            fun.min = function(z) { quantile(z,0.25) },
                                            fun.max = function(z) { quantile(z,0.75) },
                                            fun = median) + 
                               stat_summary(geom="pointrange",
                                            fun = median, 
                                            aes(col = tissue_id), position = position_dodge(width = 0.35), size = 0.25) +
                               geom_hline(yintercept = 0.475, lty = 2) + 
                               theme_AL_box_rotX() + 
                               theme(axis.text.x = element_text(size=6), 
                                     axis.text.y = element_text(size=6), 
                                     axis.title.x = element_text(size=6), 
                                     axis.title.y = element_text(size=6),
                                     legend.title = element_blank(),
                                     legend.text = element_text(size=4),
                                     legend.key.height = unit(3, 'mm'))+
                               scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits = c(0,0.5))+
                               labs(x="participant", y="allelic expression (median)")+
                               scale_color_manual(values = color_valzz)+
                               guides(col = guide_legend(ncol = 1))))
       )      
)




