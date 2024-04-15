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
suppl_table4 <- fread("supplementary_tables/Supplementary_Table4.tsv")
suppl_table4$sig <- ifelse(suppl_table4$binomP.adj < 0.01, yes = "sig", no = "ns")

# Supplementary table 5
suppl_table5 <- fread("supplementary_tables/Supplementary_Table5.tsv")

# merge individual and tissues to single-column
suppl_table4$individual_tissue <- paste0(suppl_table4$individual, "_", suppl_table4$tissue_id)

# Set factor orders
individual_tissue_order <- unique(sort(suppl_table4$individual_tissue))
new_category_order <- c("PAR","escape_across_tissues", "escape_in_single_tissue","inactive_across_tissues", "inactive_in_single_tissue",  "variable_across_tissues")
new_category_color_orderz <- c("PAR" = "#00A651", "escape_across_tissues" = "red", "escape_in_single_tissue" = "pink", "inactive_across_tissues" = "black", "inactive_in_single_tissue"="grey90", "variable_across_tissues" = "blue")

# add binary read outs for easier plotting
#suppl_table4$binary_readout_AE <- ifelse(suppl_table4$allelic_expression <= 0.4, yes = "escape", no = "inactive")
suppl_table4$binary_readout_binomP <- ifelse(suppl_table4$binomP.adj < 0.01, yes = "*", no = NA)

# p- and q-arm
suppl_table4$arm <- ifelse(suppl_table4$position < 61000001, yes = "p", no = "q")

# check the arms
armz <- unique(suppl_table4[,c("arm", "position", "gene")])

# make complete df to have missing data as grey
suppl_table4_complete <- tidyr::complete(suppl_table4, individual_tissue, gene)
#suppl_table4_complete[is.na(suppl_table4_complete$allelic_expression),]$binary_readout_AE <- "NA"

# remove superfluous columns
suppl_table4_complete_short <- unique(dplyr::select(suppl_table4, gene, new_category ))

# add placeholder to make both heatmaps equally long
suppl_table4_complete_short$placeholder <- "nmXCI_2_esophagus_musc"

# read in annotation
chrX_gencode <- subset(data.frame(import("annotation_data/gencode.v41.annotation.gtf")), seqnames == "chrX")

# add asterix to curated genes
suppl_table4_complete_short$curation <- ifelse(suppl_table4_complete_short$gene %in% suppl_table5$gene, yes = "*", no = NA)

# get order of X-linked genes for nicer plotting
chrX_gene_order <- suppl_table4
chrX_gene_order <- unique(chrX_gene_order[,c("gene", "position")][order(chrX_gene_order$position, decreasing = F),]$gene)

# Sanity check to see we haven't lost any genes along the way
intersect(suppl_table4_complete_short$gene, chrX_gene_order)



# plot
ggsave2(paste0(plot_dir, "Suppl_Figure_5.pdf"), height = 20, width = 60, limitsize = F, 
        plot_grid(ncol=1,rel_heights = c(1, -0.675, 1), align = "h",
                  plot_grid(
                    ggplot(suppl_table4_complete, 
                           aes(x=factor(gene, levels = chrX_gene_order), 
                               y=factor(individual_tissue, levels = rev(individual_tissue_order)), 
                               fill = as.numeric(allelic_expression))) + 
                      geom_tile(col="white") + 
                      #geom_tile() + 
                      geom_text(aes(label=binary_readout_binomP, y = factor(individual_tissue, levels = rev(individual_tissue_order))), nudge_y = -0.275, col = "black", size = 5) +
                      #scale_fill_gradient() +
                      #scale_fill_steps(low="beige", high = "red") +
                      scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.4) +
                      #scale_fill_manual(values = c( "NA" = "grey", "escape" = "red", "inactive" = "black")) + 
                      coord_fixed(1/1) + 
                      theme_AL_box_rotX(legend.position="top", legend.title = element_blank()) +
                      #guides(fill=guide_legend(ncol = 3))+
                      theme(axis.text.x = element_blank()) +
                      labs(x="", y=""))
                  ,
                  
                  plot_grid(NULL)
                  
                  ,
                  
                  plot_grid(
                    ggplot(suppl_table4_complete_short, 
                           aes(x=factor(gene, levels = chrX_gene_order), y=placeholder, fill = factor(new_category, levels = new_category_order))) + 
                      coord_fixed(1/1) + 
                      geom_tile(col="black") + 
                      scale_fill_manual(values=new_category_color_orderz) +
                      geom_text(aes(label=curation), col = "white") +
                      theme_AL_box_rotX(legend.position="bottom", legend.title = element_blank()) +
                      guides(fill=guide_legend(ncol =6))+
                      labs(x="", y="")
                  )
        )
)



ggsave2(paste0(plot_dir, "Figure_3E.pdf"), height = 20, width = 20, limitsize = F, 
        plot_grid(ncol=1,rel_heights = c(1, -0.5275, 1), align = "h",
                  plot_grid(
                    ggplot(suppl_table4_complete[!suppl_table4_complete$gene %in% unique(suppl_table4[suppl_table4$new_category %in% c("inactive_across_tissues", "inactive_in_single_tissue"),]$gene),], 
                           aes(x=factor(gene, levels = chrX_gene_order), 
                               y=factor(individual_tissue, levels = rev(individual_tissue_order)), 
                               fill = as.numeric(allelic_expression))) + 
                      geom_tile(col="white") + 
                      #geom_tile() + 
                      geom_text(aes(label=binary_readout_binomP, y = factor(individual_tissue, levels = rev(individual_tissue_order))), nudge_y = -0.275, col = "black", size = 5) +
                      #scale_fill_gradient() +
                      #scale_fill_steps(low="beige", high = "red") +
                      scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0.4, na.value = "grey") +
                      #scale_fill_manual(values = c( "NA" = "grey", "escape" = "red", "inactive" = "black")) + 
                      coord_fixed(1/1) + 
                      theme_AL_box_rotX(legend.position="top", legend.title = element_blank()) +
                      #guides(fill=guide_legend(ncol = 3))+
                      theme(axis.text.x = element_blank()) +
                      labs(x="", y="")
                    
                    )
                  ,
                  plot_grid(NULL)
                  ,
                  
                  plot_grid(
                    ggplot(suppl_table4_complete_short[!suppl_table4_complete_short$new_category %in% c("inactive_across_tissues", "inactive_in_single_tissue"),],
                           aes(x=factor(gene, levels = chrX_gene_order), y=placeholder, fill = factor(new_category, levels = c("PAR", "escape_across_tissues", "escape_in_single_tissue", "variable_across_tissues")))) + 
                      coord_fixed(1/1) + 
                      geom_tile(col="black") + 
                      scale_fill_manual(values=new_category_color_orderz) +
                      geom_text(aes(label=curation), col = "white") +
                      theme_AL_box_rotX(legend.position="bottom", legend.title = element_blank()) +
                      guides(fill=guide_legend(ncol =6))+
                      labs(x="", y="")
                  )
        )
)
