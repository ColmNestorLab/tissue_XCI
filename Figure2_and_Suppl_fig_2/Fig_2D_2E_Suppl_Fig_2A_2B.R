library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(liftOver)
library(rtracklayer)
library(readxl)
library(ggpubr)
library(ggrepel)
library(VennDiagram)


# plot parameters
source("r_sources/plot_parameters.R")

# read in data #
suppl_table_4 <- fread(file = "supplementary_tables/Supplementary_Table4.tsv")

####### Tukiainen UPIC data reading in and pre-processing #######
# read in Tukiainen UPIC table
tuki_UPIC <- fread("annotation_data/Suppl.Table.5.csv")
tuki_UPIC_table_13 <- fread("annotation_data/landscape.Suppl.Table.13.csv")

# change tissue names to match ours
tuki_UPIC[tuki_UPIC$Tissue == "ESPMCS",]$Tissue <- "esophagus-muco"
tuki_UPIC[tuki_UPIC$Tissue == "STMACH",]$Tissue <- "stomach"
tuki_UPIC[tuki_UPIC$Tissue == "ESPMSL",]$Tissue <- "esophagus-musc"
tuki_UPIC[tuki_UPIC$Tissue == "ARTAORT",]$Tissue <- "artery-aort"
tuki_UPIC[tuki_UPIC$Tissue == "LUNG",]$Tissue <- "lung"
tuki_UPIC[tuki_UPIC$Tissue == "ARTCRN",]$Tissue <- "artery-coro"
tuki_UPIC[tuki_UPIC$Tissue == "SKINS",]$Tissue <- "skin-lleg"
tuki_UPIC[tuki_UPIC$Tissue == "UTERUS",]$Tissue <- "uterus"
tuki_UPIC[tuki_UPIC$Tissue == "LCL",]$Tissue <- "lymphocytes"
tuki_UPIC[tuki_UPIC$Tissue == "PNCREAS",]$Tissue <- "pancreas"
tuki_UPIC[tuki_UPIC$Tissue == "KDNCTX",]$Tissue <- "kidney-cort"
tuki_UPIC[tuki_UPIC$Tissue == "VAGINA",]$Tissue <- "vagina"
tuki_UPIC[tuki_UPIC$Tissue == "CLNTRN",]$Tissue <- "colon-tran"
tuki_UPIC[tuki_UPIC$Tissue == "WHLBLD",]$Tissue <- "whole blood"
tuki_UPIC[tuki_UPIC$Tissue == "LIVER",]$Tissue <- "liver"
tuki_UPIC[tuki_UPIC$Tissue == "THYROID",]$Tissue <- "thyroid"

# change colnames to line up better with our colnames
colnames(tuki_UPIC) <- c("gene", "gene_ID", "position_hg19", "tissue_id","refCount", "altCount", "totalCount","XaCount","XiCount","Xi_to_total_expr.","P_value", "Q_value", "incomplete_XCI")

# remove 'X:' from the hetSNP position_hg19
tuki_UPIC$position_hg19 <- gsub(tuki_UPIC$position_hg19, pattern = "X:", replacement = "")

# add contig
tuki_UPIC$contig <- "chrX"

# make short df
tuki_UPIC_short <- dplyr::select(tuki_UPIC, contig, position_hg19, gene, tissue_id)

# make position_hg19 numeric
tuki_UPIC_short$position_hg19 <- as.numeric(tuki_UPIC_short$position_hg19)

# add stop column
tuki_UPIC_short$stop <- tuki_UPIC_short$position_hg19+1

# make GRanges object
tuki_UPIC_short_GR <- makeGRangesFromDataFrame(tuki_UPIC_short,
                                               keep.extra.columns=TRUE,
                                               ignore.strand=TRUE,
                                               seqnames.field=c("contig"),
                                               start.field="position_hg19",
                                               end.field=c("end", "stop"))

# chain file from http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/
chain <- import.chain("annotation_data/hg19ToHg38.over.chain")

# liftOver to hg38
tuki_UPIC_short_GR_Hg38 <- liftOver(tuki_UPIC_short_GR, chain)

# convert to data frame and keep only columns I need
tuki_UPIC_short_GR_Hg38_DF <- dplyr::select(as.data.frame(tuki_UPIC_short_GR_Hg38), start, gene, tissue_id)

# rename columns
colnames(tuki_UPIC_short_GR_Hg38_DF) <- c("position_hg38", "gene", "tissue_id")

# merge with original data frame
tuki_UPIC_fin <- merge(tuki_UPIC, tuki_UPIC_short_GR_Hg38_DF, by = c("gene", "tissue_id"))

# calculate allelic expression
tuki_UPIC_fin$allelic_expression_tuki <- abs(0.5 - tuki_UPIC_fin$refCount / (tuki_UPIC_fin$refCount + tuki_UPIC_fin$altCount))

## load gene annotations
gtf_raw <- import("annotation_data/gencode.v41.annotation.gtf")

## add gene annotations to WES/ASE df
gr.WES <- with(tuki_UPIC_fin, GRanges(seqnames=contig, IRanges(position_hg38,width = 1)))
ol.WES <- findOverlaps(gr.WES, gtf_raw[gtf_raw$type %in% c("exon", "UTR")])
nm.WES <- tapply(gtf_raw[gtf_raw$type %in% c("exon", "UTR")]$gene_name[subjectHits(ol.WES)],queryHits(ol.WES),function(x) paste0(unique(x),collapse=";") )

tuki_UPIC_fin[, gene_new := NA]
tuki_UPIC_fin$gene_new[as.numeric(names(nm.WES))] <- nm.WES

# remove overlapping genes
tuki_UPIC_fin <- tuki_UPIC_fin[!is.na(tuki_UPIC_fin$gene_new) & !grepl(tuki_UPIC_fin$gene_new, pattern = ";"),]





####### Our UPIC data reading in and pre-processing #######
#Split the data table on individual.
suppl_table5_filtered_UPIC <- data.table(suppl_table_4[suppl_table_4$individual == "UPIC",])


# overlapping sites
samez <- unique(suppl_table5_filtered_UPIC[suppl_table5_filtered_UPIC$position %in% tuki_UPIC_fin$position_hg38,]$position)
diffz <- unique(suppl_table5_filtered_UPIC[!suppl_table5_filtered_UPIC$position %in% tuki_UPIC_fin$position_hg38,]$position)

# merge our and their data table
smash_UPIC <- merge(suppl_table5_filtered_UPIC, tuki_UPIC_fin, by = c("gene", "tissue_id"), all = F)

smash_UPIC$samez <- ifelse(smash_UPIC$position == smash_UPIC$position_hg38, yes = "same_hetSNP", no = "different_hetSNP")

smash_UPIC2 <- merge(suppl_table5_filtered_UPIC, tuki_UPIC_fin, by = c("gene", "tissue_id"), all = T)



####### figure 2E ######
# Generate 3 sets
set1_snps <- unique(tuki_UPIC_fin$position_hg38)
set2_snps <- unique(suppl_table5_filtered_UPIC$position)

# Chart
fig2E1 <- venn.diagram(
  x = list(set1_snps, set2_snps),
  category.names = c("Tukiainen" , "Gylemo"),
  filename = NULL,
  imagetype="tiff",
  main = "hetSNPs",
  fill = c("red", "yellow"), alpha = c(0.33, 0.33), cex = 1, cat.fontface = 2,
  lty =2
)

set1_genes <- unique(tuki_UPIC_fin$gene_new)
set2_genes <- unique(suppl_table5_filtered_UPIC$gene)

# Chart
fig2E2 <- venn.diagram(
  x = list(set1_genes, set2_genes),
  category.names = c("Tukiainen " , "Gylemo"),
  filename = NULL,
  imagetype="tiff",
  main = "genes",
  fill = c("red", "yellow"), alpha = c(0.33, 0.33), cex = 1,cat.fontface = 2,
  lty =2
)


ggsave2(filename = paste0(plot_dir, "Fig_2E.pdf"), height = 4,  width=8,
        plot_grid(fig2E1, fig2E2)
)

# differing genes
same_genes <- intersect(set1_genes, set2_genes)
unique_Tuki_genes <- setdiff(set1_genes, set2_genes)
unique_our_genes <- setdiff(set2_genes, set1_genes)



checking <- tuki_UPIC_fin[tuki_UPIC_fin$gene_new %in% unique_Tuki_genes,]



ggsave2(filename = paste0(plot_dir, "Fig_2D_and_Suppl_Fig_2A.pdf"), units = "mm", width = 400, height = 400,
        plot_grid(ncol=1, labels = c("A", "B"),
                  ggplot(smash_UPIC, aes(x=allelic_expression, y=allelic_expression_tuki)) +
                    geom_point(data=smash_UPIC[!smash_UPIC$gene %in% c("AKAP17A", "ARHGAP6"),]) + 
                    geom_point(data=smash_UPIC[smash_UPIC$gene %in% c("AKAP17A", "ARHGAP6"),], col = "red") +
                    geom_text(data=smash_UPIC[smash_UPIC$gene %in% c("AKAP17A", "ARHGAP6"),], aes(label=gene)) + 
                    #geom_abline(slope=1, intercept=0) +
                    geom_smooth(method = "lm", level=0.95) +
                    stat_cor(method = "spearman") + 
                    theme_AL_box()+
                    labs(title = "Gene AE correlation", subtitle = "spearman", y = "allelic expression (Tukiainen analysis)", x = "allelic expression (our analysis)")+
                    geom_hline(yintercept = c(0.1,0.4)) +
                    geom_vline(xintercept = c(0.1,0.4)) +
                    theme(legend.position = "top")+
                    scale_y_continuous(breaks = c(0,0.25,0.5))+
                    scale_x_continuous(breaks = c(0,0.25,0.5)),
                  
                  ggplot(smash_UPIC, aes(x=allelic_expression, y=allelic_expression_tuki)) + 
                    geom_point() + 
                    #geom_abline(slope=1, intercept=0) +
                    geom_smooth(method = "lm", level=0.95) +
                    stat_cor(method = "spearman") + 
                    theme_AL_box()+
                    #geom_hline(yintercept = 0.4) +
                    #geom_vline(xintercept = 0.4) +
                    labs(title = "Gene AE correlation tissue split", subtitle = "spearman", y = "allelic expression (Tukiainen analysis)", x = "allelic expression (our analysis)")+
                    facet_wrap(~tissue_id)+
                    guides(col=guide_legend(title="same hetSNP"))+
                    theme(legend.position = "top") +
                    scale_y_continuous(breaks = c(0,0.25,0.5))+
                    scale_x_continuous(breaks = c(0,0.25,0.5)))
)


# Supplementary Figure 2B
genes_unique_to_Tukiainen_UPIC <- fread("supplementary_tables/Supplementary_Table9.csv", header = F)

unique(genes_unique_to_Tukiainen_UPIC$V2)

genes_unique_to_Tukiainen_UPIC$new_category <- "SNP not found"
genes_unique_to_Tukiainen_UPIC[grepl(genes_unique_to_Tukiainen_UPIC$V2, pattern = "Too_low"),]$new_category <- "low readcount"
genes_unique_to_Tukiainen_UPIC[grepl(genes_unique_to_Tukiainen_UPIC$V2, pattern = "Intronic"),]$new_category <- "SNP intronic"
genes_unique_to_Tukiainen_UPIC[grepl(genes_unique_to_Tukiainen_UPIC$V2, pattern = "transformed_lymphocytes_which_is_excluded_in_our_analysis"),]$new_category <- "found in excluded tissue (EBV-transformed lymphycytes)"

ggsave2(filename = paste0(plot_dir, "Suppl_Fig_2B.pdf"), 
ggplot(genes_unique_to_Tukiainen_UPIC, aes(y=factor(new_category, levels = rev(c("found in excluded tissue (EBV-transformed lymphycytes)", "SNP not found", "SNP intronic",  "low readcount"))))) + geom_bar(stat="count") + geom_text(stat="count", aes(label=after_stat(count))) + theme_AL_box_rotX()
)

