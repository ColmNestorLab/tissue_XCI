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


####### Figure 4B ###########

tukiainen <- fread("annotation_data/landscape.Suppl.Table.13.csv", header = T)
names(tukiainen)[names(tukiainen) == 'Gene name'] <- 'Gene_name'

# small fixes
tukiainen <- tukiainen[!(tukiainen$Gene_name == "IDS" & is.na(tukiainen$XCI_across_tissues)),]
tukiainen[tukiainen$Gene_name == "KAL1",]$Gene_name <- "ANOS1"
tukiainen[tukiainen$Gene_name == "HDHD1",]$Gene_name <- "PUDP"
tukiainen[tukiainen$Gene_name == "ARSE",]$Gene_name <- "ARSL"
tukiainen[tukiainen$Gene_name == "PHF16",]$Gene_name <- "JADE3"
tukiainen[tukiainen$Gene_name == "CXorf24",]$Gene_name <- "LINC01560"
tukiainen[tukiainen$Gene_name == "FAM155B",]$Gene_name <- "NALF2"
tukiainen[tukiainen$Gene_name == "CXorf56",]$Gene_name <- "STEEP1"
tukiainen[tukiainen$Gene_name == "FAM127C",]$Gene_name <- "RTL8B"
tukiainen[tukiainen$Gene_name == "FAM127B",]$Gene_name <- "RTL8A"
tukiainen[tukiainen$Gene_name == "LINC00087",]$Gene_name <- "SMIM10L2B"
tukiainen[tukiainen$Gene_name == "TAZ",]$Gene_name <- "TAFAZZIN"

##### Alluvial ######
# Keep only Tukiainen genes that we also find
df_classification <- dplyr::select(tukiainen[tukiainen$Gene_name %in% suppl_table_4$gene,], Gene_name, XCI_across_tissues)

# rename columns, add axes column (required for alluvial plots)
colnames(df_classification) <- c("gene", "category")
df_classification$axes <- "Tukiainen"

# change classification names
df_classification[(df_classification$category == "Partial, similar"),]$category <- "escape_across_tissues"
df_classification[(df_classification$category == "Partial, heterogeneous"),]$category <- "variable_across_tissues"
df_classification[(df_classification$category == "Full"),]$category <- "inactive_across_tissues"


# make the part of the alluvial df that's related to our data
df_gylemo <- unique(dplyr::select(suppl_table_4, gene, new_category))
df_gylemo$axes <- "Gylemo"
colnames(df_gylemo) <- c("gene", "category", "axes")

# ID genes that we ID but they don't
novelz <- data.frame(gene = setdiff(df_gylemo$gene, df_classification$gene),
                     category = "not_investigated")
colnames(novelz) <- c("gene", "category")
novelz$axes <- "Tukiainen"

# bind 'em
df_classification <- rbind(df_classification, novelz)

# bind 'em
smash_fin <- rbind(df_classification, df_gylemo)

# NA category means they either a) are not present in the Tukiainen analyis or b) are only ID'd in one tissue.
smash_fin[is.na(smash_fin$category),]$category <- "not_investigated"

library(ggalluvial)

#smash_fin <- unique(smash_fin[smash_fin$gene %in% PAR$`Approved symbol`,])

# check its in lodes form (required data format for alluvial plots)
is_lodes_form(smash_fin, key = "axes", value = "category", id = "gene")

# Change category to PAR if the gene is in the PAR reigon
smash_fin[smash_fin$axes == "Tukiainen" & smash_fin$gene %in% tukiainen[tukiainen$Region == "PAR"]$Gene_name,]$category <- "PAR"

# check categories
unique(smash_fin$category)

# level it
smash_fin$category <- factor(smash_fin$category, levels = c("PAR","not_investigated", "escape_across_tissues", "variable_across_tissues", "inactive_across_tissues", "escape_in_single_tissue", "inactive_in_single_tissue"))

# count genes per categories
asdfff <- smash_fin %>% dplyr::group_by(category, axes) %>% dplyr::count()

# change names
smash_fin[smash_fin$axes == "Tukiainen",]$axes <- "Previous"
smash_fin$axes <- factor(smash_fin$axes, levels = c("Previous", "Gylemo"))

#qqqq <- smash_fin %>% dplyr::group_by(gene, axes) %>% count()
#factor(smash_fin$category)

smash_fin_test <- data.frame(smash_fin)

ggsave2(filename = paste0(plot_dir, "Figure_4B.pdf"), width=3, height = 4,
        plot_grid(ggplot(smash_fin_test, aes(alluvium = factor(gene), x = axes, stratum = category)) + 
                    geom_flow(aes(fill=category)) +
                    geom_stratum(aes(fill=category))  + 
                    theme_AL_simple(legend.title = element_blank())+
                    theme(axis.text.x = element_text(size=6),
                          axis.text.y = element_text(size=6),
                          legend.text = element_text(size=6))+
                    labs(x="",y="gene count")+
                    geom_text(stat="count",aes(col= category, label=..count.., group = interaction(category, axes)), show.legend = F, position = "stack")+
                    scale_fill_manual(values = c("PAR" = "#00A651", "escape_across_tissues" = "red", "escape_in_single_tissue" = "pink", "inactive_across_tissues" = "black", "inactive_in_single_tissue"="grey90", "variable_across_tissues" = "blue", "not_investigated" = "yellow"))+
                    scale_color_manual(values = c("PAR" = "#00A651", "escape_across_tissues" = "red", "escape_in_single_tissue" = "pink", "inactive_across_tissues" = "black", "inactive_in_single_tissue"="grey90", "variable_across_tissues" = "blue", "not_investigated" = "yellow"))
        )
)

