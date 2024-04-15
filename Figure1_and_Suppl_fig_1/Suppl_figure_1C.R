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
library(ggbeeswarm)
source("r_sources/plot_parameters.R")


### experimental -  FIN ###
# showing chr17 expression 

tpmz <- melt(fread("supplementary_tables/Supplementary_Table8.tsv"))

gtf <- import("annotation_data/gencode.v41.annotation.gtf")

gtf_filt <- subset(data.frame(gtf), seqnames %in% c("chr17") & type == "gene" & gene_type == "protein_coding")

tpmz_filt <- tpmz[tpmz$gene %in% gtf_filt$gene_name,]

tpmz_filt$contig <- "LOL"
tpmz_filt[tpmz_filt$gene %in% gtf_filt[gtf_filt$seqnames == "chr17",]$gene_name]$contig <- "chr17"
tpmz_filt$variable <- gsub(tpmz_filt$variable, pattern = "\\.", replacement = "-")

# add metadata (tissue)
sample_meta <- fread("annotation_data/sample.tsv")

# add UPIC meta (it is not included in the above metadata)
UPIC_meta <- dplyr::select(subset(fread("annotation_data/fixed_GTEx_Data_V6_Annotations_SampleAttributesDS.txt"), PARTICIPANT == "GTEX-UPIC" & SMGEBTCHT == "TrueSeq.v1"), SAMPID, SMTSD)

#change colnames to match
colnames(UPIC_meta) <- c("entity:sample_id", "tissue_id")

# bind together
meta_fin <- rbind(dplyr::select(sample_meta, "entity:sample_id", tissue_id), UPIC_meta)

meta_fin[meta_fin$tissue_id == "Adipose_Subcutaneous",]$tissue_id <- "adipose-subc"
meta_fin[meta_fin$tissue_id == "Adipose_Visceral_Omentum",]$tissue_id <- "adipose-visc"
meta_fin[meta_fin$tissue_id == "Artery_Aorta",]$tissue_id <- "artery-aort"
meta_fin[meta_fin$tissue_id == "Artery - Aorta",]$tissue_id <- "artery-aort"
meta_fin[meta_fin$tissue_id == "Artery_Coronary",]$tissue_id <- "artery-coro"
meta_fin[meta_fin$tissue_id == "Artery - Coronary",]$tissue_id <- "artery-coro"
meta_fin[meta_fin$tissue_id == "Artery_Tibial",]$tissue_id <- "artery-tibi"
meta_fin[meta_fin$tissue_id == "Brain_Amygdala",]$tissue_id <- "brain-amyg"
meta_fin[meta_fin$tissue_id == "Brain_Anterior_cingulate_cortex_BA24",]$tissue_id <- "brain-ante"
meta_fin[meta_fin$tissue_id == "Brain_Caudate_basal_ganglia",]$tissue_id <- "brain-caud"
meta_fin[meta_fin$tissue_id == "Brain_Cerebellar_Hemisphere",]$tissue_id <- "brain-cehe"
meta_fin[meta_fin$tissue_id == "Brain_Cerebellum",]$tissue_id <- "brain-cere"
meta_fin[meta_fin$tissue_id == "Brain_Cortex",]$tissue_id <- "brain-cort"
meta_fin[meta_fin$tissue_id == "Brain_Frontal_Cortex_BA9",]$tissue_id <- "brain-frco"
meta_fin[meta_fin$tissue_id == "Brain_Hippocampus",]$tissue_id <- "brain-hipp"
meta_fin[meta_fin$tissue_id == "Brain_Hypothalamus",]$tissue_id <- "brain-hypo"
meta_fin[meta_fin$tissue_id == "Brain_Nucleus_accumbens_basal_ganglia",]$tissue_id <- "brain-nucl"
meta_fin[meta_fin$tissue_id == "Brain_Putamen_basal_ganglia",]$tissue_id <- "brain-puta"
meta_fin[meta_fin$tissue_id == "Brain_Spinal_cord_cervical_c-1",]$tissue_id <- "brain-spin"
meta_fin[meta_fin$tissue_id == "Brain_Substantia_nigra",]$tissue_id <- "brain-subs"
meta_fin[meta_fin$tissue_id == "Breast_Mammary_Tissue",]$tissue_id <- "breast"
meta_fin[meta_fin$tissue_id == "Minor_Salivary_Gland",]$tissue_id <- "salivary gland"
meta_fin[meta_fin$tissue_id == "Cervix_Ectocervix",]$tissue_id <- "cervix-ecto"
meta_fin[meta_fin$tissue_id == "Cervix_Endocervix",]$tissue_id <- "cervix-endo"
meta_fin[meta_fin$tissue_id == "Colon_Sigmoid",]$tissue_id <- "colon-sigm"
meta_fin[meta_fin$tissue_id == "Colon_Transverse",]$tissue_id <- "colon-tran"
meta_fin[meta_fin$tissue_id == "Colon - Transverse",]$tissue_id <- "colon-tran"
meta_fin[meta_fin$tissue_id == "Esophagus_Gastroesophageal_Junction",]$tissue_id <- "esophagus-gaju"
meta_fin[meta_fin$tissue_id == "Esophagus_Mucosa",]$tissue_id <- "esophagus-muco"
meta_fin[meta_fin$tissue_id == "Esophagus - Mucosa",]$tissue_id <- "esophagus-muco"
meta_fin[meta_fin$tissue_id == "Esophagus_Muscularis",]$tissue_id <- "esophagus-musc"
meta_fin[meta_fin$tissue_id == "Esophagus - Muscularis",]$tissue_id <- "esophagus-musc"
meta_fin[meta_fin$tissue_id == "Cells_Cultured_fibroblasts",]$tissue_id <- "fibroblasts"
meta_fin[meta_fin$tissue_id == "Heart_Atrial_Appendage",]$tissue_id <- "heart-atri"
meta_fin[meta_fin$tissue_id == "Heart_Left_Ventricle",]$tissue_id <- "heart-vent"
meta_fin[meta_fin$tissue_id == "Kidney_Cortex",]$tissue_id <- "kidney-cort"
meta_fin[meta_fin$tissue_id == "Kidney - Cortex",]$tissue_id <- "kidney-cort"
meta_fin[meta_fin$tissue_id == "Kidney_Medulla",]$tissue_id <- "kidney-medu"
meta_fin[meta_fin$tissue_id == "Cells_EBV-transformed_lymphocytes",]$tissue_id <- "lymphocytes"
meta_fin[meta_fin$tissue_id == "Cells - EBV-transformed lymphocytes",]$tissue_id <- "lymphocytes"
meta_fin[meta_fin$tissue_id == "Muscle_Skeletal",]$tissue_id <- "muscle"
meta_fin[meta_fin$tissue_id == "Nerve_Tibial",]$tissue_id <- "nerve"
meta_fin[meta_fin$tissue_id == "Skin_Not_Sun_Exposed_Suprapubic",]$tissue_id <- "skin-supr"
meta_fin[meta_fin$tissue_id == "Skin_Sun_Exposed_Lower_leg",]$tissue_id <- "skin-lleg"
meta_fin[meta_fin$tissue_id == "Skin - Sun Exposed (Lower leg)",]$tissue_id <- "skin-lleg"
meta_fin[meta_fin$tissue_id == "Small_Intestine_Terminal_Ileum",]$tissue_id <- "small intestine"
meta_fin[meta_fin$tissue_id == "Adrenal_Gland",]$tissue_id <- "adrenal gland"
meta_fin[meta_fin$tissue_id == "Bladder",]$tissue_id <- "bladder"
meta_fin[meta_fin$tissue_id == "Liver",]$tissue_id <- "liver"
meta_fin[meta_fin$tissue_id == "Lung",]$tissue_id <- "lung"
meta_fin[meta_fin$tissue_id == "Pancreas",]$tissue_id <- "pancreas"
meta_fin[meta_fin$tissue_id == "Pituitary",]$tissue_id <- "pituitary"
meta_fin[meta_fin$tissue_id == "Spleen",]$tissue_id <- "spleen"
meta_fin[meta_fin$tissue_id == "Stomach",]$tissue_id <- "stomach"
meta_fin[meta_fin$tissue_id == "Thymocytes",]$tissue_id <- "thymocytes"
meta_fin[meta_fin$tissue_id == "Thyroid",]$tissue_id <- "thyroid"
meta_fin[meta_fin$tissue_id == "Whole_Blood",]$tissue_id <- "whole blood"
meta_fin[meta_fin$tissue_id == "Whole Blood",]$tissue_id <- "whole blood"
meta_fin[meta_fin$tissue_id == "Uterus",]$tissue_id <- "uterus"
meta_fin[meta_fin$tissue_id == "Vagina",]$tissue_id <- "vagina"
meta_fin[meta_fin$tissue_id == "Ovary",]$tissue_id <- "ovary"
meta_fin[meta_fin$tissue_id == "Fallopian_Tube",]$tissue_id <- "fallopian tube"

meta_fin$individual <- "lolz"

meta_fin[grepl(meta_fin$`entity:sample_id`, pattern = "UPIC"),]$individual <- "UPIC"
meta_fin[grepl(meta_fin$`entity:sample_id`, pattern = "13PLJ"),]$individual <- "nmXCI-1"
meta_fin[grepl(meta_fin$`entity:sample_id`, pattern = "ZZPU"),]$individual <- "nmXCI-2"

tpmz_filt_meta <- merge(tpmz_filt, meta_fin, by.x = "variable", by.y = "entity:sample_id")


chr17_tpmz <- tpmz_filt_meta[tpmz_filt_meta$contig == "chr17" ,]

chr17_tpmz$arm <- ifelse(chr17_tpmz$gene %in% gtf_filt[gtf_filt$gene_name %in% chr17_tpmz$gene & gtf_filt$end < 25000000,]$gene_name, yes = "p", no = "q")

ggsave2(
  paste0(plot_dir, "Sfig_1C.pdf"),
  ggplot(chr17_tpmz[chr17_tpmz$individual  %in% c("nmXCI-2","UPIC"),], 
         aes(x=tissue_id, y=value, fill=arm)) + 
    geom_boxplot(outlier.shape = NA, alpha = .4) + 
    theme_AL_box_rotX() + 
    facet_grid(~individual, scales = "free_x", space = "free_x") + 
    coord_cartesian(ylim=c(0,70))
)

