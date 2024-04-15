library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(rtracklayer)
library(GenomicRanges)
library(ggplotify)
library(karyoploteR)
source("r_sources/plot_parameters.R")

##### Read in data ##### 
# Supplementary table 4
suppl_table4 <- subset(fread("supplementary_tables/Supplementary_Table4.tsv"), individual != "UPIC")

# read in gencode v41
gtf_raw <- import("annotation_data/gencode.v41.annotation.gtf")
gtf_raw <- gtf_raw[seqnames(gtf_raw) == "chrX",]

####### figure 2c ######
# make zoom
entire_X_zoom.region <- toGRanges(data.frame("chrX", 1, 160000000))

# keep only gene markers we are interested in
markerz_our <- gtf_raw[gtf_raw$gene_name %in% unique(suppl_table4$gene) & gtf_raw$type == "gene",]

# add PAR tag
markerz_our$PAR <- ifelse(markerz_our$gene_name %in% suppl_table4[suppl_table4$PAR == "PAR",]$gene, yes = "PAR", no = "nonPAR")

plot_fin_alles <- as.ggplot(expression(
  kp <- plotKaryotype(plot.type=3, chromosomes = "chrX",  main = "gene overlap", genome = "hg38", zoom = entire_X_zoom.region),
  kpAddBaseNumbers(kp),
  kpArrows(kp, chr="chrX", x0=2784257, x1=2784257, y0=0.1, y1=0, cex=0.5, data.panel = 2),
  kpArrows(kp, chr="chrX", x0=155701383, x1=155701383, y0=0.1, y1=0, cex=0.5, data.panel = 2),
  kpPoints(kp, 
           data = markerz_our[markerz_our$PAR == "PAR",],
           chr = markerz_our[markerz_our$PAR == "PAR",]$seqnames,
           x=markerz_our[markerz_our$PAR == "PAR",]$start,
           cex = 0.5,
           y=0.1,
           data.panel = 1,
           col="#00A651"
  ),
  kpPoints(kp, 
           data = markerz_our[markerz_our$PAR == "nonPAR",],
           chr = markerz_our[markerz_our$PAR == "nonPAR",]$seqnames,
           x=markerz_our[markerz_our$PAR == "nonPAR",]$start,
           cex = 0.5,
           y=0.1,
           data.panel = 1
  )
  )
)


ggsave2(filename = paste0(plot_dir, "Fig_2B.pdf"), height = 14, width = 24, units = "cm",
        plot_grid(plot_fin_alles)
)

