# read in libraries
library(QDNAseq)
library(QDNAseq.hg38)

# set bin size
bins <- getBinAnnotations(binSize=500, genome="hg38")

# read in data
readCounts <- binReadCounts(bins, path="/mnt/wwn-0x5000cca2b0dc9d33-part1/skewing_project/paper_2/data/cXCI_women/redo_SNP_calling/WGS/BAM_fin")

readCounts <- applyFilters(readCounts)
readCounts <- estimateCorrection(readCounts)
readCounts <- applyFilters(readCounts, chromosomes="Y")
copyNumbers <- correctBins(readCounts)

copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

plot(copyNumbersSegmented)

copyNumbersCalled <- callBins(copyNumbersSegmented)

# export this for supplementary figure 1B.
plot(copyNumbersCalled)

# Export these for figures 1D.
c17 <- chromosomes(copyNumbersCalled) == c(17)
plot(copyNumbersCalled[c17,])

# Manually export each plot (annoyingly enough).