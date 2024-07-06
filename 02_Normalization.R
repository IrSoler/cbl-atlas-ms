# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Manuscript: Single cell landscape of sex differences in the progression of 
# multiple sclerosis
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Loading libraries --------------------------------------------------------

library(SingleCellExperiment)
library(scater)
library(scran)
library(SummarizedExperiment)
library(GenomicRanges)
library(stats4)
library(BiocGenerics)
library(parallel)
library(dplyr)
library(BiocParallel)
library(ggplot2)
library(viridis)

# --- Loading data -------------------------------------------------------------
sce <- readRDS("results/sce_qc.rds")

################################################################################
#   Normalization by deconvolution                                             #
################################################################################

# Method: https://rdrr.io/bioc/scran/man/computeSumFactors.html
set.seed(100)
quick_clusters = quickCluster(sce, BPPARAM=MulticoreParam(10), method = "igraph")
sce = computeSumFactors(sce, clusters = quick_clusters, BPPARAM=MulticoreParam(10))
sce = logNormCounts(sce, assay.type = "counts", size.factors = sce$sizeFactor, log=TRUE, pseudo.count=1)

# Explore results
summary(sizeFactors(sce))
sce@colData$deconvolution_size_factors <- sizeFactors(sce)

hist(sizeFactors(sce), xlab="Size factor")
hist(log10(sizeFactors(sce)), xlab="log10(Size facto)")

plot(sce$nCountAF, sizeFactors(sce), xlab="Library size", ylab="Size_factor")

################################################################################
#   Plot normalized data                                                       #
################################################################################

# Colors
fu = "#D46C40"
fd = "#FFC7AF"
mu = "#497EC2"
md = "#ADCBF8"
su = "#906BBA"
sd = "#CCB8E3"

# Parsing
sce@colData$nCountnorm = Matrix::colSums(logcounts(sce))
metadata_4 = as.data.frame(sce@colData)
metadata_4$group_sex = factor(metadata_4$group_sex, levels = c("Control_female", "MS_female", "Control_male", "MS_male"))

#  Normalized data explored by group_sex (explore other covariables)
metadata_4 %>%
  ggplot(aes(x=sample, y=nCountnorm, fill=group_sex)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(1000, 4000)) +
  ggtitle("After normalization")+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom", legend.justification = "left") +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  xlab("Sample ID")+
  ylab("Normalized library size")+
  scale_fill_manual(values = c(fd, fu, md, mu), name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))
dev.off()

################################################################################
#   Save results                                                               #
################################################################################

saveRDS(sce,"results/sce_norm.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


