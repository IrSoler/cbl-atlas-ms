# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Manuscript: Single cell landscape of sex differences in the progression of 
# multiple sclerosis
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Loading libraries --------------------------------------------------------

library(scater)
library(scran)
library(SingleCellExperiment)
library(ggplot2)
library(knitr)
library(BiocParallel)

# --- Loading data -------------------------------------------------------------
sce = readRDS("results/sce_norm.rds")

# Covariable combining technical batches
sce$batch = paste(sce$Capbatch, sce$Seqbatch, sep="_")

################################################################################
#   Variance of the log-counts                                                 #
################################################################################

# Modeling Gene variation
dec.sce = modelGeneVar(sce, BPPARAM=MulticoreParam(60), assay.type = "logcounts", block=sce$batch, density.weights=FALSE) #block the variance identificated in removal confounders

# Visualizing the fit
blocked.stats <- dec.sce$per.block

for (i in colnames(blocked.stats)) {
  png(file=paste0("./figures/4.Feature_Selection/var/var_sample_",i,".png"), width=6, height=4, units="in", res=600)
  current <- blocked.stats[[i]]
  plot(current$mean, current$total, main=i, pch=16, cex=0.5,
       xlab="Mean", ylab="Variance")
  curfit <- metadata(current)
  curve(curfit$trend(x), col='#798BFF', add=TRUE, lwd=2) 
  dev.off()
}

# Ordering by most interesting genes for inspection.
dec.sce[order(dec.sce$bio, decreasing=TRUE),]

################################################################################
#   Explore coeficient of variation                                            #
################################################################################

# Model coefficient of variation
dec.cv2 = modelGeneCV2(sce, BPPARAM=MulticoreParam(60), assay.type="logcounts") 
fit.cv2 = metadata(dec.cv2)

jpeg(file="./figures/4.Feature_Selection/Coeff_of_variation.jpeg", width=6, height=4, units="in", res=300)
plot(fit.cv2$mean, fit.cv2$cv2, log="xy",col="#440154")
curve(fit.cv2$trend(x), col="#B8DE29", add=TRUE, lwd=2)
dev.off()

# Ordering by most interesting genes for inspection.
dec.cv2[order(dec.cv2$ratio, decreasing=TRUE),]

################################################################################
#   Select High Variable Genes                                                 #
################################################################################

# Select 20% of most variable genes based on variance of the log-counts.
chosen = getTopHVGs(dec.sce, prop=0.2)

# Subset
sce.hvg = sce[chosen,]
dim(sce.hvg)

################################################################################
#   Save results                                                               #
################################################################################

saveRDS(sce.hvg,"results/sce_hvg_20.rds")
saveRDS(sce,"results/sce_feature_selection_20.rds")
write(rownames(sce.hvg), file = "results/hvg_20_list.txt", sep = "\n")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
