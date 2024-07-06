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
library(cowplot)
library(ggplot2)
library(PCAtools)
library(uwot)
library(BiocParallel)
library(viridis)

# --- Loading data -------------------------------------------------------------

sce.hvg = readRDS("results/sce_hvg_20.rds")
sce.hvg$group_sex = factor(sce.hvg$group_sex, levels = c("Control_female", "MS_female", "Control_male", "MS_male"))

# Colors
fu = "#D46C40"
fd = "#FFC7AF"
mu = "#497EC2"
md = "#ADCBF8"
su = "#906BBA"
sd = "#CCB8E3"

numero_core = 20

################################################################################
#   Principal component analysis                                               #
################################################################################

## Run PCA
set.seed(100)
sce.hvg = runPCA(sce.hvg, exprs_values = "logcounts")

## Choose the number of PC using the elbow point
percent.var = attr(reducedDim(sce.hvg), "percentVar")
chosen.elbow = PCAtools::findElbowPoint(percent.var)

# --- Plots --------------------------------------------------------------------

## Percent variance
par(mar = c(5, 8, 2, 2))
plot(percent.var[0:10], axes = F, cex.lab = 2, xlab="PC", ylab="Percentage (%) of \n explained variability", 
     cex = 4, pch=21, col="#AF979B", bg="#AF979B", xlim=c(1,10), ylim=c(0,30))
axis(1, cex.axis = 2, at=0:10)
axis(2, cex.axis = 2)
axis.break(axis=1,breakpos)
abline(v=chosen.elbow, col="#891B2F", lwd = 8)

## Select PCs
reducedDim(sce.hvg, "PCA") = reducedDim(sce.hvg, "PCA")[,1:chosen.elbow]

## Plot PC1 vs PC2 (explore other covariables)
plotReducedDim(sce.hvg, dimred="PCA", colour_by="group_sex")+
  xlab("PC1")+
  ylab("PC2")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))

## Plot all selected PCs
plotReducedDim(sce.hvg, dimred="PCA", colour_by="group_sex", ncomponents = chosen.elbow)+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))

## Plot separating by group
plotReducedDim(sce.hvg, dimred="PCA", colour_by="group_sex", ncomponents = 1:2, point_alpha = 0.7) +
  facet_grid(sce.hvg$diagnosis ~ sce.hvg$sex, switch = "y", labeller = labeller(.rows = c("Control" = "Control", "MS" = "Multiple sclerosis"), .cols = c("female" = "Females", "male" = "Males"))) + 
  xlab("PC 1") + ylab("PC 2") + theme_bw(base_size = 18) +
  theme(legend.position = "bottom", legend.justification = "left")+ 
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))

################################################################################
#   UMAP                                                                       #
################################################################################

## Run UMAP
set.seed(100)
sce.hvg = runUMAP(sce.hvg, dimred="PCA", BPPARAM=MulticoreParam(numero_core), n_neighbors = 50)
reducedDim(sce.hvg, "UMAP_seed_100_n50") = reducedDim(sce.hvg, "UMAP")

# --- Plots --------------------------------------------------------------------
## Explore with other covariables

## All groups together
plotReducedDim(sce.hvg2, dimred="UMAP_seed_100_n50",colour_by="group_sex", point_alpha=0.3)+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))

## Separated groups
plotReducedDim(sce.hvg, dimred="UMAP_seed_100_n50", colour_by="group_sex", ncomponents = 1:2, point_alpha = 0.7) +
  facet_grid(sce.hvg$diagnosis ~ sce.hvg$sex, switch = "y", labeller = labeller(.rows = c("Control" = "Control", "MS" = "Multiple sclerosis"), .cols = c("female" = "Females", "male" = "Males"))) + 
  xlab("UMAP 1") + ylab("UMAP 2") + theme_bw() +
  theme(legend.position = "top")+ 
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


################################################################################
#   tSNE                                                                       #
################################################################################

## Run tSNE
set.seed(100)
sce.hvg = runTSNE(sce.hvg, dimred="PCA", perplexity = 40)
reducedDim(sce.hvg, "TSNE_seed_100_per_40") = reducedDim(sce.hvg, "TSNE")

# --- Plots --------------------------------------------------------------------
## Explore with other covariables

## All groups together
plotReducedDim(sce.hvg2, dimred="TSNE_seed_100_per_40",colour_by="group_sex", point_alpha=0.3)+
  xlab("tSNE 1")+
  ylab("tSNE 2")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))

## Separated groups
plotReducedDim(sce.hvg, dimred="TSNE_seed_100_per_40", colour_by="group_sex", ncomponents = 1:2, point_alpha = 0.7) +
  facet_grid(sce.hvg$diagnosis ~ sce.hvg$sex, switch = "y", labeller = labeller(.rows = c("Control" = "Control", "MS" = "Multiple sclerosis"), .cols = c("female" = "Females", "male" = "Males"))) + 
  xlab("UMAP 1") + ylab("UMAP 2") + theme_bw() +
  theme(legend.position = "top")+ 
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


################################################################################
#   Save results                                                               #
################################################################################

saveRDS(sce.hvg,"results/dim_reduc_20.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
