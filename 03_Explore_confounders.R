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
library(batchelor)
library(BiocParallel)
library(ggplot2)
library(viridis)
library(MetBrewer)

# --- Loading data -------------------------------------------------------------
sce = readRDS("results/sce_norm.rds")

################################################################################
#   Exploratory variables                                                      #
################################################################################

# --- Explanatory variables by gene --------------------------------------------

## Calculates the percentage of variance of each gene’s expression that is 
## explained by each variable in the colData of the SingleCellExperiment object.
vars <- scater::getVarianceExplained(sce, exprs_values = "logcounts", variables=c("sample", "age", "nFeaturesAF", "nCountnorm", "group_sex", "Capbatch", "Seqbatch","CellCycle"))
head(vars)

## Plot
plotExplanatoryVariables(vars) +
  scale_color_manual(values=met.brewer(name = "Signac", n = 8), name = "Covariates", labels = c("Sample ID", "Normalised library size", "Expressed genes", "Capbatch", "Seqbatch", "Group", "Cell cycle stage", "Age"))


# --- Explanatory variables by principal component -----------------------------

## Identify PCs that correlate strongly to certain QC or metadata values
sce <- scater::runPCA(sce, exprs_values = "logcounts")
percent.var = attr(reducedDim(sce), "percentVar")

vars_pca <- scater::getExplanatoryPCs(sce, variables=c("sample", "age", "nFeaturesAF", "nCountnorm", "group_sex", "Capbatch", "Seqbatch","CellCycle"))
vars_pca_1 <- vars_pca / 100

## Plot
plotExplanatoryPCs(vars_pca_1) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  scale_color_manual(values=met.brewer(name = "Signac", n = 8), name = "Covariates", labels = c("Sample ID", "Normalised library size", "Expressed genes", "Capbatch", "Seqbatch", "Group", "Cell cycle stage", "Age"))


## --- Explore reported batch effect -------------------------------------------

plotReducedDim(sce, dimred="PCA", colour_by="Capbatch", ncomponents = 3)+
  ggtitle("PCA capture batch")+
  scale_color_viridis(option = 'viridis', discrete = TRUE)

plotReducedDim(sce, dimred="PCA", colour_by="Seqbatch", ncomponents = 3)+
  ggtitle("PCA sequencing batch")+
  scale_color_viridis(option = 'viridis', discrete = TRUE)
dev.off()

################################################################################
#   Exploratory variables adding cell type annotation                          #
################################################################################
# Execute this code after the cell annotation step

# --- Loading data -------------------------------------------------------------
sce = readRDS("results/sce_norm.rds")
sce.hvg = readRDS("results/cell_type_20_filtered.rds")

sce = sce[,colnames(sce.hvg)]
sce$cell.type = sce.hvg$cell.type

# --- Explanatory variables by gene --------------------------------------------

## Calculates the percentage of variance of each gene’s expression that is 
## explained by each variable in the colData of the SingleCellExperiment object
vars <- scater::getVarianceExplained(sce, exprs_values = "logcounts", variables=c("cell.type", "sample", "age", "nFeaturesAF", "nCountnorm", "group_sex", "Capbatch", "Seqbatch","CellCycle"))
head(vars)

## Plot
plotExplanatoryVariables(vars) +
  scale_color_manual(values=met.brewer(name = "Signac", n = 9), name = "Covariates", labels = c("Cell type", "Sample ID", "Normalised library size", "Expressed genes", "Capbatch", "Seqbatch", "Group", "Cell cycle stage", "Age"))

# --- Explanatory variables by principal component -----------------------------

## Identify PCs that correlate strongly to certain QC or metadata values
sce <- scater::runPCA(sce, exprs_values = "logcounts")
percent.var = attr(reducedDim(sce), "percentVar")

vars_pca <- scater::getExplanatoryPCs(sce, variables=c("cell.type", "sample", "age", "nFeaturesAF", "nCountnorm", "group_sex", "Capbatch", "Seqbatch","CellCycle"))
head(vars_pca)

vars_pca_1 <- vars_pca / 100
head(vars_pca_1)

## Plot
plotExplanatoryPCs(vars_pca_1) +
  scale_x_continuous(breaks=seq(0, 10, 1)) +
  theme_bw(base_size = 18) +
  scale_color_manual(values=met.brewer(name = "Signac", n = 9), name = "Covariates", labels = c("Cell type", "Sample ID", "Normalised library size", "Expressed genes", "Capbatch", "Seqbatch", "Group", "Cell cycle stage", "Age"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

