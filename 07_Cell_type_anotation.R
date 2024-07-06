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
library(ggplot2)
library(BRETIGEA)
library(BiocParallel)
library(SingleR)
library(viridis)
library(pheatmap)
library(celda)


# --- Loading data -------------------------------------------------------------

sce.hvg = readRDS("results/clustering_20_k20.rds")


################################################################################
#   Brain samples                                                              #
################################################################################-

# --- Cell type identification -------------------------------------------------

sce.matrix = as.matrix(logcounts(sce.hvg))

ct_res = brainCells(sce.matrix, species = "human", nMarker = 60, scale = "TRUE", celltypes = c("ast", "mic", "neu", "oli", "opc"))
type = colnames(ct_res)[apply(ct_res,1,which.max)]
sce.hvg$cell.type = type

sce.hvg$cell.type = factor(sce.hvg$cell.type,  levels = c("ast","mic","neu","oli","opc"),labels = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes", "OPCs"))
                                                                                                  
# --- Filter cells not well annotated ------------------------------------------

sce.hvg$cell.type2 = factor(ifelse(sce.hvg$cluster %in% c(13, 1, 9), "Astrocytes", 
                          ifelse(sce.hvg$cluster %in% c(8,11,3,2,7,18,10,5), "Neurons",
                            ifelse(sce.hvg$cluster %in% c(6), "Oligodendrocytes",
                              ifelse(sce.hvg$cluster %in% c(4), "OPCs",
                                     ifelse(sce.hvg$cluster %in% c(14), "Microglia","N"))))))

for (i in which(sce.hvg$cell.type2 == "N")) {
  sce.hvg$cell.type2[i] = sce.hvg$cell.type[i]
}

sce.hvg$cell.type2 = factor(sce.hvg$cell.type2, levels = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes", "OPCs"))
sce.hvg2 = sce.hvg[,!(sce.hvg$cluster %in% c(12, 15, 17))]


################################################################################
#   Blood samples                                                              #
################################################################################-

# --- Cell type identification -------------------------------------------------

ref = celldex::MonacoImmuneData() # Reference cell types

pred = SingleR(test=sce.hvg, ref=ref, labels=ref$label.main)
sce.hvg$cell.type = pred$pruned.labels

# --- Filter cells not well annotated ------------------------------------------

## Remove not assigned cell types
a = which(is.na(sce.hvg$cell.type2))
sce.hvg = sce.hvg[,-a]
rm(a)

## Remove Basophils, Neutrophils, Progenitors and T cells not identified as
## CD4+ or CD8+
rm.cell = c("Basophils", "Neutrophils", "Progenitors", "T cells")
sce.hvg2 = sce.hvg[,!sce.hvg$cell.type2 %in% rm.cell]


################################################################################
#   Plots                                                                      #
################################################################################

# Brain samples
mycol = c("#957464", "#FFF175", "#A68FA6", "#2D5075", "#e8aeae")
mycell = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes", "OPCs")
mylabel = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes", "OPCs")

# Blood samples
# mycol = c("#D2A8C1", "#44959A", "#3D405B", "#E2D166", "#D7962C", "#81E8B7")
# mycell = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "Dendritic cells", "Monocytes", "NK-cells")
# mylabel = c("B cells", "CD4+ T cells", "CD8+ T cells", "Dendritic cells", "Monocytes", "NK cells")


# --- PCA ----------------------------------------------------------------------
plotReducedDim(sce.hvg2, dimred="PCA", colour_by="cell.type", point_size=0.7)+
  xlab("PCA1")+
  ylab("PCA2")+
  theme(legend.position = "top") +
  scale_color_manual(values = mycol, name = "Cell types")

# --- UMAP ---------------------------------------------------------------------
plotReducedDim(sce.hvg2, dimred="UMAP_seed_100_n50", colour_by="cell.type", point_size=0.7)+
  xlab("UMAP1")+
  ylab("UMAP2")+
  theme(legend.position = "top") +
  scale_color_manual(values = mycol, name = "Cell types")

# --- tSNE ---------------------------------------------------------------------
plotReducedDim(sce.hvg2, dimred="TSNE_seed_100_per_40", colour_by="cell.type", point_size=0.7)+
  xlab("tSNE1")+
  ylab("tSNE2")+
  theme(legend.position = "top") +
  scale_color_manual(values = mycol, name = "Cell types")



################################################################################
#   Save results                                                               #
################################################################################

saveRDS(sce.hvg2,"results/cell_type_20.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


