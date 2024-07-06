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
library(umap)
library(pheatmap)
library(igraph)
library(cowplot)
library(SingleR)
library(reshape2)
library(ggpubr)
library(viridis)
library(Matrix)
library(BiocParallel)
library(bluster)
library(celldex)
library(SingleR)
library(grid)
library(MetBrewer)
library(RColorBrewer)

# --- Loading data -------------------------------------------------------------
sce.hvg2 = readRDS("results/dim_reduc_20.rds")

################################################################################
#   SNN Clustering                                                             #
################################################################################-

## Elaborate SNN graph evaluating 20 neighbors
g = buildSNNGraph(sce.hvg2, use.dimred="PCA", type="number", k=20)

## Look for communities using random walks.
cluster = igraph::cluster_walktrap(g)$membership
sce.hvg2$cluster = factor(cluster)
table(sce.hvg$cluster)

# --- Plots --------------------------------------------------------------------
mycol = colorRampPalette(brewer.pal(name="Set3", n = 12))(20)

## PCA by cluster
plotReducedDim(sce.hvg2, dimred="PCA",colour_by="cluster", text_by="cluster", point_alpha = 0.7)+
  xlab("PC 1")+
  ylab("PC 2")+
  scale_color_manual(values = mycol, name = "Cluster")

## tSNE by cluster
plotReducedDim(sce.hvg2, dimred="TSNE_seed_100_per_40", colour_by="cluster", text_by="cluster", point_size=0.7)+
  xlab("tSNE 1")+
  ylab("tSNE 2")+
  theme_bw(base_size = 18) +
  theme(legend.position = "none") +
  scale_color_manual(values = mycol, name = "Cluster")

## UMAP by cluster
plotReducedDim(sce.hvg, dimred="UMAP_seed_100_n50", colour_by="cluster", text_by="cluster")+
  xlab("UMAP 1")+
  ylab("UMAP 2")+
  scale_color_manual(values = mycol, name = "Cluster")


################################################################################
#   Clustering evaluation                                                      #
################################################################################-

# --- Ratio of observed to expected edge: how well separated the clusters are --
ratio = pairwiseModularity(g, cluster, as.ratio=TRUE) 

## Plot
jpeg(file="6.Clustering_20/Pheatmap_ratio.jpeg", width=10, height=6, units="in", res=300)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE, col=viridis::inferno(100))
  setHook("grid.newpage", NULL, "replace")
  grid.text("Groups", y=-0.02, gp=gpar(fontsize=10))
  grid.text("Groups", x=-0.02, rot=90, gp=gpar(fontsize=10))
dev.off()


# --- Bootstrap to compute the co-assignment probability -----------------------
pcs = reducedDim(sce.hvg, "PCA")

set.seed(100)
cluster_function = function(x) {
  g = makeSNNGraph(x, type="number", k=20)
  igraph::cluster_walktrap(g)$membership
}

ass.prob = bootstrapStability(pcs, FUN=cluster_function, clusters=sce.hvg$cluster)

## Plot
jpeg(file="6.Clustering_20/Pheatmap_bootstrap.jpeg", width=10, height=6, units="in", res=300)
  setHook("grid.newpage", function() pushViewport(viewport(x=1,y=1,width=0.95, height=0.95, name="vp", just=c("right","top"))), action="prepend")
  pheatmap(ass.prob, cluster_cols=FALSE, cluster_rows=FALSE, col=viridis::inferno(100))
  setHook("grid.newpage", NULL, "replace")
  grid.text("Groups", y=-0.02, gp=gpar(fontsize=10))
  grid.text("Groups", x=-0.02, rot=90, gp=gpar(fontsize=10))
dev.off()


# --- Purity -------------------------------------------------------------------
purity <- neighborPurity(reducedDim(sce.hvg, "PCA"), clusters=sce.hvg$cluster, k=20)

pure.data <- as.data.frame(purity)
pure.data$maximum <- factor(pure.data$maximum)
pure.data$cluster <- factor(cluster)

jpeg(file="6.Clustering_20/purity.jpeg", width=10, height=6, units="in", res=300)
  ggplot(pure.data, aes(x=cluster, y=purity, colour=maximum)) +
    ggbeeswarm::geom_quasirandom(method="smiley") + 
    xlab("Groups") + 
    ylab("Purity")
dev.off()

################################################################################
#   Save results                                                               #
################################################################################

saveRDS(sce.hvg,"results/clustering_20_k20.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------




