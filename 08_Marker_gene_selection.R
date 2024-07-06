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
library(MAST)
library(tibble)
library(dplyr)
library(viridis)
library(BiocParallel)

# --- Loading data -------------------------------------------------------------

sce = readRDS("./sce_norm.rds")
sce.hvg = readRDS("./cell_type_20.rds")

sce$cell.type = sce.hvg$cell.type
sce$cluster = sce.hvg$cluster

################################################################################
#   Identification of marker genes by cluster                                  #
################################################################################

# --- Identify gene markers ----------------------------------------------------
markers_genes = findMarkers(x = sce,
                            groups = sce$cluster,
                            test.type="wilcox",
                            pval.type = "all",
                            direction = "up")

# --- Explore gene markers -----------------------------------------------------

## Data frame
data.show = lapply(names(markers_genes), function(x) { temp = markers_genes[[x]][, 1:3] ; temp$gene = rownames(markers_genes[[x]]) ; temp$url = ifelse(temp$gene=="NA",NA,paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", temp$gene,"'>",temp$gene,"</a>",sep="")) ; return(temp) } )
names(data.show) = names(markers_genes)

## Show results (execute it for each cluster)
results = as.data.frame(data.show[["1"]][,c(1:2, 5)])
sig = results$FDR < 0.05 
sig = results[sig,]
nsig = nrow(sig)
DT::datatable(results, caption = paste0(nsig, " marker genes"), escape = FALSE)

################################################################################
#   Identification of marker genes by cell type                                #
################################################################################

# --- Identify gene markers ----------------------------------------------------
markers_genes = findMarkers(x = sce,
                            groups = sce$cell.type,
                            test.type="wilcox",
                            pval.type = "all",
                            direction = "up")

# --- Explore gene markers -----------------------------------------------------

## Data frame
data.show = lapply(names(markers_genes), function(x) { temp = markers_genes[[x]][, 1:3] ; temp$gene = rownames(markers_genes[[x]]) ; temp$url = ifelse(temp$gene=="NA",NA,paste("<a href='https://www.genecards.org/cgi-bin/carddisp.pl?gene=", temp$gene,"'>",temp$gene,"</a>",sep="")) ; return(temp) } )
names(data.show) = names(markers_genes)

## Show results (execute it for each cell type)
results = as.data.frame(data.show[["CD4+ T-cells"]][,c(1:2, 5)])
sig = results$FDR < 0.05 
sig = results[sig,]
nsig = nrow(sig)
DT::datatable(results, caption = paste0(nsig, " marker genes"), escape = FALSE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------