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
library(CellChat)
library(patchwork)
library(DT)

# --- Loading data -------------------------------------------------------------
sce = readRDS("results/sce_norm.rds")
sce.hvg = readRDS("results/cell_type_20.rds")
sce$cell.type = sce.hvg$cell.type
rm(sce.hvg)

################################################################################
#   CellChat pipeline                                                          #
################################################################################

# --- Build CellChat object ----------------------------------------------------

## Control female
sce_c_F = sce[,sce$group_sex=="Control_female"]
cellchat_c_F = createCellChat(object = sce_c_F, group.by = "cell.type", assay = "logcounts")

## MS female
sce_ms_F = sce[,sce$group_sex=="MS_female"]
cellchat_ms_F = createCellChat(object = sce_ms_F, group.by = "cell.type", assay = "logcounts")

## Control male
sce_c_M = sce[,sce$group_sex=="Control_male"]
cellchat_c_M = createCellChat(object = sce_c_M, group.by = "cell.type", assay = "logcounts")

## MS male
sce_ms_M = sce[,sce$group_sex=="MS_male"]
cellchat_ms_M = createCellChat(object = sce_ms_M, group.by = "cell.type", assay = "logcounts")

## Join cellchat objects
object.list <- list(CF = cellchat_c_F, MSF = cellchat_ms_F, CM = cellchat_c_M, MSM = cellchat_ms_M)

# --- Add CellChat database ----------------------------------------------------

## Load database
CellChatDB = CellChatDB.human
CellChatDB.use = CellChatDB

## Add to each group
for (i in 1:4) {
  
  element = object.list[[i]]
  element@DB = CellChatDB.use
  object.list[[i]] = element
}

# --- Identify significant genes and interactions ------------------------------

for (i in 1:4) {
  
  cellchat = object.list[[i]]
  
  ## Subset the expression data of signaling genes for saving computation cost.
  cellchat <- subsetData(cellchat)

  ## identifyOverExpressedGenes has been modify locally to adjust p.value
  ## by the method BH
  cellchat <- identifyOverExpressedGenesLocal(cellchat, data.use = cellchat@data.signaling) # Select all genes

  ## Parsing results
  cellchat@var.features[["features.info"]] = cellchat@var.features[["features.info"]][cellchat@var.features[["features.info"]][, "pvalues.adj"] < 0.05,]
  cellchat@var.features[["features"]] = cellchat@var.features[["features.info"]][,"features"]

  ## Idenfity overexpressed interactions
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  ## Update object list
  object.list[[i]] = cellchat
}

# --- Project gene expression data onto PPI ------------------------------------

for (i in 1:4) {
  cellchat = object.list[[i]]
  cellchat <- projectData(cellchat, PPI.human)
  object.list[[i]] = cellchat
}

# --- Compute the communication probability and infer cellular communication network ---
for (i in 1:4) {
  
  cellchat = object.list[[i]]
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  object.list[[i]] = cellchat

}


# --- Extract the inferred cellular communication network ----------------------

## Control female
df.net.CF <- subsetCommunication(object.list[['CF']], thresh = 1.1)
df.net.CF$p.adjusted = stats::p.adjust(df.net.CF$pval, method = "BH")
write.csv(x = df.net.CF, file = "all_interactions_Control_female.csv")

### MS female
df.net.MSF <- subsetCommunication(object.list[['MSF']], thresh = 1.1)
df.net.MSF$p.adjusted = stats::p.adjust(df.net.MSF$pval, method = "BH")
write.csv(x = df.net.MSF, file = "all_interactions_MS_female.csv")

### Control male
df.net.CM <- subsetCommunication(object.list[['CM']], thresh = 1.1)
df.net.CM$p.adjusted = stats::p.adjust(df.net.CM$pval, method = "bonferroni")
write.csv(x = df.net.CM, file = "all_interactions_Control_male.csv")

### MS male
df.net.MSM <- subsetCommunication(object.list[['MSM']], thresh = 1.1)
df.net.MSM$p.adjusted = stats::p.adjust(df.net.MSM$pval, method = "bonferroni")
write.csv(x = df.net.MSM, file = "all_interactions_MS_male.csv")



# --- Infer the cell-cell communication at a signaling pathway level -----------

## Compute probability
for (i in 1:4) {
  
  cellchat = object.list[[i]]
  cellchat <- computeCommunProbPathway(cellchat)
  object.list[[i]] = cellchat
}

## Aggregate net
for (i in 1:4) {
  
  cellchat = object.list[[i]]
  cellchat <- aggregateNet(cellchat)
  object.list[[i]] = cellchat
}

# --- Save results -------------------------------------------------------------

cellchat <- mergeCellChat(object.list, add.names = names(object.list))
saveRDS(object.list, file = "ccc_list.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------