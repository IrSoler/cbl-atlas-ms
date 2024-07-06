# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Manuscript: Single cell landscape of sex differences in the progression of 
# multiple sclerosis
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Loading libraries --------------------------------------------------------

library(scran)
library(scater)
library(SingleCellExperiment)
library(MAST)
library(data.table)
library(edgeR)
library(dplyr)

# --- Local functions ----------------------------------------------------------

## Calculate log2(TPM + 1) counts
tpm <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

## Extract results from lrtest
extract_results = function(zz) {
  
  results = data.frame(rownames(zz))
  colnames(results) = "gene_id"
  
  results[,"lambda"] = data.frame(zz[,3,1]) #Statistic
  results[,"p.value"] = data.frame(zz[,3,3]) #p.value hurdle
  results[,"p.adjusted"] = p.adjust(results[,"p.value"], "BH") #p.adjusted BH
  
  return(results)
}

## logFC calculation
LFC.MAST = function (zz, control, case) {
  
  # Variables
  coefnames = colnames(coef(zz, 'D'))
  
  # Constrast 0
  contrast0 = setNames(rep(0, length(coefnames)), coefnames)
  contrast0[control] = 1
  
  # Contrast 1
  contrast1 = diag(length(coefnames))
  rownames(contrast1) = colnames(contrast1)<-coefnames
  
  # Calculate logFC and select test of interest
  LFC = MAST::logFC(zz, contrast0, contrast1)
  
  results = data.frame(rownames(LFC[['logFC']]))
  colnames(results) = "gene_id"
  results$logFC = as.data.frame(LFC[['logFC']])[,case]

  return(results)
}

# --- Load data ----------------------------------------------------------------

## SCE object
sce = readRDS("results/sce_norm.rds")
sce.hvg = readRDS("results/cell_type_20.rds")
sce$cell.type = sce.hvg$cell.type

################################################################################
#   zlm                                                                        #
################################################################################

# --- Normalization ------------------------------------------------------------

## Transcript length
## Data download from Biomart Ensembl Genes 106 -- Human Genes (GRCh38.p13) 05-10-2022
transcript = read.delim(file = "transcript_length.txt", header = TRUE, sep=",")

p = which(transcript$Gene.stable.ID %in% rowData(sce)[,"ENSEMBL"]) 
transcript = transcript[p,]

p2 = which(rowData(sce)[,"ENSEMBL"] %in% transcript$Gene.stable.ID)
sce = sce[p2,]

ensembl.id = rowData(sce)[,"ENSEMBL"]
transcript = transcript[match(ensembl.id, transcript$Gene.stable.ID),]
all.equal(rowData(sce)[,"ENSEMBL"], transcript$Gene.stable.ID)

## Add length data to SCE object
rowData(sce)[,"gene.start"] = transcript$Gene.start..bp.
rowData(sce)[,"gene.end"] = transcript$Gene.end..bp.
rowData(sce)[,"gene.length"] = transcript$Gene.end..bp. - transcript$Gene.start..bp. + 1

## Calculate tpm and log2tpm assays
assay(sce, "tpm") = tpm(counts = assay(sce, "counts"), len = rowData(sce)[,"gene.length"])
assay(sce, "log2tpm") = log2(assay(sce, "tpm") + 1)


# --- Subset by cell type ------------------------------------------------------

## Brain cell types
sce_astro = sce[,sce$cell.type == "Astrocytes"]
sce_microglia = sce[,colData(sce)$cell.type == "Microglia"]
sce_neurons = sce[,colData(sce)$cell.type == "Neurons"]
sce_oligo = sce[,colData(sce)$cell.type == "Oligodendrocytes"]
sce_oligoPC = sce[,colData(sce)$cell.type == "Oligodendrocyte PC"]

## Blood cell types
# sce_b_cells = sce[,sce$cell.type == "B cells"]
# sce_cd4 = sce[,colData(sce)$cell.type == "CD4+ T cells"]
# sce_cd8 = sce[,colData(sce)$cell.type == "CD8+ T cells"]
# sce_nk = sce[,colData(sce)$cell.type == "NK cells"]
# sce_dc = sce[,colData(sce)$cell.type == "Dendritic cells"]
# sce_mono = sce[,colData(sce)$cell.type == "Monocytes"]


# --- Hurdle model -------------------------------------------------------------

## Brain cell types
cell_objects = list(sce_astro, sce_microglia, sce_neurons, sce_oligo, sce_oligoPC)
cell_names = c("astro", "microglia", "neurons", "oligo", "oligoPC")

## Blood cell types
# cell_objects = list(sce_b_cells, sce_cd4, sce_cd8, sce_nk, sce_dc, sce_mono)
# cell_names = c("b_cells", "cd4", "cd8", "nk", "dc", "mono")

for (i in 1:length(cell_names)) {
  
  ## Scaling number of genes
  cdr2_cells = colSums(assay(cell_objects[[i]], "log2tpm")>0)
  colData(cell_objects[[i]])$cngeneson = scale(cdr2_cells)
  
  ## Construct sca object
  sce.assay_cells = SceToSingleCellAssay(cell_objects[[i]])
  rowData(sce.assay_cells)$primerid = rownames(sce.assay_cells)
  
  ## Model brain
  zlmCond_cells = zlm(~0+group_sex+cngeneson+sample+stage+age+region+Capbatch+Seqbatch+CellCycle, sce.assay_cells, exprs_values = 'log2tpm', method = "bayesglm", ebayes = TRUE, parallel = TRUE)
  
  ## Model blood
  # zlmCond_cells = zlm(~0+group_sex+cngeneson+donor+batch_pair+age_sampling+prev_treatments+CellCycle, sce.assay_cells, exprs_values = "log2tpm", method = "bayesglm", ebayes = TRUE, parallel = TRUE)
  saveRDS(object = zlmCond_cells, file = paste0("zlm_", cell_names[i], ".rds"))
}


################################################################################
#   DGE test                                                                   #
################################################################################

## Brain cell types
cell_objects = list(sce_astro, sce_microglia, sce_neurons, sce_oligo, sce_oligoPC)
cell_names = c("astro", "microglia", "neurons", "oligo", "oligoPC")

## Blood cell types
# cell_objects = list(sce_b_cells, sce_cd4, sce_cd8, sce_nk, sce_dc, sce_mono)
# cell_names = c("b_cells", "cd4", "cd8", "nk", "dc", "mono")

for (i in 1:length(cell_objects)) {
  
  ## Load zlm object
  zlmCond_cells = readRDS(paste0("zlm_", cell_names[i], ".rds"))
  
  
  ## Comparison Impact of Disease in Females (IDF)
  zlmCond_cellsF = lrTest(zlmCond_cells, Hypothesis("group_sexMS_female - group_sexControl_female"))
  results = extract_results(zz = zlmCond_cellsF)
  
  LFC.mF = LFC.MAST(zz = zlmCond_cells, case = "group_sexMS_female", control = "group_sexControl_female")
  results = merge(results, LFC.mF, by ="gene_id")
  
  write.table(results, file=paste0("all_", cell_names[i], "_F.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  
  
  ## Comparison Impact of Disease in Males (IDM)
  zlmCond_cellsM = lrTest(zlmCond_cells, Hypothesis("group_sexMS_male - group_sexControl_male"))
  results = extract_results(zz = zlmCond_cellsM)
  
  LFC.mM = LFC.MAST(zz = zlmCond_cells, control = "group_sexControl_male", case = "group_sexMS_male")
  results = merge(results, LFC.mM, by ="gene_id")
  
  write.table(results, file=paste0("all_", cell_names[i], "_M.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  
  
  ## Comparison Sex Differential Impact of Disease (SDID)
  
  zlmCond_cellsFM = lrTest(zlmCond_cells, Hypothesis("(group_sexPPMS_female - group_sexHI3_female) - (group_sexPPMS_male - group_sexHI3_male)"))
  results = extract_results(zz = zlmCond_cellsFM)
  
  ### SDID logFC
  coefnames = colnames(coef(zlmCond_cells, 'D'))

  contrast0 = setNames(rep(0, length(coefnames)), coefnames)
  contrast0["group_sexControl_female"] = 1
  contrast0["group_sexMS_male"] = 1

  contrast1 = diag(length(coefnames))
  rownames(contrast1) = colnames(contrast1)<-coefnames
  contrast1["group_sexControl_male", "group_sexMS_female"] = 1
  
  LFC = MAST::logFC(zlmCond_cells, contrast0, contrast1)
  p = data.frame(rownames(LFC[['logFC']]))
  colnames(p) = "gene_id"
  p$logFC = as.data.frame(LFC[['logFC']])[,"group_sexMS_female"]
  results = merge(results, p, by ="gene_id")
  
  ### Save results
  write.table(results, file=paste0("all_", cell_names[i], "_FM.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
  
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
