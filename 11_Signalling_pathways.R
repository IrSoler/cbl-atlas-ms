# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Manuscript: Single cell landscape of sex differences in the progression of 
# multiple sclerosis
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Loading libraries --------------------------------------------------------

library(hipathia)
library(SingleCellExperiment)
library(scuttle)
library(edgeR)
library(MAST)

# --- Loading data -------------------------------------------------------------
## SCE
sce.hvg = readRDS("/clinicfs/userhomes/isoler/TFM/PRJNA544731/results/cell_type_20.rds")
sce = readRDS("/clinicfs/userhomes/isoler/TFM/PRJNA544731/results/sce_norm.rds")
sce$cell.type = sce.hvg$cell.type
rm(sce.hvg)

## Pathways
pathways = load_pathways(species = "hsa")

# --- Local functions ----------------------------------------------------------
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
  colnames(results) = "pathID"
  results$logFC = as.data.frame(LFC[['logFC']])[,case]
  
  return(results)
}

################################################################################
#   Hipathia transformations                                                   #
################################################################################

## Translate data to ENTREZID
sce2 = as(sce, "SummarizedExperiment")
rownames(sce2) = rownames(sce)
rowData(sce2) = rowData(sce)
assay(sce2, "logcounts") = as.matrix(assay(sce, "logcounts"))

sce = sce2
rm(sce2)

sce = translate_data(sce, "hsa", sel_assay = 2)

## Scaled data
sce_scaled = normalize_data(sce, sel_assay = 1)

## Hipathia analysis
results = hipathia(sce_scaled, pathways, decompose = FALSE, verbose=TRUE)

## Save results
saveRDS(results, file="./results/hipathia.rds")

################################################################################
#   Differential activation analysis                                           #
################################################################################

results = readRDS("hipathia.rds")
path_vals = get_paths_data(results)
test = c("(group_sexMS_female - group_sexControl_female)", "(group_sexMS_male - group_sexControl_male)",
         "(group_sexMS_female - group_sexControl_female) - (group_sexMS_male - group_sexControl_male)")
x = 1

## Brain samples
names = c("astro", "micro", "neurons", "oligo", "OPC")
cells = c("Astrocytes", "Microglia", "Neurons", "Oligodendrocytes", "Oligodendrocyte PC")

## Blood samples
# names = c("b_cells", "cd4", "cd8", "nk", "dc", "mono")
# cells = c("B cells", "CD4+ T cells", "CD8+ T cells", "NK cells", "Dendritic cells", "Monocytes")


for (cell in cells) {
  
  ## Extract results for each cell type
  is.cell = which(sce.n$cell.type==cell)
  sce.cell = sce.n[ , is.cell]  
  cell_path = path_vals[, is.cell]

  ## Check if the order is established correctly
  print(all(colnames(sce.cell) == colnames(cell_path)))
  
  ## Build SCA object
  sce.cell_path = as(cell_path, "SingleCellExperiment")
  sca.cell = SceToSingleCellAssay(sce.cell_path)
  rowData(sca.cell)$primerid = rownames(cell_path)
  
  ## Model
  zlmCond = zlm(~0+group_sex+sample+Capbatch+Seqbatch, sca.cell)
  
  ## Statistical analysis
  for (i in c(1:3)) {
    
    ## lrTest
    zlmCond.test = lrTest(zlmCond, Hypothesis(test[i]))
    
    ## Extract results
    results = data.frame(rownames(zlmCond.test)) # pathID
    colnames(results) = "pathID"
    results[,"lambda"] = data.frame(zlmCond.test[,3,1]) # statistic
    results[,"p.value"] = data.frame(zlmCond.test[,3,3]) # p-value
    results[,"p.adjusted"] = p.adjust(results[,3], "BH") # p-adjusted by BH
    
    ## logFC
    if (i == 1) {
      contrast = "F"
      LFC.mF = LFC.MAST(zz = zlmCond, control = "group_sexControl_female", case = "group_sexMS_female")
      results = merge(results, LFC.mF, by ="pathID")
    }
    
    if (i == 2) {
      contrast = "M"
      LFC.mM = LFC.MAST(zz = zlmCond, control = "group_sexControl_male", case = "group_sexMS_male")
      results = merge(results, LFC.mM, by ="pathID")
      
    }
    
    if (i == 3) {
      contrast = "FM"
      
      # Variables
      coefnames = colnames(coef(zlmCond, 'D'))
      
      # Constrast 0
      contrast0 = setNames(rep(0, length(coefnames)), coefnames)
      contrast0["group_sexControl_female"] = 1
      contrast0["group_sexMS_male"] = 1
      
      # Contrast 1
      contrast1 = diag(length(coefnames))
      rownames(contrast1) = colnames(contrast1)<-coefnames
      contrast1["group_sexControl_male", "group_sexMS_female"] = 1
      
      # Calculate logFC and select test of interest
      LFC = MAST::logFC(zlmCond, contrast0, contrast1)
      
      p = data.frame(rownames(LFC[['logFC']]))
      colnames(p) = "pathID"
      p$logFC = as.data.frame(LFC[['logFC']])[,"group_sexMS_female"]
      
      results = merge(results, p, by ="pathID")

    }
    
    feat.name = rowData(sca.cell)
    feat.name = as.data.frame(feat.name[,c(1,2)])
    colnames(feat.name) = c("pathID", "path_name")
    
    ## Path names
    results = merge(results, feat.name, by ="pathID")
    
    ## Save results
    name = paste0( names[x], "_paths_", contrast, ".tsv")
    write.table(results[,c(1:4, 5, 6)], file = name, col.names=TRUE, row.names = FALSE, sep = "\t")
  }
  x = x + 1
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
