# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# Manuscript: Single cell landscape of sex differences in the progression of 
# multiple sclerosis
# Author: Irene Soler Sáez
# Computational Biomedicine Laboratory, CIPF (Valencia, Spain)
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# --- Loading libraries --------------------------------------------------------

library(genefilter)
library(GEOquery)
library(hgu133plus2.db)
library(topGO)
library(Hmisc)
library(limma)
library(dplyr)
library(GO.db)

# --- Loading data -------------------------------------------------------------

GO_def = suppressMessages(AnnotationDbi::select(GO.db,
                                                keys=(keys(GO.db)),
                                                columns=c("GOID","TERM","DEFINITION"),
                                                keytype="GOID"))
rownames(GO_def) = GO_def$GOID

# --- Local functions ----------------------------------------------------------

gene_filter <- function(dataset, cutoff =  0.05, fc = 0.5) {
  
  # Filtering genes
  sig = dataset[dataset$p.adjusted < cutoff,]
  down = rownames(sig[(sig$logFC < (-fc)),])
  up = rownames(sig[(sig$logFC > (fc)),])
  geneList = as.vector(dataset$p.adjusted)
  names(geneList) = row.names(dataset)
  
  # Selection of significant upregulated genes
  upList <- factor(as.integer(names(geneList) %in% up))
  names(upList) <- names(geneList)
  
  # Selection of significant downregulated genes
  downList <- factor(as.integer(names(geneList) %in% down))
  names(downList) <- names(geneList)

  return(list(upList, downList))
}

# --- General variables --------------------------------------------------------

cutoff = 0.05
fc = 0.5
org = "org.Hs.eg.db"
notation = "symbol"
met = "Rel"
alg = "personalizado"
cufoff2 = 0.8

################################################################################
#   Functional analysis                                                        #
################################################################################

## Brain cell types
cells = c("astro", "microglia", "neurons", "oligo", "oligoPC")

## Blood cell types
# cells = c("b_cells", "cd4", "cd8", "t_cells", "nk", "dc", "mono")

## Comparisons
comparisons = c("F", "M", "FM")

for (cell in cells) {
  
  for (comparison in comparisons) {
    
    # --- Load data ------------------------------------------------------------
    
    ## Read DGE results
    dataset = read.table(file = paste0("all_", cell, "_", comparison, ".tsv"),
                         sep = '\t', header = TRUE)
    rownames(dataset) = dataset$gene_id
    dataset = dataset[,-1]
    
    # Substitute NAs to 0s
    dataset = dataset %>%
      mutate(across(everything(), ~ ifelse(is.na(.), 0, .)))
    
    ## Select significant genes: p.adj < 0.05 and |logFC| > 0.5
    genes_sig_up = gene_filter(dataset = dataset)[[1]] # Genes upregulated
    genes_sig_down = gene_filter(dataset = dataset)[[2]]  #Genes downregulated
    
    # --- Functional profiling for upregulated genes ---------------------------
    if (sum(genes_sig_up == 1) > 0) {
      
      ## Build topGOdata object
      pre_up <- new("topGOdata", 
                    ontology = "BP", 
                    allGenes = genes_sig_up, 
                    geneSel = function(x)(x == 1), 
                    nodeSize = 5,
                    annot = annFUN.org,
                    mapping = org,
                    ID = notation)
      
      ## Perform the analysis
      weight01_up <- runTest(pre_up, algorithm = "weight01", statistic = "fisher", cutOff = pvalor)
      weight01_t_up = GenTable(pre_up, weight01 = weight01_up,  topNodes = length(score(weight01_up)))
      weight01_t_up$TERM = GO_def[weight01_t_up$GO.ID, "TERM"]
      weight01_t_up$p.value = sort(weight01_up@score) 
      
      ## Add gene info to the data.frame
      sigGenes <- sigGenes(pre_up)
      AnnoList <- lapply(weight01_t_up$"GO.ID", function(x) as.character(unlist(genesInTerm(object = pre_up, whichGO = x))))
      SigList <- lapply(AnnoList, function(x) intersect(x, sigGenes))
      weight01_t_up$"Genes" <- sapply(SigList, paste, collapse = ",")
      
      ## Save results
      write.table(weight01_t_up, file=paste0("all_weight01_", cell, "_", comparison, "_UP.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")

    }
    
    # --- Functional profiling for downregulated genes ---------------------------
    if (sum(genes_sig_down == 1) > 0) {
      
      ## Build topGOdata object
      pre_down <- new("topGOdata", 
                      ontology = "BP", 
                      allGenes = genes_sig_down, 
                      geneSel = function(x)(x == 1), 
                      nodeSize = 5,
                      annot = annFUN.org,
                      mapping = organismo,
                      ID = notacion)
      
      ## Perform the analysis
      weight01_down <- runTest(pre_down, algorithm = "weight01", statistic = "fisher", cutOff = pvalor)
      weight01_t_down = GenTable(pre_down, weight01 = weight01_down,  topNodes = length(score(weight01_down)))
      weight01_t_down$TERM = GO_def[weight01_t_down$GO.ID, "TERM"]
      weight01_t_down$p.value = sort(weight01_down@score)
      
      ## Add gene info to the data.frame
      sigGenes <- sigGenes(pre_down)
      AnnoList <- lapply(weight01_t_down$"GO.ID", function(x) as.character(unlist(genesInTerm(object = pre_down, whichGO = x))))
      SigList <- lapply(AnnoList, function(x) intersect(x, sigGenes))
      weight01_t_down$"Genes" <- sapply(SigList, paste, collapse = ",")
      
      ## Save results
      write.table(weight01_t_down, file=paste0("all_weight01_", cell, "_", comparison, "_DOWN.tsv"), col.names=TRUE, row.names=FALSE, sep="\t")
    }
  }  
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------