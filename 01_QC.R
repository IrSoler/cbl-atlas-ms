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
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(Matrix)
library(cowplot)
library(viridis)
library(BiocParallel)
library(scDblFinder)

# --- Loading data -------------------------------------------------------------

# COUNT MATRIX
notfiltered.counts = readMM("matrix.mtx")
class(notfiltered.counts)
dim(notfiltered.counts)

# ROW DATA
notfiltered.rowMetadata = read.delim("genes.tsv", header = F)
class(notfiltered.rowMetadata)
dim(notfiltered.rowMetadata)

rownames(notfiltered.counts) = notfiltered.rowMetadata[,1]

# COLUMN DATA
notfiltered.colData = read.delim("barcodes.tsv", header = F)
class(notfiltered.colData)
dim(notfiltered.colData)

colnames(notfiltered.counts) = notfiltered.colData[,1]

# METADATA
notfiltered.colMetadata = read.delim("meta.txt", sep = "\t", header = TRUE)
class (notfiltered.colMetadata)
dim(notfiltered.colMetadata)

# OBJECT SingleCellExperiment 
sce = SingleCellExperiment(assays = list(counts = notfiltered.counts))
rowData(sce)$SYMBOL = notfiltered.rowMetadata[,1]
rowData(sce)$ENSEMBL = notfiltered.rowMetadata[,1]

# --- Preparing metadata -------------------------------------------------------

#Metadata
cell_metadata = notfiltered.colMetadata

# Check coldata names
all.equal(cell_metadata[,1], colnames(sce))

rownames(cell_metadata) = cell_metadata[,1] #Assign cell names to row names
cell_metadata = cell_metadata[,-1] #Remove first column if cell names is the first column

# Make new covariable group_sex
group_sex <- paste(cell_metadata$diagnosis, cell_metadata$sex, sep = '_')

#Assign metadata to colData of sce object
colData(sce) = DataFrame(cell_metadata, group_sex) 
sce$group_sex = factor(sce$group_sex, levels = c("Control_female", "MS_female", "Control_male", "MS_male"))

################################################################################
#   QC metrics                                                                 #
################################################################################

# METRICS BY CELL

## Doublets identification
set.seed(123)
sce <- scDblFinder(sce, BPPARAM = MulticoreParam(60), samples = sce$Capbatch)

## Library size. Number of genes counts per cell
sce@colData$nCount = Matrix::colSums(counts(sce)) 

## Number of genes per cell
sce@colData$nFeatures = Matrix::colSums(counts(sce)>0) 

## Mitochondrial genes
mito_genes = rownames(sce)[grep("^MT-",rownames(sce))] 
sce@colData$percent_mito = Matrix::colSums(counts(sce)[mito_genes, ]) / sce@colData$nCount 

## Ribosomal genes
ribo_genes = rownames(sce)[grep("^RP[SL]",rownames(sce))] 
sce@colData$percent_ribo = Matrix::colSums(counts(sce)[ribo_genes, ]) / sce@colData$nCount 

# METRICS BY GENES

## nCells for each gene
rowData(sce)$nCells <- Matrix::rowSums(counts(sce)>0)

## pseudogenes are identify with the last consonant 'P'
pseudogenes = rownames(sce)[grep("P$",rownames(sce))]
pseudogenes1 = rownames(sce)[grep("P[0-9]+$",rownames(sce))]

## tRNA
trna = rownames(sce)[grep("^TR.-",rownames(sce))]

## Small nuclear RNAs
snrna = rownames(sce)[grep("^RNU",rownames(sce))]

## Small nucleolar RNAs
snorna = rownames(sce)[grep("^SNORD",rownames(sce))] # SNORD# for “small nucleolar RNA, C/D box” genes
snorna1 = rownames(sce)[grep("^SNORA",rownames(sce))] # SNORA# for “small nucleolar RNA, H/ACA box” genes
snorna2 = rownames(sce)[grep("^SCARNA",rownames(sce))] # SCARNA# for “small Cajal body‐specific RNA” genes

## Ribosomal RNAs
rrna = rownames(sce)[grep("^RNA[0-9]",rownames(sce))]
rprna = rownames(sce)[grep("^RP[0-9]",rownames(sce))]

#Create metadata dataframe
metadata_2 = as.data.frame(sce@colData)

################################################################################
#   QC plots before filtering                                                  #
################################################################################

# Colors
fu = "#D46C40"
fd = "#FFC7AF"
mu = "#497EC2"
md = "#ADCBF8"
su = "#906BBA"
sd = "#CCB8E3"


# Density plot number of genes
df = data.frame(nFeatures = sce$nFeatures, nCount = sce$nCount)
plot(density(sce$nFeatures), xlab="Number of expressed genes by cell", ylab="Density\n", main = "", las = 1, ylim = c(0, 0.0008))


# Number of genes (explore other covariables)
plotColData(sce,y = "nFeatures",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("Number of expressed genes")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.position = "top") +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Density plot for number of counts 
plot(density(sce$nCount), xlab="Library size", ylab="Density\n",  main = "", las = 1, ylim = c(0, 0.0004))


# Number of counts (explore other covariables)
plotColData(sce,y = "nCount",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("Library size")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.position = "top") +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Percentage of mitochondrial genes (explore other covariables)
plotColData(sce,y = "percent_mito",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("Mitochondrial \n gene ratio")+
  ggtitle("Before quality control filtering")+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.position = "bottom", legend.justification = "left") +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Number of counts vs number of genes for each cell (explore other covariables)
plotColData(sce,x = "nCount",y = "nFeatures",colour_by = "sample", point_alpha=0.3)+
  xlab("Library size")+
  ylab("Number of expressed genes")+
  ggtitle("Before quality control filtering")+
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom", legend.justification = "left", legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_colour_viridis_d(option = "C",  name = "Sample ID", direction = 1, alpha = 1)


# Mitochondrial percentage vs number of genes for each cell (explore other covariables)
plotColData(sce,x = "percent_mito",y = "nFeatures",colour_by = "group_sex")+
  xlab("Mithocondrial gene ratio")+
  ylab("Number of genes")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Plot number of cells by patient (explore other covariables)
metadata_2 %>%
  ggplot(aes(x=sample, fill=group_sex)) +
  geom_bar() +
  xlab("Sample ID")+
  ylab("Number of cells")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top") +
  scale_fill_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


################################################################################
#   QC filtering                                                               #
################################################################################

# UCSC-MS: filtering cutoff in Material and Methods
# GSE144744 Data was pre-filtered: check thresholds

# FILTERING BY CELL

discard_count = sce@colData$nCount <= 1000 # from the article
table(discard_count)

discard_ngenes = sce@colData$nFeatures <= 500 # from the article
table(discard_ngenes)

discard_mito = isOutlier(sce$percent_mito,type= "higher")
table(discard_mito)

discard_ribo = isOutlier(sce$percent_ribo, type="higher")
table(discard_ribo)

discard_doublet = sce@colData$scDblFinder.class == "doublet"
table(discard_doublet)

discard = discard_count | discard_ngenes | discard_mito | discard_ribo | discard_doublet

# FILTERING BY GENE
selected_features = rowData(sce)$nCells > 3 # from the article

a <- c(pseudogenes, pseudogenes1, trna, snrna, snorna, snorna1, snorna2, rrna, rprna)
not_selected_pseudogenes = rownames(sce) %in% a
nkeep = !not_selected_pseudogenes & selected_features

# COMPLETE FILTERING
sce.filt = sce[nkeep, !discard]

################################################################################
#   QC metrics after filtering                                                 #
################################################################################

# Library size
sce.filt@colData$nCountAF = Matrix::colSums(counts(sce.filt)) 

# Number of genes per cell
sce.filt@colData$nFeaturesAF = Matrix::colSums(counts(sce.filt)>0) 

# Mitochondrial genes
mito_genesAF = rownames(sce.filt)[grep("^MT-",rownames(sce.filt))] 
sce.filt@colData$percent_mitoAF = Matrix::colSums(counts(sce.filt)[mito_genesAF, ]) / sce.filt@colData$nCountAF 

# Ribosomal genes
ribo_genesAF = rownames(sce.filt)[grep("^RP[SL]",rownames(sce.filt))] 
sce.filt@colData$percent_riboAF = Matrix::colSums(counts(sce.filt)[ribo_genesAF, ]) / sce.filt@colData$nCountAF 

# Relative expression of each gene per cell
plotHighestExprs(sce.filt, exprs_values = "counts")

################################################################################
#   QC plot after filtering                                                    #
################################################################################

# Create metadata dataframe
metadata_3 = as.data.frame(sce.filt@colData)


# Density plot number of genes
df = data.frame(nFeaturesAF = sce.filt$nFeaturesAF, nCountAF = sce.filt$nCountAF)
summary(sce.filt$nFeaturesAF)
plot(density(sce.filt$nFeaturesAF), xlab="Number of expressed genes by cell", ylab="Density\n", main = "", las = 1, ylim = c(0, 0.0010))


# Number of genes for each covariable (explore other covariables)
plotColData(sce.filt,y = "nFeaturesAF",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("Number of expressed genes")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.position = "top") +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Density plot for number of counts
plot(density(sce.filt$nCountAF), xlab="Library size", ylab="Density\n", main = "", las = 1, ylim = c(0, 0.0005))


# Number of counts (explore other covariables)
plotColData(sce.filt,y = "nCountAF",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("Library size")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.position = "top") +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Percentage of mitochondrial genes (explore other covariables)
plotColData(sce.filt,y = "percent_mito",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("Mitochondrial \n gene ratio")+
  ggtitle("After quality control filtering")+
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),  legend.position = "bottom", legend.justification = "left") +
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Number of counts vs number of genes (explore other covariables)
plotColData(sce.filt,x = "nCountAF",y = "nFeaturesAF",colour_by = "sample", point_alpha=0.3)+
  xlab("Library size")+
  ylab("Number of expressed genes")+
  ggtitle("After quality control filtering")+
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom", legend.justification = "left", legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  scale_colour_viridis_d(option = "C",  name = "Sample ID", direction = 1, alpha = 1)


# Mitochondrial percentage vs number of genes (explore other covariables)
plotColData(sce.filt,x = "percent_mitoAF",y = "nFeaturesAF",colour_by = "group_sex")+
  xlab("Gene mithocondrial ratio")+
  ylab("Number of genes")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# Plot number of cells by patient (explore other covariables)
metadata_3 %>%
  ggplot(aes(x=sample, fill=group_sex)) +
  geom_bar() +
  xlab("Sample ID")+
  ylab("Number of cells")+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top") +
  scale_fill_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


################################################################################
#   Cell cycle                                                                 #
################################################################################

# Obtain human marker genes
hm.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

# Calculate cell cycle phase scores and state
assignments = cyclone(sce.filt, hm.pairs, gene.names=rowData(sce.filt)$ENSEMBL, verbose=TRUE, BPPARAM = MulticoreParam(60))
sce.filt$CellCycle = assignments$phases
sce.filt$G1.score = assignments$scores$G1
sce.filt$G2M.score = assignments$scores$G2M
sce.filt$S.score = assignments$scores$S

# --- Plot cell cycle scores ---------------------------------------------------

## G2M vs G1 (explore other covariables)
plotColData(sce.filt,y = "G2M.score",x = "G1.score",colour_by = "group_sex")+
  xlab("Puntuación para fase G1")+
  ylab("Puntuación para las fases G2-M")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# G2M (explore other covariables)
plotColData(sce.filt,y = "G2M.score",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("G2M score")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# G1 (explore other covariables)
plotColData(sce.filt,y = "G1.score",x = "sample",colour_by = "group_sex")+
  xlab("Sample ID")+
  ylab("G1 score")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


# S (explore other covariables)
plotColData(sce.filt,y = "S.score",x = "sample",colour_by = "group_sex")+
  xlab("Sample")+
  ylab("S score")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "top")+
  scale_colour_manual(values = c(fd, fu, md, mu),  name = "Group", labels = c("Control female", "MS female", "Control male", "MS male"))


################################################################################
#   Save results                                                               #
################################################################################

saveRDS(sce.filt,"results/sce_qc.rds")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------