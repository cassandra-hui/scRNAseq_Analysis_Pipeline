# PI:   
# Date:  
# Title: 
# Alingment Type: 

# Organism: 
# Conditions:
# Rep: 
# Design:


#############

library(Seurat)
library(dplyr)
library(ggplot2)
library(scDblFinder)
library(BiocSingular)
library(scater)
library(SingleCellExperiment)


source("/Users/cassandrahui/Documents/Projects/Scripts/scRNAseq/functions/functions.R")

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
fig_path <- "../05_figures/"
###################


dat <- readRDS("02_objects/01_seurat_obj_before_QC.rds")

slotNames(dat)
head(dat[[]])
dat


hist(dat$nCount_RNA, breaks = 400)
hist(dat$nCount_RNA, breaks = 600, xlim = range(1, 12000))  # to see lower end better

# Decide on cut off, halfway up the left side from lowest to highest

hist(dat$nFeature_RNA, breaks = 100)
hist(dat$nFeature_RNA, breaks = 600, xlim = range(1, 3000))  # to see lower end better

# Decide on cut off, halfway up the left side from lowest to highest
hist(dat$nFeature_RNA, breaks = 600, xlim = range(1, 1500))  # to see lower end better
abline(v = 40, col = "red", lwd = 2)

# Using 40

# Create a scatter plot with horizontal line at y = 7000
Seurat::FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_hline(yintercept = 7000, color = "red", linetype = "dashed")

# Using 7000



# Add in mitochondrial genes to see deteriorated cells
# Can skip if working with single nurclei
dat <- Seurat::PercentageFeatureSet(dat, 
                                    pattern = "^mt-", # mouse
                                    #pattern = "^MT-", # human
                                    col.name = "percent.mt")


plot <- Seurat::VlnPlot(dat, features = c("percent.mito",
                                  "nCount_RNA",
                                  "nFeature_RNA"))
plot
SaveFigure(plot, "01_vln_QC", width = 18, height = 12)


Seurat::FeatureScatter(dat, 
                       feature1 = "percent.mt", 
                       feature2 = "nFeature_RNA")


# Number of cells before QC
table(dat$orig.ident)

# sample_1 sample_2 sample_3 sample_4 
# 752     1227     2110      622   

# Trim low quality data and outliers on high end
dat <- subset(dat, 
              subset = nFeature_RNA < 7000 & percent.mt < 25 & nFeature_RNA > 40)


table(dat$orig.ident)
# Number of cells after filtering

# sample_1 sample_2 sample_3 sample_4 
# 717     1148     2070      591 

plot <- RidgePlot(dat, features=c("nFeature_RNA","nCount_RNA", "percent.mt"), ncol = 2)
plot
SaveFigure(plot, "02_ridge_QC_persample_postfilter", width = 12, height = 10)


### Remove Doublets


# Convert Seurat object to SingleCellExperiment with metadata
SCE <- SingleCellExperiment(
  assays = list(counts = GetAssayData(dat, layer = "counts")),
  colData = DataFrame(dat@meta.data),
  rowData = DataFrame(gene_names = rownames(dat))
)


# Compute Doublets
set.seed(42)
SCE$doublet_score <- computeDoubletDensity(SCE)

# Visualize
ggcells(SCE, aes(x = sample, y = doublet_score)) +
  geom_violin(fill = 'blue')



# Arbitrary threshold
doublet_score_threshold <- 2
SCE$probable_doublet <- (SCE$doublet_score > doublet_score_threshold)

plotColData(SCE, "doublet_score", colour_by = "probable_doublet")
sum(SCE$probable_doublet) #84


SCE_qc <- SCE[, !SCE$probable_doublet]
dim(SCE_qc)
# 25923  4442 

# 4442 cells now

dim(SCE) # Number of cells before
# 25923  4526



# Subset the original Seurat object
cells_to_keep <- colnames(SCE_qc)
dat_qc <- subset(dat, cells = cells_to_keep)
dat_qc

table(dat_qc$orig.ident)
# Number of cells after doublet removal

# sample_1 sample_2 sample_3 sample_4 
# 711     1131     2039      561 

saveRDS(dat_qc, "02_objects/Seurat_QC.rds")

