# PI:   
# Date:  
# Title: 
# Alingment Type: 

# Organism: 
# Conditions:
# Rep: 
# Design:


#########################################################################
#########################################################################
# libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(kableExtra)
library(biomaRt)
library(knitr)
library(limma)
library(topGO)
library(scMRMA)
library(DoubletFinder)
library(org.Mm.eg.db)
#rm(list = ls())

source('functions.R')

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
data_path <- "../01_data/DGE_filtered"
fig_path <- "../05_figures"
###################

# Read in Data
mat_path <- data_path
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

# Create object
dat <- CreateSeuratObject(mat, min.cells = 100,
                           names.field = 0, meta.data = cell_meta)



# Rename orig.ident by samples
dat@meta.data$orig.ident <- dat@meta.data$sample
Idents(dat) <- dat@meta.data$orig.ident


SaveObject(dat, "01_seurat_obj_before_QC")
#dat <- ReadObject("01_seurat_obj_before_QC")

slotNames(dat)
head(dat[[]])

set.seed(12345)
kable(do.call("cbind", tapply(dat$nFeature_RNA, 
                              Idents(dat),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of Genes/Cell by Sample") %>% kable_styling()


kable(do.call("cbind", tapply(dat$nCount_RNA, 
                              Idents(dat),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of UMI/Cell by Sample") %>% kable_styling()


# View Mitochondrial percentatge (MT- for Human)
dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^mt-")
plot <- VlnPlot(dat, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot
SaveFigure(plot, "01_vln_QC", width = 18, height = 12)


# View 
plot1 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 +plot2
SaveFigure((plot1 + plot2),"01_scatter_beforeQC", width = 13, height = 6, res = 200)

# Number of cells before QC
table(dat$orig.ident)

# sample_1 sample_2 sample_3 sample_4 
# 200      452      851      338     total: 1841



## Cut to see data better
dat <- subset(dat, subset = nFeature_RNA < 6000 & percent.mt < 40)

# Plot again
# Edit lines to determine best cut
plot <- FeatureScatter(
  dat, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5, shuffle = TRUE)  + geom_vline(xintercept = c(200,15000)) + geom_hline(yintercept = c(150, 3000))
plot
SaveFigure((plot),"02_QC_cuts", width = 7, height = 6, res = 200)
table(dat$orig.ident)

# Number of cells after extreme cuts
#
# sample_1 sample_2 sample_3 sample_4 
# 752     1227     2110      622     total: 4711


plot <- RidgePlot(dat, features=c("nFeature_RNA","nCount_RNA", "percent.mt"))
plot
SaveFigure(plot, "01_ridge_beforeQC_persample", width = 15, height = 10)
plot(sort(Matrix::rowSums(GetAssayData(dat) >= 3), decreasing = TRUE) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")

# Perform the filtering
dat <- subset(dat, subset = nCount_RNA > 200 & nCount_RNA < 15000 & nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 25)

table(dat$orig.ident)
# Number of cells after filtering

# sample_1 sample_2 sample_3 sample_4 
# 200      452      851      338     total: 1841

plot <- RidgePlot(dat, features=c("nFeature_RNA","nCount_RNA", "percent.mt"), ncol = 2)
plot
SaveFigure(plot, "02_ridge_QC_persample_postfilter", width = 12, height = 10)

####################################################
#Normalize
dat <- NormalizeData(dat, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Features
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dat), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
SaveFigure((plot1 + plot2), "03_var_features", width = 12, height = 6)

SaveObject(dat, "02_normalized_dat")
#dat <- ReadObject("02_normalized_dat")

####################################################
#Scale and Run PCA
dat <- ScaleData(dat)

dat <- RunPCA(dat, npcs = 100)
SaveObject(dat, "03_seurat_obj_after_PCA")
#dat <- ReadObject("03_seurat_obj_after_PCA")

####################################################
# Examine and visualize PCA results a few different ways
print(dat[["pca"]], dims = 1:5, nfeatures = 5)

plot <- DimPlot(object = dat, reduction = "pca")
plot
SaveFigure(plot, "03_PCA", width = 8, height = 7.5)
plot <- VizDimLoadings(dat, dims = 1:2, reduction = "pca")
plot
SaveFigure(plot, "03_viz_PCA_loadings", width = 10, height = 8)



# Image doesn't save as png unless fast = FALSE
plot <- DimHeatmap(dat, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
plot
SaveFigure(plot, "03_dim_heatmap1", width = 6, height = 6)


plot <- DimHeatmap(dat, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
plot
SaveFigure(plot, "03_dim_heatmap1_15", width = 12, height = 18)


# View to determine best number of PCs or Run Jackstraw 
plot <- ElbowPlot(dat,ndims = 50)
plot
SaveFigure(plot, "03_PC_elbow_plot", width = 8, height = 8)

####################################
# This will take a long time with large data

## Jackstraw
############

dat <- JackStraw(dat, dims = 100) # Check this is right
#dat <- JackStraw(dat, num.replicate = 100) 
dat <- ScoreJackStraw(dat, dims = 1:100)

plot <- JackStrawPlot(object = dat, dims = 1:100) + theme(legend.position="bottom")
plot
SaveFigure(plot, "03_Jack_Straw", width = 12, height = 8)

SaveObject(dat, "04_seurat_obj_after_Jackstraw")
#dat <- ReadObject("04_seurat_obj_after_Jackstraw")
############

# Select pcs
use.pcs <- 1:14
dat <- FindNeighbors(dat, dims = use.pcs)
# Create numtiple resolutions to view
dat <- FindClusters(dat, resolution = seq(0.25, 4, 0.5), verbose = FALSE)


dat <- RunTSNE(
  object = dat,
  reduction.use = "pca",
  dims = use.pcs,
  do.fast = TRUE)



plot <- DimPlot(object = dat,
                group.by=grep("res",colnames(dat@meta.data),
                              value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "tsne", label = T)
plot
SaveFigure(plot, "04_Resolutions", width = 14, height = 12)


# Select resolution to use
Idents(dat) <- "RNA_snn_res.0.25"
dat <- RunUMAP(dat, dims = use.pcs)
plot <- DimPlot(dat, reduction = "umap")
plot
SaveFigure(plot, "05_inital_UMAP", width = 6.5, height = 6)


# Merging and Breaking Clusters 
#######################################################################
# Merging cluster 1 into 2
####
# dat.m = dat
# 
# dat.m <- RenameIdents(
#   object = dat.m,
#   '2' = '1'
# )
# table(Idents(dat.m))
# plot <- DimPlot(dat.m, reduction = "umap")
# plot

# Check features to determine cell types

# FeaturePlot(dat, features = "Kit")
# FeaturePlot(dat, features = "Ano1")
# FeaturePlot(dat, features = "Cd34")

# Breaking down cluster 1 ino 0.75 Res
####
# DimPlot(object = dat,
#         group.by=grep("res",colnames(dat@meta.data),
#                       value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "umap", label = T)
# 
# 
# 
# newIdent = as.character(Idents(dat))
# newIdent[newIdent == '1'] = paste0("R",as.character(dat$RNA_snn_res.0.75[newIdent == '1']))
# 
# Idents(dat) <- as.factor(newIdent)
# table(Idents(dat))
# plot <- DimPlot(object = dat, pt.size=0.5, label = F, reduction = "umap")
# plot
# SaveFigure(plot, "10_clust1_broken", width = 14, height = 12)


# SaveObject(dat.m, "05_seurat_obj_after_merging_clusters")
# #dat <- ReadObject("05_seurat_obj_after_merging_clusters")

#######################################################################
 # Reorder cluster
dat <- BuildClusterTree(dat, reorder = TRUE, reorder.numeric = TRUE)

plot <- PlotClusterTree(dat)
plot
SaveFigure(plot, "06_cluster_tree", width = 5, height = 5)

plot <- DimPlot(dat, reduction = "umap")
          # + NoAxes() + NoLegend()
plot
SaveFigure(plot, "07_final_UMAP", width = 6.2, height = 6)


## UMAP by samples
plot <- DimPlot(object = dat, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
plot
SaveFigure(plot, "07_umap_by_sample", width = 8, height = 8)

# Individual Samples
sample1 <- subset(dat, sample == "sample_1")
plot <- DimPlot(sample1, reduction = "umap")
plot
SaveFigure(plot, "sample_1", width = 8, height = 7)


sample2 <- subset(dat, sample == "sample_2")
plot <- DimPlot(sample2, reduction = "umap")
plot
SaveFigure(plot, "sample_2", width = 8, height = 7)

#########
# See distribution of genes and mito

plot <- FeaturePlot(dat, features = c('nCount_RNA'), pt.size=0.5)
plot
SaveFigure(plot, "07_umap_counts", width = 8, height = 8)
plot <- FeaturePlot(dat, features = c('nFeature_RNA'), pt.size=0.5)
plot
SaveFigure(plot, "07_umap_features", width = 8, height = 8)
plot <- FeaturePlot(dat, features = c('percent.mt'), pt.size=0.5)
plot
SaveFigure(plot, "07_umap_mito", width = 8, height = 8)

SaveObject(dat, "05_seurat_obj_after_umap")
#dat <- ReadObject("05_seurat_obj_after_umap")

#######################################################
# Find Markers (This also takes a while)

dat_markers <- FindAllMarkers(dat, min.pct = 0.25, logfc.threshold = 0.25)
SaveObject(dat_markers, "00_markers")
#dat_markers <- ReadObject("00_markers")

# Heatmap
top20 <- dat_markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
plot <- DoHeatmap(
  object = dat,
  features = top20$gene)
plot
SaveFigure(plot, "08_heatmap_top20_markers", width = 12, height = 20)

# Dot plot
top5 <- dat_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)
plot <- DotPlot(dat, features = to_plot, group.by = "tree.ident") + coord_flip()
plot
SaveFigure(plot, "08_dplot_top5_markers", width = 9, height = 20)


# Plot individual expression of genes
plot1 <- VlnPlot(dat, features = c("Ano1"), group.by = "tree.ident")
plot1
plot2 <- VlnPlot(dat, features = c("Kit"), group.by = "tree.ident")
plot2
SaveFigure(plot1 + plot2, "08_vln_exp_markers", width = 10, height = 8)


## Manual Markers 

FeaturePlot(dat, features = "Kit")

# ICC
markers1 <- c(
  "Kit", "Ano1", "Cd34" # Last one negative
)

ICC <- FeaturePlot(dat, features = markers1)
ICC 
SaveFigure(ICC, "ICC", width = 8, height = 7)


###############################################
## Cell Type Identification 

##### scMRMA

no_cluster_result<-scMRMA(input = dat,
                          species = "Mm",
                          db="panglaodb")

# # Select higher resolutions for more cell types to come through
high_res_result<-scMRMA(input = dat,
                        species = "Mm",
                        db='panglaodb',
                        selfClusters = dat$RNA_snn_res.2.25)



dat<- AddMetaData(dat, no_cluster_result$uniformR$annotationResult, col.name = "CellType")
table(dat$CellType, dat$orig.ident)

#                                   sample_1 sample_2 sample_3 sample_4
# Smooth muscle cells                    39      149      309      113
# Fibroblasts                            72      173      325      131
# Keratinocytes                          43       45       44        8
# Endothelial cells                      13       31       46       33
# Macrophages                            12       28       37       12
# Enterocytes                             9       14       22       15
# Pulmonary alveolar type II cells        4        7       27       17
# Schwann cells                           8        4       27        7
# Interneurons                            0        1       14        2




plot <- DimPlot(dat, group.by = "CellType", label=T)
plot
SaveFigure(plot, "09_umap_CellType", width = 8, height = 6)

table(no_cluster_result$uniformR$annotationResult)
table(dat$CellType, dat$seurat_clusters)
table(dat$seurat_clusters, dat$CellType)




a<-DimPlot(dat, group.by = "CellType", label = T, repel=T,combine=T) +
  ggtitle("Cell Type")
a

b<-DimPlot(object = dat, pt.size=0.5, reduction = "umap", label = T)+
  ggtitle("Clusters")
b
c<-DimPlot(object = dat, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
c
a+b
SaveFigure(a + b, "09_umap_Cluster_CellType", width = 16, height = 6)

SaveObject(dat, "06_dat_after_scSMRA")
#SaveObject(dat, "06_dat_after_scSMRA")



