### Dimensionality Reduction
# PI:   
# Date:  
# Title: 

library(Seurat)
library(clustree)

source('functions.R')

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
fig_path <- "../05_figures/"
###################


dat <- readRDS("02_objects/Seurat_Norm.rds")
Idents(dat) <- dat$orig.ident

######## Run PCA

dat <- RunPCA(dat)


print(dat[["pca"]], dims = 1:5, nfeatures = 5)

plot <- DimPlot(object = dat, reduction = "pca")
plot
SaveFigure(plot, "04_PCA", width = 8, height = 7.5)


# View to determine best number of PCs pick the elbow
# Too few will give less defined clusters, too many will be computationally heavy
# You can select a few different options to see at what point the TSNE or UMAP are not effected by more PCs

plot <- ElbowPlot(dat,ndims = 50)
plot
SaveFigure(plot, "04_PC_elbow_plot", width = 8, height = 8)

# Using 25
use.pcs <- 1:25

############# Run TSNE

dat <- RunTSNE(dat, dims = use.pcs)

plot <- DimPlot(dat, reduction = "tsne")
plot


dat <- FindNeighbors(dat, dims = use.pcs)

dat <- FindClusters(dat, resolution = seq(0.1, 0.8, by=0.1))
head(dat@meta.data)


# Too see association of clusters
clustree(dat@meta.data[,grep("SCT_snn_res", colnames(dat@meta.data))],
                   prefix = "SCT_snn_res.")


# Visualize TSNE
DimPlot(dat, group.by = "SCT_snn_res.0.1")
DimPlot(dat, group.by = "SCT_snn_res.0.3")


################ Run UMAP

dat <- RunUMAP(dat, dims = use.pcs)

DimPlot(dat, group.by = "SCT_snn_res.0.1")
DimPlot(dat, group.by = "SCT_snn_res.0.3")

# Set resolution to use
Idents(dat) <- "SCT_snn_res.0.1"


plot <- DimPlot(dat, reduction = "umap")
plot
SaveFigure(plot, "05_inital_UMAP", width = 6.5, height = 6)


## UMAP by samples
plot <- DimPlot(object = dat, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
plot
SaveFigure(plot, "05_umap_by_sample", width = 8, height = 8)


saveRDS(dat, "02_objects/Seurat_DR.rds")

