######################## 
# Intergration
########################

# !!!

## Only use if there is a batch effect, different datasets, 
## or different species you want to align cell types between



# PI:   
# Date:  
# Title: 

library(Seurat)


source("/Users/cassandrahui/Documents/Projects/Scripts/scRNAseq/functions/functions.R")

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
fig_path <- "../05_figures/"
###################

use.pcs <- 1:25
dat <- readRDS("02_objects/Seurat_DR.rds")

# Confirm there is a clear sepparation in samples
DimPlot(dat, reduction = "umap", 
        group.by = "orig.ident")


# There are two ways to do this
# If the first option errors use the second



# Option 1

####################################################
### Split by samples and integrate, then rejoin

Assays(dat)
#DefaultAssay(dat) <- "RNA"

# Use "SCT" for SCTransform or "RNA" for manual Normalization, Find VarFeatures, and Scale
dat[["SCT"]] <- split(dat[["SCT"]], f = dat$orig.ident)

features <- VariableFeatures(dat)

dat <- IntegrateLayers(object = dat, 
                       method = CCAIntegration, 
                       orig.reduction = "pca", 
                       new.reduction = "integrated.cca", 
                       features = features, 
                       verbose = FALSE)

# re-join layers after integration
dat[["SCT"]] <- JoinLayers(dat[["SCT"]])


### Run UMAP again with appropriate "reduction"


dat <- RunUMAP(dat, dims = use.pcs, reduction = "integrated.cca")

DimPlot(dat, reduction = "umap", group.by = "orig.ident")



# If needed, set identities to your samples
Idents(dat) <- "orig.ident"

# Then you can plot without group.by
DimPlot(dat, reduction = "umap")

saveRDS(dat, "Seurat_Intergrated.rds")

##########################################



# Option 2

##########################################
## Integrating multiple objects

# Split your dataset of you do not have multiple objects

batch1 <- subset(dat, orig.ident == "sample_1" | orig.ident == "sample_2")
batch2 <- subset(dat, orig.ident == "sample_3" | orig.ident == "sample_4")

table(batch1$orig.ident)
table(batch2$orig.ident)


my_datasets <- list(batch1, batch2) # You can have more than 2


## Find Anchors
common_variable_features <- SelectIntegrationFeatures(my_datasets) # nfeatures = 2000 by default

head(common_variable_features)
length(common_variable_features)


mutual_nearest_neighbors <- FindIntegrationAnchors(my_datasets, 
                                                   anchor.features = common_variable_features) 
# reduction = "cca" by default

mutual_nearest_neighbors
head(mutual_nearest_neighbors@anchors)

## Intergrate

my_integrated_dataset <- IntegrateData(mutual_nearest_neighbors)
my_integrated_dataset
my_integrated_dataset@assays

# Select new assay
DefaultAssay(my_integrated_dataset) <- "integrated"

# Rerun Scale, PCA, and UMAP
my_integrated_dataset <- ScaleData(my_integrated_dataset, verbose = FALSE)
my_integrated_dataset <- RunPCA(my_integrated_dataset, verbose = FALSE)

plot <- DimPlot(object = my_integrated_dataset, reduction = "pca")
plot
SaveFigure(plot, "04_PCA_intergration", width = 8, height = 7.5)

plot <- ElbowPlot(my_integrated_dataset, ndims = 50)
plot
SaveFigure(plot, "04_PC_elbow_plot_intergration", width = 8, height = 8)

# Using 25
use.pcs <- 1:25


set.seed(42)
my_integrated_dataset <- RunUMAP(my_integrated_dataset, dims = use.pcs, verbose = FALSE)
head(my_integrated_dataset@meta.data)

plot <- DimPlot(my_integrated_dataset, alpha = 0.5)
plot
SaveFigure(plot, "05_Intergration_UMAP", width = 6, height = 6)


plot <- DimPlot(my_integrated_dataset, group.by = "SCT_snn_res.0.1")
plot
SaveFigure(plot, "05_UMAP_intergration", width = 6.5, height = 6)


saveRDS(my_integrated_dataset, "Seurat_Intergrated.rds")

