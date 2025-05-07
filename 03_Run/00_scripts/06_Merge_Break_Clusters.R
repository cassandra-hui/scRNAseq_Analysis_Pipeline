# Merging and Breaking Clusters
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

dat <- readRDS("02_objects/Seurat_CellType.rds")




# Merging and Breaking Clusters 
#######################################################################
# Merging cluster 1 into 2
####

dat <- RenameIdents(
  object = dat,
  '2' = '1'
)
table(Idents(dat))
plot <- DimPlot(dat, reduction = "umap")
plot

# Check features to determine cell types

FeaturePlot(dat, features = "Kit")
FeaturePlot(dat, features = "Ano1")
FeaturePlot(dat, features = "Cd34")

# Breaking down cluster 1 into 0.8 Res
####
DimPlot(object = dat,
        group.by=grep("res",colnames(dat@meta.data),
                      value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "umap", label = T)


# Rename cluster 1 to R1
newIdent = as.character(Idents(dat))
newIdent[newIdent == '1'] = paste0("R",as.character(dat$SCT_snn_res.0.8[newIdent == '1']))

Idents(dat) <- as.factor(newIdent)
table(Idents(dat))

# Merge other clusters accidentally added

dat <- RenameIdents(
  object = dat,
  'R3' = 'R1',
  'R9' = 'R1'
)

plot <- DimPlot(object = dat, pt.size=0.5, label = F, reduction = "umap")
plot
SaveFigure(plot, "10_clust1_broken", width = 14, height = 12)


saveRDS(dat, "02_objects/Seurat_Merge.rds")

