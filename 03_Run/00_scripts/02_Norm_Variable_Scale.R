## ONLY IF NOT USING SCTRANFORM
# PI:   
# Date:  
# Title: 

library(Seurat)

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
fig_path <- "../05_figures/"
###################

dat <- readRDS("02_objects/Seurat_QC.rds")

## Normalize, Find Variable Features, and Scale if not using SCTransfrom

# Rename the assay from "originalexp" to "RNA"
DefaultAssay(dat)
dat[["RNA"]] <- dat[["originalexp"]]
Assays(dat)
DefaultAssay(dat) <- "RNA"
dat[["originalexp"]] <- NULL  # Remove the old assay

# Now check the assay name
Assays(dat)  # Should be "RNA"


# Log Normailze
dat <- NormalizeData(dat,
                     normalization.method = "LogNormalize",
                     scale.factor = 10000)

# Find Variable Features
dat <- FindVariableFeatures(dat, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dat), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(dat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
SaveFigure((plot2), "03_var_features", width = 6, height = 6)


# Scale data
dat <- ScaleData(dat)


saveRDS(dat, "02_objects/Seurat_Norm.rds")