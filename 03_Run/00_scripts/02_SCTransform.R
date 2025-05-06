## SCTransfrom
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

dat <- readRDS("02_objects/Seurat_QC.rds")

### Using SCTransform to Normalize, Find Variable Features, and Scale
# This is supposed to work better for scRNA data, 
# but can use the regular Normalize and Scale script instead

# Rename the assay from "originalexp" to "RNA"
DefaultAssay(dat)
dat[["RNA"]] <- dat[["originalexp"]]
Assays(dat)
DefaultAssay(dat) <- "RNA"
dat[["originalexp"]] <- NULL  # Remove the old assay

# Now check the assay name
Assays(dat)  # Should be "RNA"


# May error if data is too large. 
dat <- SCTransform(dat)


saveRDS(dat, "02_objects/Seurat_Norm.rds")
