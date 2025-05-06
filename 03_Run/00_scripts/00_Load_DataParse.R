# PI:   
# Date:  
# Title: 

library(Seurat)

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
dat <- CreateSeuratObject(mat, 
                          min.cells = 3,
                          min.features = 100,
                          names.field = 0, 
                          meta.data = cell_meta)



# Rename orig.ident by samples
dat@meta.data$orig.ident <- dat@meta.data$sample
Idents(dat) <- dat@meta.data$orig.ident


saveRDS(dat, "02_objects/01_seurat_obj_before_QC.rds")
#dat <- loadRDS("02_objects/01_seurat_obj_before_QC.rds")
