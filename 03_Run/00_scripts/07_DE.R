## Differential Expression
# PI:   
# Date:  
# Title: 


####################

# Load required packages
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)

source("/Users/cassandrahui/Documents/Projects/Scripts/scRNAseq/functions/functions.R")
source("/Users/cassandrahui/Documents/Projects/Scripts/scRNAseq/functions/analyze_subset.R")

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
fig_path <- "../05_figures/"
###################

#dat <- readRDS("Seurat_CellType.rds")
dat <- readRDS("Seurat_Merge.rds")



# Differential expression preformed with Seurat FindMarkers

# View metadata
head(dat@meta.data)

# Assign groups to compare
# Or use orig.ident for individual samples

dat$group <- NA
dat$group[dat$orig.ident == "sample_1" | dat$orig.ident == "sample_2"] <- "control"
dat$group[dat$orig.ident == "sample_3" | dat$orig.ident == "sample_4"] <- "treatment"
head(dat@meta.data)

# Run function with object, column to compare, type to compare, and specific
results <- analyze_subset(dat,
                         group_col = "group",
                         subset_col = "CellType",
                         subset_value = "Fibroblasts")

# Visualize
head(results$comparisons[[1]])

# Save single comparison
comparison_name <- names(results$comparisons)[1]  # Get the name of first comparison
write.csv(results$comparisons[[comparison_name]], 
          file = paste0("../04_outputs/DE_results_Fibroblasts_", comparison_name, ".csv"))

# Or to save all comparisons as separate files
for (name in names(results$comparisons)) {
  write.csv(results$comparisons[[name]], 
            file = paste0("../04_outputs/DE_results_Fibroblasts_", name, ".csv"))
  cat("Saved:", name, "\n")
}

