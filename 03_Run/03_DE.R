# PI:   
# Date:  
# Title: 
# Alingment Type: 

# Organism: 
# Conditions:
# Rep: 
# Design:

####################

# Load required packages
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)

source("/Users/cassandrahui/Documents/Projects/Scripts/scRNAseq/functions/analyze_subset.R")

table(dat$orig.ident)
res <- analyze_subset(dat, "orig.ident", "tree.ident", "2")
head(res$comparisons)


