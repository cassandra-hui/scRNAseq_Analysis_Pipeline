# PI:   
# Date:  
# Title: 

library(Seurat)
library(dplyr)
library(kableExtra)
library(ggplot2)

source('functions.R')


###################
# Set for functions to work
#setwd("~/Documents/Projects/")
data_path <- "../01_data/DGE_filtered"
fig_path <- "../05_figures"
###################

# Identifying
experiment_name = "Misty's AstyMax scRNAseq"
dataset_loc <- "../01_data"
ids <- c("M", "P", "S", "T")

experiment.metrics <- read.csv(file.path(dataset_loc, ids[1], "outs/metrics_summary.csv"))
sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))
rownames(sequencing.metrics) <- gsub("\\.", " ", rownames(sequencing.metrics))
colnames(sequencing.metrics) <- "All samples"
sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")


d10x.metrics <- lapply(ids, function(i){
  metrics <- read.csv(file.path(dataset_loc,paste0(i,"/outs"),"metrics_summary.csv"), colClasses = "character")
})
experiment.metrics <- do.call("rbind", d10x.metrics)
rownames(experiment.metrics) <- ids

sequencing.metrics <- data.frame(t(experiment.metrics[,c(1:19)]))

row.names(sequencing.metrics) <- gsub("\\."," ", rownames(sequencing.metrics))

sequencing.metrics %>%
  kable(caption = 'Cell Ranger Results') %>%
  pack_rows("Overview", 1, 3, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Sequencing Characteristics", 4, 9, label_row_css = "background-color: #666; color: #fff;") %>%
  pack_rows("Mapping Characteristics", 10, 19, label_row_css = "background-color: #666; color: #fff;") %>%
  kable_styling("striped")


d10x.data <- lapply(ids, function(i){
  d10x <- Read10X_h5(file.path(dataset_loc, i, "/outs","filtered_feature_bc_matrix.h5"))
  colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="_")
  d10x
})
names(d10x.data) <- ids

str(d10x.data)



plot_cellranger_cells(1)
plot_cellranger_cells(2)
plot_cellranger_cells(3)
plot_cellranger_cells(4)


experiment.data <- do.call("cbind", d10x.data)

dat <- CreateSeuratObject(
  experiment.data,
  project = experiment_name,
  min.cells = 3,
  min.features = 100,
  names.field = 2,
  names.delim = "\\_")

dat
str(dat)


saveRDS(dat, "02_objects/01_seurat_obj_before_QC.rds")
#dat <- loadRDS("02_objects/01_seurat_obj_before_QC.rds")

