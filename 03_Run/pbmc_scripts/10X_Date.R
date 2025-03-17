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



plot_cellranger_cells <- function(ind){
  xbreaks = c(1,1e1,1e2,1e3,1e4,1e5,1e6)
  xlabels = c("1","10","100","1000","10k","100K","1M")
  ybreaks = c(1,2,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000,100000,200000,500000,1000000)
  ylabels = c("1","2","5","10","2","5","100","2","5","1000","2","5","10k","2","5","100K","2","5","1M")
  
  pl1 <- data.frame(index=seq.int(1,ncol(d10x.data[[ind]])),
                    nCount_RNA = sort(Matrix:::colSums(d10x.data[[ind]])+1,decreasing=T),
                    nFeature_RNA = sort(Matrix:::colSums(d10x.data[[ind]]>0)+1,decreasing=T)) %>%
    ggplot() +
    scale_color_manual(values=c("red2","blue4"), labels=c("Features", "UMI"), name=NULL) +
    ggtitle(paste("CellRanger filltered cells:",ids[ind],sep=" ")) + xlab("Barcodes") + ylab("counts (UMI or Features") +
    scale_x_continuous(trans = 'log2', breaks=xbreaks, labels = xlabels) +
    scale_y_continuous(trans = 'log2', breaks=ybreaks, labels = ylabels) +
    geom_line(aes(x=index, y=nCount_RNA, color = "UMI"), size=1.75) +
    geom_line(aes(x=index, y=nFeature_RNA, color = "Features"), size=1.25)
  
  return(pl1)
}

plot_cellranger_cells(1)
plot_cellranger_cells(2)
plot_cellranger_cells(3)
plot_cellranger_cells(4)


experiment.data <- do.call("cbind", d10x.data)

pbmc <- CreateSeuratObject(
  experiment.data,
  project = experiment_name,
  min.cells = 10,
  min.features = 100,
  names.field = 2,
  names.delim = "\\_")

pbmc
str(pbmc)


SaveObject(pbmc, "01_seurat_obj_before_QC")
#pbmc <- ReadObject("01_seurat_obj_before_QC")

slotNames(pbmc)
head(pbmc[[]])

set.seed(12345)
kable(do.call("cbind", tapply(pbmc$nFeature_RNA, 
                              Idents(pbmc),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of Genes/Cell by Sample") %>% kable_styling()


kable(do.call("cbind", tapply(pbmc$nCount_RNA, 
                              Idents(pbmc),quantile,probs=seq(0,1,0.05))),
      caption = "5% Quantiles of UMI/Cell by Sample") %>% kable_styling()


# View Mitochondrial percentatge (MT- for Human)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
plot <- VlnPlot(pbmc, pt.size = 0.10,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot
SaveFigure(plot, "01_vln_QC", width = 18, height = 12)


# View 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 +plot2
SaveFigure((plot1 + plot2),"01_scatter_beforeQC", width = 13, height = 6, res = 200)

# Number of cells before QC
table(pbmc$orig.ident)

# sample_1 sample_2 sample_3 sample_4 
# 200      452      851      338     total: 1841



## Cut to see data better
pbmc <- subset(pbmc, subset = nFeature_RNA < 6000 & percent.mt < 40)

# Plot again
# Edit lines to determine best cut
plot <- FeatureScatter(
  pbmc, "nCount_RNA", "nFeature_RNA",
  pt.size = 0.5, shuffle = TRUE)  + geom_vline(xintercept = c(200,15000)) + geom_hline(yintercept = c(150, 3000))
plot
SaveFigure((plot),"02_QC_cuts", width = 7, height = 6, res = 200)
table(pbmc$orig.ident)

# Number of cells after extreme cuts
#
# sample_1 sample_2 sample_3 sample_4 
# 752     1227     2110      622     total: 4711


plot <- RidgePlot(pbmc, features=c("nFeature_RNA","nCount_RNA", "percent.mt"))
plot
SaveFigure(plot, "01_ridge_beforeQC_persample", width = 15, height = 10)
plot(sort(Matrix::rowSums(GetAssayData(pbmc) >= 3), decreasing = TRUE) , xlab="gene rank", ylab="number of cells", main="Cells per genes (reads/gene >= 3 )")

# Perform the filtering
pbmc <- subset(pbmc, subset = nCount_RNA > 200 & nCount_RNA < 15000 & nFeature_RNA > 150 & nFeature_RNA < 3000 & percent.mt < 25)

table(pbmc$orig.ident)
# Number of cells after filtering

# sample_1 sample_2 sample_3 sample_4 
# 200      452      851      338     total: 1841

plot <- RidgePlot(pbmc, features=c("nFeature_RNA","nCount_RNA", "percent.mt"), ncol = 2)
plot
SaveFigure(plot, "02_ridge_QC_persample_postfilter", width = 12, height = 10)

####################################################
#Normalize
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Find Variable Features
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
SaveFigure((plot1 + plot2), "03_var_features", width = 12, height = 6)

SaveObject(pbmc, "02_normalized_pbmc")
#pbmc <- ReadObject("02_normalized_pbmc")

####################################################
#Scale and Run PCA
pbmc <- ScaleData(pbmc)

pbmc <- RunPCA(pbmc, npcs = 100)
SaveObject(pbmc, "03_seurat_obj_after_PCA")
#pbmc <- ReadObject("03_seurat_obj_after_PCA")

####################################################
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

plot <- DimPlot(object = pbmc, reduction = "pca")
plot
SaveFigure(plot, "03_PCA", width = 8, height = 7.5)
plot <- VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
plot
SaveFigure(plot, "03_viz_PCA_loadings", width = 10, height = 8)



# Image doesn't save as png unless fast = FALSE
plot <- DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
plot
SaveFigure(plot, "03_dim_heatmap1", width = 6, height = 6)


plot <- DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
plot
SaveFigure(plot, "03_dim_heatmap1_15", width = 12, height = 18)


# View to determine best number of PCs or Run Jackstraw 
plot <- ElbowPlot(pbmc,ndims = 50)
plot
SaveFigure(plot, "03_PC_elbow_plot", width = 8, height = 8)

####################################
# This will take a long time with large data

## Jackstraw
############

pbmc <- JackStraw(pbmc, dims = 100) # Check this is right
#pbmc <- JackStraw(pbmc, num.replicate = 100) 
pbmc <- ScoreJackStraw(pbmc, dims = 1:100)

plot <- JackStrawPlot(object = pbmc, dims = 1:100) + theme(legend.position="bottom")
plot
SaveFigure(plot, "03_Jack_Straw", width = 12, height = 8)

SaveObject(pbmc, "04_seurat_obj_after_Jackstraw")
#pbmc <- ReadObject("04_seurat_obj_after_Jackstraw")
############

# Select pcs
use.pcs <- 1:14
pbmc <- FindNeighbors(pbmc, dims = use.pcs)
# Create numtiple resolutions to view
pbmc <- FindClusters(pbmc, resolution = seq(0.25, 4, 0.5), verbose = FALSE)


pbmc <- RunTSNE(
  object = pbmc,
  reduction.use = "pca",
  dims = use.pcs,
  do.fast = TRUE)



plot <- DimPlot(object = pbmc,
                group.by=grep("res",colnames(pbmc@meta.data),
                              value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "tsne", label = T)
plot
SaveFigure(plot, "04_Resolutions", width = 14, height = 12)


# Select resolution to use
Idents(pbmc) <- "RNA_snn_res.0.25"
pbmc <- RunUMAP(pbmc, dims = use.pcs)
plot <- DimPlot(pbmc, reduction = "umap")
plot
SaveFigure(plot, "05_inital_UMAP", width = 6.5, height = 6)


# Merging and Breaking Clusters 
#######################################################################
# Merging cluster 1 into 2
####
# pbmc.m = pbmc
# 
# pbmc.m <- RenameIdents(
#   object = pbmc.m,
#   '2' = '1'
# )
# table(Idents(pbmc.m))
# plot <- DimPlot(pbmc.m, reduction = "umap")
# plot

# Check features to determine cell types

# FeaturePlot(pbmc, features = "Kit")
# FeaturePlot(pbmc, features = "Ano1")
# FeaturePlot(pbmc, features = "Cd34")

# Breaking down cluster 1 ino 0.75 Res
####
# DimPlot(object = pbmc,
#         group.by=grep("res",colnames(pbmc@meta.data),
#                       value = TRUE)[1:4], ncol=2 , pt.size=3.0, reduction = "umap", label = T)
# 
# 
# 
# newIdent = as.character(Idents(pbmc))
# newIdent[newIdent == '1'] = paste0("R",as.character(pbmc$RNA_snn_res.0.75[newIdent == '1']))
# 
# Idents(pbmc) <- as.factor(newIdent)
# table(Idents(pbmc))
# plot <- DimPlot(object = pbmc, pt.size=0.5, label = F, reduction = "umap")
# plot
# SaveFigure(plot, "10_clust1_broken", width = 14, height = 12)


# SaveObject(pbmc.m, "05_seurat_obj_after_merging_clusters")
# #pbmc <- ReadObject("05_seurat_obj_after_merging_clusters")

#######################################################################
# Reorder cluster
pbmc <- BuildClusterTree(pbmc, reorder = TRUE, reorder.numeric = TRUE)

plot <- PlotClusterTree(pbmc)
plot
SaveFigure(plot, "06_cluster_tree", width = 5, height = 5)

plot <- DimPlot(pbmc, reduction = "umap")
# + NoAxes() + NoLegend()
plot
SaveFigure(plot, "07_final_UMAP", width = 6.2, height = 6)


## UMAP by samples
plot <- DimPlot(object = pbmc, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
plot
SaveFigure(plot, "07_umap_by_sample", width = 8, height = 8)

# Individual Samples
sample1 <- subset(pbmc, sample == "sample_1")
plot <- DimPlot(sample1, reduction = "umap")
plot
SaveFigure(plot, "sample_1", width = 8, height = 7)


sample2 <- subset(pbmc, sample == "sample_2")
plot <- DimPlot(sample2, reduction = "umap")
plot
SaveFigure(plot, "sample_2", width = 8, height = 7)

#########
# See distribution of genes and mito

plot <- FeaturePlot(pbmc, features = c('nCount_RNA'), pt.size=0.5)
plot
SaveFigure(plot, "07_umap_counts", width = 8, height = 8)
plot <- FeaturePlot(pbmc, features = c('nFeature_RNA'), pt.size=0.5)
plot
SaveFigure(plot, "07_umap_features", width = 8, height = 8)
plot <- FeaturePlot(pbmc, features = c('percent.mt'), pt.size=0.5)
plot
SaveFigure(plot, "07_umap_mito", width = 8, height = 8)

SaveObject(pbmc, "05_seurat_obj_after_umap")
#pbmc <- ReadObject("05_seurat_obj_after_umap")

#######################################################
# Find Markers (This also takes a while)

pbmc_markers <- FindAllMarkers(pbmc, min.pct = 0.25, logfc.threshold = 0.25)
SaveObject(pbmc_markers, "00_markers")
#pbmc_markers <- ReadObject("00_markers")

# Heatmap
top20 <- pbmc_markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
plot <- DoHeatmap(
  object = pbmc,
  features = top20$gene)
plot
SaveFigure(plot, "08_heatmap_top20_markers", width = 12, height = 20)

# Dot plot
top5 <- pbmc_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)
plot <- DotPlot(pbmc, features = to_plot, group.by = "tree.ident") + coord_flip()
plot
SaveFigure(plot, "08_dplot_top5_markers", width = 9, height = 20)


# Plot individual expression of genes
plot1 <- VlnPlot(pbmc, features = c("Ano1"), group.by = "tree.ident")
plot1
plot2 <- VlnPlot(pbmc, features = c("Kit"), group.by = "tree.ident")
plot2
SaveFigure(plot1 + plot2, "08_vln_exp_markers", width = 10, height = 8)


## Manual Markers 

FeaturePlot(pbmc, features = "Kit")

# ICC
markers1 <- c(
  "Kit", "Ano1", "Cd34" # Last one negative
)

ICC <- FeaturePlot(pbmc, features = markers1)
ICC 
SaveFigure(ICC, "ICC", width = 8, height = 7)


###############################################
## Cell Type Identification 

##### scMRMA

no_cluster_result<-scMRMA(input = pbmc,
                          species = "Mm",
                          db="panglaodb")

# # Select higher resolutions for more cell types to come through
high_res_result<-scMRMA(input = pbmc,
                        species = "Mm",
                        db='panglaodb',
                        selfClusters = pbmc$RNA_snn_res.2.25)



pbmc<- AddMetaData(pbmc, no_cluster_result$uniformR$annotationResult, col.name = "CellType")
table(pbmc$CellType, pbmc$orig.ident)

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




plot <- DimPlot(pbmc, group.by = "CellType", label=T)
plot
SaveFigure(plot, "09_umap_CellType", width = 8, height = 6)

table(no_cluster_result$uniformR$annotationResult)
table(pbmc$CellType, pbmc$seurat_clusters)
table(pbmc$seurat_clusters, pbmc$CellType)




a<-DimPlot(pbmc, group.by = "CellType", label = T, repel=T,combine=T) +
  ggtitle("Cell Type")
a

b<-DimPlot(object = pbmc, pt.size=0.5, reduction = "umap", label = T)+
  ggtitle("Clusters")
b
c<-DimPlot(object = pbmc, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
c
a+b
SaveFigure(a + b, "09_umap_Cluster_CellType", width = 16, height = 6)

SaveObject(pbmc, "06_pbmc_after_scSMRA")
#SaveObject(pbmc, "06_pbmc_after_scSMRA")

##############################################################
## Gene Ontology
#############################################################

# Gene Ontology

# View current sample names
unique(pbmc$orig.ident)

# Identifiy a cluster or cell type
###################
# Cluster 5: Fibroblasts
cluster5 <- subset(pbmc, idents = '5')
expr <- as.matrix(GetAssayData(cluster5))

# # Cell type: Interneuron
# ICC <- subset(pbmc, subset = CellType == "Interneurons")
# expr <- as.matrix(GetAssayData(ICC))


# Filter out genes that are 0 for every cell in this cluster

summary(rowSums(expr))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.985   48.537   85.471  146.666  164.874 4718.396 

# If there are none, DO NOT run this code 
# bad <- which(rowSums(expr) == 0)
# expr <- expr[-bad,]

# Select genes that are expressed > 0 in at least half of cells
n.gt.0 <- apply(expr, 1, function(x)length(which(x > 0)))
expressed.genes <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.5)]
all.genes <- rownames(expr)

# Define geneList as 1 if gene is in expressed.genes, 0 otherwise
geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
names(geneList) <- all.genes

#########
# Create topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP", # use biological process ontology
              allGenes = geneList,
              geneSelectionFun = function(x)(x == 1),
              annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")

# Test for enrichment using Fisher's Exact Test
resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
clust5_BP <- data.frame(GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60))
# Save as a CSV file
write.csv(clust5_BP, "../07_reports/clust5_BP.csv", row.names = FALSE)

############################################
# Differential Expression
###########################################
## DE for Fibroblasts

# Filter genes to those expressed in at least 10% of cells
keep <- rownames(expr)[which(n.gt.0/ncol(expr) >= 0.1)]
expr2 <- expr[keep,]

# Set up "design matrix" with statistical model
cluster5$proper.ident <- make.names(cluster5$orig.ident)
mm <- model.matrix(~0 + proper.ident + percent.mt + nFeature_RNA, data = cluster5[[]])
head(mm)

# Fit model in limma
fit <- lmFit(expr2, mm)
head(coef(fit))

# Test 'sample_1' - 'sample_2'
contr <- makeContrasts(proper.identsample_1 - proper.identsample_2, levels = colnames(coef(fit)))
contr

fit2 <- contrasts.fit(fit, contrasts = contr)
fit2 <- eBayes(fit2)
out <- topTable(fit2, n = Inf, sort.by = "P")
head(out, 30)
write.csv(out, "../07_reports/DE_between_samps1_2_in_clust5_fibroblasts.csv", row.names = TRUE)

## View marker genes 
FeaturePlot(sample1, features = c("Rn18s-rs5", "Fos"))
FeaturePlot(sample2, features = c("Rn18s-rs5", "Fos"))

FeaturePlot(sample1, features = c("Nr4a1", "Mir6236"))
FeaturePlot(sample2, features = c("Nr4a1", "Mir6236"))






