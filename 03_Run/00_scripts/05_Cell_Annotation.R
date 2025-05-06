## Cell Annotation
# PI:   
# Date:  
# Title: 

library(Seurat)
library(celldex)
library(dittoSeq)
library(SingleR)
library(scMRMA)


source("/Users/cassandrahui/Documents/Projects/Scripts/scRNAseq/functions/functions.R")

###################
# Set for functions to work
#setwd("~/Documents/Projects/")
fig_path <- "../05_figures/"
###################


#dat <- readRDS("02_objects/Seurat_DR.rds")
dat <- readRDS("02_objects/Seurat_Intergrated.rds")


###################################
# You can decide on manual or automatic cell annotations

# SingleR uses references you select
# Made by Seurat with Human and Mouse data (https://azimuth.hubmapconsortium.org/)

# With mannual annotation you can use, Seurat AddModuleScore()
# But if you want to incorporate negative genes you need to use UCell (on Bioconductor, https://github.com/carmonalab/UCell)

###################################


#########################
### Cell Cycle Scoring
########################


# Extract built in genes
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes


dat <- CellCycleScoring(dat,
                        s.features = s.genes,
                        g2m.features = g2m.genes)

head(dat@meta.data)

plot <- DimPlot(dat, group.by = "Phase")
plot
SaveFigure(plot, "06_CellCylce_UMAP", width = 5, height = 5)


###########################
## Manual Identification ##
###########################


# To help with manual identification you can provide a list of marker genes for each cluster
#############
# Find Markers (This takes a while)

dat_markers <- FindAllMarkers(dat, min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(dat_markers, "02_objects/markers")

# Heatmap
top20 <- dat_markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
plot <- DoHeatmap(
  object = dat,
  features = top20$gene)
plot
SaveFigure(plot, "07_heatmap_top20_markers", width = 12, height = 20)

# Dot plot
top5 <- dat_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
to_plot <- unique(top5$gene)
plot <- DotPlot(dat, features = to_plot, group.by = "SCT_snn_res.0.1") + coord_flip()
plot
SaveFigure(plot, "07_dplot_top5_markers", width = 9, height = 20)


#########


# Set default clusters
dat <- SetIdent(dat, value = dat$SCT_snn_res.0.1)

# Use original count data not integrated or SCT data
DefaultAssay(dat) <- "RNA"

# To view one gene
FeaturePlot(dat, "C3")
VlnPlot(dat, features = c("Ano1"))


## Manually define cell types

tcell_genes <- c("Il7r", "Ltb", "Trac", "Cd3d")
monocyte_genes <- c("Cd14", "Cst3", "Cd68", "Ctss")
FeaturePlot(dat, tcell_genes, ncol=2)
FeaturePlot(dat, monocyte_genes, ncol=2)

VlnPlot(dat,
        features = tcell_genes,
        ncol = 2)

VlnPlot(dat,
        features = monocyte_genes,
        ncol = 2)


# Calculate the score
dat <- AddModuleScore(dat,
                      features = list(monocyte_genes),
                      name = "monocyte_genes")

# A column was added to seu@meta.data.
head(dat@meta.data)

plot1 <- FeaturePlot(dat, "monocyte_genes1")
plot2<- VlnPlot(dat, "monocyte_genes1")
plot1 + plot2
SaveFigure(plot1 + plot2, "08_monocytes", width = 9, height = 20)

# Rename Clusters
#########

# First examine your current clusters
table(dat$SCT_snn_res.0.1)

# Create a new column for cell type annotations
dat$cell_type <- "Unassigned" # Default value

# Assign cell types based on clusters (adjust cluster numbers based on your data)
dat$cell_type[dat$SCT_snn_res.0.1 == 0] <- "T Cells"
dat$cell_type[dat$SCT_snn_res.0.1 == 1] <- "Monocytes"
dat$cell_type[dat$SCT_snn_res.0.1 == 2] <- "B Cells"
# Continue for other clusters...

# Visualize with the new annotations
plot <- DimPlot(dat, group.by = "cell_type", label = TRUE) + NoLegend()
plot
SaveFigure(plot, "09_CellType_UMAP", width = 5, height = 5)

saveRDS(dat, "Seurat_CellType.rds")

##############################
## Automated Identification ##
##############################

# Options 1: SingleR
######################

## Selcet a reference dataset from celldex
# https://bioconductor.org/packages/devel/data/experiment/manuals/celldex/man/celldex.pdf

# Using standard Mouse here

ref <- MouseRNAseqData()
class(ref)
# Show cell types
table(ref$label.main)


dat_SingleR <- SingleR::SingleR(test = Seurat::GetAssayData(dat),
                                ref = ref,
                                labels = ref$label.main)

head(dat_SingleR)

SingleR::plotScoreHeatmap(dat_SingleR)
SingleR::plotDeltaDistribution(dat_SingleR)


## Remove cell types with less than 10 cells so not to clog our plots
singleR_labels <- dat_SingleR$labels
t <- table(singleR_labels)
other <- names(t)[t < 10]
singleR_labels[singleR_labels %in% other] <- "none"

## Add to object
dat$SingleR_annot <- singleR_labels


dittoSeq::dittoDimPlot(dat, "SingleR_annot", size = 0.7)
dittoSeq::dittoBarPlot(dat, var = "SingleR_annot", group.by = "orig.ident")

dittoSeq::dittoBarPlot(dat, 
                       var = "SingleR_annot", 
                       group.by = "SCT_snn_res.0.1")

saveRDS(dat, "Seurat_CellType.rds")


# Options 2: scMRMA
######################

# Based on genes only
no_cluster_result<-scMRMA(input = dat,
                          species = "Mm", # Hs for Humans
                          db="panglaodb")

# Select higher resolutions for more cell types to come through
high_res_result<-scMRMA(input = dat,
                        species = "Mm",
                        db='panglaodb',
                        selfClusters = dat$SCT_snn_res.0.8)



dat<- AddMetaData(dat, no_cluster_result$uniformR$annotationResult, col.name = "CellType")
table(dat$CellType, dat$orig.ident)

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


dat<- AddMetaData(dat, high_res_result$uniformR$annotationResult, col.name = "CellType2")
table(dat$CellType2, dat$orig.ident)


plot <- DimPlot(dat, alpha = 0.8, group.by = "CellType", label=T)
plot
SaveFigure(plot, "09_umap_CellType", width = 8, height = 6)

table(no_cluster_result$uniformR$annotationResult)
table(dat$CellType, dat$seurat_clusters)
table(dat$seurat_clusters, dat$CellType)




a<-DimPlot(dat, group.by = "CellType", label = T, repel=T,combine=T) +
  ggtitle("Cell Type")
a

b<-DimPlot(object = dat, pt.size=0.5, reduction = "umap", label = T)+
  ggtitle("Clusters")
b
c<-DimPlot(object = dat, group.by="orig.ident", pt.size=0.5, reduction = "umap", shuffle = TRUE)
c
a+b
SaveFigure(a + b, "10_umap_Cluster_CellType", width = 14, height = 6)


saveRDS(dat, "Seurat_CellType.rds")
