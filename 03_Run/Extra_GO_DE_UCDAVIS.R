##############################################################
## Gene Ontology
#############################################################

# Gene Ontology

# View current sample names
unique(dat$orig.ident)

# Identifiy a cluster or cell type
###################
# Cluster 5: Fibroblasts
cluster5 <- subset(dat, idents = '5')
expr <- as.matrix(GetAssayData(cluster5))

# # Cell type: Interneuron
# ICC <- subset(dat, subset = CellType == "Interneurons")
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


