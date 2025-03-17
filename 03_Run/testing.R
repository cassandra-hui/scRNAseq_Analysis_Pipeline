# Load required libraries
library(Seurat)
library(SeuratObject)
library(dplyr)



experiment.aggregate <- readRDS("~/Documents/Projects/Sal Baker/obj.RDS")


results <- analyze_subset(experiment.aggregate, 
                          subset_col = "CellType", 
                          subset_value = "Mesothelial cells")


results <- analyze_subset(experiment.aggregate, "Mesothelial cells")

# For cell types
results <- analyze_subset(experiment.aggregate, "CellType", "Mesothelial cells")


# Test in your R console to see how metadata is accessed
head(experiment.aggregate[[]])  # Should show metadata
table(experiment.aggregate[["CellType"]])  # Should show cell types


column = "CellType"
cell_type = "Mesothelial cells"

# Debug steps
# 1. Create test subset
test_cells <- subset(experiment.aggregate, column == cell_type)

# Method 2: Using subsetting with a variable
column <- "CellType"
value <- "Mesothelial cells"
cells_2 <- subset(experiment.aggregate, eval(parse(text = paste0(column, "=='", value, "'"))))
print(dim(cells_2))

# 2. Look at the group column directly from meta.data
print("Checking group structure:")
str(test_cells@meta.data$group)

# 3. Get unique groups using factor levels
print("Group levels:")
test_groups <- levels(test_cells@meta.data$group)
print(test_groups)

# 4. Check group distribution
print("Group distribution:")
print(table(test_cells@meta.data$group))

# 5. Try creating pairs from the levels
test_pairs <- combn(test_groups, 2, simplify = FALSE)
print("Pair combinations:")
print(test_pairs)

# 6. Test a single comparison
if(length(test_groups) >= 2) {
  print("Testing FindMarkers with first two groups:")
  test_de <- FindMarkers(test_cells,
                         group.by = "group",
                         ident.1 = test_groups[1],
                         ident.2 = test_groups[2],
                         min.pct = 0.25)
  print(head(test_de))
}





# Load required libraries
library(Seurat)
library(dplyr)

# Test different subsetting methods
print("Testing Seurat subsetting methods:")

# Method 1: Direct subsetting (known working)
cells_1 <- subset(experiment.aggregate, CellType == "Mesothelial cells")
print("Method 1 dimensions:")
print(dim(cells_1))

# Method 2: Using subset() with a string column name
column <- "CellType"
value <- "Mesothelial cells"
cells_2 <- subset(experiment.aggregate, subset = get(column) == value)
print("Method 2 dimensions:")
print(dim(cells_2))

# Method 3: Using WhichCells
cells_3 <- subset(experiment.aggregate, 
                  cells = WhichCells(experiment.aggregate, 
                                     expression = CellType == "Mesothelial cells"))
print("Method 3 dimensions:")
print(dim(cells_3))

# Method 4: Using metadata directly
meta_data <- experiment.aggregate[[]]
selected_cells <- rownames(meta_data)[meta_data$CellType == "Mesothelial cells"]
cells_4 <- subset(experiment.aggregate, cells = selected_cells)
print("Method 4 dimensions:")
print(dim(cells_4))


# Or with different columns
results <- analyze_subset(experiment.aggregate, "tree.ident", "3")

str(results)
View(results$Sham_vs_Partial)

experiment <- readRDS("~/Documents/Projects/Sal Baker/sham.RDS")

# Or with different columns
results <- analyze_subset(experiment, "group", "tree.ident", "3")

# Compare by tissue
results <- analyze_subset(experiment.aggregate, 
                          subset_col = "CellType", 
                          subset_value = "Mesothelial cells")


res <- analyze_subset(experiment, "tissue", "CellType", "Mesothelial cells")

