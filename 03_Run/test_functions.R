


analyze_subset <- function(seurat_obj, cell_type) {
  # Subset to specific cell type
  cells <- subset(seurat_obj, CellType == cell_type)
  
  # Print basic info
  print(paste("Analyzing:", cell_type))
  print(table(cells$group))
  
  # Get groups
  groups <- unique(cells$group)
  
  # Initialize results list
  de_results <- list()
  
  # Get all pairwise comparisons
  pairs <- combn(groups, 2, simplify = FALSE)
  
  # Run comparisons
  for (pair in pairs) {
    group1 <- pair[1]
    group2 <- pair[2]
    name <- paste(group1, "vs", group2, sep="_")
    
    print(paste("Running comparison:", name))
    
    de_results[[name]] <- FindMarkers(cells,
                                      group.by = "group",
                                      ident.1 = group1,
                                      ident.2 = group2,
                                      min.pct = 0.25)
  }
  
  return(de_results)
}

analyze_subset <- function(seurat_obj, subset_col, subset_value, min_cells = 3) {
  # Get metadata and subset cells
  meta_data <- seurat_obj[[]]
  selected_cells <- rownames(meta_data)[meta_data[[subset_col]] == subset_value]
  cells <- subset(seurat_obj, cells = selected_cells)
  
  # Print basic info
  print(paste("Analyzing", subset_col, ":", subset_value))
  print("Group distribution:")
  group_counts <- table(cells@meta.data$group)
  print(group_counts)
  
  # Check cell numbers
  if (any(group_counts < min_cells)) {
    print(paste("Skipping - some groups have fewer than", min_cells, "cells"))
    return(NULL)
  }
  
  # Get groups from metadata
  groups <- levels(cells@meta.data$group)
  
  # Initialize results list
  de_results <- list()
  
  # Get all pairwise comparisons
  pairs <- combn(groups, 2, simplify = FALSE)
  
  # Run comparisons
  for (pair in pairs) {
    group1 <- pair[1]
    group2 <- pair[2]
    name <- paste(group1, "vs", group2, sep="_")
    
    print(paste("Running comparison:", name))
    
    de_results[[name]] <- FindMarkers(cells,
                                      group.by = "group",
                                      ident.1 = group1,
                                      ident.2 = group2,
                                      min.pct = 0.25)
  }
  
  return(de_results)
}



analyze_subset <- function(seurat_obj, group_col = "group", subset_col, subset_value, min_cells = 3) {
  # Get metadata and subset cells
  meta_data <- seurat_obj[[]]
  selected_cells <- rownames(meta_data)[meta_data[[subset_col]] == subset_value]
  cells <- subset(seurat_obj, cells = selected_cells)
  
  # Print basic info
  print(paste("Analyzing", subset_col, ":", subset_value))
  print(paste("Comparing groups in column:", group_col))
  
  # Get group distribution using the specified group column
  group_counts <- table(cells@meta.data[[group_col]])
  print("Group distribution:")
  print(group_counts)
  
  # Check cell numbers
  if (any(group_counts < min_cells)) {
    print(paste("Skipping - some groups have fewer than", min_cells, "cells"))
    return(NULL)
  }
  
  # Get groups from metadata using the specified group column
  groups <- levels(cells@meta.data[[group_col]])
  if (is.null(groups)) {
    groups <- unique(cells@meta.data[[group_col]])
  }
  
  # Initialize results list
  de_results <- list()
  
  # Get all pairwise comparisons
  pairs <- combn(groups, 2, simplify = FALSE)
  
  # Run comparisons
  for (pair in pairs) {
    group1 <- pair[1]
    group2 <- pair[2]
    name <- paste(group1, "vs", group2, sep="_")
    
    print(paste("Running comparison:", name))
    
    de_results[[name]] <- FindMarkers(cells,
                                      group.by = group_col,
                                      ident.1 = group1,
                                      ident.2 = group2,
                                      min.pct = 0.25)
  }
  
  return(de_results)
}
