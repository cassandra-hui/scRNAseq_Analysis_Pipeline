#' Perform differential expression analysis on subsets of a Seurat object
#'
#' @param seurat_obj A Seurat object containing single-cell RNA sequencing data
#' @param group_col Character string specifying the column name for group comparisons (default: "group")
#' @param subset_col Character string specifying the column name for subsetting
#' @param subset_value Character string specifying the value to subset on
#' @param min_cells Integer specifying the minimum number of cells required per group (default: 3)
#'
#' @return A list containing:
#'   \itemize{
#'     \item metadata: Analysis information including timestamps and cell counts
#'     \item comparisons: Results of pairwise differential expression analyses
#'   }
#'   Or NULL if minimum cell number requirements are not met
#'
#' @examples
#' \dontrun{
#' results <- analyze_subset(seurat_obj,
#'                          group_col = "group",
#'                          subset_col = "CellType",
#'                          subset_value = "Mesothelial cells")
#' }
#'
#' @importFrom Seurat FindMarkers
#' @importFrom SeuratObject subset
#' @importFrom methods inherits
#' @import dplyr
#'
#' @export
#' 
#' @import Seurat
#' @import SeuratObject
#' @import dplyr
#' @importFrom methods inherits
NULL

# Load required packages
library(Seurat)
library(SeuratObject)
library(dplyr)


analyze_subset <- function(seurat_obj, group_col = "orig.ident", subset_col, subset_value, min_cells = 3) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  # Check if columns exist in metadata
  meta_cols <- colnames(seurat_obj@meta.data)
  if (!subset_col %in% meta_cols) {
    stop(sprintf("Column '%s' not found in metadata. Available columns: %s", 
                 subset_col, paste(meta_cols, collapse = ", ")))
  }
  if (!group_col %in% meta_cols) {
    stop(sprintf("Column '%s' not found in metadata. Available columns: %s", 
                 group_col, paste(meta_cols, collapse = ", ")))
  }
  
  
  # Get metadata and subset cells
  meta_data <- seurat_obj[[]]
  selected_cells <- rownames(meta_data)[meta_data[[subset_col]] == subset_value]
  cells <- subset(seurat_obj, cells = selected_cells)
  
  # Verify we have cells after subsetting
  if (ncol(cells) == 0) {
    stop("No cells remained after subsetting")
  }
  
  # Print basic info
  print(paste("Analyzing", subset_col, ":", subset_value))
  print(paste("Comparing groups in column:", group_col))
  
  # Get group distribution using the specified group column
  group_counts <- table(cells@meta.data[[group_col]])
  print("Group distribution:")
  print(group_counts)
  
  # Check minimum cell numbers
  if (any(group_counts < min_cells)) {
    insufficient_groups <- names(group_counts)[group_counts < min_cells]
    warning_msg <- sprintf(
      "Skipping - insufficient cells in groups: %s\nCounts: %s", 
      paste(insufficient_groups, collapse = ", "),
      paste(sprintf("%s (%d)", names(group_counts), group_counts), collapse = ", ")
    )
    warning(warning_msg)
    return(NULL)
  }
  
  # Get groups from metadata using the specified group column
  groups <- levels(cells@meta.data[[group_col]])
  if (is.null(groups)) {
    groups <- unique(cells@meta.data[[group_col]])
  }
  
  # Check if we have enough groups for comparison
  if (length(groups) < 2) {
    stop(sprintf("Need at least 2 groups in column '%s' for comparison. Found: %s", 
                 group_col, paste(groups, collapse = ", ")))
  }
  
  # Initialize results list with metadata
  de_results <- list(
    metadata = list(
      subset_col = subset_col,
      subset_value = subset_value,
      group_col = group_col,
      total_cells = ncol(cells),
      group_counts = group_counts,
      timestamp = Sys.time()
    ),
    comparisons = list()
  )
  
  # Get all pairwise comparisons
  pairs <- combn(groups, 2, simplify = FALSE)
  
  # Run comparisons with progress reporting
  total_comparisons <- length(pairs)
  for (i in seq_along(pairs)) {
    pair <- pairs[[i]]
    group1 <- pair[1]
    group2 <- pair[2]
    name <- paste(group1, "vs", group2, sep="_")
    
    print(sprintf("Running comparison %d/%d: %s", i, total_comparisons, name))
    
    tryCatch({
      de_results$comparisons[[name]] <- FindMarkers(cells,
                                                    group.by = group_col,
                                                    ident.1 = group1,
                                                    ident.2 = group2,
                                                    min.pct = 0.25)
    }, error = function(e) {
      warning(sprintf("Error in comparison %s: %s", name, e$message))
      return(NULL)
    })
  }
  
  return(de_results)
}
