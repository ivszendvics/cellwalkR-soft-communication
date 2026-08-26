# 02_run_cellwalker.R
# Runs CellWalkR on the parsed expression matrix + marker gene sets to get
# a cells x cell-types probability matrix (the "soft labels" everything
# downstream depends on).
#
# Install the cellwalker2 branch, not the CRAN release, to get the features
# this project relies on:
#   devtools::install_github("PFPrzytycki/CellWalkR", ref = "cellwalker2")

library(CellWalkR)
library(Seurat)

parsed <- readRDS("data/gse103322_parsed.rds")
expr <- parsed$expr
cell_meta <- parsed$meta

# --- Marker gene sets ---
# Start from Puram et al.'s supplementary marker tables (Table S1/S5 in the
# paper) for malignant, fibroblast, T cell, B cell, myeloid, endothelial,
# and mast cell populations. Cross-check ambiguous ones against PanglaoDB.
# Format CellWalkR expects: a named list, one character vector of marker
# genes per cell type.
#
# marker_list <- list(
#   Malignant   = c("EPCAM", "SFN", "KRT17", ...),
#   Fibroblast  = c("COL1A1", "COL1A2", "PDGFRB", ...),
#   Tcell       = c("CD3D", "CD3E", "CD2", ...),
#   Bcell       = c("CD79A", "MS4A1", "CD19", ...),
#   Myeloid     = c("CD68", "LYZ", "AIF1", ...),
#   Endothelial = c("PECAM1", "VWF", "CDH5", ...),
#   Mast        = c("TPSAB1", "TPSB2", "CPA3", ...)
# )
marker_list <- readRDS("data/marker_genes.rds")  # build this from the paper's supplement

# --- Build label edges from markers, per CellWalkR's expected input ---
labelEdges <- computeTypeLabelEdges(expr, marker_list)  # confirm exact fn name/signature against installed version

# --- Run the walk ---
# walkCells() (or the current cellwalker2 equivalent -- check
# ?CellWalkR::walkCells / package NEWS since the API shifted between v1 and
# v2) returns per-cell scores against each label node.
walk_result <- walkCells(
  labelEdges = labelEdges,
  cellEdges  = computeCellSimilarity(expr)  # confirm exact fn name in installed version
)

# --- Convert to a proper probability matrix ---
# walk_result$cellLabels (name may differ -- inspect str(walk_result)) is a
# cells x types score matrix from the random walk. Row-normalize to turn
# diffusion scores into a probability distribution per cell.
raw_scores <- walk_result$cellLabels
type_probs <- sweep(raw_scores, 1, rowSums(raw_scores), "/")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(type_probs, "data/processed/type_probs.rds")

message("Saved cells x types probability matrix to data/processed/type_probs.rds")
message("Dimensions: ", paste(dim(type_probs), collapse = " x "))
