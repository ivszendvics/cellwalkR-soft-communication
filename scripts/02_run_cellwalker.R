# 02_run_cellwalker.R
# Runs CellWalkR on the parsed expression matrix + marker gene sets to get
# a cells x cell-types probability matrix (the "soft labels" everything
# downstream depends on).
#
# --- Install notes (confirmed against the installed cellwalker2 branch) ---
# devtools::install_github("PFPrzytycki/CellWalkR", ref = "cellwalker2")
#
# CellWalkR's DESCRIPTION pins Seurat to (>= 4.3, < 5.0). If you have Seurat
# 5.x installed (the CRAN default as of 2026), processRNASeq() below will
# fail inside Seurat::RunPCA() with "subscript out of bounds", because
# CellWalkR calls RunPCA() without explicit `features =`, relying on Seurat
# 4.x's default feature-resolution behavior. Fix by pinning BOTH packages to
# the last Seurat 4.x release pair (Seurat 5.x requires SeuratObject 5.x,
# which removes the `slot=` argument Seurat 4.x's internals use):
#   remotes::install_version("Seurat", version = "4.4.0")
#   remotes::install_version("SeuratObject", version = "4.1.4")

library(CellWalkR)
devtools::load_all(".")  # for markers_from_list()

parsed <- readRDS("data/gse103322_parsed.rds")
expr <- parsed$expr
cell_meta <- parsed$meta

# --- Marker gene sets ---
# Start from Puram et al.'s supplementary marker tables (Table S1/S5 in the
# paper) for malignant, fibroblast, T cell, B cell, myeloid, endothelial,
# and mast cell populations. Cross-check ambiguous ones against PanglaoDB.
# Format: a named list, one character vector of marker genes per cell type.
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

# --- Preprocess with Seurat (normalize, build cell-cell KNN graph) ---
# do.findMarkers = FALSE / group.col = NULL: we supply curated markers
# ourselves rather than deriving them from a Seurat clustering pass.
dataset <- processRNASeq(expr, do.findMarkers = FALSE, computeKNN = TRUE)

# --- Cell-to-label edges from curated markers ---
# computeTypeEdges() expects a markers data.frame with columns gene, cluster,
# avg_log2FC. markers_from_list() (R/utils.R) builds that from marker_list,
# assigning a placeholder avg_log2FC of 1 since curated markers don't carry
# an empirical fold-change -- log2FC.cutoff = 0 below stops that placeholder
# from filtering any markers out.
markers <- markers_from_list(marker_list)
labelEdges <- computeTypeEdges(dataset$expr_norm, markers, log2FC.cutoff = 0)

# --- Run the walk ---
# annotateCells() returns list(cellWalk, weight1); weight1 = NULL means
# CellWalkR tunes the label/cell edge-weight ratio itself via a subsample of
# `sampleDepth` cells.
walk_result <- annotateCells(dataset$cellGraph, labelEdges, weight1 = NULL,
                              sampleDepth = 3000, labelEdgeOpts = 10^seq(-5, 3, 1))
cellWalk <- walk_result[[1]]

# --- Extract a proper probability matrix ---
# cellWalk$normMat is NOT what we want here: it's a per-column Z-score-like
# normalization (see normalizeInfluence() -- (x - mean)/sd, then /max), and
# goes negative, so it can't be row-normalized into probabilities.
# The raw, non-negative random-walk-with-restart scores live in the
# label-to-cell block of the combined (cells + labels) influence matrix:
# cellWalk$infMat is (cells + labels) x (cells + labels), with the first
# `n_labels` rows/columns corresponding to labels.
n_labels <- ncol(labelEdges)
raw_scores <- cellWalk$infMat[-(seq_len(n_labels)), seq_len(n_labels)]
type_probs <- sweep(raw_scores, 1, rowSums(raw_scores), "/")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(type_probs, "data/processed/type_probs.rds")

message("Saved cells x types probability matrix to data/processed/type_probs.rds")
message("Dimensions: ", paste(dim(type_probs), collapse = " x "))
